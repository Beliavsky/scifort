! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! RBFInterpolator for scattered data in any positive dimension.
! Independent implementation of the augmented radial-basis collocation system.
! SciPy is a numerical reference, not a source of implementation code.
module scifort_interpolate_rbf
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only: dp
    use scifort_interpolate_kdtree, only: kd_tree, kd_tree_type, kdtree_success
    implicit none
    private

    integer, parameter, public :: rbf_status_success = 0
    integer, parameter, public :: rbf_status_invalid_input = 1
    integer, parameter, public :: rbf_status_singular = 2

    type :: rbf_core
        type(kd_tree_type) :: spatial_index ! local RBF nearest-neighbor queries
        real(dp), allocatable :: y(:,:)        ! (observations, coordinates)
        complex(dp), allocatable :: d(:,:)     ! (observations, channels)
        real(dp), allocatable :: smooth(:)     ! individual regularization parameters
        complex(dp), allocatable :: coeff(:,:) ! (observations + polynomial terms, channels)
        real(dp), allocatable :: shift(:), scale(:)
        integer, allocatable :: powers(:,:)    ! (coordinates, polynomial terms)
        character(len=32) :: kernel = 'thin_plate_spline'
        real(dp) :: epsilon_scale = 1.0_dp
        integer :: degree = 1
        integer :: neighbors = 0               ! 0 => one global fit
    end type rbf_core

    type, public :: rbf_interpolator_type
        private
        type(rbf_core) :: core
    contains
        procedure :: valid => valid_real
        procedure :: dimensions => dimensions_real
        procedure :: channels => channels_real
        procedure :: value => value_real
        procedure, private :: evaluate_one_real, evaluate_many_real
        generic :: evaluate => evaluate_one_real, evaluate_many_real
    end type rbf_interpolator_type

    type, public :: complex_rbf_interpolator_type
        private
        type(rbf_core) :: core
    contains
        procedure :: valid => valid_complex
        procedure :: dimensions => dimensions_complex
        procedure :: channels => channels_complex
        procedure :: value => value_complex
        procedure, private :: evaluate_one_complex, evaluate_many_complex
        generic :: evaluate => evaluate_one_complex, evaluate_many_complex
    end type complex_rbf_interpolator_type

    public :: rbf_interpolator
    interface rbf_interpolator
        module procedure create_real_vector, create_real_matrix
        module procedure create_complex_vector, create_complex_matrix
    end interface rbf_interpolator

contains

    pure real(dp) function nan() result(v)
        v = ieee_value(0.0_dp,ieee_quiet_nan)
    end function nan

    pure logical function core_valid(p) result(ok)
        type(rbf_core), intent(in) :: p !! interpolation data
        ok = allocated(p%y) .and. allocated(p%d) .and. allocated(p%smooth)
        if (.not. ok) return
        ok = (p%neighbors > 0 .and. p%spatial_index%valid()) .or. allocated(p%coeff)
    end function core_valid

    pure logical function valid_real(self) result(ok)
        class(rbf_interpolator_type), intent(in) :: self !! interpolator
        ok = core_valid(self%core)
    end function valid_real

    pure logical function valid_complex(self) result(ok)
        class(complex_rbf_interpolator_type), intent(in) :: self !! interpolator
        ok = core_valid(self%core)
    end function valid_complex

    pure integer function dimensions_real(self) result(n)
        class(rbf_interpolator_type), intent(in) :: self !! interpolator
        n = 0
        if (allocated(self%core%y)) n = size(self%core%y,2)
    end function dimensions_real

    pure integer function dimensions_complex(self) result(n)
        class(complex_rbf_interpolator_type), intent(in) :: self !! interpolator
        n = 0
        if (allocated(self%core%y)) n = size(self%core%y,2)
    end function dimensions_complex

    pure integer function channels_real(self) result(n)
        class(rbf_interpolator_type), intent(in) :: self !! interpolator
        n = 0
        if (allocated(self%core%d)) n = size(self%core%d,2)
    end function channels_real

    pure integer function channels_complex(self) result(n)
        class(complex_rbf_interpolator_type), intent(in) :: self !! interpolator
        n = 0
        if (allocated(self%core%d)) n = size(self%core%d,2)
    end function channels_complex

    function create_real_vector(y,d,kernel,epsilon,degree,smoothing,smoothing_values,neighbors,status) result(p)
        real(dp), intent(in) :: y(:,:)             !! (observations, coordinates)
        real(dp), intent(in) :: d(:)               !! observed values
        character(len=*), intent(in), optional :: kernel !! radial kernel name
        real(dp), intent(in), optional :: epsilon  !! radial scale, mandatory for nonhomogeneous kernels
        integer, intent(in), optional :: degree    !! polynomial degree, -1 disables polynomial
        real(dp), intent(in), optional :: smoothing !! common nonnegative smoothing parameter
        real(dp), intent(in), optional :: smoothing_values(:) !! per-node smoothing
        integer, intent(in), optional :: neighbors !! number of local nearest data sites
        integer, intent(out), optional :: status   !! 0 success, 1 invalid, 2 singular
        type(rbf_interpolator_type) :: p
        p = create_real_matrix(y,reshape(d,[size(d),1]),kernel,epsilon,degree, &
            smoothing,smoothing_values,neighbors,status)
    end function create_real_vector

    function create_real_matrix(y,d,kernel,epsilon,degree,smoothing,smoothing_values,neighbors,status) result(p)
        real(dp), intent(in) :: y(:,:)            !! (observations, coordinates)
        real(dp), intent(in) :: d(:,:)            !! (observations, channels)
        character(len=*), intent(in), optional :: kernel !! name of RBF kernel
        real(dp), intent(in), optional :: epsilon !! kernel shape parameter
        integer, intent(in), optional :: degree   !! polynomial degree, -1 to suppress
        real(dp), intent(in), optional :: smoothing !! scalar smoothing
        real(dp), intent(in), optional :: smoothing_values(:) !! per-site smoothing
        integer, intent(in), optional :: neighbors !! local neighbor count
        integer, intent(out), optional :: status !! return status
        type(rbf_interpolator_type) :: p
        integer :: code
        call initialize(p%core,y,cmplx(d,0.0_dp,kind=dp),kernel,epsilon,degree, &
            smoothing,smoothing_values,neighbors,code)
        if (present(status)) status=code
    end function create_real_matrix

    function create_complex_vector(y,d,kernel,epsilon,degree,smoothing,smoothing_values,neighbors,status) result(p)
        real(dp), intent(in) :: y(:,:) !! (observations, coordinates)
        complex(dp), intent(in) :: d(:) !! complex observations
        character(len=*), intent(in), optional :: kernel !! radial kernel name
        real(dp), intent(in), optional :: epsilon !! radial scale
        integer, intent(in), optional :: degree !! polynomial degree
        real(dp), intent(in), optional :: smoothing !! shared smoothing
        real(dp), intent(in), optional :: smoothing_values(:) !! vector smoothing
        integer, intent(in), optional :: neighbors !! nearest sites per query
        integer, intent(out), optional :: status !! error code
        type(complex_rbf_interpolator_type) :: p
        p=create_complex_matrix(y,reshape(d,[size(d),1]),kernel,epsilon,degree, &
            smoothing,smoothing_values,neighbors,status)
    end function create_complex_vector

    function create_complex_matrix(y,d,kernel,epsilon,degree,smoothing,smoothing_values,neighbors,status) result(p)
        real(dp), intent(in) :: y(:,:) !! coordinates
        complex(dp), intent(in) :: d(:,:) !! (observations, channels)
        character(len=*), intent(in), optional :: kernel !! radial kernel name
        real(dp), intent(in), optional :: epsilon !! radial scale
        integer, intent(in), optional :: degree !! polynomial degree
        real(dp), intent(in), optional :: smoothing !! scalar smoothing
        real(dp), intent(in), optional :: smoothing_values(:) !! smoothing per site
        integer, intent(in), optional :: neighbors !! number of nearest points
        integer, intent(out), optional :: status !! error status
        type(complex_rbf_interpolator_type) :: p
        integer :: code
        call initialize(p%core,y,d,kernel,epsilon,degree,smoothing,smoothing_values,neighbors,code)
        if (present(status)) status=code
    end function create_complex_matrix

    pure integer function min_degree(kernel) result(k)
        character(len=*), intent(in) :: kernel !! kernel name
        select case(kernel)
        case('linear','multiquadric')
            k=0
        case('thin_plate_spline','cubic')
            k=1
        case('quintic')
            k=2
        case('inverse_multiquadric','inverse_quadratic','gaussian', &
             'matern1_2','matern3_2','matern5_2')
            k=0
        case default
            k=-99
        end select
    end function min_degree

    pure logical function homogeneous_kernel(kernel) result(ok)
        character(len=*), intent(in) :: kernel !! kernel name
        ok = kernel=='linear' .or. kernel=='thin_plate_spline' .or. &
             kernel=='cubic' .or. kernel=='quintic'
    end function homogeneous_kernel

    subroutine initialize(p,y,d,kernel,epsilon,degree,smoothing,smoothing_values,neighbors,code)
        type(rbf_core), intent(out) :: p !! fitted interpolator
        real(dp), intent(in) :: y(:,:) !! (observation, dimension)
        complex(dp), intent(in) :: d(:,:) !! observed values, may have multiple channels
        character(len=*), intent(in), optional :: kernel !! RBF kernel name
        real(dp), intent(in), optional :: epsilon !! shape parameter
        integer, intent(in), optional :: degree !! polynomial degree
        real(dp), intent(in), optional :: smoothing !! scalar nonnegative smoothing
        real(dp), intent(in), optional :: smoothing_values(:) !! site-specific smoothing
        integer, intent(in), optional :: neighbors !! local interpolation count
        integer, intent(out) :: code !! 0 success, 1 invalid or 2 singular
        integer :: n,m,nb,nvar
        real(dp) :: lambda
        character(len=32) :: name
        code=rbf_status_invalid_input
        n=size(y,1)
        m=size(y,2)
        if (n<1 .or. m<1 .or. size(d,1)/=n .or. size(d,2)<1) return
        if (.not.all(ieee_is_finite(y))) return
        if (.not.all(ieee_is_finite(real(d,dp)))) return
        if (.not.all(ieee_is_finite(aimag(d)))) return
        name='thin_plate_spline'
        if (present(kernel)) name=kernel
        if (min_degree(trim(name)) == -99) return
        p%kernel=name
        p%degree=min_degree(trim(name))
        if (present(degree)) p%degree=degree
        if (p%degree < -1 .or. p%degree>5) return
        p%epsilon_scale=1.0_dp
        if (.not.homogeneous_kernel(trim(name)) .and. .not.present(epsilon)) return
        if (present(epsilon)) p%epsilon_scale=epsilon
        if (.not.ieee_is_finite(p%epsilon_scale) .or. p%epsilon_scale<=0.0_dp) return
        if (present(smoothing) .and. present(smoothing_values)) return
        lambda=0.0_dp
        if (present(smoothing)) lambda=smoothing
        if (.not.ieee_is_finite(lambda) .or. lambda<0.0_dp) return
        allocate(p%smooth(n))
        p%smooth=lambda
        if (present(smoothing_values)) then
            if (size(smoothing_values)/=n) return
            if (.not.all(ieee_is_finite(smoothing_values))) return
            if (any(smoothing_values<0.0_dp)) return
            p%smooth=smoothing_values
        end if
        p%neighbors=0
        if (present(neighbors)) then
            if (neighbors<1) return
            p%neighbors=min(n,neighbors)
        end if
        p%y=y
        p%d=d
        call make_powers(m,p%degree,p%powers)
        nb=size(p%powers,2)
        if (p%degree>=0 .and. nb==0) then
            call invalidate(p)
            return
        end if
        if (n<nb) then
            call invalidate(p)
            return
        end if
        if (p%neighbors>0) then
            if (p%neighbors<nb) then
                call invalidate(p)
                return
            end if
            p%spatial_index=kd_tree(y,status=code)
            if (code/=kdtree_success) then
                call invalidate(p)
                return
            end if
            code=rbf_status_success
        else
            call solve_subset(p,[(nvar,nvar=1,n)],p%coeff,p%shift,p%scale,code)
            if (code/=rbf_status_success) call invalidate(p)
        end if
    end subroutine initialize

    subroutine invalidate(p)
        type(rbf_core), intent(inout) :: p !! interpolation object to mark invalid
        if (allocated(p%y)) deallocate(p%y)
        if (allocated(p%d)) deallocate(p%d)
    end subroutine invalidate

    subroutine make_powers(dim,deg,powers)
        integer, intent(in) :: dim !! coordinate dimension
        integer, intent(in) :: deg !! maximum total polynomial degree
        integer, allocatable, intent(out) :: powers(:,:) !! multi-indices
        integer, allocatable :: candidate(:), out(:,:)
        integer :: n, index, radix, j, d, total, count
        if (deg<0) then
            allocate(powers(dim,0))
            return
        end if
        ! The public degree cap is five, and high dimensions are rejected
        ! if the number of monomials would overflow bounded storage.
        n=1
        do j=1,dim
            if (n>100000/(deg+1)) then
                allocate(powers(dim,0))
                return
            end if
            n=n*(deg+1)
        end do
        allocate(out(dim,n), candidate(dim))
        count=0
        do index=0,n-1
            radix=index
            total=0
            do j=1,dim
                d=mod(radix,deg+1)
                candidate(j)=d
                total=total+d
                radix=radix/(deg+1)
            end do
            if (total>deg) cycle
            count=count+1
            out(:,count)=candidate
        end do
        powers=out(:,:count)
    end subroutine make_powers

    pure real(dp) function kernel_value(r,kernel) result(v)
        real(dp), intent(in) :: r !! nonnegative scaled Euclidean distance
        character(len=*), intent(in) :: kernel !! kernel name
        select case(kernel)
        case('linear')
            v=-r
        case('thin_plate_spline')
            v=0.0_dp
            if (r>0.0_dp) v=r*r*log(r)
        case('cubic')
            v=r**3
        case('quintic')
            v=-r**5
        case('multiquadric')
            v=-sqrt(1.0_dp+r*r)
        case('inverse_multiquadric')
            v=1.0_dp/sqrt(1.0_dp+r*r)
        case('inverse_quadratic')
            v=1.0_dp/(1.0_dp+r*r)
        case('gaussian')
            v=exp(-r*r)
        case('matern1_2')
            v=exp(-r)
        case('matern3_2')
            v=(1.0_dp+sqrt(3.0_dp)*r)*exp(-sqrt(3.0_dp)*r)
        case('matern5_2')
            v=(1.0_dp+sqrt(5.0_dp)*r+5.0_dp*r*r/3.0_dp)*exp(-sqrt(5.0_dp)*r)
        case default
            v=0.0_dp
        end select
    end function kernel_value

    pure real(dp) function basis_term(x,powers) result(v)
        real(dp), intent(in) :: x(:) !! normalized coordinate
        integer, intent(in) :: powers(:) !! nonnegative exponents
        integer :: j
        v=1.0_dp
        do j=1,size(x)
            if (powers(j)>0) v=v*x(j)**powers(j)
        end do
    end function basis_term

    subroutine solve_subset(p,ids,coeff,shift,scale,code)
        type(rbf_core), intent(in) :: p !! RBF data and hyperparameters
        integer, intent(in) :: ids(:) !! one-based observation IDs
        complex(dp), allocatable, intent(out) :: coeff(:,:) !! fit coefficients
        real(dp), allocatable, intent(out) :: shift(:),scale(:) !! polynomial normalization
        integer, intent(out) :: code !! error flag
        integer :: n,dim,nb,m,i,j,k,c
        real(dp), allocatable :: mat(:,:),rhs(:,:),tmp(:)
        real(dp) :: r,alpha, threshold
        n=size(ids)
        dim=size(p%y,2)
        nb=size(p%powers,2)
        m=n+nb
        code=rbf_status_singular
        allocate(shift(dim),scale(dim),mat(m,m),rhs(m,2*size(p%d,2)),tmp(max(m,2*size(p%d,2))))
        shift=0.5_dp*(maxval(p%y(ids,:),dim=1)+minval(p%y(ids,:),dim=1))
        scale=0.5_dp*(maxval(p%y(ids,:),dim=1)-minval(p%y(ids,:),dim=1))
        where (scale==0.0_dp) scale=1.0_dp
        mat=0.0_dp
        rhs=0.0_dp
        do i=1,n
            rhs(i,1:size(p%d,2))=real(p%d(ids(i),:),dp)
            rhs(i,size(p%d,2)+1:)=aimag(p%d(ids(i),:))
            do j=1,i
                r=p%epsilon_scale*sqrt(sum((p%y(ids(i),:)-p%y(ids(j),:))**2))
                mat(i,j)=kernel_value(r,trim(p%kernel))
                mat(j,i)=mat(i,j)
            end do
            mat(i,i)=mat(i,i)+p%smooth(ids(i))
            do j=1,nb
                alpha=basis_term((p%y(ids(i),:)-shift)/scale,p%powers(:,j))
                mat(i,n+j)=alpha
                mat(n+j,i)=alpha
            end do
        end do
        threshold=32.0_dp*epsilon(1.0_dp)*maxval(abs(mat))
        if (threshold<=0.0_dp) return
        do k=1,m
            j=k-1+maxloc(abs(mat(k:m,k)),dim=1)
            if (abs(mat(j,k))<=threshold) return
            if (j/=k) then
                tmp(:m)=mat(k,:)
                mat(k,:)=mat(j,:)
                mat(j,:)=tmp(:m)
                tmp(:size(rhs,2))=rhs(k,:)
                rhs(k,:)=rhs(j,:)
                rhs(j,:)=tmp(:size(rhs,2))
            end if
            do i=k+1,m
                alpha=mat(i,k)/mat(k,k)
                mat(i,k)=0.0_dp
                mat(i,k+1:m)=mat(i,k+1:m)-alpha*mat(k,k+1:m)
                rhs(i,:)=rhs(i,:)-alpha*rhs(k,:)
            end do
        end do
        do i=m,1,-1
            do j=i+1,m
                rhs(i,:)=rhs(i,:)-mat(i,j)*rhs(j,:)
            end do
            rhs(i,:)=rhs(i,:)/mat(i,i)
        end do
        c=size(p%d,2)
        allocate(coeff(m,c))
        coeff=cmplx(rhs(:,1:c),rhs(:,c+1:),dp)
        code=rbf_status_success
    end subroutine solve_subset

    ! Strict total ordering for stable nearest-neighbor ties. A larger
    ! distance (or later input index when distances tie) is "worse".
    subroutine nearest_ids(p,x,ids)
        type(rbf_core), intent(in) :: p !! RBF sites and constructed spatial index
        real(dp), intent(in) :: x(:) !! query point
        integer, allocatable, intent(out) :: ids(:) !! neighbors, increasing distance and input ID
        call p%spatial_index%query(x,ids,k=p%neighbors)
    end subroutine nearest_ids

    function core_evaluate(p,x,code) result(v)
        type(rbf_core), intent(in) :: p !! global or local fitted interpolator
        real(dp), intent(in) :: x(:) !! query coordinates
        integer, intent(out) :: code !! zero success, 1 invalid, 2 singular local system
        complex(dp), allocatable :: v(:)
        complex(dp), allocatable :: coeff(:,:)
        real(dp), allocatable :: shift(:), scale(:)
        integer, allocatable :: ids(:)
        integer :: n,nb,j
        real(dp) :: r,phi
        code=rbf_status_invalid_input
        if (.not.core_valid(p)) then
            allocate(v(0))
            return
        end if
        allocate(v(size(p%d,2)))
        v=cmplx(nan(),nan(),dp)
        if (size(x)/=size(p%y,2)) return
        if (.not.all(ieee_is_finite(x))) return
        if (p%neighbors>0) then
            call nearest_ids(p,x,ids)
            call solve_subset(p,ids,coeff,shift,scale,code)
            if (code/=rbf_status_success) return
        else
            ids=[(j,j=1,size(p%y,1))]
            coeff=p%coeff
            shift=p%shift
            scale=p%scale
        end if
        v=(0.0_dp,0.0_dp)
        n=size(ids)
        nb=size(p%powers,2)
        do j=1,n
            r=p%epsilon_scale*sqrt(sum((p%y(ids(j),:)-x)**2))
            phi=kernel_value(r,trim(p%kernel))
            v=v+phi*coeff(j,:)
        end do
        do j=1,nb
            phi=basis_term((x-shift)/scale,p%powers(:,j))
            v=v+phi*coeff(n+j,:)
        end do
        code=rbf_status_success
    end function core_evaluate

    function value_real(self,x,channel,status) result(v)
        class(rbf_interpolator_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x(:) !! query coordinates
        integer, intent(in), optional :: channel !! output channel, default one
        integer, intent(out), optional :: status !! returned status
        real(dp) :: v
        complex(dp), allocatable :: ans(:)
        integer :: c,code
        ans=core_evaluate(self%core,x,code)
        c=1
        if (present(channel)) c=channel
        v=nan()
        if (code==0) then
            if (c<1 .or. c>size(ans)) then
                code=rbf_status_invalid_input
            else
                v=real(ans(c),dp)
            end if
        end if
        if (present(status)) status=code
    end function value_real

    function value_complex(self,x,channel,status) result(v)
        class(complex_rbf_interpolator_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x(:) !! query coordinates
        integer, intent(in), optional :: channel !! output channel, default one
        integer, intent(out), optional :: status !! return code
        complex(dp) :: v
        complex(dp), allocatable :: ans(:)
        integer :: c,code
        ans=core_evaluate(self%core,x,code)
        c=1
        if (present(channel)) c=channel
        v=cmplx(nan(),nan(),dp)
        if (code==0) then
            if (c<1 .or. c>size(ans)) then
                code=rbf_status_invalid_input
            else
                v=ans(c)
            end if
        end if
        if (present(status)) status=code
    end function value_complex

    function evaluate_one_real(self,x,status) result(v)
        class(rbf_interpolator_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x(:) !! query coordinates
        integer, intent(out), optional :: status !! return status
        real(dp), allocatable :: v(:)
        complex(dp), allocatable :: ans(:)
        integer :: code
        ans=core_evaluate(self%core,x,code)
        v=real(ans,dp)
        if (present(status)) status=code
    end function evaluate_one_real

    function evaluate_many_real(self,x,status) result(v)
        class(rbf_interpolator_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x(:,:) !! (queries, coordinates)
        integer, intent(out), optional :: status !! first query error, if any
        real(dp), allocatable :: v(:,:)
        real(dp), allocatable :: one(:)
        integer :: i,code,first_error
        allocate(v(size(x,1),self%channels()))
        first_error=0
        do i=1,size(x,1)
            one=self%evaluate_one_real(x(i,:),code)
            if (size(one)==size(v,2)) v(i,:)=one
            if (size(one)/=size(v,2)) v(i,:)=nan()
            if (first_error==0) first_error=code
        end do
        if (present(status)) status=first_error
    end function evaluate_many_real

    function evaluate_one_complex(self,x,status) result(v)
        class(complex_rbf_interpolator_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x(:) !! query coordinates
        integer, intent(out), optional :: status !! return status
        complex(dp), allocatable :: v(:)
        integer :: code
        v=core_evaluate(self%core,x,code)
        if (present(status)) status=code
    end function evaluate_one_complex

    function evaluate_many_complex(self,x,status) result(v)
        class(complex_rbf_interpolator_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x(:,:) !! (queries, coordinates)
        integer, intent(out), optional :: status !! first failure
        complex(dp), allocatable :: v(:,:)
        complex(dp), allocatable :: one(:)
        integer :: i,code,first_error
        allocate(v(size(x,1),self%channels()))
        first_error=0
        do i=1,size(x,1)
            one=self%evaluate_one_complex(x(i,:),code)
            if (size(one)==size(v,2)) v(i,:)=one
            if (size(one)/=size(v,2)) v(i,:)=cmplx(nan(),nan(),dp)
            if (first_error==0) first_error=code
        end do
        if (present(status)) status=first_error
    end function evaluate_many_complex

end module scifort_interpolate_rbf
