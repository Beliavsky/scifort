! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Legacy scipy.interpolate.interp1d compatibility for scalar and multichannel
! ordinates (Python axis=0). Higher-degree kinds reuse the native B-splines.
module scifort_interpolate_legacy
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_interpolate_bspline, only : bspline_array_interpolator_type, make_interp_spline
    implicit none
    private

    integer, parameter, public :: interp1d_success=0
    integer, parameter, public :: interp1d_invalid_input=1
    integer, parameter, public :: interp1d_out_of_bounds=2

    type :: interp1d_core
        real(dp), allocatable :: x(:)           ! sorted, ascending knots
        complex(dp), allocatable :: y(:,:)      ! (knots, channels)
        type(bspline_array_interpolator_type) :: spline_re, spline_im
        character(len=12) :: kind='linear'
        logical :: bounds_error=.true.
        logical :: extrapolate=.false.
        complex(dp) :: below=(0.0_dp,0.0_dp),above=(0.0_dp,0.0_dp)
    end type interp1d_core

    type, public :: interp1d_type
        private
        type(interp1d_core) :: core
    contains
        procedure :: valid => real_valid
        procedure :: channels => real_channels
        procedure :: value => real_value
        procedure, private :: real_one,real_many
        generic :: evaluate => real_one,real_many
    end type interp1d_type

    type, public :: complex_interp1d_type
        private
        type(interp1d_core) :: core
    contains
        procedure :: valid => complex_valid
        procedure :: channels => complex_channels
        procedure :: value => complex_value
        procedure, private :: complex_one,complex_many
        generic :: evaluate => complex_one,complex_many
    end type complex_interp1d_type

    public :: interp1d
    interface interp1d
        module procedure create_real_vector,create_real_matrix
        module procedure create_complex_vector,create_complex_matrix
    end interface

contains

    pure real(dp) function nan() result(v)
        v=ieee_value(0.0_dp,ieee_quiet_nan)
    end function nan

    pure logical function real_valid(self) result(ok)
        class(interp1d_type), intent(in) :: self !! interpolator
        ok=allocated(self%core%x) .and. allocated(self%core%y)
    end function real_valid

    pure logical function complex_valid(self) result(ok)
        class(complex_interp1d_type), intent(in) :: self !! interpolator
        ok=allocated(self%core%x) .and. allocated(self%core%y)
    end function complex_valid

    pure integer function real_channels(self) result(n)
        class(interp1d_type), intent(in) :: self !! interpolator
        n=0
        if (allocated(self%core%y)) n=size(self%core%y,2)
    end function real_channels

    pure integer function complex_channels(self) result(n)
        class(complex_interp1d_type), intent(in) :: self !! interpolator
        n=0
        if (allocated(self%core%y)) n=size(self%core%y,2)
    end function complex_channels

    function create_real_vector(x,y,kind,bounds_error,fill_value,fill_below,fill_above, &
        extrapolate,assume_sorted,status) result(f)
        real(dp), intent(in) :: x(:) !! input abscissas, may be unsorted
        real(dp), intent(in) :: y(:) !! observed values
        character(len=*), intent(in), optional :: kind !! linear, nearest, previous, next, ...
        logical, intent(in), optional :: bounds_error !! default true, unless extrapolating
        real(dp), intent(in), optional :: fill_value !! common fill value; default NaN
        real(dp), intent(in), optional :: fill_below !! asymmetric fill for low values
        real(dp), intent(in), optional :: fill_above !! asymmetric fill for high values
        logical, intent(in), optional :: extrapolate !! use endpoint extrapolation
        logical, intent(in), optional :: assume_sorted !! require sorted input instead of sorting
        integer, intent(out), optional :: status !! 0 success, 1 invalid data
        type(interp1d_type) :: f
        f=create_real_matrix(x,reshape(y,[size(y),1]),kind,bounds_error,fill_value, &
            fill_below,fill_above,extrapolate,assume_sorted,status)
    end function create_real_vector

    function create_real_matrix(x,y,kind,bounds_error,fill_value,fill_below,fill_above, &
        extrapolate,assume_sorted,status) result(f)
        real(dp), intent(in) :: x(:) !! input abscissas
        real(dp), intent(in) :: y(:,:) !! ordinates, shape (sample, channel)
        character(len=*), intent(in), optional :: kind !! interpolation mode
        logical, intent(in), optional :: bounds_error !! reject points outside domain
        real(dp), intent(in), optional :: fill_value !! fallback outside domain
        real(dp), intent(in), optional :: fill_below !! fallback below domain
        real(dp), intent(in), optional :: fill_above !! fallback above domain
        logical, intent(in), optional :: extrapolate !! extend interpolation beyond domain
        logical, intent(in), optional :: assume_sorted !! require x sorted, else sort input
        integer, intent(out), optional :: status !! status flag
        type(interp1d_type) :: f
        integer :: code
        call initialize(f%core,x,cmplx(y,0.0_dp,dp),kind,bounds_error, &
            cmplx_fill(fill_value),cmplx_fill(fill_below),cmplx_fill(fill_above), &
            present(fill_value),present(fill_below),present(fill_above),extrapolate,assume_sorted,code)
        if (present(status)) status=code
    end function create_real_matrix

    function create_complex_vector(x,y,kind,bounds_error,fill_value,fill_below,fill_above, &
        extrapolate,assume_sorted,status) result(f)
        real(dp), intent(in) :: x(:) !! sample abscissas
        complex(dp), intent(in) :: y(:) !! complex samples
        character(len=*), intent(in), optional :: kind !! interpolation mode
        logical, intent(in), optional :: bounds_error !! reject out-of-bounds queries
        complex(dp), intent(in), optional :: fill_value !! common out-of-bounds fill
        complex(dp), intent(in), optional :: fill_below !! low fill
        complex(dp), intent(in), optional :: fill_above !! high fill
        logical, intent(in), optional :: extrapolate !! enable extrapolation
        logical, intent(in), optional :: assume_sorted !! require sorted nodes
        integer, intent(out), optional :: status !! status code
        type(complex_interp1d_type) :: f
        f=create_complex_matrix(x,reshape(y,[size(y),1]),kind,bounds_error, &
            fill_value,fill_below,fill_above,extrapolate,assume_sorted,status)
    end function create_complex_vector

    function create_complex_matrix(x,y,kind,bounds_error,fill_value,fill_below,fill_above, &
        extrapolate,assume_sorted,status) result(f)
        real(dp), intent(in) :: x(:) !! sample abscissas
        complex(dp), intent(in) :: y(:,:) !! complex ordinates
        character(len=*), intent(in), optional :: kind !! interpolation mode
        logical, intent(in), optional :: bounds_error !! reject out-of-bounds queries
        complex(dp), intent(in), optional :: fill_value !! common fill
        complex(dp), intent(in), optional :: fill_below !! below-domain fill
        complex(dp), intent(in), optional :: fill_above !! above-domain fill
        logical, intent(in), optional :: extrapolate !! interpolate outside data range
        logical, intent(in), optional :: assume_sorted !! require already sorted nodes
        integer, intent(out), optional :: status !! error status
        type(complex_interp1d_type) :: f
        integer :: code
        call initialize(f%core,x,y,kind,bounds_error, &
            complex_fill(fill_value),complex_fill(fill_below),complex_fill(fill_above), &
            present(fill_value),present(fill_below),present(fill_above),extrapolate,assume_sorted,code)
        if (present(status)) status=code
    end function create_complex_matrix

    pure complex(dp) function cmplx_fill(v) result(z)
        real(dp), intent(in), optional :: v !! optional real fill value
        z=cmplx(nan(),0.0_dp,dp)
        if (present(v)) z=cmplx(v,0.0_dp,dp)
    end function cmplx_fill

    pure complex(dp) function complex_fill(v) result(z)
        complex(dp), intent(in), optional :: v !! optional complex fill value
        z=cmplx(nan(),0.0_dp,dp)
        if (present(v)) z=v
    end function complex_fill

    subroutine initialize(p,x,y,kind,bounds_error,fill_value,fill_below,fill_above, &
        got_value,got_below,got_above,extrapolate,assume_sorted,code)
        type(interp1d_core), intent(out) :: p !! fitted legacy interpolator
        real(dp), intent(in) :: x(:) !! sample coordinates
        complex(dp), intent(in) :: y(:,:) !! samples, (knot, channel)
        character(len=*), intent(in), optional :: kind !! interpolation mode
        logical, intent(in), optional :: bounds_error !! reject outside domain
        complex(dp), intent(in) :: fill_value !! default fill if provided
        complex(dp), intent(in) :: fill_below !! low-side fill value
        complex(dp), intent(in) :: fill_above !! high-side fill value
        logical, intent(in) :: got_value !! whether fill_value was passed
        logical, intent(in) :: got_below !! whether low fill was passed
        logical, intent(in) :: got_above !! whether high fill was passed
        logical, intent(in), optional :: extrapolate !! extrapolate at ends
        logical, intent(in), optional :: assume_sorted !! require sorted input
        integer, intent(out) :: code !! error code
        integer :: n,degree,i,j,k,idx
        real(dp) :: tx
        complex(dp), allocatable :: ty(:)
        logical :: sorted
        character(len=12) :: mode
        code=interp1d_invalid_input
        n=size(x)
        if (n<2 .or. size(y,1)/=n .or. size(y,2)<1) return
        if (.not.all(ieee_is_finite(x))) return
        if (.not.all(ieee_is_finite(real(y,dp))) .or. .not.all(ieee_is_finite(aimag(y)))) return
        mode='linear'
        if (present(kind)) mode=kind
        select case(mode)
        case('linear','nearest','nearest-up','previous','next','zero','slinear','quadratic','cubic')
        case default
            return
        end select
        degree=-1
        select case(mode)
        case('slinear')
            degree=1
        case('quadratic')
            degree=2
        case('cubic')
            degree=3
        end select
        if (degree>=0 .and. n<degree+1) return
        sorted=.false.
        if (present(assume_sorted)) sorted=assume_sorted
        if (sorted .and. any(x(2:)<=x(:n-1))) return
        p%x=x
        p%y=y
        if (.not.sorted) then
            allocate(ty(size(y,2)))
            do i=2,n
                tx=p%x(i)
                ty=p%y(i,:)
                j=i-1
                do while (j>=1)
                    if (p%x(j)<=tx) exit
                    p%x(j+1)=p%x(j)
                    p%y(j+1,:)=p%y(j,:)
                    j=j-1
                end do
                p%x(j+1)=tx
                p%y(j+1,:)=ty
            end do
        end if
        if (any(p%x(2:)<=p%x(:n-1))) then
            deallocate(p%x,p%y)
            return
        end if
        p%kind=mode
        p%extrapolate=.false.
        if (present(extrapolate)) p%extrapolate=extrapolate
        p%bounds_error=.not.p%extrapolate
        if (present(bounds_error)) p%bounds_error=bounds_error
        if (p%bounds_error .and. p%extrapolate) then
            deallocate(p%x,p%y)
            return
        end if
        p%below=cmplx(nan(),0.0_dp,dp)
        p%above=cmplx(nan(),0.0_dp,dp)
        if (got_value) then
            p%below=fill_value
            p%above=fill_value
        end if
        if (got_below) p%below=fill_below
        if (got_above) p%above=fill_above
        if (degree>=0) then
            p%spline_re=make_interp_spline(p%x,real(p%y,dp),k=degree,status=code)
            if (code/=0) then
                deallocate(p%x,p%y)
                return
            end if
            p%spline_im=make_interp_spline(p%x,aimag(p%y),k=degree,status=code)
            if (code/=0) then
                deallocate(p%x,p%y)
                return
            end if
        end if
        code=interp1d_success
    end subroutine initialize

    function evaluate_core(p,q,code) result(v)
        type(interp1d_core), intent(in) :: p !! constructed interpolator
        real(dp), intent(in) :: q !! query coordinate
        integer, intent(out) :: code !! status
        complex(dp), allocatable :: v(:)
        integer :: i,n,lo,hi,m
        real(dp) :: t
        code=interp1d_invalid_input
        if (.not.allocated(p%x)) then
            allocate(v(0))
            return
        end if
        allocate(v(size(p%y,2)))
        v=cmplx(nan(),0.0_dp,dp)
        if (.not.ieee_is_finite(q)) return
        n=size(p%x)
        code=interp1d_success
        if (q<p%x(1)) then
            if (p%bounds_error) then
                code=interp1d_out_of_bounds
                return
            end if
            if (.not.p%extrapolate) then
                v=p%below
                return
            end if
            if (p%kind=='previous') return
        else if (q>p%x(n)) then
            if (p%bounds_error) then
                code=interp1d_out_of_bounds
                return
            end if
            if (.not.p%extrapolate) then
                v=p%above
                return
            end if
            if (p%kind=='next') return
        end if
        lo=1
        hi=n
        do while (hi-lo>1)
            m=(lo+hi)/2
            if (p%x(m)<=q) then
                lo=m
            else
                hi=m
            end if
        end do
        if (q>=p%x(n)) lo=n
        select case(p%kind)
        case('nearest','nearest-up')
            lo=min(n-1,lo)
            t=0.5_dp*(p%x(lo)+p%x(lo+1))
            if (p%kind=='nearest') then
                if (q>t) lo=lo+1
            else
                if (q>=t) lo=lo+1
            end if
            v=p%y(lo,:)
        case('previous','zero')
            v=p%y(lo,:)
        case('next')
            if (q>p%x(lo) .and. lo<n) lo=lo+1
            v=p%y(lo,:)
        case('linear')
            lo=min(n-1,lo)
            t=(q-p%x(lo))/(p%x(lo+1)-p%x(lo))
            v=(1.0_dp-t)*p%y(lo,:)+t*p%y(lo+1,:)
        case('slinear','quadratic','cubic')
            v=cmplx(p%spline_re%evaluate(q),p%spline_im%evaluate(q),dp)
        end select
    end function evaluate_core

    function real_value(self,x,channel,status) result(v)
        class(interp1d_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x !! query coordinate
        integer, intent(in), optional :: channel !! output channel, default 1
        integer, intent(out), optional :: status !! error status
        real(dp) :: v
        complex(dp), allocatable :: out(:)
        integer :: code,c
        out=evaluate_core(self%core,x,code)
        c=1
        if (present(channel)) c=channel
        v=nan()
        if (code==0) then
            if (c<1 .or. c>size(out)) then
                code=interp1d_invalid_input
            else
                v=real(out(c),dp)
            end if
        end if
        if (present(status)) status=code
    end function real_value

    function complex_value(self,x,channel,status) result(v)
        class(complex_interp1d_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x !! query coordinate
        integer, intent(in), optional :: channel !! output channel, default 1
        integer, intent(out), optional :: status !! return status
        complex(dp) :: v
        complex(dp), allocatable :: out(:)
        integer :: code,c
        out=evaluate_core(self%core,x,code)
        c=1
        if (present(channel)) c=channel
        v=cmplx(nan(),nan(),dp)
        if (code==0) then
            if (c<1 .or. c>size(out)) then
                code=interp1d_invalid_input
            else
                v=out(c)
            end if
        end if
        if (present(status)) status=code
    end function complex_value

    function real_one(self,x,status) result(v)
        class(interp1d_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x !! query coordinate
        integer, intent(out), optional :: status !! return status
        real(dp), allocatable :: v(:)
        complex(dp), allocatable :: ans(:)
        integer :: code
        ans=evaluate_core(self%core,x,code)
        v=real(ans,dp)
        if (present(status)) status=code
    end function real_one

    function complex_one(self,x,status) result(v)
        class(complex_interp1d_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x !! query coordinate
        integer, intent(out), optional :: status !! return status
        complex(dp), allocatable :: v(:)
        integer :: code
        v=evaluate_core(self%core,x,code)
        if (present(status)) status=code
    end function complex_one

    function real_many(self,x,status) result(v)
        class(interp1d_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x(:) !! query coordinates
        integer, intent(out), optional :: status !! first error if any
        real(dp), allocatable :: v(:,:)
        complex(dp), allocatable :: ans(:)
        integer :: i,code,first_error
        allocate(v(size(x),self%channels()))
        first_error=0
        do i=1,size(x)
            ans=evaluate_core(self%core,x(i),code)
            if (size(ans)==size(v,2)) v(i,:)=real(ans,dp)
            if (size(ans)/=size(v,2)) v(i,:)=nan()
            if (first_error==0) first_error=code
        end do
        if (present(status)) status=first_error
    end function real_many

    function complex_many(self,x,status) result(v)
        class(complex_interp1d_type), intent(in) :: self !! interpolator
        real(dp), intent(in) :: x(:) !! query coordinates
        integer, intent(out), optional :: status !! first error if any
        complex(dp), allocatable :: v(:,:)
        complex(dp), allocatable :: ans(:)
        integer :: i,code,first_error
        allocate(v(size(x),self%channels()))
        first_error=0
        do i=1,size(x)
            ans=evaluate_core(self%core,x(i),code)
            if (size(ans)==size(v,2)) v(i,:)=ans
            if (size(ans)/=size(v,2)) v(i,:)=cmplx(nan(),nan(),dp)
            if (first_error==0) first_error=code
        end do
        if (present(status)) status=first_error
    end function complex_many

end module scifort_interpolate_legacy
