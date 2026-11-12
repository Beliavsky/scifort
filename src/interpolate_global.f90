! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Confluent Newton/Hermite (Krogh) and Floater--Hormann interpolation.
! Original implementation from the published divided-difference and rational
! barycentric formulas; SciPy is used only as an independent test reference.
module scifort_interpolate_global
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: global_interp_success = 0
    integer, parameter, public :: global_interp_invalid_input = 1

    type :: krogh_core
        real(dp), allocatable :: nodes(:)
        complex(dp), allocatable :: coefficients(:,:)
    end type krogh_core

    type :: floater_core
        real(dp), allocatable :: nodes(:), weights(:)
        complex(dp), allocatable :: ordinates(:,:)
        integer :: degree = -1
    end type floater_core

    type, public :: krogh_interpolator_type
        private
        type(krogh_core) :: core
    contains
        procedure :: valid => real_krogh_valid
        procedure :: channels => real_krogh_channels
        procedure :: nodes => real_krogh_nodes
        procedure :: value => real_krogh_value
        procedure :: values => real_krogh_values
        procedure :: derivative => real_krogh_derivative
        procedure :: derivatives => real_krogh_derivatives
    end type krogh_interpolator_type

    type, public :: complex_krogh_interpolator_type
        private
        type(krogh_core) :: core
    contains
        procedure :: valid => complex_krogh_valid
        procedure :: channels => complex_krogh_channels
        procedure :: nodes => complex_krogh_nodes
        procedure :: value => complex_krogh_value
        procedure :: values => complex_krogh_values
        procedure :: derivative => complex_krogh_derivative
        procedure :: derivatives => complex_krogh_derivatives
    end type complex_krogh_interpolator_type

    type, public :: floater_hormann_interpolator_type
        private
        type(floater_core) :: core
    contains
        procedure :: valid => real_floater_valid
        procedure :: channels => real_floater_channels
        procedure :: nodes => real_floater_nodes
        procedure :: weights => real_floater_weights
        procedure :: value => real_floater_value
        procedure :: values => real_floater_values
        procedure :: degree => real_floater_degree
    end type floater_hormann_interpolator_type

    type, public :: complex_floater_hormann_interpolator_type
        private
        type(floater_core) :: core
    contains
        procedure :: valid => complex_floater_valid
        procedure :: channels => complex_floater_channels
        procedure :: nodes => complex_floater_nodes
        procedure :: weights => complex_floater_weights
        procedure :: value => complex_floater_value
        procedure :: values => complex_floater_values
        procedure :: degree => complex_floater_degree
    end type complex_floater_hormann_interpolator_type

    public :: krogh_interpolator, krogh_interpolate, floater_hormann_interpolator
    interface krogh_interpolator
        module procedure real_krogh_vector_constructor
        module procedure real_krogh_matrix_constructor
        module procedure complex_krogh_vector_constructor
        module procedure complex_krogh_matrix_constructor
    end interface
    interface floater_hormann_interpolator
        module procedure real_floater_vector_constructor
        module procedure real_floater_matrix_constructor
        module procedure complex_floater_vector_constructor
        module procedure complex_floater_matrix_constructor
    end interface
    interface krogh_interpolate
        module procedure real_krogh_convenience
        module procedure complex_krogh_convenience
    end interface
contains

    pure logical function finite_complex(y)
        complex(dp), intent(in) :: y(:,:) !! ordinates to validate
        finite_complex = all(ieee_is_finite(real(y,dp))) .and. all(ieee_is_finite(aimag(y)))
    end function finite_complex

    ! The repeated abscissae for a Hermite condition must appear in one
    ! contiguous group, with y entries holding successive derivatives.
    ! Different groups may occur in arbitrary order.
    subroutine build_krogh(x,y,p,status)
        real(dp), intent(in) :: x(:) !! real nodes, repeated consecutively for derivatives
        complex(dp), intent(in) :: y(:,:) !! successive derivative values at repeated nodes
        type(krogh_core), intent(out) :: p !! fitted confluent Newton representation
        integer, intent(out) :: status !! zero if the input is valid
        complex(dp), allocatable :: tab(:,:,:)
        integer, allocatable :: starts(:)
        integer :: n, i, j, k, nchan
        real(dp) :: factorial
        status = global_interp_invalid_input
        n = size(x)
        nchan = size(y,2)
        if (n < 1 .or. nchan < 1 .or. size(y,1) /= n) return
        if (.not. all(ieee_is_finite(x)) .or. .not. finite_complex(y)) return
        allocate(starts(n))
        do i=1,n
            starts(i) = i
            if (i > 1) then
                if (x(i) == x(i-1)) starts(i) = starts(i-1)
            end if
            ! A node may not appear again after its derivative block.
            if (starts(i) == i) then
                if (any(x(:i-1) == x(i))) return
            end if
        end do
        allocate(tab(n,n,nchan))
        tab = (0.0_dp,0.0_dp)
        do i=1,n
            tab(i,1,:) = y(starts(i),:)
        end do
        factorial = 1.0_dp
        do j=2,n
            factorial = factorial * real(j-1,dp)
            do i=j,n
                if (x(i) == x(i-j+1)) then
                    ! All nodes in this divided-difference block coincide.
                    tab(i,j,:) = y(starts(i)+j-1,:)/factorial
                else
                    tab(i,j,:) = (tab(i,j-1,:)-tab(i-1,j-1,:))/(x(i)-x(i-j+1))
                end if
            end do
        end do
        p%nodes=x
        allocate(p%coefficients(n,nchan))
        do i=1,n
            p%coefficients(i,:)=tab(i,i,:)
        end do
        status = global_interp_success
    end subroutine build_krogh

    ! Evaluate the Newton representation and an arbitrary derivative by a
    ! differentiated Horner recurrence, never expanding monomial coefficients.
    function krogh_at(p,q,nu) result(ans)
        type(krogh_core), intent(in) :: p !! fitted Newton coefficients
        real(dp), intent(in) :: q !! evaluation coordinate
        integer, intent(in) :: nu !! nonnegative derivative order
        complex(dp) :: ans(size(p%coefficients,2))
        complex(dp) :: v(0:min(max(0,nu),size(p%nodes)-1),size(p%coefficients,2))
        integer :: i,j,n
        ans=(0.0_dp,0.0_dp)
        n=size(p%nodes)
        if (nu<0 .or. nu>=n .or. .not.ieee_is_finite(q)) return
        v=(0.0_dp,0.0_dp)
        v(0,:)=p%coefficients(n,:)
        do i=n-1,1,-1
            do j=nu,1,-1
                v(j,:)=(q-p%nodes(i))*v(j,:) + real(j,dp)*v(j-1,:)
            end do
            v(0,:)=(q-p%nodes(i))*v(0,:)+p%coefficients(i,:)
        end do
        ans=v(nu,:)
    end function krogh_at

    ! Produce derivative orders 0 .. count-1 in one Horner traversal.
    ! The entries for orders at least the polynomial degree+1 are zero.
    function krogh_all_at(p,q,count) result(v)
        type(krogh_core), intent(in) :: p !! fitted Newton coefficients
        real(dp), intent(in) :: q !! evaluation coordinate
        integer, intent(in) :: count !! number of orders, including order zero
        complex(dp) :: v(0:count-1,size(p%coefficients,2))
        integer :: i,j,n,limit
        v=(0.0_dp,0.0_dp)
        if (.not.ieee_is_finite(q)) return
        n=size(p%nodes)
        limit=min(count-1,n-1)
        v(0,:)=p%coefficients(n,:)
        do i=n-1,1,-1
            do j=limit,1,-1
                v(j,:)=(q-p%nodes(i))*v(j,:) + real(j,dp)*v(j-1,:)
            end do
            v(0,:)=(q-p%nodes(i))*v(0,:)+p%coefficients(i,:)
        end do
    end function krogh_all_at

    ! Magnitudes of the Floater--Hormann weights are positive sums of
    ! inverse distance products. Form the terms in log space, then globally
    ! normalize. The alternating sign follows the node index. Multiplying
    ! all weights by a common constant does not change the interpolant.
    subroutine build_floater(x,y,d,p,status)
        real(dp), intent(in) :: x(:) !! strictly monotonic distinct real nodes
        complex(dp), intent(in) :: y(:,:) !! finite node values (nodes,channels)
        integer, intent(in) :: d !! local polynomial degree, 0 <= d < nodes
        type(floater_core), intent(out) :: p !! fitted rational interpolator
        integer, intent(out) :: status !! zero on successful construction
        integer :: n,i,j,k,lo,hi
        real(dp) :: span, logsmax, vmax, s
        real(dp), allocatable :: log_weights(:), terms(:)
        status=global_interp_invalid_input
        n=size(x)
        if (n<1 .or. d<0 .or. d>=n .or. size(y,1)/=n .or. size(y,2)<1) return
        if (.not. all(ieee_is_finite(x)) .or. .not.finite_complex(y)) return
        if (n>1) then
            if (.not.(all(x(2:)>x(:n-1)) .or. all(x(2:)<x(:n-1)))) return
        end if
        allocate(log_weights(n))
        span = maxval(x)-minval(x)
        if (n==1) span=1.0_dp
        if (.not.ieee_is_finite(span) .or. span<=0.0_dp) return
        do i=1,n
            lo=max(1,i-d)
            hi=min(i,n-d)
            allocate(terms(hi-lo+1))
            do k=lo,hi
                s=0.0_dp
                do j=k,k+d
                    if (j==i) cycle
                    s=s-log(abs((x(i)-x(j))/span))
                end do
                terms(k-lo+1)=s
            end do
            vmax=maxval(terms)
            log_weights(i)=vmax+log(sum(exp(terms-vmax)))
            deallocate(terms)
        end do
        logsmax=maxval(log_weights)
        allocate(p%weights(n))
        p%weights=exp(log_weights-logsmax)
        if (any(p%weights==0.0_dp)) return
        do i=2,n,2
            p%weights(i)=-p%weights(i)
        end do
        p%nodes=x
        p%ordinates=y
        p%degree=d
        status=global_interp_success
    end subroutine build_floater

    function floater_at(p,q) result(ans)
        type(floater_core), intent(in) :: p !! fitted Floater--Hormann object
        real(dp), intent(in) :: q !! finite query coordinate
        complex(dp) :: ans(size(p%ordinates,2))
        real(dp) :: terms(size(p%nodes)), scale, denominator
        integer :: i
        ans=(0.0_dp,0.0_dp)
        if (.not.ieee_is_finite(q)) return
        do i=1,size(p%nodes)
            if (q==p%nodes(i)) then
                ans=p%ordinates(i,:)
                return
            end if
        end do
        scale=max(1.0_dp,abs(q),maxval(abs(p%nodes)))
        terms=p%weights/((q-p%nodes)/scale)
        if (maxval(abs(terms))==0.0_dp) return
        terms=terms/maxval(abs(terms))
        denominator=sum(terms)
        if (denominator==0.0_dp) return
        do i=1,size(p%nodes)
            ans=ans+terms(i)*p%ordinates(i,:)
        end do
        ans=ans/denominator
    end function floater_at
    function real_krogh_vector_constructor(x,y,status) result(p)
        real(dp), intent(in) :: x(:) !! real interpolation nodes
        real(dp), intent(in) :: y(:) !! finite nodal values, derivatives for repeated Krogh nodes
        integer, intent(out), optional :: status !! zero for a valid interpolator
        type(krogh_interpolator_type) :: p
        integer :: code
        code=global_interp_invalid_input
        if (size(x)==size(y)) then
            call build_krogh(x,reshape(cmplx(y,0.0_dp,dp),[size(y),1]),p%core,code)
        end if
        if (present(status)) status=code
    end function real_krogh_vector_constructor

    function real_krogh_matrix_constructor(x,y,status) result(p)
        real(dp), intent(in) :: x(:) !! real interpolation nodes
        real(dp), intent(in) :: y(:,:) !! finite nodal values, derivatives for repeated Krogh nodes
        integer, intent(out), optional :: status !! zero for a valid interpolator
        type(krogh_interpolator_type) :: p
        integer :: code
        code=global_interp_invalid_input
        if (size(x)==size(y,1)) then
            call build_krogh(x,cmplx(y,0.0_dp,dp),p%core,code)
        end if
        if (present(status)) status=code
    end function real_krogh_matrix_constructor

    pure logical function real_krogh_valid(self)
        class(krogh_interpolator_type), intent(in) :: self !! object to validate
        real_krogh_valid=allocated(self%core%nodes) .and. allocated(self%core%coefficients)
    end function real_krogh_valid

    pure integer function real_krogh_channels(self)
        class(krogh_interpolator_type), intent(in) :: self !! interpolator to inspect
        real_krogh_channels=0
        if (self%valid()) real_krogh_channels=size(self%core%coefficients,2)
    end function real_krogh_channels

    function real_krogh_nodes(self) result(v)
        class(krogh_interpolator_type), intent(in) :: self !! fitted interpolation object
        real(dp), allocatable :: v(:)
        if (self%valid()) then
            v=self%core%nodes
        else
            allocate(v(0))
        end if
    end function real_krogh_nodes

    function real_krogh_value(self,q,channel,nu) result(v)
        class(krogh_interpolator_type), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: q !! query coordinate
        integer, intent(in), optional :: channel !! channel index, default one
        integer, intent(in), optional :: nu !! derivative order, default zero
        real(dp) :: v
        integer :: index, order
        v=0.0_dp
        if (.not.self%valid()) return
        index=1
        if (present(channel)) index=channel
        if (index<1 .or. index>self%channels()) return
        order=0
        if (present(nu)) order=nu
        block
            complex(dp) :: work(self%channels())
            work=krogh_at(self%core,q,order)
            v=real(work(index),dp)
        end block
    end function real_krogh_value

    function real_krogh_values(self,q,nu) result(v)
        class(krogh_interpolator_type), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: q(:) !! query coordinates
        integer, intent(in), optional :: nu !! nonnegative derivative order (default zero)
        real(dp), allocatable :: v(:,:)
        integer :: i, order
        if (.not.self%valid()) then
            allocate(v(0,0))
            return
        end if
        allocate(v(size(q),self%channels()))
        order=0
        if (present(nu)) order=nu
        do i=1,size(q)
            block
                complex(dp) :: work(self%channels())
                work=krogh_at(self%core,q(i),order)
                v(i,:)=real(work,dp)
            end block
        end do
    end function real_krogh_values

    function real_krogh_derivative(self,q,nu,channel) result(v)
        class(krogh_interpolator_type), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: q !! query coordinate
        integer, intent(in), optional :: nu !! derivative order (default 1)
        integer, intent(in), optional :: channel !! output channel, default one
        real(dp) :: v
        integer :: order
        order=1
        if (present(nu)) order=nu
        v=self%value(q,channel=channel,nu=order)
    end function real_krogh_derivative

    function complex_krogh_vector_constructor(x,y,status) result(p)
        real(dp), intent(in) :: x(:) !! real interpolation nodes
        complex(dp), intent(in) :: y(:) !! finite nodal values, derivatives for repeated Krogh nodes
        integer, intent(out), optional :: status !! zero for a valid interpolator
        type(complex_krogh_interpolator_type) :: p
        integer :: code
        code=global_interp_invalid_input
        if (size(x)==size(y)) then
            call build_krogh(x,reshape(y,[size(y),1]),p%core,code)
        end if
        if (present(status)) status=code
    end function complex_krogh_vector_constructor

    function complex_krogh_matrix_constructor(x,y,status) result(p)
        real(dp), intent(in) :: x(:) !! real interpolation nodes
        complex(dp), intent(in) :: y(:,:) !! finite nodal values, derivatives for repeated Krogh nodes
        integer, intent(out), optional :: status !! zero for a valid interpolator
        type(complex_krogh_interpolator_type) :: p
        integer :: code
        code=global_interp_invalid_input
        if (size(x)==size(y,1)) then
            call build_krogh(x,y,p%core,code)
        end if
        if (present(status)) status=code
    end function complex_krogh_matrix_constructor

    pure logical function complex_krogh_valid(self)
        class(complex_krogh_interpolator_type), intent(in) :: self !! object to validate
        complex_krogh_valid=allocated(self%core%nodes) .and. allocated(self%core%coefficients)
    end function complex_krogh_valid

    pure integer function complex_krogh_channels(self)
        class(complex_krogh_interpolator_type), intent(in) :: self !! interpolator to inspect
        complex_krogh_channels=0
        if (self%valid()) complex_krogh_channels=size(self%core%coefficients,2)
    end function complex_krogh_channels

    function complex_krogh_nodes(self) result(v)
        class(complex_krogh_interpolator_type), intent(in) :: self !! fitted interpolation object
        real(dp), allocatable :: v(:)
        if (self%valid()) then
            v=self%core%nodes
        else
            allocate(v(0))
        end if
    end function complex_krogh_nodes

    function complex_krogh_value(self,q,channel,nu) result(v)
        class(complex_krogh_interpolator_type), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: q !! query coordinate
        integer, intent(in), optional :: channel !! channel index, default one
        integer, intent(in), optional :: nu !! derivative order, default zero
        complex(dp) :: v
        integer :: index, order
        v=0.0_dp
        if (.not.self%valid()) return
        index=1
        if (present(channel)) index=channel
        if (index<1 .or. index>self%channels()) return
        order=0
        if (present(nu)) order=nu
        block
            complex(dp) :: work(self%channels())
            work=krogh_at(self%core,q,order)
            v=work(index)
        end block
    end function complex_krogh_value

    function complex_krogh_values(self,q,nu) result(v)
        class(complex_krogh_interpolator_type), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: q(:) !! query coordinates
        integer, intent(in), optional :: nu !! nonnegative derivative order (default zero)
        complex(dp), allocatable :: v(:,:)
        integer :: i, order
        if (.not.self%valid()) then
            allocate(v(0,0))
            return
        end if
        allocate(v(size(q),self%channels()))
        order=0
        if (present(nu)) order=nu
        do i=1,size(q)
            block
                complex(dp) :: work(self%channels())
                work=krogh_at(self%core,q(i),order)
                v(i,:)=work
            end block
        end do
    end function complex_krogh_values

    function complex_krogh_derivative(self,q,nu,channel) result(v)
        class(complex_krogh_interpolator_type), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: q !! query coordinate
        integer, intent(in), optional :: nu !! derivative order (default 1)
        integer, intent(in), optional :: channel !! output channel, default one
        complex(dp) :: v
        integer :: order
        order=1
        if (present(nu)) order=nu
        v=self%value(q,channel=channel,nu=order)
    end function complex_krogh_derivative

    function real_floater_vector_constructor(x,y,d,status) result(p)
        real(dp), intent(in) :: x(:) !! real interpolation nodes
        real(dp), intent(in) :: y(:) !! finite nodal values, derivatives for repeated Krogh nodes
        integer, intent(in), optional :: d !! polynomial order (default 3), 0 <= d < number of nodes
        integer, intent(out), optional :: status !! zero for a valid interpolator
        type(floater_hormann_interpolator_type) :: p
        integer :: code
        integer :: order
        order=3
        if (present(d)) order=d
        code=global_interp_invalid_input
        if (size(x)==size(y)) then
            call build_floater(x,reshape(cmplx(y,0.0_dp,dp),[size(y),1]),order,p%core,code)
        end if
        if (present(status)) status=code
    end function real_floater_vector_constructor

    function real_floater_matrix_constructor(x,y,d,status) result(p)
        real(dp), intent(in) :: x(:) !! real interpolation nodes
        real(dp), intent(in) :: y(:,:) !! finite nodal values, derivatives for repeated Krogh nodes
        integer, intent(in), optional :: d !! polynomial order (default 3), 0 <= d < number of nodes
        integer, intent(out), optional :: status !! zero for a valid interpolator
        type(floater_hormann_interpolator_type) :: p
        integer :: code
        integer :: order
        order=3
        if (present(d)) order=d
        code=global_interp_invalid_input
        if (size(x)==size(y,1)) then
            call build_floater(x,cmplx(y,0.0_dp,dp),order,p%core,code)
        end if
        if (present(status)) status=code
    end function real_floater_matrix_constructor

    pure logical function real_floater_valid(self)
        class(floater_hormann_interpolator_type), intent(in) :: self !! object to validate
        real_floater_valid=allocated(self%core%nodes) .and. allocated(self%core%ordinates)
    end function real_floater_valid

    pure integer function real_floater_channels(self)
        class(floater_hormann_interpolator_type), intent(in) :: self !! interpolator to inspect
        real_floater_channels=0
        if (self%valid()) real_floater_channels=size(self%core%ordinates,2)
    end function real_floater_channels

    function real_floater_nodes(self) result(v)
        class(floater_hormann_interpolator_type), intent(in) :: self !! fitted interpolation object
        real(dp), allocatable :: v(:)
        if (self%valid()) then
            v=self%core%nodes
        else
            allocate(v(0))
        end if
    end function real_floater_nodes

    function real_floater_weights(self) result(v)
        class(floater_hormann_interpolator_type), intent(in) :: self !! rational interpolation object
        real(dp), allocatable :: v(:)
        if (self%valid()) then
            v=self%core%weights
        else
            allocate(v(0))
        end if
    end function real_floater_weights

    pure integer function real_floater_degree(self)
        class(floater_hormann_interpolator_type), intent(in) :: self !! rational interpolation object
        real_floater_degree=self%core%degree
    end function real_floater_degree

    function real_floater_value(self,q,channel) result(v)
        class(floater_hormann_interpolator_type), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: q !! query coordinate
        integer, intent(in), optional :: channel !! channel index, default one
        real(dp) :: v
        integer :: index
        v=0.0_dp
        if (.not.self%valid()) return
        index=1
        if (present(channel)) index=channel
        if (index<1 .or. index>self%channels()) return
        block
            complex(dp) :: work(self%channels())
            work=floater_at(self%core,q)
            v=real(work(index),dp)
        end block
    end function real_floater_value

    function real_floater_values(self,q) result(v)
        class(floater_hormann_interpolator_type), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: q(:) !! query coordinates
        real(dp), allocatable :: v(:,:)
        integer :: i
        if (.not.self%valid()) then
            allocate(v(0,0))
            return
        end if
        allocate(v(size(q),self%channels()))
        do i=1,size(q)
            block
                complex(dp) :: work(self%channels())
                work=floater_at(self%core,q(i))
                v(i,:)=real(work,dp)
            end block
        end do
    end function real_floater_values

    function complex_floater_vector_constructor(x,y,d,status) result(p)
        real(dp), intent(in) :: x(:) !! real interpolation nodes
        complex(dp), intent(in) :: y(:) !! finite nodal values, derivatives for repeated Krogh nodes
        integer, intent(in), optional :: d !! polynomial order (default 3), 0 <= d < number of nodes
        integer, intent(out), optional :: status !! zero for a valid interpolator
        type(complex_floater_hormann_interpolator_type) :: p
        integer :: code
        integer :: order
        order=3
        if (present(d)) order=d
        code=global_interp_invalid_input
        if (size(x)==size(y)) then
            call build_floater(x,reshape(y,[size(y),1]),order,p%core,code)
        end if
        if (present(status)) status=code
    end function complex_floater_vector_constructor

    function complex_floater_matrix_constructor(x,y,d,status) result(p)
        real(dp), intent(in) :: x(:) !! real interpolation nodes
        complex(dp), intent(in) :: y(:,:) !! finite nodal values, derivatives for repeated Krogh nodes
        integer, intent(in), optional :: d !! polynomial order (default 3), 0 <= d < number of nodes
        integer, intent(out), optional :: status !! zero for a valid interpolator
        type(complex_floater_hormann_interpolator_type) :: p
        integer :: code
        integer :: order
        order=3
        if (present(d)) order=d
        code=global_interp_invalid_input
        if (size(x)==size(y,1)) then
            call build_floater(x,y,order,p%core,code)
        end if
        if (present(status)) status=code
    end function complex_floater_matrix_constructor

    pure logical function complex_floater_valid(self)
        class(complex_floater_hormann_interpolator_type), intent(in) :: self !! object to validate
        complex_floater_valid=allocated(self%core%nodes) .and. allocated(self%core%ordinates)
    end function complex_floater_valid

    pure integer function complex_floater_channels(self)
        class(complex_floater_hormann_interpolator_type), intent(in) :: self !! interpolator to inspect
        complex_floater_channels=0
        if (self%valid()) complex_floater_channels=size(self%core%ordinates,2)
    end function complex_floater_channels

    function complex_floater_nodes(self) result(v)
        class(complex_floater_hormann_interpolator_type), intent(in) :: self !! fitted interpolation object
        real(dp), allocatable :: v(:)
        if (self%valid()) then
            v=self%core%nodes
        else
            allocate(v(0))
        end if
    end function complex_floater_nodes

    function complex_floater_weights(self) result(v)
        class(complex_floater_hormann_interpolator_type), intent(in) :: self !! rational interpolation object
        real(dp), allocatable :: v(:)
        if (self%valid()) then
            v=self%core%weights
        else
            allocate(v(0))
        end if
    end function complex_floater_weights

    pure integer function complex_floater_degree(self)
        class(complex_floater_hormann_interpolator_type), intent(in) :: self !! rational interpolation object
        complex_floater_degree=self%core%degree
    end function complex_floater_degree

    function complex_floater_value(self,q,channel) result(v)
        class(complex_floater_hormann_interpolator_type), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: q !! query coordinate
        integer, intent(in), optional :: channel !! channel index, default one
        complex(dp) :: v
        integer :: index
        v=0.0_dp
        if (.not.self%valid()) return
        index=1
        if (present(channel)) index=channel
        if (index<1 .or. index>self%channels()) return
        block
            complex(dp) :: work(self%channels())
            work=floater_at(self%core,q)
            v=work(index)
        end block
    end function complex_floater_value

    function complex_floater_values(self,q) result(v)
        class(complex_floater_hormann_interpolator_type), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: q(:) !! query coordinates
        complex(dp), allocatable :: v(:,:)
        integer :: i
        if (.not.self%valid()) then
            allocate(v(0,0))
            return
        end if
        allocate(v(size(q),self%channels()))
        do i=1,size(q)
            block
                complex(dp) :: work(self%channels())
                work=floater_at(self%core,q(i))
                v(i,:)=work
            end block
        end do
    end function complex_floater_values

    function real_krogh_derivatives(self,q,der) result(v)
        class(krogh_interpolator_type), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: q(:) !! query coordinates
        integer, intent(in), optional :: der !! number of derivative orders, default node count
        real(dp), allocatable :: v(:,:,:)
        complex(dp), allocatable :: tmp(:,:)
        integer :: i,count
        if (.not.self%valid()) then
            allocate(v(0,0,0))
            return
        end if
        count=size(self%core%nodes)
        if (present(der)) count=der
        if (count<0) then
            allocate(v(0,0,0))
            return
        end if
        allocate(v(count,size(q),self%channels()))
        do i=1,size(q)
            if (count==0) exit
            tmp=krogh_all_at(self%core,q(i),count)
            v(:,i,:)=real(tmp,dp)
        end do
    end function real_krogh_derivatives

    function complex_krogh_derivatives(self,q,der) result(v)
        class(complex_krogh_interpolator_type), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: q(:) !! query coordinates
        integer, intent(in), optional :: der !! number of derivative orders, default node count
        complex(dp), allocatable :: v(:,:,:)
        complex(dp), allocatable :: tmp(:,:)
        integer :: i,count
        if (.not.self%valid()) then
            allocate(v(0,0,0))
            return
        end if
        count=size(self%core%nodes)
        if (present(der)) count=der
        if (count<0) then
            allocate(v(0,0,0))
            return
        end if
        allocate(v(count,size(q),self%channels()))
        do i=1,size(q)
            if (count==0) exit
            tmp=krogh_all_at(self%core,q(i),count)
            v(:,i,:)=tmp
        end do
    end function complex_krogh_derivatives

    function real_krogh_convenience(x,y,q,nu,status) result(v)
        real(dp), intent(in) :: x(:) !! knots, repeated for derivative constraints
        real(dp), intent(in) :: y(:) !! function and consecutive derivative values
        real(dp), intent(in) :: q(:) !! evaluation coordinates
        integer, intent(in), optional :: nu !! derivative order (default zero)
        integer, intent(out), optional :: status !! zero for valid input
        real(dp), allocatable :: v(:)
        real(dp), allocatable :: work(:,:)
        type(krogh_interpolator_type) :: p
        integer :: code
        p=krogh_interpolator(x,y,status=code)
        if (code==global_interp_success) then
            work=p%values(q,nu=nu)
            v=work(:,1)
        else
            allocate(v(0))
        end if
        if (present(status)) status=code
    end function real_krogh_convenience

    function complex_krogh_convenience(x,y,q,nu,status) result(v)
        real(dp), intent(in) :: x(:) !! knots, repeated for derivative constraints
        complex(dp), intent(in) :: y(:) !! function and consecutive derivative values
        real(dp), intent(in) :: q(:) !! evaluation coordinates
        integer, intent(in), optional :: nu !! derivative order (default zero)
        integer, intent(out), optional :: status !! zero for valid input
        complex(dp), allocatable :: v(:)
        complex(dp), allocatable :: work(:,:)
        type(complex_krogh_interpolator_type) :: p
        integer :: code
        p=krogh_interpolator(x,y,status=code)
        if (code==global_interp_success) then
            work=p%values(q,nu=nu)
            v=work(:,1)
        else
            allocate(v(0))
        end if
        if (present(status)) status=code
    end function complex_krogh_convenience

end module scifort_interpolate_global
