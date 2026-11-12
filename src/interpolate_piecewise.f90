! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Real-valued piecewise polynomials in local power (PPoly) and Bernstein
! (BPoly) bases, independently implemented from the defining polynomials.
! Coefficients follow SciPy's (degree+1, intervals, channels) layout.
module scifort_interpolate_piecewise
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: piecewise_status_success = 0
    integer, parameter, public :: piecewise_status_invalid_input = 1

    ! Optional per-segment real polynomial solver. The module remains fully
    ! usable without third-party root-finding libraries. Optional adapters
    ! (such as polyroots-fortran) implement this exact callback signature.
    abstract interface
        function piecewise_root_solver(coeff,lower,upper) result(roots)
            import :: dp
            real(dp), intent(in) :: coeff(:) !! descending segment coefficients
            real(dp), intent(in) :: lower !! inclusive lower local coordinate
            real(dp), intent(in) :: upper !! inclusive upper local coordinate
            real(dp), allocatable :: roots(:) !! finite real roots in interval
        end function piecewise_root_solver
    end interface
    public :: piecewise_native_real_roots

    type, public :: ppoly_type
        private
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: c(:,:,:)
        logical :: allow_extrapolation = .true.
        logical :: periodic_extension = .false.
    contains
        procedure :: valid => pp_valid
        procedure :: knots => pp_knots
        procedure :: coefficients => pp_coefficients
        procedure :: channels => pp_channels
        procedure, private :: pp_extend_scalar, pp_extend_channels
        generic :: extend => pp_extend_scalar, pp_extend_channels
        procedure, private :: pp_evaluate_scalar, pp_evaluate_array
        generic :: evaluate => pp_evaluate_scalar, pp_evaluate_array
        procedure :: value => pp_value
        procedure :: values => pp_values
        procedure :: derivative => pp_derivative
        procedure :: antiderivative => pp_antiderivative
        procedure :: integrate => pp_integrate
        procedure :: integral => pp_integral
        procedure :: solve => pp_solve
        procedure :: roots => pp_roots
    end type ppoly_type

    type, public :: bpoly_type
        private
        real(dp), allocatable :: x(:)
        real(dp), allocatable :: c(:,:,:)
        logical :: allow_extrapolation = .true.
        logical :: periodic_extension = .false.
    contains
        procedure :: valid => bp_valid
        procedure :: knots => bp_knots
        procedure :: coefficients => bp_coefficients
        procedure :: channels => bp_channels
        procedure, private :: bp_extend_scalar, bp_extend_channels
        generic :: extend => bp_extend_scalar, bp_extend_channels
        procedure, private :: bp_evaluate_scalar, bp_evaluate_array
        generic :: evaluate => bp_evaluate_scalar, bp_evaluate_array
        procedure :: value => bp_value
        procedure :: values => bp_values
        procedure :: derivative => bp_derivative
        procedure :: antiderivative => bp_antiderivative
        procedure :: integrate => bp_integrate
        procedure :: integral => bp_integral
    end type bpoly_type


    ! Complex-valued counterparts keep the same (degree+1, intervals, channels)
    ! coefficient ordering, while the breakpoint coordinates remain real.
    type, public :: complex_ppoly_type
        private
        real(dp), allocatable :: x(:)
        complex(dp), allocatable :: c(:,:,:)
        logical :: allow_extrapolation = .true.
        logical :: periodic_extension = .false.
    contains
        procedure :: valid => cpp_valid
        procedure :: knots => cpp_knots
        procedure :: coefficients => cpp_coefficients
        procedure :: channels => cpp_channels
        procedure, private :: cpp_extend_scalar, cpp_extend_channels
        generic :: extend => cpp_extend_scalar, cpp_extend_channels
        procedure, private :: cpp_evaluate_scalar, cpp_evaluate_array
        generic :: evaluate => cpp_evaluate_scalar, cpp_evaluate_array
        procedure :: value => cpp_value
        procedure :: values => cpp_values
        procedure :: derivative => cpp_derivative
        procedure :: antiderivative => cpp_antiderivative
        procedure :: integrate => cpp_integrate
        procedure :: integral => cpp_integral
    end type complex_ppoly_type

    type, public :: complex_bpoly_type
        private
        real(dp), allocatable :: x(:)
        complex(dp), allocatable :: c(:,:,:)
        logical :: allow_extrapolation = .true.
        logical :: periodic_extension = .false.
    contains
        procedure :: valid => cbp_valid
        procedure :: knots => cbp_knots
        procedure :: coefficients => cbp_coefficients
        procedure :: channels => cbp_channels
        procedure, private :: cbp_extend_scalar, cbp_extend_channels
        generic :: extend => cbp_extend_scalar, cbp_extend_channels
        procedure, private :: cbp_evaluate_scalar, cbp_evaluate_array
        generic :: evaluate => cbp_evaluate_scalar, cbp_evaluate_array
        procedure :: value => cbp_value
        procedure :: values => cbp_values
        procedure :: derivative => cbp_derivative
        procedure :: antiderivative => cbp_antiderivative
        procedure :: integrate => cbp_integrate
        procedure :: integral => cbp_integral
    end type complex_bpoly_type

    public :: complex_ppoly_from_bernstein, complex_bpoly_from_power
    public :: ppoly, bpoly, ppoly_from_bernstein, bpoly_from_power, bpoly_from_derivatives
    interface ppoly
        module procedure ppoly_scalar_coeffs, ppoly_channel_coeffs
        module procedure cppoly_scalar_coeffs, cppoly_channel_coeffs
    end interface
    interface bpoly
        module procedure bpoly_scalar_coeffs, bpoly_channel_coeffs
        module procedure cbpoly_scalar_coeffs, cbpoly_channel_coeffs
    end interface
    interface bpoly_from_derivatives
        module procedure bpoly_derivatives_scalar, bpoly_derivatives_channels
    end interface
contains

    function ppoly_scalar_coeffs(c, x, extrapolate, periodic, status) result(p)
        real(dp), intent(in) :: c(:,:) !! descending power coefficients (degree+1, intervals)
        real(dp), intent(in) :: x(:) !! nondecreasing or nonincreasing breakpoints
        logical, intent(in), optional :: extrapolate !! allow out-of-range evaluation (default true)
        logical, intent(in), optional :: periodic !! periodically extend outside the interval (default false)
        integer, intent(out), optional :: status !! zero for success, nonzero for invalid dimensions/values
        type(ppoly_type) :: p
        real(dp), allocatable :: cm(:,:,:)
        allocate(cm(size(c,1),size(c,2),1))
        cm(:,:,1) = c
        p = ppoly_channel_coeffs(cm, x, extrapolate, periodic, status)
    end function ppoly_scalar_coeffs

    function ppoly_channel_coeffs(c, x, extrapolate, periodic, status) result(p)
        real(dp), intent(in) :: c(:,:,:) !! descending coefficients (degree+1, intervals, channels)
        real(dp), intent(in) :: x(:) !! breakpoints, either monotone direction
        logical, intent(in), optional :: extrapolate !! extrapolate outside original knot span
        logical, intent(in), optional :: periodic !! map out-of-range queries into the original period
        integer, intent(out), optional :: status !! result status (nonzero for invalid input)
        type(ppoly_type) :: p
        integer :: code
        code = piecewise_status_invalid_input
        if (piecewise_valid_input(c,x)) then
            p%x = x
            p%c = c
            if (present(extrapolate)) p%allow_extrapolation = extrapolate
            if (present(periodic)) p%periodic_extension = periodic
            if (p%periodic_extension .and. x(size(x)) == x(1)) then
                deallocate(p%x, p%c)
            else
                code = piecewise_status_success
            end if
        end if
        if (present(status)) status = code
    end function ppoly_channel_coeffs

    function bpoly_scalar_coeffs(c, x, extrapolate, periodic, status) result(p)
        real(dp), intent(in) :: c(:,:) !! Bernstein coefficients (degree+1, intervals)
        real(dp), intent(in) :: x(:) !! monotone breakpoints
        logical, intent(in), optional :: extrapolate !! extrapolation policy
        logical, intent(in), optional :: periodic !! periodic extrapolation policy
        integer, intent(out), optional :: status !! status code
        type(bpoly_type) :: p
        real(dp), allocatable :: cm(:,:,:)
        allocate(cm(size(c,1),size(c,2),1))
        cm(:,:,1) = c
        p = bpoly_channel_coeffs(cm,x,extrapolate,periodic,status)
    end function bpoly_scalar_coeffs

    function bpoly_channel_coeffs(c, x, extrapolate, periodic, status) result(p)
        real(dp), intent(in) :: c(:,:,:) !! Bernstein coefficients (degree+1, intervals, channels)
        real(dp), intent(in) :: x(:) !! monotone breakpoints
        logical, intent(in), optional :: extrapolate !! extrapolation policy
        logical, intent(in), optional :: periodic !! periodic extrapolation policy
        integer, intent(out), optional :: status !! status code
        type(bpoly_type) :: p
        integer :: code
        code = piecewise_status_invalid_input
        if (piecewise_valid_input(c,x)) then
            p%x = x
            p%c = c
            if (present(extrapolate)) p%allow_extrapolation = extrapolate
            if (present(periodic)) p%periodic_extension = periodic
            if (p%periodic_extension .and. x(size(x)) == x(1)) then
                deallocate(p%x, p%c)
            else
                code = piecewise_status_success
            end if
        end if
        if (present(status)) status = code
    end function bpoly_channel_coeffs

    pure logical function piecewise_valid_input(c,x) result(valid)
        real(dp), intent(in) :: c(:,:,:) !! coefficient cube to validate
        real(dp), intent(in) :: x(:) !! candidate breakpoints
        integer :: n
        valid = .false.
        n = size(x)
        if (n < 2 .or. size(c,1) < 1 .or. size(c,2) /= n-1 .or. size(c,3) < 1) return
        if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(c))) return
        if (.not. (all(x(2:n)>=x(1:n-1)) .or. all(x(2:n)<=x(1:n-1)))) return
        if (x(1) == x(n)) return
        valid = .true.
    end function piecewise_valid_input

    pure logical function pp_valid(self)
        class(ppoly_type), intent(in) :: self !! power-basis polynomial
        pp_valid = allocated(self%x) .and. allocated(self%c)
    end function pp_valid

    pure logical function bp_valid(self)
        class(bpoly_type), intent(in) :: self !! Bernstein-basis polynomial
        bp_valid = allocated(self%x) .and. allocated(self%c)
    end function bp_valid

    function pp_knots(self) result(x)
        class(ppoly_type), intent(in) :: self !! power-basis polynomial
        real(dp), allocatable :: x(:)
        if (self%valid()) x = self%x
    end function pp_knots

    function bp_knots(self) result(x)
        class(bpoly_type), intent(in) :: self !! Bernstein-basis polynomial
        real(dp), allocatable :: x(:)
        if (self%valid()) x = self%x
    end function bp_knots

    function pp_coefficients(self) result(c)
        class(ppoly_type), intent(in) :: self !! source polynomial
        real(dp), allocatable :: c(:,:,:)
        if (self%valid()) c = self%c
    end function pp_coefficients

    function bp_coefficients(self) result(c)
        class(bpoly_type), intent(in) :: self !! source polynomial
        real(dp), allocatable :: c(:,:,:)
        if (self%valid()) c = self%c
    end function bp_coefficients

    pure integer function pp_channels(self) result(n)
        class(ppoly_type), intent(in) :: self !! polynomial object
        n = 0
        if (self%valid()) n = size(self%c,3)
    end function pp_channels

    pure integer function bp_channels(self) result(n)
        class(bpoly_type), intent(in) :: self !! polynomial object
        n = 0
        if (self%valid()) n = size(self%c,3)
    end function bp_channels

    pure real(dp) function nan_value() result(v)
        v = ieee_value(0.0_dp, ieee_quiet_nan)
    end function nan_value

    ! Select the interval on the RIGHT at an internal knot, in either
    ! breakpoint direction, as SciPy's PPoly/BPoly evaluator does.
    pure integer function find_interval(x, query) result(j)
        real(dp), intent(in) :: x(:) !! monotone knot vector
        real(dp), intent(in) :: query !! evaluation coordinate
        integer :: lo, hi, mid
        logical :: ascending
        lo = 1
        hi = size(x)
        ascending = x(hi) > x(lo)
        do while (hi-lo>1)
            mid = (lo+hi)/2
            if ((ascending .and. query >= x(mid)) .or. (.not. ascending .and. query <= x(mid))) then
                lo = mid
            else
                hi = mid
            end if
        end do
        j = min(lo,size(x)-1)
    end function find_interval

    pure subroutine adjust_query(x, q, extrapolate, periodic, mapped, allowed)
        real(dp), intent(in) :: x(:) !! breakpoints
        real(dp), intent(in) :: q !! original query
        logical, intent(in) :: extrapolate !! allow unconstrained continuation
        logical, intent(in) :: periodic !! use periodic wrapping
        real(dp), intent(out) :: mapped !! effective query coordinate
        logical, intent(out) :: allowed !! false if query is out of bounds or nonfinite
        real(dp) :: lo, hi, width
        mapped = q
        allowed = .false.
        if (.not. ieee_is_finite(q)) return
        lo = min(x(1),x(size(x)))
        hi = max(x(1),x(size(x)))
        if (periodic) then
            width = x(size(x))-x(1)
            mapped = x(1) + modulo(q-x(1),width)
        else if (.not. extrapolate) then
            if (q < lo .or. q > hi) return
        end if
        allowed = .true.
    end subroutine adjust_query

    pure real(dp) function evaluate_power(c, h, nu) result(v)
        real(dp), intent(in) :: c(:) !! high-to-low power coefficients
        real(dp), intent(in) :: h !! displacement from interval's first endpoint
        integer, intent(in) :: nu !! derivative order, nonnegative
        integer :: k, i, d
        real(dp) :: fac
        k = size(c)-1
        if (nu > k) then
            v = 0.0_dp
            return
        end if
        v = 0.0_dp
        do i = 1,k-nu+1
            fac = 1.0_dp
            do d = 1,nu
                fac = fac*real(k-(i-1)-(d-1),dp)
            end do
            v = v*h + fac*c(i)
        end do
    end function evaluate_power

    pure real(dp) function evaluate_bernstein(c, h, width, nu) result(v)
        real(dp), intent(in) :: c(:) !! Bernstein coefficients of a single segment
        real(dp), intent(in) :: h !! offset from interval start
        real(dp), intent(in) :: width !! interval width (signed)
        integer, intent(in) :: nu !! derivative order
        real(dp) :: scratch(size(c)), tmp(size(c)), t, factor
        integer :: k, j, r, p, d
        k = size(c)-1
        if (nu > k) then
            v = 0.0_dp
            return
        end if
        scratch = c
        if (nu > 0) then
            if (width == 0.0_dp) then
                v = nan_value()
                return
            end if
            factor = 1.0_dp
            do d = 1,nu
                factor = factor*real(k-d+1,dp)/width
                do j = 1,k-d+1
                    tmp(j) = scratch(j+1)-scratch(j)
                end do
                scratch(1:k-d+1) = tmp(1:k-d+1)
            end do
            scratch(1:k-nu+1) = scratch(1:k-nu+1)*factor
        end if
        p = k-nu
        if (width == 0.0_dp) then
            t = 0.0_dp
        else
            t = h/width
        end if
        do r = 1,p
            do j = 1,p-r+1
                scratch(j) = (1.0_dp-t)*scratch(j)+t*scratch(j+1)
            end do
        end do
        v = scratch(1)
    end function evaluate_bernstein

    function pp_evaluate_scalar(self, t, nu, extrapolate) result(y)
        class(ppoly_type), intent(in) :: self !! power-basis polynomial
        real(dp), intent(in) :: t !! evaluation location
        integer, intent(in), optional :: nu !! derivative order (default zero)
        logical, intent(in), optional :: extrapolate !! override nonperiodic extrapolation
        real(dp), allocatable :: y(:)
        real(dp) :: mapped
        integer :: j, r, deriv
        logical :: allowed, ext
        if (.not. self%valid()) then
            allocate(y(0))
            return
        end if
        allocate(y(self%channels()))
        deriv = 0
        if (present(nu)) deriv = nu
        if (deriv<0) then
            y = nan_value()
            return
        end if
        ext = self%allow_extrapolation
        if (present(extrapolate)) ext = extrapolate
        call adjust_query(self%x,t,ext,self%periodic_extension .and. .not. present(extrapolate),mapped,allowed)
        if (.not. allowed) then
            y = nan_value()
            return
        end if
        j = find_interval(self%x,mapped)
        do r = 1,self%channels()
            y(r) = evaluate_power(self%c(:,j,r),mapped-self%x(j),deriv)
        end do
    end function pp_evaluate_scalar

    function bp_evaluate_scalar(self, t, nu, extrapolate) result(y)
        class(bpoly_type), intent(in) :: self !! Bernstein-basis polynomial
        real(dp), intent(in) :: t !! evaluation location
        integer, intent(in), optional :: nu !! derivative order (default zero)
        logical, intent(in), optional :: extrapolate !! override nonperiodic extrapolation
        real(dp), allocatable :: y(:)
        real(dp) :: mapped
        integer :: j, r, deriv
        logical :: allowed, ext
        if (.not. self%valid()) then
            allocate(y(0))
            return
        end if
        allocate(y(self%channels()))
        deriv = 0
        if (present(nu)) deriv = nu
        if (deriv<0) then
            y = nan_value()
            return
        end if
        ext = self%allow_extrapolation
        if (present(extrapolate)) ext = extrapolate
        call adjust_query(self%x,t,ext,self%periodic_extension .and. .not. present(extrapolate),mapped,allowed)
        if (.not. allowed) then
            y = nan_value()
            return
        end if
        j = find_interval(self%x,mapped)
        do r = 1,self%channels()
            y(r) = evaluate_bernstein(self%c(:,j,r),mapped-self%x(j),self%x(j+1)-self%x(j),deriv)
        end do
    end function bp_evaluate_scalar

    function pp_evaluate_array(self,t,nu,extrapolate) result(y)
        class(ppoly_type), intent(in) :: self !! power-basis polynomial
        real(dp), intent(in) :: t(:) !! query coordinates
        integer, intent(in), optional :: nu !! derivative order
        logical, intent(in), optional :: extrapolate !! extrapolation override
        real(dp), allocatable :: y(:,:)
        integer :: i
        allocate(y(size(t),self%channels()))
        do i=1,size(t)
            y(i,:) = self%evaluate(t(i),nu,extrapolate)
        end do
    end function pp_evaluate_array

    function bp_evaluate_array(self,t,nu,extrapolate) result(y)
        class(bpoly_type), intent(in) :: self !! Bernstein-basis polynomial
        real(dp), intent(in) :: t(:) !! query coordinates
        integer, intent(in), optional :: nu !! derivative order
        logical, intent(in), optional :: extrapolate !! extrapolation override
        real(dp), allocatable :: y(:,:)
        integer :: i
        allocate(y(size(t),self%channels()))
        do i=1,size(t)
            y(i,:) = self%evaluate(t(i),nu,extrapolate)
        end do
    end function bp_evaluate_array

    function pp_value(self,t,nu,extrapolate) result(y)
        class(ppoly_type), intent(in) :: self !! one-channel polynomial
        real(dp), intent(in) :: t !! query coordinate
        integer, intent(in), optional :: nu !! derivative order
        logical, intent(in), optional :: extrapolate !! extrapolation override
        real(dp) :: y
        real(dp), allocatable :: a(:)
        a = self%evaluate(t,nu,extrapolate)
        y = nan_value()
        if (size(a)==1) y=a(1)
    end function pp_value

    function bp_value(self,t,nu,extrapolate) result(y)
        class(bpoly_type), intent(in) :: self !! one-channel polynomial
        real(dp), intent(in) :: t !! query coordinate
        integer, intent(in), optional :: nu !! derivative order
        logical, intent(in), optional :: extrapolate !! extrapolation override
        real(dp) :: y
        real(dp), allocatable :: a(:)
        a = self%evaluate(t,nu,extrapolate)
        y = nan_value()
        if (size(a)==1) y=a(1)
    end function bp_value

    function pp_values(self,t,nu,extrapolate) result(y)
        class(ppoly_type), intent(in) :: self !! one-channel polynomial
        real(dp), intent(in) :: t(:) !! query coordinates
        integer, intent(in), optional :: nu !! derivative order
        logical, intent(in), optional :: extrapolate !! extrapolation override
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: a(:,:)
        a = self%evaluate(t,nu,extrapolate)
        allocate(y(size(t)))
        y = nan_value()
        if (size(a,2)==1) y=a(:,1)
    end function pp_values

    function bp_values(self,t,nu,extrapolate) result(y)
        class(bpoly_type), intent(in) :: self !! one-channel polynomial
        real(dp), intent(in) :: t(:) !! query coordinates
        integer, intent(in), optional :: nu !! derivative order
        logical, intent(in), optional :: extrapolate !! extrapolation override
        real(dp), allocatable :: y(:)
        real(dp), allocatable :: a(:,:)
        a = self%evaluate(t,nu,extrapolate)
        allocate(y(size(t)))
        y = nan_value()
        if (size(a,2)==1) y=a(:,1)
    end function bp_values

    recursive function pp_derivative(self,nu) result(out)
        class(ppoly_type), intent(in) :: self !! source polynomial
        integer, intent(in), optional :: nu !! differentiation order (default one; negative integrates)
        type(ppoly_type) :: out
        real(dp), allocatable :: cnew(:,:,:)
        integer :: order, k, i, d
        real(dp) :: scale
        if (.not. self%valid()) return
        order=1
        if (present(nu)) order=nu
        if (order<0) then
            out=self%antiderivative(-order)
            return
        end if
        if (order==0) then
            out=self
            return
        end if
        k=size(self%c,1)-1
        allocate(cnew(max(1,k+1-order),size(self%c,2),size(self%c,3)))
        cnew=0.0_dp
        if (order<=k) then
            do i=1,k+1-order
                scale=1.0_dp
                do d=1,order
                    scale=scale*real(k-i+2-d,dp)
                end do
                cnew(i,:,:)=self%c(i,:,:)*scale
            end do
        end if
        out=ppoly(cnew,self%x,extrapolate=self%allow_extrapolation,periodic=self%periodic_extension)
    end function pp_derivative

    recursive function pp_antiderivative(self,nu) result(out)
        class(ppoly_type), intent(in) :: self !! source polynomial
        integer, intent(in), optional :: nu !! integration order (default one; negative differentiates)
        type(ppoly_type) :: out
        real(dp), allocatable :: cnew(:,:,:)
        real(dp) :: dx, scale, val
        integer :: k, order, i, j, r, q
        if (.not. self%valid()) return
        order=1
        if (present(nu)) order=nu
        if (order<0) then
            out=self%derivative(-order)
            return
        end if
        if (order==0) then
            out=self
            return
        end if
        k=size(self%c,1)-1
        allocate(cnew(k+order+1,size(self%c,2),size(self%c,3)))
        cnew=0.0_dp
        do i=1,k+1
            scale=1.0_dp
            do q=1,order
                scale=scale/real(k-i+1+q,dp)
            end do
            cnew(i,:,:)=self%c(i,:,:)*scale
        end do
        ! Integrate repeatedly to ensure continuity of every intermediate
        ! antiderivative, rather than correcting only the final constant.
        out=ppoly(cnew,self%x,extrapolate=self%allow_extrapolation,periodic=.false.)
        if (order>1) then
            out=self%antiderivative(1)
            do q=2,order
                out=out%antiderivative(1)
            end do
            return
        end if
        do j=2,size(self%x)-1
            dx=self%x(j)-self%x(j-1)
            do r=1,size(self%c,3)
                val=evaluate_power(out%c(:,j-1,r),dx,0)
                out%c(k+2,j,r)=val
            end do
        end do
    end function pp_antiderivative

    recursive function bp_derivative(self,nu) result(out)
        class(bpoly_type), intent(in) :: self !! source Bernstein polynomial
        integer, intent(in), optional :: nu !! derivative order; negative integrates
        type(bpoly_type) :: out
        real(dp), allocatable :: cnew(:,:,:)
        integer :: order,k,j,m
        real(dp) :: width
        if (.not. self%valid()) return
        order=1
        if (present(nu)) order=nu
        if (order<0) then
            out=self%antiderivative(-order)
            return
        end if
        if (order==0) then
            out=self
            return
        end if
        if (order>1) then
            out=self%derivative(1)
            do j=2,order
                out=out%derivative(1)
            end do
            return
        end if
        k=size(self%c,1)-1
        allocate(cnew(max(1,k),size(self%c,2),size(self%c,3)))
        cnew=0.0_dp
        if (k>0) then
            do m=1,size(self%c,2)
                width=self%x(m+1)-self%x(m)
                if (width==0.0_dp) then
                    cnew(:,m,:)=nan_value()
                else
                    cnew(:,m,:)=real(k,dp)/width*(self%c(2:k+1,m,:)-self%c(1:k,m,:))
                end if
            end do
        end if
        ! Internal constructor, not input validation: zero-width intervals
        ! may legitimately have nonfinite differentiated coefficients.
        out%x=self%x
        out%c=cnew
        out%allow_extrapolation=self%allow_extrapolation
        out%periodic_extension=self%periodic_extension
    end function bp_derivative

    recursive function bp_antiderivative(self,nu) result(out)
        class(bpoly_type), intent(in) :: self !! source Bernstein polynomial
        integer, intent(in), optional :: nu !! integration order; negative differentiates
        type(bpoly_type) :: out
        real(dp), allocatable :: cnew(:,:,:)
        real(dp) :: width
        integer :: order,k,i,j
        if (.not. self%valid()) return
        order=1
        if (present(nu)) order=nu
        if (order<0) then
            out=self%derivative(-order)
            return
        end if
        if (order==0) then
            out=self
            return
        end if
        if (order>1) then
            out=self%antiderivative(1)
            do i=2,order
                out=out%antiderivative(1)
            end do
            return
        end if
        k=size(self%c,1)-1
        allocate(cnew(k+2,size(self%c,2),size(self%c,3)))
        cnew=0.0_dp
        do j=1,size(self%c,2)
            if (j>1) cnew(1,j,:)=cnew(k+2,j-1,:)
            width=(self%x(j+1)-self%x(j))/real(k+1,dp)
            do i=1,k+1
                cnew(i+1,j,:)=cnew(i,j,:)+width*self%c(i,j,:)
            end do
        end do
        out=bpoly(cnew,self%x,extrapolate=self%allow_extrapolation,periodic=.false.)
    end function bp_antiderivative

    function pp_integrate(self,a,b,extrapolate) result(y)
        class(ppoly_type), intent(in) :: self !! source polynomial
        real(dp), intent(in) :: a !! lower integration limit
        real(dp), intent(in) :: b !! upper integration limit
        logical, intent(in), optional :: extrapolate !! override nonperiodic extrapolation
        real(dp), allocatable :: y(:)
        type(ppoly_type) :: primitive
        if (.not. self%valid()) then
            allocate(y(0))
            return
        end if
        primitive=self%antiderivative()
        y=integrate_via_primitive(primitive,self%x,a,b,self%allow_extrapolation,&
            self%periodic_extension,extrapolate)
    end function pp_integrate

    function bp_integrate(self,a,b,extrapolate) result(y)
        class(bpoly_type), intent(in) :: self !! source Bernstein polynomial
        real(dp), intent(in) :: a !! lower integration limit
        real(dp), intent(in) :: b !! upper integration limit
        logical, intent(in), optional :: extrapolate !! override nonperiodic extrapolation
        real(dp), allocatable :: y(:)
        type(bpoly_type) :: primitive
        logical :: periodic, ext
        real(dp) :: lo, hi, width, qa, qb
        if (.not. self%valid()) then
            allocate(y(0))
            return
        end if
        primitive=self%antiderivative()
        ext=self%allow_extrapolation
        if (present(extrapolate)) ext=extrapolate
        periodic=self%periodic_extension .and. .not. present(extrapolate)
        lo=minval(self%x)
        hi=maxval(self%x)
        width=hi-lo
        allocate(y(self%channels()))
        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
            y=nan_value()
        else if (periodic) then
            qa=lo+modulo(a-lo,width)
            qb=lo+modulo(b-lo,width)
            y=primitive%evaluate(qb)-primitive%evaluate(qa)+&
                (floor((b-lo)/width)-floor((a-lo)/width))* &
                (primitive%evaluate(hi)-primitive%evaluate(lo))
        else if (.not. ext .and. (a<lo .or. a>hi .or. b<lo .or. b>hi)) then
            y=nan_value()
        else
            y=primitive%evaluate(b)-primitive%evaluate(a)
        end if
    end function bp_integrate

    function integrate_via_primitive(primitive,x,a,b,allow_extrapolation,periodic_extension,override) result(y)
        type(ppoly_type), intent(in) :: primitive !! nonperiodic integral polynomial
        real(dp), intent(in) :: x(:) !! source breakpoints
        real(dp), intent(in) :: a !! first integration bound
        real(dp), intent(in) :: b !! second integration bound
        logical, intent(in) :: allow_extrapolation !! source extrapolation rule
        logical, intent(in) :: periodic_extension !! source periodic rule
        logical, intent(in), optional :: override !! caller extrapolation override
        real(dp), allocatable :: y(:)
        real(dp) :: lo,hi,width,qa,qb
        logical :: ext,periodic
        ext=allow_extrapolation
        if (present(override)) ext=override
        periodic=periodic_extension .and. .not. present(override)
        lo=minval(x)
        hi=maxval(x)
        width=hi-lo
        allocate(y(primitive%channels()))
        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
            y=nan_value()
        else if (periodic) then
            qa=lo+modulo(a-lo,width)
            qb=lo+modulo(b-lo,width)
            y=primitive%evaluate(qb)-primitive%evaluate(qa)+&
                (floor((b-lo)/width)-floor((a-lo)/width))* &
                (primitive%evaluate(hi)-primitive%evaluate(lo))
        else if (.not. ext .and. (a<lo .or. a>hi .or. b<lo .or. b>hi)) then
            y=nan_value()
        else
            y=primitive%evaluate(b)-primitive%evaluate(a)
        end if
    end function integrate_via_primitive

    function pp_integral(self,a,b,extrapolate) result(y)
        class(ppoly_type), intent(in) :: self !! single-channel source polynomial
        real(dp), intent(in) :: a !! first integration bound
        real(dp), intent(in) :: b !! second integration bound
        logical, intent(in), optional :: extrapolate !! override extrapolation
        real(dp) :: y
        real(dp), allocatable :: vals(:)
        vals=self%integrate(a,b,extrapolate)
        y=nan_value()
        if (size(vals)==1) y=vals(1)
    end function pp_integral

    function bp_integral(self,a,b,extrapolate) result(y)
        class(bpoly_type), intent(in) :: self !! single-channel source polynomial
        real(dp), intent(in) :: a !! first integration bound
        real(dp), intent(in) :: b !! second integration bound
        logical, intent(in), optional :: extrapolate !! override extrapolation
        real(dp) :: y
        real(dp), allocatable :: vals(:)
        vals=self%integrate(a,b,extrapolate)
        y=nan_value()
        if (size(vals)==1) y=vals(1)
    end function bp_integral

    pure real(dp) function binomial_coefficient(n,k) result(v)
        integer, intent(in) :: n !! upper binomial argument
        integer, intent(in) :: k !! lower binomial argument
        integer :: j
        v=1.0_dp
        if (k<0 .or. k>n) then
            v=0.0_dp
            return
        end if
        do j=1,min(k,n-k)
            v=v*real(n-j+1,dp)/real(j,dp)
        end do
    end function binomial_coefficient

    function bpoly_from_power(pp,extrapolate,periodic) result(bp)
        type(ppoly_type), intent(in) :: pp !! source PPoly object
        logical, intent(in), optional :: extrapolate !! override extrapolation policy
        logical, intent(in), optional :: periodic !! override periodic extension
        type(bpoly_type) :: bp
        real(dp), allocatable :: c(:,:,:)
        real(dp) :: width, v
        integer :: k,i,j,r,s
        logical :: ext,per
        if (.not. pp%valid()) return
        k=size(pp%c,1)-1
        allocate(c(k+1,size(pp%c,2),pp%channels()))
        do j=1,size(pp%c,2)
            width=pp%x(j+1)-pp%x(j)
            do r=1,pp%channels()
                do i=0,k
                    v=0.0_dp
                    do s=0,i
                        v=v+pp%c(k+1-s,j,r)*width**s* &
                            binomial_coefficient(i,s)/binomial_coefficient(k,s)
                    end do
                    c(i+1,j,r)=v
                end do
            end do
        end do
        ext=pp%allow_extrapolation
        per=pp%periodic_extension
        if (present(extrapolate)) ext=extrapolate
        if (present(periodic)) per=periodic
        bp=bpoly(c,pp%x,extrapolate=ext,periodic=per)
    end function bpoly_from_power

    function ppoly_from_bernstein(bp,extrapolate,periodic) result(pp)
        type(bpoly_type), intent(in) :: bp !! source BPoly object
        logical, intent(in), optional :: extrapolate !! override extrapolation policy
        logical, intent(in), optional :: periodic !! override periodic extension
        type(ppoly_type) :: pp
        real(dp), allocatable :: c(:,:,:)
        real(dp) :: width, v
        integer :: k,i,j,r,s
        logical :: ext,per
        if (.not. bp%valid()) return
        k=size(bp%c,1)-1
        allocate(c(k+1,size(bp%c,2),bp%channels()))
        do j=1,size(bp%c,2)
            width=bp%x(j+1)-bp%x(j)
            do r=1,bp%channels()
                do i=0,k
                    v=0.0_dp
                    do s=0,i
                        v=v+(-1.0_dp)**(i-s)*binomial_coefficient(i,s)*bp%c(s+1,j,r)
                    end do
                    if (width==0.0_dp .and. i>0) then
                        c(k+1-i,j,r)=nan_value()
                    else
                        c(k+1-i,j,r)=v*binomial_coefficient(k,i)/width**i
                    end if
                end do
            end do
        end do
        ext=bp%allow_extrapolation
        per=bp%periodic_extension
        if (present(extrapolate)) ext=extrapolate
        if (present(periodic)) per=periodic
        ! Construct internally because conversions on repeated knots may have
        ! nonfinite, mathematically undefined power coefficients.
        pp%x=bp%x
        pp%c=c
        pp%allow_extrapolation=ext
        pp%periodic_extension=per
    end function ppoly_from_bernstein

    ! Hermite interpolation from arbitrary numbers of endpoint derivatives.
    ! Derivatives at a given node must be supplied from order zero upwards.
    ! Each interval uses its two endpoint constraints; lower-degree intervals
    ! are elevated exactly to the maximum degree using the Bernstein identity.
    function bpoly_derivatives_scalar(x, y, counts, extrapolate, status, orders, order) result(bp)
        real(dp), intent(in) :: x(:) !! strictly increasing interpolation knots
        real(dp), intent(in) :: y(:,:) !! derivatives (knot, derivative order+1)
        integer, intent(in) :: counts(:) !! available derivative orders at each knot, >=1
        logical, intent(in), optional :: extrapolate !! allow extrapolation beyond first/last knot
        integer, intent(out), optional :: status !! zero on success, invalid-input otherwise
        integer, intent(in), optional :: orders(:) !! per-interval polynomial degrees (>=1), as in SciPy
        integer, intent(in), optional :: order !! one polynomial degree for all intervals (Fortran scalar convenience)
        type(bpoly_type) :: bp
        real(dp), allocatable :: yc(:,:,:)
        allocate(yc(size(y,1),size(y,2),1))
        yc(:,:,1)=y
        bp=bpoly_derivatives_channels(x,yc,counts,extrapolate,status,orders,order)
    end function bpoly_derivatives_scalar

    function bpoly_derivatives_channels(x,y,counts,extrapolate,status,orders,order) result(bp)
        real(dp), intent(in) :: x(:) !! strictly increasing knots, at least two
        real(dp), intent(in) :: y(:,:,:) !! derivatives (knot, derivative order+1, channel)
        integer, intent(in) :: counts(:) !! number of derivatives provided at each knot, >=1
        logical, intent(in), optional :: extrapolate !! allow extrapolation (default true)
        integer, intent(out), optional :: status !! zero on success, invalid-input otherwise
        integer, intent(in), optional :: orders(:) !! per-interval polynomial degrees (>=1), as in SciPy
        integer, intent(in), optional :: order !! one polynomial degree for all intervals (Fortran scalar convenience)
        type(bpoly_type) :: bp
        real(dp), allocatable :: c(:,:,:), local(:,:)
        real(dp) :: width, multiplier, ratio
        integer :: n, nchan, highest, degree, i, r, j, ch, current, left_count, right_count, requested
        integer :: code
        code=piecewise_status_invalid_input
        n=size(x)
        if (n<2 .or. size(y,1)/=n .or. size(counts)/=n .or. size(y,3)<1) goto 900
        if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y))) goto 900
        if (.not. all(x(2:n)>x(1:n-1))) goto 900
        if (any(counts<1) .or. any(counts>size(y,2))) goto 900
        if (present(orders) .and. present(order)) goto 900
        if (present(orders)) then
            if (size(orders)/=n-1) goto 900
            if (any(orders<1)) goto 900
        end if
        if (present(order)) then
            if (order<1) goto 900
        end if
        highest=maxval(counts(1:n-1)+counts(2:n)-1)
        if (present(orders)) highest=max(highest,maxval(orders))
        if (present(order)) highest=max(highest,order)
        nchan=size(y,3)
        allocate(c(highest+1,n-1,nchan),local(highest+1,nchan))
        c=0.0_dp
        do i=1,n-1
            width=x(i+1)-x(i)
            left_count=counts(i)
            right_count=counts(i+1)
            requested=-1
            if (present(orders)) requested=orders(i)
            if (present(order)) requested=order
            if (requested>=1) then
                if (requested+1>left_count+right_count) goto 900
                ! Retain as many constraints as possible from both ends,
                ! preferring the right endpoint when the total is odd.
                left_count=min((requested+1)/2,left_count)
                right_count=min(requested+1-left_count,right_count)
                left_count=min(requested+1-right_count,counts(i))
                if (left_count+right_count/=requested+1) goto 900
            end if
            degree=left_count+right_count-1
            local=0.0_dp
            do ch=1,nchan
                do r=0,left_count-1
                    do j=0,r
                        multiplier=1.0_dp
                        do current=0,j-1
                            multiplier=multiplier*width/real(degree-current,dp)
                        end do
                        local(r+1,ch)=local(r+1,ch)+&
                            binomial_coefficient(r,j)*multiplier*y(i,j+1,ch)
                    end do
                end do
                do r=0,right_count-1
                    do j=0,r
                        multiplier=1.0_dp
                        do current=0,j-1
                            multiplier=multiplier*(-width)/real(degree-current,dp)
                        end do
                        local(degree-r+1,ch)=local(degree-r+1,ch)+&
                            binomial_coefficient(r,j)*multiplier*y(i+1,j+1,ch)
                    end do
                end do
            end do
            do current=degree,highest-1
                local(current+2,:)=local(current+1,:)
                do j=current+1,2,-1
                    ratio=real(j-1,dp)/real(current+1,dp)
                    local(j,:)=ratio*local(j-1,:)+(1.0_dp-ratio)*local(j,:)
                end do
            end do
            c(:,i,:)=local
        end do
        bp=bpoly(c,x,extrapolate=extrapolate)
        if (bp%valid()) code=piecewise_status_success
900     if (present(status)) status=code
    end function bpoly_derivatives_channels

    ! Isolate every distinct real root using the real zeros of the derivative
    ! to divide each interval into monotone subintervals. This also detects
    ! tangential roots, which a simple sign-change search misses. It never
    ! requires a complex companion matrix or a system LAPACK installation.
    recursive function real_polynomial_roots(coeff, lower, upper) result(roots)
        real(dp), intent(in) :: coeff(:) !! descending real power coefficients
        real(dp), intent(in) :: lower !! smallest local coordinate, finite
        real(dp), intent(in) :: upper !! largest local coordinate, finite
        real(dp), allocatable :: roots(:)
        real(dp), allocatable :: critical(:), derivative_coeff(:), nodes(:)
        real(dp) :: fleft, fright, middle, lo, hi, fm, root
        integer :: first, degree, j, i, iteration

        allocate(roots(0))
        first = 1
        do while (first < size(coeff))
            if (coeff(first) /= 0.0_dp) exit
            first = first + 1
        end do
        degree = size(coeff) - first
        if (degree == 0 .or. lower > upper) return
        if (degree == 1) then
            root = -coeff(first+1)/coeff(first)
            if (ieee_is_finite(root)) then
                if (root >= lower - 32.0_dp*epsilon(root)*max(1.0_dp,abs(lower)) .and. &
                    root <= upper + 32.0_dp*epsilon(root)*max(1.0_dp,abs(upper))) then
                    roots = [max(lower,min(upper,root))]
                end if
            end if
            return
        end if

        allocate(derivative_coeff(degree))
        do i = 1, degree
            derivative_coeff(i) = coeff(first+i-1)*real(degree-i+1,dp)
        end do
        critical = real_polynomial_roots(derivative_coeff,lower,upper)
        nodes = [lower,critical,upper]
        do j = 1, size(nodes)
            if (polynomial_is_zero(coeff(first:),nodes(j))) &
                call add_distinct_root(roots,nodes(j))
        end do
        do j = 1, size(nodes)-1
            lo = nodes(j)
            hi = nodes(j+1)
            if (hi <= lo) cycle
            fleft = evaluate_power(coeff(first:),lo,0)
            fright = evaluate_power(coeff(first:),hi,0)
            if (fleft == 0.0_dp .or. fright == 0.0_dp) cycle
            if (sign(1.0_dp,fleft) == sign(1.0_dp,fright)) cycle
            ! Bisection on a monotone segment; the derivative recursively
            ! locates extrema even when a root has even multiplicity.
            do iteration = 1, 110
                middle = lo + 0.5_dp*(hi-lo)
                if (middle == lo .or. middle == hi) exit
                fm = evaluate_power(coeff(first:),middle,0)
                if (fm == 0.0_dp) then
                    lo = middle
                    hi = middle
                    exit
                else if (sign(1.0_dp,fm) == sign(1.0_dp,fleft)) then
                    lo = middle
                    fleft = fm
                else
                    hi = middle
                end if
            end do
            call add_distinct_root(roots,lo+0.5_dp*(hi-lo))
        end do
        call sort_real_roots(roots)
    end function real_polynomial_roots

    function piecewise_native_real_roots(coeff,lower,upper) result(roots)
        real(dp), intent(in) :: coeff(:) !! descending polynomial coefficients
        real(dp), intent(in) :: lower !! finite inclusive minimum local coordinate
        real(dp), intent(in) :: upper !! finite inclusive maximum local coordinate
        real(dp), allocatable :: roots(:)
        roots = real_polynomial_roots(coeff,lower,upper)
    end function piecewise_native_real_roots

    pure logical function polynomial_is_zero(coeff, x) result(is_zero)
        real(dp), intent(in) :: coeff(:) !! descending local coefficients
        real(dp), intent(in) :: x !! local coordinate
        real(dp) :: bound, residual
        integer :: j
        residual = evaluate_power(coeff,x,0)
        bound = abs(coeff(1))
        do j = 2,size(coeff)
            bound = bound*abs(x)+abs(coeff(j))
        end do
        is_zero = abs(residual) <= 16.0_dp*epsilon(1.0_dp)*bound
    end function polynomial_is_zero

    subroutine add_distinct_root(roots, root)
        real(dp), allocatable, intent(inout) :: roots(:) !! unordered accumulated distinct roots
        real(dp), intent(in) :: root !! candidate real root
        integer :: j
        if (.not. ieee_is_finite(root)) return
        do j = 1,size(roots)
            if (abs(roots(j)-root) <= 64.0_dp*epsilon(root)* &
                max(1.0_dp,abs(root),abs(roots(j)))) return
        end do
        roots = [roots,root]
    end subroutine add_distinct_root

    subroutine sort_real_roots(roots)
        real(dp), intent(inout) :: roots(:) !! real roots to sort in increasing order
        real(dp) :: value
        integer :: i, j
        do i = 2,size(roots)
            value = roots(i)
            j = i-1
            do while (j >= 1)
                if (roots(j) <= value) exit
                roots(j+1) = roots(j)
                j = j-1
            end do
            roots(j+1) = value
        end do
    end subroutine sort_real_roots

    pure real(dp) function polynomial_root_bound(coeff, width) result(bound)
        real(dp), intent(in) :: coeff(:) !! descending local coefficients, nonzero polynomial
        real(dp), intent(in) :: width !! magnitude of the original segment width
        real(dp) :: log_ratio
        integer :: i, first
        bound = max(1.0_dp,2.0_dp*abs(width))
        first = 1
        do while (first < size(coeff))
            if (coeff(first) /= 0.0_dp) exit
            first = first+1
        end do
        if (coeff(first) == 0.0_dp) return
        do i = first+1,size(coeff)
            if (coeff(i) == 0.0_dp) cycle
            log_ratio = (log(abs(coeff(i)))-log(abs(coeff(first))))/real(i-first,dp)
            bound = max(bound,2.0_dp*exp(min(log_ratio,log(huge(bound))-log(4.0_dp))))
        end do
    end function polynomial_root_bound

    ! SciPy-compatible real equation solving: polynomial roots, tangencies,
    ! discontinuity sign changes, and the [start, NaN] zero-section marker.
    ! For multiple channels, select a channel explicitly (default: first).
    function pp_solve(self,y,discontinuity,extrapolate,channel,status,root_solver) result(roots)
        class(ppoly_type), intent(in) :: self !! piecewise real polynomial
        real(dp), intent(in), optional :: y !! solve p(x)=y; default zero
        logical, intent(in), optional :: discontinuity !! include jumps crossing y; default true
        logical, intent(in), optional :: extrapolate !! extend first/last pieces; default source setting
        integer, intent(in), optional :: channel !! selected output channel (1-based, default 1)
        integer, intent(out), optional :: status !! 0 success, 1 invalid source/inputs
        procedure(piecewise_root_solver), optional :: root_solver !! optional real-root backend
        real(dp), allocatable :: roots(:)
        real(dp), allocatable :: coeff(:), local(:), vals(:), buffer(:)
        real(dp) :: target, dx, lo, hi, bound, left_value, right_value, root
        integer :: m, k, ch, j, i, n, r, code
        logical :: ext, disc, ascending

        allocate(roots(0))
        code = piecewise_status_invalid_input
        if (.not. self%valid()) goto 900
        ch = 1
        if (present(channel)) ch = channel
        if (ch < 1 .or. ch > self%channels()) goto 900
        target = 0.0_dp
        if (present(y)) target = y
        if (.not. ieee_is_finite(target)) goto 900
        ext = self%allow_extrapolation
        if (self%periodic_extension) ext = .true.
        if (present(extrapolate)) ext = extrapolate
        disc = .true.
        if (present(discontinuity)) disc = discontinuity
        m = size(self%c,2)
        k = size(self%c,1)-1
        ascending = self%x(m+1) > self%x(1)
        allocate(buffer(m*(k+5)+2),coeff(k+1))
        n = 0
        do j = 1,m
            coeff = self%c(:,j,ch)
            coeff(k+1) = coeff(k+1)-target
            dx = self%x(j+1)-self%x(j)

            if (disc .and. j > 1) then
                left_value = evaluate_power(self%c(:,j-1,ch),self%x(j)-self%x(j-1),0)-target
                right_value = coeff(k+1)
                if (left_value /= 0.0_dp .and. right_value /= 0.0_dp) then
                    if (sign(1.0_dp,left_value) /= sign(1.0_dp,right_value)) &
                        call append_root(buffer,n,self%x(j),.true.)
                end if
            end if

            if (all(coeff == 0.0_dp)) then
                call append_root(buffer,n,self%x(j),.false.)
                ! A degenerate zero-length segment has just a point root;
                ! SciPy omits the NaN marker reserved for open zero spans.
                if (dx /= 0.0_dp) call append_root(buffer,n,nan_value(),.false.)
                cycle
            end if

            lo = min(0.0_dp,dx)
            hi = max(0.0_dp,dx)
            if (ext .and. (j == 1 .or. j == m)) then
                bound = polynomial_root_bound(coeff,dx)
                if (j == 1) then
                    if (ascending) then
                        lo = -bound
                    else
                        hi = bound
                    end if
                end if
                if (j == m) then
                    if (ascending) then
                        hi = bound
                    else
                        lo = -bound
                    end if
                end if
            end if
            if (present(root_solver)) then
                local = root_solver(coeff,lo,hi)
                if (any(.not.ieee_is_finite(local))) local = real_polynomial_roots(coeff,lo,hi)
                ! A third-party solver can give roots in arbitrary order and
                ! slightly outside the requested closed interval. Validate
                ! and order its output without changing the piecewise policy.
                do r = 1,size(local)
                    if (local(r) < lo-64.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(lo)) .or. &
                        local(r) > hi+64.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(hi))) then
                        local(r) = nan_value()
                    else
                        local(r) = min(hi,max(lo,local(r)))
                    end if
                end do
                local = pack(local,ieee_is_finite(local))
                call sort_real_roots(local)
            else
                local = real_polynomial_roots(coeff,lo,hi)
            end if
            ! SciPy orders roots within each polynomial section by increasing
            ! coordinate, even if the original breakpoint array descends.
            do i = 1,size(local)
                root = self%x(j)+local(i)
                call append_root(buffer,n,root,.true.)
            end do
        end do
        roots = buffer(:n)
        code = piecewise_status_success
900     if (present(status)) status = code
    end function pp_solve

    function pp_roots(self,discontinuity,extrapolate,channel,status,root_solver) result(roots)
        class(ppoly_type), intent(in) :: self !! real piecewise polynomial
        logical, intent(in), optional :: discontinuity !! include jump sign changes, default true
        logical, intent(in), optional :: extrapolate !! endpoint extrapolation override
        integer, intent(in), optional :: channel !! selected output channel (default first)
        integer, intent(out), optional :: status !! zero for success, nonzero for invalid inputs
        procedure(piecewise_root_solver), optional :: root_solver !! optional per-segment root backend
        real(dp), allocatable :: roots(:)
        roots = self%solve(y=0.0_dp,discontinuity=discontinuity, &
            extrapolate=extrapolate,channel=channel,status=status,root_solver=root_solver)
    end function pp_roots

    subroutine append_root(buffer,n,root,deduplicate)
        real(dp), intent(inout) :: buffer(:) !! allocated output workspace
        integer, intent(inout) :: n !! current root count
        real(dp), intent(in) :: root !! root or NaN zero-interval marker
        logical, intent(in) :: deduplicate !! suppress adjacent near-equal finite roots
        if (deduplicate .and. n>0 .and. ieee_is_finite(root)) then
            if (ieee_is_finite(buffer(n))) then
                if (abs(buffer(n)-root) <= 64.0_dp*epsilon(root)* &
                    max(1.0_dp,abs(root),abs(buffer(n)))) return
            end if
        end if
        if (n >= size(buffer)) error stop 'ppoly roots: workspace exhausted'
        n = n+1
        buffer(n) = root
    end subroutine append_root


    pure logical function complex_valid_input(c,x) result(ok)
        complex(dp), intent(in) :: c(:,:,:) !! candidate complex coefficient array
        real(dp), intent(in) :: x(:) !! real monotone knot array
        integer :: n
        n = size(x)
        ok = .false.
        if (n<2 .or. size(c,1)<1 .or. size(c,2)/=n-1 .or. size(c,3)<1) return
        if (.not. all(ieee_is_finite(x))) return
        if (.not. all(ieee_is_finite(real(c,dp))) .or. .not. all(ieee_is_finite(aimag(c)))) return
        if (.not. (all(x(2:n)>=x(1:n-1)) .or. all(x(2:n)<=x(1:n-1)))) return
        if (x(1)==x(n)) return
        ok = .true.
    end function complex_valid_input

    function cppoly_scalar_coeffs(c,x,extrapolate,periodic,status) result(p)
        complex(dp), intent(in) :: c(:,:) !! power coefficients (degree+1,intervals)
        real(dp), intent(in) :: x(:) !! real breakpoints, either monotone direction
        logical, intent(in), optional :: extrapolate !! allow extrapolation (default true)
        logical, intent(in), optional :: periodic !! wrap out-of-bounds coordinates (default false)
        integer, intent(out), optional :: status !! 0 success; 1 invalid input
        type(complex_ppoly_type) :: p
        complex(dp), allocatable :: cm(:,:,:)
        allocate(cm(size(c,1),size(c,2),1))
        cm(:,:,1)=c
        p=cppoly_channel_coeffs(cm,x,extrapolate,periodic,status)
    end function cppoly_scalar_coeffs

    function cppoly_channel_coeffs(c,x,extrapolate,periodic,status) result(p)
        complex(dp), intent(in) :: c(:,:,:) !! power coefficients (degree+1,intervals,channels)
        real(dp), intent(in) :: x(:) !! monotone real breakpoints
        logical, intent(in), optional :: extrapolate !! allow extrapolation (default true)
        logical, intent(in), optional :: periodic !! periodic extension (default false)
        integer, intent(out), optional :: status !! 0 success; 1 invalid input
        type(complex_ppoly_type) :: p
        integer :: code
        code=piecewise_status_invalid_input
        if (complex_valid_input(c,x)) then
            p%x=x
            p%c=c
            if (present(extrapolate)) p%allow_extrapolation=extrapolate
            if (present(periodic)) p%periodic_extension=periodic
            code=piecewise_status_success
        end if
        if (present(status)) status=code
    end function cppoly_channel_coeffs

    function cbpoly_scalar_coeffs(c,x,extrapolate,periodic,status) result(p)
        complex(dp), intent(in) :: c(:,:) !! Bernstein coefficients (degree+1,intervals)
        real(dp), intent(in) :: x(:) !! monotone real breakpoints
        logical, intent(in), optional :: extrapolate !! allow extrapolation (default true)
        logical, intent(in), optional :: periodic !! periodic extension (default false)
        integer, intent(out), optional :: status !! 0 success; 1 invalid input
        type(complex_bpoly_type) :: p
        complex(dp), allocatable :: cm(:,:,:)
        allocate(cm(size(c,1),size(c,2),1))
        cm(:,:,1)=c
        p=cbpoly_channel_coeffs(cm,x,extrapolate,periodic,status)
    end function cbpoly_scalar_coeffs

    function cbpoly_channel_coeffs(c,x,extrapolate,periodic,status) result(p)
        complex(dp), intent(in) :: c(:,:,:) !! Bernstein coefficients (degree+1,intervals,channels)
        real(dp), intent(in) :: x(:) !! monotone real breakpoints
        logical, intent(in), optional :: extrapolate !! allow extrapolation (default true)
        logical, intent(in), optional :: periodic !! periodic extension (default false)
        integer, intent(out), optional :: status !! 0 success; 1 invalid input
        type(complex_bpoly_type) :: p
        integer :: code
        code=piecewise_status_invalid_input
        if (complex_valid_input(c,x)) then
            p%x=x
            p%c=c
            if (present(extrapolate)) p%allow_extrapolation=extrapolate
            if (present(periodic)) p%periodic_extension=periodic
            code=piecewise_status_success
        end if
        if (present(status)) status=code
    end function cbpoly_channel_coeffs

    pure logical function cpp_valid(self)
        class(complex_ppoly_type), intent(in) :: self !! complex PPoly object
        cpp_valid=allocated(self%x) .and. allocated(self%c)
    end function cpp_valid

    pure logical function cbp_valid(self)
        class(complex_bpoly_type), intent(in) :: self !! complex BPoly object
        cbp_valid=allocated(self%x) .and. allocated(self%c)
    end function cbp_valid

    function cpp_knots(self) result(x)
        class(complex_ppoly_type), intent(in) :: self !! complex PPoly object
        real(dp), allocatable :: x(:)
        if (self%valid()) x=self%x
    end function cpp_knots

    function cbp_knots(self) result(x)
        class(complex_bpoly_type), intent(in) :: self !! complex BPoly object
        real(dp), allocatable :: x(:)
        if (self%valid()) x=self%x
    end function cbp_knots

    function cpp_coefficients(self) result(c)
        class(complex_ppoly_type), intent(in) :: self !! source polynomial
        complex(dp), allocatable :: c(:,:,:)
        if (self%valid()) c=self%c
    end function cpp_coefficients

    function cbp_coefficients(self) result(c)
        class(complex_bpoly_type), intent(in) :: self !! source polynomial
        complex(dp), allocatable :: c(:,:,:)
        if (self%valid()) c=self%c
    end function cbp_coefficients

    pure integer function cpp_channels(self) result(n)
        class(complex_ppoly_type), intent(in) :: self !! source polynomial
        n=0
        if (self%valid()) n=size(self%c,3)
    end function cpp_channels

    pure integer function cbp_channels(self) result(n)
        class(complex_bpoly_type), intent(in) :: self !! source polynomial
        n=0
        if (self%valid()) n=size(self%c,3)
    end function cbp_channels

    pure complex(dp) function complex_nan_value() result(v)
        v=cmplx(nan_value(),nan_value(),dp)
    end function complex_nan_value

    pure complex(dp) function evaluate_complex_power(c,h,nu) result(v)
        complex(dp), intent(in) :: c(:) !! descending power coefficients
        real(dp), intent(in) :: h !! local real displacement
        integer, intent(in) :: nu !! derivative order, >=0
        integer :: k, i, j
        real(dp) :: factor
        k=size(c)-1
        v=cmplx(0.0_dp,0.0_dp,dp)
        if (nu>k) return
        do i=1,k-nu+1
            factor=1.0_dp
            do j=1,nu
                factor=factor*real(k-i+2-j,dp)
            end do
            v=v*h+factor*c(i)
        end do
    end function evaluate_complex_power

    pure complex(dp) function evaluate_complex_bernstein(c,h,width,nu) result(v)
        complex(dp), intent(in) :: c(:) !! Bernstein control points
        real(dp), intent(in) :: h !! local displacement from knot
        real(dp), intent(in) :: width !! signed interval width
        integer, intent(in) :: nu !! derivative order, >=0
        complex(dp) :: scratch(size(c)), tmp(size(c))
        real(dp) :: t, factor
        integer :: degree, j, d, r, p
        degree=size(c)-1
        v=cmplx(0.0_dp,0.0_dp,dp)
        if (nu>degree) return
        scratch=c
        if (nu>0) then
            if (width==0.0_dp) then
                v=complex_nan_value()
                return
            end if
            factor=1.0_dp
            do d=1,nu
                factor=factor*real(degree-d+1,dp)/width
                do j=1,degree-d+1
                    tmp(j)=scratch(j+1)-scratch(j)
                end do
                scratch(1:degree-d+1)=tmp(1:degree-d+1)
            end do
            scratch(1:degree-nu+1)=scratch(1:degree-nu+1)*factor
        end if
        p=degree-nu
        t=0.0_dp
        if (width/=0.0_dp) t=h/width
        do r=1,p
            do j=1,p-r+1
                scratch(j)=(1.0_dp-t)*scratch(j)+t*scratch(j+1)
            end do
        end do
        v=scratch(1)
    end function evaluate_complex_bernstein

    function cpp_evaluate_scalar(self,t,nu,extrapolate) result(y)
        class(complex_ppoly_type), intent(in) :: self !! source polynomial
        real(dp), intent(in) :: t !! real evaluation coordinate
        integer, intent(in), optional :: nu !! derivative order, default zero
        logical, intent(in), optional :: extrapolate !! override nonperiodic extrapolation
        complex(dp), allocatable :: y(:)
        real(dp) :: mapped
        integer :: deriv,j,r
        logical :: ext,allowed
        allocate(y(self%channels()))
        if (.not. self%valid()) return
        deriv=0
        if (present(nu)) deriv=nu
        if (deriv<0) then
            y=complex_nan_value()
            return
        end if
        ext=self%allow_extrapolation
        if (present(extrapolate)) ext=extrapolate
        call adjust_query(self%x,t,ext,self%periodic_extension .and. .not. present(extrapolate),mapped,allowed)
        if (.not. allowed) then
            y=complex_nan_value()
            return
        end if
        j=find_interval(self%x,mapped)
        do r=1,self%channels()
            y(r)=evaluate_complex_power(self%c(:,j,r),mapped-self%x(j),deriv)
        end do
    end function cpp_evaluate_scalar

    function cbp_evaluate_scalar(self,t,nu,extrapolate) result(y)
        class(complex_bpoly_type), intent(in) :: self !! source polynomial
        real(dp), intent(in) :: t !! real evaluation coordinate
        integer, intent(in), optional :: nu !! derivative order, default zero
        logical, intent(in), optional :: extrapolate !! override nonperiodic extrapolation
        complex(dp), allocatable :: y(:)
        real(dp) :: mapped
        integer :: deriv,j,r
        logical :: ext,allowed
        allocate(y(self%channels()))
        if (.not. self%valid()) return
        deriv=0
        if (present(nu)) deriv=nu
        if (deriv<0) then
            y=complex_nan_value()
            return
        end if
        ext=self%allow_extrapolation
        if (present(extrapolate)) ext=extrapolate
        call adjust_query(self%x,t,ext,self%periodic_extension .and. .not. present(extrapolate),mapped,allowed)
        if (.not. allowed) then
            y=complex_nan_value()
            return
        end if
        j=find_interval(self%x,mapped)
        do r=1,self%channels()
            y(r)=evaluate_complex_bernstein(self%c(:,j,r),mapped-self%x(j),&
                self%x(j+1)-self%x(j),deriv)
        end do
    end function cbp_evaluate_scalar

    function cpp_evaluate_array(self,t,nu,extrapolate) result(y)
        class(complex_ppoly_type), intent(in) :: self !! source polynomial
        real(dp), intent(in) :: t(:) !! query coordinates
        integer, intent(in), optional :: nu !! derivative order, default zero
        logical, intent(in), optional :: extrapolate !! extrapolation override
        complex(dp), allocatable :: y(:,:)
        integer :: i
        allocate(y(size(t),self%channels()))
        do i=1,size(t)
            y(i,:)=self%evaluate(t(i),nu,extrapolate)
        end do
    end function cpp_evaluate_array

    function cbp_evaluate_array(self,t,nu,extrapolate) result(y)
        class(complex_bpoly_type), intent(in) :: self !! source polynomial
        real(dp), intent(in) :: t(:) !! query coordinates
        integer, intent(in), optional :: nu !! derivative order, default zero
        logical, intent(in), optional :: extrapolate !! extrapolation override
        complex(dp), allocatable :: y(:,:)
        integer :: i
        allocate(y(size(t),self%channels()))
        do i=1,size(t)
            y(i,:)=self%evaluate(t(i),nu,extrapolate)
        end do
    end function cbp_evaluate_array

    function cpp_value(self,t,nu,extrapolate) result(y)
        class(complex_ppoly_type), intent(in) :: self !! single-channel polynomial
        real(dp), intent(in) :: t !! real query coordinate
        integer, intent(in), optional :: nu !! derivative order
        logical, intent(in), optional :: extrapolate !! extrapolation override
        complex(dp) :: y
        complex(dp), allocatable :: values(:)
        y=complex_nan_value()
        values=self%evaluate(t,nu,extrapolate)
        if (size(values)==1) y=values(1)
    end function cpp_value

    function cbp_value(self,t,nu,extrapolate) result(y)
        class(complex_bpoly_type), intent(in) :: self !! single-channel polynomial
        real(dp), intent(in) :: t !! real query coordinate
        integer, intent(in), optional :: nu !! derivative order
        logical, intent(in), optional :: extrapolate !! extrapolation override
        complex(dp) :: y
        complex(dp), allocatable :: values(:)
        y=complex_nan_value()
        values=self%evaluate(t,nu,extrapolate)
        if (size(values)==1) y=values(1)
    end function cbp_value

    function cpp_values(self,t,nu,extrapolate) result(y)
        class(complex_ppoly_type), intent(in) :: self !! single-channel polynomial
        real(dp), intent(in) :: t(:) !! real query coordinates
        integer, intent(in), optional :: nu !! derivative order
        logical, intent(in), optional :: extrapolate !! extrapolation override
        complex(dp), allocatable :: y(:), values(:,:)
        allocate(y(size(t)))
        y=complex_nan_value()
        values=self%evaluate(t,nu,extrapolate)
        if (size(values,2)==1) y=values(:,1)
    end function cpp_values

    function cbp_values(self,t,nu,extrapolate) result(y)
        class(complex_bpoly_type), intent(in) :: self !! single-channel polynomial
        real(dp), intent(in) :: t(:) !! real query coordinates
        integer, intent(in), optional :: nu !! derivative order
        logical, intent(in), optional :: extrapolate !! extrapolation override
        complex(dp), allocatable :: y(:), values(:,:)
        allocate(y(size(t)))
        y=complex_nan_value()
        values=self%evaluate(t,nu,extrapolate)
        if (size(values,2)==1) y=values(:,1)
    end function cbp_values

    recursive function cpp_derivative(self,nu) result(out)
        class(complex_ppoly_type), intent(in) :: self !! complex polynomial
        integer, intent(in), optional :: nu !! order (default one; negative integrates)
        type(complex_ppoly_type) :: out
        complex(dp), allocatable :: cnew(:,:,:)
        real(dp) :: factor
        integer :: order,degree,i,d
        if (.not. self%valid()) return
        order=1
        if (present(nu)) order=nu
        if (order<0) then
            out=self%antiderivative(-order)
            return
        end if
        if (order==0) then
            out=self
            return
        end if
        degree=size(self%c,1)-1
        allocate(cnew(max(1,degree+1-order),size(self%c,2),self%channels()))
        cnew=cmplx(0.0_dp,0.0_dp,dp)
        if (order<=degree) then
            do i=1,degree+1-order
                factor=1.0_dp
                do d=1,order
                    factor=factor*real(degree-i+2-d,dp)
                end do
                cnew(i,:,:)=self%c(i,:,:)*factor
            end do
        end if
        out=ppoly(cnew,self%x,extrapolate=self%allow_extrapolation,periodic=self%periodic_extension)
    end function cpp_derivative

    recursive function cpp_antiderivative(self,nu) result(out)
        class(complex_ppoly_type), intent(in) :: self !! complex polynomial
        integer, intent(in), optional :: nu !! order (default one; negative differentiates)
        type(complex_ppoly_type) :: out
        complex(dp), allocatable :: cnew(:,:,:)
        integer :: order,degree,i,j,r
        real(dp) :: width
        if (.not. self%valid()) return
        order=1
        if (present(nu)) order=nu
        if (order<0) then
            out=self%derivative(-order)
            return
        end if
        if (order==0) then
            out=self
            return
        end if
        if (order>1) then
            out=self%antiderivative(1)
            do i=2,order
                out=out%antiderivative(1)
            end do
            return
        end if
        degree=size(self%c,1)-1
        allocate(cnew(degree+2,size(self%c,2),self%channels()))
        cnew=cmplx(0.0_dp,0.0_dp,dp)
        do i=1,degree+1
            cnew(i,:,:)=self%c(i,:,:)/real(degree-i+2,dp)
        end do
        do j=2,size(self%x)-1
            width=self%x(j)-self%x(j-1)
            do r=1,self%channels()
                cnew(degree+2,j,r)=evaluate_complex_power(cnew(:,j-1,r),width,0)
            end do
        end do
        out=ppoly(cnew,self%x,extrapolate=self%allow_extrapolation,periodic=.false.)
    end function cpp_antiderivative

    recursive function cbp_derivative(self,nu) result(out)
        class(complex_bpoly_type), intent(in) :: self !! complex Bernstein polynomial
        integer, intent(in), optional :: nu !! order (default one; negative integrates)
        type(complex_bpoly_type) :: out
        complex(dp), allocatable :: cnew(:,:,:)
        integer :: order,degree,j,i
        real(dp) :: width
        if (.not. self%valid()) return
        order=1
        if (present(nu)) order=nu
        if (order<0) then
            out=self%antiderivative(-order)
            return
        end if
        if (order==0) then
            out=self
            return
        end if
        if (order>1) then
            out=self%derivative(1)
            do j=2,order
                out=out%derivative(1)
            end do
            return
        end if
        degree=size(self%c,1)-1
        allocate(cnew(max(1,degree),size(self%c,2),self%channels()))
        cnew=cmplx(0.0_dp,0.0_dp,dp)
        if (degree>0) then
            do j=1,size(self%c,2)
                width=self%x(j+1)-self%x(j)
                if (width==0.0_dp) then
                    cnew(:,j,:)=complex_nan_value()
                else
                    cnew(:,j,:)=real(degree,dp)/width*(self%c(2:degree+1,j,:)-self%c(1:degree,j,:))
                end if
            end do
        end if
        ! Degenerate intervals may produce nonfinite derivatives. Preserve
        ! their representations internally, just as the real type does.
        out%x=self%x
        out%c=cnew
        out%allow_extrapolation=self%allow_extrapolation
        out%periodic_extension=self%periodic_extension
    end function cbp_derivative

    recursive function cbp_antiderivative(self,nu) result(out)
        class(complex_bpoly_type), intent(in) :: self !! complex Bernstein polynomial
        integer, intent(in), optional :: nu !! order (default one; negative differentiates)
        type(complex_bpoly_type) :: out
        complex(dp), allocatable :: cnew(:,:,:)
        integer :: order,degree,j,i
        real(dp) :: width
        if (.not. self%valid()) return
        order=1
        if (present(nu)) order=nu
        if (order<0) then
            out=self%derivative(-order)
            return
        end if
        if (order==0) then
            out=self
            return
        end if
        if (order>1) then
            out=self%antiderivative(1)
            do i=2,order
                out=out%antiderivative(1)
            end do
            return
        end if
        degree=size(self%c,1)-1
        allocate(cnew(degree+2,size(self%c,2),self%channels()))
        cnew=cmplx(0.0_dp,0.0_dp,dp)
        do j=1,size(self%c,2)
            if (j>1) cnew(1,j,:)=cnew(degree+2,j-1,:)
            width=(self%x(j+1)-self%x(j))/real(degree+1,dp)
            do i=1,degree+1
                cnew(i+1,j,:)=cnew(i,j,:)+width*self%c(i,j,:)
            end do
        end do
        out=bpoly(cnew,self%x,extrapolate=self%allow_extrapolation,periodic=.false.)
    end function cbp_antiderivative

    function cpp_integrate(self,a,b,extrapolate) result(y)
        class(complex_ppoly_type), intent(in) :: self !! complex polynomial
        real(dp), intent(in) :: a !! first integration bound
        real(dp), intent(in) :: b !! second integration bound
        logical, intent(in), optional :: extrapolate !! override extrapolation
        complex(dp), allocatable :: y(:)
        type(complex_ppoly_type) :: primitive
        logical :: ext,per
        real(dp) :: lo,hi,width,qa,qb
        allocate(y(self%channels()))
        if (.not. self%valid()) return
        ext=self%allow_extrapolation
        if (present(extrapolate)) ext=extrapolate
        per=self%periodic_extension .and. .not. present(extrapolate)
        lo=minval(self%x)
        hi=maxval(self%x)
        width=hi-lo
        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
            y=complex_nan_value()
        else if (per) then
            primitive=self%antiderivative()
            qa=lo+modulo(a-lo,width)
            qb=lo+modulo(b-lo,width)
            y=primitive%evaluate(qb)-primitive%evaluate(qa)+&
                real(floor((b-lo)/width,kind=8)-floor((a-lo)/width,kind=8),dp)*&
                (primitive%evaluate(hi)-primitive%evaluate(lo))
        else if (.not. ext .and. (a<lo .or. a>hi .or. b<lo .or. b>hi)) then
            y=complex_nan_value()
        else
            primitive=self%antiderivative()
            y=primitive%evaluate(b,extrapolate=.true.)-primitive%evaluate(a,extrapolate=.true.)
        end if
    end function cpp_integrate

    function cbp_integrate(self,a,b,extrapolate) result(y)
        class(complex_bpoly_type), intent(in) :: self !! complex polynomial
        real(dp), intent(in) :: a !! first integration bound
        real(dp), intent(in) :: b !! second integration bound
        logical, intent(in), optional :: extrapolate !! override extrapolation
        complex(dp), allocatable :: y(:)
        type(complex_bpoly_type) :: primitive
        logical :: ext,per
        real(dp) :: lo,hi,width,qa,qb
        allocate(y(self%channels()))
        if (.not. self%valid()) return
        ext=self%allow_extrapolation
        if (present(extrapolate)) ext=extrapolate
        per=self%periodic_extension .and. .not. present(extrapolate)
        lo=minval(self%x)
        hi=maxval(self%x)
        width=hi-lo
        if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) then
            y=complex_nan_value()
        else if (per) then
            primitive=self%antiderivative()
            qa=lo+modulo(a-lo,width)
            qb=lo+modulo(b-lo,width)
            y=primitive%evaluate(qb)-primitive%evaluate(qa)+&
                real(floor((b-lo)/width,kind=8)-floor((a-lo)/width,kind=8),dp)*&
                (primitive%evaluate(hi)-primitive%evaluate(lo))
        else if (.not. ext .and. (a<lo .or. a>hi .or. b<lo .or. b>hi)) then
            y=complex_nan_value()
        else
            primitive=self%antiderivative()
            y=primitive%evaluate(b,extrapolate=.true.)-primitive%evaluate(a,extrapolate=.true.)
        end if
    end function cbp_integrate

    function cpp_integral(self,a,b,extrapolate) result(y)
        class(complex_ppoly_type), intent(in) :: self !! single-channel complex polynomial
        real(dp), intent(in) :: a !! first integration bound
        real(dp), intent(in) :: b !! second integration bound
        logical, intent(in), optional :: extrapolate !! extrapolation override
        complex(dp) :: y
        complex(dp), allocatable :: values(:)
        y=complex_nan_value()
        values=self%integrate(a,b,extrapolate)
        if (size(values)==1) y=values(1)
    end function cpp_integral

    function cbp_integral(self,a,b,extrapolate) result(y)
        class(complex_bpoly_type), intent(in) :: self !! single-channel complex polynomial
        real(dp), intent(in) :: a !! first integration bound
        real(dp), intent(in) :: b !! second integration bound
        logical, intent(in), optional :: extrapolate !! extrapolation override
        complex(dp) :: y
        complex(dp), allocatable :: values(:)
        y=complex_nan_value()
        values=self%integrate(a,b,extrapolate)
        if (size(values)==1) y=values(1)
    end function cbp_integral

    function complex_bpoly_from_power(pp,extrapolate,periodic) result(bp)
        type(complex_ppoly_type), intent(in) :: pp !! complex power-basis polynomial
        logical, intent(in), optional :: extrapolate !! override extrapolation policy
        logical, intent(in), optional :: periodic !! override periodic extension
        type(complex_bpoly_type) :: bp
        complex(dp), allocatable :: c(:,:,:)
        complex(dp) :: value
        real(dp) :: width
        integer :: k,i,j,r,s
        logical :: ext,per
        if (.not. pp%valid()) return
        k=size(pp%c,1)-1
        allocate(c(k+1,size(pp%c,2),pp%channels()))
        do j=1,size(pp%c,2)
            width=pp%x(j+1)-pp%x(j)
            do r=1,pp%channels()
                do i=0,k
                    value=cmplx(0.0_dp,0.0_dp,dp)
                    do s=0,i
                        value=value+pp%c(k+1-s,j,r)*width**s*&
                            binomial_coefficient(i,s)/binomial_coefficient(k,s)
                    end do
                    c(i+1,j,r)=value
                end do
            end do
        end do
        ext=pp%allow_extrapolation
        per=pp%periodic_extension
        if (present(extrapolate)) ext=extrapolate
        if (present(periodic)) per=periodic
        bp=bpoly(c,pp%x,extrapolate=ext,periodic=per)
    end function complex_bpoly_from_power

    function complex_ppoly_from_bernstein(bp,extrapolate,periodic) result(pp)
        type(complex_bpoly_type), intent(in) :: bp !! complex Bernstein polynomial
        logical, intent(in), optional :: extrapolate !! override extrapolation policy
        logical, intent(in), optional :: periodic !! override periodic extension
        type(complex_ppoly_type) :: pp
        complex(dp), allocatable :: c(:,:,:)
        complex(dp) :: value
        real(dp) :: width
        integer :: k,i,j,r,s
        logical :: ext,per
        if (.not. bp%valid()) return
        k=size(bp%c,1)-1
        allocate(c(k+1,size(bp%c,2),bp%channels()))
        do j=1,size(bp%c,2)
            width=bp%x(j+1)-bp%x(j)
            do r=1,bp%channels()
                do i=0,k
                    value=cmplx(0.0_dp,0.0_dp,dp)
                    do s=0,i
                        value=value+(-1.0_dp)**(i-s)*binomial_coefficient(i,s)*bp%c(s+1,j,r)
                    end do
                    if (width==0.0_dp .and. i>0) then
                        c(k+1-i,j,r)=complex_nan_value()
                    else
                        c(k+1-i,j,r)=value*binomial_coefficient(k,i)/width**i
                    end if
                end do
            end do
        end do
        ext=bp%allow_extrapolation
        per=bp%periodic_extension
        if (present(extrapolate)) ext=extrapolate
        if (present(periodic)) per=periodic
        ! Repeated knots can make the local power derivatives undefined;
        ! retain the resulting coefficients internally rather than rejecting.
        pp%x=bp%x
        pp%c=c
        pp%allow_extrapolation=ext
        pp%periodic_extension=per
    end function complex_ppoly_from_bernstein

    ! Append or prepend intervals, keeping the receiver unchanged on error.
    ! Incoming x consists of new breakpoints, not a repeated boundary knot.
    subroutine pp_extend_scalar(self, c, x, status)
        class(ppoly_type), intent(inout) :: self !! polynomial to extend in place
        real(dp), intent(in) :: c(:,:) !! coefficients (degree+1, new intervals)
        real(dp), intent(in) :: x(:) !! new breakpoints, same orientation as existing knots
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        real(dp), allocatable :: cube(:,:,:)
        allocate(cube(size(c,1),size(c,2),1))
        cube(:,:,1)=c
        call pp_extend_channels(self,cube,x,status)
    end subroutine pp_extend_scalar

    subroutine pp_extend_channels(self, c, x, status)
        class(ppoly_type), intent(inout) :: self !! polynomial receiver, modified on success only
        real(dp), intent(in) :: c(:,:,:) !! coefficients (degree+1, new intervals, channels)
        real(dp), intent(in) :: x(:) !! new breakpoints, sorted with same orientation
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        real(dp), allocatable :: merged(:,:,:), original(:,:,:), incoming(:,:,:)
        real(dp), allocatable :: extended_x(:)
        integer :: code, old_n, added, new_degree
        logical :: append, ascending
        code=piecewise_status_invalid_input
        if (.not. self%valid()) goto 900
        old_n=size(self%c,2)
        added=size(c,2)
        if (size(x)/=added .or. size(c,1)<1 .or. size(c,3)/=size(self%c,3)) goto 900
        if (added==0) then
            code=piecewise_status_success
            goto 900
        end if
        if (.not. all(ieee_is_finite(x))) goto 900
        if (.not. (all(x(2:added)>=x(1:added-1)) .or. &
            all(x(2:added)<=x(1:added-1)))) goto 900
        ascending=self%x(size(self%x))>=self%x(1)
        if (ascending) then
            if (x(added)<x(1)) goto 900
            if (x(1)>=self%x(old_n+1)) then
                append=.true.
            else if (x(added)<=self%x(1)) then
                append=.false.
            else
                goto 900
            end if
        else
            if (x(added)>x(1)) goto 900
            if (x(1)<=self%x(old_n+1)) then
                append=.true.
            else if (x(added)>=self%x(1)) then
                append=.false.
            else
                goto 900
            end if
        end if
        if (.not. all(ieee_is_finite(c))) goto 900
        new_degree=max(size(self%c,1),size(c,1))
        allocate(merged(new_degree,old_n+added,size(c,3)),extended_x(old_n+added+1))
        merged=0.0_dp
        if (append) then
            merged(new_degree-size(self%c,1)+1:new_degree,1:old_n,:)=self%c
            merged(new_degree-size(c,1)+1:new_degree,old_n+1:,:)=c
            extended_x(:old_n+1)=self%x
            extended_x(old_n+2:)=x
        else
            merged(new_degree-size(c,1)+1:new_degree,1:added,:)=c
            merged(new_degree-size(self%c,1)+1:new_degree,added+1:,:)=self%c
            extended_x(:added)=x
            extended_x(added+1:)=self%x
        end if
        call move_alloc(extended_x,self%x)
        call move_alloc(merged,self%c)
        code=piecewise_status_success
900     if (present(status)) status=code
    end subroutine pp_extend_channels

    ! Append or prepend intervals, keeping the receiver unchanged on error.
    ! Incoming x consists of new breakpoints, not a repeated boundary knot.
    subroutine bp_extend_scalar(self, c, x, status)
        class(bpoly_type), intent(inout) :: self !! polynomial to extend in place
        real(dp), intent(in) :: c(:,:) !! coefficients (degree+1, new intervals)
        real(dp), intent(in) :: x(:) !! new breakpoints, same orientation as existing knots
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        real(dp), allocatable :: cube(:,:,:)
        allocate(cube(size(c,1),size(c,2),1))
        cube(:,:,1)=c
        call bp_extend_channels(self,cube,x,status)
    end subroutine bp_extend_scalar

    subroutine bp_extend_channels(self, c, x, status)
        class(bpoly_type), intent(inout) :: self !! polynomial receiver, modified on success only
        real(dp), intent(in) :: c(:,:,:) !! coefficients (degree+1, new intervals, channels)
        real(dp), intent(in) :: x(:) !! new breakpoints, sorted with same orientation
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        real(dp), allocatable :: merged(:,:,:), original(:,:,:), incoming(:,:,:)
        real(dp), allocatable :: extended_x(:)
        integer :: code, old_n, added, new_degree
        logical :: append, ascending
        code=piecewise_status_invalid_input
        if (.not. self%valid()) goto 900
        old_n=size(self%c,2)
        added=size(c,2)
        if (size(x)/=added .or. size(c,1)<1 .or. size(c,3)/=size(self%c,3)) goto 900
        if (added==0) then
            code=piecewise_status_success
            goto 900
        end if
        if (.not. all(ieee_is_finite(x))) goto 900
        if (.not. (all(x(2:added)>=x(1:added-1)) .or. &
            all(x(2:added)<=x(1:added-1)))) goto 900
        ascending=self%x(size(self%x))>=self%x(1)
        if (ascending) then
            if (x(added)<x(1)) goto 900
            if (x(1)>=self%x(old_n+1)) then
                append=.true.
            else if (x(added)<=self%x(1)) then
                append=.false.
            else
                goto 900
            end if
        else
            if (x(added)>x(1)) goto 900
            if (x(1)<=self%x(old_n+1)) then
                append=.true.
            else if (x(added)>=self%x(1)) then
                append=.false.
            else
                goto 900
            end if
        end if
        if (.not. all(ieee_is_finite(c))) goto 900
        new_degree=max(size(self%c,1),size(c,1))
        allocate(merged(new_degree,old_n+added,size(c,3)),extended_x(old_n+added+1))
        merged=0.0_dp
        original=raise_bernstein_real(self%c,new_degree)
        incoming=raise_bernstein_real(c,new_degree)
        if (append) then
            merged(new_degree-size(original,1)+1:new_degree,1:old_n,:)=original
            merged(new_degree-size(incoming,1)+1:new_degree,old_n+1:,:)=incoming
            extended_x(:old_n+1)=self%x
            extended_x(old_n+2:)=x
        else
            merged(new_degree-size(incoming,1)+1:new_degree,1:added,:)=incoming
            merged(new_degree-size(original,1)+1:new_degree,added+1:,:)=original
            extended_x(:added)=x
            extended_x(added+1:)=self%x
        end if
        call move_alloc(extended_x,self%x)
        call move_alloc(merged,self%c)
        code=piecewise_status_success
900     if (present(status)) status=code
    end subroutine bp_extend_channels

    ! Append or prepend intervals, keeping the receiver unchanged on error.
    ! Incoming x consists of new breakpoints, not a repeated boundary knot.
    subroutine cpp_extend_scalar(self, c, x, status)
        class(complex_ppoly_type), intent(inout) :: self !! polynomial to extend in place
        complex(dp), intent(in) :: c(:,:) !! coefficients (degree+1, new intervals)
        real(dp), intent(in) :: x(:) !! new breakpoints, same orientation as existing knots
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        complex(dp), allocatable :: cube(:,:,:)
        allocate(cube(size(c,1),size(c,2),1))
        cube(:,:,1)=c
        call cpp_extend_channels(self,cube,x,status)
    end subroutine cpp_extend_scalar

    subroutine cpp_extend_channels(self, c, x, status)
        class(complex_ppoly_type), intent(inout) :: self !! polynomial receiver, modified on success only
        complex(dp), intent(in) :: c(:,:,:) !! coefficients (degree+1, new intervals, channels)
        real(dp), intent(in) :: x(:) !! new breakpoints, sorted with same orientation
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        complex(dp), allocatable :: merged(:,:,:), original(:,:,:), incoming(:,:,:)
        real(dp), allocatable :: extended_x(:)
        integer :: code, old_n, added, new_degree
        logical :: append, ascending
        code=piecewise_status_invalid_input
        if (.not. self%valid()) goto 900
        old_n=size(self%c,2)
        added=size(c,2)
        if (size(x)/=added .or. size(c,1)<1 .or. size(c,3)/=size(self%c,3)) goto 900
        if (added==0) then
            code=piecewise_status_success
            goto 900
        end if
        if (.not. all(ieee_is_finite(x))) goto 900
        if (.not. (all(x(2:added)>=x(1:added-1)) .or. &
            all(x(2:added)<=x(1:added-1)))) goto 900
        ascending=self%x(size(self%x))>=self%x(1)
        if (ascending) then
            if (x(added)<x(1)) goto 900
            if (x(1)>=self%x(old_n+1)) then
                append=.true.
            else if (x(added)<=self%x(1)) then
                append=.false.
            else
                goto 900
            end if
        else
            if (x(added)>x(1)) goto 900
            if (x(1)<=self%x(old_n+1)) then
                append=.true.
            else if (x(added)>=self%x(1)) then
                append=.false.
            else
                goto 900
            end if
        end if
        if (.not. all(ieee_is_finite(real(c,dp)))) goto 900
        if (.not. all(ieee_is_finite(aimag(c)))) goto 900
        new_degree=max(size(self%c,1),size(c,1))
        allocate(merged(new_degree,old_n+added,size(c,3)),extended_x(old_n+added+1))
        merged=cmplx(0.0_dp,0.0_dp,dp)
        if (append) then
            merged(new_degree-size(self%c,1)+1:new_degree,1:old_n,:)=self%c
            merged(new_degree-size(c,1)+1:new_degree,old_n+1:,:)=c
            extended_x(:old_n+1)=self%x
            extended_x(old_n+2:)=x
        else
            merged(new_degree-size(c,1)+1:new_degree,1:added,:)=c
            merged(new_degree-size(self%c,1)+1:new_degree,added+1:,:)=self%c
            extended_x(:added)=x
            extended_x(added+1:)=self%x
        end if
        call move_alloc(extended_x,self%x)
        call move_alloc(merged,self%c)
        code=piecewise_status_success
900     if (present(status)) status=code
    end subroutine cpp_extend_channels

    ! Append or prepend intervals, keeping the receiver unchanged on error.
    ! Incoming x consists of new breakpoints, not a repeated boundary knot.
    subroutine cbp_extend_scalar(self, c, x, status)
        class(complex_bpoly_type), intent(inout) :: self !! polynomial to extend in place
        complex(dp), intent(in) :: c(:,:) !! coefficients (degree+1, new intervals)
        real(dp), intent(in) :: x(:) !! new breakpoints, same orientation as existing knots
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        complex(dp), allocatable :: cube(:,:,:)
        allocate(cube(size(c,1),size(c,2),1))
        cube(:,:,1)=c
        call cbp_extend_channels(self,cube,x,status)
    end subroutine cbp_extend_scalar

    subroutine cbp_extend_channels(self, c, x, status)
        class(complex_bpoly_type), intent(inout) :: self !! polynomial receiver, modified on success only
        complex(dp), intent(in) :: c(:,:,:) !! coefficients (degree+1, new intervals, channels)
        real(dp), intent(in) :: x(:) !! new breakpoints, sorted with same orientation
        integer, intent(out), optional :: status !! 0 success, 1 invalid input
        complex(dp), allocatable :: merged(:,:,:), original(:,:,:), incoming(:,:,:)
        real(dp), allocatable :: extended_x(:)
        integer :: code, old_n, added, new_degree
        logical :: append, ascending
        code=piecewise_status_invalid_input
        if (.not. self%valid()) goto 900
        old_n=size(self%c,2)
        added=size(c,2)
        if (size(x)/=added .or. size(c,1)<1 .or. size(c,3)/=size(self%c,3)) goto 900
        if (added==0) then
            code=piecewise_status_success
            goto 900
        end if
        if (.not. all(ieee_is_finite(x))) goto 900
        if (.not. (all(x(2:added)>=x(1:added-1)) .or. &
            all(x(2:added)<=x(1:added-1)))) goto 900
        ascending=self%x(size(self%x))>=self%x(1)
        if (ascending) then
            if (x(added)<x(1)) goto 900
            if (x(1)>=self%x(old_n+1)) then
                append=.true.
            else if (x(added)<=self%x(1)) then
                append=.false.
            else
                goto 900
            end if
        else
            if (x(added)>x(1)) goto 900
            if (x(1)<=self%x(old_n+1)) then
                append=.true.
            else if (x(added)>=self%x(1)) then
                append=.false.
            else
                goto 900
            end if
        end if
        if (.not. all(ieee_is_finite(real(c,dp)))) goto 900
        if (.not. all(ieee_is_finite(aimag(c)))) goto 900
        new_degree=max(size(self%c,1),size(c,1))
        allocate(merged(new_degree,old_n+added,size(c,3)),extended_x(old_n+added+1))
        merged=cmplx(0.0_dp,0.0_dp,dp)
        original=raise_bernstein_complex(self%c,new_degree)
        incoming=raise_bernstein_complex(c,new_degree)
        if (append) then
            merged(new_degree-size(original,1)+1:new_degree,1:old_n,:)=original
            merged(new_degree-size(incoming,1)+1:new_degree,old_n+1:,:)=incoming
            extended_x(:old_n+1)=self%x
            extended_x(old_n+2:)=x
        else
            merged(new_degree-size(incoming,1)+1:new_degree,1:added,:)=incoming
            merged(new_degree-size(original,1)+1:new_degree,added+1:,:)=original
            extended_x(:added)=x
            extended_x(added+1:)=self%x
        end if
        call move_alloc(extended_x,self%x)
        call move_alloc(merged,self%c)
        code=piecewise_status_success
900     if (present(status)) status=code
    end subroutine cbp_extend_channels

    ! Stable degree elevation, performed one Bernstein degree at a time.
    function raise_bernstein_real(c, new_order) result(elevated)
        real(dp), intent(in) :: c(:,:,:) !! original Bernstein coefficients
        integer, intent(in) :: new_order !! output number of coefficients >= input count
        real(dp), allocatable :: elevated(:,:,:), updated(:,:,:)
        real(dp) :: alpha
        integer :: current, j
        elevated=c
        do current=size(c,1),new_order-1
            allocate(updated(current+1,size(c,2),size(c,3)))
            updated(1,:,:)=elevated(1,:,:)
            updated(current+1,:,:)=elevated(current,:,:)
            do j=2,current
                alpha=real(j-1,dp)/real(current,dp)
                updated(j,:,:)=alpha*elevated(j-1,:,:)+(1.0_dp-alpha)*elevated(j,:,:)
            end do
            call move_alloc(updated,elevated)
        end do
    end function raise_bernstein_real

    ! Stable degree elevation, performed one Bernstein degree at a time.
    function raise_bernstein_complex(c, new_order) result(elevated)
        complex(dp), intent(in) :: c(:,:,:) !! original Bernstein coefficients
        integer, intent(in) :: new_order !! output number of coefficients >= input count
        complex(dp), allocatable :: elevated(:,:,:), updated(:,:,:)
        real(dp) :: alpha
        integer :: current, j
        elevated=c
        do current=size(c,1),new_order-1
            allocate(updated(current+1,size(c,2),size(c,3)))
            updated(1,:,:)=elevated(1,:,:)
            updated(current+1,:,:)=elevated(current,:,:)
            do j=2,current
                alpha=real(j-1,dp)/real(current,dp)
                updated(j,:,:)=alpha*elevated(j-1,:,:)+(1.0_dp-alpha)*elevated(j,:,:)
            end do
            call move_alloc(updated,elevated)
        end do
    end function raise_bernstein_complex

end module scifort_interpolate_piecewise
