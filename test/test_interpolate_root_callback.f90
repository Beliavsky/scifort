! SPDX-License-Identifier: MIT
program test_interpolate_root_callback
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_interpolate, only : ppoly, ppoly_type, piecewise_native_real_roots
    implicit none
    type(ppoly_type) :: pp
    real(dp), allocatable :: roots(:), usual(:)
    integer :: calls
    calls=0
    pp=ppoly(reshape([1.0_dp,-1.0_dp,0.0_dp],[3,1]),[0.0_dp,2.0_dp],extrapolate=.false.)
    roots=pp%roots(root_solver=tracking_roots)
    usual=pp%roots()
    if (calls/=1 .or. size(roots)/=size(usual)) error stop 'callback invocation'
    if (maxval(abs(roots-usual))>1.0e-13_dp) error stop 'callback numerical roots'
    pp=ppoly(reshape([0.0_dp,1.0_dp],[1,2]),[0.0_dp,1.0_dp,2.0_dp],extrapolate=.false.)
    roots=pp%roots(root_solver=tracking_roots)
    if (size(roots)/=2 .or. roots(1)/=0.0_dp) error stop 'zero section behavior'
    if (.not.ieee_is_nan(roots(2))) error stop 'zero section NaN marker'
    print '(a)', 'optional PPoly root callback regression PASS'
contains
    function tracking_roots(coeff,lower,upper) result(r)
        real(dp), intent(in) :: coeff(:) !! descending polynomial coefficients
        real(dp), intent(in) :: lower !! local interval minimum
        real(dp), intent(in) :: upper !! local interval maximum
        real(dp), allocatable :: r(:)
        calls=calls+1
        r=piecewise_native_real_roots(coeff,lower,upper)
    end function tracking_roots
end program test_interpolate_root_callback
