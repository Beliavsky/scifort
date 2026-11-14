! SPDX-License-Identifier: MIT
program test_interpolate_fitpack_tck
    use scifort_kinds, only : dp
    use scifort_interpolate, only : fitpack_tck_type, splrep, splev, splint, &
        splder, splantider, spalde, fitpack_tck_success, &
        fitpack_tck_out_of_bounds
    implicit none
    integer, parameter :: n = 14
    real(dp) :: x(n), y(n), w(n), d, val, expected
    real(dp), allocatable :: ders(:)
    type(fitpack_tck_type) :: tck, derivative_tck, antiderivative_tck, lsq
    integer :: i, status

    do i=1,n
        x(i) = 0.25_dp*real(i-1,dp)
        y(i) = 2.0_dp*x(i)**3-3.0_dp*x(i)**2+1.0_dp
        w(i) = 1.0_dp+0.1_dp*real(i,dp)
    end do
    tck=splrep(x,y,k=3,s=0.0_dp,status=status)
    if (status /= fitpack_tck_success .or. .not. tck%valid()) error stop 'splrep exact'
    if (size(tck%c) /= size(tck%t)) error stop 'FITPACK padded tuple'
    do i=1,n
        if (abs(splev(x(i),tck)-y(i)) > 1.0e-10_dp) error stop 'exact interpolation'
        if (abs(splev(x(i),tck,der=1)-(6.0_dp*x(i)**2-6.0_dp*x(i))) > 1.0e-9_dp) &
            error stop 'spline first derivative'
    end do
    val=splev(-0.75_dp,tck,ext=1)
    if (abs(val) > 1.0e-14_dp) error stop 'zero outside'
    val=splev(-0.75_dp,tck,ext=3)
    if (abs(val-splev(x(1),tck)) > 1.0e-14_dp) error stop 'clamp outside'
    val=splev(-0.75_dp,tck,ext=2,status=status)
    if (status /= fitpack_tck_out_of_bounds) error stop 'out-of-bounds error'
    ders=spalde(0.75_dp,tck,status=status)
    if (status /= 0 .or. size(ders) /= 4) error stop 'spalde size'
    if (abs(ders(3)-(12.0_dp*0.75_dp-6.0_dp)) > 1.0e-9_dp) error stop 'spalde second'
    derivative_tck=splder(tck,status=status)
    if (status /= 0) error stop 'splder'
    antiderivative_tck=splantider(tck,status=status)
    if (status /= 0) error stop 'splantider'
    d=0.75_dp
    if (abs(splev(d,derivative_tck)-(6.0_dp*d*d-6.0_dp*d)) > 1.0e-9_dp) &
        error stop 'derived spline'
    if (abs(splev(d,antiderivative_tck)-splev(0.0_dp,antiderivative_tck)- &
        splint(0.0_dp,d,tck)) > 1.0e-9_dp) error stop 'antiderivative'
    expected=(0.5_dp*x(n)**4-x(n)**3+x(n))
    if (abs(splint(-100.0_dp,100.0_dp,tck)-expected) > 1.0e-9_dp) error stop 'clipped integral'
    lsq=splrep(x,y,w=w,k=3,task=-1,t=[1.0_dp,2.0_dp],status=status)
    if (status /= 0) error stop 'fixed-knot LSQ'
    if (maxval(abs(splev(x,lsq)-y)) > 1.0e-10_dp) error stop 'LSQ cubic reproduction'
    lsq=splrep(x,y,s=1.0_dp,status=status)
    if (status /= 0 .or. .not. lsq%valid()) error stop 'native positive smoothing'
    if (size(lsq%t) /= 8) error stop 'polynomial smoothing limit'
    print *, 'FITPACK TCK compatibility PASS'
end program test_interpolate_fitpack_tck
