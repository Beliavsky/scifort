! SPDX-License-Identifier: MIT
program test_interpolate_legacy
    use scifort_kinds, only:dp
    use scifort_interpolate, only:interp1d,interp1d_type,complex_interp1d_type, &
        interp1d_success,interp1d_invalid_input,interp1d_out_of_bounds
    implicit none
    type(interp1d_type) :: f
    type(complex_interp1d_type) :: g
    real(dp) :: x(5),y(5,2),q(6),poly
    complex(dp) :: z(5,2)
    real(dp), allocatable :: a(:,:),one(:)
    complex(dp), allocatable :: c(:,:)
    integer :: i,code
    x=[2.0_dp,0.0_dp,4.0_dp,1.0_dp,3.0_dp]
    y(:,1)=x**2
    y(:,2)=3.0_dp+x**2
    g=interp1d(x,cmplx(y,2*y,dp),kind='quadratic',status=code)
    if (code/=interp1d_success .or. .not.g%valid()) error stop 'complex interp1d create'
    q=[0.0_dp,0.5_dp,1.5_dp,2.3_dp,3.5_dp,4.0_dp]
    c=g%evaluate(q,status=code)
    if (code/=0) error stop 'complex interp1d evaluation status'
    do i=1,size(q)
        if (abs(c(i,1)-cmplx(q(i)**2,2*q(i)**2,dp))>3.0e-11_dp) &
            error stop 'quadratic interp1d exact polynomial'
        if (abs(c(i,2)-cmplx(3.0_dp+q(i)**2,2*(3.0_dp+q(i)**2),dp))>3.0e-11_dp) &
            error stop 'multichannel quadratic interp1d'
    end do
    f=interp1d([0.0_dp,1.0_dp,2.0_dp],[10.0_dp,20.0_dp,40.0_dp], &
        kind='nearest',extrapolate=.true.,status=code)
    if (code/=0 .or. abs(f%value(0.5_dp)-10.0_dp)>1.0e-12_dp) &
        error stop 'nearest tie down'
    f=interp1d([0.0_dp,1.0_dp,2.0_dp],[10.0_dp,20.0_dp,40.0_dp], &
        kind='nearest-up',extrapolate=.true.,status=code)
    if (code/=0 .or. abs(f%value(0.5_dp)-20.0_dp)>1.0e-12_dp) &
        error stop 'nearest tie up'
    f=interp1d([0.0_dp,1.0_dp,2.0_dp],[10.0_dp,20.0_dp,40.0_dp], &
        kind='previous',extrapolate=.true.,status=code)
    if (code/=0 .or. abs(f%value(4.0_dp)-40.0_dp)>1.0e-12_dp) &
        error stop 'previous high extrapolation'
    f=interp1d([0.0_dp,1.0_dp,2.0_dp],[10.0_dp,20.0_dp,40.0_dp], &
        kind='next',extrapolate=.true.,status=code)
    if (code/=0 .or. abs(f%value(-1.0_dp)-10.0_dp)>1.0e-12_dp) &
        error stop 'next low extrapolation'
    f=interp1d([2.0_dp,0.0_dp,1.0_dp],[4.0_dp,0.0_dp,1.0_dp], &
        kind='linear',bounds_error=.false.,fill_below=-9.0_dp,fill_above=99.0_dp,status=code)
    if (code/=0 .or. f%channels()/=1) error stop 'real interp1d metadata'
    a=f%evaluate([-1.0_dp,0.5_dp,4.0_dp],status=code)
    if (code/=0 .or. any(abs(a(:,1)-[-9.0_dp,0.5_dp,99.0_dp])>1.0e-12_dp)) &
        error stop 'real interp1d fills and unsorted x'
    f=interp1d(x,y,kind='cubic',extrapolate=.true.,status=code)
    if (code/=0) error stop 'cubic interpolation fit'
    if (abs(f%value(4.5_dp)-4.5_dp**2)>1.0e-10_dp) error stop 'cubic extrapolation'
    f=interp1d(x,y,kind='cubic',status=code)
    poly=f%value(-1.0_dp,status=code)
    if (code/=interp1d_out_of_bounds) error stop 'bounds_error default'
    f=interp1d([0.0_dp,0.0_dp,1.0_dp],[1.0_dp,2.0_dp,3.0_dp],status=code)
    if (code/=interp1d_invalid_input .or. f%valid()) error stop 'duplicate knots accepted'
    f=interp1d([2.0_dp,0.0_dp,1.0_dp],[4.0_dp,0.0_dp,1.0_dp], &
        assume_sorted=.true.,status=code)
    if (code/=interp1d_invalid_input .or. f%valid()) error stop 'assume_sorted error'
    print '(a)', 'interp1d legacy modes, fills, bounds, real/complex PASS'
end program test_interpolate_legacy
