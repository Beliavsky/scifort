! SPDX-License-Identifier: MIT
! Exercise dependency-free SciPy-compatible parametric splprep TCK support.
program test_interpolate_fitpack_parametric
    use scifort_kinds, only : dp
    use scifort_interpolate, only : splprep, splev, fitpack_parametric_tck_type, &
        fitpack_tck_success, fitpack_tck_invalid_input, fitpack_tck_unsupported_smoothing, &
        fitpack_tck_out_of_bounds
    implicit none
    integer, parameter :: n=25
    real(dp) :: u(n), x(n,3), w(n), q(17), trueval(3), derivative(3)
    real(dp), allocatable :: values(:,:), v(:), t(:), fitted_u(:)
    type(fitpack_parametric_tck_type) :: curve, fixed, bad
    integer :: i, code

    do i=1,n
        u(i) = real(i-1,dp)/real(n-1,dp)
        x(i,1) = 1.0_dp+2.0_dp*u(i)+0.4_dp*u(i)**3
        x(i,2) = 2.0_dp-3.0_dp*u(i)+0.7_dp*u(i)**2
        x(i,3) = -1.0_dp+0.5_dp*u(i)**3
        w(i) = 0.5_dp+0.1_dp*real(i,dp)
    end do
    q = [(real(i-2,dp)/real(size(q)-3,dp),i=1,size(q))]
    curve = splprep(x,u=u,k=3,s=0.0_dp,status=code)
    if (code /= fitpack_tck_success .or. .not. curve%valid()) error stop 'splprep exact fit'
    if (curve%channels() /= 3) error stop 'splprep channels'
    if (size(curve%c,1) /= n .or. size(curve%t) /= n+4) error stop 'splprep TCK shape'
    if (maxval(abs(curve%u-u)) > 1.0e-15_dp) error stop 'stored explicit parameter'
    values = splev(q,curve,status=code)
    if (code /= 0 .or. any(shape(values) /= [size(q),3])) error stop 'vector evaluation shape'
    do i=1,size(q)
        trueval = [1.0_dp+2.0_dp*q(i)+0.4_dp*q(i)**3, &
            2.0_dp-3.0_dp*q(i)+0.7_dp*q(i)**2, -1.0_dp+0.5_dp*q(i)**3]
        if (maxval(abs(values(i,:)-trueval)) > 5.0e-11_dp) error stop 'parametric cubic reproduction'
        derivative = [2.0_dp+1.2_dp*q(i)**2, -3.0_dp+1.4_dp*q(i),1.5_dp*q(i)**2]
        v=splev(q(i),curve,der=1,status=code)
        if (code /= 0 .or. maxval(abs(v-derivative)) > 5.0e-11_dp) &
            error stop 'parametric derivative'
    end do
    v=splev(-0.2_dp,curve,ext=1,status=code)
    if (code /= 0 .or. maxval(abs(v)) > 1.0e-15_dp) error stop 'zero outside'
    v=splev(-0.2_dp,curve,ext=3,status=code)
    if (code /= 0 .or. maxval(abs(v-x(1,:))) > 5.0e-13_dp) error stop 'clamp outside'
    v=splev(-0.2_dp,curve,ext=2,status=code)
    if (code /= fitpack_tck_out_of_bounds) error stop 'out-of-domain failure'

    ! FITPACK's task=-1 takes a complete knot vector, not only the interior knots.
    t = [spread(u(1),1,4),0.25_dp,0.5_dp,0.75_dp,spread(u(n),1,4)]
    fixed = splprep(x,u=u,w=w,k=3,task=-1,t=t,status=code)
    if (code /= 0 .or. .not. fixed%valid()) error stop 'fixed knot fit'
    if (maxval(abs(fixed%t-t)) > 1.0e-14_dp) error stop 'fixed knot vector'
    values = splev(u,fixed)
    if (maxval(abs(values-x)) > 2.0e-11_dp) error stop 'weighted cubic reproduction'

    ! Default parameterization is accumulated Euclidean chord length.
    x(:,1)=u
    x(:,2)=2.0_dp*u
    x(:,3)=-u
    curve = splprep(x,k=3,status=code)
    if (code /= 0) error stop 'chord length fit'
    if (maxval(abs(curve%u-u)) > 1.0e-13_dp) error stop 'chord length parameter'
    v = splev(0.4_dp,curve)
    if (maxval(abs(v-[0.4_dp,0.8_dp,-0.4_dp])) > 1.0e-12_dp) &
        error stop 'chord length eval'

    bad=splprep(x,u=u,s=1.0_dp,status=code)
    if (code /= fitpack_tck_unsupported_smoothing .or. bad%valid()) &
        error stop 'positive smoothing must not be approximated'
    bad=splprep(x,u=u,per=.true.,status=code)
    if (code /= fitpack_tck_unsupported_smoothing .or. bad%valid()) &
        error stop 'periodic fit must not be approximated'
    bad=splprep(x,u=u,k=6,status=code)
    if (code /= fitpack_tck_invalid_input .or. bad%valid()) error stop 'invalid degree'
    bad=splprep(x,u=u,w=w(:n-1),status=code)
    if (code /= fitpack_tck_invalid_input) error stop 'invalid weights'
    fitted_u=u
    fitted_u(8)=fitted_u(7)
    bad=splprep(x,u=fitted_u,status=code)
    if (code /= fitpack_tck_invalid_input) error stop 'duplicate parameter'
    x(6,:)=x(5,:)
    bad=splprep(x,status=code)
    if (code /= fitpack_tck_invalid_input) error stop 'zero chord length'
    print *, 'FITPACK parametric splprep PASS'
end program test_interpolate_fitpack_parametric
