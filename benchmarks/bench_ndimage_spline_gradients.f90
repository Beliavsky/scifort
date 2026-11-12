! SPDX-License-Identifier: MIT
! Compare analytic cached-spline derivatives with coordinate finite differences.
program bench_ndimage_spline_gradients
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_interpolation_plan, &
        prepare_spline_plan_nd, evaluate_spline_plan_nd, &
        evaluate_spline_plan_with_gradient_nd
    implicit none
    integer, parameter :: npts=900,nruns=4
    real(dp) :: x(64,60),coords(2,npts),perturbed(2,npts)
    real(dp) :: values(npts),grad(2,npts),fd(2,npts),yplus(npts),yminus(npts)
    real(dp) :: t0,t1,analytic_seconds,finite_seconds,step,err
    type(spline_interpolation_plan) :: plan
    integer :: i,j,k,n,axis,code

    do j=1,60
        do i=1,64
            x(i,j)=sin(0.17_dp*i+0.24_dp*j)+0.003_dp*i*j
        end do
    end do
    do k=1,npts
        coords(:,k)=[1.3_dp+62.4_dp*real(mod(37*k,900),dp)/900.0_dp, &
                     1.4_dp+57.2_dp*real(mod(101*k,901),dp)/901.0_dp]
    end do
    call prepare_spline_plan_nd(x,plan,order=5,mode='mirror',status=code)
    if (code/=0) error stop 'prepare failed'
    step=1e-5_dp
    call cpu_time(t0)
    do n=1,nruns
        call evaluate_spline_plan_with_gradient_nd(plan,coords,values,grad,code)
        if (code/=0) error stop 'analytic failed'
    end do
    call cpu_time(t1)
    analytic_seconds=t1-t0
    call cpu_time(t0)
    do n=1,nruns
        do axis=1,2
            perturbed=coords
            perturbed(axis,:)=coords(axis,:)+step
            call evaluate_spline_plan_nd(plan,perturbed,yplus,code)
            if (code/=0) error stop 'positive sample failed'
            perturbed(axis,:)=coords(axis,:)-step
            call evaluate_spline_plan_nd(plan,perturbed,yminus,code)
            if (code/=0) error stop 'negative sample failed'
            fd(axis,:)=(yplus-yminus)/(2.0_dp*step)
        end do
    end do
    call cpu_time(t1)
    finite_seconds=t1-t0
    err=maxval(abs(grad-fd))
    if (err>1e-5_dp) error stop 'derivative comparison failed'
    print '(a,es12.4)', 'Analytic seconds: ',analytic_seconds
    print '(a,es12.4)', 'Finite-difference seconds: ',finite_seconds
    print '(a,f9.3)', 'Speedup: ',finite_seconds/max(analytic_seconds,tiny(1.0_dp))
    print '(a,es12.4)', 'Maximum gradient error: ',err
end program bench_ndimage_spline_gradients
