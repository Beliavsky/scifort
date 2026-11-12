! SPDX-License-Identifier: MIT
! Same source can benchmark the previous and the indexed Clough-Tocher builds.
program benchmark_clough_tocher_driver
    use scifort_kinds, only: dp
    use scifort_interpolate_clough_tocher, only: clough_tocher_2d_interpolator, &
        clough_tocher_2d_interpolator_type, ct_status_success
    implicit none
    integer, parameter :: n=450, nq=9000, nch=4, repetitions=4
    real(dp) :: x(n,2), y(n,nch), q(nq,2), t0,t1,t2,build_time,checksum
    real(dp), allocatable :: ans(:,:)
    integer, allocatable :: triangles(:,:)
    integer :: i,j,t,status,r
    type(clough_tocher_2d_interpolator_type) :: interp
    do i=1,n
        x(i,1)=modulo(sin(12.9898_dp*real(i,dp))*43758.5453_dp,1.0_dp)
        x(i,2)=modulo(sin(78.233_dp*real(i,dp))*12475.231_dp,1.0_dp)
        do j=1,nch
            y(i,j)=sin(x(i,1)*real(j,dp)*2.2_dp+2.0_dp*x(i,2))+cos(0.15_dp*real(i*j,dp))
        end do
    end do
    call cpu_time(t0)
    interp=clough_tocher_2d_interpolator(x,y,tol=1.0e-10_dp, &
        maxiter=600,fill_value=-3.125_dp,status=status)
    call cpu_time(t1)
    build_time=t1-t0
    if (status/=ct_status_success .or. .not.interp%valid()) error stop 'benchmark construction'
    triangles=interp%simplices()
    do i=1,nq
        if (mod(i,7)==0) then
            q(i,:)=[1.5_dp,-0.5_dp]
        else
            t=mod(37*i+31,size(triangles,2))+1
            q(i,:)=0.2_dp*x(triangles(1,t),:)+0.3_dp*x(triangles(2,t),:)+ &
                0.5_dp*x(triangles(3,t),:)
        end if
    end do
    call cpu_time(t2)
    do r=1,repetitions
        ans=interp%evaluate(q,status=status)
        if (status/=ct_status_success) error stop 'benchmark evaluation'
    end do
    call cpu_time(t0)
    checksum=sum(ans)
    print '(a,es24.16)', 'build_sec=',build_time
    print '(a,es24.16)', 'evaluate_sec=',(t0-t2)/real(repetitions,dp)
    print '(a,es24.16)', 'checksum=',checksum
    print '(a,i0)', 'triangles=',size(triangles,2)
end program benchmark_clough_tocher_driver
