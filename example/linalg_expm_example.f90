! SPDX-License-Identifier: MIT
! Build with kinds.f90, linalg.f90, linalg_dense.f90, linalg_expm.f90.
program linalg_expm_example
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_expm, only : expm
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: a(2,2), x(2,2), batches(2,2,2), results(2,2,2)
    complex(dp) :: z(2,2), y(2,2)
    integer :: status, squarings

    a = 0.0_dp
    a(1,2) = -pi/2.0_dp
    a(2,1) = pi/2.0_dp
    call expm(a,x,status,scalings=squarings)
    if (status /= linalg_status_success) error stop 'expm real failed'
    print '(a)', 'exp([0, -pi/2; pi/2, 0]):'
    print '(2f12.7)', x(1,:)
    print '(2f12.7)', x(2,:)
    print '(a,i0)', 'Scaling steps: ',squarings

    z = cmplx(a,0.4_dp*a,dp)
    call expm(z,y,status)
    if (status /= linalg_status_success) error stop 'expm complex failed'
    print '(a,2es12.4)', 'complex exp(A)(1,1): ',real(y(1,1),dp),aimag(y(1,1))

    batches(1,:,:) = a
    batches(2,:,:) = -a
    call expm(batches,results,status)
    if (status /= linalg_status_success) error stop 'expm batch failed'
    print '(a,f12.7)', 'batch exp(-A)(1,2): ',results(2,1,2)
end program linalg_expm_example
