! SPDX-License-Identifier: MIT
program test_signal_sos_responses
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_discrete_simulation, only : dstep, dimpulse, dstep_sos, dimpulse_sos, &
        discrete_response_result, discrete_sim_success, discrete_sim_invalid_input
    use scifort_signal_sos, only : sosfilt
    implicit none
    real(dp), parameter :: tol = 2.0e-13_dp
    real(dp) :: sos(2,6), bad(2,6), x(19), zi(2,2)
    complex(dp) :: csos(2,6), czi(2,2)
    real(dp), allocatable :: yr(:)
    complex(dp), allocatable :: yc(:)
    type(discrete_response_result) :: a, b, impulse, step
    integer :: j

    sos(1,:) = [0.24_dp,0.30_dp,0.10_dp,1.0_dp,-0.4_dp,0.12_dp]
    sos(2,:) = [0.75_dp,-0.04_dp,0.0_dp,1.0_dp,0.2_dp,0.05_dp]
    zi(1,:) = [0.03_dp,-0.02_dp]
    zi(2,:) = [-0.01_dp,0.02_dp]
    a = dstep_sos(sos,0.2_dp,n=19)
    b = dstep(sos,0.2_dp,n=19)
    call check(a%status==discrete_sim_success,'real step status')
    call check(b%status==discrete_sim_success,'generic real step')
    call check(all(shape(a%y)==[19,1,1]),'shape')
    call check(maxval(abs(a%t-[(0.2_dp*real(j-1,dp),j=1,19)]))<tol,'times')
    call check(maxval(abs(a%y-b%y))<tol,'generic alias')
    yr = sosfilt(sos,[(1.0_dp,j=1,19)])
    call check(maxval(abs(a%y(:,1,1)-yr))<tol,'real step vs filtering')
    impulse = dimpulse_sos(sos,1.0_dp,n=19)
    step = dstep_sos(sos,1.0_dp,n=19)
    call check(impulse%status==discrete_sim_success,'real impulse')
    call check(maxval(abs(impulse%y(:,1,1)-sosfilt(sos, &
        [1.0_dp,(0.0_dp,j=2,19)])))<tol,'real impulse vs filtering')
    call check(maxval(abs(step%y(:,1,1)-sum_prefix(impulse%y(:,1,1))))<tol, &
        'step is cumulative impulse')
    a = dstep_sos(sos(1:2,:),0.1_dp,n=19,zi=zi)
    yr = sosfilt(sos,[(1.0_dp,j=1,19)],zi=zi)
    call check(maxval(abs(a%y(:,1,1)-yr))<tol,'explicit real state')

    csos = cmplx(sos,0.0_dp,dp)
    csos(1,2) = cmplx(0.30_dp,0.07_dp,dp)
    csos(2,5) = cmplx(0.20_dp,-0.15_dp,dp)
    czi = cmplx(zi,0.04_dp,dp)
    a = dstep(csos,0.125_dp,n=19,zi=czi)
    yc = sosfilt(csos,[(cmplx(1.0_dp,0.0_dp,dp),j=1,19)],zi=czi)
    call check(a%status==discrete_sim_success,'complex step status')
    call check(maxval(abs(a%y(:,1,1)-yc))<tol,'complex step vs filtering')
    b = dimpulse(csos,0.125_dp,n=19,zi=czi)
    yc = sosfilt(csos,[cmplx(1.0_dp,0.0_dp,dp), &
        (cmplx(0.0_dp,0.0_dp,dp),j=2,19)],zi=czi)
    call check(b%status==discrete_sim_success,'complex impulse status')
    call check(maxval(abs(b%y(:,1,1)-yc))<tol,'complex impulse vs filtering')

    a = dstep_sos(sos,-1.0_dp,n=12)
    call check(a%status==discrete_sim_invalid_input,'negative dt')
    call check(.not. allocated(a%y),'negative dt unallocated')
    a = dimpulse_sos(sos,0.1_dp,n=0)
    call check(a%status==discrete_sim_invalid_input,'nonpositive n')
    a = dstep_sos(sos,ieee_value(0.0_dp,ieee_quiet_nan))
    call check(a%status==discrete_sim_invalid_input,'nan dt')
    bad = sos
    bad(2,4) = 2.0_dp
    a = dstep_sos(bad,0.5_dp,n=4)
    call check(a%status==discrete_sim_invalid_input,'unnormalized a0')
    a = dimpulse_sos(sos(:,:5),0.5_dp,n=4)
    call check(a%status==discrete_sim_invalid_input,'wrong columns')
    a = dstep_sos(sos,0.5_dp,n=7,zi=zi(:1,:))
    call check(a%status==discrete_sim_invalid_input,'wrong real state')
    a = dimpulse_sos(csos,0.5_dp,n=7,zi=czi(:1,:))
    call check(a%status==discrete_sim_invalid_input,'wrong complex state')
    call check(.not. allocated(a%y),'invalid states preserve empty result')
    print '(a)', 'PASS native SOS step and impulse responses'
contains
    subroutine check(ok,message)
        logical, intent(in) :: ok !! assertion result
        character(*), intent(in) :: message !! assertion description
        if (.not. ok) then
            print '(2a)', 'FAIL: ',message
            error stop 1
        end if
    end subroutine check
    function sum_prefix(v) result(out)
        complex(dp), intent(in) :: v(:) !! impulse response
        complex(dp) :: out(size(v))
        integer :: k
        out(1) = v(1)
        do k=2,size(v)
            out(k) = out(k-1)+v(k)
        end do
    end function sum_prefix
end program test_signal_sos_responses
