! SPDX-License-Identifier: MIT
! Balanced truncation for complex MIMO stable systems.
program test_signal_dlti_balanced
    use scifort_signal, only : dp, dlti_system, dlti_ss, dlti_balanced_result, &
        dlti_success, dlti_invalid_input, dlti_numerical_failure, &
        mimo_frequency_response
    implicit none
    complex(dp) :: a(3,3), b(3,3), c(3,3), d(3,3), &
        mixing(3,3), ag(3,3), bg(3,3), cg(3,3)
    real(dp) :: hsv_true(3), w(15), pi, maxerr
    type(dlti_system) :: sys, bad, zero
    type(dlti_balanced_result) :: reduced, complete, no_state
    type(mimo_frequency_response) :: fr_full, fr_small
    integer :: i,j

    pi = acos(-1.0_dp)
    a = (0.0_dp,0.0_dp)
    b = (0.0_dp,0.0_dp)
    c = (0.0_dp,0.0_dp)
    d = (0.0_dp,0.0_dp)
    a(1,1) = cmplx(0.15_dp,0.10_dp,dp)
    a(2,2) = cmplx(0.40_dp,-0.20_dp,dp)
    a(3,3) = cmplx(0.65_dp,0.0_dp,dp)
    b(1,1) = (1.2_dp,0.0_dp)
    b(2,2) = (0.7_dp,0.0_dp)
    b(3,3) = (0.2_dp,0.0_dp)
    c(1,1) = (0.5_dp,0.0_dp)
    c(2,2) = (0.8_dp,0.0_dp)
    c(3,3) = (0.2_dp,0.0_dp)
    d(1,1) = (0.2_dp,0.1_dp)
    d(2,2) = (0.0_dp,-0.2_dp)
    hsv_true = [0.6_dp/(1.0_dp-abs(a(1,1))**2), &
        0.56_dp/(1.0_dp-abs(a(2,2))**2), &
        0.04_dp/(1.0_dp-abs(a(3,3))**2)]
    hsv_true = [hsv_true(2),hsv_true(1),hsv_true(3)]
    ! Mix internal states by a unitary DFT; Hankel singular values are invariant.
    do j=1,3
        do i=1,3
            mixing(i,j) = exp(cmplx(0.0_dp,2.0_dp*pi*real((i-1)*(j-1),dp)/3.0_dp,dp))/sqrt(3.0_dp)
        end do
    end do
    ag = matmul(matmul(mixing,a),transpose(conjg(mixing)))
    bg = matmul(mixing,b)
    cg = matmul(c,transpose(conjg(mixing)))
    sys = dlti_ss(ag,bg,cg,d,0.25_dp)
    call check(sys%status == dlti_success,'complex system constructor')
    reduced = sys%balanced_truncate(2)
    call check(reduced%status == dlti_success,'rank 2 balanced model')
    call check(reduced%model%nstates() == 2,'retained exactly two states')
    call check(reduced%model%ninputs() == 3,'preserved input channels')
    call check(reduced%model%noutputs() == 3,'preserved output channels')
    call check(reduced%model%dt == 0.25_dp,'preserved dt')
    call check(maxval(abs(reduced%hsv-hsv_true)) < 1.0e-10_dp,'analytic singular values')
    call check(abs(reduced%error_bound-2.0_dp*hsv_true(3)) < 1.0e-10_dp,'error bound')
    call check(maxval(abs(reduced%model%d-d)) < 1.0e-14_dp,'feedthrough preserved')
    do i=1,size(w)
        w(i) = pi*real(i,dp)/real(size(w)+1,dp)
    end do
    fr_full = sys%freqresp_mimo(w=w,solver='direct')
    fr_small = reduced%model%freqresp_mimo(w=w,solver='direct')
    call check(fr_full%status == 0 .and. fr_small%status == 0,'response status')
    maxerr = maxval(abs(fr_small%response-fr_full%response))
    call check(maxerr <= reduced%error_bound+1.0e-9_dp,'discrete H-infinity bound')
    complete = sys%balanced_truncate(3)
    call check(complete%status == dlti_success,'full balancing')
    fr_small = complete%model%freqresp_mimo(w=w,solver='direct')
    call check(maxval(abs(fr_small%response-fr_full%response)) < 1.0e-10_dp, &
        'full-order input-output equivalence')
    no_state = sys%balanced_truncate(0)
    call check(no_state%status == dlti_success,'zero-state reduction')
    call check(no_state%model%nstates() == 0,'zero-state result')
    call check(maxval(abs(no_state%model%d-d)) < 1.0e-15_dp,'static feedthrough')
    complete = sys%balanced_truncate(4)
    call check(complete%status == dlti_invalid_input,'invalid high order')
    complete = sys%balanced_truncate(-1)
    call check(complete%status == dlti_invalid_input,'invalid negative order')
    complete = sys%balanced_truncate(1,tol=-1.0_dp)
    call check(complete%status == dlti_invalid_input,'invalid tolerance')
    a(3,3) = (1.1_dp,0.0_dp)
    bad = dlti_ss(a,b,c,d)
    complete = bad%balanced_truncate(1)
    call check(complete%status == dlti_numerical_failure,'unstable rejected')
    zero = dlti_ss(ag(:0,:0),bg(:0,:),cg(:,:0),d)
    complete = zero%balanced_truncate(0)
    call check(complete%status == dlti_success,'empty state')
    call check(size(complete%hsv) == 0,'empty Hankel spectrum')
    print '(a)', 'PASS: balanced truncation, MIMO, complex, stability, ranks'
contains
    subroutine check(condition,label)
        logical, intent(in) :: condition !! assertion to test
        character(*), intent(in) :: label !! short diagnostic
        if (.not. condition) then
            print '(a,1x,a)', 'FAIL:', label
            error stop 1
        end if
    end subroutine check
end program test_signal_dlti_balanced
