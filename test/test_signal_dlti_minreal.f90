! SPDX-License-Identifier: MIT
! Controllability/observability state reduction, including complex similarity.
program test_signal_dlti_minreal
    use scifort_kinds, only : dp
    use scifort_signal_dlti, only : dlti_system, dlti_ss, dlti_gramian_result, &
        dlti_success, dlti_invalid_input, dlti_numerical_failure
    use scifort_signal_discrete_simulation, only : discrete_simulation_result
    use scifort_signal_lti_frequency, only : mimo_frequency_response
    implicit none
    type(dlti_system) :: original,reduced,bad,zero
    type(dlti_gramian_result) :: gram
    type(discrete_simulation_result) :: simulation_a,simulation_b
    type(mimo_frequency_response) :: fr_a,fr_b
    complex(dp) :: a(4,4),b(4,2),c(2,4),d(2,2), &
        unitary(4,4),diag(4,4),base_b(4,2),base_c(2,4)
    complex(dp) :: u(32,2)
    real(dp) :: pi,w(8)
    integer :: i,j,k

    pi = acos(-1.0_dp)
    unitary = (0.0_dp,0.0_dp)
    do j=1,4
        do i=1,4
            unitary(i,j) = exp(cmplx(0.0_dp,2.0_dp*pi*real((i-1)*(j-1),dp)/4.0_dp,dp))/2.0_dp
        end do
    end do
    diag = (0.0_dp,0.0_dp)
    diag(1,1) = cmplx(0.4_dp,0.05_dp,dp)
    diag(2,2) = cmplx(-0.2_dp,-0.1_dp,dp)
    diag(3,3) = cmplx(1.05_dp,0.01_dp,dp)  ! unreachable, potentially observable
    diag(4,4) = cmplx(0.7_dp,-0.04_dp,dp) ! reachable, unobservable
    base_b = (0.0_dp,0.0_dp)
    base_b(1,:) = [cmplx(1.0_dp,0.2_dp,dp),cmplx(0.1_dp,0.0_dp,dp)]
    base_b(2,:) = [cmplx(0.3_dp,-0.1_dp,dp),cmplx(0.8_dp,0.1_dp,dp)]
    base_b(4,:) = [cmplx(0.3_dp,0.0_dp,dp),cmplx(0.4_dp,0.0_dp,dp)]
    base_c = (0.0_dp,0.0_dp)
    base_c(:,1) = [cmplx(1.0_dp,0.0_dp,dp),cmplx(0.5_dp,-0.1_dp,dp)]
    base_c(:,2) = [cmplx(-0.2_dp,0.1_dp,dp),cmplx(0.7_dp,-0.2_dp,dp)]
    base_c(:,3) = [cmplx(0.1_dp,0.0_dp,dp),cmplx(0.3_dp,0.2_dp,dp)]
    a = matmul(matmul(unitary,diag),transpose(conjg(unitary)))
    b = matmul(unitary,base_b)
    c = matmul(base_c,transpose(conjg(unitary)))
    d = reshape([cmplx(0.2_dp,0.1_dp,dp),cmplx(0.1_dp,0.0_dp,dp), &
        cmplx(0.0_dp,0.3_dp,dp),cmplx(-0.1_dp,0.0_dp,dp)],[2,2])
    original = dlti_ss(a,b,c,d,dt=0.125_dp)
    call check(original%status == dlti_success,'original')
    reduced = original%minreal()
    call check(reduced%status == dlti_success,'minreal status')
    call check(reduced%nstates() == 2,'removed unreachable and unobservable modes')
    call check(reduced%ninputs() == 2 .and. reduced%noutputs() == 2,'preserved MIMO shape')
    call check(maxval(abs(reduced%d-d)) < 1.0e-15_dp,'preserved feedthrough')
    do i=1,size(w)
        w(i) = 0.06_dp*i
    end do
    fr_a = original%freqresp_mimo(w=w,solver='direct')
    fr_b = reduced%freqresp_mimo(w=w,solver='direct')
    call check(fr_a%status == 0 .and. fr_b%status == 0,'frequency statuses')
    call check(maxval(abs(fr_a%response-fr_b%response)) < 3.0e-12_dp, &
        'preserve MIMO transfer matrix')
    do k=1,size(u,1)
        u(k,:) = [cmplx(cos(0.12_dp*k),sin(0.1_dp*k),dp), &
            cmplx(sin(0.18_dp*k),cos(0.2_dp*k),dp)]
    end do
    simulation_a = original%output(u)
    simulation_b = reduced%output(u)
    call check(maxval(abs(simulation_a%y-simulation_b%y)) < 1.0e-11_dp, &
        'preserve MIMO zero-state response')
    gram = original%gramians()
    call check(gram%status == dlti_numerical_failure,'original contains hidden unstable mode')
    gram = reduced%gramians()
    call check(gram%status == dlti_success,'reduced Gramians converge')
    bad = original%minreal(tol=-1.0_dp)
    call check(bad%status == dlti_invalid_input,'invalid negative rank tolerance')
    bad = original%minreal(tol=1.0_dp)
    call check(bad%status == dlti_invalid_input,'invalid large rank tolerance')
    zero = dlti_ss(reshape([cmplx(0.0_dp,0.0_dp,dp)],[0,0]), &
        reshape([cmplx(0.0_dp,0.0_dp,dp)],[0,1]), &
        reshape([cmplx(0.0_dp,0.0_dp,dp)],[1,0]), &
        reshape([cmplx(0.4_dp,0.2_dp,dp)],[1,1]))
    reduced = zero%minreal()
    call check(reduced%status == dlti_success .and. reduced%nstates() == 0,'zero-state minreal')
    print '(a)', 'PASS: numerical minimal realization, MIMO response and hidden modes'
contains
    subroutine check(ok,label)
        logical, intent(in) :: ok !! assertion
        character(*), intent(in) :: label !! failure context
        if (.not. ok) then
            print '(a,1x,a)', 'FAIL:',label
            error stop 1
        end if
    end subroutine check
end program test_signal_dlti_minreal
