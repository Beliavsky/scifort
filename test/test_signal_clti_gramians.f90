! SPDX-License-Identifier: MIT
! Continuous-time Lyapunov Gramians and balanced state-space reduction.
program test_signal_clti_gramians
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_signal, only : dp, clti_system, clti_ss, clti_gramian_result, &
        clti_balanced_result, clti_success, clti_invalid_input, clti_numerical_failure, &
        mimo_frequency_response
    implicit none
    type(clti_system) :: sys, unstable, static, dissipative
    type(clti_gramian_result) :: gr, gr_d, gr_bad
    type(clti_balanced_result) :: reduced, complete, empty, bad
    type(mimo_frequency_response) :: response, smaller
    complex(dp) :: a(3,3),b(3,3),c(3,3),d(3,3),mix(3,3)
    complex(dp) :: am(3,3),bm(3,3),cm(3,3)
    complex(dp), allocatable :: ctl(:,:),obs(:,:)
    real(dp) :: hsv_true(3),h2_true,w(40),pi,err
    integer :: i,j

    pi = acos(-1.0_dp)
    a = (0.0_dp,0.0_dp)
    b = (0.0_dp,0.0_dp)
    c = (0.0_dp,0.0_dp)
    d = (0.0_dp,0.0_dp)
    a(1,1) = cmplx(-0.6_dp,0.8_dp,dp)
    a(2,2) = cmplx(-1.5_dp,-0.5_dp,dp)
    a(3,3) = cmplx(-2.0_dp,0.0_dp,dp)
    b(1,1) = (1.2_dp,0.0_dp)
    b(2,2) = (1.0_dp,0.0_dp)
    b(3,3) = (0.5_dp,0.0_dp)
    c(1,1) = (0.6_dp,0.0_dp)
    c(2,2) = (1.3_dp,0.0_dp)
    c(3,3) = (0.2_dp,0.0_dp)
    hsv_true = [0.6_dp,1.3_dp/3.0_dp,0.025_dp]
    h2_true = sqrt(0.72_dp**2/1.2_dp+1.3_dp**2/3.0_dp+0.1_dp**2/4.0_dp)
    do j=1,3
        do i=1,3
            mix(i,j) = exp(cmplx(0.0_dp,2.0_dp*pi*real((i-1)*(j-1),dp)/3.0_dp,dp))/sqrt(3.0_dp)
        end do
    end do
    am = matmul(matmul(mix,a),transpose(conjg(mix)))
    bm = matmul(mix,b)
    cm = matmul(c,transpose(conjg(mix)))
    sys = clti_ss(am,bm,cm,d)
    call check(sys%status() == clti_success,'MIMO continuous system')
    ctl = sys%controllability()
    obs = sys%observability()
    call check(all(shape(ctl) == [3,9]),'controllability shape')
    call check(all(shape(obs) == [9,3]),'observability shape')
    gr = sys%gramians(tol=1.0e-13_dp)
    call check(gr%status == clti_success,'Cayley Gramians')
    call check(maxval(abs(matmul(am,gr%wc)+matmul(gr%wc,transpose(conjg(am)))+ &
        matmul(bm,transpose(conjg(bm))))) < 5.0e-11_dp,'controllability Lyapunov residual')
    call check(maxval(abs(matmul(transpose(conjg(am)),gr%wo)+matmul(gr%wo,am)+ &
        matmul(transpose(conjg(cm)),cm))) < 5.0e-11_dp,'observability Lyapunov residual')
    call check(abs(gr%h2_norm-h2_true) < 5.0e-12_dp,'continuous-time H2 norm')

    reduced = sys%balanced_truncate(2)
    call check(reduced%status == clti_success,'balanced truncation status')
    call check(reduced%model%nstates() == 2,'retained state count')
    call check(maxval(abs(reduced%hsv-hsv_true)) < 1.0e-10_dp,'Hankel singular values')
    call check(abs(reduced%error_bound-2.0_dp*hsv_true(3)) < 1.0e-10_dp,'balanced error bound')
    do i=1,size(w)
        w(i) = 5.0_dp*real(i-1,dp)/real(size(w)-1,dp)
    end do
    response = sys%freqresp_mimo(w=w,solver='direct')
    smaller = reduced%model%freqresp_mimo(w=w,solver='direct')
    call check(response%status == 0 .and. smaller%status == 0,'analog frequency responses')
    err = maxval(abs(response%response-smaller%response))
    call check(err <= reduced%error_bound+1.0e-10_dp,'analog H-infinity bound')
    complete = sys%balanced_truncate(3)
    call check(complete%status == clti_success,'full rank balanced transform')
    smaller = complete%model%freqresp_mimo(w=w,solver='direct')
    call check(maxval(abs(response%response-smaller%response)) < 2.0e-11_dp, &
        'full-order continuous transfer function')
    empty = sys%balanced_truncate(0)
    call check(empty%status == clti_success,'order zero reduction')
    call check(empty%model%nstates() == 0,'zero-state model')
    ! A nonzero continuous-time D gives infinite H2 energy (Dirac impulse).
    d(1,1) = (0.1_dp,-0.2_dp)
    dissipative = clti_ss(am,bm,cm,d)
    gr_d = dissipative%gramians()
    call check(gr_d%status == clti_success,'feedthrough does not alter Gramians')
    call check(.not. ieee_is_finite(gr_d%h2_norm),'non-strictly-proper H2 infinite')
    reduced = dissipative%balanced_truncate(2)
    call check(reduced%status == clti_success,'feedthrough balanced status')
    smaller = reduced%model%freqresp_mimo(w=w,solver='direct')
    call check(abs(smaller%response(1,1,1)-response%response(1,1,1)) < &
        reduced%error_bound+abs(d(1,1))+1.0e-10_dp,'feedthrough frequency test')
    bad = sys%balanced_truncate(-1)
    call check(bad%status == clti_invalid_input,'invalid negative order')
    bad = sys%balanced_truncate(4)
    call check(bad%status == clti_invalid_input,'invalid large order')
    bad = sys%balanced_truncate(1,tol=-0.1_dp)
    call check(bad%status == clti_invalid_input,'invalid rank tolerance')
    a(3,3) = cmplx(0.2_dp,0.0_dp,dp)
    unstable = clti_ss(a,b,c,d)
    gr_bad = unstable%gramians()
    call check(gr_bad%status == clti_numerical_failure,'unstable Gramian rejected')
    bad = unstable%balanced_truncate(1)
    call check(bad%status == clti_numerical_failure,'unstable balance rejected')
    static = clti_ss(am(:0,:0),bm(:0,:),cm(:,:0),d)
    gr_d = static%gramians()
    call check(gr_d%status == clti_success,'zero-state Gramians')
    call check(size(gr_d%wc) == 0 .and. size(gr_d%wo) == 0,'zero-state arrays')
    empty = static%balanced_truncate(0)
    call check(empty%status == clti_success,'zero-state balancing')
    print '(a)', 'PASS: continuous Gramians, H2 norm, balancing, invalid cases'
contains
    subroutine check(ok,label)
        logical, intent(in) :: ok !! tested condition
        character(*), intent(in) :: label !! test label
        if (.not. ok) then
            print '(a,1x,a)', 'FAIL:',label
            error stop 1
        end if
    end subroutine check
end program test_signal_clti_gramians
