! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_missing_algorithms
    use scifort_kinds, only : dp
    use scifort_signal_stft_dual, only : closest_stft_dual_window, &
        stft_dual_status_success, stft_dual_status_noninvertible
    use scifort_signal_czt, only : czt_points, czt_status_success, czt_status_invalid_input
    use scifort_signal_analog, only : findfreqs
    use scifort_signal_symiir, only : symiirorder1, symiirorder2, &
        symiir_status_success, symiir_status_invalid_input, symiir_status_no_convergence
    use scifort_signal_fir2d, only : firwin_2d
    use scifort_signal_unique_roots, only : unique_roots, &
        unique_roots_status_success, unique_roots_status_invalid_input
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), allocatable :: dual(:), freqs(:), rows(:,:), ra(:), smoothed(:), smoothed2(:,:)
    complex(dp), allocatable :: cdual(:), cz(:), cp(:), cr(:), csmoothed(:)
    complex(dp) :: alpha_complex, input_complex(6)
    integer, allocatable :: mult(:)
    real(dp) :: alpha_real, data(8), samples(120), matrix(3,120)
    integer :: i, k, st

    data = 1.0_dp
    call closest_stft_dual_window(data,2,dual,alpha_real,status=st)
    if (st /= stft_dual_status_success) error stop 'constant dual status'
    if (maxval(abs(dual-0.25_dp)) > 1e-13_dp) error stop 'constant dual samples'
    if (abs(alpha_real-0.25_dp) > 1e-13_dp) error stop 'constant dual scale'
    call closest_stft_dual_window(data,2,dual,alpha_real,scaled=.false.,status=st)
    if (st /= 0 .or. abs(alpha_real-1.0_dp) > 1e-14_dp) error stop 'unscaled dual status'
    if (maxval(abs(dual-0.25_dp)) > 1e-13_dp) error stop 'unscaled dual values'
    data = 0.0_dp
    data(1::2) = 1.0_dp
    call closest_stft_dual_window(data,2,dual,status=st)
    if (st /= stft_dual_status_noninvertible .or. size(dual) /= 0) error stop 'singular dual'
    input_complex = [(cmplx(real(i,dp),0.125_dp*real(i-2,dp),dp),i=1,6)]
    call closest_stft_dual_window(input_complex,2,cdual,alpha_complex,scaled=.false.,status=st)
    if (st /= 0 .or. abs(alpha_complex-1.0_dp) > 1.0e-14_dp) error stop 'complex dual status'
    do k = 1,2
        if (abs(sum(conjg(input_complex(k::2))*cdual(k::2))-1.0_dp) > 1.0e-12_dp) &
            error stop 'complex perfect reconstruction'
    end do

    cz = czt_points(8,status=st)
    if (st /= czt_status_success .or. size(cz) /= 8) error stop 'CZT points status'
    do i = 1,8
        if (abs(cz(i)-exp(cmplx(0.0_dp,2.0_dp*pi*real(i-1,dp)/8.0_dp,dp))) > 1.0e-13_dp) &
            error stop 'CZT unit circle'
    end do
    cz = czt_points(5,w=cmplx(0.8_dp,-0.2_dp,dp),a=cmplx(1.3_dp,0.4_dp,dp))
    do i = 1,5
        if (abs(cz(i)-cmplx(1.3_dp,0.4_dp,dp)/cmplx(0.8_dp,-0.2_dp,dp)**(i-1)) > 1.0e-12_dp) &
            error stop 'CZT logarithmic spiral'
    end do
    cz = czt_points(0,status=st)
    if (st /= czt_status_invalid_input .or. size(cz) /= 0) error stop 'invalid CZT points'

    freqs = findfreqs([1.0_dp,0.0_dp],[1.0_dp,8.0_dp,25.0_dp],9,status=st)
    if (st /= 0 .or. size(freqs) /= 9) error stop 'findfreqs real status'
    do i = 1,9
        if (abs(log10(freqs(i))-(-2.0_dp+0.5_dp*real(i-1,dp))) > 1.0e-13_dp) &
            error stop 'findfreqs real values'
    end do
    freqs = findfreqs([cmplx(0.0_dp,0.0_dp,dp)], &
        [cmplx(-4.0_dp,3.0_dp,dp),cmplx(-4.0_dp,-3.0_dp,dp)],10,kind='zp',status=st)
    if (st /= 0 .or. any(freqs <= 0.0_dp)) error stop 'findfreqs complex zpk'
    freqs = findfreqs([1.0_dp],[1.0_dp],0,status=st)
    if (st == 0 .or. size(freqs) /= 0) error stop 'invalid findfreqs count'

    call unique_roots([1.0_dp,1.009_dp,1.018_dp,5.0_dp],ra,mult, &
        tol=0.01_dp,rtype='avg',status=st)
    if (st /= unique_roots_status_success .or. size(ra) /= 3) error stop 'grouping status'
    if (any(mult /= [2,1,1])) error stop 'greedy root grouping count'
    if (abs(ra(1)-1.0045_dp) > 1.0e-12_dp) error stop 'average root'
    cp = [cmplx(1.0_dp,0.2_dp,dp),cmplx(1.0_dp,0.21_dp,dp), &
          cmplx(2.0_dp,0.0_dp,dp)]
    call unique_roots(cp,cr,mult,tol=0.02_dp,rtype='max',status=st)
    if (st /= 0 .or. any(mult /= [2,1])) error stop 'complex root grouping'
    if (abs(cr(1)-cp(2)) > 1.0e-14_dp) error stop 'complex lexicographic max'
    call unique_roots(cp,cr,mult,rtype='bad',status=st)
    if (st /= unique_roots_status_invalid_input .or. size(cr) /= 0) error stop 'invalid grouping'

    rows = firwin_2d([5,7],'hamming',0.1_dp,status=st)
    if (st /= 0 .or. any(shape(rows) /= [5,7])) error stop 'separable FIR shape'
    if (abs(sum(rows)-1.0_dp) > 1.0e-13_dp) error stop 'separable FIR gain'
    if (maxval(abs(rows-rows(5:1:-1,7:1:-1))) > 1.0e-13_dp) error stop 'separable FIR symmetry'
    rows = firwin_2d([5,7],'hamming',0.1_dp,circular=.true.,status=st)
    if (st /= 0 .or. any(shape(rows) /= [7,5])) error stop 'circular FIR shape'
    if (maxval(abs(rows-rows(7:1:-1,5:1:-1))) > 1.0e-13_dp) error stop 'circular FIR symmetry'
    rows = firwin_2d([0,4],'hamming',0.1_dp,status=st)
    if (st == 0 .or. size(rows) /= 0) error stop 'invalid FIR shape'
    samples = [(sin(0.1_dp*real(i,dp))+0.2_dp*cos(0.37_dp*real(i,dp)),i=1,120)]
    smoothed = symiirorder1(samples,0.25_dp,0.5_dp,precision=1.0e-12_dp,status=st)
    if (st /= symiir_status_success) error stop 'symiirorder1 status'
    csmoothed = symiirorder1(cmplx(samples,0.25_dp*samples,dp), &
        0.25_dp,0.5_dp,precision=1.0e-12_dp,status=st)
    if (st /= 0 .or. maxval(abs(csmoothed-cmplx(smoothed,0.25_dp*smoothed,dp))) > 1.0e-13_dp) &
        error stop 'symiirorder1 complex linearity'
    smoothed = symiirorder2(samples,0.5_dp,0.8_dp,precision=1.0e-12_dp,status=st)
    if (st /= symiir_status_success .or. size(smoothed) /= 120) error stop 'symiirorder2 status'
    do i = 1,3
        matrix(i,:) = samples*real(i,dp)
    end do
    smoothed2 = symiirorder2(matrix,0.5_dp,0.8_dp,precision=1.0e-12_dp,status=st)
    if (st /= 0 .or. any(shape(smoothed2) /= [3,120])) error stop 'symiirorder2 batch status'
    if (maxval(abs(smoothed2(1,:)-smoothed)) > 1.0e-12_dp) error stop 'symiirorder2 batch values'
    smoothed = symiirorder1(samples,1.0_dp,1.0_dp,status=st)
    if (st /= symiir_status_invalid_input .or. size(smoothed) /= 0) error stop 'invalid symmetric pole'
    smoothed = symiirorder2(samples,0.999_dp,0.6_dp,status=st)
    if (st /= symiir_status_no_convergence .or. size(smoothed) /= 0) error stop 'short symmetric tail'
    print '(a)', 'PASS: seven signal algorithms (dual STFT, FIR2D, symmetric IIR, roots and CZT)' 
end program test_signal_missing_algorithms
