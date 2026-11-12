! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_missing_algorithms_public
    use scifort_signal, only : dp, czt_points, findfreqs, unique_roots, &
        firwin_2d, closest_stft_dual_window, symiirorder1, symiirorder2, &
        czt_status_success, design_status_success, stft_dual_status_success, &
        unique_roots_status_success, symiir_status_success
    implicit none
    real(dp), allocatable :: frequencies(:), grouped(:), dual(:), fir(:,:), smooth(:)
    complex(dp), allocatable :: z(:)
    integer, allocatable :: counts(:)
    real(dp) :: data(100)
    integer :: st, i

    z = czt_points(4,status=st)
    if (st /= czt_status_success .or. abs(z(3)+1.0_dp) > 1.0e-12_dp) error stop 'public CZT points'
    frequencies = findfreqs([1.0_dp],[1.0_dp,1.0_dp],10,status=st)
    if (st /= design_status_success .or. size(frequencies) /= 10) error stop 'public findfreqs'
    call unique_roots([1.0_dp,1.001_dp,3.0_dp],grouped,counts,tol=0.01_dp,status=st)
    if (st /= unique_roots_status_success .or. any(counts /= [2,1])) error stop 'public unique roots'
    fir = firwin_2d([7,9],'hamming',0.2_dp,status=st)
    if (st /= design_status_success .or. any(shape(fir) /= [7,9])) error stop 'public FIR2D'
    call closest_stft_dual_window([1.0_dp,1.0_dp,1.0_dp,1.0_dp],2,dual,status=st)
    if (st /= stft_dual_status_success .or. maxval(abs(dual-0.5_dp)) > 1.0e-12_dp) &
        error stop 'public dual STFT'
    data = [(sin(0.12_dp*real(i,dp)),i=1,100)]
    smooth = symiirorder1(data,0.25_dp,0.5_dp,status=st)
    if (st /= symiir_status_success .or. size(smooth) /= 100) error stop 'public IIR order1'
    smooth = symiirorder2(data,0.5_dp,0.7_dp,status=st)
    if (st /= symiir_status_success .or. size(smooth) /= 100) error stop 'public IIR order2'
    print '(a)', 'PASS: seven new signal algorithms via public scifort_signal'
end program test_signal_missing_algorithms_public
