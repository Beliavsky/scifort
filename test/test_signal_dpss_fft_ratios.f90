! SPDX-License-Identifier: MIT
! Compare exact Toeplitz concentration with zero-padded FFT autocorrelation.
program test_signal_dpss_fft_ratios
    use scifort_kinds, only : dp
    use scifort_signal_dpss, only : dpss, dpss_result, &
        dpss_status_success, dpss_status_invalid_input
    implicit none
    integer, parameter :: lengths(*) = [2,3,4,7,8,15,16,17,31,32,33, &
        63,64,65,127,128,129,255,256,257,511,512,513,1024,2048]
    type(dpss_result) :: direct, fourier, automatic
    integer :: i, m, k, ncases, isym, iband
    real(dp) :: nw, error, worst
    character(12) :: solver

    ncases = 0
    worst = 0.0_dp
    do i = 1, size(lengths)
        m = lengths(i)
        k = min(5,m)
        solver = 'ql'
        if (m >= 32) solver = 'partial'
        do isym = 0, 1
            do iband = 1, 3
                select case (iband)
                case (1)
                    nw = min(0.3_dp,0.25_dp*real(m,dp))
                case (2)
                    nw = min(3.5_dp,0.35_dp*real(m,dp))
                case default
                    nw = 0.44_dp*real(m,dp)
                end select
                direct = dpss(m,nw,kmax=k,sym=isym==1, &
                    solver=trim(solver),ratio_method='direct')
                fourier = dpss(m,nw,kmax=k,sym=isym==1, &
                    solver=trim(solver),ratio_method='fft')
                automatic = dpss(m,nw,kmax=k,sym=isym==1, &
                    solver=trim(solver))
                if (direct%status /= dpss_status_success) error stop 'direct DPSS failed'
                if (fourier%status /= dpss_status_success) error stop 'FFT DPSS failed'
                if (automatic%status /= dpss_status_success) error stop 'auto DPSS failed'
                if (maxval(abs(direct%tapers-fourier%tapers)) /= 0.0_dp) then
                    error stop 'FFT ratio option changed tapers'
                end if
                error = maxval(abs(direct%ratios-fourier%ratios))
                worst = max(worst,error)
                if (error > 4.0e-11_dp) then
                    print *, 'bad ratio',m,nw,isym,error
                    error stop 'FFT and direct concentration disagree'
                end if
                if (m >= 256 .or. (m == 255 .and. isym == 0)) then
                    if (any(automatic%ratios /= fourier%ratios)) then
                        error stop 'automatic method did not select FFT'
                    end if
                else
                    if (any(automatic%ratios /= direct%ratios)) then
                        error stop 'automatic method did not select direct sum'
                    end if
                end if
                if (any(fourier%ratios < 0.0_dp) .or. &
                    any(fourier%ratios > 1.0_dp)) error stop 'ratios out of range'
                ncases = ncases+1
            end do
        end do
    end do
    ! Explicit FFT handling of degenerate size and invalid options.
    fourier = dpss(1,0.1_dp,ratio_method='fft')
    if (fourier%status /= dpss_status_success) error stop 'length-one FFT ratio'
    if (abs(fourier%ratios(1)-1.0_dp) > epsilon(1.0_dp)) error stop 'length-one ratio'
    fourier = dpss(64,3.5_dp,kmax=3,ratio_method='unknown')
    if (fourier%status /= dpss_status_invalid_input) error stop 'bad method accepted'
    fourier = dpss(64,3.5_dp,kmax=3,ratio_method=' FFT')
    if (fourier%status /= dpss_status_invalid_input) error stop 'invalid whitespace accepted'
    fourier = dpss(0,0.5_dp,ratio_method='fft')
    if (fourier%status /= dpss_status_success) error stop 'empty DPSS changed'

    print '(a,i0,a,es12.4)', 'test_signal_dpss_fft_ratios passed: ', &
        ncases, ' method comparison cases; max difference=',worst
end program test_signal_dpss_fft_ratios
