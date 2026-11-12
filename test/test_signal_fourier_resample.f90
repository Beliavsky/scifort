! SPDX-License-Identifier: MIT
program test_signal_fourier_resample
    use scifort_signal, only : dp, resample, &
        fourier_resample_status_success, fourier_resample_status_invalid_input
    use scifort_signal_spectral, only : forward_fft
    implicit none
    real(dp), allocatable :: y(:), times(:)
    complex(dp), allocatable :: z(:), fx(:)
    real(dp), parameter :: tol=1.0e-12_dp
    integer :: st, k

    y=resample([1.0_dp,2.0_dp,3.0_dp,4.0_dp],2,status=st)
    if (st /= fourier_resample_status_success) error stop 'resample short output'
    call check(y,[1.5_dp,3.5_dp],tol)
    y=resample([1.0_dp,2.0_dp,3.0_dp,4.0_dp],7, &
        t=[1.0_dp,1.5_dp],output_times=times,status=st)
    call check(y,[1.0_dp,1.2059391826513937_dp,2.1980774557257_dp, &
        2.655340227855494_dp,3.52310770609061_dp, &
        4.147933280089347_dp,2.7696021475874533_dp],tol)
    if (size(times) /= 7) error stop 'resample timestamps length'
    if (abs(times(7)-(1.0_dp+6.0_dp*2.0_dp/7.0_dp)) > tol) &
        error stop 'resample timestamps value'
    z=resample([cmplx(1.0_dp,1.0_dp,dp),cmplx(2.0_dp,-1.0_dp,dp), &
        cmplx(3.0_dp,2.0_dp,dp)],6,status=st)
    if (st /= fourier_resample_status_success) error stop 'resample complex status'
    if (maxval(abs(z-[cmplx(1.0_dp,1.0_dp,dp), &
        cmplx(1.0_dp,-2.0_dp/3.0_dp,dp),cmplx(2.0_dp,-1.0_dp,dp), &
        cmplx(3.0_dp,1.0_dp/3.0_dp,dp),cmplx(3.0_dp,2.0_dp,dp), &
        cmplx(2.0_dp,7.0_dp/3.0_dp,dp)]))>tol) error stop 'complex Nyquist handling'
    fx=[cmplx(1.0_dp,1.0_dp,dp),cmplx(2.0_dp,-1.0_dp,dp),cmplx(3.0_dp,2.0_dp,dp)]
    call forward_fft(fx)
    z=resample(fx,6,domain='freq',status=st)
    if (st /= fourier_resample_status_success) error stop 'frequency-domain resample status'
    if (maxval(abs(z-[cmplx(1.0_dp,1.0_dp,dp), &
        cmplx(1.0_dp,-2.0_dp/3.0_dp,dp),cmplx(2.0_dp,-1.0_dp,dp), &
        cmplx(3.0_dp,1.0_dp/3.0_dp,dp),cmplx(3.0_dp,2.0_dp,dp), &
        cmplx(2.0_dp,7.0_dp/3.0_dp,dp)]))>tol) error stop 'frequency-domain values'
    y=resample([1.0_dp,2.0_dp],0,status=st)
    if (st /= fourier_resample_status_invalid_input .or. size(y) /= 0) error stop 'bad num'
    y=resample([1.0_dp,2.0_dp],4,window='unavailable',status=st)
    if (st /= fourier_resample_status_invalid_input .or. size(y) /= 0) error stop 'bad window'
    z=resample([cmplx(1.0_dp,0.0_dp,dp)],3,domain='unknown',status=st)
    if (st /= fourier_resample_status_invalid_input .or. size(z) /= 0) error stop 'bad domain'
    print '(a)', 'test_signal_fourier_resample passed'
contains
    subroutine check(a,b,atol)
        real(dp), intent(in) :: a(:) !! actual samples
        real(dp), intent(in) :: b(:) !! SciPy reference values
        real(dp), intent(in) :: atol !! absolute error bound
        if (size(a) /= size(b)) error stop 'resample shape mismatch'
        if (size(a)>0) then
            if (maxval(abs(a-b))>atol) error stop 'resample value mismatch'
        end if
    end subroutine check
end program test_signal_fourier_resample
