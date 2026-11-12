! SPDX-License-Identifier: MIT
! Whole-record and streaming multitaper detrending, including precedence and
! synchronization of multichannel streamed spectra.
program test_signal_multitaper_whole_stream_detrend
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : spectral_result, cross_spectral_result, &
        spectral_status_success, spectral_status_invalid_input
    use scifort_signal_multitaper, only : multitaper_psd, multitaper_psd_ci, &
        multitaper_spectral_matrix, multitaper_csd, &
        multitaper_coherence, multitaper_harmonic_ftest, multitaper_harmonic_scan, &
        multitaper_stream, multitaper_spectral_matrix_spectrogram, &
        multitaper_matrix_result, multitaper_confidence_result, multitaper_harmonic_result, &
        multitaper_harmonic_scan_result, &
        multitaper_spectral_matrix_spectrogram_result
    implicit none
    integer, parameter :: n = 116, nc = 3, seg = 28
    real(dp) :: data(n,nc), linear(n), noise(n,nc), x(n), y(n), pi
    complex(dp) :: z(n,nc), zz(n)
    type(multitaper_stream) :: stream
    type(multitaper_matrix_result) :: matrix, transposed, bad_matrix, avg
    type(multitaper_spectral_matrix_spectrogram_result) :: batch, part
    type(multitaper_confidence_result) :: ci
    type(spectral_result) :: p, q, coherence, harm
    type(cross_spectral_result) :: csd
    complex(dp), allocatable :: combined(:,:,:,:)
    integer :: i, c, j, status, first, last, got, m, nf, iw, it
    integer, parameter :: chunks(12) = [1,0,8,19,3,7,11,1,18,0,9,39]
    character(8), parameter :: weightings(3) = [character(8) :: 'equal','eigen','adaptive']
    real(dp), parameter :: tol = 3.0e-10_dp

    pi = acos(-1.0_dp)
    do i = 1,n
        linear(i) = 15.0_dp + 0.42_dp*real(i,dp)
        do c = 1,nc
            noise(i,c) = cos(0.13_dp*real(i*c,dp)) + 0.31_dp*sin(0.37_dp*real(i+c,dp))
            data(i,c) = noise(i,c)+linear(i)*real(c,dp)
            z(i,c) = cmplx(data(i,c),0.6_dp*linear(i)+0.5_dp*noise(i,c),dp)
        end do
    end do
    x=data(:,1)
    y=data(:,2)
    zz=z(:,1)

    ! An exactly affine record must be eliminated, even in complex arithmetic.
    p = multitaper_psd(linear,3.0_dp,kmax=4,detrend_mode='linear')
    if (p%status /= spectral_status_success .or. maxval(p%power) > 1.0e-20_dp) &
        error stop 'real whole-record affine rejection'
    p = multitaper_psd(cmplx(linear,-1.4_dp*linear+3.0_dp,dp),3.0_dp, &
        kmax=4,detrend_mode='linear')
    if (p%status /= spectral_status_success .or. maxval(p%power) > 1.0e-19_dp) &
        error stop 'complex whole-record affine rejection'

    do iw=1,3
        p=multitaper_psd(x,3.0_dp,kmax=4,nfft=132,weighting=trim(weightings(iw)), &
            detrend=.false.,detrend_mode='linear')
        q=multitaper_psd(x,3.0_dp,kmax=4,nfft=132,weighting=trim(weightings(iw)), &
            detrend_mode='linear')
        if (p%status /= spectral_status_success .or. q%status /= spectral_status_success) &
            error stop 'whole-record PSD status'
        if (maxval(abs(p%power-q%power)) > tol) error stop 'detrend_mode precedence'
        p=multitaper_psd(zz,3.0_dp,kmax=4,nfft=132,weighting=trim(weightings(iw)), &
            detrend_mode='linear')
        if (p%status /= spectral_status_success) error stop 'complex PSD status'
        matrix=multitaper_spectral_matrix(data,3.0_dp,kmax=4,nfft=132, &
            weighting=trim(weightings(iw)),detrend_mode='linear')
        if (matrix%status /= spectral_status_success) error stop 'real matrix status'
        transposed=multitaper_spectral_matrix(transpose(data),3.0_dp,kmax=4,nfft=132, &
            weighting=trim(weightings(iw)),detrend_mode='linear',axis=2)
        if (maxval(abs(matrix%power-transposed%power)) > tol) error stop 'whole-record axis=2'
        q=multitaper_psd(x,3.0_dp,kmax=4,nfft=132,weighting=trim(weightings(iw)), &
            detrend_mode='linear')
        if (maxval(abs(real(matrix%power(:,1,1),dp)-q%power)) > tol) &
            error stop 'matrix/PSD parity'
        if (maxval(abs(matrix%power(:,1,2)-conjg(matrix%power(:,2,1)))) > tol) &
            error stop 'Hermitian whole record matrix'
        csd=multitaper_csd(x,y,3.0_dp,kmax=4,nfft=132,weighting=trim(weightings(iw)), &
            detrend_mode='linear')
        if (csd%status /= spectral_status_success .or. &
            maxval(abs(csd%power-matrix%power(:,1,2))) > tol) error stop 'whole-record CSD'
        coherence=multitaper_coherence(x,y,3.0_dp,kmax=4,nfft=132, &
            weighting=trim(weightings(iw)),detrend_mode='linear')
        if (coherence%status /= spectral_status_success) error stop 'whole-record coherence'
        if (any(coherence%power < 0.0_dp) .or. any(coherence%power > 1.0_dp)) &
            error stop 'coherence bounds'
        ci=multitaper_psd_ci(x,3.0_dp,kmax=4,nfft=132, &
            weighting=trim(weightings(iw)),detrend_mode='linear')
        if (ci%status /= spectral_status_success) error stop 'CI status'
        if (maxval(abs(ci%power-q%power)) > tol) error stop 'CI/PSD agreement'
    end do

    p=multitaper_psd(x,3.0_dp,detrend_mode='constant')
    q=multitaper_psd(x,3.0_dp,detrend=.true.)
    if (maxval(abs(p%power-q%power)) > tol) error stop 'default PSD regression'
    p=multitaper_psd(x,3.0_dp,detrend_mode='none')
    q=multitaper_psd(x,3.0_dp,detrend=.false.)
    if (maxval(abs(p%power-q%power)) > tol) error stop 'none PSD regression'
    p=multitaper_psd(x,3.0_dp,detrend_mode='quadratic')
    if (p%status /= spectral_status_invalid_input .or. size(p%power) /= 0) &
        error stop 'bad PSD detrend mode'
    bad_matrix=multitaper_spectral_matrix(data,3.0_dp,detrend_mode='bogus')
    if (bad_matrix%status /= spectral_status_invalid_input .or. size(bad_matrix%power) /= 0) &
        error stop 'bad matrix detrend mode'
    ci=multitaper_psd_ci(x,3.0_dp,kmax=4,detrend_mode='unknown')
    if (ci%status /= spectral_status_invalid_input) error stop 'bad CI detrend mode'
    harm=multitaper_psd(x,3.0_dp,detrend_mode='linear')
    if (harm%status /= spectral_status_success) error stop 'harmonic prerequisite'
    block
        type(multitaper_harmonic_result) :: h1, h2
        type(multitaper_harmonic_scan_result) :: scan
        h1=multitaper_harmonic_ftest(linear,3.0_dp,kmax=5,detrend_mode='linear')
        if (h1%status /= spectral_status_success) error stop 'harmonic linear status'
        h2=multitaper_harmonic_ftest(x,3.0_dp,kmax=5,detrend_mode='linear')
        if (h2%status /= spectral_status_success) error stop 'harmonic F status'
        scan=multitaper_harmonic_scan(x,3.0_dp,kmax=5,detrend_mode='linear')
        if (scan%status /= spectral_status_success) error stop 'harmonic scan status'
        if (maxval(abs(scan%f_statistic-h2%f_statistic)) > tol) &
            error stop 'harmonic F scan parity'
    end block

    ! Compare streamed matrices for irregular chunks against batch frames.
    do iw=1,3
        batch=multitaper_spectral_matrix_spectrogram(data,2.5_dp,seg,noverlap=17, &
            kmax=4,nfft=33,weighting=trim(weightings(iw)),detrend_mode='linear')
        call stream%initialize(2.5_dp,seg,nc,status,noverlap=17,kmax=4,nfft=33, &
            weighting=trim(weightings(iw)),detrend=.false.,detrend_mode='linear')
        if (status /= spectral_status_success) error stop 'stream init'
        m=size(batch%times)
        allocate(combined(size(batch%power,1),m,nc,nc))
        got=0
        first=1
        do j=1,size(chunks)
            if (first > n) exit
            last=min(n,first+chunks(j)-1)
            part=stream%push(data(first:last,:))
            if (part%status /= spectral_status_success) error stop 'stream push'
            nf=size(part%times)
            if (nf>0) combined(:,got+1:got+nf,:,:)=part%power
            got=got+nf
            first=last+1
        end do
        if (first <= n) then
            part=stream%push(data(first:n,:))
            if (part%status /= spectral_status_success) error stop 'final push'
            nf=size(part%times)
            if (nf>0) combined(:,got+1:got+nf,:,:)=part%power
            got=got+nf
        end if
        if (got/=m .or. stream%window_count()/=int(m,int64) .or. &
            stream%sample_count()/=int(n,int64)) error stop 'stream count'
        if (maxval(abs(combined-batch%power)) > tol) error stop 'stream/batch linear mismatch'
        avg=stream%mean_matrix()
        if (maxval(abs(avg%power-sum(batch%power,dim=2)/real(m,dp))) > tol) &
            error stop 'stream mean mismatch'
        call stream%reset()
        if (stream%sample_count()/=0_int64 .or. stream%window_count()/=0_int64) &
            error stop 'stream reset counters'
        part=stream%push(data)
        if (maxval(abs(part%power-batch%power)) > tol) error stop 'stream reset trend lost'
        deallocate(combined)
    end do

    batch=multitaper_spectral_matrix_spectrogram(z,2.5_dp,seg,noverlap=13,kmax=4, &
        nfft=32,weighting='adaptive',detrend_mode='linear')
    call stream%initialize(2.5_dp,seg,nc,status,noverlap=13,kmax=4, &
        nfft=32,weighting='adaptive',complex_input=.true.,detrend_mode='linear')
    if (status /= spectral_status_success) error stop 'complex stream init'
    part=stream%push(transpose(z),axis=2)
    if (part%status /= spectral_status_success .or. &
        maxval(abs(part%power-batch%power)) > tol) error stop 'complex stream linear axis=2'
    i=int(stream%sample_count())
    part=stream%push(transpose(z(1:4,:)),axis=3)
    if (part%status /= spectral_status_invalid_input .or. &
        stream%sample_count()/=int(i,int64)) error stop 'invalid push changed state'
    call stream%initialize(2.5_dp,seg,nc,status,detrend_mode='quadratic')
    if (status /= spectral_status_invalid_input) error stop 'bad stream detrend mode'

    print '(a)', 'Whole-record and streaming multitaper detrending passed'
end program test_signal_multitaper_whole_stream_detrend
