! SPDX-License-Identifier: MIT
! Regression tests for missing scipy.signal normalization, smoother and FFT dispatch.
program test_signal_whittaker_normalize_conv
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_whittaker, only : whittaker_henderson, whittaker_result, &
        whittaker_status_success, whittaker_status_invalid_input, whittaker_status_singular
    use scifort_signal_tf_normalize, only : normalize, normalized_tf_result, &
        normalize_status_success, normalize_status_invalid_input
    use scifort_signal_abcd_normalize, only : abcd_normalize
    use scifort_signal_state_space, only : state_space_result, state_space_success, state_space_invalid_input
    use scifort_signal_convolve, only : convolve, correlate, choose_conv_method, &
        convolution_status_success, convolution_status_invalid_input
    implicit none
    type(whittaker_result) :: wh
    type(normalized_tf_result) :: tf
    type(state_space_result) :: sys
    real(dp), allocatable :: v(:), yy(:)
    complex(dp), allocatable :: zz(:)
    real(dp) :: t(2), nan
    integer :: ierr, i
    character(len=6) :: method

    wh = whittaker_henderson([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp],150.0_dp)
    if (wh%status/=whittaker_status_success .or. maxval(abs(wh%x-[1,2,3,4,5]))>1.e-11_dp) &
        error stop 'WH linear polynomial exactness'
    wh = whittaker_henderson([1.0_dp,2.0_dp,0.0_dp,4.0_dp,5.0_dp],10.0_dp, &
        weights=[1.0_dp,1.0_dp,0.0_dp,1.0_dp,1.0_dp])
    if (wh%status/=whittaker_status_success .or. maxval(abs(wh%x-[1,2,3,4,5]))>1.e-11_dp) &
        error stop 'WH interpolation'
    wh = whittaker_henderson([2.0_dp,-4.0_dp,1.0_dp],0.0_dp,order=1)
    if (wh%status/=whittaker_status_success .or. any(wh%x/=[2.0_dp,-4.0_dp,1.0_dp])) &
        error stop 'WH lambda zero'
    nan = ieee_value(0.0_dp,ieee_quiet_nan)
    wh = whittaker_henderson([1.0_dp,nan,3.0_dp,4.0_dp],1.0_dp, &
        weights=[1.0_dp,0.0_dp,1.0_dp,1.0_dp])
    if (wh%status/=whittaker_status_success) error stop 'WH masked NaN'
    wh = whittaker_henderson([1.0_dp,nan,3.0_dp],1.0_dp)
    if (wh%status/=whittaker_status_invalid_input) error stop 'WH unmasked NaN'
    wh = whittaker_henderson([1.0_dp,2.0_dp,3.0_dp],-1.0_dp)
    if (wh%status/=whittaker_status_invalid_input) error stop 'WH invalid penalty'
    wh = whittaker_henderson([1.0_dp,2.0_dp,3.0_dp],1.0_dp,weights=[0.0_dp,0.0_dp,0.0_dp])
    if (wh%status/=whittaker_status_singular) error stop 'WH singular system'

    tf = normalize([0.0_dp,3.0_dp,6.0_dp],[0.0_dp,2.0_dp,-5.0_dp,4.0_dp])
    if (tf%status/=normalize_status_success .or. tf%trimmed/=1) error stop 'TF numerator trimming'
    if (maxval(abs(tf%b(1,:)-cmplx([1.5_dp,3.0_dp],0.0_dp,dp)))>1.e-14_dp) error stop 'TF numerator'
    if (maxval(abs(tf%a-cmplx([1.0_dp,-2.5_dp,2.0_dp],0.0_dp,dp)))>1.e-14_dp) error stop 'TF denominator'
    tf = normalize(reshape([0.0_dp,0.0_dp,1.0_dp,3.0_dp,2.0_dp,4.0_dp],[2,3]),[2.0_dp,1.0_dp])
    if (tf%status/=normalize_status_success .or. tf%trimmed/=1) error stop 'TF multi-output trimming'
    if (any(shape(tf%b)/=[2,2])) error stop 'TF multi-output shape'
    tf = normalize([1.0_dp,2.0_dp],cmplx([1.0_dp,2.0_dp],[1.0_dp,-1.0_dp],dp))
    if (tf%status/=normalize_status_success) error stop 'TF real numerator, complex denominator'
    tf = normalize(cmplx([1.0_dp,2.0_dp],[1.0_dp,2.0_dp],dp),[2.0_dp,1.0_dp])
    if (tf%status/=normalize_status_success) error stop 'TF complex numerator, real denominator'
    if (maxval(abs(tf%b(1,:) - [cmplx(.5_dp,.5_dp,dp), cmplx(1.0_dp,1.0_dp,dp)])) &
        >1.e-14_dp) error stop 'TF mixed numerical coefficients'
    tf = normalize([0.0_dp,0.0_dp],[0.0_dp,0.0_dp])
    if (tf%status/=normalize_status_invalid_input) error stop 'TF invalid denominator'

    sys = abcd_normalize(a=reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp],[2,2]), &
        b=reshape([-1.0_dp,5.0_dp],[2,1]),d=reshape([2.5_dp],[1,1]))
    if (sys%status/=state_space_success) error stop 'ABCD normalize missing C'
    if (any(shape(sys%a)/=[2,2]) .or. any(shape(sys%b)/=[2,1]) .or. &
        any(shape(sys%c)/=[1,2]) .or. any(shape(sys%d)/=[1,1])) error stop 'ABCD shapes'
    if (any(abs(sys%c)>1.e-15_dp)) error stop 'ABCD zeros C'
    sys = abcd_normalize(a=reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp],[2,2]), &
        b=reshape([-1.0_dp,5.0_dp],[2,1]),d=reshape([2.5_dp,6.0_dp],[2,1]))
    if (sys%status/=state_space_success .or. any(shape(sys%c)/=[2,2])) error stop 'ABCD two outputs'
    sys = abcd_normalize(a=reshape([1.0_dp],[1,1]),b=reshape([1.0_dp,2.0_dp],[2,1]), &
        c=reshape([1.0_dp],[1,1]))
    if (sys%status/=state_space_invalid_input) error stop 'ABCD invalid shape'
    sys = abcd_normalize(b=reshape([1,2],[2,1]),c=reshape([3,4],[1,2]))
    if (sys%status/=state_space_success .or. any(shape(sys%d)/=[1,1])) &
        error stop 'ABCD mixed integer matrices'
    if (abs(sys%b(2,1) - cmplx(2.0_dp,0.0_dp,dp))>1.e-14_dp) &
        error stop 'ABCD integer conversion'
    sys = abcd_normalize(a=reshape([real(dp) :: ],[0,0]), &
        d=reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp],[2,3]))
    if (sys%status/=state_space_success .or. any(shape(sys%b)/=[0,3]) .or. &
        any(shape(sys%c)/=[2,0])) error stop 'ABCD zero-state system'
    sys = abcd_normalize()
    if (sys%status/=state_space_invalid_input) error stop 'ABCD missing all'

    v = convolve([1.0_dp,2.0_dp,3.0_dp],[4.0_dp,5.0_dp],method='fft',status=ierr)
    if (ierr/=convolution_status_success .or. any(abs(v-[4,13,22,15])>1.e-12_dp)) &
        error stop 'FFT convolution'
    zz = correlate(cmplx([1.0_dp,2.0_dp],[1.0_dp,-1.0_dp],dp), &
        cmplx([1.0_dp,3.0_dp],[2.0_dp,1.0_dp],dp),method='fft',status=ierr)
    if (ierr/=convolution_status_success) error stop 'FFT correlation status'
    if (maxval(abs(zz - correlate(cmplx([1.0_dp,2.0_dp],[1.0_dp,-1.0_dp],dp), &
        cmplx([1.0_dp,3.0_dp],[2.0_dp,1.0_dp],dp),method='direct')))>1.e-12_dp) &
        error stop 'FFT complex correlation'
    method = choose_conv_method([1.0_dp,2.0_dp],[3.0_dp,4.0_dp],status=ierr)
    if (method/='direct' .or. ierr/=convolution_status_success) error stop 'small choose method'
    allocate(yy(4096))
    yy = [(sin(real(i,dp)*0.1_dp),i=1,4096)]
    method = choose_conv_method(yy,yy)
    if (method/='fft') error stop 'large choose method'
    method = choose_conv_method(yy(1:128),yy(1:100),measure=.true.,times=t)
    if (.not. (method=='direct' .or. method=='fft') .or. any(t<0.0_dp)) error stop 'measured method'
    v = convolve(yy(1:512),yy(1:400),method='auto')
    if (maxval(abs(v-convolve(yy(1:512),yy(1:400),method='direct')))>1.e-10_dp) &
        error stop 'auto convolution'
    v = convolve([1.0_dp,2.0_dp],[3.0_dp],method='unsupported',status=ierr)
    if (ierr/=convolution_status_invalid_input .or. size(v)/=0) error stop 'invalid method'

    print '(a)', 'signal WH, TF/ABCD normalization and convolution-method tests passed'
end program test_signal_whittaker_normalize_conv
