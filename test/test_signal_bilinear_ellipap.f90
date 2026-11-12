! SPDX-License-Identifier: MIT
! Native regression for bilinear and elliptic analog prototype functions.
program test_signal_bilinear_ellipap
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_bilinear, only : bilinear, bilinear_zpk
    use scifort_signal_ellipap, only : ellipap
    use scifort_signal_analog_prototypes, only : analog_prototype_result, &
        prototype_success, prototype_invalid_input, prototype_too_expensive, besselap
    use scifort_signal_frequency_transforms, only : analog_tf_transform_result, &
        analog_zpk_transform_result, frequency_transform_success, &
        frequency_transform_invalid_input, frequency_transform_too_expensive
    implicit none
    type(analog_tf_transform_result) :: tf
    type(analog_zpk_transform_result) :: zpk
    type(analog_prototype_result) :: proto
    complex(dp) :: z(1), p(2), s, h1, h2, gain
    integer :: j, n

    tf = bilinear([1.0_dp],[1.0_dp,1.0_dp],fs=2.0_dp)
    call check_tf(tf, [0.2_dp,0.2_dp], [1.0_dp,-0.6_dp])
    tf = bilinear([1.0_dp,2.0_dp],[1.0_dp,3.0_dp],fs=1.0_dp)
    call check_tf(tf, [4.0_dp/5.0_dp,0.0_dp],[1.0_dp,1.0_dp/5.0_dp])
    tf = bilinear([cmplx(1.0_dp,2.0_dp,dp)], &
        [cmplx(1.0_dp,0.0_dp,dp),cmplx(1.0_dp,-1.0_dp,dp)],fs=2.0_dp)
    if (tf%status /= frequency_transform_success) error stop 'complex bilinear'
    if (abs(tf%b(1)-cmplx(1.0_dp,2.0_dp,dp)/cmplx(5.0_dp,-1.0_dp,dp)) > 2.0e-14_dp) &
        error stop 'complex bilinear numerator'
    tf = bilinear([1.0_dp],[1.0_dp],fs=1.0_dp)
    call check_tf(tf,[1.0_dp],[1.0_dp])
    tf = bilinear([0.0_dp],[1.0_dp,1.0_dp])
    if (tf%status /= frequency_transform_success) error stop 'zero numerator'
    if (maxval(abs(tf%b)) /= 0.0_dp) error stop 'zero numerator coefficients'
    tf = bilinear([1.0_dp],[1.0_dp,1.0_dp],fs=0.0_dp)
    if (tf%status /= frequency_transform_invalid_input) error stop 'invalid sampling rate'
    tf = bilinear([1.0_dp],[1.0_dp,1.0_dp], &
        fs=ieee_value(0.0_dp,ieee_quiet_nan))
    if (tf%status /= frequency_transform_invalid_input) error stop 'NaN rate'
    tf = bilinear([1.0_dp],[1.0_dp,0.0_dp,1.0_dp],fs=-1.0_dp)
    if (tf%status /= frequency_transform_invalid_input) error stop 'negative rate'
    tf = bilinear([1.0_dp],[real(dp) ::],fs=1.0_dp)
    if (tf%status /= frequency_transform_invalid_input) error stop 'empty denominator'
    tf = bilinear([1.0_dp],[(1.0_dp,j=1,34)])
    if (tf%status /= frequency_transform_too_expensive) error stop 'large polynomial'

    z = [cmplx(-2.0_dp,0.0_dp,dp)]
    p = [cmplx(-1.0_dp,0.5_dp,dp),cmplx(-1.0_dp,-0.5_dp,dp)]
    gain = cmplx(2.0_dp,0.0_dp,dp)
    zpk = bilinear_zpk(z,p,gain,10.0_dp)
    if (zpk%status /= frequency_transform_success) error stop 'ZPK status'
    if (size(zpk%zeros) /= 2 .or. size(zpk%poles) /= 2) error stop 'ZPK lengths'
    if (abs(zpk%zeros(2)+1.0_dp) > 1.0e-15_dp) error stop 'Nyquist zero'
    do j = 1, 24
        s = cmplx(cos(0.11_dp*j),sin(0.11_dp*j),dp)
        h1 = eval_zpk(z,p,gain,20.0_dp*(s-1.0_dp)/(s+1.0_dp))
        h2 = eval_zpk(zpk%zeros,zpk%poles,zpk%gain,s)
        if (abs(h1-h2) > 1.0e-12_dp*max(1.0_dp,abs(h1))) &
            error stop 'ZPK frequency response'
    end do
    zpk = bilinear_zpk(z,p,cmplx(2.0_dp,0.3_dp,dp),10.0_dp)
    if (zpk%status /= frequency_transform_success) error stop 'complex ZPK gain'
    zpk = bilinear_zpk([p,z],p,2.0_dp,10.0_dp)
    if (zpk%status /= frequency_transform_invalid_input) error stop 'improper ZPK'
    zpk = bilinear_zpk(z,[cmplx(20.0_dp,0.0_dp,dp),p],2.0_dp,10.0_dp)
    if (zpk%status /= frequency_transform_invalid_input) error stop 'ZPK singular transform'
    zpk = bilinear_zpk(z,p,2.0_dp,0.0_dp)
    if (zpk%status /= frequency_transform_invalid_input) error stop 'ZPK invalid rate'

    proto = ellipap(0,1.0_dp,40.0_dp)
    if (proto%status /= prototype_success) error stop 'elliptic order zero'
    if (size(proto%zeros) /= 0 .or. size(proto%poles) /= 0) error stop 'order-zero roots'
    if (abs(proto%gain-10.0_dp**(-0.05_dp)) > 1.0e-14_dp) error stop 'elliptic order-zero gain'
    do n = 1, 12
        proto = ellipap(n,0.5_dp,45.0_dp)
        if (proto%status /= prototype_success) error stop 'elliptic prototype status'
        if (size(proto%poles) /= n .or. size(proto%zeros) /= n-mod(n,2)) &
            error stop 'elliptic prototype lengths'
        if (any(real(proto%poles,dp) >= 0.0_dp)) error stop 'unstable elliptic poles'
        if (abs(abs(eval_zpk(proto%zeros,proto%poles, &
            cmplx(proto%gain,0.0_dp,dp),cmplx(0.0_dp,0.0_dp,dp)))- &
            10.0_dp**(-0.5_dp/20.0_dp*real(1-mod(n,2),dp))) > 1.0e-9_dp) &
            error stop 'elliptic DC normalization'
    end do
    proto = ellipap(-1,1.0_dp,40.0_dp)
    if (proto%status /= prototype_invalid_input) error stop 'negative elliptic order'
    proto = ellipap(21,1.0_dp,40.0_dp)
    if (proto%status /= prototype_too_expensive) error stop 'large elliptic order'
    proto = ellipap(4,0.0_dp,40.0_dp)
    if (proto%status /= prototype_invalid_input) error stop 'zero elliptic ripple'
    proto = ellipap(4,2.0_dp,1.0_dp)
    if (proto%status /= prototype_invalid_input) error stop 'elliptic invalid attenuation'
    proto = besselap(13,'phase')
    if (proto%status /= prototype_success) error stop 'Bessel order 13'
    proto = besselap(20,'phase')
    if (proto%status /= prototype_success) error stop 'Bessel order 20'
    proto = besselap(21,'phase')
    if (proto%status /= prototype_too_expensive) error stop 'Bessel order 21 guard'
    print '(a)', 'test_signal_bilinear_ellipap passed'
contains
    subroutine check_tf(out, b, a)
        type(analog_tf_transform_result), intent(in) :: out !! bilinear result
        real(dp), intent(in) :: b(:) !! expected numerator coefficients
        real(dp), intent(in) :: a(:) !! expected denominator coefficients
        if (out%status /= frequency_transform_success) error stop 'bilinear TF status'
        if (size(out%b) /= size(b) .or. size(out%a) /= size(a)) error stop 'bilinear TF sizes'
        if (any(abs(out%b-b) > 1.0e-13_dp)) error stop 'bilinear TF numerator'
        if (any(abs(out%a-a) > 1.0e-13_dp)) error stop 'bilinear TF denominator'
    end subroutine check_tf

    pure function eval_zpk(zv,pv,k,s) result(h)
        complex(dp), intent(in) :: zv(:) !! zeros
        complex(dp), intent(in) :: pv(:) !! poles
        complex(dp), intent(in) :: k !! gain
        complex(dp), intent(in) :: s !! evaluation coordinate
        complex(dp) :: h
        integer :: i
        h = k
        do i = 1, size(zv)
            h = h*(s-zv(i))
        end do
        do i = 1, size(pv)
            h = h/(s-pv(i))
        end do
    end function eval_zpk
end program test_signal_bilinear_ellipap
