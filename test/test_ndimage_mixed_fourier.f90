! SPDX-License-Identifier: MIT
! Per-axis boundary modes, complex callback and FFT-domain numerical checks.
program test_ndimage_mixed_fourier
    use scifort_kinds, only : dp
    use scifort_ndimage, only : gaussian_filter_nd, gaussian_filter1d_nd, uniform_filter_nd, &
        generic_filter_nd, minimum_filter_nd, maximum_filter_nd, correlate_nd, convolve_nd, &
        fourier_gaussian_nd, fourier_uniform_nd, fourier_shift_nd, &
        ndimage_status_invalid_input
    implicit none
    real(dp) :: a(3,4),b(3,4),c(3,4),work(3,4),h15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: r15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    complex(dp) :: z(3,4),w(3,4),cz15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    complex(dp) :: cw15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    complex(dp) :: spectrum(4),shifted(4),gaussian(4),uniform(4)
    character(len=16) :: modes(2),many_modes(15)
    integer :: status,i,j
    real(dp) :: pi

    pi=acos(-1.0_dp)
    do j=1,4
        do i=1,3
            a(i,j)=real(i+2*j,dp)
        end do
    end do
    z=cmplx(a,-2.0_dp*a,dp)
    modes=['nearest         ','constant        ']
    call gaussian_filter_nd(a,[0.8_dp,1.3_dp],b,status=status,modes=modes,cval=-0.2_dp)
    if (status/=0) error stop 'mixed Gaussian status'
    call gaussian_filter1d_nd(a,0.8_dp,1,work,status=status,mode='nearest',cval=-0.2_dp)
    if (status/=0) error stop 'axis one Gaussian'
    call gaussian_filter1d_nd(work,1.3_dp,2,c,status=status,mode='constant',cval=-0.2_dp)
    if (status/=0 .or. any(abs(b-c)>1e-12_dp)) error stop 'mixed Gaussian composition'
    call gaussian_filter_nd(a,0.0_dp,b,status=status,modes=modes)
    if (status/=0 .or. any(b/=a)) error stop 'scalar Gaussian mixed identity'
    call gaussian_filter_nd(z,[0.8_dp,1.3_dp],w,status=status,modes=modes)
    if (status/=0) error stop 'mixed complex Gaussian'
    call uniform_filter_nd(a,[1,1],b,status=status,modes=modes)
    if (status/=0 .or. any(abs(b-a)>1e-14_dp)) error stop 'mixed uniform identity'
    call uniform_filter_nd(z,[1,1],w,status=status,modes=modes)
    if (status/=0 .or. any(abs(w-z)>1e-14_dp)) error stop 'mixed complex uniform identity'
    call minimum_filter_nd(a,[1,1],b,status=status,modes=modes)
    if (status/=0 .or. any(b/=a)) error stop 'mixed minimum identity'
    call maximum_filter_nd(a,[1,1],b,status=status,modes=modes)
    if (status/=0 .or. any(b/=a)) error stop 'mixed maximum identity'
    call correlate_nd(a,reshape([1.0_dp],[1,1]),b,status=status,modes=modes)
    if (status/=0 .or. any(b/=a)) error stop 'mixed correlate identity'
    call convolve_nd(z,reshape([(1.0_dp,0.0_dp)],[1,1]),w,status=status,modes=modes)
    if (status/=0 .or. any(w/=z)) error stop 'mixed convolve identity'
    call generic_filter_nd(a,mean_real,b,status,sizes=[1,1],modes=modes)
    if (status/=0 .or. any(abs(b-a)>1e-14_dp)) error stop 'mixed generic real identity'
    call generic_filter_nd(z,mean_complex,w,status,sizes=[1,1],modes=modes)
    if (status/=0 .or. any(abs(w-z)>1e-14_dp)) error stop 'complex generic identity'
    call generic_filter_nd(z,mean_complex,w,status,sizes=[1,1],mode='bad')
    if (status/=ndimage_status_invalid_input) error stop 'complex callback invalid mode'
    call generic_filter_nd(z,mean_complex,w,status,sizes=[1,1],modes=['bad             ','nearest         '])
    if (status/=ndimage_status_invalid_input) error stop 'complex callback invalid axis mode'
    call generic_filter_nd(z,mean_complex,w,status,sizes=[1,1],mode='nearest',modes=modes)
    if (status/=ndimage_status_invalid_input) error stop 'complex callback duplicate mode'
    call gaussian_filter_nd(a,1.0_dp,b,status=status,modes=['nearest         '])
    if (status/=ndimage_status_invalid_input) error stop 'Gaussian wrong number of modes'
    call gaussian_filter_nd(a,1.0_dp,b,status=status,mode='nearest',modes=modes)
    if (status/=ndimage_status_invalid_input) error stop 'Gaussian duplicate mode'
    call uniform_filter_nd(z,3,w,status=status,modes=['nearest         ','invalid         '])
    if (status/=ndimage_status_invalid_input) error stop 'uniform invalid mode'
    call uniform_filter_nd(a,3,b,status=status,mode='nearest',modes=modes)
    if (status/=ndimage_status_invalid_input) error stop 'uniform duplicate mode'
    call maximum_filter_nd(a,2,b,status=status,mode='wrap',modes=modes)
    if (status/=ndimage_status_invalid_input) error stop 'maximum duplicate mode'
    call minimum_filter_nd(a,2,b,status=status,modes=['nearest         ','bad             '])
    if (status/=ndimage_status_invalid_input) error stop 'minimum invalid mode'
    call correlate_nd(a,reshape([1.0_dp],[1,1]),b,status=status,mode='nearest',modes=modes)
    if (status/=ndimage_status_invalid_input) error stop 'correlate duplicate mode'
    call convolve_nd(z,reshape([(1.0_dp,0.0_dp)],[1,1]),w,status=status, &
        modes=['nearest         ','invalid         '])
    if (status/=ndimage_status_invalid_input) error stop 'convolve invalid mode' 

    spectrum=(1.0_dp,0.0_dp)
    call fourier_shift_nd(spectrum,1.0_dp,shifted,status)
    if (status/=0) error stop 'Fourier shift status'
    if (abs(shifted(1)-(1.0_dp,0.0_dp))>1e-13_dp) error stop 'DC Fourier shift'
    if (abs(shifted(2)-cmplx(0.0_dp,-1.0_dp,dp))>1e-13_dp) error stop 'first harmonic shift'
    if (abs(shifted(3)-cmplx(-1.0_dp,0.0_dp,dp))>1e-13_dp) error stop 'Nyquist shift'
    if (abs(shifted(4)-cmplx(0.0_dp,1.0_dp,dp))>1e-13_dp) error stop 'negative frequency shift'
    call fourier_gaussian_nd(spectrum,1.0_dp,gaussian,status)
    if (status/=0) error stop 'Fourier Gaussian status'
    if (abs(gaussian(2)-cmplx(exp(-2*pi*pi/16.0_dp),0.0_dp,dp))>1e-13_dp) &
        error stop 'Gaussian frequency factor'
    call fourier_uniform_nd(spectrum,2.0_dp,uniform,status)
    if (status/=0) error stop 'Fourier uniform status'
    if (abs(uniform(2)-cmplx(2.0_dp/pi,0.0_dp,dp))>1e-13_dp) &
        error stop 'uniform sinc factor'
    call fourier_shift_nd(spectrum,[1.0_dp],shifted,status,n=8)
    if (status/=ndimage_status_invalid_input) error stop 'invalid RFFT spectrum size'
    call fourier_gaussian_nd(spectrum,1.0_dp,gaussian,status,n=0)
    if (status/=ndimage_status_invalid_input) error stop 'invalid FFT n=0'
    call fourier_uniform_nd(spectrum,1.0_dp,uniform,status,axis=2)
    if (status/=ndimage_status_invalid_input) error stop 'invalid FFT axis'
    call fourier_gaussian_nd(spectrum,-1.0_dp,gaussian,status)
    if (status/=ndimage_status_invalid_input) error stop 'invalid Gaussian sigma'

    h15=5.0_dp
    many_modes='nearest'
    many_modes(1)='wrap'
    call gaussian_filter_nd(h15,spread(0.0_dp,1,15),r15,status=status,modes=many_modes)
    if (status/=0 .or. any(r15/=h15)) error stop 'rank 15 Gaussian identity'
    call uniform_filter_nd(h15,spread(1,1,15),r15,status=status,modes=many_modes)
    if (status/=0 .or. any(r15/=h15)) error stop 'rank 15 uniform identity'
    cz15=cmplx(1.0_dp,-3.0_dp,dp)
    call generic_filter_nd(cz15,mean_complex,cw15,status, &
        sizes=spread(1,1,15),modes=many_modes)
    if (status/=0 .or. any(abs(cw15-cz15)>1e-14_dp)) error stop 'rank 15 complex generic'
    call fourier_shift_nd(cz15,0.0_dp,cw15,status)
    if (status/=0 .or. any(abs(cw15-cz15)>1e-14_dp)) error stop 'rank 15 Fourier shift'
    call fourier_gaussian_nd(cz15,0.0_dp,cw15,status)
    if (status/=0 .or. any(abs(cw15-cz15)>1e-14_dp)) error stop 'rank 15 Fourier Gaussian'
    print '(a)', 'PASS: mixed boundary modes, complex callbacks and Fourier transforms'
contains
    function mean_real(values) result(value)
        real(dp), intent(in) :: values(:)
        real(dp) :: value
        value=sum(values)/real(size(values),dp)
    end function mean_real

    function mean_complex(values) result(value)
        complex(dp), intent(in) :: values(:)
        complex(dp) :: value
        value=sum(values)/real(size(values),dp)
    end function mean_complex
end program test_ndimage_mixed_fourier
