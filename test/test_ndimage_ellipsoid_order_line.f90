! SPDX-License-Identifier: MIT
! Native API and error checks for Fourier ellipsoids, mixed ranks, line callbacks.
program test_ndimage_ellipsoid_order_line
    use scifort_kinds, only : dp
    use scifort_ndimage, only : fourier_ellipsoid_nd, median_filter_nd, &
        rank_filter_nd, percentile_filter_nd, generic_filter1d_nd, &
        ndimage_status_invalid_input, ndimage_status_too_expensive
    implicit none
    complex(dp) :: frequencies(4),result(4)
    complex(dp) :: spectrum15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    complex(dp) :: filtered15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: a(3,4),b(3,4),filtered(3,4)
    real(dp) :: image15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: result15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    character(len=16) :: modes(2),high_modes(15)
    integer :: i,j,s

    frequencies=(1.0_dp,0.0_dp)
    call fourier_ellipsoid_nd(frequencies,2.0_dp,result,status=s)
    if (s/=0) error stop '1D ellipsoid status'
    if (abs(result(1)-cmplx(1.0_dp,0.0_dp,dp))>1.0e-14_dp) &
        error stop 'ellipsoid must preserve DC'
    if (abs(real(result(2),dp)-2.0_dp/acos(-1.0_dp))>1.0e-13_dp) &
        error stop '1D ellipsoid sine response'
    call fourier_ellipsoid_nd(frequencies,0.0_dp,result,status=s)
    if (s/=0 .or. any(result/=frequencies)) error stop 'zero ellipsoid size identity'
    call fourier_ellipsoid_nd(frequencies,-2.0_dp,result,status=s)
    if (s/=0 .or. abs(real(result(2),dp)-2.0_dp/acos(-1.0_dp))>1.0e-13_dp) &
        error stop 'negative diameter symmetry'
    call fourier_ellipsoid_nd(frequencies,[1.0_dp,2.0_dp],result,status=s)
    if (s/=ndimage_status_invalid_input) error stop 'ellipsoid parameter length'
    call fourier_ellipsoid_nd(frequencies,2.0_dp,result,status=s,n=4)
    if (s/=ndimage_status_invalid_input) error stop 'packed FFT shape'
    call fourier_ellipsoid_nd(frequencies,2.0_dp,result,status=s,axis=2)
    if (s/=ndimage_status_invalid_input) error stop 'FFT axis bounds'
    call fourier_ellipsoid_nd(frequencies,2.0_dp,result,status=s,n=0)
    if (s/=ndimage_status_invalid_input) error stop 'packed FFT length zero'

    modes=['nearest         ','constant        ']
    do j=1,4
        do i=1,3
            a(i,j)=real(2*i-5*j,dp)
        end do
    end do
    call median_filter_nd(a,[1,1],b,status=s,modes=modes)
    if (s/=0 .or. any(a/=b)) error stop 'mixed median identity'
    call rank_filter_nd(a,-1,[1,1],b,status=s,modes=modes)
    if (s/=0 .or. any(a/=b)) error stop 'mixed rank identity'
    call percentile_filter_nd(a,50.0_dp,[1,1],b,status=s,modes=modes)
    if (s/=0 .or. any(a/=b)) error stop 'mixed percentile identity'
    call median_filter_nd(a,[2,3],b,status=s,modes=modes,origin=[-1,1])
    if (s/=0) error stop 'median asymmetric mixed stencil'
    call percentile_filter_nd(a,75.0_dp,[2,3],b,status=s,modes=modes,origin=[-1,1])
    if (s/=0) error stop 'percentile asymmetric mixed stencil'
    call rank_filter_nd(a,2,[2,3],b,status=s,modes=modes,origin=[-1,1])
    if (s/=0) error stop 'rank asymmetric mixed stencil'
    call median_filter_nd(a,[2,3],b,status=s,modes=['invalid         ','nearest         '])
    if (s/=ndimage_status_invalid_input) error stop 'invalid per-axis median mode'
    call median_filter_nd(a,[2,3],b,status=s,mode='nearest',modes=modes)
    if (s/=ndimage_status_invalid_input) error stop 'median conflicting modes'
    call rank_filter_nd(a,1,[2,3],b,status=s,modes=['nearest         '])
    if (s/=ndimage_status_invalid_input) error stop 'rank wrong number of modes'
    call percentile_filter_nd(a,50.0_dp,[2,3],b,status=s,mode='reflect',modes=modes)
    if (s/=ndimage_status_invalid_input) error stop 'percentile conflicting modes'

    call generic_filter1d_nd(a,copy_axis,1,1,b,status=s)
    if (s/=0 .or. any(b/=a)) error stop '1D callback identity'
    call generic_filter1d_nd(a,copy_axis,2,2,b,status=s,mode='constant',cval=-9.0_dp,origin=0)
    if (s/=0) error stop '1D callback padded line'
    if (any(b(:,1)/= -9.0_dp) .or. any(b(:,2:)/=a(:,:3))) &
        error stop '1D callback constant left edge'
    call generic_filter1d_nd(a,copy_axis,1,3,b,status=s)
    if (s/=ndimage_status_invalid_input) error stop 'invalid callback axis'
    call generic_filter1d_nd(a,copy_axis,0,1,b,status=s)
    if (s/=ndimage_status_invalid_input) error stop 'invalid callback filter size'
    call generic_filter1d_nd(a,copy_axis,2,1,b,status=s,origin=1)
    if (s/=ndimage_status_invalid_input) error stop 'invalid callback origin'
    call generic_filter1d_nd(a,copy_axis,2,1,b,status=s,mode='unknown')
    if (s/=ndimage_status_invalid_input) error stop 'invalid callback mode'
    call generic_filter1d_nd(a,copy_axis,100001,1,b,status=s)
    if (s/=ndimage_status_too_expensive) error stop 'oversized callback stencil guard'

    image15=4.0_dp
    high_modes='nearest'
    high_modes(1)='wrap'
    call median_filter_nd(image15,spread(1,1,15),result15,status=s,modes=high_modes)
    if (s/=0 .or. any(result15/=image15)) error stop 'rank15 median mixed'
    call rank_filter_nd(image15,0,spread(1,1,15),result15,status=s,modes=high_modes)
    if (s/=0 .or. any(result15/=image15)) error stop 'rank15 rank mixed'
    call percentile_filter_nd(image15,80.0_dp,spread(1,1,15),result15,status=s,modes=high_modes)
    if (s/=0 .or. any(result15/=image15)) error stop 'rank15 percentile mixed'
    call generic_filter1d_nd(image15,copy_axis,1,1,result15,status=s)
    if (s/=0 .or. any(result15/=image15)) error stop 'rank15 line callback identity'
    spectrum15=cmplx(2.0_dp,-3.0_dp,dp)
    call fourier_ellipsoid_nd(spectrum15,1.5_dp,filtered15,status=s)
    if (s/=0) error stop 'rank15 ellipsoid status'
    if (abs(filtered15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)-spectrum15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)) &
        >1.0e-14_dp) error stop 'rank15 Fourier DC'
    if (abs(filtered15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1))>=abs(spectrum15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1))) &
        error stop 'rank15 ellipsoid must attenuate nonzero frequency'

    print '(a)', 'PASS: Fourier ellipsoid, per-axis order statistics and line callback regression'
contains
    subroutine copy_axis(input_line,output_line)
        real(dp), intent(in) :: input_line(:)
        real(dp), intent(out) :: output_line(:)
        output_line=input_line(1:size(output_line))
    end subroutine copy_axis
end program test_ndimage_ellipsoid_order_line
