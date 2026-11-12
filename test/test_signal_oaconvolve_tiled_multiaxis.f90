! SPDX-License-Identifier: MIT
! Check reusable-spectrum overlap-add with multiple tiled dimensions.
program test_signal_oaconvolve_tiled_multiaxis
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_signal_fftconvolve, only : oaconvolve_nd, fftconvolve_nd
    use scifort_signal_convolve, only : convolution_status_success, &
        convolution_status_invalid_input
    implicit none
    real(dp) :: a(129,65), b(3,5), full(131,69), same(129,65)
    real(dp) :: valid(127,61), ref(131,69), ref_same(129,65)
    real(dp) :: a3(129,2,65), b3(3,1,5), out3(131,2,69), ref3(131,2,69)
    real(dp) :: swapped_a(129,3), swapped_b(5,131)
    real(dp) :: swapped_out(133,133), swapped_ref(133,133)
    complex(dp) :: ac(129,65), bc(3,5), complex_out(131,69)
    complex(dp) :: complex_ref(131,69)
    real(dp) :: rank15(1,1,1,1,129,1,1,1,1,1,1,1,1,1,65)
    real(dp) :: kernel15(1,1,1,1,3,1,1,1,1,1,1,1,1,1,5)
    real(dp) :: output15(1,1,1,1,131,1,1,1,1,1,1,1,1,1,69)
    integer :: code, i, j

    a = reshape([(sin(0.087_dp*i)+cos(0.11_dp*i), i=1,size(a))],shape(a))
    b = reshape([(cos(0.23_dp*i)-0.3_dp*sin(0.31_dp*i), i=1,size(b))],shape(b))
    call oaconvolve_nd(a,b,full,axes=[1,2],status=code)
    if (code /= convolution_status_success) error stop 'multiple tiled axes failed'
    call fftconvolve_nd(a,b,ref,axes=[1,2],status=code)
    if (code /= convolution_status_success) error stop 'reference failed'
    if (maxval(abs(full-ref)) > 1.0e-11_dp) error stop 'multiple tiled axes mismatch'
    call oaconvolve_nd(a,b,same,mode='same',axes=[1,2],status=code)
    if (code /= convolution_status_success) error stop 'same tiled axes failed'
    call fftconvolve_nd(a,b,ref_same,mode='same',axes=[1,2],status=code)
    if (maxval(abs(same-ref_same)) > 1.0e-11_dp) error stop 'same tiled axes mismatch'
    call oaconvolve_nd(a,b,valid,mode='valid',axes=[1,2],status=code)
    if (code /= convolution_status_success) error stop 'valid tiled axes failed'
    if (maxval(abs(valid-ref(3:129,5:65))) > 1.0e-11_dp) error stop 'valid tiled axes mismatch'
    call oaconvolve_nd(b,a,full,axes=[1,2],status=code)
    if (code /= convolution_status_success) error stop 'reversed tiled axes failed'
    if (maxval(abs(full-ref)) > 1.0e-11_dp) error stop 'reversed tiled axes mismatch'

    do j=1,2
        a3(:,j,:) = a*real(j,dp)
    end do
    b3(:,1,:) = b
    call oaconvolve_nd(a3,b3,out3,axes=[1,3],status=code)
    if (code /= convolution_status_success) error stop 'nonconvolved broadcast failed'
    call fftconvolve_nd(a3,b3,ref3,axes=[1,3],status=code)
    if (maxval(abs(out3-ref3)) > 1.0e-11_dp) error stop 'nonconvolved broadcast mismatch'

    ! Opposite operands are longer along different transformed dimensions:
    ! keep the kernel spectrum fixed and tile the most asymmetric side.
    swapped_a = a(:,1:3)
    swapped_b = reshape([(0.15_dp*sin(0.039_dp*i),i=1,size(swapped_b))],shape(swapped_b))
    call oaconvolve_nd(swapped_a,swapped_b,swapped_out,status=code)
    if (code /= convolution_status_success) error stop 'opposite long directions failed'
    call fftconvolve_nd(swapped_a,swapped_b,swapped_ref,status=code)
    if (maxval(abs(swapped_out-swapped_ref)) > 1.0e-11_dp) error stop 'opposite long directions mismatch'

    ac = cmplx(a,0.2_dp*a,dp)
    bc = cmplx(b,-0.4_dp*b,dp)
    call oaconvolve_nd(ac,bc,complex_out,status=code)
    if (code /= convolution_status_success) error stop 'complex axes failed'
    call fftconvolve_nd(ac,bc,complex_ref,status=code)
    if (maxval(abs(complex_out-complex_ref)) > 2.0e-11_dp) error stop 'complex axes mismatch'
    call oaconvolve_nd(a,bc,complex_out,status=code)
    if (code /= convolution_status_success) error stop 'mixed real complex failed'
    call fftconvolve_nd(a,bc,complex_ref,status=code)
    if (maxval(abs(complex_out-complex_ref)) > 2.0e-11_dp) error stop 'mixed real complex mismatch'

    rank15 = reshape(a,shape(rank15))
    kernel15 = reshape(b,shape(kernel15))
    call oaconvolve_nd(rank15,kernel15,output15,axes=[5,15],status=code)
    if (code /= convolution_status_success) error stop 'rank15 multi-axis failed'
    if (maxval(abs(reshape(output15,[size(full)])-reshape(ref,[size(ref)]))) > &
        1.0e-11_dp) error stop 'rank15 multi-axis mismatch'
    call oaconvolve_nd(a,b,full,axes=[1,1],status=code)
    if (code /= convolution_status_invalid_input) error stop 'invalid axes accepted'
    if (.not. all(ieee_is_nan(full))) error stop 'invalid axes not NaN'

    print '(a)', 'test_signal_oaconvolve_tiled_multiaxis passed'
end program test_signal_oaconvolve_tiled_multiaxis
