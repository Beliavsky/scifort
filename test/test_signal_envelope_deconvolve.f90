! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_envelope_deconvolve
    use scifort_kinds, only : dp
    use scifort_signal_envelope, only : envelope, envelope_status_success, envelope_status_invalid_input
    use scifort_signal_deconvolve, only : deconvolve, &
        deconvolution_status_success, deconvolution_status_invalid_input
    implicit none
    integer, parameter :: n = 32
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: tol = 2.0e-11_dp
    integer :: j, k, i, status
    real(dp) :: x(n), x2(n,2), x3(2,3,n), amplitude
    real(dp), allocatable :: q(:), rem(:), env(:), residue(:)
    real(dp), allocatable :: env2(:,:), residue2(:,:), env3(:,:,:), residue3(:,:,:)
    complex(dp) :: z(n), divisor(2), signal(4), known_q(3)
    complex(dp), allocatable :: cq(:), crem(:), cres(:), zenv(:)

    call deconvolve([1.0_dp,3.0_dp,3.0_dp,1.0_dp], [1.0_dp,1.0_dp], q, rem, status)
    if (status /= deconvolution_status_success) error stop 'deconvolve real status'
    if (maxval(abs(q-[1.0_dp,2.0_dp,1.0_dp])) > tol) error stop 'deconvolve quotient'
    if (maxval(abs(rem)) > tol) error stop 'deconvolve remainder'

    call deconvolve([1.0_dp,2.0_dp], [1.0_dp,2.0_dp,3.0_dp], q, rem, status)
    if (status /= deconvolution_status_success) error stop 'deconvolve overlong divisor'
    if (size(q) /= 0 .or. size(rem) /= 2) error stop 'overlong divisor shape'
    if (maxval(abs(rem-[1.0_dp,2.0_dp])) > tol) error stop 'overlong divisor remainder'

    divisor = [cmplx(1.0_dp,0.5_dp,dp), cmplx(-0.3_dp,0.8_dp,dp)]
    known_q = [cmplx(2.0_dp,1.0_dp,dp), cmplx(-1.0_dp,2.0_dp,dp), &
        cmplx(0.5_dp,-0.7_dp,dp)]
    signal = cmplx(0.0_dp,0.0_dp,dp)
    do j = 1, size(known_q)
        do k = 1, size(divisor)
            signal(j+k-1) = signal(j+k-1) + known_q(j)*divisor(k)
        end do
    end do
    call deconvolve(signal, divisor, cq, crem, status)
    if (status /= deconvolution_status_success) error stop 'deconvolve complex status'
    if (maxval(abs(cq-known_q)) > tol) error stop 'deconvolve complex quotient'
    if (maxval(abs(crem)) > tol) error stop 'deconvolve complex remainder'

    call deconvolve([1.0_dp,2.0_dp,3.0_dp], divisor, cq, crem, status)
    if (status /= deconvolution_status_success) error stop 'deconvolve mixed real-complex'
    call deconvolve(signal, [1.0_dp,2.0_dp], cq, crem, status)
    if (status /= deconvolution_status_success) error stop 'deconvolve mixed complex-real'

    call deconvolve([1.0_dp,2.0_dp], [0.0_dp,1.0_dp], q, rem, status)
    if (status /= deconvolution_status_invalid_input) error stop 'leading zero divisor'
    if (size(q) /= 0 .or. size(rem) /= 0) error stop 'invalid deconvolution shape'
    call deconvolve([1.0_dp,2.0_dp], [real(dp) ::], q, rem, status)
    if (status /= deconvolution_status_invalid_input) error stop 'empty deconvolution divisor'

    do j = 1, n
        x(j) = 2.0_dp + 3.0_dp*cos(2.0_dp*pi*3.0_dp*real(j-1,dp)/real(n,dp))
        z(j) = 2.0_dp*cmplx(cos(2.0_dp*pi*5.0_dp*real(j-1,dp)/real(n,dp)), &
            sin(2.0_dp*pi*5.0_dp*real(j-1,dp)/real(n,dp)),dp)
    end do
    env = envelope(x, remainder=residue, status=status)
    if (status /= envelope_status_success) error stop 'real envelope status'
    if (maxval(abs(env-3.0_dp)) > tol) error stop 'real carrier envelope'
    if (maxval(abs(residue-2.0_dp)) > tol) error stop 'real lowpass residue'
    env = envelope(x, n_out=16, squared=.true., residual='none', remainder=residue, status=status)
    if (status /= envelope_status_success .or. size(env) /= 16) error stop 'envelope downsampling'
    if (size(residue) /= 0) error stop 'no residual mode rank1'
    if (maxval(abs(env-9.0_dp)) > tol) error stop 'squared envelope'
    env = envelope(z, remainder=cres, residual='all', status=status)
    if (status /= envelope_status_success) error stop 'complex envelope status'
    if (maxval(abs(env-2.0_dp)) > tol) error stop 'complex analytic envelope'
    if (maxval(abs(cres)) > tol) error stop 'complex all residual'

    do j = 1, n
        z(j) = 3.0_dp*cmplx(cos(-2.0_dp*pi*5.0_dp*real(j-1,dp)/real(n,dp)), &
            sin(-2.0_dp*pi*5.0_dp*real(j-1,dp)/real(n,dp)),dp)
    end do
    env = envelope(z, bp_low=-6, bp_high=-4, squared=.true., residual='all', &
        remainder=cres, status=status)
    if (status /= envelope_status_success) error stop 'negative band status'
    if (maxval(abs(env-9.0_dp)) > tol) error stop 'negative frequency band'
    if (maxval(abs(cres)) > tol) error stop 'negative band residual'

    x2(:,1) = x
    x2(:,2) = 2.0_dp*x
    env2 = envelope(x2, axis=1, remainder=residue2, status=status)
    if (status /= envelope_status_success) error stop '2D axis1 status'
    if (any(shape(env2) /= [n,2])) error stop '2D axis1 shape'
    if (maxval(abs(env2(:,1)-3.0_dp)) > tol) error stop '2D axis1 first channel'
    if (maxval(abs(env2(:,2)-6.0_dp)) > tol) error stop '2D axis1 second channel'
    if (maxval(abs(residue2(:,2)-4.0_dp)) > tol) error stop '2D axis1 residual'
    env2 = envelope(transpose(x2), axis=2, n_out=17, residual='none', &
        remainder=residue2, status=status)
    if (status /= envelope_status_success) error stop '2D axis2 status'
    if (any(shape(env2) /= [2,17])) error stop '2D resample shape'
    if (size(residue2) /= 0) error stop '2D none residual shape'
    if (maxval(abs(env2(1,:)-3.0_dp)) > tol) error stop '2D resampled envelope'

    do j=1,n
        do k=1,3
            do i=1,2
                amplitude = real(i+k,dp)
                x3(i,k,j) = amplitude*x(j)
            end do
        end do
    end do
    env3 = envelope(x3, n_out=16, axis=3, remainder=residue3, status=status)
    if (status /= envelope_status_success) error stop '3D envelope status'
    if (any(shape(env3) /= [2,3,16])) error stop '3D envelope shape'
    if (maxval(abs(env3(2,3,:)-15.0_dp)) > tol) error stop '3D envelope channel'
    if (maxval(abs(residue3(2,3,:)-10.0_dp)) > tol) error stop '3D residual channel'
    env3 = envelope(x3, axis=3, residual='none', remainder=residue3, status=status)
    if (status /= envelope_status_success .or. size(residue3) /= 0) &
        error stop '3D no residual'
    env3 = envelope(x3(1:2:2,:,:), axis=3, status=status)
    if (status /= envelope_status_success .or. size(env3,1) /= 1) &
        error stop '3D strided view'

    env = envelope(x, bp_low=16, bp_high=16, status=status)
    if (status /= envelope_status_invalid_input .or. size(env) /= 0) &
        error stop 'invalid envelope band'
    env2 = envelope(x2, axis=0, status=status)
    if (status /= envelope_status_invalid_input .or. size(env2) /= 0) &
        error stop 'invalid envelope axis'
    env3 = envelope(x3, n_out=0, status=status)
    if (status /= envelope_status_invalid_input .or. size(env3) /= 0) &
        error stop 'invalid envelope resample length'
    env2 = envelope(x2, residual='badmode', status=status)
    if (status /= envelope_status_invalid_input .or. size(env2) /= 0) &
        error stop 'invalid envelope residual option'
    env2 = envelope(reshape([real(dp) :: ],[0,32]),axis=2,bp_low=20,status=status)
    if (status /= envelope_status_invalid_input) error stop 'empty batch invalid band'
    env2 = envelope(reshape([real(dp) :: ],[0,32]),axis=2,status=status)
    if (status /= envelope_status_success .or. any(shape(env2)/=[0,32])) &
        error stop 'zero-channel 2D envelope'

    print *, 'PASS test_signal_envelope_deconvolve'
end program test_signal_envelope_deconvolve
