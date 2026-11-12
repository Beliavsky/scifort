! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_envelope_rank15_into
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal, only : envelope, envelope_into, &
        envelope_status_success, envelope_status_invalid_input
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp), tol = 1.0e-10_dp
    integer, parameter :: n = 32
    real(dp) :: x4(2,n,3,2), y4(2,16,3,2), r4(2,16,3,2), xs(2,n,3,2)
    complex(dp) :: z4(2,n,3,2), cr4(2,16,3,2)
    real(dp) :: x15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,n)
    real(dp) :: y15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,16)
    real(dp) :: r15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,16)
    complex(dp) :: z15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,n)
    complex(dp) :: cr15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,16)
    real(dp) :: zero_x(0,n,1,1), zero_y(0,16,1,1)
    real(dp), allocatable :: expected(:)
    integer :: i, j, k, t, st
    real(dp) :: amp, baseline, phi, bad_scalar_residual

    do k=1,2
        do j=1,3
            do i=1,2
                amp = real(i+j+k,dp)
                baseline = real(j+2*k,dp)
                do t=1,n
                    phi = 2.0_dp*pi*5.0_dp*real(t-1,dp)/real(n,dp)
                    x4(i,t,j,k) = baseline + amp*cos(phi)
                    z4(i,t,j,k) = cmplx(baseline,0.0_dp,dp) + &
                        cmplx(amp*cos(phi),amp*sin(phi),dp)
                end do
            end do
        end do
    end do
    y4 = -111.0_dp
    r4 = -222.0_dp
    call envelope_into(x4,y4,axis=2,n_out=16,remainder=r4,status=st)
    if (st /= envelope_status_success) error stop 'rank4 real status'
    do k=1,2
        do j=1,3
            do i=1,2
                if (maxval(abs(y4(i,:,j,k)-real(i+j+k,dp))) > tol) &
                    error stop 'rank4 real envelope values'
                if (maxval(abs(r4(i,:,j,k)-real(j+2*k,dp))) > tol) &
                    error stop 'rank4 real residual values'
            end do
        end do
    end do
    call envelope_into(z4,y4,axis=-3,n_out=16,residual='all',remainder=cr4,status=st)
    if (st /= envelope_status_success) error stop 'rank4 complex status'
    do k=1,2
        do j=1,3
            do i=1,2
                if (maxval(abs(y4(i,:,j,k)-real(i+j+k,dp))) > tol) &
                    error stop 'rank4 complex envelope values'
                if (maxval(abs(cr4(i,:,j,k)-cmplx(real(j+2*k,dp),0.0_dp,dp))) > tol) &
                    error stop 'rank4 complex residual values'
            end do
        end do
    end do
    r4 = -31.0_dp
    call envelope_into(x4,y4,axis=2,n_out=16,residual='none',remainder=r4,status=st)
    if (st /= envelope_status_success .or. any(r4 /= -31.0_dp)) &
        error stop 'none mode residual preservation'

    do i=1,2
        x15(i,1,1,1,1,1,1,1,1,1,1,1,1,1,:) = x4(i,:,1,1)
        z15(i,1,1,1,1,1,1,1,1,1,1,1,1,1,:) = z4(i,:,1,1)
    end do
    y15 = -11.0_dp
    r15 = -12.0_dp
    call envelope_into(x15,y15,axis=-1,n_out=16,remainder=r15,status=st)
    if (st /= envelope_status_success) error stop 'rank15 status'
    do i=1,2
        expected = envelope(x4(i,:,1,1),n_out=16,status=st)
        if (st /= envelope_status_success) error stop 'vector envelope status'
        if (maxval(abs(y15(i,1,1,1,1,1,1,1,1,1,1,1,1,1,:)-expected)) > tol) &
            error stop 'rank15 versus vector'
    end do
    y15 = -10.0_dp
    cr15 = cmplx(-7.0_dp,-8.0_dp,dp)
    call envelope_into(z15,y15,axis=-1,n_out=16,residual='all',remainder=cr15,status=st)
    if (st /= envelope_status_success) error stop 'rank15 complex status'
    do i=1,2
        if (maxval(abs(y15(i,1,1,1,1,1,1,1,1,1,1,1,1,1,:)-real(i+2,dp))) > tol) &
            error stop 'rank15 complex values'
        if (maxval(abs(cr15(i,1,1,1,1,1,1,1,1,1,1,1,1,1,:)- &
            cmplx(3.0_dp,0.0_dp,dp))) > tol) error stop 'rank15 complex residual'
    end do
    zero_y = -1.0_dp
    call envelope_into(zero_x,zero_y,axis=2,n_out=16,status=st)
    if (st /= envelope_status_success) error stop 'empty batch status'
    y15 = -13.0_dp
    call envelope_into(x15(1:2:2,:,:,:,:,:,:,:,:,:,:,:,:,:,:), &
        y15(1:2:2,:,:,:,:,:,:,:,:,:,:,:,:,:,:), n_out=16,status=st)
    if (st /= envelope_status_success) error stop 'rank15 strided status'
    if (maxval(abs(y15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,:)-3.0_dp)) > tol) &
        error stop 'rank15 strided values'
    if (any(y15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,:) /= -13.0_dp)) &
        error stop 'rank15 untouched second channel'

    xs = x4
    xs(2,n,3,2) = ieee_value(0.0_dp,ieee_quiet_nan)
    y4 = -9.0_dp
    r4 = -8.0_dp
    call envelope_into(xs,y4,axis=2,n_out=16,remainder=r4,status=st)
    if (st /= envelope_status_invalid_input) error stop 'invalid later channel status'
    if (any(y4 /= -9.0_dp) .or. any(r4 /= -8.0_dp)) &
        error stop 'later channel modified outputs'
    call envelope_into(x4,y4,axis=2,n_out=16,bp_low=17,status=st)
    if (st /= envelope_status_invalid_input .or. any(y4 /= -9.0_dp)) &
        error stop 'invalid band changed output'
    call envelope_into(x4,y4,axis=0,n_out=16,status=st)
    if (st /= envelope_status_invalid_input .or. any(y4 /= -9.0_dp)) &
        error stop 'invalid axis changed output'
    call envelope_into(x4,y4,axis=2,n_out=8,status=st)
    if (st /= envelope_status_invalid_input .or. any(y4 /= -9.0_dp)) &
        error stop 'invalid output shape changed output'
    call envelope_into(x4,y4,axis=2,n_out=16,remainder=r15,status=st)
    if (st /= envelope_status_invalid_input .or. any(y4 /= -9.0_dp)) &
        error stop 'invalid residual shape changed output'
    bad_scalar_residual = -45.0_dp
    call envelope_into(x4,y4,axis=2,n_out=16,remainder=bad_scalar_residual,status=st)
    if (st /= envelope_status_invalid_input .or. any(y4 /= -9.0_dp)) &
        error stop 'scalar residual rank accepted'
    if (bad_scalar_residual /= -45.0_dp) error stop 'scalar residual modified'
    print *, 'PASS test_signal_envelope_rank15_into'
end program test_signal_envelope_rank15_into
