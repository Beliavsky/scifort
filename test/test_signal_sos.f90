! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_sos
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_signal, only : dp, sosfilt, sosfilt_zi, lfilter, &
        sos_status_success, sos_status_invalid_input, sos_status_no_steady_state
    implicit none
    real(dp) :: sos(2,6), x(45), bad(2,6)
    real(dp), allocatable :: y(:), z(:), y1(:), y2(:), zi(:,:), zf(:,:), zi_next(:,:)
    complex(dp) :: csos(2,6), cx(45)
    complex(dp), allocatable :: cy(:), cz(:), czi(:,:), czf(:,:), cy2(:)
    integer :: i, status

    sos(1,:) = [0.2_dp,0.3_dp,0.1_dp,1.0_dp,-0.25_dp,0.08_dp]
    sos(2,:) = [1.0_dp,-0.1_dp,0.2_dp,1.0_dp,-0.3_dp,0.1_dp]
    x = [(0.05_dp*i + sin(0.24_dp*i), i=1,45)]
    y = sosfilt(sos,x,status=status)
    if (status /= sos_status_success) error stop 'sos real status'
    z = lfilter(sos(1,1:3),sos(1,4:6),x)
    z = lfilter(sos(2,1:3),sos(2,4:6),z)
    call assert_real(y,z,1.0e-13_dp)
    y1 = sosfilt(sos,x(:17),zf=zf)
    y2 = sosfilt(sos,x(18:),zi=zf,zf=zi_next)
    call assert_real(y,[y1,y2],1.0e-13_dp)
    if (size(zi_next,1) /= 2 .or. size(zi_next,2) /= 2) error stop 'sos state shape'

    zi = sosfilt_zi(sos,status=status)
    if (status /= sos_status_success) error stop 'sos step status'
    z = sosfilt(sos, spread(1.0_dp,1,45), zi=zi)
    if (any(abs(z-product(sum(sos(:,1:3),dim=2)/sum(sos(:,4:6),dim=2)))>1.0e-12_dp)) &
        error stop 'sos step steady state'

    csos = cmplx(sos,0.0_dp,dp)
    csos(1,2) = cmplx(0.3_dp,0.07_dp,dp)
    csos(1,5) = cmplx(-0.25_dp,0.03_dp,dp)
    cx = cmplx(x,0.3_dp*cos(x),dp)
    cy = sosfilt(csos,cx)
    cz = lfilter(csos(1,1:3),csos(1,4:6),cx)
    cz = lfilter(csos(2,1:3),csos(2,4:6),cz)
    call assert_complex(cy,cz,1.0e-13_dp)
    cy2 = sosfilt(csos,cx(:20),zf=czf)
    cz = sosfilt(csos,cx(21:),zi=czf)
    call assert_complex(cy,[cy2,cz],1.0e-13_dp)
    czi = sosfilt_zi(csos,status=status)
    if (status /= sos_status_success) error stop 'complex step state status'
    cy = sosfilt(csos, spread(cmplx(1.0_dp,0.0_dp,dp),1,45),zi=czi)
    if (any(abs(cy-product(sum(csos(:,1:3),dim=2)/sum(csos(:,4:6),dim=2)))>1.0e-12_dp)) &
        error stop 'complex step steady state'

    bad = sos
    bad(1,4) = 2.0_dp
    y = sosfilt(bad,x,status=status)
    if (status /= sos_status_invalid_input .or. .not. all(ieee_is_nan(y))) error stop 'invalid a0'
    y = sosfilt(sos,x,zi=reshape([1.0_dp,2.0_dp],[1,2]),status=status)
    if (status /= sos_status_invalid_input) error stop 'invalid zi dimensions'
    bad = sos
    bad(1,5) = -1.0_dp
    bad(1,6) = 0.0_dp
    zi = sosfilt_zi(bad,status=status)
    if (status /= sos_status_no_steady_state) error stop 'unit pole'
    y = sosfilt(sos,[real(dp)::],status=status)
    if (status /= sos_status_success .or. size(y) /= 0) error stop 'empty samples'
    print '(a)', 'test_signal_sos passed'
contains
    subroutine assert_real(v,w,tol)
        real(dp), intent(in) :: v(:) !! actual real array
        real(dp), intent(in) :: w(:) !! expected real array
        real(dp), intent(in) :: tol !! maximum absolute error
        if (size(v) /= size(w)) error stop 'real size'
        if (any(abs(v-w)>tol)) error stop 'real mismatch'
    end subroutine
    subroutine assert_complex(v,w,tol)
        complex(dp), intent(in) :: v(:) !! actual complex array
        complex(dp), intent(in) :: w(:) !! expected complex array
        real(dp), intent(in) :: tol !! maximum absolute error
        if (size(v) /= size(w)) error stop 'complex size'
        if (any(abs(v-w)>tol)) error stop 'complex mismatch'
    end subroutine
end program
