! SPDX-License-Identifier: MIT
! Native regression tests for real-to-complex Schur conversion.
program test_linalg_rsf2csf
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, &
        linalg_status_invalid_shape, linalg_status_no_convergence
    use scifort_linalg_schur, only : schur, rsf2csf
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none

    real(dp) :: t2(2,2), z2(2,2), a4(4,4), t4(4,4), z4(4,4)
    real(dp) :: saved2(2,2), no_matrix(0,0), wrong(2,3), one(1,1)
    real(dp) :: scaled(4,4)
    complex(dp) :: tc2(2,2), zc2(2,2), tc4(4,4), zc4(4,4)
    complex(dp) :: empty_tc(0,0), empty_zc(0,0), wrong_tc(2,3)
    integer :: status, k

    z2 = 0.0_dp
    z2(1,1) = 1.0_dp
    z2(2,2) = 1.0_dp
    t2(1,:) = [0.0_dp, -2.0_dp]
    t2(2,:) = [2.0_dp, 0.0_dp]
    saved2 = t2
    call rsf2csf(t2,z2,tc2,zc2,status)
    if (status /= linalg_status_success) error stop '2x2 rotation conversion'
    if (any(t2 /= saved2)) error stop 'input was modified'
    call verify(cmplx(t2,0.0_dp,dp),tc2,zc2)
    if (abs(aimag(tc2(1,1))) < 1.9_dp) error stop 'complex pair lost'
    t2 = saved2*1.0e-260_dp
    call rsf2csf(t2,z2,tc2,zc2,status)
    if (status /= linalg_status_success) error stop 'tiny real Schur block'
    call verify(cmplx(t2,0.0_dp,dp),tc2,zc2)
    if (abs(aimag(tc2(1,1))/1.0e-260_dp) < 1.9_dp) &
        error stop 'tiny conjugate eigenvalues lost'
    t2 = saved2*1.0e260_dp
    call rsf2csf(t2,z2,tc2,zc2,status)
    if (status /= linalg_status_success) error stop 'large real Schur block'
    call verify(cmplx(t2,0.0_dp,dp),tc2,zc2)
    t2 = 0.0_dp
    call rsf2csf(t2,z2,tc2,zc2,status)
    if (status /= linalg_status_success .or. any(abs(tc2) > 0.0_dp)) &
        error stop 'zero real Schur form'
    t2(1,:) = [1.0_dp, 2.0_dp]
    t2(2,:) = [3.0_dp, 4.0_dp]
    call rsf2csf(t2,z2,tc2,zc2,status)
    if (status /= linalg_status_success) error stop '2x2 real eigenvalues'
    call verify(cmplx(t2,0.0_dp,dp),tc2,zc2)
    t2 = saved2

    a4 = 0.0_dp
    a4(1,1) = -3.0_dp
    a4(2,2) = 0.2_dp
    a4(2,3) = -4.0_dp
    a4(3,2) = 4.0_dp
    a4(3,3) = 0.2_dp
    a4(4,4) = 8.0_dp
    a4(1,2) = 0.7_dp
    a4(1,3) = 0.6_dp
    a4(2,4) = -2.1_dp
    a4(3,4) = 1.9_dp
    z4 = 0.0_dp
    do k = 1, 4
        z4(k,k) = 1.0_dp
    end do
    call rsf2csf(a4,z4,tc4,zc4,status)
    if (status /= linalg_status_success) error stop '4x4 block conversion'
    call verify(cmplx(a4,0.0_dp,dp),tc4,zc4)

    ! Exercise the real Schur -> complex Schur pipeline directly.
    a4(4,1) = 0.13_dp
    call schur(a4,t4,z4,status)
    if (status /= linalg_status_success) error stop 'real schur prerequisite'
    call rsf2csf(t4,z4,tc4,zc4,status)
    if (status /= linalg_status_success) error stop 'schur conversion pipeline'
    call verify(cmplx(a4,0.0_dp,dp),tc4,zc4)

    scaled = a4*1.0e145_dp
    call schur(scaled,t4,z4,status)
    if (status /= linalg_status_success) error stop 'scaled real schur'
    call rsf2csf(t4,z4,tc4,zc4,status)
    if (status /= linalg_status_success) error stop 'scaled conversion'
    call verify(cmplx(scaled,0.0_dp,dp),tc4,zc4)

    call rsf2csf(no_matrix,no_matrix,empty_tc,empty_zc,status)
    if (status /= linalg_status_success) error stop 'empty input'
    call rsf2csf(t2,wrong,tc2,zc2,status)
    if (status /= linalg_status_invalid_shape) error stop 'wrong basis accepted'
    call rsf2csf(t2,z2,wrong_tc,zc2,status)
    if (status /= linalg_status_invalid_shape) error stop 'wrong output accepted'

    t2(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    call rsf2csf(t2,z2,tc2,zc2,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite input accepted'
    t2 = saved2
    a4(3,1) = 0.2_dp
    call rsf2csf(a4,z4,tc4,zc4,status)
    if (status /= linalg_status_no_convergence) error stop 'non-Schur matrix accepted'
    a4 = 0.0_dp
    a4(2,1) = 1.0_dp
    a4(3,2) = 1.0_dp
    call rsf2csf(a4,z4,tc4,zc4,status)
    if (status /= linalg_status_no_convergence) error stop 'overlapping blocks accepted'

    one(1,1) = 4.0_dp
    call rsf2csf(one,reshape([1.0_dp],[1,1]),tc2,zc2,status)
    if (status /= linalg_status_invalid_shape) error stop 'mismatched size accepted'
    print '(a)', 'rsf2csf native tests passed'
contains
    subroutine verify(a,t,z)
        complex(dp), intent(in) :: a(:, :) !! original input matrix
        complex(dp), intent(in) :: t(:, :) !! output complex triangular factor
        complex(dp), intent(in) :: z(:, :) !! output unitary Schur vectors
        complex(dp) :: eye(size(a,1),size(a,1))
        real(dp) :: magnitude, rec_err, ortho_err
        integer :: n, j

        n = size(a,1)
        eye = cmplx(0.0_dp,0.0_dp,dp)
        do j = 1, n
            eye(j,j) = cmplx(1.0_dp,0.0_dp,dp)
        end do
        magnitude = max(tiny(1.0_dp),maxval(abs(a)))
        rec_err = maxval(abs(matmul(matmul(z,t/magnitude),transpose(conjg(z)))-a/magnitude))
        ortho_err = maxval(abs(matmul(transpose(conjg(z)),z)-eye))
        if (rec_err > 1.0e-8_dp) error stop 'reconstruction failed'
        if (ortho_err > 1.0e-8_dp) error stop 'unitarity failed'
        do j = 1, n-1
            if (any(abs(t(j+1:n,j)/magnitude) > 1.0e-8_dp)) error stop 'not triangular'
        end do
    end subroutine verify
end program test_linalg_rsf2csf
