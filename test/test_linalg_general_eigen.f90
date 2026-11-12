! SPDX-License-Identifier: MIT
! Native checks for scipy.linalg.eig/eigvals-style ordinary eigenproblems.
program test_linalg_general_eigen
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_general_eigen, only : eig, eigvals
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none

    real(dp) :: a(4,4), r(2,2), bad(2,3), too_large(41,41)
    real(dp) :: znorm, scale
    complex(dp) :: c(3,3), w(4), vr(4,4), vl(4,4), wc(3), vc(3,3), lc(3,3)
    complex(dp) :: roots(2), wshort(1), huge_w(41)
    complex(dp) :: emptyw(0), emptya(0,0), emptyv(0,0)
    integer :: i, status

    ! Real matrix with a conjugate pair; verify left and right eigenvectors.
    r = reshape([0.0_dp, 1.0_dp, -1.0_dp, 0.0_dp],[2,2])
    call eig(r, roots, status)
    if (status /= linalg_status_success) error stop 'rotation: status'
    if (abs(roots(1)-cmplx(0.0_dp,1.0_dp,dp)) > 1e-11_dp .and. &
        abs(roots(2)-cmplx(0.0_dp,1.0_dp,dp)) > 1e-11_dp) error stop 'rotation: +i'
    if (abs(roots(1)-cmplx(0.0_dp,-1.0_dp,dp)) > 1e-11_dp .and. &
        abs(roots(2)-cmplx(0.0_dp,-1.0_dp,dp)) > 1e-11_dp) error stop 'rotation: -i'

    ! Real lower triangular, deliberately nonnormal with repeated values.
    a = 0.0_dp
    do i = 1, 4
        a(i,i) = 1.0_dp + real(i,dp)/3.0_dp
    end do
    a(2,1) = 2.0_dp
    a(3,1) = -4.0_dp
    a(3,2) = 3.0_dp
    a(4,3) = 0.5_dp
    call eig(a, w, status, vr=vr, vl=vl)
    if (status /= linalg_status_success) error stop 'lower triangular: status'
    do i = 1, 4
        if (abs(w(i)-cmplx(a(5-i,5-i),0.0_dp,dp)) > 1e-12_dp) &
            error stop 'lower triangular: eigenvalues'
        scale = max(1.0_dp,maxval(abs(a)))
        znorm = maxval(abs(matmul(cmplx(a,kind=dp),vr(:,i))-w(i)*vr(:,i)))/scale
        if (znorm > 1e-10_dp) error stop 'right eigenvector residual'
        znorm = maxval(abs(matmul(transpose(cmplx(a,kind=dp)),vl(:,i))-conjg(w(i))*vl(:,i)))/scale
        if (znorm > 1e-10_dp) error stop 'left eigenvector residual'
        if (abs(sum(abs(vr(:,i))**2)-1.0_dp) > 1e-12_dp) error stop 'right normalization'
        if (abs(sum(abs(vl(:,i))**2)-1.0_dp) > 1e-12_dp) error stop 'left normalization'
    end do
    call eigvals(a,w,status)
    if (status /= linalg_status_success) error stop 'eigvals'

    ! Complex nonnormal case and zero matrix.
    c = cmplx(0.0_dp,0.0_dp,dp)
    c(1,1) = cmplx(1.0_dp,2.0_dp,dp)
    c(2,2) = cmplx(-2.0_dp,1.0_dp,dp)
    c(3,3) = cmplx(3.0_dp,-1.0_dp,dp)
    c(1,2) = cmplx(2.0_dp,-0.2_dp,dp)
    c(2,3) = cmplx(-1.0_dp,1.0_dp,dp)
    c(3,1) = cmplx(0.2_dp,0.3_dp,dp)
    call eig(c,wc,status,vr=vc,vl=lc)
    if (status /= linalg_status_success) error stop 'complex matrix: status'
    do i = 1, 3
        if (maxval(abs(matmul(c,vc(:,i))-wc(i)*vc(:,i))) > 1e-9_dp) &
            error stop 'complex right eigenvectors'
        if (maxval(abs(matmul(transpose(conjg(c)),lc(:,i))-conjg(wc(i))*lc(:,i))) > 1e-9_dp) &
            error stop 'complex left eigenvectors'
    end do
    c = cmplx(0.0_dp,0.0_dp,dp)
    call eig(c,wc,status,vr=vc,vl=lc)
    if (status /= linalg_status_success .or. any(abs(wc) > 0.0_dp)) &
        error stop 'zero matrix eigenvalues'
    do i=1,3
        if (abs(vc(i,i)-cmplx(1.0_dp,0.0_dp,dp)) > 0.0_dp) error stop 'zero matrix vectors'
    end do
    call eig(emptya,emptyw,status,vr=emptyv)
    if (status /= linalg_status_success) error stop 'zero-order matrices'

    ! Invalid shapes and nonfinite inputs fail without invalid arithmetic.
    bad = 0.0_dp
    call eigvals(bad,roots,status)
    if (status /= linalg_status_invalid_shape) error stop 'nonsquare accepted'
    call eigvals(r,wshort,status)
    if (status /= linalg_status_invalid_shape) error stop 'short outputs accepted'
    too_large = 0.0_dp
    call eigvals(too_large,huge_w,status)
    if (status /= linalg_status_invalid_shape) error stop 'too large accepted'
    r(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    call eigvals(r,roots,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite accepted'

    print '(a)', 'general eig native tests passed'
end program test_linalg_general_eigen
