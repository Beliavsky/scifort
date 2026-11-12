! SPDX-License-Identifier: MIT
program test_linalg_banded
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_singular, linalg_status_not_positive_definite
    use scifort_linalg_banded, only : solve_banded, solveh_banded, cholesky_banded, cho_solve_banded
    implicit none
    integer :: status, i, j, n
    real(dp) :: a(4,4), ab(3,4), xtrue(4), b(4), x(4)
    real(dp) :: br(4,2), xr(4,2), ur(3,4), lr(3,4), cb_r(3,4)
    real(dp), allocatable :: ab_big(:,:), b_big(:), x_big(:)
    complex(dp) :: ac(4,4), abc(3,4), bc(4), xc(4), xtruth(4), lc(4,4)
    complex(dp) :: bc2(4,2), xc2(4,2), cb_c(3,4), ab_l(3,4), ab_u(3,4)

    a = 0.0_dp
    a(2,1)=2.0_dp
    a(1,2)=1.0_dp
    a(2,2)=4.0_dp
    a(3,2)=1.0_dp
    a(2,3)=-0.5_dp
    a(3,3)=5.0_dp
    a(4,3)=0.2_dp
    a(3,4)=0.6_dp
    a(4,4)=3.0_dp
    call pack_real(a,1,1,ab)
    xtrue = [2.0_dp, -1.0_dp, 0.3_dp, 0.8_dp]
    b = matmul(a,xtrue)
    call solve_banded(ab,1,1,b,x,status)
    if (status /= linalg_status_success) error stop 'real pivoting status'
    if (maxval(abs(x-xtrue)) > 1.0e-13_dp) error stop 'real pivoting accuracy'
    br(:,1) = b
    br(:,2) = matmul(a,[1.0_dp,2.0_dp,3.0_dp,4.0_dp])
    call solve_banded(ab,1,1,br,xr,status)
    if (status /= linalg_status_success) error stop 'real multirhs status'
    if (maxval(abs(xr(:,1)-xtrue)) > 1.0e-13_dp .or. &
            maxval(abs(xr(:,2)-[1.0_dp,2.0_dp,3.0_dp,4.0_dp])) > 1.0e-13_dp) &
        error stop 'real multirhs accuracy'

    ac = cmplx(a, kind=dp)
    ac(2,1) = cmplx(2.0_dp,0.2_dp,dp)
    ac(2,3) = cmplx(-0.5_dp,0.3_dp,dp)
    ac(3,4) = cmplx(0.6_dp,-0.2_dp,dp)
    call pack_complex(ac,1,1,abc)
    xtruth = [cmplx(2.0_dp,0.1_dp,dp), cmplx(-1.0_dp,0.3_dp,dp), &
        cmplx(0.3_dp,-0.7_dp,dp), cmplx(0.8_dp,0.1_dp,dp)]
    bc = matmul(ac,xtruth)
    call solve_banded(abc,1,1,bc,xc,status)
    if (status /= linalg_status_success .or. maxval(abs(xc-xtruth)) > 2.0e-13_dp) &
        error stop 'complex pivoting'

    ! Real SPD band: both SciPy lower and upper compact storage formats.
    a = 0.0_dp
    do i=1,4
        a(i,i) = 4.0_dp
        if (i < 4) then
            a(i+1,i) = -0.5_dp
            a(i,i+1) = -0.5_dp
        end if
    end do
    call pack_real_positive(a,2,.true.,lr)
    call pack_real_positive(a,2,.false.,ur)
    br(:,1) = matmul(a,xtrue)
    br(:,2) = matmul(a,[1.0_dp,2.0_dp,3.0_dp,4.0_dp])
    call solveh_banded(lr,br,xr,status,lower=.true.)
    if (status /= linalg_status_success .or. maxval(abs(xr(:,1)-xtrue)) > 2.0e-13_dp) &
        error stop 'real lower solveh'
    call solveh_banded(ur,br,xr,status)
    if (status /= linalg_status_success .or. maxval(abs(xr(:,1)-xtrue)) > 2.0e-13_dp) &
        error stop 'real upper solveh'
    call cholesky_banded(ur,cb_r,status)
    if (status /= linalg_status_success) error stop 'real cholesky factor'
    call cho_solve_banded(cb_r,br,xr,status)
    if (status /= linalg_status_success .or. maxval(abs(xr(:,1)-xtrue)) > 2.0e-13_dp) &
        error stop 'real cholesky reuse'

    ! Hermitian positive definite complex band generated from L L^H.
    lc = cmplx(0.0_dp,0.0_dp,dp)
    do i=1,4
        lc(i,i) = cmplx(2.0_dp+0.2_dp*i,0.0_dp,dp)
        if (i < 4) lc(i+1,i) = cmplx(0.17_dp,-0.12_dp,dp)
        if (i < 3) lc(i+2,i) = cmplx(-0.05_dp,0.16_dp,dp)
    end do
    ac = matmul(lc,transpose(conjg(lc)))
    call pack_complex_positive(ac,2,.true.,ab_l)
    call pack_complex_positive(ac,2,.false.,ab_u)
    bc = matmul(ac,xtruth)
    bc2(:,1) = bc
    bc2(:,2) = matmul(ac,[cmplx(1.0_dp,0.3_dp,dp), &
        cmplx(2.0_dp,0.6_dp,dp),cmplx(3.0_dp,-0.3_dp,dp),cmplx(4.0_dp,0.2_dp,dp)])
    call solveh_banded(ab_l,bc,xc,status,lower=.true.)
    if (status /= linalg_status_success .or. maxval(abs(xc-xtruth)) > 3.0e-13_dp) &
        error stop 'complex lower solveh'
    call solveh_banded(ab_u,bc2,xc2,status,lower=.false.)
    if (status /= linalg_status_success .or. maxval(abs(xc2(:,1)-xtruth)) > 3.0e-13_dp) &
        error stop 'complex upper multirhs solveh'
    call cholesky_banded(ab_l,cb_c,status,lower=.true.)
    if (status /= linalg_status_success) error stop 'complex cholesky lower factor'
    do j=1,4
        do i=j,min(4,j+2)
            if (abs(cb_c(1+i-j,j)-lc(i,j)) > 2.0e-13_dp) &
                error stop 'complex lower Cholesky values'
        end do
    end do
    call cho_solve_banded(cb_c,bc,xc,status,lower=.true.)
    if (status /= linalg_status_success .or. maxval(abs(xc-xtruth)) > 3.0e-13_dp) &
        error stop 'complex cholesky lower reuse'
    call cholesky_banded(ab_u,cb_c,status,lower=.false.)
    if (status /= linalg_status_success) error stop 'complex cholesky upper factor'
    do j=1,4
        do i=max(1,j-2),j
            if (abs(cb_c(3+i-j,j)-conjg(lc(j,i))) > 2.0e-13_dp) &
                error stop 'complex upper Cholesky values'
        end do
    end do
    call cho_solve_banded(cb_c,bc,xc,status,lower=.false.)
    if (status /= linalg_status_success .or. maxval(abs(xc-xtruth)) > 3.0e-13_dp) &
        error stop 'complex cholesky upper reuse'

    ! Failure tests: invalid bandwidth, singular pivot, and non-PD diagonal.
    call solve_banded(ab,3,1,b,x,status)
    if (status /= linalg_status_invalid_shape .or. any(x /= 0.0_dp)) &
        error stop 'bad band shape'
    ab = 0.0_dp
    call solve_banded(ab,1,1,b,x,status)
    if (status /= linalg_status_singular .or. any(x /= 0.0_dp)) &
        error stop 'singular band'
    lr = 0.0_dp
    lr(1,:) = [1.0_dp,1.0_dp,-1.0_dp,1.0_dp]
    call cholesky_banded(lr,cb_r,status,lower=.true.)
    if (status /= linalg_status_not_positive_definite) error stop 'not positive definite'
    ab_l(1,1) = cmplx(1.0_dp,1.0_dp,dp)
    call solveh_banded(ab_l,bc,xc,status,lower=.true.)
    if (status /= linalg_status_not_positive_definite) error stop 'non-Hermitian diagonal'

    ! Banded storage should scale linearly: a 50,000 x 50,000 tri-diagonal
    ! system needs O(n) storage, not a dense 20 GB coefficient matrix.
    n = 50000
    allocate(ab_big(3,n),b_big(n),x_big(n))
    ab_big = 0.0_dp
    ab_big(2,:) = 4.0_dp
    ab_big(1,2:n) = -0.3_dp
    ab_big(3,1:n-1) = -0.2_dp
    b_big = 3.5_dp
    b_big(1) = 3.7_dp
    b_big(n) = 3.8_dp
    call solve_banded(ab_big,1,1,b_big,x_big,status)
    if (status /= linalg_status_success) error stop 'large band solve status'
    if (maxval(abs(x_big-1.0_dp)) > 2.0e-13_dp) error stop 'large band solve accuracy'
    print '(A)', 'banded real/complex, pivoted/HPD, shape failures, and n=50000 tests passed'
contains
    subroutine pack_real(m,lower,upper,packed)
        real(dp), intent(in) :: m(:,:) !! dense test matrix
        integer, intent(in) :: lower !! subdiagonals
        integer, intent(in) :: upper !! superdiagonals
        real(dp), intent(out) :: packed(:,:) !! packed output
        integer :: ii,jj,nn
        nn = size(m,1)
        packed = 0.0_dp
        do jj=1,nn
            do ii=max(1,jj-upper),min(nn,jj+lower)
                packed(upper+1+ii-jj,jj) = m(ii,jj)
            end do
        end do
    end subroutine pack_real

    subroutine pack_complex(m,lower,upper,packed)
        complex(dp), intent(in) :: m(:,:) !! dense test matrix
        integer, intent(in) :: lower !! subdiagonals
        integer, intent(in) :: upper !! superdiagonals
        complex(dp), intent(out) :: packed(:,:) !! packed output
        integer :: ii,jj,nn
        nn = size(m,1)
        packed = cmplx(0.0_dp,0.0_dp,dp)
        do jj=1,nn
            do ii=max(1,jj-upper),min(nn,jj+lower)
                packed(upper+1+ii-jj,jj) = m(ii,jj)
            end do
        end do
    end subroutine pack_complex

    subroutine pack_real_positive(m,bw,lower,packed)
        real(dp), intent(in) :: m(:,:) !! symmetric test matrix
        integer, intent(in) :: bw !! half-bandwidth
        logical, intent(in) :: lower !! use lower or upper storage
        real(dp), intent(out) :: packed(:,:) !! packed output
        integer :: ii,jj,nn
        nn = size(m,1)
        packed = 0.0_dp
        do jj=1,nn
            if (lower) then
                do ii=jj,min(nn,jj+bw)
                    packed(1+ii-jj,jj) = m(ii,jj)
                end do
            else
                do ii=max(1,jj-bw),jj
                    packed(bw+1+ii-jj,jj) = m(ii,jj)
                end do
            end if
        end do
    end subroutine pack_real_positive

    subroutine pack_complex_positive(m,bw,lower,packed)
        complex(dp), intent(in) :: m(:,:) !! Hermitian test matrix
        integer, intent(in) :: bw !! half-bandwidth
        logical, intent(in) :: lower !! use lower or upper storage
        complex(dp), intent(out) :: packed(:,:) !! packed output
        integer :: ii,jj,nn
        nn = size(m,1)
        packed = cmplx(0.0_dp,0.0_dp,dp)
        do jj=1,nn
            if (lower) then
                do ii=jj,min(nn,jj+bw)
                    packed(1+ii-jj,jj) = m(ii,jj)
                end do
            else
                do ii=max(1,jj-bw),jj
                    packed(bw+1+ii-jj,jj) = m(ii,jj)
                end do
            end if
        end do
    end subroutine pack_complex_positive
end program test_linalg_banded
