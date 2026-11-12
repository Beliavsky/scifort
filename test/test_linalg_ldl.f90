! SPDX-License-Identifier: MIT
program test_linalg_ldl
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_ldl, only : ldl
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none

    real(dp) :: a(6,6), lu(6,6), d(6,6), nanval
    real(dp) :: small(2,2), lu2(2,2), d2(2,2)
    complex(dp) :: c(6,6), lcu(6,6), dc(6,6)
    integer :: perm(6), perm2(2), status, i, j, n, k
    complex(dp) :: cz(0,0), zeros(0,0), dzero(0,0)
    integer :: pz(0)

    ! A zero diagonal may require a 2x2 block, not a singularity error.
    small(1,:) = [0.0_dp,1.0_dp]
    small(2,:) = [1.0_dp,0.0_dp]
    call ldl(small,lu2,d2,perm2,status)
    if (status /= linalg_status_success) error stop 'indefinite 2x2 failed'
    if (abs(d2(1,2)) < 0.9_dp) error stop 'missing 2x2 pivot'
    call check_real(small)
    call check_real(small,lower=.false.)

    ! Singular diagonal (including an all-zero matrix) is factorable.
    small = 0.0_dp
    small(1,1) = 2.0_dp
    call check_real(small)
    small = 0.0_dp
    call check_real(small)

    ! Deliberate distant 2x2 swap with a preceding completed pivot.
    a = 0.0_dp
    do i = 1, 6
        a(i,i) = real(i,dp)
    end do
    a(1,1) = 12.0_dp
    a(2,2) = 0.0_dp
    a(5,5) = 0.0_dp
    a(2,5) = 7.0_dp
    a(5,2) = 7.0_dp
    a(3,2) = -0.4_dp
    a(2,3) = -0.4_dp
    call ldl(a,lu,d,perm,status)
    if (status /= linalg_status_success) error stop 'distant 2x2 swap failed'
    if (maxval(abs(d - diag_only(d))) < 1.0_dp) error stop 'missing off-diagonal pivot'
    call check_real(a)
    call check_real(a,lower=.false.)

    ! Complex Hermitian case with non-real off-diagonal entries.
    c = cmplx(a,0.0_dp,dp)
    c(2,5) = cmplx(3.0_dp,4.0_dp,dp)
    c(5,2) = conjg(c(2,5))
    c(4,1) = cmplx(0.3_dp,1.2_dp,dp)
    c(1,4) = conjg(c(4,1))
    call check_complex(c,.true.,.true.)
    call check_complex(c,.true.,.false.)

    ! Complex symmetric (not Hermitian) with imaginary diagonal entries.
    c = cmplx(a,0.0_dp,dp)
    c(2,5) = cmplx(2.0_dp,3.0_dp,dp)
    c(5,2) = c(2,5)
    c(1,1) = cmplx(12.0_dp,3.0_dp,dp)
    c(4,1) = cmplx(0.3_dp,1.2_dp,dp)
    c(1,4) = c(4,1)
    call check_complex(c,.false.,.true.)
    call check_complex(c,.false.,.false.)
    call ldl(c,lcu,dc,perm,status)
    if (status /= linalg_status_invalid_shape) error stop 'non-Hermitian diagonal accepted'

    ! Reproducible matrix family, varied dimensions and pivot patterns.
    do n = 1, 15
        do k = 1, 5
            block
                real(dp) :: r(n,n)
                complex(dp) :: z(n,n)
                do j = 1, n
                    do i = 1, n
                        r(i,j) = sin(real(17*i+29*j+7*k,dp)*0.371_dp)
                        z(i,j) = cmplx(r(i,j),cos(real(11*i+19*j+9*k,dp)*0.419_dp),dp)
                    end do
                end do
                r = 0.5_dp*(r+transpose(r))
                z = z+transpose(conjg(z))
                if (mod(k,2) == 0) then
                    do i = 1, n
                        r(i,i) = 0.0_dp
                        z(i,i) = cmplx(0.0_dp,0.0_dp,dp)
                    end do
                end if
                call check_real(r)
                call check_real(r,lower=.false.)
                call check_complex(z,.true.,.true.)
                call check_complex(z,.true.,.false.)
                z = cmplx(r,0.0_dp,dp)+cmplx(0.0_dp,1.0_dp,dp)* &
                    (z+transpose(z))*0.37_dp
                call check_complex(z,.false.,.true.)
                call check_complex(z,.false.,.false.)
            end block
        end do
    end do

    ! Test that input data in the unused triangle are never referenced.
    nanval = ieee_value(0.0_dp,ieee_quiet_nan)
    a = 0.0_dp
    do i = 1, 6
        a(i,i) = real(i,dp)
        do j = 1, i-1
            a(i,j) = 0.05_dp*real(i+j,dp)
            a(j,i) = nanval
        end do
    end do
    call check_real(a)
    a = transpose(a)
    call check_real(a,lower=.false.)
    small = 0.0_dp
    small(2,1) = nanval
    call ldl(small,lu2,d2,perm2,status)
    if (status /= linalg_status_no_convergence) error stop 'NaN not rejected'
    call ldl(small,lu2(1:1,:),d2,perm2,status)
    if (status /= linalg_status_invalid_shape) error stop 'invalid output shape'

    call ldl(cz,zeros,dzero,pz,status)
    if (status /= linalg_status_success) error stop 'empty array failed'
    print '(a)', 'ldl symmetric/Hermitian tests passed'

contains
    function diag_only(x) result(y)
        real(dp), intent(in) :: x(:,:) !! square matrix to diagonalize
        real(dp) :: y(size(x,1),size(x,2))
        integer :: ii
        y = 0.0_dp
        do ii = 1, size(x,1)
            y(ii,ii) = x(ii,ii)
        end do
    end function diag_only

    subroutine check_real(x,lower)
        real(dp), intent(in) :: x(:,:) !! input whose selected triangle is authoritative
        logical, intent(in), optional :: lower !! use lower part, otherwise upper
        real(dp) :: ll(size(x,1),size(x,1)), dd(size(x,1),size(x,1))
        real(dp) :: full(size(x,1),size(x,1)), t(size(x,1),size(x,1))
        integer :: p(size(x,1)), s, nn, ii, jj
        logical :: lo
        nn = size(x,1)
        lo = .true.
        if (present(lower)) lo = lower
        call ldl(x,ll,dd,p,s,lower=lo)
        if (s /= linalg_status_success) then
            print *, 'real size, status: ', nn, s
            error stop 'real LDL failed'
        end if
        full = 0.0_dp
        do jj = 1, nn
            do ii = 1, nn
                if (lo) then
                    full(ii,jj) = x(max(ii,jj),min(ii,jj))
                else
                    full(ii,jj) = x(min(ii,jj),max(ii,jj))
                end if
            end do
        end do
        t = matmul(ll,matmul(dd,transpose(ll)))
        if (maxval(abs(full-t)) > 2.0e-10_dp*max(1.0_dp,maxval(abs(full)))) &
            error stop 'real reconstruction failed'
        call check_shape(cmplx(ll,kind=dp),cmplx(dd,kind=dp),p,lo)
    end subroutine check_real

    subroutine check_complex(x,herm,lower)
        complex(dp), intent(in) :: x(:,:) !! input whose selected triangle is authoritative
        logical, intent(in) :: herm !! true for Hermitian; false for complex symmetric
        logical, intent(in) :: lower !! select lower or upper triangle
        complex(dp) :: ll(size(x,1),size(x,1)), dd(size(x,1),size(x,1))
        complex(dp) :: full(size(x,1),size(x,1)), t(size(x,1),size(x,1))
        integer :: p(size(x,1)), s, nn, ii, jj
        nn = size(x,1)
        call ldl(x,ll,dd,p,s,lower=lower,hermitian=herm)
        if (s /= linalg_status_success) then
            print *, 'complex size, hermitian, lower, status: ', nn, herm, lower, s
            error stop 'complex LDL failed'
        end if
        full = cmplx(0.0_dp,0.0_dp,dp)
        do jj = 1, nn
            do ii = 1, nn
                if (lower) then
                    if (ii >= jj) then
                        full(ii,jj) = x(ii,jj)
                    else
                        full(ii,jj) = x(jj,ii)
                        if (herm) full(ii,jj) = conjg(full(ii,jj))
                    end if
                else
                    if (ii <= jj) then
                        full(ii,jj) = x(ii,jj)
                    else
                        full(ii,jj) = x(jj,ii)
                        if (herm) full(ii,jj) = conjg(full(ii,jj))
                    end if
                end if
            end do
        end do
        if (herm) then
            t = matmul(ll,matmul(dd,transpose(conjg(ll))))
        else
            t = matmul(ll,matmul(dd,transpose(ll)))
        end if
        if (maxval(abs(full-t)) > 2.0e-10_dp*max(1.0_dp,maxval(abs(full)))) &
            error stop 'complex reconstruction failed'
        call check_shape(ll,dd,p,lower)
    end subroutine check_complex

    subroutine check_shape(ll,dd,p,lower)
        complex(dp), intent(in) :: ll(:,:) !! permuted outer factor
        complex(dp), intent(in) :: dd(:,:) !! block diagonal factor
        integer, intent(in) :: p(:) !! one-based row order for unit triangular form
        logical, intent(in) :: lower !! check lower or upper triangular form
        integer :: nn, ii, jj
        nn = size(ll,1)
        if (any(p < 1) .or. any(p > nn)) error stop 'invalid permutation entry'
        do ii = 1, nn
            if (count(p == ii) /= 1) error stop 'non-bijective permutation'
            if (abs(ll(p(ii),ii)-cmplx(1.0_dp,0.0_dp,dp)) > 1.0e-11_dp) &
                error stop 'outer diagonal not unit'
            do jj = 1, nn
                if (lower .and. jj > ii .and. abs(ll(p(ii),jj)) > 1.0e-11_dp) &
                    error stop 'outer factor not lower triangular'
                if (.not. lower .and. jj < ii .and. abs(ll(p(ii),jj)) > 1.0e-11_dp) &
                    error stop 'outer factor not upper triangular'
                if (abs(ii-jj) > 1 .and. abs(dd(ii,jj)) > 1.0e-11_dp) &
                    error stop 'non-block-diagonal D'
            end do
        end do
    end subroutine check_shape
end program test_linalg_ldl
