! SPDX-License-Identifier: MIT
program test_linalg_indefinite_solve
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_singular, &
        linalg_status_invalid_shape, linalg_status_no_convergence
    use scifort_linalg_ldl, only : ldl
    use scifort_linalg_indefinite_solve, only : ldl_solve, solve_indefinite
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none

    real(dp) :: a(2,2), lu(2,2), d(2,2), b(2), x(2), xx(2,2), aa(2,2), nanval
    complex(dp) :: c(2,2), cl(2,2), cd(2,2), cb(2), cx(2)
    real(dp) :: empty(0,0), emptyb(0,3), emptyx(0,3), none(2,0), nonex(2,0)
    integer :: perm(2), status, n, trial, i, j

    ! Pure 2-by-2 pivot, including exact zeros on both diagonal entries.
    a(1,:) = [0.0_dp, 3.0_dp]
    a(2,:) = [3.0_dp, 0.0_dp]
    b = [6.0_dp, 9.0_dp]
    call solve_indefinite(a,b,x,status)
    call assert_success(status,'2x2 real indefinite')
    if (maxval(abs(x-[3.0_dp,2.0_dp])) > 1.0e-13_dp) error stop '2x2 real solution'
    call solve_indefinite(a,b,x,status,lower=.false.)
    call assert_success(status,'2x2 real upper')
    if (maxval(abs(x-[3.0_dp,2.0_dp])) > 1.0e-13_dp) error stop '2x2 upper solution'

    ! Off-diagonal complex 2x2 pivot, conjugation vs ordinary transpose.
    c = cmplx(0.0_dp,0.0_dp,dp)
    c(1,2) = cmplx(2.0_dp,3.0_dp,dp)
    c(2,1) = conjg(c(1,2))
    call check_complex(c,.true.,.true.)
    call check_complex(c,.true.,.false.)
    c(2,1) = c(1,2)
    call check_complex(c,.false.,.true.)
    call check_complex(c,.false.,.false.)

    ! Multiple completed pivots followed by a remote symmetric swap.
    block
        real(dp) :: swapped(6,6)
        integer :: k
        swapped = 0.0_dp
        do k = 1, 6
            swapped(k,k) = real(k,dp)
        end do
        swapped(1,1) = 12.0_dp
        swapped(2,2) = 0.0_dp
        swapped(5,5) = 0.0_dp
        swapped(2,5) = 7.0_dp
        swapped(5,2) = 7.0_dp
        swapped(2,3) = -0.4_dp
        swapped(3,2) = -0.4_dp
        call check_real(swapped,.true.)
        call check_real(swapped,.false.)
    end block

    ! Permutation, previous blocks, and repeated RHS solve.
    do n = 1, 24
        do trial = 1, 6
            block
                real(dp) :: m(n,n)
                complex(dp) :: z(n,n)
                do j = 1, n
                    do i = 1, n
                        m(i,j) = sin(0.197_dp*real(9*i+17*j+3*trial,dp))
                        z(i,j) = cmplx(m(i,j), &
                            cos(0.127_dp*real(3*i+11*j+19*trial,dp)),dp)
                    end do
                end do
                m = 0.5_dp*(m+transpose(m))
                z = z+transpose(conjg(z))
                do i = 1, n
                    if (mod(i+trial,2) == 0) then
                        m(i,i) = m(i,i)+real(n+2,dp)
                        z(i,i) = z(i,i)+real(2*n+4,dp)
                    else
                        m(i,i) = m(i,i)-real(n+2,dp)
                        z(i,i) = z(i,i)-real(2*n+4,dp)
                    end if
                end do
                call check_real(m,.true.)
                call check_real(m,.false.)
                call check_complex(z,.true.,.true.)
                call check_complex(z,.true.,.false.)
                z = cmplx(m,kind=dp)+cmplx(0.0_dp,1.0_dp,dp)* &
                    (z+transpose(z))*0.03_dp
                call check_complex(z,.false.,.true.)
                call check_complex(z,.false.,.false.)
            end block
        end do
    end do

    ! Missing or duplicate permutation entries must be rejected.
    call ldl(a,lu,d,perm,status)
    call assert_success(status,'2x2 factor')
    perm = [1,1]
    call ldl_solve(lu,d,perm,b,x,status)
    if (status /= linalg_status_invalid_shape .or. any(x /= 0.0_dp)) &
        error stop 'duplicate permutation accepted'
    perm = [0,2]
    call ldl_solve(lu,d,perm,b,x,status)
    if (status /= linalg_status_invalid_shape) error stop 'bad permutation accepted'

    ! A singular 1-by-1 pivot is not a solvable linear system.
    a = 0.0_dp
    a(1,1) = 2.0_dp
    call solve_indefinite(a,b,x,status)
    if (status /= linalg_status_singular .or. any(x /= 0.0_dp)) &
        error stop 'singular real system accepted'
    call ldl(a,lu,d,perm,status)
    call assert_success(status,'factorable singular real')
    call ldl_solve(lu,d,perm,b,x,status)
    if (status /= linalg_status_singular) error stop 'singular LDL solve accepted'

    ! Malformed D and numerically singular 2-by-2 D are not inverted.
    lu = 0.0_dp
    lu(1,1) = 1.0_dp
    lu(2,2) = 1.0_dp
    perm = [1,2]
    d = 1.0_dp
    call ldl_solve(lu,d,perm,b,x,status)
    if (status /= linalg_status_singular) error stop 'rank-one 2x2 D accepted'
    d = 0.0_dp
    d(1,2) = 1.0_dp
    call ldl_solve(lu,d,perm,b,x,status)
    if (status /= linalg_status_invalid_shape) error stop 'asymmetric D accepted'

    ! Real and complex argument shape checks.
    aa = 0.0_dp
    call solve_indefinite(aa,[1.0_dp],x,status)
    if (status /= linalg_status_invalid_shape) error stop 'bad real shape accepted'
    call solve_indefinite(aa,reshape([1.0_dp,2.0_dp],[1,2]),xx,status)
    if (status /= linalg_status_invalid_shape) error stop 'bad RHS matrix shape accepted'

    ! Hermitian rejects imaginary diagonals, complex symmetric does not.
    c = cmplx(0.0_dp,0.0_dp,dp)
    c(1,1) = cmplx(1.0_dp,1.0_dp,dp)
    c(2,2) = cmplx(1.0_dp,0.0_dp,dp)
    cb = cmplx([1.0_dp,2.0_dp],kind=dp)
    call solve_indefinite(c,cb,cx,status)
    if (status /= linalg_status_invalid_shape) error stop 'Hermitian imaginary diagonal accepted'
    call solve_indefinite(c,cb,cx,status,hermitian=.false.)
    call assert_success(status,'complex symmetric imaginary diagonal')
    if (maxval(abs(matmul(c,cx)-cb)) > 1.0e-12_dp) &
        error stop 'complex symmetric diagonal solution failed'

    nanval = ieee_value(0.0_dp,ieee_quiet_nan)
    c = cmplx(0.0_dp,0.0_dp,dp)
    c(1,1) = cmplx(2.0_dp,0.0_dp,dp)
    c(2,2) = cmplx(3.0_dp,0.0_dp,dp)
    call ldl(c,cl,cd,perm,status)
    call assert_success(status,'complex finite factors')
    cb(1) = cmplx(nanval,0.0_dp,dp)
    call ldl_solve(cl,cd,perm,cb,cx,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite RHS accepted'

    ! Empty system and zero RHS columns are both legal.
    call solve_indefinite(empty,emptyb,emptyx,status)
    call assert_success(status,'n=0 real')
    aa = 0.0_dp
    aa(1,1) = 2.0_dp
    aa(2,2) = 3.0_dp
    call solve_indefinite(aa,none,nonex,status)
    call assert_success(status,'nrhs=0 real')

    print *, 'indefinite LDL solve native tests passed'
contains

    subroutine assert_success(s, where)
        integer, intent(in) :: s !! observed status
        character(*), intent(in) :: where !! case name
        if (s /= linalg_status_success) then
            print *, where, s
            error stop 'unexpected non-success status'
        end if
    end subroutine assert_success

    subroutine check_real(m,lower)
        real(dp), intent(in) :: m(:, :) !! nonsingular symmetric coefficient matrix
        logical, intent(in) :: lower !! selected input triangle
        real(dp) :: fac(size(m,1),size(m,1)), blockd(size(m,1),size(m,1))
        real(dp) :: rhs(size(m,1),3), truex(size(m,1),3)
        real(dp) :: solution(size(m,1),3), via_factor(size(m,1),3)
        real(dp) :: one(size(m,1)), one_x(size(m,1)), altered(size(m,1),size(m,1))
        integer :: p(size(m,1)), s, nn, ii, jj
        nn = size(m,1)
        do jj = 1, 3
            do ii = 1, nn
                truex(ii,jj) = sin(real(3*ii+17*jj,dp)*0.128_dp)
            end do
        end do
        rhs = matmul(m,truex)
        altered = m
        do jj = 1, nn
            do ii = 1, nn
                if ((lower .and. ii < jj) .or. (.not. lower .and. ii > jj)) &
                    altered(ii,jj) = -12.0_dp
            end do
        end do
        call ldl(altered,fac,blockd,p,s,lower=lower)
        call assert_success(s,'real ldl')
        call ldl_solve(fac,blockd,p,rhs,via_factor,s,lower=lower)
        call assert_success(s,'real factor solve')
        call solve_indefinite(altered,rhs,solution,s,lower=lower)
        call assert_success(s,'real direct solve')
        if (maxval(abs(solution-truex)) > 2.0e-10_dp) error stop 'real direct residual'
        if (maxval(abs(via_factor-truex)) > 2.0e-10_dp) error stop 'real factor residual'
        one = rhs(:,1)
        call ldl_solve(fac,blockd,p,one,one_x,s,lower=lower)
        call assert_success(s,'real vector factor solve')
        if (maxval(abs(one_x-truex(:,1))) > 2.0e-10_dp) error stop 'real vector factor residual'
        call solve_indefinite(altered,one,one_x,s,lower=lower)
        call assert_success(s,'real vector direct solve')
        if (maxval(abs(one_x-truex(:,1))) > 2.0e-10_dp) error stop 'real vector direct residual'
    end subroutine check_real

    subroutine check_complex(m,hermitian,lower)
        complex(dp), intent(in) :: m(:, :) !! nonsingular symmetric or Hermitian matrix
        logical, intent(in) :: hermitian !! conjugate transpose or transpose
        logical, intent(in) :: lower !! selected input triangle
        complex(dp) :: fac(size(m,1),size(m,1)), blockd(size(m,1),size(m,1))
        complex(dp) :: rhs(size(m,1),3), truex(size(m,1),3)
        complex(dp) :: solution(size(m,1),3), via_factor(size(m,1),3)
        complex(dp) :: one(size(m,1)), one_x(size(m,1)), altered(size(m,1),size(m,1))
        integer :: p(size(m,1)), s, nn, ii, jj
        nn = size(m,1)
        do jj = 1, 3
            do ii = 1, nn
                truex(ii,jj) = cmplx(sin(real(ii+17*jj,dp)*0.13_dp), &
                    cos(real(3*ii+7*jj,dp)*0.251_dp),dp)
            end do
        end do
        rhs = matmul(m,truex)
        altered = m
        do jj = 1, nn
            do ii = 1, nn
                if ((lower .and. ii < jj) .or. (.not. lower .and. ii > jj)) &
                    altered(ii,jj) = cmplx(-12.0_dp,7.0_dp,dp)
            end do
        end do
        call ldl(altered,fac,blockd,p,s,lower=lower,hermitian=hermitian)
        call assert_success(s,'complex ldl')
        call ldl_solve(fac,blockd,p,rhs,via_factor,s,lower=lower,hermitian=hermitian)
        call assert_success(s,'complex factor solve')
        call solve_indefinite(altered,rhs,solution,s,lower=lower,hermitian=hermitian)
        call assert_success(s,'complex direct solve')
        if (maxval(abs(solution-truex)) > 2.0e-10_dp) error stop 'complex direct residual'
        if (maxval(abs(via_factor-truex)) > 2.0e-10_dp) error stop 'complex factor residual'
        one = rhs(:,1)
        call ldl_solve(fac,blockd,p,one,one_x,s,lower=lower,hermitian=hermitian)
        call assert_success(s,'complex vector factor solve')
        if (maxval(abs(one_x-truex(:,1))) > 2.0e-10_dp) error stop 'complex vector factor residual'
        call solve_indefinite(altered,one,one_x,s,lower=lower,hermitian=hermitian)
        call assert_success(s,'complex vector direct solve')
        if (maxval(abs(one_x-truex(:,1))) > 2.0e-10_dp) error stop 'complex vector direct residual'
    end subroutine check_complex
end program test_linalg_indefinite_solve
