! SPDX-License-Identifier: MIT
program test_linalg_toeplitz
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence, linalg_status_singular
    use scifort_linalg_toeplitz, only : toeplitz, solve_toeplitz
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    implicit none
    integer :: n, k, status, trial
    logical :: used_dense
    real(dp) :: c0(0), b0(0), x0(0), t0(0,0)
    real(dp) :: cn(2), rn(2), bn(2), xn(2), t(2,2), nanval
    real(dp) :: zeros(2,0), no_rhs(2,0)
    complex(dp) :: cz(2), bz(2), xz(2), tz(2,2)

    call solve_toeplitz(c0,b0,x0,status)
    call success(status,'empty solve')
    call toeplitz(c0,t0,status)
    call success(status,'empty constructor')
    cn = [2.0_dp,0.125_dp]
    call solve_toeplitz(cn,zeros,no_rhs,status)
    call success(status,'zero RHS columns')

    cn = [0.0_dp,2.0_dp]
    rn = [99.0_dp,3.0_dp]
    call toeplitz(cn,t,status,rn)
    call success(status,'non-Hermitian 2x2 constructor')
    if (maxval(abs(t-reshape([0.0_dp,2.0_dp,3.0_dp,0.0_dp],[2,2]))) > 1.e-14_dp) &
        error stop 'r(1) must be ignored'
    bn = matmul(t,[4.0_dp,-3.0_dp])
    call solve_toeplitz(cn,bn,xn,status,rn,used_dense=used_dense)
    call success(status,'zero diagonal fallback')
    if (.not. used_dense) error stop 'expected dense fallback'
    if (maxval(abs(xn-[4.0_dp,-3.0_dp])) > 1.e-13_dp) error stop 'fallback result'

    cn = [0.0_dp,0.0_dp]
    call solve_toeplitz(cn,bn,xn,status)
    if (status /= linalg_status_singular .or. any(xn /= 0.0_dp)) &
        error stop 'zero Toeplitz matrix must be singular'
    ! An all-zero leading diagonal is not by itself proof of singularity.
    ! Here the whole Toeplitz matrix really is singular, and the pivoted
    ! fallback must report it while leaving the solution fully initialized.
    cn = [0.0_dp,1.0_dp]
    rn = [10.0_dp,0.0_dp]
    call solve_toeplitz(cn,bn,xn,status,rn,used_dense=used_dense)
    if (.not. used_dense .or. status /= linalg_status_singular) &
        error stop 'singular dense fallback status'
    if (any(xn /= 0.0_dp)) error stop 'singular dense fallback output'
    cn = [2.0_dp,1.0_dp]
    call solve_toeplitz(cn,bn,xn,status,r=[1.0_dp])
    if (status /= linalg_status_invalid_shape .or. any(xn /= 0.0_dp)) &
        error stop 'invalid first row length'
    call toeplitz(cn,t0,status)
    if (status /= linalg_status_invalid_shape) error stop 'bad Toeplitz constructor shape'

    nanval = ieee_value(0.0_dp,ieee_quiet_nan)
    cn = [2.0_dp,0.5_dp]
    bn = [1.0_dp,nanval]
    call solve_toeplitz(cn,bn,xn,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite RHS accepted'
    bn = [1.0_dp,2.0_dp]
    call solve_toeplitz(cn,bn,xn,status,r=[nanval,0.5_dp])
    if (status /= linalg_status_no_convergence) error stop 'nonfinite first row accepted'

    cz = [cmplx(3.0_dp,2.0_dp,dp),cmplx(0.2_dp,-0.6_dp,dp)]
    call toeplitz(cz,tz,status)
    call success(status,'default complex row')
    if (abs(tz(1,2)-conjg(cz(2))) > 0.0_dp) error stop 'complex row conjugation'
    if (abs(tz(1,1)-cz(1)) > 0.0_dp) error stop 'complex diagonal altered'
    bz = matmul(tz,[cmplx(0.5_dp,1.1_dp,dp),cmplx(-2.0_dp,0.2_dp,dp)])
    call solve_toeplitz(cz,bz,xz,status,used_dense=used_dense)
    call success(status,'complex diagonal default row')
    if (used_dense) error stop 'unexpected complex default fallback'
    if (maxval(abs(xz-[cmplx(0.5_dp,1.1_dp,dp),cmplx(-2.0_dp,0.2_dp,dp)])) &
        > 1.0e-13_dp) error stop 'complex 2x2 result'

    ! Zero leading complex diagonal, but nonsingular final Toeplitz matrix.
    cz = [cmplx(0.0_dp,0.0_dp,dp),cmplx(2.0_dp,1.0_dp,dp)]
    call toeplitz(cz,tz,status,[cmplx(999.0_dp,0.0_dp,dp), &
        cmplx(3.0_dp,-1.0_dp,dp)])
    call success(status,'complex zero pivot constructor')
    bz = matmul(tz,[cmplx(1.0_dp,2.0_dp,dp),cmplx(3.0_dp,-1.0_dp,dp)])
    call solve_toeplitz(cz,bz,xz,status, &
        r=[cmplx(999.0_dp,0.0_dp,dp),cmplx(3.0_dp,-1.0_dp,dp)],used_dense=used_dense)
    call success(status,'complex zero pivot solve')
    if (.not. used_dense) error stop 'expected complex dense fallback'
    if (maxval(abs(xz-[cmplx(1.0_dp,2.0_dp,dp),cmplx(3.0_dp,-1.0_dp,dp)])) &
        > 1.0e-13_dp) error stop 'complex fallback result'

    ! General unsymmetric real and complex systems, including repeated RHS.
    do n = 1, 40
        do trial = 1, 4
            call exercise_real(n,trial)
            call exercise_complex(n,trial)
        end do
    end do

    ! Rectangular Toeplitz construction, with an intentionally ignored r(1).
    block
        real(dp) :: c(3), r(5), a(3,5)
        c = [1.0_dp,2.0_dp,3.0_dp]
        r = [-100.0_dp,4.0_dp,5.0_dp,6.0_dp,7.0_dp]
        call toeplitz(c,a,status,r)
        call success(status,'rectangular Toeplitz')
        if (any(abs(a(1,:)-[1.0_dp,4.0_dp,5.0_dp,6.0_dp,7.0_dp]) &
            > 1.0e-14_dp)) error stop 'rectangular first row'
        if (any(abs(a(:,1)-c) > 1.0e-14_dp)) error stop 'rectangular first column'
    end block

    ! Solver must remain scale-invariant over a wide physical-unit range.
    do k = -80, 80, 80
        block
            real(dp) :: c(5), r(5), tt(5,5), x_ref(5), b(5), x(5), factor
            factor = 10.0_dp**k
            c = factor*[3.0_dp,0.1_dp,-0.25_dp,0.08_dp,0.12_dp]
            r = factor*[55.0_dp,-0.13_dp,0.2_dp,-0.05_dp,0.07_dp]
            call toeplitz(c,tt,status,r)
            call success(status,'scaled Toeplitz construction')
            x_ref = [0.3_dp,-0.2_dp,0.5_dp,1.1_dp,-0.4_dp]
            b = matmul(tt,x_ref)
            call solve_toeplitz(c,b,x,status,r)
            call success(status,'scaled Toeplitz solve')
            if (maxval(abs(x-x_ref)) > 2.0e-12_dp) error stop 'scaled solution'
        end block
    end do

    print '(a)', 'test_linalg_toeplitz: PASS'

contains

    subroutine success(code, label)
        integer, intent(in) :: code !! solver status to check
        character(len=*), intent(in) :: label !! case identification for a failure
        if (code /= linalg_status_success) then
            print '(a,1x,a,1x,i0)', 'FAILED',label,code
            error stop 'Toeplitz operation failed'
        end if
    end subroutine success

    subroutine exercise_real(nn, seed)
        integer, intent(in) :: nn !! system order, positive
        integer, intent(in) :: seed !! deterministic case index
        real(dp) :: c(nn), r(nn), a(nn,nn), y(nn,3), rhs(nn,3), xx(nn,3)
        real(dp) :: vv(nn), bb(nn), sol(nn), symmetry(nn,nn)
        integer :: i, q, code
        logical :: dense

        c(1) = real(nn+4,dp)
        r(1) = -10.0_dp
        do i = 2, nn
            c(i) = 0.2_dp*sin(real(11*i+3*seed,dp))
            r(i) = 0.3_dp*cos(real(7*i-4*seed,dp))
        end do
        call toeplitz(c,a,code,r)
        call success(code,'real constructor')
        do q = 1, 3
            do i = 1, nn
                y(i,q) = sin(real(3*i+q+seed,dp))
            end do
        end do
        rhs = matmul(a,y)
        call solve_toeplitz(c,rhs,xx,code,r,used_dense=dense)
        call success(code,'real multiple RHS')
        if (dense) error stop 'unexpected real dense fallback'
        if (maxval(abs(y-xx)) > 5.0e-12_dp) error stop 'real multiple RHS mismatch'
        call solve_toeplitz(c,rhs(:,2),sol,code,r)
        call success(code,'real vector RHS')
        if (maxval(abs(sol-y(:,2))) > 5.0e-12_dp) error stop 'real vector RHS mismatch'

        ! SciPy's default row is the conjugate of c, which is c for real input.
        call toeplitz(c,symmetry,code)
        call success(code,'symmetric constructor')
        bb = matmul(symmetry,y(:,1))
        call solve_toeplitz(c,bb,vv,code)
        call success(code,'symmetric default solve')
        if (maxval(abs(vv-y(:,1))) > 5.0e-12_dp) error stop 'default row real solve'

        block
            complex(dp) :: zrhs(nn,3), zsol(nn,3)
            zrhs = cmplx(rhs,matmul(a,y)*0.2_dp,dp)
            call solve_toeplitz(c,zrhs,zsol,code,r)
            call success(code,'real coefficients complex RHS')
            if (maxval(abs(zsol-cmplx(y,0.2_dp*y,dp))) > 5.0e-12_dp) &
                error stop 'real coefficients complex RHS mismatch'
        end block
    end subroutine exercise_real

    subroutine exercise_complex(nn, seed)
        integer, intent(in) :: nn !! system order, positive
        integer, intent(in) :: seed !! deterministic case index
        complex(dp) :: c(nn), r(nn), a(nn,nn), y(nn,2), rhs(nn,2), xx(nn,2)
        complex(dp) :: bb(nn), sol(nn), hh(nn,nn)
        real(dp) :: real_rhs(nn,2)
        integer :: i, q, code
        logical :: dense

        c(1) = cmplx(real(nn+4,dp),0.7_dp,dp)
        r(1) = cmplx(8.0_dp,-22.0_dp,dp)
        do i = 2, nn
            c(i) = cmplx(0.25_dp*sin(real(11*i+seed,dp)), &
                0.15_dp*cos(real(13*i+seed,dp)),dp)
            r(i) = cmplx(0.14_dp*cos(real(7*i+seed,dp)), &
                0.17_dp*sin(real(9*i+seed,dp)),dp)
        end do
        call toeplitz(c,a,code,r)
        call success(code,'complex constructor')
        do q = 1, 2
            do i = 1, nn
                y(i,q) = cmplx(sin(real(3*i+q+seed,dp)), &
                    cos(real(5*i-2*q+seed,dp)),dp)
            end do
        end do
        rhs = matmul(a,y)
        call solve_toeplitz(c,rhs,xx,code,r,used_dense=dense)
        call success(code,'complex multiple RHS')
        if (dense) error stop 'unexpected complex dense fallback'
        if (maxval(abs(y-xx)) > 5.0e-12_dp) error stop 'complex multiple RHS mismatch'
        call solve_toeplitz(c,rhs(:,2),sol,code,r)
        call success(code,'complex vector RHS')
        if (maxval(abs(sol-y(:,2))) > 5.0e-12_dp) error stop 'complex vector RHS mismatch'

        call toeplitz(c,hh,code)
        call success(code,'conjugate-row constructor')
        bb = matmul(hh,y(:,1))
        call solve_toeplitz(c,bb,sol,code)
        call success(code,'conjugate-row solve')
        if (maxval(abs(sol-y(:,1))) > 5.0e-12_dp) error stop 'default row complex solve'

        real_rhs = real(rhs,dp)
        call solve_toeplitz(c,real_rhs,xx,code,r)
        call success(code,'complex coefficients real RHS')
        if (maxval(abs(matmul(a,xx)-cmplx(real_rhs,kind=dp))) > 5.0e-12_dp) &
            error stop 'complex coefficients real RHS mismatch'
    end subroutine exercise_complex
end program test_linalg_toeplitz
