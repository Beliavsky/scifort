! SPDX-License-Identifier: MIT
program test_linalg_interpolative_svd_operator
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_interpolative_svd_operator, only : randomized_interpolative_svd_operator, &
        randomized_interpolative_svd_operator_complex
    implicit none
    real(dp), allocatable :: a(:, :), u(:, :), s(:), v(:, :), ud(:, :), sd(:), vd(:, :)
    real(dp), allocatable :: a_recon(:, :), us(:, :)
    complex(dp), allocatable :: z(:, :), uc(:, :), vc(:, :), z_recon(:, :), ucs(:, :)
    complex(dp), allocatable :: ucd(:, :), vcd(:, :)
    integer :: m, n, i, j, k, ka, variant, status, checks, fwd, adj
    logical :: fail_forward, fail_adjoint, emit_nan
    real(dp) :: err, scale, eps, est
    integer(int64), parameter :: seed = 20261110_int64

    checks = 0
    fail_forward = .false.
    fail_adjoint = .false.
    emit_nan = .false.

    do m = 0, 11
        do n = 0, 13
            allocate(a(m,n),z(m,n),a_recon(m,n),z_recon(m,n))
            do variant = 1, 6
                do j = 1, n
                    do i = 1, m
                        a(i,j) = cos(real(5*i+2*j+variant,dp)*0.23_dp) + &
                            sin(real(2*i-3*j+variant,dp)*0.39_dp)
                        z(i,j) = cmplx(a(i,j),sin(real(6*i+5*j+variant,dp)*0.27_dp),dp)
                    end do
                end do
                select case (variant)
                case (2)
                    do j = 2, n
                        a(:,j) = a(:,1)*real(j,dp)
                        z(:,j) = z(:,1)*cmplx(real(j,dp),real(j,dp)*0.2_dp,dp)
                    end do
                case (3)
                    a = 0.0_dp
                    z = cmplx(0.0_dp,0.0_dp,dp)
                case (4)
                    a = a*1.0e-80_dp
                    z = z*1.0e-80_dp
                case (5)
                    a = a*1.0e80_dp
                    z = z*1.0e80_dp
                case (6)
                    if (n > 2) then
                        do j = 3, n
                            a(:,j) = a(:,1)*real(j,dp) + a(:,2)/real(j,dp)
                            z(:,j) = z(:,1)*real(j,dp) + z(:,2)/real(j,dp)
                        end do
                    end if
                end select
                do ka = 1, 4
                    select case(ka)
                    case(1)
                        k = 0
                    case(2)
                        k = min(m,n)
                    case(3)
                        k = min(1,min(m,n))
                    case(4)
                        k = min(3,min(m,n))
                    end select
                    allocate(u(m,k),s(k),v(n,k),uc(m,k),vc(n,k),us(m,k),ucs(m,k))
                    call randomized_interpolative_svd_operator(m,n,k,real_mv,real_rmv,u,s,v,status, &
                        seed=seed,oversampling=5)
                    call check(status == linalg_status_success,'real fixed status')
                    call check(all(s >= 0.0_dp),'real singular values positive')
                    call check_sorted(s)
                    call test_real_factors(k)
                    call randomized_interpolative_svd_operator_complex(m,n,k,complex_mv,complex_rmv, &
                        uc,s,vc,status,seed=seed,oversampling=5)
                    call check(status == linalg_status_success,'complex fixed status')
                    call check(all(s >= 0.0_dp),'complex singular values positive')
                    call check_sorted(s)
                    call test_complex_factors(k)
                    deallocate(u,s,v,uc,vc,us,ucs)
                end do
                if (variant == 2 .or. variant == 3 .or. variant == 6) then
                    eps = 1.0e-5_dp
                    call randomized_interpolative_svd_operator(m,n,eps,real_mv,real_rmv,k, &
                        ud,sd,vd,status,seed=seed,verification_probes=20,estimated_error=est)
                    call check(status == linalg_status_success,'real adaptive status')
                    call check(size(ud,1) == m .and. size(ud,2) == k,'real adaptive shape')
                    call check(est <= eps,'real adaptive estimate')
                    call check_true_real(ud,sd,vd,k,1.0e-6_dp)
                    call randomized_interpolative_svd_operator_complex(m,n,eps,complex_mv, &
                        complex_rmv,k,ucd,sd,vcd,status,seed=seed,verification_probes=20, &
                        estimated_error=est)
                    call check(status == linalg_status_success,'complex adaptive status')
                    call check(size(ucd,1) == m .and. size(ucd,2) == k,'complex adaptive shape')
                    call check(est <= eps,'complex adaptive estimate')
                    call check_true_complex(ucd,sd,vcd,k,1.0e-6_dp)
                end if
            end do
            deallocate(a,z,a_recon,z_recon)
        end do
    end do

    allocate(a(4,3),z(4,3),u(4,2),s(2),v(3,2),uc(4,2),vc(3,2))
    a = 1.0_dp
    z = cmplx(1.0_dp,0.2_dp,dp)
    call randomized_interpolative_svd_operator(4,3,2,real_mv,real_rmv,u,s,v,status, &
        oversampling=-1)
    call check(status == linalg_status_invalid_shape,'negative oversampling')
    call randomized_interpolative_svd_operator(4,3,2,real_mv,real_rmv,u(:3,:),s,v,status)
    call check(status == linalg_status_invalid_shape,'incorrect factor shape')
    fail_forward = .true.
    call randomized_interpolative_svd_operator(4,3,2,real_mv,real_rmv,u,s,v,status)
    call check(status == linalg_status_no_convergence,'failed real forward callback')
    fail_forward = .false.
    fail_adjoint = .true.
    call randomized_interpolative_svd_operator_complex(4,3,2,complex_mv,complex_rmv,uc,s,vc,status)
    call check(status == linalg_status_no_convergence,'failed complex adjoint callback')
    fail_adjoint = .false.
    emit_nan = .true.
    call randomized_interpolative_svd_operator(4,3,2,real_mv,real_rmv,u,s,v,status)
    call check(status == linalg_status_no_convergence,'nonfinite real callback')
    emit_nan = .false.
    call randomized_interpolative_svd_operator(4,3,-1,real_mv,real_rmv,u,s,v,status)
    call check(status == linalg_status_invalid_shape,'negative rank')
    call randomized_interpolative_svd_operator(4,3,0.0_dp,real_mv,real_rmv,k, &
        ud,sd,vd,status,seed=seed)
    call check(status == linalg_status_success,'adaptive zero tolerance')
    call randomized_interpolative_svd_operator(4,3,-0.1_dp,real_mv,real_rmv,k, &
        ud,sd,vd,status,seed=seed)
    call check(status == linalg_status_invalid_shape,'invalid tolerance')
    call check(k == 0 .and. size(sd) == 0,'error leaves empty results')
    call randomized_interpolative_svd_operator(4,3,0.1_dp,real_mv,real_rmv,k, &
        ud,sd,vd,status,verification_probes=0)
    call check(status == linalg_status_invalid_shape,'invalid probe count')
    print '(a,i0)', 'Matrix-free interpolative SVD native checks passed: ', checks

contains
    subroutine check(ok, message)
        logical, intent(in) :: ok !! required condition
        character(*), intent(in) :: message !! diagnostic description
        if (.not. ok) then
            print '(a)', message
            error stop 'matrix-free SVD check failed'
        end if
        checks = checks + 1
    end subroutine check
    subroutine check_sorted(vals)
        real(dp), intent(in) :: vals(:) !! descending singular values
        integer :: j
        do j = 2, size(vals)
            call check(vals(j) <= vals(j-1)*(1.0_dp+1.0e-12_dp) + tiny(1.0_dp), &
                'singular values descending')
        end do
    end subroutine check_sorted
    subroutine test_real_factors(rank)
        integer, intent(in) :: rank !! requested factor count
        real(dp) :: denom
        integer :: j
        us = u
        do j = 1, rank
            us(:,j) = us(:,j)*s(j)
        end do
        a_recon = matmul(us,transpose(v))
        denom = 0.0_dp
        if (size(a) > 0) denom = maxval(abs(a))
        if (denom > 0.0_dp) then
            err = maxval(abs(a-a_recon))/denom
            if (rank == min(m,n)) call check(err < 1.0e-8_dp,'real full reconstruction')
            if (variant == 2 .and. rank >= 1) call check(err < 1.0e-8_dp,'real low-rank reconstruction')
        end if
        if (rank > 0) then
            err = maxval(abs(matmul(transpose(u),u)-identity_real(rank)))
            call check(err < 1.0e-8_dp,'real orthogonality U')
            err = maxval(abs(matmul(transpose(v),v)-identity_real(rank)))
            call check(err < 1.0e-8_dp,'real orthogonality V')
        end if
    end subroutine test_real_factors
    subroutine test_complex_factors(rank)
        integer, intent(in) :: rank !! requested factor count
        real(dp) :: denom
        integer :: j
        ucs = uc
        do j = 1, rank
            ucs(:,j) = ucs(:,j)*s(j)
        end do
        z_recon = matmul(ucs,transpose(conjg(vc)))
        denom = 0.0_dp
        if (size(z) > 0) denom = maxval(abs(z))
        if (denom > 0.0_dp) then
            err = maxval(abs(z-z_recon))/denom
            if (rank == min(m,n)) call check(err < 1.0e-8_dp,'complex full reconstruction')
            if (variant == 2 .and. rank >= 1) call check(err < 1.0e-8_dp,'complex low-rank reconstruction')
        end if
        if (rank > 0) then
            err = maxval(abs(matmul(transpose(conjg(uc)),uc)-identity_complex(rank)))
            call check(err < 1.0e-8_dp,'complex orthogonality U')
            err = maxval(abs(matmul(transpose(conjg(vc)),vc)-identity_complex(rank)))
            call check(err < 1.0e-8_dp,'complex orthogonality V')
        end if
    end subroutine test_complex_factors
    subroutine check_true_real(left,values,right,rank,limit)
        real(dp), intent(in) :: left(:, :) !! m-by-k left vectors
        real(dp), intent(in) :: values(:) !! singular values
        real(dp), intent(in) :: right(:, :) !! n-by-k right vectors
        integer, intent(in) :: rank !! numerical rank
        real(dp), intent(in) :: limit !! relative error tolerance
        real(dp), allocatable :: scaled(:, :)
        real(dp) :: denom
        integer :: j
        scaled = left
        do j = 1, rank
            scaled(:,j) = scaled(:,j)*values(j)
        end do
        a_recon = matmul(scaled,transpose(right))
        denom = 0.0_dp
        if (size(a) > 0) denom = maxval(abs(a))
        if (denom > 0.0_dp) call check(maxval(abs(a-a_recon))/denom <= limit, &
            'real actual residual')
    end subroutine check_true_real
    subroutine check_true_complex(left,values,right,rank,limit)
        complex(dp), intent(in) :: left(:, :) !! m-by-k left vectors
        real(dp), intent(in) :: values(:) !! singular values
        complex(dp), intent(in) :: right(:, :) !! n-by-k right vectors
        integer, intent(in) :: rank !! numerical rank
        real(dp), intent(in) :: limit !! relative error tolerance
        complex(dp), allocatable :: scaled(:, :)
        real(dp) :: denom
        integer :: j
        scaled = left
        do j = 1, rank
            scaled(:,j) = scaled(:,j)*values(j)
        end do
        z_recon = matmul(scaled,transpose(conjg(right)))
        denom = 0.0_dp
        if (size(z) > 0) denom = maxval(abs(z))
        if (denom > 0.0_dp) call check(maxval(abs(z-z_recon))/denom <= limit, &
            'complex actual residual')
    end subroutine check_true_complex
    function identity_real(d) result(out)
        integer, intent(in) :: d !! square order
        real(dp) :: out(d,d)
        integer :: j
        out = 0.0_dp
        do j = 1,d
            out(j,j) = 1.0_dp
        end do
    end function identity_real
    function identity_complex(d) result(out)
        integer, intent(in) :: d !! square order
        complex(dp) :: out(d,d)
        integer :: j
        out = cmplx(0.0_dp,0.0_dp,dp)
        do j = 1,d
            out(j,j) = cmplx(1.0_dp,0.0_dp,dp)
        end do
    end function identity_complex
    subroutine real_mv(x,y,istat)
        real(dp), intent(in) :: x(:) !! n-vector
        real(dp), intent(out) :: y(:) !! m-vector
        integer, intent(out) :: istat !! callback outcome
        fwd = fwd + 1
        istat = 31
        if (fail_forward) return
        y = matmul(a,x)
        if (emit_nan .and. size(y)>0) y(1) = ieee_value(0.0_dp,ieee_quiet_nan)
        istat = 0
    end subroutine real_mv
    subroutine real_rmv(x,y,istat)
        real(dp), intent(in) :: x(:) !! m-vector
        real(dp), intent(out) :: y(:) !! n-vector
        integer, intent(out) :: istat !! callback outcome
        adj = adj + 1
        istat = 32
        if (fail_adjoint) return
        y = matmul(transpose(a),x)
        if (emit_nan .and. size(y)>0) y(1) = ieee_value(0.0_dp,ieee_quiet_nan)
        istat = 0
    end subroutine real_rmv
    subroutine complex_mv(x,y,istat)
        complex(dp), intent(in) :: x(:) !! n-vector
        complex(dp), intent(out) :: y(:) !! m-vector
        integer, intent(out) :: istat !! callback outcome
        fwd = fwd + 1
        istat = 33
        if (fail_forward) return
        y = matmul(z,x)
        istat = 0
    end subroutine complex_mv
    subroutine complex_rmv(x,y,istat)
        complex(dp), intent(in) :: x(:) !! m-vector
        complex(dp), intent(out) :: y(:) !! n-vector
        integer, intent(out) :: istat !! callback outcome
        adj = adj + 1
        istat = 34
        if (fail_adjoint) return
        y = matmul(transpose(conjg(z)),x)
        istat = 0
    end subroutine complex_rmv
end program test_linalg_interpolative_svd_operator
