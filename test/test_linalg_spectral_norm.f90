! SPDX-License-Identifier: MIT
program test_linalg_spectral_norm
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_spectral_norm, only : estimate_spectral_norm, &
        estimate_spectral_norm_diff, estimate_spectral_norm_operator, &
        estimate_spectral_norm_operator_complex, estimate_spectral_norm_diff_operator, &
        estimate_spectral_norm_diff_operator_complex
    implicit none
    real(dp), allocatable :: a(:, :), b(:, :), x(:)
    complex(dp), allocatable :: c(:, :), d(:, :), xc(:)
    real(dp) :: e1, e2, scale, nan
    integer :: i, j, m, n, variant, status, checks
    checks = 0
    do m = 0, 12
        do n = 0, 13
            allocate(a(m,n),b(m,n),c(m,n),d(m,n),x(n),xc(n))
            do variant = 1, 5
                do j = 1, n
                    x(j) = cos(real(j,dp)) + 0.2_dp
                    xc(j) = cmplx(x(j),sin(real(j,dp)),dp)
                    do i = 1, m
                        a(i,j) = sin(real(i*7+j*11+variant,dp)) + &
                            0.2_dp*cos(real(i*j+variant,dp))
                        b(i,j) = 0.2_dp*cos(real(i*13-j*5,dp))
                        c(i,j) = cmplx(a(i,j),sin(real(i*3-j*8,dp)),dp)
                        d(i,j) = cmplx(b(i,j),0.3_dp*cos(real(i*j,dp)),dp)
                    end do
                end do
                select case (variant)
                case (2)
                    a = 0.0_dp
                    b = 0.0_dp
                    c = cmplx(0.0_dp,0.0_dp,dp)
                    d = cmplx(0.0_dp,0.0_dp,dp)
                case (3)
                    if (m*n > 0) then
                        do j = 1, n
                            a(:,j) = a(:,1) * real(j,dp)
                            c(:,j) = c(:,1) * cmplx(real(j,dp),0.2_dp*real(j,dp),dp)
                        end do
                    end if
                case (4)
                    a = a * 1.0e100_dp
                    b = b * 1.0e100_dp
                    c = c * 1.0e100_dp
                    d = d * 1.0e100_dp
                case (5)
                    a = a * 1.0e-100_dp
                    b = b * 1.0e-100_dp
                    c = c * 1.0e-100_dp
                    d = d * 1.0e-100_dp
                end select
                ! The matrix-vector operator and dense implementation must
                ! yield identical power iterates from the same initial vector.
                call estimate_spectral_norm(a,e1,status,its=25,start=x)
                if (status /= linalg_status_success) error stop 'real dense'
                call estimate_spectral_norm_operator(m,n,mv_a,rmv_a,e2,status,its=25,start=x)
                if (status /= linalg_status_success) error stop 'real operator'
                scale = max(1.0e-120_dp,maxval(abs(a)))
                if (abs(e1-e2) > 2.0e-12_dp * scale) then
                    print *, 'real comparison', m,n,variant,e1,e2
                    error stop 'real dense vs operator'
                end if
                call estimate_spectral_norm(c,e1,status,its=25,start=xc)
                if (status /= linalg_status_success) error stop 'complex dense'
                call estimate_spectral_norm_operator_complex(m,n,mv_c,rmv_c,e2,status,its=25,start=xc)
                if (status /= linalg_status_success) error stop 'complex operator'
                scale = max(1.0e-120_dp,maxval(abs(c)))
                if (abs(e1-e2) > 2.0e-12_dp*scale) then
                    print *, 'complex comparison',m,n,variant,e1,e2
                    error stop 'complex dense vs operator'
                end if
                call estimate_spectral_norm_diff(a,b,e1,status,its=25,start=x)
                if (status /= linalg_status_success) error stop 'real difference dense'
                call estimate_spectral_norm_diff_operator( &
                    m,n,mv_a,rmv_a,mv_b,rmv_b,e2,status,its=25,start=x)
                if (status /= linalg_status_success) error stop 'real difference operator'
                scale = max(1.0e-120_dp,max(maxval(abs(a)),maxval(abs(b))))
                if (abs(e1-e2) > 2.0e-12_dp*scale) error stop 'real difference mismatch'
                call estimate_spectral_norm_diff(c,d,e1,status,its=25,start=xc)
                if (status /= linalg_status_success) error stop 'complex difference dense'
                call estimate_spectral_norm_diff_operator_complex( &
                    m,n,mv_c,rmv_c,mv_d,rmv_d,e2,status,its=25,start=xc)
                if (status /= linalg_status_success) error stop 'complex difference operator'
                scale = max(1.0e-120_dp,max(maxval(abs(c)),maxval(abs(d))))
                if (abs(e1-e2) > 2.0e-12_dp*scale) error stop 'complex difference mismatch'
                checks = checks + 8
            end do
            deallocate(a,b,c,d,x,xc)
        end do
    end do
    ! Rank-one, known spectral norm exactly sqrt(125).
    allocate(a(2,2),b(2,2),c(2,2),d(2,2),x(2),xc(2))
    a(1,:) = [3.0_dp,4.0_dp]
    a(2,:) = [6.0_dp,8.0_dp]
    call estimate_spectral_norm(a,e1,status,its=3,seed=42_int64)
    if (status /= 0 .or. abs(e1-sqrt(125.0_dp)) > 1.0e-13_dp) error stop 'rank one real'
    c = cmplx(a,0.0_dp,dp) * cmplx(0.0_dp,1.0_dp,dp)
    call estimate_spectral_norm(c,e1,status,its=3,seed=42_int64)
    if (status /= 0 .or. abs(e1-sqrt(125.0_dp)) > 1.0e-13_dp) error stop 'rank one complex'
    checks = checks + 2
    ! Optional seeds give deterministic results independent of random state.
    call estimate_spectral_norm(a,e1,status,its=0,seed=17_int64)
    call estimate_spectral_norm(a,e2,status,its=0,seed=17_int64)
    if (e1 /= e2) error stop 'repeatable seed'
    call estimate_spectral_norm(c,e1,status,its=0,seed=17_int64)
    call estimate_spectral_norm(c,e2,status,its=0,seed=17_int64)
    if (e1 /= e2) error stop 'repeatable complex seed'
    checks = checks+2
    ! Error handling, including shape/iteration and invalid callback.
    call estimate_spectral_norm(a,e1,status,its=-1)
    if (status /= linalg_status_invalid_shape) error stop 'negative iterations'
    call estimate_spectral_norm(a,e1,status,start=[1.0_dp])
    if (status /= linalg_status_invalid_shape) error stop 'wrong start length'
    call estimate_spectral_norm(a,e1,status,start=[0.0_dp,0.0_dp])
    if (status /= linalg_status_invalid_shape) error stop 'zero initial vector'
    call estimate_spectral_norm_diff(a,a(:,1:1),e1,status)
    if (status /= linalg_status_invalid_shape) error stop 'wrong difference shape'
    call estimate_spectral_norm_operator(2,2,bad_mv,rmv_a,e1,status)
    if (status /= linalg_status_no_convergence) error stop 'failed callback'
    nan = ieee_value(0.0_dp,ieee_quiet_nan)
    a(1,1) = nan
    call estimate_spectral_norm(a,e1,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite input'
    call estimate_spectral_norm_diff(a,b,e1,status)
    if (status /= linalg_status_no_convergence) error stop 'nonfinite difference'
    checks = checks+7
    ! Nearly equal, large entries: subtraction must precede scaling.
    a = 0.0_dp
    b = 0.0_dp
    a(1,1) = 1.0e100_dp + 1.0e88_dp
    b(1,1) = 1.0e100_dp
    call estimate_spectral_norm_diff(a,b,e1,status,its=3,seed=9_int64)
    if (status /= linalg_status_success) error stop 'large cancellation status'
    if (abs(e1-(a(1,1)-b(1,1))) > abs(a(1,1)-b(1,1))*1.0e-13_dp) &
        error stop 'large cancellation result'
    checks = checks + 1
    ! An unrepresentable norm must produce a non-success status, not Inf.
    a = 0.0_dp
    b = 0.0_dp
    a(1,1) = 0.6_dp*huge(1.0_dp)
    b(1,1) = -a(1,1)
    call estimate_spectral_norm_diff(a,b,e1,status,its=3,seed=9_int64)
    if (status /= linalg_status_no_convergence) error stop 'overflow status'
    checks = checks + 1
    print *, 'spectral norm tests passed:', checks
contains
    subroutine mv_a(v,y,st)
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: y(:)
        integer, intent(out) :: st
        y = matmul(a,v)
        st = 0
    end subroutine mv_a
    subroutine rmv_a(v,y,st)
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: y(:)
        integer, intent(out) :: st
        y = matmul(transpose(a),v)
        st = 0
    end subroutine rmv_a
    subroutine mv_b(v,y,st)
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: y(:)
        integer, intent(out) :: st
        y = matmul(b,v)
        st = 0
    end subroutine mv_b
    subroutine rmv_b(v,y,st)
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: y(:)
        integer, intent(out) :: st
        y = matmul(transpose(b),v)
        st = 0
    end subroutine rmv_b
    subroutine mv_c(v,y,st)
        complex(dp), intent(in) :: v(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: st
        y = matmul(c,v)
        st = 0
    end subroutine mv_c
    subroutine rmv_c(v,y,st)
        complex(dp), intent(in) :: v(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: st
        y = matmul(transpose(conjg(c)),v)
        st = 0
    end subroutine rmv_c
    subroutine mv_d(v,y,st)
        complex(dp), intent(in) :: v(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: st
        y = matmul(d,v)
        st = 0
    end subroutine mv_d
    subroutine rmv_d(v,y,st)
        complex(dp), intent(in) :: v(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: st
        y = matmul(transpose(conjg(d)),v)
        st = 0
    end subroutine rmv_d
    subroutine bad_mv(v,y,st)
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: y(:)
        integer, intent(out) :: st
        y = 0.0_dp
        st = 3
    end subroutine bad_mv
end program test_linalg_spectral_norm
