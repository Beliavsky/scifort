! SPDX-License-Identifier: MIT
! Complex sparse matrix-function actions: Arnoldi/Lanczos, restart and callbacks.
program test_sparse_complex_funm
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        csr_complex_from_dense, csr_complex_to_csc
    use scifort_sparse_linalg_complex_funm, only : &
        funm_multiply_krylov_complex, funm_multiply_krylov_complex_operator, &
        complex_funm_invalid_input, complex_funm_breakdown, complex_funm_callback_failure, &
        complex_funm_no_convergence
    use scifort_sparse_linalg_complex_expm, only : expm_multiply_complex
    implicit none

    type(csr_complex_matrix) :: a
    type(csc_complex_matrix) :: ac
    complex(dp) :: dense(12,12), rhs(12), rhsblock(12,2), t, z(2)
    complex(dp), allocatable :: y(:), yref(:), yo(:), yc(:), yblock(:,:), expectedblock(:,:)
    real(dp) :: update, nanval
    integer :: i, info, steps, products, k

    dense = cmplx(0.0_dp,0.0_dp,dp)
    dense(1,1) = cmplx(-0.4_dp,0.2_dp,dp)
    dense(2,2) = cmplx(0.1_dp,-0.8_dp,dp)
    dense(1,2) = cmplx(0.3_dp,-1.0_dp,dp)
    dense(2,1) = cmplx(-0.5_dp,0.9_dp,dp)
    a = csr_complex_from_dense(dense(:2,:2))
    ac = csr_complex_to_csc(a)
    z = [cmplx(1.0_dp,1.0_dp,dp),cmplx(0.2_dp,-1.2_dp,dp)]
    t = cmplx(0.35_dp,-0.19_dp,dp)
    call funm_multiply_krylov_complex(a,z,poly_f,y,info,t=t,cycles=steps,products=products)
    call require(info == 0 .and. steps == 1 .and. products == 2,'full Arnoldi polynomial status')
    yref = matmul(t*dense(:2,:2),matmul(t*dense(:2,:2),z)) + &
           cmplx(0.2_dp,-0.3_dp,dp)*matmul(t*dense(:2,:2),z) + &
           cmplx(0.75_dp,0.1_dp,dp)*z
    call close(y,yref,1.0e-12_dp,'full Arnoldi polynomial')
    call funm_multiply_krylov_complex(ac,z,poly_f,yc,info,t=t)
    call require(info == 0,'CSC matrix function')
    call close(yc,yref,1.0e-12_dp,'CSC polynomial')
    call funm_multiply_krylov_complex_operator(2,z,small_operator,poly_f,yo,info,t=t)
    call require(info == 0,'matrix-free full Arnoldi')
    call close(yo,yref,1.0e-12_dp,'matrix-free full Arnoldi')

    ! Independent action of the exponential on a nonsymmetric complex matrix.
    dense = cmplx(0.0_dp,0.0_dp,dp)
    do i = 1,12
        dense(i,i) = cmplx(-0.9_dp-0.07_dp*real(i,dp),0.1_dp*real(i,dp),dp)
        if (i < 12) dense(i,i+1) = cmplx(0.45_dp,0.3_dp,dp)
        if (i > 1) dense(i,i-1) = cmplx(-0.21_dp,0.12_dp,dp)
        rhs(i) = cmplx(sin(0.25_dp*real(i,dp)),cos(0.21_dp*real(i,dp)),dp)
    end do
    a = csr_complex_from_dense(dense)
    ac = csr_complex_to_csc(a)
    t = cmplx(0.33_dp,0.0_dp,dp)
    call funm_multiply_krylov_complex(a,rhs,exp_f,y,info,t=t, &
        restart_every_m=4,max_restarts=6,rtol=1.0e-10_dp, &
        cycles=steps,products=products,update_norm=update)
    call require(info == 0 .and. steps > 1 .and. products >= steps*4,'restarted Arnoldi status')
    call expm_multiply_complex(a,rhs,yref,info,t=0.33_dp)
    call require(info == 0,'complex exponential reference')
    call close(y,yref,1.0e-8_dp,'restarted Arnoldi exponential')
    call funm_multiply_krylov_complex(ac,rhs,exp_f,yc,info,t=t, &
        restart_every_m=4,max_restarts=6,rtol=1.0e-10_dp)
    call require(info == 0,'restarted CSC')
    call close(yc,yref,1.0e-8_dp,'restarted CSC exp')

    ! Multiple RHS are evaluated without a dense input matrix-function evaluation.
    rhsblock(:,1) = rhs
    rhsblock(:,2) = 0.25_dp*rhs + cmplx(0.0_dp,0.8_dp,dp)
    call funm_multiply_krylov_complex(a,rhsblock,exp_f,yblock,info,t=t, &
        restart_every_m=4,max_restarts=6,rtol=1.0e-10_dp, &
        cycles=steps,products=products)
    call require(info == 0 .and. steps > 2 .and. products >= steps*4,'block CSR status')
    call expm_multiply_complex(a,rhsblock,expectedblock,info,t=0.33_dp)
    call require(info == 0,'block exponential reference')
    call close(reshape(yblock,[24]),reshape(expectedblock,[24]),1.0e-8_dp,'block CSR exp')
    call funm_multiply_krylov_complex(ac,rhsblock,exp_f,yblock,info,t=t, &
        restart_every_m=4,max_restarts=6,rtol=1.0e-10_dp)
    call require(info == 0,'block CSC status')
    call close(reshape(yblock,[24]),reshape(expectedblock,[24]),1.0e-8_dp,'block CSC exp')

    ! Complex Hermitian operator, using the real-tridiagonal Lanczos recurrence.
    dense = cmplx(0.0_dp,0.0_dp,dp)
    do i = 1,12
        dense(i,i) = cmplx(-0.8_dp+0.11_dp*real(i,dp),0.0_dp,dp)
        if (i < 12) then
            dense(i+1,i) = cmplx(0.38_dp,0.21_dp,dp)
            dense(i,i+1) = conjg(dense(i+1,i))
        end if
    end do
    a = csr_complex_from_dense(dense)
    t = cmplx(0.0_dp,0.27_dp,dp)
    call funm_multiply_krylov_complex(a,rhs,exp_f,y,info,t=t, &
        restart_every_m=5,max_restarts=7,rtol=1.0e-11_dp,assume_hermitian=.true.)
    call require(info == 0,'complex Hermitian Lanczos status')
    call close(y,exp_action(t*dense,rhs),2.0e-9_dp,'complex Lanczos')

    ! Explicitly report unverified truncation with one short Arnoldi cycle.
    call funm_multiply_krylov_complex(a,rhs,exp_f,y,info,restart_every_m=2,max_restarts=1)
    call require(info == complex_funm_no_convergence,'finite iteration budget')

    ! Zero-column RHS blocks still validate options and perform no callbacks.
    call funm_multiply_krylov_complex(a,rhsblock(:,:0),broken_f,yblock,info)
    call require(info == 0 .and. size(yblock,1) == 12 .and. &
        size(yblock,2) == 0,'zero-column block')
    call funm_multiply_krylov_complex(a,rhsblock(:,:0),broken_f,yblock,info,max_restarts=0)
    call require(info == complex_funm_invalid_input .and. size(yblock) == 0,'zero block bad options')

    ! An oversized basis is clamped to n, without creating a huge projection.
    call funm_multiply_krylov_complex(a,rhs,poly_f,y,info,restart_every_m=huge(0))
    call require(info == 0,'oversized basis clamp')

    ! Zero input must not invoke the dense function or matrix-free operator.
    call funm_multiply_krylov_complex(a,0.0_dp*rhs,broken_f,y,info)
    call require(info == 0 .and. all(y == cmplx(0.0_dp,0.0_dp,dp)),'zero RHS')

    ! Reject invalid shapes, nonfinite entries, callback failure, and invalid options.
    call funm_multiply_krylov_complex(a,rhs(:11),exp_f,y,info)
    call require(info == complex_funm_invalid_input .and. size(y) == 0,'wrong RHS length')
    call funm_multiply_krylov_complex(a,rhs,exp_f,y,info,restart_every_m=0)
    call require(info == complex_funm_invalid_input,'zero basis')
    call funm_multiply_krylov_complex(a,rhs,exp_f,y,info,rtol=-1.0_dp)
    call require(info == complex_funm_invalid_input,'negative tolerance')
    nanval = ieee_value(0.0_dp,ieee_quiet_nan)
    call funm_multiply_krylov_complex(a,rhs,exp_f,y,info,t=cmplx(nanval,0.0_dp,dp))
    call require(info == complex_funm_invalid_input,'nonfinite time')
    call funm_multiply_krylov_complex(a,rhs,broken_f,y,info)
    call require(info == complex_funm_callback_failure,'failed dense evaluator')
    call funm_multiply_krylov_complex(a,rhs,nan_f,y,info)
    call require(info == complex_funm_breakdown,'nonfinite dense evaluator')
    call funm_multiply_krylov_complex_operator(12,rhs,broken_op,exp_f,y,info)
    call require(info == complex_funm_callback_failure,'failed matvec callback')
    call funm_multiply_krylov_complex_operator(12,rhs,nan_op,exp_f,y,info)
    call require(info == complex_funm_breakdown,'nonfinite matvec callback')
    call funm_multiply_krylov_complex_operator(0,[complex(dp)::],broken_op,broken_f,y,info)
    call require(info == 0 .and. size(y) == 0,'zero order operator')

    ! Large matrix-free complex operator with an exactly linear matrix function.
    call test_large_matrix_free()
    print '(a)', 'test_sparse_complex_funm: PASS'

contains
    subroutine require(condition,message)
        logical, intent(in) :: condition !! Expected invariant.
        character(len=*), intent(in) :: message !! Failure label.
        if (.not. condition) then
            print '(3a)', 'FAIL: ',message,''
            error stop 1
        end if
    end subroutine require

    subroutine close(got, reference, tolerance, message)
        complex(dp), intent(in) :: got(:) !! Computed values.
        complex(dp), intent(in) :: reference(:) !! Expected values.
        real(dp), intent(in) :: tolerance !! Absolute tolerance.
        character(len=*), intent(in) :: message !! Failure label.
        call require(size(got) == size(reference),message//' shape')
        if (size(got) > 0) then
            if (maxval(abs(got-reference)) > tolerance) then
                print *, message, ' error: ',maxval(abs(got-reference))
                error stop 2
            end if
        end if
    end subroutine close

    subroutine poly_f(h,f,stat)
        complex(dp), intent(in) :: h(:,:) !! Small projected matrix.
        complex(dp), intent(out) :: f(:,:) !! Polynomial applied to matrix.
        integer, intent(out) :: stat !! Zero on success.
        integer :: j
        f = matmul(h,h) + cmplx(0.2_dp,-0.3_dp,dp)*h
        do j=1,size(h,1)
            f(j,j) = f(j,j) + cmplx(0.75_dp,0.1_dp,dp)
        end do
        stat = 0
    end subroutine poly_f

    subroutine exp_f(h,f,stat)
        complex(dp), intent(in) :: h(:,:) !! Projected matrix.
        complex(dp), intent(out) :: f(:,:) !! exp(h).
        integer, intent(out) :: stat !! Status.
        complex(dp), allocatable :: term(:,:), small(:,:)
        integer :: j, k, scale_steps, n
        real(dp) :: norm_bound
        n = size(h,1)
        allocate(term(n,n),small(n,n))
        norm_bound = 0.0_dp
        do j=1,n
            norm_bound = max(norm_bound,sum(abs(h(j,:))))
        end do
        scale_steps = 0
        if (norm_bound > 0.5_dp) scale_steps = max(0,ceiling(log(norm_bound/0.5_dp)/log(2.0_dp)))
        scale_steps = min(scale_steps,30)
        small = h*2.0_dp**(-scale_steps)
        f = cmplx(0.0_dp,0.0_dp,dp)
        term = cmplx(0.0_dp,0.0_dp,dp)
        do j=1,n
            f(j,j) = cmplx(1.0_dp,0.0_dp,dp)
            term(j,j) = cmplx(1.0_dp,0.0_dp,dp)
        end do
        do k=1,100
            term = matmul(term,small)/real(k,dp)
            f = f+term
            if (maxval(abs(term)) <= 2.0e-16_dp*max(1.0_dp,maxval(abs(f)))) exit
        end do
        do j=1,scale_steps
            f = matmul(f,f)
        end do
        stat=0
    end subroutine exp_f

    function exp_action(h,b) result(out)
        complex(dp), intent(in) :: h(:,:) !! Dense input matrix.
        complex(dp), intent(in) :: b(:) !! Vector to multiply.
        complex(dp) :: out(size(b)), fh(size(b),size(b))
        integer :: ierr
        call exp_f(h,fh,ierr)
        out = matmul(fh,b)
    end function exp_action

    subroutine small_operator(x,ax,stat)
        complex(dp), intent(in) :: x(:) !! Vector length 2.
        complex(dp), intent(out) :: ax(:) !! Product with initial 2x2 matrix.
        integer, intent(out) :: stat !! Zero on success.
        complex(dp), parameter :: matrix(2,2) = reshape([ &
            cmplx(-0.4_dp,0.2_dp,dp),cmplx(-0.5_dp,0.9_dp,dp), &
            cmplx(0.3_dp,-1.0_dp,dp),cmplx(0.1_dp,-0.8_dp,dp)],[2,2])
        ax = matmul(matrix,x)
        stat = 0
    end subroutine small_operator

    subroutine broken_f(h,f,stat)
        complex(dp), intent(in) :: h(:,:) !! Ignored.
        complex(dp), intent(out) :: f(:,:) !! Ignored.
        integer, intent(out) :: stat !! Failure status.
        f = cmplx(0.0_dp,0.0_dp,dp)
        stat = 1
    end subroutine broken_f

    subroutine nan_f(h,f,stat)
        complex(dp), intent(in) :: h(:,:) !! Ignored.
        complex(dp), intent(out) :: f(:,:) !! Nonfinite output.
        integer, intent(out) :: stat !! Zero but invalid numerical result.
        real(dp) :: qnan
        qnan = ieee_value(0.0_dp,ieee_quiet_nan)
        f = cmplx(qnan,0.0_dp,dp)
        stat = 0
    end subroutine nan_f

    subroutine broken_op(x,ax,stat)
        complex(dp), intent(in) :: x(:) !! Ignored.
        complex(dp), intent(out) :: ax(:) !! Output vector.
        integer, intent(out) :: stat !! Failure status.
        ax = cmplx(0.0_dp,0.0_dp,dp)
        stat = 1
    end subroutine broken_op

    subroutine nan_op(x,ax,stat)
        complex(dp), intent(in) :: x(:) !! Ignored.
        complex(dp), intent(out) :: ax(:) !! Invalid vector.
        integer, intent(out) :: stat !! Zero status.
        real(dp) :: qnan
        qnan = ieee_value(0.0_dp,ieee_quiet_nan)
        ax = cmplx(qnan,0.0_dp,dp)
        stat = 0
    end subroutine nan_op

    subroutine test_large_matrix_free()
        integer, parameter :: nlarge=1024
        complex(dp) :: b(nlarge), expected(nlarge), ax(nlarge)
        complex(dp), allocatable :: computed(:)
        integer :: j, stat, calls, cycles_done
        do j=1,nlarge
            b(j) = cmplx(sin(real(j,dp)*0.013_dp),cos(real(j,dp)*0.037_dp),dp)
        end do
        call tridiag_op(b,ax,stat)
        expected = b+cmplx(0.2_dp,-0.4_dp,dp)*ax
        call funm_multiply_krylov_complex_operator(nlarge,b,tridiag_op,linear_f, &
            computed,stat,restart_every_m=4,max_restarts=3,rtol=1.0e-11_dp, &
            cycles=cycles_done,products=calls)
        call require(stat == 0 .and. calls <= 12 .and. cycles_done >= 2,'large matrix free status')
        call close(computed,expected,5.0e-10_dp,'large matrix free affine')
    end subroutine test_large_matrix_free

    subroutine tridiag_op(x,ax,stat)
        complex(dp), intent(in) :: x(:) !! Complex input.
        complex(dp), intent(out) :: ax(:) !! Tridiagonal action.
        integer, intent(out) :: stat !! Zero on success.
        integer :: j
        ax = cmplx(-0.7_dp,0.1_dp,dp)*x
        do j=2,size(x)
            ax(j) = ax(j) + cmplx(0.03_dp,-0.02_dp,dp)*x(j-1)
            ax(j-1) = ax(j-1) + cmplx(-0.04_dp,-0.01_dp,dp)*x(j)
        end do
        stat = 0
    end subroutine tridiag_op

    subroutine linear_f(h,f,stat)
        complex(dp), intent(in) :: h(:,:) !! Projected matrix.
        complex(dp), intent(out) :: f(:,:) !! I+(0.2-0.4i)*h.
        integer, intent(out) :: stat !! Success flag.
        integer :: j
        f = cmplx(0.2_dp,-0.4_dp,dp)*h
        do j=1,size(h,1)
            f(j,j)=f(j,j)+cmplx(1.0_dp,0.0_dp,dp)
        end do
        stat = 0
    end subroutine linear_f
end program test_sparse_complex_funm
