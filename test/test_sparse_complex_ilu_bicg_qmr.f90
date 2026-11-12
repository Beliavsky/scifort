! SPDX-License-Identifier: MIT
program test_sparse_complex_ilu_bicg_qmr
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        csr_complex_from_dense, csr_complex_from_coo, csr_complex_to_csc
    use scifort_sparse_complex_ilu, only : complex_ilu0_factor, ilu0_complex, &
        valid_ilu0_complex, ilu0_complex_solve, &
        complex_ilu_bad_diagonal, complex_ilu_invalid_input
    use scifort_sparse_linalg_complex, only : bicg_complex, qmr_complex, &
        bicgstab_complex, gmres_complex, sparse_complex_linalg_invalid_input
    implicit none
    type(csr_complex_matrix) :: a, bad
    type(csc_complex_matrix) :: csc
    type(complex_ilu0_factor) :: ilu
    complex(dp) :: d(5,5), full(5,5), target(5), rhs(5), one(5), z(2,2), empty(0,0)
    complex(dp), allocatable :: x(:), y(:)
    real(dp) :: residual
    integer :: info, steps, i

    d = cmplx(0.0_dp,0.0_dp,dp)
    do i = 1, 5
        d(i,i) = cmplx(8.0_dp+i,0.3_dp*i,dp)
        if (i < 5) then
            d(i,i+1) = cmplx(0.5_dp,-0.2_dp,dp)
            d(i+1,i) = cmplx(-0.7_dp,0.6_dp,dp)
        end if
    end do
    d(1,4) = cmplx(0.3_dp,0.4_dp,dp)
    d(5,2) = cmplx(-0.1_dp,0.5_dp,dp)
    target = [cmplx(1.0_dp,0.5_dp,dp),cmplx(-2.0_dp,0.2_dp,dp), &
        cmplx(3.0_dp,-0.5_dp,dp),cmplx(-0.3_dp,1.2_dp,dp),cmplx(-1.0_dp,-0.2_dp,dp)]
    rhs = matmul(d,target)
    a = csr_complex_from_dense(d,info)
    call require(info == 0,'CSR construction')
    csc = csr_complex_to_csc(a)
    call ilu0_complex(a,ilu,info)
    call require(info == 0 .and. valid_ilu0_complex(ilu),'ILU(0) factorization')
    call ilu0_complex_solve(ilu,rhs,x,info)
    call require(info == 0,'ILU(0) forward solve')
    ! With fill dropping, the approximate inverse need not solve Ax=b exactly.
    call require(all(abs(x) < 100.0_dp),'ILU(0) reasonable solution')
    one = matmul(transpose(conjg(d)),target)
    call ilu0_complex_solve(ilu,one,y,info,adjoint=.true.)
    call require(info == 0 .and. all(abs(y) < 100.0_dp),'ILU(0) adjoint solve')

    call bicg_complex(a,rhs,x,info,rtol=1.0e-12_dp,maxiter=30, &
        jacobi=.true.,iterations=steps,residual_norm=residual)
    call check_solution('BiCG CSR Jacobi')
    call bicg_complex(csc,rhs,x,info,rtol=1.0e-12_dp,maxiter=30, &
        ilu=ilu,iterations=steps,residual_norm=residual)
    call check_solution('BiCG CSC ILU')
    call bicg_complex(a,rhs,x,info,rtol=1.0e-12_dp,maxiter=30, &
        x0=target,iterations=steps,residual_norm=residual)
    call require(info == 0 .and. steps == 0,'BiCG exact initial guess')

    call qmr_complex(a,rhs,x,info,rtol=1.0e-12_dp,maxiter=30, &
        iterations=steps,residual_norm=residual)
    call check_solution('QMR CSR')
    call qmr_complex(csc,rhs,x,info,rtol=1.0e-12_dp,maxiter=30, &
        ilu=ilu,iterations=steps,residual_norm=residual)
    call check_solution('QMR CSC ILU')
    call qmr_complex(a,rhs,x,info,rtol=1.0e-12_dp,maxiter=30, &
        jacobi=.true.,iterations=steps,residual_norm=residual)
    call check_solution('QMR CSR Jacobi')

    call bicgstab_complex(a,rhs,x,info,rtol=1.0e-12_dp,maxiter=30, &
        ilu=ilu,iterations=steps,residual_norm=residual)
    call check_solution('BiCGSTAB CSR ILU')
    call gmres_complex(csc,rhs,x,info,rtol=1.0e-12_dp,maxiter=30,restart=5, &
        ilu=ilu,iterations=steps,residual_norm=residual)
    call check_solution('GMRES CSC ILU')

    call bicg_complex(a,rhs,x,info,ilu=ilu,jacobi=.true.)
    call require(info == sparse_complex_linalg_invalid_input,'preconditioners are exclusive')
    call qmr_complex(a,rhs,x,info,ilu=ilu,jacobi=.true.)
    call require(info == sparse_complex_linalg_invalid_input,'QMR preconditioners exclusive')

    z = reshape([cmplx(0.0_dp,0.0_dp,dp),cmplx(1.0_dp,1.0_dp,dp), &
        cmplx(0.0_dp,0.0_dp,dp),cmplx(2.0_dp,0.0_dp,dp)],[2,2])
    bad = csr_complex_from_dense(z)
    call ilu0_complex(bad,ilu,info)
    call require(info == complex_ilu_bad_diagonal,'detect missing diagonal')
    call ilu0_complex(bad,ilu,info,diagonal_shift=1.0_dp)
    call require(info == complex_ilu_bad_diagonal,'shift does not invent a missing diagonal')

    ! For a complete sparsity pattern, no nonzeros are dropped by ILU(0):
    ! its factorization and its conjugate transpose must solve exactly.
    full = d+cmplx(0.125_dp,0.03_dp,dp)
    a = csr_complex_from_dense(full)
    call ilu0_complex(a,ilu,info)
    call require(info == 0,'full-pattern ILU factor')
    rhs = matmul(full,target)
    call ilu0_complex_solve(ilu,rhs,x,info)
    call require(info == 0 .and. maxval(abs(x-target)) < 1.0e-12_dp, &
        'exact full-pattern LU solve')
    rhs = matmul(transpose(conjg(full)),target)
    call ilu0_complex_solve(ilu,rhs,x,info,adjoint=.true.)
    call require(info == 0 .and. maxval(abs(x-target)) < 1.0e-12_dp, &
        'exact full-pattern LU adjoint solve')

    ! An explicitly stored zero pivot is different from a missing diagonal.
    bad = csr_complex_from_coo(2,2,[1,1,2,2],[1,2,1,2], &
        [cmplx(0.0_dp,0.0_dp,dp),cmplx(1.0_dp,1.0_dp,dp), &
         cmplx(1.0_dp,0.0_dp,dp),cmplx(2.0_dp,0.0_dp,dp)],info)
    call require(info == 0,'explicit zero diagonal preserved')
    call ilu0_complex(bad,ilu,info)
    call require(info == complex_ilu_bad_diagonal,'zero pivot detected')
    call ilu0_complex(bad,ilu,info,diagonal_shift=1.0_dp)
    call require(info == 0 .and. valid_ilu0_complex(ilu),'diagonal shift resolves zero pivot')
    call ilu0_complex(bad,ilu,info,diagonal_shift=-1.0_dp)
    call require(info == complex_ilu_invalid_input,'negative diagonal shift rejected')

    a = csr_complex_from_dense(empty,info)
    call require(info == 0,'empty complex CSR input')
    call ilu0_complex(a,ilu,info)
    call require(info == 0 .and. valid_ilu0_complex(ilu),'empty ILU factor')
    call ilu0_complex_solve(ilu,[complex(dp) ::],x,info)
    call require(info == 0 .and. size(x) == 0,'empty ILU solve')
    call bicg_complex(a,[complex(dp) ::],x,info,iterations=steps)
    call require(info == 0 .and. steps == 0,'empty BiCG')
    call qmr_complex(a,[complex(dp) ::],x,info,iterations=steps)
    call require(info == 0 .and. steps == 0,'empty QMR')

    print *, 'PASS complex ILU(0), BiCG, QMR' 
contains
    subroutine check_solution(message)
        character(*), intent(in) :: message !! Test description.
        if (info /= 0 .or. steps == 0 .or. residual > 2.0e-9_dp) then
            print *, trim(message), ' info/steps/residual:',info,steps,residual
            error stop 'complex iterative failure'
        end if
        if (maxval(abs(x-target)) > 2.0e-9_dp) then
            print *, trim(message), ' solution error:', maxval(abs(x-target))
            error stop 'wrong complex solution'
        end if
    end subroutine check_solution

    subroutine require(ok, message)
        logical, intent(in) :: ok !! Assertion condition.
        character(*), intent(in) :: message !! Failure explanation.
        if (.not. ok) then
            print *, trim(message)
            error stop 'complex regression failure'
        end if
    end subroutine require
end program test_sparse_complex_ilu_bicg_qmr
