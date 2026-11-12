! SPDX-License-Identifier: MIT
program test_sparse_arnoldi
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense, csr_from_coo
    use scifort_sparse_linalg, only : eigs_arnoldi, sparse_linalg_invalid_input
    implicit none
    integer, parameter :: n = 64
    real(dp) :: dense(n,n), x, y, err
    type(csr_matrix) :: a
    complex(dp), allocatable :: values(:), vectors(:,:)
    real(dp), allocatable :: residuals(:)
    integer :: i, j, info, cycles, locked, operations, mode
    character(len=2) :: target

    do mode = 1,5
        dense = 0.0_dp
        do i = 1,n
            dense(i,i) = 0.4_dp*cos(real(i,dp))
            if (i < n) then
                dense(i,i+1) = 0.01_dp
                dense(i+1,i) = -0.01_dp
            end if
        end do
        select case(mode)
        case(1,2)
            target = 'LM'
            if (mode == 2) target = 'LR'
            dense(1,1) = 8.0_dp
            dense(2,2) = 7.0_dp
        case(3)
            target = 'SR'
            dense(1,1) = -8.0_dp
            dense(2,2) = -7.0_dp
        case(4,5)
            target = 'LI'
            if (mode == 5) target = 'SI'
            dense(1:4,1:4) = 0.0_dp
            dense(1,2) = -8.0_dp
            dense(2,1) = 8.0_dp
            dense(3,4) = -7.0_dp
            dense(4,3) = 7.0_dp
        end select
        a = csr_from_dense(dense)
        call eigs_arnoldi(a,2,values,vectors,info,which=target,tol=1.0e-8_dp, &
            ncv=6,maxiter=60,iterations=cycles,locked_count=locked, &
            matvecs=operations,residual_norms=residuals)
        call require(info == 0, 'restarted Arnoldi convergence '//target)
        call require(cycles > 1 .and. locked == 2, 'soft locking after restart')
        call require(operations > 0 .and. size(residuals) == 2, 'solver diagnostics')
        call require(maxval(residuals) < 5.0e-7_dp, 'true Ritz residual '//target)
        do j = 1,2
            x = sqrt(sum(abs(matmul(dense,vectors(:,j))-values(j)*vectors(:,j))**2))
            call require(abs(x-residuals(j)) < 1.0e-10_dp, 'independent residual')
            y = sqrt(sum(abs(vectors(:,j))**2))
            call require(abs(y-1.0_dp)<1.0e-10_dp, 'unit eigenvector')
        end do
        select case(mode)
        case(1,2)
            call require(abs(real(values(1),dp)-8.0_dp) < 0.01_dp, 'largest eigenvalue')
        case(3)
            call require(abs(real(values(1),dp)+8.0_dp) < 0.01_dp, 'smallest eigenvalue')
        case(4)
            call require(abs(aimag(values(1))-8.0_dp) < 0.01_dp, 'positive complex eigenvalue')
        case(5)
            call require(abs(aimag(values(1))+8.0_dp) < 0.01_dp, 'negative complex eigenvalue')
        end select
    end do

    ! A large sparse operator, without ever forming a dense n-by-n matrix.
    block
        integer, parameter :: nn = 1200
        integer :: rows(3*nn-2), cols(3*nn-2), t
        real(dp) :: entries(3*nn-2)
        do t = 1,nn
            rows(t) = t
            cols(t) = t
            entries(t) = 0.3_dp*cos(real(t,dp))
        end do
        entries(1) = 9.0_dp
        entries(2) = 8.0_dp
        do t = 1,nn-1
            rows(nn+t) = t
            cols(nn+t) = t+1
            entries(nn+t) = 0.015_dp
            rows(2*nn-1+t) = t+1
            cols(2*nn-1+t) = t
            entries(2*nn-1+t) = -0.012_dp
        end do
        a = csr_from_coo(nn,nn,rows,cols,entries)
        call eigs_arnoldi(a,2,values,vectors,info,which='LR',tol=1.0e-8_dp, &
            ncv=10,maxiter=70,residual_norms=residuals)
        call require(info == 0, '1200-row sparse Arnoldi convergence')
        call require(maxval(residuals) < 1.0e-6_dp, 'large sparse true residual')
        call require(abs(real(values(1),dp)-9.0_dp) < 0.01_dp, 'large sparse eigenvalue')
    end block

    dense = 0.0_dp
    a = csr_from_dense(dense)
    call eigs_arnoldi(a,2,values,vectors,info,ncv=8,maxiter=8,residual_norms=residuals)
    call require(info == 0 .and. all(residuals < 1.0e-12_dp), 'zero matrix')
    call eigs_arnoldi(a,2,values,vectors,info,maxiter=0,iterations=cycles)
    call require(info > 0 .and. cycles == 0, 'zero iteration cap')
    call eigs_arnoldi(a,2,values,vectors,info,which='XX')
    call require(info == sparse_linalg_invalid_input, 'invalid selection')
    call eigs_arnoldi(a,2,values,vectors,info,which='SM',ncv=8)
    call require(info == sparse_linalg_invalid_input, 'reduced smallest magnitude')
    call eigs_arnoldi(a,2,values,vectors,info,ncv=2)
    call require(info == sparse_linalg_invalid_input, 'invalid Krylov width')
    call eigs_arnoldi(a,2,values,vectors,info,tol=-0.5_dp)
    call require(info == sparse_linalg_invalid_input, 'negative tolerance')
    call eigs_arnoldi(a,2,values,vectors,info,v0= &
        [(cmplx(0.0_dp,0.0_dp,dp),i=1,n)])
    call require(info == sparse_linalg_invalid_input, 'zero starting vector')
    call eigs_arnoldi(a,2,values,vectors,info,v0= &
        [(cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp),i=1,n)])
    call require(info == sparse_linalg_invalid_input, 'nonfinite starting vector')
    a = csr_matrix()
    call eigs_arnoldi(a,1,values,vectors,info)
    call require(info == sparse_linalg_invalid_input, 'invalid CSR')
    print '(a)', 'test_sparse_arnoldi: PASS'
contains
    subroutine require(condition, message)
        logical, intent(in) :: condition !! Required condition for regression.
        character(len=*), intent(in) :: message !! Description of the assertion.
        if (.not. condition) then
            print '(a)', 'FAILED: '//message
            error stop 1
        end if
    end subroutine require
end program test_sparse_arnoldi
