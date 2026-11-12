! SPDX-License-Identifier: MIT
program test_sparse_lanczos
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense, csr_from_coo, csr_matvec
    use scifort_sparse_linalg, only : eigsh_lanczos, eigsh, eigs, &
        sparse_linalg_invalid_input
    implicit none
    integer, parameter :: n = 80
    type(csr_matrix) :: a
    real(dp) :: matrix(6,6), sample(6), nan_value
    real(dp), allocatable :: values(:), vectors(:,:), errors(:), y(:)
    complex(dp), allocatable :: cvalues(:), cvectors(:,:)
    integer :: i, j, info, cycles, locked, products
    character(len=2) :: choice

    sample = [-7.0_dp,-2.0_dp,0.0_dp,1.0_dp,3.0_dp,8.0_dp]
    matrix = 0.0_dp
    do i = 1,6
        matrix(i,i) = sample(i)
    end do
    a = csr_from_dense(matrix)
    do j = 1,4
        select case (j)
        case (1)
            choice = 'LM'
        case (2)
            choice = 'LA'
        case (3)
            choice = 'SA'
        case (4)
            choice = 'SM'
        end select
        call eigsh_lanczos(a,2,values,vectors,info,which=choice,ncv=6, &
            tol=1.0e-10_dp,maxiter=10,iterations=cycles, &
            locked_count=locked,matvecs=products,residual_norms=errors)
        call require(info == 0, 'Lanczos full-subspace '//choice)
        call require(locked == 2 .and. products > 0 .and. cycles > 0, 'diagnostics')
        select case (choice)
        case ('LM')
            call close(values,[-7.0_dp,8.0_dp],1.0e-9_dp,'LM')
        case ('LA')
            call close(values,[3.0_dp,8.0_dp],1.0e-9_dp,'LA')
        case ('SA')
            call close(values,[-7.0_dp,-2.0_dp],1.0e-9_dp,'SA')
        case ('SM')
            call close(values,[0.0_dp,1.0_dp],1.0e-9_dp,'SM')
        end select
        call require(maxval(errors) < 1.0e-9_dp, 'full-space residuals')
        call require(abs(dot_product(vectors(:,1),vectors(:,2))) < 1.0e-10_dp, &
            'orthogonal Ritz vectors')
    end do

    ! Nontrivial sparse symmetric spectrum, reduced Krylov width and restart.
    block
        integer :: rows(3*n-2), columns(3*n-2)
        real(dp) :: entries(3*n-2)
        do i = 1,n
            rows(i) = i
            columns(i) = i
            entries(i) = -0.3_dp + 0.04_dp*cos(real(i,dp))
        end do
        entries(n) = 9.0_dp
        entries(n-1) = 8.0_dp
        entries(1) = -10.0_dp
        entries(2) = -9.0_dp
        do i = 1,n-1
            rows(n+i) = i
            columns(n+i) = i+1
            entries(n+i) = 0.012_dp
            rows(2*n-1+i) = i+1
            columns(2*n-1+i) = i
            entries(2*n-1+i) = 0.012_dp
        end do
        a = csr_from_coo(n,n,rows,columns,entries)
    end block
    do j = 1,3
        select case (j)
        case (1)
            choice = 'LA'
        case (2)
            choice = 'SA'
        case (3)
            choice = 'LM'
        end select
        call eigsh_lanczos(a,2,values,vectors,info,which=choice, &
            ncv=12,tol=1.0e-8_dp,maxiter=60,iterations=cycles, &
            locked_count=locked,residual_norms=errors,matvecs=products)
        call require(info == 0, 'restarted sparse Lanczos '//choice)
        call require(locked == 2 .and. cycles >= 1, 'Ritz convergence and soft locking')
        call require(products >= 12 .and. maxval(errors) < 2.0e-7_dp, &
            'Lanczos products and residuals')
        do i = 1,2
            y = csr_matvec(a,vectors(:,i))
            call require(abs(norm2(y-values(i)*vectors(:,i))-errors(i)) < 1.0e-10_dp, &
                'independently verified true residual')
            call require(abs(norm2(vectors(:,i))-1.0_dp) < 1.0e-9_dp, &
                'unit Ritz vectors')
        end do
    end do

    ! A slowly varying spectrum demands multiple restarts, rather than
    ! apparently converging in a single block expansion.
    block
        integer, parameter :: nn = 100
        integer :: rows(3*nn-2), cols(3*nn-2)
        real(dp) :: entries(3*nn-2)
        do i = 1,nn
            rows(i) = i
            cols(i) = i
            entries(i) = 2.0_dp+real(i,dp)/30.0_dp
        end do
        do i = 1,nn-1
            rows(nn+i) = i
            cols(nn+i) = i+1
            entries(nn+i) = -0.1_dp
            rows(2*nn-1+i) = i+1
            cols(2*nn-1+i) = i
            entries(2*nn-1+i) = -0.1_dp
        end do
        a = csr_from_coo(nn,nn,rows,cols,entries)
        call eigsh_lanczos(a,2,values,vectors,info,which='LA',ncv=20, &
            maxiter=60,tol=1.0e-8_dp,iterations=cycles, &
            locked_count=locked,residual_norms=errors)
        call require(info == 0 .and. cycles > 1 .and. locked == 2, &
            'multiple restart cycles for clustered spectral edge')
        call require(maxval(errors) < 6.0e-8_dp, 'clustered edge residuals')
        call eigsh_lanczos(a,2,values,vectors,info,which='LA',ncv=16, &
            maxiter=1,tol=1.0e-12_dp,residual_norms=errors)
        call require(info > 0 .and. maxval(errors) > 1.0e-9_dp, &
            'true nonconvergence reported')
    end block

    ! The optional method dispatcher must match direct Lanczos results.
    call eigsh(a,2,values,vectors,info,method='lanczos',which='LA', &
        ncv=20,maxiter=60,tol=1.0e-8_dp,residual_norms=errors)
    call require(info == 0 .and. maxval(errors) < 2.0e-7_dp, &
        'eigsh Lanczos method dispatch')
    call eigsh(a,2,values,vectors,info,method='nonsense')
    call require(info == sparse_linalg_invalid_input, 'invalid eigsh method')
    ! Shift-invert now works through the method dispatcher as well.
    matrix = 0.0_dp
    do i=1,6
        matrix(i,i) = sample(i)
    end do
    a = csr_from_dense(matrix)
    call eigsh(a,2,values,vectors,info,method='lanczos',sigma=2.4_dp, &
        which='LM',ncv=6,maxiter=10,tol=1.0e-8_dp)
    call require(info == 0, 'shifted Lanczos method dispatch')
    call close(values,[1.0_dp,3.0_dp],1.0e-8_dp,'shifted Lanczos values')
    call eigsh_lanczos(a,2,values,vectors,info,which='SM',ncv=12)
    call require(info == sparse_linalg_invalid_input, 'smallest magnitude requires full subspace')
    call eigsh_lanczos(a,2,values,vectors,info,ncv=2)
    call require(info == sparse_linalg_invalid_input, 'invalid Krylov width')
    call eigsh_lanczos(a,2,values,vectors,info,maxiter=0,iterations=cycles)
    call require(info > 0 .and. cycles == 0, 'zero iteration budget')
    call eigsh_lanczos(a,2,values,vectors,info,v0=[(0.0_dp,i=1,n)])
    call require(info == sparse_linalg_invalid_input, 'zero starting vector')
    nan_value = ieee_value(0.0_dp,ieee_quiet_nan)
    call eigsh_lanczos(a,2,values,vectors,info,tol=nan_value)
    call require(info == sparse_linalg_invalid_input, 'nonfinite tolerance')
    call eigsh_lanczos(a,2,values,vectors,info,which='ZZ')
    call require(info == sparse_linalg_invalid_input, 'unknown selection')
    a = csr_from_dense(reshape([1.0_dp,2.0_dp,0.0_dp,3.0_dp],[2,2]))
    call eigsh_lanczos(a,1,values,vectors,info)
    call require(info == sparse_linalg_invalid_input, 'asymmetric input')
    a = csr_matrix()
    call eigsh_lanczos(a,1,values,vectors,info)
    call require(info == sparse_linalg_invalid_input, 'invalid CSR')

    ! Existing nonsymmetric eigs can now opt in to restarted Arnoldi.
    matrix = 0.0_dp
    do i = 1,6
        matrix(i,i) = sample(i)
    end do
    a = csr_from_dense(matrix)
    call eigs(a,2,cvalues,cvectors,info,method='arnoldi',which='LM',ncv=6, &
        maxiter=10)
    call require(info == 0, 'eigs Arnoldi method dispatch')
    call eigs(a,2,cvalues,cvectors,info,method='bogus')
    call require(info == sparse_linalg_invalid_input, 'invalid eigs method')
    call eigs(a,1,cvalues,cvectors,info,method='arnoldi', &
        sigma=cmplx(2.6_dp,0.0_dp,dp),ncv=6,maxiter=10,tol=1.0e-8_dp)
    call require(info == 0, 'shifted Arnoldi method dispatch')
    call require(abs(cvalues(1)-cmplx(3.0_dp,0.0_dp,dp))<1.0e-8_dp, &
        'shifted Arnoldi values')

    print '(a)', 'test_sparse_lanczos: PASS'
contains
    subroutine require(condition, message)
        logical, intent(in) :: condition !! Condition that must be satisfied.
        character(len=*), intent(in) :: message !! Error message on failure.
        if (.not. condition) then
            print '(a)', 'FAILED: '//message
            error stop 1
        end if
    end subroutine require

    subroutine close(actual, expected, threshold, message)
        real(dp), intent(in) :: actual(:) !! Computed eigenvalues.
        real(dp), intent(in) :: expected(:) !! Expected eigenvalues.
        real(dp), intent(in) :: threshold !! Maximum permitted difference.
        character(len=*), intent(in) :: message !! Test-case name.
        call require(size(actual) == size(expected), 'shape '//message)
        call require(maxval(abs(actual-expected)) <= threshold, message)
    end subroutine close
end program test_sparse_lanczos
