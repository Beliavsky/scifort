! SPDX-License-Identifier: MIT
program test_sparse_complex_direct
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        csr_complex_from_dense, csr_complex_to_csc, csr_complex_from_coo
    use scifort_sparse_linalg_complex_direct, only : splu_complex_factor, &
        splu_complex, valid_splu_complex, splu_complex_solve, spsolve_complex, &
        complex_direct_invalid_input, complex_direct_singular
    implicit none

    type(csr_complex_matrix) :: a, bad
    type(csc_complex_matrix) :: csc
    type(splu_complex_factor) :: fac
    complex(dp) :: d(4,4), rhs(4,3), known(4,3), zero(0,0), singular(2,2), tiny(1,1)
    complex(dp), allocatable :: x(:), y(:,:), b(:,:), matrix(:,:)
    integer :: ierr, i, j, k, n, count
    integer, allocatable :: rows(:), cols(:)
    complex(dp), allocatable :: values(:)
    character(len=1) :: trans

    d = cmplx(0.0_dp,0.0_dp,dp)
    d(1,:) = [cmplx(0.0_dp,0.0_dp,dp), cmplx(2.0_dp,1.0_dp,dp), &
              cmplx(0.0_dp,0.0_dp,dp), cmplx(1.0_dp,0.0_dp,dp)]
    d(2,:) = [cmplx(3.0_dp,-1.0_dp,dp), cmplx(5.0_dp,1.0_dp,dp), &
              cmplx(0.5_dp,0.0_dp,dp), cmplx(0.0_dp,0.0_dp,dp)]
    d(3,:) = [cmplx(0.0_dp,0.0_dp,dp), cmplx(-1.0_dp,0.3_dp,dp), &
              cmplx(4.0_dp,0.5_dp,dp), cmplx(2.0_dp,0.0_dp,dp)]
    d(4,:) = [cmplx(0.2_dp,-0.1_dp,dp), cmplx(0.0_dp,0.0_dp,dp), &
              cmplx(1.0_dp,1.0_dp,dp), cmplx(6.0_dp,0.5_dp,dp)]
    do j = 1,3
        do i = 1,4
            known(i,j) = cmplx(real(i*j,dp)/3.0_dp,real(i-j,dp)/7.0_dp,dp)
        end do
    end do
    a = csr_complex_from_dense(d,ierr)
    call require(ierr == 0,'dense to complex CSR')
    csc = csr_complex_to_csc(a,ierr)
    call require(ierr == 0,'CSR to CSC')
    call splu_complex(a,fac,ierr)
    call require(ierr == 0 .and. valid_splu_complex(fac),'factor CSR')
    call require(any(fac%row_perm /= [1,2,3,4]),'row pivot required')
    do k = 1,3
        select case (k)
        case (1)
            trans = 'N'
            rhs = matmul(d,known)
        case (2)
            trans = 'T'
            rhs = matmul(transpose(d),known)
        case (3)
            trans = 'H'
            rhs = matmul(transpose(conjg(d)),known)
        end select
        call splu_complex_solve(fac,rhs,y,ierr,trans)
        call require(ierr == 0,'factor multi RHS solve')
        call require(maxval(abs(y-known)) < 2.0e-12_dp,'factor multi RHS values')
        call splu_complex_solve(fac,rhs(:,1),x,ierr,trans)
        call require(ierr == 0,'factor vector solve')
        call require(maxval(abs(x-known(:,1))) < 2.0e-12_dp,'factor vector values')
    end do
    call splu_complex(csc,fac,ierr)
    call require(ierr == 0 .and. valid_splu_complex(fac),'factor CSC')
    call splu_complex_solve(fac,rhs(:,1),x,ierr,trans='H')
    call require(ierr == 0 .and. maxval(abs(x-known(:,1))) < 2.0e-12_dp,'CSC adjoint')
    rhs = matmul(d,known)
    call spsolve_complex(a,rhs,y,ierr)
    call require(ierr == 0 .and. maxval(abs(y-known)) < 2.0e-12_dp,'CSR spsolve block')
    call spsolve_complex(csc,rhs,y,ierr)
    call require(ierr == 0 .and. maxval(abs(y-known)) < 2.0e-12_dp,'CSC spsolve block')
    call spsolve_complex(a,rhs(:,1),x,ierr)
    call require(ierr == 0 .and. maxval(abs(x-known(:,1))) < 2.0e-12_dp,'CSR spsolve vector')
    call spsolve_complex(csc,rhs(:,1),x,ierr)
    call require(ierr == 0 .and. maxval(abs(x-known(:,1))) < 2.0e-12_dp,'CSC spsolve vector')

    call splu_complex_solve(fac,rhs(:,1),x,ierr,trans='X')
    call require(ierr == complex_direct_invalid_input .and. size(x) == 0,'reject bad transpose')
    call splu_complex_solve(fac,rhs(:,1),x,ierr,trans='nn')
    call require(ierr == complex_direct_invalid_input,'reject overlong transpose')
    call splu_complex_solve(fac,rhs(:3,1),x,ierr)
    call require(ierr == complex_direct_invalid_input,'reject RHS shape')
    rhs(1,1) = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
    call splu_complex_solve(fac,rhs(:,1),x,ierr)
    call require(ierr == complex_direct_invalid_input,'reject NaN RHS')
    rhs = matmul(d,known)
    call splu_complex(a,fac,ierr)
    call require(ierr == 0,'factor for edge-case solves')
    allocate(b(4,0))
    call splu_complex_solve(fac,b,y,ierr,trans='T')
    call require(ierr == 0 .and. all(shape(y) == [4,0]),'empty RHS block transpose')
    deallocate(b)
    fac%row_perm(1) = fac%row_perm(2)
    call require(.not. valid_splu_complex(fac),'reject invalid permutation')
    call splu_complex_solve(fac,rhs(:,1),x,ierr)
    call require(ierr == complex_direct_invalid_input .and. size(x) == 0, &
        'reject invalid factor')
    call splu_complex(a,fac,ierr,pivot_tol=-1.0_dp)
    call require(ierr == complex_direct_invalid_input .and. .not. valid_splu_complex(fac), &
        'reject negative pivot tolerance')
    call splu_complex(a,fac,ierr,pivot_tol=10.0_dp)
    call require(ierr == complex_direct_singular .and. .not. valid_splu_complex(fac), &
        'report pivot threshold')

    singular = cmplx(0.0_dp,0.0_dp,dp)
    singular(1,:) = [cmplx(0.0_dp,0.0_dp,dp),cmplx(1.0_dp,1.0_dp,dp)]
    singular(2,:) = [cmplx(0.0_dp,0.0_dp,dp),cmplx(2.0_dp,1.0_dp,dp)]
    bad = csr_complex_from_dense(singular)
    call splu_complex(bad,fac,ierr)
    call require(ierr == complex_direct_singular,'missing pivot singular')
    call spsolve_complex(bad,[cmplx(1.0_dp,0.0_dp,dp),cmplx(2.0_dp,0.0_dp,dp)],x,ierr)
    call require(ierr == complex_direct_singular .and. size(x) == 0,'singular spsolve')
    bad = csr_complex_from_coo(2,2,[1,1,2,2],[1,2,1,2], &
        [cmplx(0.0_dp,0.0_dp,dp),cmplx(2.0_dp,1.0_dp,dp), &
         cmplx(1.0_dp,-1.0_dp,dp),cmplx(3.0_dp,0.0_dp,dp)],ierr)
    call require(ierr == 0,'explicit zero CSR')
    call splu_complex(bad,fac,ierr)
    call require(ierr == 0,'pivot on explicit zero')
    call splu_complex_solve(fac,[cmplx(2.0_dp,1.0_dp,dp),cmplx(5.0_dp,-1.0_dp,dp)],x,ierr)
    call require(ierr == 0,'explicit zero solve')

    bad = csr_complex_from_dense(zero,ierr)
    call require(ierr == 0,'empty CSR construction')
    call splu_complex(bad,fac,ierr)
    call require(ierr == 0 .and. valid_splu_complex(fac),'empty matrix factor')
    call splu_complex_solve(fac,[complex(dp)::],x,ierr)
    call require(ierr == 0 .and. size(x) == 0,'empty RHS solve')
    call spsolve_complex(bad,zero,y,ierr)
    call require(ierr == 0 .and. size(y,1) == 0,'empty block solve')

    tiny(1,1) = cmplx(2.0_dp,3.0_dp,dp)
    a = csr_complex_from_dense(tiny)
    call spsolve_complex(a,[cmplx(4.0_dp,-2.0_dp,dp)],x,ierr)
    call require(ierr == 0,'singleton solve')
    call require(abs(x(1)*tiny(1,1)-cmplx(4.0_dp,-2.0_dp,dp)) < 1.0e-14_dp, &
        'singleton accuracy')

    ! A sparse tridiagonal example exercises fill storage without dense LU.
    n = 256
    allocate(rows(3*n-2),cols(3*n-2),values(3*n-2))
    count = 0
    do i = 1,n
        count = count+1
        rows(count) = i
        cols(count) = i
        values(count) = cmplx(4.0_dp,0.25_dp,dp)
        if (i < n) then
            count = count+1
            rows(count) = i
            cols(count) = i+1
            values(count) = cmplx(-1.0_dp,-0.2_dp,dp)
            count = count+1
            rows(count) = i+1
            cols(count) = i
            values(count) = cmplx(-1.0_dp,0.3_dp,dp)
        end if
    end do
    a = csr_complex_from_coo(n,n,rows,cols,values,ierr)
    call require(ierr == 0,'banded sparse assembly')
    call splu_complex(a,fac,ierr)
    call require(ierr == 0,'banded sparse factor')
    call require(size(fac%lu%data) < 6*n,'banded LU remains sparse')
    allocate(matrix(n,2))
    matrix(:,1) = cmplx(1.0_dp,0.3_dp,dp)
    matrix(:,2) = cmplx(-1.0_dp,0.7_dp,dp)
    allocate(b(n,2))
    ! Form b using sparse row products, not a dense n-by-n matrix.
    b = cmplx(0.0_dp,0.0_dp,dp)
    do i = 1,n
        do k = a%indptr(i),a%indptr(i+1)-1
            b(i,:) = b(i,:) + a%data(k)*matrix(a%indices(k),:)
        end do
    end do
    call splu_complex_solve(fac,b,y,ierr)
    call require(ierr == 0 .and. maxval(abs(y-matrix)) < 1.0e-12_dp,'banded sparse solve')
    print '(a)', 'PASS complex sparse direct LU/spsolve, transpose/adjoint, multi-RHS and failures'
contains
    subroutine require(condition,label)
        logical, intent(in) :: condition !! Assertion condition.
        character(len=*), intent(in) :: label !! Failure label.
        if (.not. condition) then
            print '(a)', 'FAIL '//label
            error stop 1
        end if
    end subroutine require
end program test_sparse_complex_direct
