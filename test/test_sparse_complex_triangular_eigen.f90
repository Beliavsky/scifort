! SPDX-License-Identifier: MIT
program test_sparse_complex_triangular_eigen
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        csr_complex_from_dense, csr_complex_from_coo, csr_complex_to_csc, csr_complex_matvec
    use scifort_sparse_linalg_complex_triangular, only : spsolve_triangular_complex, &
        complex_triangular_invalid_input, complex_triangular_bad_diagonal
    use scifort_sparse_linalg_complex_eigen, only : eigsh_complex, complex_eigen_invalid_input
    implicit none
    type(csr_complex_matrix) :: a, bad, hermitian
    type(csc_complex_matrix) :: ac
    complex(dp) :: d(4,4), rhs(4,2), truth(4,2), hh(2,2), empty(0,0), repeated(6,6)
    complex(dp), allocatable :: sol(:,:), vec(:), vectors(:,:), mv(:)
    real(dp), allocatable :: values(:), residuals(:)
    integer :: info, i, j, mode, steps
    character(len=1) :: trans

    d = cmplx(0.0_dp,0.0_dp,dp)
    d(1,1) = cmplx(2.0_dp,1.0_dp,dp)
    d(2,2) = cmplx(-3.0_dp,0.4_dp,dp)
    d(3,3) = cmplx(4.0_dp,-0.5_dp,dp)
    d(4,4) = cmplx(2.5_dp,0.0_dp,dp)
    d(2,1) = cmplx(1.0_dp,0.75_dp,dp)
    d(3,1) = cmplx(-0.5_dp,1.5_dp,dp)
    d(3,2) = cmplx(2.0_dp,-0.2_dp,dp)
    d(4,1) = cmplx(0.25_dp,0.0_dp,dp)
    d(4,3) = cmplx(1.0_dp,2.0_dp,dp)
    do j = 1,2
        do i = 1,4
            truth(i,j) = cmplx(real(i+j,dp),real(2*i-j,dp)/3.0_dp,dp)
        end do
    end do
    do mode = 1,2
        if (mode == 2) d = transpose(d)
        a = csr_complex_from_dense(d,info)
        call require(info == 0,'triangular csr constructor')
        ac = csr_complex_to_csc(a,info)
        call require(info == 0,'triangular csc constructor')
        do j = 1,3
            select case(j)
            case(1)
                trans = 'N'
                rhs = matmul(d,truth)
            case(2)
                trans = 'T'
                rhs = matmul(transpose(d),truth)
            case(3)
                trans = 'H'
                rhs = matmul(conjg(transpose(d)),truth)
            end select
            call spsolve_triangular_complex(a,rhs,sol,info,lower=(mode==1),trans=trans)
            call require(info == 0,'csr triangular block')
            call require(maxval(abs(sol-truth)) < 3.0e-13_dp,'csr triangular solve accuracy')
            call spsolve_triangular_complex(ac,rhs,sol,info,lower=(mode==1),trans=trans)
            call require(info == 0,'csc triangular block')
            call require(maxval(abs(sol-truth)) < 3.0e-13_dp,'csc triangular solve accuracy')
            call spsolve_triangular_complex(a,rhs(:,1),vec,info,lower=(mode==1),trans=trans)
            call require(info == 0,'csr triangular vector')
            call require(maxval(abs(vec-truth(:,1))) < 3.0e-13_dp,'triangular vector accuracy')
        end do
    end do
    d = cmplx(0.0_dp,0.0_dp,dp)
    d(2,1) = cmplx(1.0_dp,0.5_dp,dp)
    d(4,2) = cmplx(0.25_dp,-1.0_dp,dp)
    a = csr_complex_from_dense(d)
    rhs = matmul(d+identity4(),truth)
    call spsolve_triangular_complex(a,rhs,sol,info,unit_diagonal=.true.)
    call require(info == 0,'unit diagonal lower')
    call require(maxval(abs(sol-truth)) < 1.0e-13_dp,'unit diagonal accuracy')
    call spsolve_triangular_complex(a,rhs,sol,info)
    call require(info == complex_triangular_bad_diagonal,'missing zero pivot rejected')
    call require(size(sol) == 0,'triangular output empty on pivot failure')
    d(1,4) = cmplx(0.25_dp,0.0_dp,dp)
    a = csr_complex_from_dense(d)
    call spsolve_triangular_complex(a,rhs,sol,info,unit_diagonal=.true.)
    call require(info == complex_triangular_invalid_input,'nontriangular rejected')
    call spsolve_triangular_complex(a,rhs(:,1),vec,info,trans='X')
    call require(info == complex_triangular_invalid_input,'invalid operation rejected')
    a = csr_complex_from_dense(empty)
    call spsolve_triangular_complex(a,empty,sol,info)
    call require(info == 0 .and. size(sol) == 0,'empty triangular solve')

    hh(1,:) = [cmplx(2.0_dp,0.0_dp,dp),cmplx(1.0_dp,1.0_dp,dp)]
    hh(2,:) = [cmplx(1.0_dp,-1.0_dp,dp),cmplx(3.0_dp,0.0_dp,dp)]
    hermitian = csr_complex_from_dense(hh)
    ac = csr_complex_to_csc(hermitian)
    call eigsh_complex(hermitian,1,values,vectors,info,which='LA',ncv=2, &
        tol=1.0e-12_dp,residual_norms=residuals,iterations=steps)
    call require(info == 0 .and. steps == 2,'Hermitian LA convergence')
    call require(abs(values(1)-4.0_dp) < 1.0e-12_dp,'largest complex eigenvalue')
    call require(residuals(1) < 1.0e-12_dp,'Hermitian LA residual')
    call eigsh_complex(ac,1,values,vectors,info,which='SA',ncv=2,tol=1.0e-12_dp)
    call require(info == 0,'CSC Hermitian SA convergence')
    call require(abs(values(1)-1.0_dp) < 1.0e-12_dp,'smallest complex eigenvalue')
    hh(1,2) = cmplx(1.0_dp,1.25_dp,dp)
    hermitian = csr_complex_from_dense(hh)
    call eigsh_complex(hermitian,1,values,vectors,info)
    call require(info == complex_eigen_invalid_input,'nonhermitian rejected')
    call eigsh_complex(ac,0,values,vectors,info)
    call require(info == complex_eigen_invalid_input,'bad k rejected')
    call eigsh_complex(ac,1,values,vectors,info,ncv=1)
    call require(info == complex_eigen_invalid_input,'bad ncv rejected')
    call eigsh_complex(ac,1,values,vectors,info,which='BE')
    call require(info == complex_eigen_invalid_input,'unsupported mode rejected')

    repeated = cmplx(0.0_dp,0.0_dp,dp)
    do i = 1,6
        if (i <= 2) then
            repeated(i,i) = cmplx(3.0_dp,0.0_dp,dp)
        else if (i <= 4) then
            repeated(i,i) = cmplx(0.0_dp,0.0_dp,dp)
        else
            repeated(i,i) = cmplx(-2.0_dp,0.0_dp,dp)
        end if
    end do
    hermitian = csr_complex_from_dense(repeated)
    call eigsh_complex(hermitian,2,values,vectors,info,which='LA',ncv=6, &
        tol=1.0e-11_dp,residual_norms=residuals)
    call require(info == 0,'repeated top eigenvalues')
    call require(maxval(abs(values-3.0_dp)) < 1.0e-10_dp,'repeat multiplicity')
    call require(abs(dot_product(vectors(:,1),vectors(:,2))) < 1.0e-10_dp,'repeat orthogonality')
    call eigsh_complex(hermitian,2,values,vectors,info,which='SM',ncv=6,tol=1.0e-11_dp)
    call require(info == 0,'repeated smallest-magnitude eigenvalues')
    call require(maxval(abs(values)) < 1.0e-10_dp,'repeated zero eigenvalues')

    block
        integer, parameter :: n = 80
        integer :: r(3*n-2), c(3*n-2), p
        complex(dp) :: dval(3*n-2)
        p = 0
        do i = 1,n
            p = p+1
            r(p) = i
            c(p) = i
            dval(p) = cmplx(2.0_dp+real(i,dp)/40.0_dp,0.0_dp,dp)
            if (i < n) then
                p = p+1
                r(p) = i
                c(p) = i+1
                dval(p) = cmplx(-0.1_dp,0.04_dp,dp)
                p = p+1
                r(p) = i+1
                c(p) = i
                dval(p) = cmplx(-0.1_dp,-0.04_dp,dp)
            end if
        end do
        hermitian = csr_complex_from_coo(n,n,r,c,dval,info)
        call require(info == 0,'moderate Hermitian construction')
        call eigsh_complex(hermitian,3,values,vectors,info,which='LA',ncv=n, &
            tol=1.0e-9_dp,residual_norms=residuals)
        call require(info == 0,'moderate Hermitian convergence')
        call require(maxval(residuals) < 1.0e-8_dp,'moderate true residuals')
        do j = 1,3
            mv = csr_complex_matvec(hermitian,vectors(:,j))
            call require(norm2(abs(mv-values(j)*vectors(:,j))) < 1.0e-8_dp, &
                'moderate Hermitian eigenpairs')
        end do
    end block

    print *, 'PASS complex triangular and Hermitian eigsh'
contains
    function identity4() result(a)
        complex(dp) :: a(4,4)
        integer :: k
        a = cmplx(0.0_dp,0.0_dp,dp)
        do k = 1,4
            a(k,k) = cmplx(1.0_dp,0.0_dp,dp)
        end do
    end function identity4

    subroutine require(condition,message)
        logical, intent(in) :: condition !! Condition expected true.
        character(len=*), intent(in) :: message !! Diagnostic assertion label.
        if (.not. condition) then
            print *, 'FAIL: ',message
            error stop 1
        end if
    end subroutine require
end program test_sparse_complex_triangular_eigen
