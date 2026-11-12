! SPDX-License-Identifier: MIT
program test_sparse_eigs
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_dense, csr_from_coo
    use scifort_sparse_linalg, only : eigs, sparse_linalg_invalid_input
    implicit none
    type(csr_matrix) :: a
    real(dp) :: matrix(6,6), r, expected
    complex(dp), allocatable :: values(:), vectors(:,:)
    real(dp), allocatable :: errors(:)
    integer :: i, mode, info, steps, n, j
    character(len=2) :: target

    ! A real matrix with two complex-conjugate pairs and two real eigenvalues.
    matrix=0.0_dp
    matrix(1,2)=-3.0_dp
    matrix(2,1)=3.0_dp
    matrix(3,3)=8.0_dp
    matrix(4,4)=-6.0_dp
    matrix(5,6)=-1.0_dp
    matrix(6,5)=1.0_dp
    a=csr_from_dense(matrix)
    do mode=1,6
        select case(mode)
        case(1)
            target='LM'
        case(2)
            target='LR'
        case(3)
            target='SR'
        case(4)
            target='LI'
        case(5)
            target='SI'
        case(6)
            target='SM'
        end select
        call eigs(a,2,values,vectors,info,which=target,tol=1.0e-10_dp, &
            ncv=6,maxiter=10,residual_norms=errors,iterations=steps)
        call require(info==0, 'small eigs '//target)
        call require(size(values)==2 .and. size(vectors,1)==6 .and. &
            size(vectors,2)==2 .and. steps==1,'small result shape')
        call require(maxval(errors)<1.0e-9_dp,'small residual '//target)
        do i=1,2
            r=sqrt(sum(abs(matmul(matrix,vectors(:,i))-values(i)*vectors(:,i))**2))
            call require(r<1.0e-9_dp,'small true residual')
            call require(abs(sqrt(sum(abs(vectors(:,i))**2))-1.0_dp)<1.0e-10_dp, &
                'unit complex eigenvector')
        end do
    end do
    ! Explicit validation and zero-update behavior.
    call eigs(a,2,values,vectors,info,which='XX')
    call require(info==sparse_linalg_invalid_input,'invalid which')
    call eigs(a,0,values,vectors,info)
    call require(info==sparse_linalg_invalid_input,'zero eigenvalue request')
    call eigs(a,6,values,vectors,info)
    call require(info==sparse_linalg_invalid_input,'k greater than n')
    call eigs(a,2,values,vectors,info,tol=-1.0_dp)
    call require(info==sparse_linalg_invalid_input,'negative tolerance')
    call eigs(a,2,values,vectors,info,ncv=2)
    call require(info==sparse_linalg_invalid_input,'ncv too small')
    call eigs(a,2,values,vectors,info,v0=[cmplx(1.0_dp,0.0_dp,dp)])
    call require(info==sparse_linalg_invalid_input,'short initial vector')
    call eigs(a,2,values,vectors,info,v0= &
        [(cmplx(0.0_dp,0.0_dp,dp),i=1,6)])
    call require(info==sparse_linalg_invalid_input,'zero initial vector')
    call eigs(a,2,values,vectors,info,v0= &
        [(cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp),i=1,6)])
    call require(info==sparse_linalg_invalid_input,'nonfinite initial vector')
    call eigs(a,2,values,vectors,info,which='SM',ncv=5,maxiter=160, &
        tol=1.0e-7_dp,residual_norms=errors)
    call require(info==0,'SM reduced subspace shift-invert')
    call require(maxval(errors)<1.0e-5_dp,'SM reduced subspace residuals')
    call require(maxval(abs(abs(values)-1.0_dp))<1.0e-5_dp, &
        'SM reduced subspace nearest pair')
    call eigs(a,2,values,vectors,info,maxiter=0)
    call require(info>0 .and. size(values)==2,'zero update budget')
    a=csr_matrix()
    call eigs(a,1,values,vectors,info)
    call require(info==sparse_linalg_invalid_input,'invalid sparse matrix')
    matrix=0.0_dp
    a=csr_from_dense(matrix)
    call eigs(a,3,values,vectors,info,which='SM',ncv=6, &
        tol=1.0e-12_dp,residual_norms=errors)
    call require(info==0,'zero nonsymmetric matrix')
    call require(all(abs(values)<1.0e-11_dp) .and. &
        all(errors<1.0e-11_dp),'zero matrix eigenpairs')

    ! Reduced-subspace tests on a sparse 60-by-60 nonnormal matrix.
    ! Four dominant real/complex eigenvalues are separated from background.
    n=60
    do mode=1,5
        select case(mode)
        case(1)
            target='LM'
        case(2)
            target='LR'
        case(3)
            target='SR'
        case(4)
            target='LI'
        case(5)
            target='SI'
        end select
        block
            real(dp) :: diag(n), off(n-1), lower(n-1)
            diag=[(0.45_dp*sin(real(i,dp)),i=1,n)]
            off=0.01_dp
            lower=-0.01_dp
            if (mode<=2) then
                diag(1)=8.0_dp
                diag(2)=7.0_dp
            else if (mode==3) then
                diag(1)=-8.0_dp
                diag(2)=-7.0_dp
            else
                diag(1:4)=0.0_dp
                off(1)=-8.0_dp
                lower(1)=8.0_dp
                off(3)=-7.0_dp
                lower(3)=7.0_dp
            end if
            a=csr_from_coo(n,n,[(i,i=1,n),(i,i=1,n-1),(i,i=2,n)], &
                [(i,i=1,n),(i+1,i=1,n-1),(i-1,i=2,n)], &
                [diag,off,lower])
        end block
        call eigs(a,2,values,vectors,info,which=target, &
            ncv=20,tol=2.0e-7_dp,maxiter=200,residual_norms=errors,iterations=steps)
        call require(info==0,'reduced-subspace eigs '//target)
        call require(maxval(errors)<2.0e-6_dp,'reduced-subspace true residual')
        if (mode<=2) then
            call require(abs(real(values(1),dp)-8.0_dp)<0.01_dp,'dominant real eigs')
        else if (mode==3) then
            call require(abs(real(values(1),dp)+8.0_dp)<0.01_dp,'smallest real eigs')
        else if (mode==4) then
            call require(abs(aimag(values(1))-8.0_dp)<0.01_dp,'largest imaginary eigs')
        else
            call require(abs(aimag(values(1))+8.0_dp)<0.01_dp,'smallest imaginary eigs')
        end if
    end do
    print '(a)','test_sparse_eigs: PASS'
contains
    subroutine require(condition, message)
        logical, intent(in) :: condition !! Required invariant.
        character(len=*), intent(in) :: message !! Error explanation.
        if (.not. condition) then
            print '(a)', 'FAILED: '//message
            error stop 1
        end if
    end subroutine require
end program test_sparse_eigs
