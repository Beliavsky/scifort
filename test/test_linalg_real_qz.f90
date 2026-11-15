! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_linalg_real_qz
    use scifort_kinds, only: dp
    use scifort_linalg, only: linalg_status_invalid_shape,linalg_status_singular
    use scifort_linalg_qz, only: qz
    use scifort_linalg_ordqz, only: ordqz,ordqz_factors
    use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
    implicit none
    integer, parameter :: n=6
    real(dp) :: a(n,n),b(n,n),aa(n,n),bb(n,n),left(n,n),right(n,n),beta(n)
    real(dp) :: base_a(n,n),base_b(n,n),saved(n,n),angle,c,s,column(n),nan,inf
    real(dp) :: empty(0,0),eaa(0,0),ebb(0,0),eq(0,0),ez(0,0)
    complex(dp) :: alpha(n)
    integer :: i,j,code,sdim,mode,checks
    character(3), parameter :: sorts(4)=[character(3)::'lhp','rhp','iuc','ouc']
    integer, parameter :: counts(4)=[1,4,1,5]
    checks=0
    base_a=0.0_dp
    base_b=0.0_dp
    do i=1,n
        base_b(i,i)=1.0_dp
    end do
    base_a(1,1)=1.0_dp
    base_a(2,2)=1.0_dp
    base_a(1,2)=-2.0_dp
    base_a(2,1)=2.0_dp
    base_a(3,3)=-3.0_dp
    base_b(3,3)=2.0_dp
    base_a(4,4)=0.25_dp
    base_a(5,5)=2.0_dp
    base_b(5,5)=0.0_dp
    base_a(6,6)=7.0_dp
    left=base_b*0.0_dp
    right=left
    do i=1,n
        left(i,i)=1.0_dp
        right(i,i)=1.0_dp
    end do
    do i=1,n-1
        angle=0.17_dp*i
        c=cos(angle)
        s=sin(angle)
        column=left(:,i)
        left(:,i)=c*column+s*left(:,i+1)
        left(:,i+1)=-s*column+c*left(:,i+1)
        column=right(:,i)
        right(:,i)=c*column-s*right(:,i+1)
        right(:,i+1)=s*column+c*right(:,i+1)
    end do
    a=matmul(left,matmul(base_a,transpose(right)))
    b=matmul(left,matmul(base_b,transpose(right)))
    call qz(a,b,aa,bb,left,right,code,alpha,beta)
    call require(code==0,'real qz')
    call verify(a,b,aa,bb,left,right)
    call require(count(beta==0.0_dp)==1,'one infinite eigenvalue')
    call require(count(abs(aimag(alpha))>1e-12_dp)==2,'one conjugate pair')
    do mode=1,4
        call ordqz(a,b,aa,bb,left,right,code,alpha,beta,sort=trim(sorts(mode)),sdim=sdim)
        call require(code==0,'real ordered qz')
        call require(sdim==counts(mode),'whole-block selection count')
        call verify(a,b,aa,bb,left,right)
    end do
    call ordqz(a,b,aa,bb,left,right,code,alpha,beta,selector=positive_imag,sdim=sdim)
    call require(code==0 .and. sdim==2,'custom predicate closes conjugate pair')
    call require(aa(2,1)/=0.0_dp,'selected real two by two block')
    call verify(a,b,aa,bb,left,right)
    call qz(a,b,aa,bb,left,right,code,method='iterative')
    call require(code==0,'real iterative method')
    call ordqz_factors(aa,bb,left,right,code,sort='lhp',sdim=sdim)
    call require(code==0 .and. sdim==1,'real factor reordering')
    call verify(a,b,aa,bb,left,right)
    aa(n,1)=1.0_dp
    saved=aa
    call ordqz_factors(aa,bb,left,right,code)
    call require(code/=0 .and. all(aa==saved),'failed real reordering transactional')
    call qz(empty,empty,eaa,ebb,eq,ez,code)
    call require(code==0,'empty real pencil')
    call qz(a,b,aa,bb,left,right,code,method='inverse')
    call require(code==linalg_status_invalid_shape,'unsupported real method rejected')
    call qz(a(:,:n-1),b,aa,bb,left,right,code)
    call require(code==linalg_status_invalid_shape,'nonsquare real input rejected')
    call qz(a,b,aa,bb,left,right,code,alpha(:n-1))
    call require(code==linalg_status_invalid_shape,'wrong homogeneous shape rejected')
    call ordqz(a,b,aa,bb,left,right,code,sort='bad')
    call require(code==linalg_status_invalid_shape,'invalid real sort rejected')
    nan=ieee_value(0.0_dp,ieee_quiet_nan)
    inf=ieee_value(0.0_dp,ieee_positive_inf)
    saved=a
    a(1,1)=nan
    call qz(a,b,aa,bb,left,right,code)
    call require(code/=0,'NaN input rejected')
    a=saved
    a(1,1)=inf
    call qz(a,b,aa,bb,left,right,code)
    call require(code/=0,'infinite input rejected')
    a=0.0_dp
    b=0.0_dp
    call qz(a,b,aa,bb,left,right,code)
    call require(code==linalg_status_singular,'indeterminate pencil rejected')
    do i=1,n
        a(i,i)=real(i,dp)
    end do
    call qz(a,b,aa,bb,left,right,code,alpha,beta)
    call require(code==0 .and. all(beta==0.0_dp),'regular all-infinite pencil')
    call verify(a,b,aa,bb,left,right)
    a=0.0_dp
    b=0.0_dp
    do i=1,n
        b(i,i)=1.0_dp
    end do
    call qz(a,b,aa,bb,left,right,code,alpha,beta)
    call require(code==0 .and. all(abs(alpha)==0.0_dp),'regular all-zero eigenvalues')
    call verify(a,b,aa,bb,left,right)
    call ordqz(a,b,aa,bb,left,right,code,sort='iuc',sdim=sdim)
    call require(code==0 .and. sdim==n,'zero eigenvalues inside unit circle')
    call qz(a,b,aa,bb,left,right,code)
    left=1e200_dp
    saved=aa
    call ordqz_factors(aa,bb,left,right,code)
    call require(code/=0 .and. all(aa==saved),'nonorthogonal large basis rejected safely')
    print '(a,i0,a)', 'Real QZ and ordered QZ: ',checks,' checks passed'
contains
    subroutine require(ok,label)
        implicit none
        logical, intent(in) :: ok !! Required numerical or API condition.
        character(*), intent(in) :: label !! Diagnostic printed on failed assertion.
        checks=checks+1
        if (.not.ok) then
            print '(a)',label
            error stop 1
        end if
    end subroutine require
    logical function positive_imag(ar,be) result(selected)
        implicit none
        complex(dp), intent(in) :: ar !! Homogeneous eigenvalue numerator.
        complex(dp), intent(in) :: be !! Homogeneous denominator, unused by this predicate.
        selected=aimag(ar)>0.0_dp
    end function positive_imag
    subroutine verify(original_a,original_b,t,u,q,z)
        implicit none
        real(dp), intent(in) :: original_a(:,:) !! First original matrix.
        real(dp), intent(in) :: original_b(:,:) !! Second original matrix.
        real(dp), intent(in) :: t(:,:) !! Candidate real quasi-triangular first factor.
        real(dp), intent(in) :: u(:,:) !! Candidate upper-triangular second factor.
        real(dp), intent(in) :: q(:,:) !! Candidate left orthogonal basis.
        real(dp), intent(in) :: z(:,:) !! Candidate right orthogonal basis.
        real(dp) :: identity(n,n),sa,sb
        integer :: row,col
        identity=0.0_dp
        do row=1,n
            identity(row,row)=1.0_dp
        end do
        sa=max(1.0_dp,maxval(abs(original_a)))
        sb=max(1.0_dp,maxval(abs(original_b)))
        call require(maxval(abs(matmul(q,matmul(t/sa,transpose(z)))-original_a/sa))<1e-12_dp, &
            'first real reconstruction')
        call require(maxval(abs(matmul(q,matmul(u/sb,transpose(z)))-original_b/sb))<1e-12_dp, &
            'second real reconstruction')
        call require(maxval(abs(matmul(transpose(q),q)-identity))<1e-12_dp,'left orthogonality')
        call require(maxval(abs(matmul(transpose(z),z)-identity))<1e-12_dp,'right orthogonality')
        do col=1,n
            call require(u(col,col)>=0.0_dp,'nonnegative beta diagonal')
            do row=col+1,n
                call require(u(row,col)==0.0_dp,'triangular second factor')
                if (row>col+1) call require(t(row,col)==0.0_dp,'quasi-triangular first factor')
            end do
        end do
    end subroutine verify
end program test_linalg_real_qz
