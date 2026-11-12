! SPDX-License-Identifier: MIT
program optimize_linear_least_squares_reference
    use scifort_kinds, only : dp
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_quiet_nan
    use scifort_optimize_linear_least_squares
    implicit none
    real(dp) :: a(3,2), b(3), x(2), rnorm, inf
    real(dp) :: rank1(3,2), bounds_lo(2), bounds_hi(2)
    integer :: status, checks
    type(linear_lsq_result) :: res
    checks=0
    a(:,1)=[1.0_dp,1.0_dp,0.0_dp]
    a(:,2)=[0.0_dp,0.0_dp,1.0_dp]
    b=[2.0_dp,1.0_dp,1.0_dp]
    call nnls(a,b,x,rnorm,status)
    call must(status==0,'nnls status')
    call must(maxval(abs(x-[1.5_dp,1.0_dp]))<1.0e-12_dp,'nnls documented vector')
    call must(abs(rnorm-sqrt(0.5_dp))<1.0e-12_dp,'nnls residual')
    b = [-1.0_dp,-1.0_dp,-1.0_dp]
    call nnls(a,b,x,rnorm,status)
    call must(status==0,'nnls negative b status')
    call must(maxval(abs(x))==0.0_dp,'nnls all zero')
    call must(abs(rnorm-sqrt(3.0_dp))<1.0e-12_dp,'nnls negative residual')
    inf=ieee_value(1.0_dp,ieee_positive_inf)
    b=[2.0_dp,1.0_dp,1.0_dp]
    bounds_lo=[0.0_dp,0.0_dp]
    bounds_hi=[1.0_dp,inf]
    res=lsq_linear(a,b,bounds_lo,bounds_hi)
    call must(res%success,'bvls converges')
    call must(maxval(abs(res%x-[1.0_dp,1.0_dp]))<1.0e-12_dp,'upper active')
    call must(res%active_mask(1)==1,'upper mask')
    call must(res%optimality<1.0e-9_dp,'small KKT')
    bounds_lo=[1.75_dp,0.2_dp]
    bounds_hi=[3.0_dp,0.6_dp]
    res=lsq_linear(a,b,bounds_lo,bounds_hi)
    call must(res%success,'bvls two bounds converges')
    call must(maxval(abs(res%x-[1.75_dp,0.6_dp]))<1.0e-12_dp,'two active bounds')
    rank1(:,1)=[1.0_dp,2.0_dp,3.0_dp]
    rank1(:,2)=2.0_dp*rank1(:,1)
    b=rank1(:,1)
    res=lsq_linear(rank1,b)
    call must(res%success,'rank deficient unbounded')
    call must(norm2(res%fun)<1.0e-11_dp,'rank deficient residual')
    res=lsq_linear(rank1,b,lower=[0.0_dp,0.0_dp])
    call must(res%success,'rank deficient bounded')
    call must(all(res%x>=0.0_dp),'rank deficient positive')
    res=lsq_linear(a,b,method='trf')
    call must(res%status==linear_lsq_invalid,'unsupported trf rejection')
    res=lsq_linear(a,b,lower=[2.0_dp,0.0_dp],upper=[1.0_dp,1.0_dp])
    call must(res%status==linear_lsq_invalid,'reversed bounds reject')
    res=lsq_linear(a,b,lower=[0.0_dp])
    call must(res%status==linear_lsq_invalid,'wrong bound size reject')
    res=lsq_linear(a,b,lower=[0.0_dp,0.0_dp],maxiter=-1)
    call must(res%status==linear_lsq_invalid,'negative maxiter reject')
    b(1)=ieee_value(1.0_dp,ieee_quiet_nan)
    res=lsq_linear(a,b)
    call must(res%status==linear_lsq_invalid,'nan input rejection')
    res=lsq_linear(a,[2.0_dp,1.0_dp,1.0_dp], &
        lower=[1.2_dp,-1.0_dp], upper=[1.2_dp,2.0_dp])
    call must(res%success,'fixed variable optimum')
    call must(abs(res%x(1)-1.2_dp)<1.0e-12_dp,'fixed variable exact')
    res=lsq_linear(a,[2.0_dp,1.0_dp,1.0_dp], &
        lower=[0.0_dp,0.0_dp],upper=[1.0_dp,1.0_dp],maxiter=0)
    call must(res%status==linear_lsq_maxiter .or. res%success, &
        'iteration limit reported or already optimal')
    res=lsq_linear(a,[2.0_dp,1.0_dp,1.0_dp],tol=-1.0_dp)
    call must(res%status==linear_lsq_invalid,'negative tolerance rejected')
    res=lsq_linear(a,[2.0_dp,1.0_dp,1.0_dp], &
        upper=[0.0_dp,ieee_value(1.0_dp,ieee_positive_inf)])
    call must(res%success .and. abs(res%x(1))<1.0e-12_dp, 'zero upper bound')
    res=lsq_linear(a, [2.0_dp,1.0_dp,1.0_dp], &
        lower=[-ieee_value(1.0_dp,ieee_positive_inf), -1.0_dp])
    call must(res%success,'minus infinity lower bound')
    res=lsq_linear(reshape([real(dp)::],[0,2]),[real(dp)::], &
        lower=[0.0_dp,1.0_dp],upper=[1.0_dp,2.0_dp])
    call must(res%success .and. size(res%fun)==0,'empty rows')
    call must(maxval(abs(res%x-[0.0_dp,1.0_dp]))<1.0e-12_dp,'empty-row feasible zero')
    res=lsq_linear(reshape([real(dp)::],[3,0]),[1.0_dp,2.0_dp,3.0_dp])
    call must(res%success .and. size(res%x)==0,'empty columns')
    call must(abs(res%cost-7.0_dp)<1.0e-12_dp,'empty-column cost')
    res=lsq_linear(a,[2.0_dp,1.0_dp,1.0_dp],upper= &
        [1.0_dp,ieee_value(1.0_dp,ieee_quiet_nan)])
    call must(res%status==linear_lsq_invalid,'nan bound rejected')
    write(*,'(a,1x,i0)') 'PASSED',checks
contains
    subroutine must(condition,description)
        logical, intent(in) :: condition !! test condition
        character(len=*), intent(in) :: description !! assertion label
        if (.not.condition) then
            print *, 'FAILED: ',description
            error stop 1
        end if
        checks=checks+1
    end subroutine must
end program optimize_linear_least_squares_reference
