! SPDX-License-Identifier: MIT
program test_spatial_kdtree_weighted_count
    use, intrinsic :: iso_fortran_env, only: int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use scifort_kinds, only: dp
    use scifort_spatial, only: kd_tree, kd_tree_type, kdtree_success, kdtree_invalid_input
    implicit none
    type(kd_tree_type) :: a, b, invalid
    real(dp), allocatable :: x(:,:), y(:,:), box(:), wa(:), wb(:), radii(:), weighted(:)
    real(dp), allocatable :: expected_weight(:), rv(:)
    integer(int64), allocatable :: counts(:), expected_count(:)
    real(dp) :: p, d, scalar_weight
    integer(int64) :: scalar_count
    integer :: dim,i,j,h,k,n,m,code,visited,nseed
    integer, allocatable :: seeds(:)
    logical :: cumulative
    real(dp), parameter :: orders(4)=[1.0_dp,1.5_dp,2.0_dp,3.0_dp]

    call random_seed(size=nseed)
    allocate(seeds(nseed))
    seeds=20261102
    call random_seed(put=seeds)
    do dim=1,7
        n=8+2*dim
        m=6+3*dim
        allocate(x(n,dim),y(m,dim),box(dim),wa(n),wb(m))
        call random_number(x)
        call random_number(y)
        call random_number(wa)
        call random_number(wb)
        wa=4.0_dp*wa-2.0_dp
        wb=3.0_dp*wb-1.5_dp
        x(2,:)=x(1,:)
        y(1,:)=x(1,:)
        box=0.0_dp
        do i=1,dim
            if (mod(i+dim,3)/=0) box(i)=1.0_dp
        end do
        a=kd_tree(x,boxsize=box,status=code)
        call require(code==kdtree_success,'construct source tree')
        b=kd_tree(y,boxsize=box,status=code)
        call require(code==kdtree_success,'construct target tree')
        do h=1,5
            if (h==5) then
                p=ieee_value(0.0_dp,ieee_positive_inf)
            else
                p=orders(h)
            end if
            do k=1,2
                cumulative=k==1
                if (cumulative) then
                    radii=[0.6_dp,0.0_dp,1.2_dp,0.31_dp,0.6_dp,4.5_dp,0.11_dp]
                else
                    radii=[0.0_dp,0.11_dp,0.31_dp,0.6_dp,0.6_dp,1.2_dp,4.5_dp]
                end if
                allocate(expected_count(size(radii)),expected_weight(size(radii)))
                expected_count=0_int64
                expected_weight=0.0_dp
                do i=1,n
                    do j=1,m
                        d=brute_lp(x(i,:),y(j,:),box,p)
                        call update_expected(d,wa(i)*wb(j),radii,cumulative,expected_count,expected_weight)
                    end do
                end do
                counts=a%count_neighbors(b,radii,p=p,cumulative=cumulative,status=code,visited=visited)
                call require(code==kdtree_success,'vector count status')
                call require(all(counts==expected_count),'vector integer counts')
                weighted=a%weighted_count_neighbors(b,radii,weights_self=wa,weights_other=wb, &
                    p=p,cumulative=cumulative,status=code,visited=visited)
                call require(code==kdtree_success,'vector weighted count status')
                call require(all(abs(weighted-expected_weight)<2.0e-11_dp),'vector weighted count values')
                ! One-sided weights are allowed; the other side defaults to unit weights.
                rv=a%weighted_count_neighbors(b,radii,weights_self=wa,p=p,cumulative=cumulative,status=code)
                call require(code==kdtree_success,'one-sided weighted count status')
                call require(all(abs(rv-brute_weights(x,y,box,p,radii,cumulative,wa=wa))<2.0e-11_dp), &
                    'one-sided source weighted count')
                rv=a%weighted_count_neighbors(b,radii,weights_other=wb,p=p,cumulative=cumulative,status=code)
                call require(code==kdtree_success,'one-sided target weighted count status')
                call require(all(abs(rv-brute_weights(x,y,box,p,radii,cumulative,wb=wb))<2.0e-11_dp), &
                    'one-sided destination weighted count')
                rv=a%weighted_count_neighbors(b,radii,p=p,cumulative=cumulative,status=code)
                call require(code==kdtree_success,'unit-weight count status')
                call require(all(abs(rv-real(expected_count,dp))<2.0e-11_dp),'unit-weight count')
                if (cumulative) then
                    do i=1,size(radii)
                        scalar_count=a%count_neighbors(b,radii(i),p=p,status=code)
                        call require(code==kdtree_success .and. scalar_count==expected_count(i), &
                            'scalar integer backwards compatibility')
                        scalar_weight=a%weighted_count_neighbors(b,radii(i),weights_self=wa, &
                            weights_other=wb,p=p,status=code)
                        call require(code==kdtree_success,'scalar weighted count status')
                        call require(abs(scalar_weight-expected_weight(i))<2.0e-11_dp, &
                            'scalar weighted counts')
                    end do
                end if
                deallocate(expected_count,expected_weight,radii)
            end do
        end do
        deallocate(x,y,box,wa,wb)
    end do

    a=kd_tree(reshape([0.0_dp,0.0_dp,1.0_dp],[3,1]))
    b=kd_tree(reshape([0.0_dp,2.0_dp],[2,1]))
    counts=a%count_neighbors(b,[0.0_dp,2.0_dp,2.0_dp],cumulative=.false.,status=code)
    call require(code==0 .and. all(counts==[2_int64,4_int64,0_int64]),'exact closed-ball boundaries')
    weighted=a%weighted_count_neighbors(b,[0.0_dp,2.0_dp],weights_self=[2.0_dp,-2.0_dp,1.0_dp], &
        weights_other=[3.0_dp,5.0_dp],status=code)
    call require(code==0 .and. all(abs(weighted-[0.0_dp,8.0_dp])<1.0e-13_dp), &
        'signed weights and cancellation')
    counts=a%count_neighbors(b,[real(dp)::],status=code)
    call require(code==0 .and. size(counts)==0,'empty radii')
    counts=a%count_neighbors(b,[1.0_dp,0.0_dp],cumulative=.false.,status=code)
    call require(code==kdtree_invalid_input .and. size(counts)==0,'reject unsorted histogram radii')
    weighted=a%weighted_count_neighbors(b,[0.0_dp,1.0_dp],weights_self=[1.0_dp],status=code)
    call require(code==kdtree_invalid_input .and. size(weighted)==0,'reject wrong weight dimension')
    weighted=a%weighted_count_neighbors(b,[0.0_dp,1.0_dp], &
        weights_other=[ieee_value(0.0_dp,ieee_quiet_nan),1.0_dp],status=code)
    call require(code==kdtree_invalid_input .and. size(weighted)==0,'reject nonfinite weights')
    scalar_weight=a%weighted_count_neighbors(b,-1.0_dp,status=code)
    call require(code==kdtree_invalid_input .and. scalar_weight==0.0_dp,'reject negative scalar radius')
    counts=a%count_neighbors(b,[0.0_dp,ieee_value(0.0_dp,ieee_quiet_nan)],status=code)
    call require(code==kdtree_invalid_input .and. size(counts)==0,'reject nonfinite radii')
    counts=a%count_neighbors(b,[0.0_dp,1.0_dp],p=0.5_dp,status=code)
    call require(code==kdtree_invalid_input .and. size(counts)==0,'reject p below one')
    invalid=kd_tree(reshape([0.0_dp,2.0_dp],[2,1]),boxsize=3.0_dp)
    counts=a%count_neighbors(invalid,[0.0_dp,1.0_dp],status=code)
    call require(code==kdtree_invalid_input .and. size(counts)==0,'reject mixed topologies')
    ! All-pairs pruning should avoid site distance evaluations.
    counts=a%count_neighbors(b,[1.0e6_dp,2.0e6_dp],visited=visited,status=code)
    call require(code==0 .and. visited==0,'full-subtree fast path')
    call require(all(counts==6_int64),'full-subtree count')
    ! For a tree queried against itself, include diagonal and both pair orders.
    a=kd_tree(reshape([0.0_dp,1.0_dp,2.0_dp],[3,1]))
    counts=a%count_neighbors(a,[0.0_dp,1.0_dp,2.0_dp],status=code)
    call require(code==0 .and. all(counts==[3_int64,7_int64,9_int64]), &
        'self-neighbor ordered counts')
    weighted=a%weighted_count_neighbors(a,[0.0_dp,1.0_dp,2.0_dp], &
        weights_self=[2.0_dp,3.0_dp,4.0_dp],weights_other=[2.0_dp,3.0_dp,4.0_dp],status=code)
    call require(code==0 .and. all(abs(weighted-[29.0_dp,65.0_dp,81.0_dp])<1.0e-13_dp), &
        'self-neighbor weighted diagonal and ordered pairs')
    ! Both the integer and weighted paths must avoid squaring tiny distances.
    a=kd_tree(reshape([0.0_dp,2.0e-200_dp],[2,1]))
    b=kd_tree(reshape([0.0_dp],[1,1]))
    counts=a%count_neighbors(b,[0.0_dp,1.0e-200_dp,3.0e-200_dp],status=code)
    call require(code==0 .and. all(counts==[1_int64,1_int64,2_int64]),'tiny-distance thresholds')
    print '(a)', 'PASS: vector and weighted KD-tree count_neighbors native regressions'
contains
    pure real(dp) function brute_lp(v,u,box,p) result(distance)
        real(dp), intent(in) :: v(:),u(:),box(:),p
        real(dp) :: gaps(size(v)),gap
        integer :: idx
        do idx=1,size(v)
            gap=abs(v(idx)-u(idx))
            if (box(idx)>0.0_dp) then
                gap=abs(modulo(v(idx),box(idx))-modulo(u(idx),box(idx)))
                gap=min(gap,box(idx)-gap)
            end if
            gaps(idx)=gap
        end do
        if (p==ieee_value(0.0_dp,ieee_positive_inf)) then
            distance=maxval(gaps)
        else
            distance=sum(gaps**p)**(1.0_dp/p)
        end if
    end function brute_lp

    subroutine update_expected(distance,pair_weight,thresholds,cumulative,counts,weights)
        real(dp), intent(in) :: distance,pair_weight,thresholds(:)
        logical, intent(in) :: cumulative
        integer(int64), intent(inout) :: counts(:)
        real(dp), intent(inout) :: weights(:)
        integer :: idx
        do idx=1,size(thresholds)
            if (distance>thresholds(idx)) cycle
            if (.not.cumulative) then
                if (idx>1) then
                    if (distance<=thresholds(idx-1)) cycle
                end if
            end if
            counts(idx)=counts(idx)+1_int64
            weights(idx)=weights(idx)+pair_weight
        end do
    end subroutine update_expected

    function brute_weights(x,y,box,p,radii,cumulative,wa,wb) result(result_weights)
        real(dp), intent(in) :: x(:,:),y(:,:),box(:),p,radii(:)
        logical, intent(in) :: cumulative
        real(dp), intent(in), optional :: wa(:),wb(:)
        real(dp) :: result_weights(size(radii)),weight,d
        integer(int64) :: unused(size(radii))
        integer :: i,j
        result_weights=0.0_dp
        unused=0_int64
        do i=1,size(x,1)
            do j=1,size(y,1)
                weight=1.0_dp
                if (present(wa)) weight=weight*wa(i)
                if (present(wb)) weight=weight*wb(j)
                d=brute_lp(x(i,:),y(j,:),box,p)
                call update_expected(d,weight,radii,cumulative,unused,result_weights)
            end do
        end do
    end function brute_weights

    subroutine require(ok,message)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: message
        if (.not.ok) then
            print *, 'FAIL: ',message
            error stop 1
        end if
    end subroutine require
end program test_spatial_kdtree_weighted_count
