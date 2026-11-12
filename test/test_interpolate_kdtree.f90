! SPDX-License-Identifier: MIT
program test_interpolate_kdtree
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use scifort_kinds, only: dp
    use scifort_interpolate_kdtree, only: kd_tree_type, kd_tree, kdtree_success, kdtree_invalid_input
    implicit none
    type(kd_tree_type) :: tree, bad
    integer, allocatable :: ids(:), brute_ids(:), ball(:), expected_ball(:), seed(:)
    real(dp), allocatable :: sites(:,:), queries(:,:), distances(:), brute_dist(:)
    integer :: n,dim,q,k,i,j,h,code,visits,total_visits,nseed
    real(dp) :: r2,rad

    bad=kd_tree(reshape([0.0_dp,ieee_value(0.0_dp,ieee_quiet_nan)],[2,1]),status=code)
    call require(code==kdtree_invalid_input .and. .not.bad%valid(),'reject nonfinite nodes')
    bad=kd_tree(reshape([0.0_dp],[0,1]),status=code)
    call require(code==kdtree_invalid_input .and. .not.bad%valid(),'reject empty input')
    tree=kd_tree(reshape([0.0_dp,0.0_dp, 0.0_dp,0.0_dp, &
        0.0_dp,0.0_dp, 3.0_dp,4.0_dp],[4,2]))
    call tree%query([0.0_dp,0.0_dp],ids,distance2=distances,k=4,status=code)
    call require(code==0 .and. all(ids==[1,2,3,4]),'all duplicate sites, stable index tie')
    call tree%query_ball_point([0.0_dp,0.0_dp],5.0_dp,ball,status=code)
    call require(code==0 .and. all(ball==[1,2,3,4]),'inclusive radius boundary')
    call tree%query([0.0_dp,0.0_dp],ids,k=5,status=code)
    call require(code==kdtree_invalid_input,'reject k exceeding point count')
    call random_seed(size=nseed)
    allocate(seed(nseed))
    seed=20261031
    call random_seed(put=seed)
    do dim=1,7
        n=129+11*dim
        allocate(sites(n,dim),queries(35,dim),brute_dist(n),brute_ids(n))
        call random_number(sites)
        call random_number(queries)
        sites(1,:)=0.5_dp
        sites(2,:)=0.5_dp ! verify deterministic tie handling on duplicates
        tree=kd_tree(sites,status=code)
        call require(tree%valid().and.code==kdtree_success,'valid construction')
        call require(tree%dimensions()==dim .and. tree%size()==n,'index dimensions')
        do q=1,size(queries,1)
            if (q==1) queries(q,:)=0.5_dp
            do i=1,n
                brute_ids(i)=i
                brute_dist(i)=sum((sites(i,:)-queries(q,:))**2)
            end do
            do i=2,n
                j=i
                do while (j>1)
                    if (.not.is_worse(brute_ids(j-1),brute_ids(j),brute_dist)) exit
                    h=brute_ids(j)
                    brute_ids(j)=brute_ids(j-1)
                    brute_ids(j-1)=h
                    j=j-1
                end do
            end do
            do k=1,5
                select case(k)
                case(1);h=1
                case(2);h=2
                case(3);h=11
                case(4);h=min(50,n)
                case(5);h=n
                end select
                call tree%query(queries(q,:),ids,distance2=distances,k=h,status=code,visited=visits)
                call require(code==kdtree_success,'k-neighbor query status')
                call require(size(ids)==h .and. all(ids==brute_ids(:h)),'k-nearest indices')
                do i=1,h
                    call require(abs(distances(i)-brute_dist(ids(i)))<1.0e-14_dp,'k-nearest distances')
                end do
                call require(visits<=n .and. visits>=h,'query visit accounting')
            end do
            call require(tree%nearest(queries(q,:))==brute_ids(1),'nearest call')
            rad=0.18_dp
            r2=rad*rad
            call tree%query_ball_point(queries(q,:),rad,ball,status=code,visited=visits)
            call require(code==kdtree_success,'radius query status')
            expected_ball=pack([(i,i=1,n)],brute_dist<=r2)
            call require(all(ball==expected_ball),'radius query matches brute force')
        end do
        deallocate(sites,queries,brute_dist,brute_ids)
    end do
    allocate(sites(4096,2),queries(256,2))
    call random_number(sites)
    call random_number(queries)
    tree=kd_tree(sites)
    total_visits=0
    do i=1,size(queries,1)
        call tree%query(queries(i,:),ids,k=8,visited=visits)
        total_visits=total_visits+visits
    end do
    call require(real(total_visits,dp)/size(queries,1)<0.12_dp*size(sites,1), 'search must prune subtrees')
    call tree%query([0.5_dp],ids,status=code)
    call require(code==kdtree_invalid_input .and. size(ids)==0,'reject wrong query dimension')
    call tree%query([0.5_dp,0.5_dp],ids,k=0,status=code)
    call require(code==kdtree_invalid_input .and. size(ids)==0,'reject zero neighbors')
    call tree%query_ball_point([0.5_dp,0.5_dp],-1.0_dp,ball,status=code)
    call require(code==kdtree_invalid_input .and. size(ball)==0,'reject negative radius')
    print '(a,i0)', 'PASS: balanced KD-tree native regressions; avg nodes visited=',total_visits/256
contains
    pure logical function is_worse(a,b,d) result(worse)
        integer, intent(in) :: a !! first site index
        integer, intent(in) :: b !! second site index
        real(dp), intent(in) :: d(:) !! squared distance by site index
        worse=d(a)>d(b) .or. (d(a)==d(b) .and. a>b)
    end function is_worse
    subroutine require(ok,reason)
        logical, intent(in) :: ok !! predicate
        character(len=*), intent(in) :: reason !! diagnostic
        if (.not.ok) then
            print *, 'FAIL: ',reason
            error stop 1
        end if
    end subroutine require
end program test_interpolate_kdtree
