! SPDX-License-Identifier: MIT
program test_spatial_kdtree_periodic
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_spatial, only : kd_tree_type, kd_tree, kdtree_success, kdtree_invalid_input
    implicit none
    type(kd_tree_type) :: a, b, invalid
    real(dp), allocatable :: x(:,:), y(:,:), box(:), query_point(:), dist(:), d2(:), expected(:), vals(:)
    real(dp), parameter :: orders(4)=[1.0_dp,1.5_dp,2.0_dp,3.0_dp]
    integer, allocatable :: ids(:), ball(:), pairs(:,:), rp(:), columns(:), rows(:), cols(:), order(:)
    integer, allocatable :: seed(:)
    integer(int64) :: count, brute_count
    integer :: nseed,dim,i,j,k,h,t,code,m,n,expected_pairs,nb,idx
    real(dp) :: p,radius

    call random_seed(size=nseed)
    allocate(seed(nseed))
    seed=20261101
    call random_seed(put=seed)
    do dim=1,6
        n=28+dim
        m=15+dim
        allocate(x(n,dim),y(m,dim),box(dim),query_point(dim),order(n),expected(n))
        call random_number(x)
        call random_number(y)
        box=0.0_dp
        do i=1,dim
            if (mod(i+dim,3)/=0) box(i)=1.0_dp
        end do
        x(2,:)=x(1,:)
        y(1,:)=x(1,:)
        a=kd_tree(x,boxsize=box,status=code)
        call require(code==kdtree_success .and. a%valid(),'periodic constructor')
        b=kd_tree(y,boxsize=box,status=code)
        call require(code==kdtree_success,'periodic other constructor')
        do h=1,5
            if (h<=4) then
                p=orders(h)
            else
                p=ieee_value(0.0_dp,ieee_positive_inf)
            end if
            call random_number(query_point)
            query_point=query_point*3.0_dp-1.0_dp ! outside periodic box must wrap
            do i=1,n
                expected(i)=brute_lp(query_point,x(i,:),box,p)
                order(i)=i
            end do
            do i=2,n
                j=i
                do while(j>1)
                    if (expected(order(j-1))<expected(order(j))) exit
                    if (expected(order(j-1))==expected(order(j)) .and. order(j-1)<order(j)) exit
                    idx=order(j)
                    order(j)=order(j-1)
                    order(j-1)=idx
                    j=j-1
                end do
            end do
            call a%query(query_point,ids,k=6,p=p,distance=dist,distance2=d2,status=code)
            call require(code==0 .and. all(ids==order(:6)),'periodic Minkowski nearest indices')
            call require(all(abs(dist-expected(ids))<2.0e-14_dp),'periodic Minkowski nearest distances')
            call require(all(abs(d2-dist**2)<2.0e-14_dp),'squared-distance backward compatibility')
            call require(a%nearest(query_point,p=p)==order(1),'periodic nearest')
            do t=1,3
                select case(t)
                case(1); radius=0.0_dp
                case(2); radius=0.31_dp
                case(3); radius=real(dim,dp)+1.0_dp
                end select
                call a%query_ball_point(query_point,radius,ball,p=p,status=code)
                call require(code==0,'periodic ball point status')
                nb=0
                do i=1,n
                    if (expected(i)<=radius) then
                        nb=nb+1
                        call require(ball(nb)==i,'periodic ball point indices')
                    end if
                end do
                call require(size(ball)==nb,'periodic ball point size')
                call a%query_pairs(radius,pairs,p=p,status=code)
                call require(code==0,'periodic query_pairs status')
                expected_pairs=0
                do i=1,n-1
                    do j=i+1,n
                        if (brute_lp(x(i,:),x(j,:),box,p)>radius) cycle
                        expected_pairs=expected_pairs+1
                        call require(all(pairs(expected_pairs,:)==[i,j]),'periodic pairs ordered')
                    end do
                end do
                call require(size(pairs,1)==expected_pairs,'periodic pairs count')
                call a%query_ball_tree(b,radius,rp,columns,p=p,status=code)
                call require(code==0 .and. size(rp)==n+1,'periodic cross query')
                brute_count=0_int64
                do i=1,n
                    nb=0
                    do j=1,m
                        if (brute_lp(x(i,:),y(j,:),box,p)>radius) cycle
                        nb=nb+1
                        call require(columns(rp(i)+nb-1)==j,'periodic cross column')
                    end do
                    call require(rp(i+1)-rp(i)==nb,'periodic cross row')
                    brute_count=brute_count+int(nb,int64)
                end do
                count=a%count_neighbors(b,radius,p=p,status=code)
                call require(code==0 .and. count==brute_count,'periodic Minkowski count')
                call a%sparse_distance_matrix(b,radius,rows,cols,vals,p=p,status=code)
                call require(code==0 .and. size(vals)==size(columns),'periodic sparse size')
                k=0
                do i=1,n
                    do j=rp(i),rp(i+1)-1
                        k=k+1
                        call require(rows(k)==i .and. cols(k)==columns(j),'periodic sparse ordering')
                        call require(abs(vals(k)-brute_lp(x(i,:),y(cols(k),:),box,p))<2.0e-14_dp, &
                            'periodic sparse distance')
                    end do
                end do
            end do
        end do
        deallocate(x,y,box,query_point,order,expected)
    end do
    a=kd_tree(reshape([0.05_dp,0.95_dp],[2,1]),boxsize=[1.0_dp],status=code)
    call a%query_ball_point([1.02_dp],0.08_dp,ball,status=code)
    call require(code==0 .and. all(ball==[1,2]),'periodic wrapped boundary radius')
    call a%query_pairs(0.099_dp,pairs,status=code)
    call require(code==0 .and. size(pairs,1)==0,'closed radius against floating rounding')
    call a%query_pairs(0.101_dp,pairs,status=code)
    call require(code==0 .and. size(pairs,1)==1,'periodic across edge pair')
    call a%query([2.98_dp],ids,k=2,distance=dist,status=code)
    call require(code==0 .and. all(ids==[2,1]),'periodic query outside box')
    ! The scalar constructor must broadcast to all coordinate dimensions.
    a=kd_tree(reshape([0.1_dp,0.9_dp,0.9_dp,0.1_dp],[2,2]),boxsize=1.0_dp,status=code)
    call require(code==0 .and. a%dimensions()==2,'scalar boxsize overload')
    call a%query_pairs(0.29_dp,pairs,status=code)
    call require(code==0 .and. size(pairs,1)==1,'scalar periodic box across two axes')
    invalid=kd_tree(reshape([0.0_dp,1.0_dp],[2,1]),boxsize=[1.0_dp],status=code)
    call require(code==kdtree_invalid_input .and. .not.invalid%valid(),'reject coordinate at L')
    invalid=kd_tree(reshape([-0.1_dp,0.9_dp],[2,1]),boxsize=[1.0_dp],status=code)
    call require(code==kdtree_invalid_input,'reject negative periodic coordinate')
    invalid=kd_tree(reshape([0.1_dp,0.9_dp],[2,1]),boxsize=[-1.0_dp],status=code)
    call require(code==kdtree_invalid_input,'reject negative box length')
    invalid=kd_tree(reshape([0.1_dp,0.9_dp],[2,1]),boxsize=[1.0_dp,1.0_dp],status=code)
    call require(code==kdtree_invalid_input,'reject box shape mismatch')
    invalid=kd_tree(reshape([0.1_dp,0.9_dp],[2,1]),boxsize=[ieee_value(0.0_dp,ieee_quiet_nan)],status=code)
    call require(code==kdtree_invalid_input,'reject nonfinite box length')
    a=kd_tree(reshape([0.05_dp,0.95_dp],[2,1]),boxsize=1.0_dp)
    b=kd_tree(reshape([0.05_dp],[1,1]))
    count=a%count_neighbors(b,0.1_dp,status=code)
    call require(code==kdtree_invalid_input,'reject mismatched periodic boxes')
    call a%query_pairs(1.0_dp,pairs,p=0.5_dp,status=code)
    call require(code==kdtree_invalid_input .and. size(pairs,1)==0,'reject Minkowski p below 1')
    call a%query([0.5_dp],ids,p=ieee_value(0.0_dp,ieee_quiet_nan),status=code)
    call require(code==kdtree_invalid_input .and. size(ids)==0,'reject NaN Minkowski p')
    call a%query_ball_tree(b,1.0_dp,rp,columns,p=0.5_dp,status=code)
    call require(code==kdtree_invalid_input,'reject bad p and cross topology')
    a=kd_tree(reshape([0.0_dp,2.0e-200_dp],[2,1]))
    call a%query([0.0_dp],ids,k=2,distance=dist,status=code)
    call require(code==0 .and. all(ids==[1,2]),'stable default nearest ordering')
    call require(dist(2)==2.0e-200_dp,'avoid underflow in unsquared nearest distance')
    a=kd_tree(reshape([0.0_dp,8.0e200_dp],[2,1]))
    call a%query([0.0_dp],ids,k=2,distance=dist,status=code)
    call require(code==0 .and. dist(2)==8.0e200_dp,'avoid overflow in unsquared nearest distance')
    print '(a)', 'PASS: periodic/mixed-topology Minkowski KD-tree native regressions'
contains
    pure real(dp) function brute_lp(a,b,box,p) result(result_distance)
        real(dp), intent(in) :: a(:),b(:),box(:),p
        real(dp) :: gaps(size(a)),v
        integer :: i
        do i=1,size(a)
            v=abs(a(i)-b(i))
            if (box(i)>0.0_dp) then
                v=abs(modulo(a(i),box(i))-modulo(b(i),box(i)))
                v=min(v,box(i)-v)
            end if
            gaps(i)=v
        end do
        if (p==ieee_value(0.0_dp,ieee_positive_inf)) then
            result_distance=maxval(gaps)
        else
            result_distance=sum(gaps**p)**(1.0_dp/p)
        end if
    end function brute_lp
    subroutine require(ok,msg)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: msg
        if (.not.ok) then
            print *, 'FAIL: ',msg
            error stop 1
        end if
    end subroutine require
end program test_spatial_kdtree_periodic
