! SPDX-License-Identifier: MIT
program test_spatial_kdtree_pairs
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_spatial, only : kd_tree_type, kd_tree, kdtree_success, kdtree_invalid_input
    implicit none
    type(kd_tree_type) :: first, second, bad
    real(dp), allocatable :: x(:,:), y(:,:), distances(:)
    real(dp) :: radius, d
    integer, allocatable :: pairs(:,:), row_ptr(:), indices(:), rows(:), columns(:), reference_pairs(:,:)
    integer, allocatable :: seed(:)
    integer(int64) :: actual_count, expected_count
    integer :: dim, n, m, i, j, k, h, nseed, matched, code, visits, pairs_seen

    call random_seed(size=nseed)
    allocate(seed(nseed))
    seed=20261101
    call random_seed(put=seed)
    do dim=1,7
        n=34+dim*7
        m=28+dim*5
        allocate(x(n,dim),y(m,dim),reference_pairs(n*(n-1)/2,2))
        call random_number(x)
        call random_number(y)
        x(2,:)=x(1,:) ! identical sites must remain distinct entries
        y(1,:)=x(1,:)
        first=kd_tree(x,status=code)
        call require(code==kdtree_success,'first tree construction')
        second=kd_tree(y,status=code)
        call require(code==kdtree_success,'second tree construction')
        do h=1,5
            select case (h)
            case(1); radius=0.0_dp
            case(2); radius=0.12_dp
            case(3); radius=0.37_dp
            case(4); radius=0.8_dp
            case(5); radius=real(dim,dp)**0.5_dp+0.5_dp
            end select
            call first%query_pairs(radius,pairs,status=code,visited=visits)
            call require(code==kdtree_success,'query_pairs status')
            pairs_seen=0
            do i=1,n-1
                do j=i+1,n
                    d=norm2(x(i,:)-x(j,:))
                    if (d<=radius) then
                        pairs_seen=pairs_seen+1
                        reference_pairs(pairs_seen,:)=[i,j]
                    end if
                end do
            end do
            call require(size(pairs,1)==pairs_seen .and. size(pairs,2)==2,'query_pairs cardinality')
            call require(all(pairs==reference_pairs(:pairs_seen,:)),'query_pairs exact sorted indices')
            call first%query_ball_tree(second,radius,row_ptr,indices,status=code,visited=visits)
            call require(code==kdtree_success,'query_ball_tree status')
            call require(size(row_ptr)==n+1 .and. row_ptr(1)==1,'CSR row pointer shape')
            call require(row_ptr(n+1)==size(indices)+1,'CSR final pointer')
            expected_count=0_int64
            do i=1,n
                matched=0
                do j=1,m
                    if (norm2(x(i,:)-y(j,:))<=radius) then
                        matched=matched+1
                        call require(indices(row_ptr(i)+matched-1)==j,'query_ball_tree CSR indices')
                    end if
                end do
                call require(row_ptr(i+1)-row_ptr(i)==matched,'query_ball_tree CSR row count')
                expected_count=expected_count+int(matched,int64)
            end do
            actual_count=first%count_neighbors(second,radius,status=code,visited=visits)
            call require(code==kdtree_success .and. actual_count==expected_count, 'count_neighbors cross trees')
            call first%sparse_distance_matrix(second,radius,rows,columns,distances,status=code)
            call require(code==kdtree_success,'sparse_distance_matrix status')
            call require(size(rows)==size(indices) .and. size(columns)==size(indices),'sparse matrix cardinality')
            call require(size(distances)==size(indices),'sparse matrix distance count')
            k=0
            do i=1,n
                do j=row_ptr(i),row_ptr(i+1)-1
                    k=k+1
                    call require(rows(k)==i .and. columns(k)==indices(j),'sparse matrix row-major order')
                    d=norm2(x(i,:)-y(columns(k),:))
                    call require(abs(distances(k)-d)<5.0e-15_dp,'sparse matrix ordinary distances')
                end do
            end do
        end do
        actual_count=first%count_neighbors(first,0.0_dp,status=code)
        call require(code==kdtree_success .and. actual_count>=int(n+2,int64),'self includes diagonal and duplicates')
        actual_count=first%count_neighbors(first,real(dim,dp)**0.5_dp+1.0_dp,status=code,visited=visits)
        call require(code==kdtree_success .and. actual_count==int(n,int64)**2,'self counts ordered pairs')
        call require(visits==0,'box-wide count is pruned as a whole')
        deallocate(x,y,reference_pairs)
    end do

    first=kd_tree(reshape([0.0_dp,1.0_dp,2.0_dp,3.0_dp],[4,1]))
    call first%query_pairs(1.0_dp,pairs,status=code)
    call require(code==0 .and. size(pairs,1)==3,'query_pairs boundary-inclusive')
    call require(all(pairs(:,1)==[1,2,3]) .and. all(pairs(:,2)==[2,3,4]), 'query_pairs sorted')
    call first%query_pairs(-1.0_dp,pairs,status=code)
    call require(code==kdtree_invalid_input .and. size(pairs,1)==0,'pairs reject negative radius')
    call first%query_pairs(ieee_value(0.0_dp,ieee_quiet_nan),pairs,status=code)
    call require(code==kdtree_invalid_input .and. size(pairs,1)==0,'pairs reject NaN radius')
    second=kd_tree(reshape([0.0_dp,0.0_dp],[1,2]))
    call first%query_ball_tree(second,1.0_dp,row_ptr,indices,status=code)
    call require(code==kdtree_invalid_input,'cross reject dimensional mismatch')
    call require(size(row_ptr)==0 .and. size(indices)==0,'cross empty arrays on error')
    actual_count=first%count_neighbors(second,1.0_dp,status=code)
    call require(code==kdtree_invalid_input .and. actual_count==0_int64,'count reject mismatched dims')
    call first%sparse_distance_matrix(second,1.0_dp,rows,columns,distances,status=code)
    call require(code==kdtree_invalid_input .and. size(rows)==0,'sparse reject mismatched dims')
    ! Regression: squaring tiny radii underflows and squaring huge radii
    ! overflows, so closed-ball comparisons must use a scaled Euclidean norm.
    first=kd_tree(reshape([0.0_dp,2.0e-200_dp],[2,1]))
    call first%query_pairs(1.0e-200_dp,pairs,status=code)
    call require(code==0 .and. size(pairs,1)==0,'avoid underflow in tiny query_pairs')
    actual_count=first%count_neighbors(first,1.0e-200_dp,status=code)
    call require(code==0 .and. actual_count==2_int64,'avoid underflow in tiny pair counting')
    first=kd_tree(reshape([0.0_dp,8.0e200_dp],[2,1]))
    call first%query_pairs(9.0e200_dp,pairs,status=code)
    call require(code==0 .and. size(pairs,1)==1,'avoid overflow in large query_pairs')
    actual_count=first%count_neighbors(first,9.0e200_dp,status=code)
    call require(code==0 .and. actual_count==4_int64,'avoid overflow in large counting')
    call first%query_pairs(7.0e200_dp,pairs,status=code)
    call require(code==0 .and. size(pairs,1)==0,'exclude outside large radius')
    ! Each row of an identity-like cross query must produce the expected diagonal.
    second=kd_tree(reshape([0.0_dp,8.0e200_dp],[2,1]))
    call first%query_ball_tree(second,0.0_dp,row_ptr,indices,status=code)
    call require(code==0 .and. all(row_ptr==[1,2,3]) .and. all(indices==[1,2]),'zero-radius CSR diagonal')
    call first%sparse_distance_matrix(second,0.0_dp,rows,columns,distances,status=code)
    call require(code==0 .and. all(distances==0.0_dp),'sparse explicit zeros')
    call first%query_pairs(ieee_value(0.0_dp,ieee_quiet_nan),pairs,status=code)
    call require(code==kdtree_invalid_input,'reject NaN radius')
    bad=kd_tree(reshape([0.0_dp],[0,1]))
    actual_count=bad%count_neighbors(first,1.0_dp,status=code)
    call require(code==kdtree_invalid_input,'count reject invalid tree')
    print '(a)', 'PASS: spatial KD-tree paired and cross-tree native regressions'
contains
    subroutine require(ok,message)
        logical, intent(in) :: ok !! desired invariant
        character(len=*), intent(in) :: message !! diagnostic label
        if (.not.ok) then
            print '(a,a)', 'FAIL: ',message
            error stop 1
        end if
    end subroutine require
end program test_spatial_kdtree_pairs
