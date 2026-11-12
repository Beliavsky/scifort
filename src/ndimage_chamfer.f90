! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Discrete city-block and chessboard distances to background pixels.
! Independent integer separable transforms; SciPy used only for validation.
module scifort_ndimage_chamfer
    use, intrinsic :: iso_fortran_env, only : int64
    implicit none
    private
    public :: distance_transform_cdt_nd, distance_transform_weighted_cdt_nd
    integer, parameter :: inf_distance=huge(1)/8
contains
    subroutine pack_mask(x,flat)
        logical, intent(in) :: x(..) !! logical foreground mask, rank 1 through 15
        logical, allocatable, intent(out) :: flat(:) !! Fortran-order flattened mask
        select rank(x)
        rank(1)
            flat=reshape(x,[size(x)])
        rank(2)
            flat=reshape(x,[size(x)])
        rank(3)
            flat=reshape(x,[size(x)])
        rank(4)
            flat=reshape(x,[size(x)])
        rank(5)
            flat=reshape(x,[size(x)])
        rank(6)
            flat=reshape(x,[size(x)])
        rank(7)
            flat=reshape(x,[size(x)])
        rank(8)
            flat=reshape(x,[size(x)])
        rank(9)
            flat=reshape(x,[size(x)])
        rank(10)
            flat=reshape(x,[size(x)])
        rank(11)
            flat=reshape(x,[size(x)])
        rank(12)
            flat=reshape(x,[size(x)])
        rank(13)
            flat=reshape(x,[size(x)])
        rank(14)
            flat=reshape(x,[size(x)])
        rank(15)
            flat=reshape(x,[size(x)])
        rank default
            allocate(flat(0))
        end select
    end subroutine pack_mask

    subroutine scatter_integer(flat,y)
        integer, intent(in) :: flat(:) !! packed integer distances
        integer, intent(out) :: y(..) !! output with shape of mask
        select rank(y)
        rank(1)
            y=reshape(flat,shape(y))
        rank(2)
            y=reshape(flat,shape(y))
        rank(3)
            y=reshape(flat,shape(y))
        rank(4)
            y=reshape(flat,shape(y))
        rank(5)
            y=reshape(flat,shape(y))
        rank(6)
            y=reshape(flat,shape(y))
        rank(7)
            y=reshape(flat,shape(y))
        rank(8)
            y=reshape(flat,shape(y))
        rank(9)
            y=reshape(flat,shape(y))
        rank(10)
            y=reshape(flat,shape(y))
        rank(11)
            y=reshape(flat,shape(y))
        rank(12)
            y=reshape(flat,shape(y))
        rank(13)
            y=reshape(flat,shape(y))
        rank(14)
            y=reshape(flat,shape(y))
        rank(15)
            y=reshape(flat,shape(y))
        rank default
        end select
    end subroutine scatter_integer

    subroutine taxicab_line(f,g,winner,step)
        integer, intent(in) :: f(:) !! initial distance to a background site
        integer, intent(out) :: g(:) !! one-axis L1 distances
        integer, intent(out) :: winner(:) !! original position of best site
        integer, intent(in) :: step !! strictly positive axis step cost
        integer :: i,n
        n=size(f)
        g=f
        do i=1,n
            winner(i)=i
            if (f(i)==inf_distance) winner(i)=0
        end do
        do i=2,n
            if (g(i-1)>=inf_distance-step) cycle
            if (g(i-1)+step<g(i)) then
                g(i)=g(i-1)+step
                winner(i)=winner(i-1)
            end if
        end do
        do i=n-1,1,-1
            if (g(i+1)>=inf_distance-step) cycle
            if (g(i+1)+step<g(i)) then
                g(i)=g(i+1)+step
                winner(i)=winner(i+1)
            end if
        end do
    end subroutine taxicab_line

    pure subroutine better_node(a,ia,b,ib,v,iv)
        integer, intent(in) :: a,ia,b,ib !! pair of minimum costs and source indices
        integer, intent(out) :: v,iv !! least cost; lowest index on a tie
        if (a<b .or. (a==b .and. ia<=ib)) then
            v=a
            iv=ia
        else
            v=b
            iv=ib
        end if
    end subroutine better_node

    subroutine range_min(tree,site,base,lo,hi,v,source)
        integer, intent(in) :: tree(:),site(:) !! binary minimum segment tree
        integer, intent(in) :: base,lo,hi !! one-based inclusive query range
        integer, intent(out) :: v,source !! cost and winning source position
        integer :: l,r,z,iz
        v=inf_distance
        source=huge(1)
        l=base+lo-1
        r=base+hi-1
        do while(l<=r)
            if (modulo(l,2)==1) then
                call better_node(v,source,tree(l),site(l),z,iz)
                v=z
                source=iz
                l=l+1
            end if
            if (modulo(r,2)==0) then
                call better_node(v,source,tree(r),site(r),z,iz)
                v=z
                source=iz
                r=r-1
            end if
            l=l/2
            r=r/2
        end do
        if (v==inf_distance) source=0
    end subroutine range_min

    subroutine chessboard_line(f,g,winner)
        integer, intent(in) :: f(:) !! prior L-infinity distances in other axes
        integer, intent(out) :: g(:) !! transformed L-infinity distances
        integer, intent(out) :: winner(:) !! winning coordinate along the axis
        integer, allocatable :: tree(:),site(:)
        integer :: base,n,i,level,v,ix,low,high,mid,maxfinite,z,iz

        n=size(f)
        if (all(f==inf_distance)) then
            g=inf_distance
            winner=0
            return
        end if
        base=1
        do while(base<n)
            base=base*2
        end do
        allocate(tree(2*base),site(2*base))
        tree=inf_distance
        site=huge(1)
        do i=1,n
            tree(base+i-1)=f(i)
            site(base+i-1)=i
        end do
        do level=base-1,1,-1
            call better_node(tree(2*level),site(2*level), &
                tree(2*level+1),site(2*level+1),z,iz)
            tree(level)=z
            site(level)=iz
        end do
        maxfinite=maxval(f,mask=f/=inf_distance)
        do i=1,n
            low=-1
            high=max(maxfinite,n-1)
            do while(high-low>1)
                mid=low+(high-low)/2
                call range_min(tree,site,base,max(1,i-mid),min(n,i+mid),v,ix)
                if (v<=mid) then
                    high=mid
                else
                    low=mid
                end if
            end do
            g(i)=high
            call range_min(tree,site,base,max(1,i-high),min(n,i+high),v,ix)
            winner(i)=ix
        end do
    end subroutine chessboard_line

    subroutine pack_weight_stencil(x,flat)
        integer, intent(in) :: x(..) !! N-dimensional 3-wide nonnegative cost stencil
        integer, allocatable, intent(out) :: flat(:) !! Fortran-order stencil weights
        select rank(x)
        rank(1)
            flat=reshape(x,[size(x)])
        rank(2)
            flat=reshape(x,[size(x)])
        rank(3)
            flat=reshape(x,[size(x)])
        rank(4)
            flat=reshape(x,[size(x)])
        rank(5)
            flat=reshape(x,[size(x)])
        rank(6)
            flat=reshape(x,[size(x)])
        rank(7)
            flat=reshape(x,[size(x)])
        rank(8)
            flat=reshape(x,[size(x)])
        rank(9)
            flat=reshape(x,[size(x)])
        rank(10)
            flat=reshape(x,[size(x)])
        rank(11)
            flat=reshape(x,[size(x)])
        rank(12)
            flat=reshape(x,[size(x)])
        rank(13)
            flat=reshape(x,[size(x)])
        rank(14)
            flat=reshape(x,[size(x)])
        rank(15)
            flat=reshape(x,[size(x)])
        rank default
            allocate(flat(0))
        end select
    end subroutine pack_weight_stencil

    pure logical function heap_precedes(a,b,dist,source) result(before)
        integer, intent(in) :: a,b !! packed linear pixel IDs
        integer(int64), intent(in) :: dist(:) !! shortest path distances
        integer, intent(in) :: source(:) !! winning background sites
        before=dist(a)<dist(b)
        if (dist(a)==dist(b)) then
            before=source(a)<source(b)
            if (source(a)==source(b)) before=a<b
        end if
    end function heap_precedes

    subroutine heap_decrease(node,heap,locations,used,dist,source)
        integer, intent(in) :: node !! inserted/decreased node
        integer, intent(inout) :: heap(:),locations(:),used !! binary min-heap
        integer(int64), intent(in) :: dist(:) !! priority distances
        integer, intent(in) :: source(:) !! priority tie-break sources
        integer :: i,parent,other

        i=locations(node)
        if (i==0) then
            used=used+1
            i=used
        end if
        do while(i>1)
            parent=i/2
            other=heap(parent)
            if (.not.heap_precedes(node,other,dist,source)) exit
            heap(i)=other
            locations(other)=i
            i=parent
        end do
        heap(i)=node
        locations(node)=i
    end subroutine heap_decrease

    subroutine heap_pop(node,heap,locations,used,dist,source)
        integer, intent(out) :: node !! removed lowest-priority pixel
        integer, intent(inout) :: heap(:),locations(:),used !! binary heap
        integer(int64), intent(in) :: dist(:) !! distances
        integer, intent(in) :: source(:) !! winning source IDs
        integer :: last,i,left,child

        node=heap(1)
        locations(node)=0
        last=heap(used)
        used=used-1
        if (used==0) return
        i=1
        do
            left=2*i
            if (left>used) exit
            child=left
            if (left<used) then
                if (heap_precedes(heap(left+1),heap(left),dist,source)) child=left+1
            end if
            if (.not.heap_precedes(heap(child),last,dist,source)) exit
            heap(i)=heap(child)
            locations(heap(i))=i
            i=child
        end do
        heap(i)=last
        locations(last)=i
    end subroutine heap_pop

    ! Exact weighted shortest paths to background pixels with a 3^rank stencil.
    ! Stencil coordinate (2,...,2)+offset stores the positive cost of moving
    ! from a pixel to its neighbor at +offset; zero forbids that move.
    ! Unlike SciPy's binary custom chamfer mask, weights are actual integer costs.
    subroutine distance_transform_weighted_cdt_nd(x,weights,y,status,feature_indices)
        logical, intent(in) :: x(..) !! foreground mask, ranks 1 through 15
        integer, intent(in) :: weights(..) !! nonnegative 3-by-...-by-3 step costs
        integer, intent(out) :: y(..) !! exact shortest weighted distances to background
        integer, optional, intent(out) :: status !! 0 success, 1 bad input, 2 cost limit
        integer, optional, intent(out) :: feature_indices(:,:) !! (rank,size), one-based sites
        integer, allocatable :: w(:),steps(:),moves(:,:),heap(:),locations(:)
        integer, allocatable :: source(:),out(:)
        integer(int64), allocatable :: dist(:),strides(:),step_linear(:)
        logical, allocatable :: foreground(:)
        integer :: dims(rank(x)),coords(rank(x)),nd,n,nonzero,used,i,j,d,node,next,k,site
        integer(int64) :: total,rem,p,candidate,neighbor
        integer(int64), parameter :: infinity=huge(1_int64)/8_int64
        logical :: inside

        if (present(status)) status=1
        nd=rank(x)
        if (nd<1 .or. nd>15 .or. rank(y)/=nd .or. rank(weights)/=nd) return
        if (any(shape(x)<1) .or. any(shape(y)/=shape(x))) return
        if (any(shape(weights)/=3)) return
        total=size(x,kind=int64)
        if (present(feature_indices)) then
            if (size(feature_indices,1)/=nd .or. size(feature_indices,2)/=total) return
        end if
        if (total>int(inf_distance,int64) .or. size(weights,kind=int64)>50000_int64) then
            if (present(status)) status=2
            return
        end if
        call pack_weight_stencil(weights,w)
        if (any(w<0)) return
        i=1
        do d=1,nd
            i=i+3**(d-1)
        end do
        if (w(i)/=0) return
        if (int(maxval(w),int64)*(total-1_int64)>=int(huge(1),int64)) then
            if (present(status)) status=2
            return
        end if
        nonzero=count(w>0)
        if (total*int(nonzero,int64)>50000000_int64) then
            if (present(status)) status=2
            return
        end if
        dims=shape(x)
        allocate(strides(nd))
        strides(1)=1_int64
        do d=2,nd
            strides(d)=strides(d-1)*int(dims(d-1),int64)
        end do
        allocate(moves(nd,nonzero),steps(nonzero),step_linear(nonzero))
        k=0
        do j=1,size(w)
            if (w(j)==0) cycle
            rem=int(j-1,int64)
            k=k+1
            steps(k)=w(j)
            step_linear(k)=0_int64
            do d=1,nd
                moves(d,k)=int(modulo(rem,3_int64))-1
                rem=rem/3_int64
                step_linear(k)=step_linear(k)+int(moves(d,k),int64)*strides(d)
            end do
        end do
        call pack_mask(x,foreground)
        n=int(total)
        allocate(dist(n),source(n),heap(n),locations(n),out(n))
        dist=infinity
        source=0
        locations=0
        used=0
        do i=1,n
            if (foreground(i)) cycle
            dist(i)=0_int64
            source(i)=i
            call heap_decrease(i,heap,locations,used,dist,source)
        end do
        do while(used>0)
            call heap_pop(node,heap,locations,used,dist,source)
            rem=int(node-1,int64)
            do d=1,nd
                coords(d)=int(modulo(rem,int(dims(d),int64)))
                rem=rem/int(dims(d),int64)
            end do
            do k=1,nonzero
                inside=.true.
                do d=1,nd
                    ! Relax the predecessor pixel. The metric permits p -> q
                    ! at +moves(:,k), so p=q-moves(:,k).
                    if (coords(d)-moves(d,k)<0 .or. coords(d)-moves(d,k)>=dims(d)) then
                        inside=.false.
                        exit
                    end if
                end do
                if (.not.inside) cycle
                neighbor=int(node,int64)-step_linear(k)
                next=int(neighbor)
                candidate=dist(node)+int(steps(k),int64)
                if (candidate>int(huge(1),int64)) cycle
                if (candidate<dist(next) .or. &
                    (candidate==dist(next) .and. source(node)<source(next))) then
                    dist(next)=candidate
                    source(next)=source(node)
                    call heap_decrease(next,heap,locations,used,dist,source)
                end if
            end do
        end do
        out=-1
        do i=1,n
            if (dist(i)/=infinity) out(i)=int(dist(i))
        end do
        call scatter_integer(out,y)
        if (present(feature_indices)) then
            feature_indices=0
            do i=1,n
                site=source(i)
                if (site==0) cycle
                rem=int(site-1,int64)
                do d=1,nd
                    feature_indices(d,i)=int(modulo(rem,int(dims(d),int64)))+1
                    rem=rem/int(dims(d),int64)
                end do
            end do
        end if
        if (present(status)) status=0
    end subroutine distance_transform_weighted_cdt_nd

    subroutine distance_transform_cdt_nd(x,y,status,metric,feature_indices,axis_weights)
        logical, intent(in) :: x(..) !! true at foreground pixels; ranks 1..15
        integer, intent(out) :: y(..) !! city-block or chessboard distance to background
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 size limit
        character(len=*), optional, intent(in) :: metric !! taxicab, cityblock, or chessboard
        integer, optional, intent(out) :: feature_indices(:,:) !! one-based (rank,total size) sites
        integer, optional, intent(in) :: axis_weights(:) !! positive weighted taxicab costs per axis
        integer :: dims(rank(x)),costs(rank(x)),n,d,i,nd
        integer(int64) :: stride,outer,block,o,j,pos,rem
        logical, allocatable :: flat(:)
        integer, allocatable :: distances(:),f(:),g(:),winner(:),features(:,:),line_features(:,:)
        character(len=:), allocatable :: kind

        if (present(status)) status=1
        nd=rank(x)
        if (nd<1 .or. nd>15 .or. rank(y)/=nd) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (size(x,kind=int64)>=int(inf_distance,int64)) then
            if (present(status)) status=2
            return
        end if
        if (present(feature_indices)) then
            if (size(feature_indices,1)/=nd .or. size(feature_indices,2)/=size(x)) return
        end if
        kind='taxicab'
        if (present(metric)) kind=trim(metric)
        if (kind/='taxicab' .and. kind/='cityblock' .and. &
            kind/='manhattan' .and. kind/='chessboard') return
        costs=1
        if (present(axis_weights)) then
            if (size(axis_weights)/=nd .or. kind=='chessboard') return
            if (any(axis_weights<1)) return
            costs=axis_weights
        end if
        dims=shape(x)
        if (sum(int(dims-1,int64)*int(costs,int64))>=int(inf_distance,int64)) then
            if (present(status)) status=2
            return
        end if
        call pack_mask(x,flat)
        allocate(distances(size(flat)))
        where(flat)
            distances=inf_distance
        elsewhere
            distances=0
        end where
        if (present(feature_indices)) then
            allocate(features(nd,size(flat)))
            do pos=1_int64,size(flat,kind=int64)
                rem=pos-1_int64
                do d=1,nd
                    features(d,pos)=int(modulo(rem,int(dims(d),int64)))+1
                    rem=rem/int(dims(d),int64)
                end do
            end do
        end if
        if (any(.not.flat)) then
            do d=1,nd
                n=dims(d)
                allocate(f(n),g(n),winner(n))
                if (present(feature_indices)) allocate(line_features(nd,n))
                stride=product(int(dims(:d-1),int64))
                outer=product(int(dims(d+1:),int64))
                block=stride*int(n,int64)
                do o=0_int64,outer-1_int64
                    do j=0_int64,stride-1_int64
                        do i=0,n-1
                            pos=o*block+j+int(i,int64)*stride+1_int64
                            f(i+1)=distances(pos)
                            if (present(feature_indices)) line_features(:,i+1)=features(:,pos)
                        end do
                        if (kind=='chessboard') then
                            call chessboard_line(f,g,winner)
                        else
                            call taxicab_line(f,g,winner,costs(d))
                        end if
                        do i=0,n-1
                            pos=o*block+j+int(i,int64)*stride+1_int64
                            distances(pos)=g(i+1)
                            if (present(feature_indices)) then
                                if (winner(i+1)>0) features(:,pos)=line_features(:,winner(i+1))
                            end if
                        end do
                    end do
                end do
                deallocate(f,g,winner)
                if (present(feature_indices)) deallocate(line_features)
            end do
        end if
        where(distances==inf_distance) distances=-1
        call scatter_integer(distances,y)
        if (present(feature_indices)) feature_indices=features
        if (present(status)) status=0
    end subroutine distance_transform_cdt_nd
end module scifort_ndimage_chamfer
