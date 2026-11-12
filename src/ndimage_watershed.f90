! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Integer image foresting transform: minimax adjacent-gradient paths,
! bucket priorities, and priority for positive markers on ties.
module scifort_ndimage_watershed
    use, intrinsic :: iso_fortran_env, only : int64
    implicit none
    private
    public :: watershed_ift_nd
    integer, parameter :: max_footprint=50000
contains
    subroutine pack_integer(x,flat)
        integer, intent(in) :: x(..) !! arbitrary-rank source image
        integer, allocatable, intent(out) :: flat(:) !! flattened source
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
    end subroutine pack_integer

    subroutine pack_logical(x,flat)
        logical, intent(in) :: x(..) !! arbitrary-rank stencil
        logical, allocatable, intent(out) :: flat(:) !! flattened stencil
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
    end subroutine pack_logical

    subroutine scatter_integer(flat,y)
        integer, intent(in) :: flat(:) !! flattened labels
        integer, intent(out) :: y(..) !! output with input's shape
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

    subroutine watershed_stencil(nd,structure,offsets,status)
        integer, intent(in) :: nd !! source image rank
        logical, optional, intent(in) :: structure(..) !! 3-by-...-by-3 stencil
        integer, allocatable, intent(out) :: offsets(:,:) !! active neighbors in C order
        integer, intent(out) :: status !! status code
        logical, allocatable :: flat(:)
        integer :: d,n,k,i,coord(nd),stride(nd),j,loc,nactive
        integer(int64) :: total,rem

        status=1
        if (.not.present(structure)) then
            allocate(offsets(nd,2*nd))
            offsets=0
            do d=1,nd
                offsets(d,d)=-1
                offsets(d,2*nd-d+1)=1
            end do
            status=0
            return
        end if
        if (rank(structure)/=nd) return
        if (any(shape(structure)/=3)) return
        total=1_int64
        do d=1,nd
            if (total>max_footprint/3_int64) then
                status=2
                return
            end if
            total=total*3_int64
        end do
        call pack_logical(structure,flat)
        nactive=count(flat)
        n=0
        allocate(offsets(nd,max(0,nactive)))
        stride(1)=1
        do d=2,nd
            stride(d)=3*stride(d-1)
        end do
        do i=0,int(total)-1
            rem=int(i,int64)
            loc=1
            do d=nd,1,-1
                coord(d)=int(modulo(rem,3_int64))
                rem=rem/3_int64
                loc=loc+coord(d)*stride(d)
            end do
            if (.not.flat(loc)) cycle
            if (all(coord==1)) cycle
            n=n+1
            offsets(:,n)=coord-1
        end do
        if (n/=nactive) offsets=offsets(:,:n)
        status=0
    end subroutine watershed_stencil

    subroutine watershed_ift_nd(x,markers,y,status,structure,legacy_linear_bounds)
        integer, intent(in) :: x(..) !! nonnegative intensity image, values 0..65535
        integer, intent(in) :: markers(..) !! nonzero seeds; negative labels are background
        integer, intent(out) :: y(..) !! propagated integer labels
        integer, optional, intent(out) :: status !! 0 success, 1 invalid, 2 too large
        logical, optional, intent(in) :: structure(..) !! neighborhood shape (3,...,3)
        logical, optional, intent(in) :: legacy_linear_bounds !! reproduce SciPy 1.17 edge wrap
        integer, allocatable :: values(:),labels(:),offsets(:,:)
        integer, allocatable :: cost(:),head(:),tail(:),next(:),prev(:)
        logical, allocatable :: queued(:),done(:)
        integer :: dims(rank(x)),loc(rank(x)),nd,n,d,k,nmax,idx,oldcost,newcost
        integer :: p,q,priority,code,marker
        integer(int64) :: stride(rank(x)),stride_c(rank(x)),rem,position,index_c,position_c
        logical :: legacy
        logical :: inside

        if (present(status)) status=1
        nd=rank(x)
        if (nd<1 .or. nd>15 .or. rank(markers)/=nd .or. rank(y)/=nd) return
        if (any(shape(x)/=shape(markers)) .or. any(shape(x)/=shape(y))) return
        if (any(shape(x)<1) .or. size(x,kind=int64)>int(huge(1),int64)) return
        call watershed_stencil(nd,structure,offsets,code)
        if (code/=0) then
            if (present(status)) status=code
            return
        end if
        call pack_integer(x,values)
        if (any(values<0) .or. any(values>65535)) return
        call pack_integer(markers,labels)
        n=size(values)
        nmax=maxval(values)
        dims=shape(x)
        stride(1)=1_int64
        do d=2,nd
            stride(d)=stride(d-1)*int(dims(d-1),int64)
        end do
        stride_c(nd)=1_int64
        do d=nd-1,1,-1
            stride_c(d)=stride_c(d+1)*int(dims(d+1),int64)
        end do
        legacy=.false.
        if (present(legacy_linear_bounds)) legacy=legacy_linear_bounds
        allocate(cost(n),next(n),prev(n),queued(n),done(n))
        allocate(head(0:nmax),tail(0:nmax))
        head=0
        tail=0
        next=0
        prev=0
        queued=.false.
        done=.false.
        cost=nmax+1
        ! Marker initialization follows C-order traversal to agree with
        ! SciPy's tie handling, even though the image is stored Fortran-order.
        do index_c=0_int64,int(n-1,int64)
            rem=index_c
            position=1_int64
            do d=nd,1,-1
                position=position+modulo(rem,int(dims(d),int64))*stride(d)
                rem=rem/int(dims(d),int64)
            end do
            p=int(position)
            marker=labels(p)
            if (marker==0) cycle
            cost(p)=0
            if (marker>0) then
                call enqueue_front(p,0)
            else
                call enqueue_back(p,0)
            end if
        end do
        do priority=0,nmax
            do while(head(priority)/=0)
                p=head(priority)
                call remove_entry(p,priority)
                done(p)=.true.
                rem=int(p-1,int64)
                do d=1,nd
                    loc(d)=int(modulo(rem,int(dims(d),int64)))
                    rem=rem/int(dims(d),int64)
                end do
                do k=1,size(offsets,2)
                    if (legacy) then
                        ! Compatibility with SciPy's current C implementation:
                        ! it tests the linearized neighbor index rather than
                        ! every original axis, so some edge offsets wrap.
                        position_c=1_int64
                        do d=1,nd
                            position_c=position_c+int(loc(d)+offsets(d,k),int64)*stride_c(d)
                        end do
                        if (position_c<1_int64 .or. position_c>int(n,int64)) cycle
                        rem=position_c-1_int64
                        position=1_int64
                        do d=1,nd
                            position=position+int(rem/stride_c(d),int64)*stride(d)
                            rem=modulo(rem,stride_c(d))
                        end do
                    else
                        inside=.true.
                        position=int(p,int64)
                        do d=1,nd
                            if (loc(d)+offsets(d,k)<0 .or. loc(d)+offsets(d,k)>=dims(d)) then
                                inside=.false.
                                exit
                            end if
                            position=position+int(offsets(d,k),int64)*stride(d)
                        end do
                        if (.not.inside) cycle
                    end if
                    q=int(position)
                    if (done(q)) cycle
                    newcost=max(cost(p),abs(values(p)-values(q)))
                    if (newcost>=cost(q)) cycle
                    oldcost=cost(q)
                    if (queued(q)) call remove_entry(q,oldcost)
                    cost(q)=newcost
                    labels(q)=labels(p)
                    if (labels(q)>0) then
                        call enqueue_front(q,newcost)
                    else
                        call enqueue_back(q,newcost)
                    end if
                end do
            end do
        end do
        call scatter_integer(labels,y)
        if (present(status)) status=0
    contains
        subroutine remove_entry(node,bucket)
            integer, intent(in) :: node !! node index
            integer, intent(in) :: bucket !! bucket priority
            integer :: a,b
            a=prev(node)
            b=next(node)
            if (a==0) then
                head(bucket)=b
            else
                next(a)=b
            end if
            if (b==0) then
                tail(bucket)=a
            else
                prev(b)=a
            end if
            next(node)=0
            prev(node)=0
            queued(node)=.false.
        end subroutine remove_entry

        subroutine enqueue_front(node,bucket)
            integer, intent(in) :: node !! node index
            integer, intent(in) :: bucket !! bucket priority
            prev(node)=0
            next(node)=head(bucket)
            if (head(bucket)/=0) prev(head(bucket))=node
            head(bucket)=node
            if (tail(bucket)==0) tail(bucket)=node
            queued(node)=.true.
        end subroutine enqueue_front

        subroutine enqueue_back(node,bucket)
            integer, intent(in) :: node !! node index
            integer, intent(in) :: bucket !! bucket priority
            next(node)=0
            prev(node)=tail(bucket)
            if (tail(bucket)/=0) next(tail(bucket))=node
            tail(bucket)=node
            if (head(bucket)==0) head(bucket)=node
            queued(node)=.true.
        end subroutine enqueue_back
    end subroutine watershed_ift_nd
end module scifort_ndimage_watershed
