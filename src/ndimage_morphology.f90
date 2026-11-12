! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Binary morphology and connected-component labeling in 1 to 15 dimensions.
! Neighborhoods are radius-one Boolean stencils with L1 connectivity.
module scifort_ndimage_morphology
    use, intrinsic :: iso_fortran_env, only : int64
    implicit none
    private
    integer, parameter :: dp=kind(1.0d0)
    integer, parameter, public :: morphology_status_success=0
    integer, parameter, public :: morphology_status_invalid_input=1
    integer, parameter, public :: morphology_status_too_expensive=2
    integer, parameter :: max_neighbors=50000
    public :: binary_dilation_nd, binary_erosion_nd, label_nd
    public :: binary_opening_nd, binary_closing_nd, binary_fill_holes_nd
    public :: binary_propagation_nd, binary_hit_or_miss_nd
    public :: generate_binary_structure_nd, iterate_structure_nd
    interface label_nd
        module procedure label_logical_nd, label_integer_nd, label_real_nd
    end interface
contains
    subroutine generate_neighbors(nd,connectivity,offsets,status)
        integer, intent(in) :: nd !! rank, 1 through 15
        integer, intent(in) :: connectivity !! maximum nonzero offset components
        integer, allocatable, intent(out) :: offsets(:,:) !! (rank, neighbors)
        integer, intent(out) :: status !! success or too large
        integer, allocatable :: tmp(:,:),vec(:)
        integer :: count
        status=morphology_status_invalid_input
        if (nd<1 .or. nd>15 .or. connectivity<1 .or. connectivity>nd) return
        allocate(tmp(nd,max_neighbors),vec(nd))
        count=0
        call enumerate(1,0)
        if (count>max_neighbors) then
            status=morphology_status_too_expensive
            return
        end if
        offsets=tmp(:,:count)
        status=morphology_status_success
    contains
        recursive subroutine enumerate(axis,used)
            integer, intent(in) :: axis !! current dimension
            integer, intent(in) :: used !! existing nonzero coordinate count
            integer :: delta
            if (count>max_neighbors) return
            if (axis>nd) then
                count=count+1
                if (count<=max_neighbors) tmp(:,count)=vec
                return
            end if
            do delta=-1,1
                if (used+abs(delta)>connectivity) cycle
                vec(axis)=delta
                call enumerate(axis+1,used+abs(delta))
                if (count>max_neighbors) exit
            end do
        end subroutine enumerate
    end subroutine generate_neighbors

    subroutine morphological_kernel(samples,dims,offsets,iterations,border,dilate,result,mask,status)
        logical, intent(in) :: samples(:) !! column-major foreground mask
        integer, intent(in) :: dims(:) !! shape
        integer, intent(in) :: offsets(:,:) !! neighboring offsets
        integer, intent(in) :: iterations !! positive passes; <=0 means until stable
        logical, intent(in) :: border !! exterior value
        logical, intent(in) :: dilate !! dilation rather than erosion
        logical, intent(out) :: result(:) !! processed foreground mask
        logical, optional, intent(in) :: mask(:) !! writable locations on each iteration
        integer, optional, intent(out) :: status !! success or iteration limit
        logical, allocatable :: source(:)
        integer(int64) :: stride(size(dims)),pos,q,remainder
        integer :: loc(size(dims)),d,k,limit
        logical :: inside,hit,stable

        stride(1)=1_int64
        do d=2,size(dims)
            stride(d)=stride(d-1)*int(dims(d-1),int64)
        end do
        source=samples
        if (present(status)) status=morphology_status_success
        limit=iterations
        if (iterations<=0) limit=size(samples)+1
        do k=1,limit
            do pos=1_int64,size(samples,kind=int64)
                if (present(mask)) then
                    if (.not.mask(pos)) then
                        result(pos)=source(pos)
                        cycle
                    end if
                end if
                remainder=pos-1_int64
                do d=1,size(dims)
                    loc(d)=int(modulo(remainder,int(dims(d),int64)))
                    remainder=remainder/int(dims(d),int64)
                end do
                result(pos)=.not.dilate
                block
                    integer :: j
                    do j=1,size(offsets,2)
                        inside=.true.
                        q=pos
                        do d=1,size(dims)
                            if (loc(d)+offsets(d,j)<0 .or. &
                                loc(d)+offsets(d,j)>=dims(d)) then
                                inside=.false.
                                exit
                            end if
                            q=q+int(offsets(d,j),int64)*stride(d)
                        end do
                        hit=border
                        if (inside) hit=source(q)
                        if (dilate .and. hit) then
                            result(pos)=.true.
                            exit
                        else if (.not.dilate .and. .not.hit) then
                            result(pos)=.false.
                            exit
                        end if
                    end do
                end block
            end do
            stable=all(result.eqv.source)
            if (stable) exit
            if (iterations<=0 .and. k==limit) then
                if (present(status)) status=morphology_status_too_expensive
                exit
            end if
            source=result
        end do
    end subroutine morphological_kernel

    integer function root(parent,i) result(r)
        integer, intent(inout) :: parent(:) !! union-find parent array
        integer, intent(in) :: i !! foreground node index
        integer :: j,p
        r=i
        do while(parent(r)/=r)
            r=parent(r)
        end do
        j=i
        do while(parent(j)/=r)
            p=parent(j)
            parent(j)=r
            j=p
        end do
    end function root

    subroutine label_kernel(samples,dims,offsets,labels,num_features)
        logical, intent(in) :: samples(:) !! foreground mask
        integer, intent(in) :: dims(:) !! extents
        integer, intent(in) :: offsets(:,:) !! connectivity neighborhood
        integer, intent(out) :: labels(:) !! contiguous labels, background zero
        integer, intent(out) :: num_features !! number of connected components
        integer, allocatable :: parent(:),map(:)
        integer(int64) :: stride(size(dims)),pos,q,remainder
        integer :: loc(size(dims)),d,k,a,b,j
        logical :: inside

        stride(1)=1_int64
        do d=2,size(dims)
            stride(d)=stride(d-1)*int(dims(d-1),int64)
        end do
        allocate(parent(size(samples)),map(size(samples)))
        parent=0
        map=0
        labels=0
        num_features=0
        do pos=1_int64,size(samples,kind=int64)
            if (.not.samples(pos)) cycle
            parent(pos)=int(pos)
            remainder=pos-1_int64
            do d=1,size(dims)
                loc(d)=int(modulo(remainder,int(dims(d),int64)))
                remainder=remainder/int(dims(d),int64)
            end do
            do k=1,size(offsets,2)
                inside=.true.
                q=pos
                do d=1,size(dims)
                    if (loc(d)+offsets(d,k)<0 .or. &
                        loc(d)+offsets(d,k)>=dims(d)) then
                        inside=.false.
                        exit
                    end if
                    q=q+int(offsets(d,k),int64)*stride(d)
                end do
                if (.not.inside .or. q>=pos) cycle
                if (.not.samples(q)) cycle
                a=root(parent,int(pos))
                b=root(parent,int(q))
                if (a/=b) parent(max(a,b))=min(a,b)
            end do
        end do
        do pos=1_int64,size(samples,kind=int64)
            if (.not.samples(pos)) cycle
            j=root(parent,int(pos))
            if (map(j)==0) then
                num_features=num_features+1
                map(j)=num_features
            end if
            labels(pos)=map(j)
        end do
    end subroutine label_kernel

    subroutine pack_logical(x,buf)
        logical, intent(in) :: x(..) !! arbitrary-rank source mask
        logical, allocatable, intent(out) :: buf(:) !! packed mask
        select rank(x)
        rank(1)
            buf=reshape(x,[size(x)])
        rank(2)
            buf=reshape(x,[size(x)])
        rank(3)
            buf=reshape(x,[size(x)])
        rank(4)
            buf=reshape(x,[size(x)])
        rank(5)
            buf=reshape(x,[size(x)])
        rank(6)
            buf=reshape(x,[size(x)])
        rank(7)
            buf=reshape(x,[size(x)])
        rank(8)
            buf=reshape(x,[size(x)])
        rank(9)
            buf=reshape(x,[size(x)])
        rank(10)
            buf=reshape(x,[size(x)])
        rank(11)
            buf=reshape(x,[size(x)])
        rank(12)
            buf=reshape(x,[size(x)])
        rank(13)
            buf=reshape(x,[size(x)])
        rank(14)
            buf=reshape(x,[size(x)])
        rank(15)
            buf=reshape(x,[size(x)])
        rank default
            allocate(buf(0))
        end select
    end subroutine pack_logical


    subroutine prepare_offsets(nd,connectivity,structure,offsets,status,complement,origin)
        integer, intent(in) :: nd !! dimensionality of the image
        integer, intent(in) :: connectivity !! connectivity if no custom structure
        logical, optional, intent(in) :: structure(..) !! custom structuring element
        integer, allocatable, intent(out) :: offsets(:,:) !! active offsets
        integer, intent(out) :: status !! success, invalid input, or resource limit
        logical, optional, intent(in) :: complement !! invert selected footprint entries
        logical, allocatable :: active(:)
        integer, optional, intent(in) :: origin(:) !! per-axis origin displacements
        integer :: dims(nd),loc(nd),origins(nd),j,d,k,n
        integer(int64) :: rem

        origins=0
        if (present(origin)) then
            if (size(origin)/=nd) then
                status=morphology_status_invalid_input
                return
            end if
            origins=origin
        end if
        if (.not.present(structure)) then
            if (any(origins<(-1)) .or. any(origins>1)) then
                status=morphology_status_invalid_input
                return
            end if
            call generate_neighbors(nd,connectivity,offsets,status)
            if (status==morphology_status_success) then
                do j=1,size(offsets,2)
                    offsets(:,j)=offsets(:,j)-origins
                end do
            end if
            return
        end if
        status=morphology_status_invalid_input
        if (rank(structure)/=nd) return
        if (any(shape(structure)<1)) return
        if (size(structure,kind=int64)>int(max_neighbors,int64)) then
            status=morphology_status_too_expensive
            return
        end if
        dims=shape(structure)
        if (any(origins<-(dims/2)) .or. any(origins>(dims-1)/2)) return
        call pack_logical(structure,active)
        if (present(complement)) then
            if (complement) active=.not.active
        end if
        n=count(active)
        allocate(offsets(nd,n))
        k=0
        do j=1,size(active)
            if (.not.active(j)) cycle
            rem=int(j-1,int64)
            do d=1,nd
                loc(d)=int(modulo(rem,int(dims(d),int64)))-dims(d)/2-origins(d)
                rem=rem/int(dims(d),int64)
            end do
            k=k+1
            offsets(:,k)=loc
        end do
        status=morphology_status_success
    end subroutine prepare_offsets

    ! A compact structure spans selected axes in ascending image-axis order,
    ! while origin(:) follows the caller's axes(:) order, like scipy.ndimage.
    ! Embed the compact offsets in full image coordinates without allocating
    ! a potentially enormous rank-15 Boolean footprint.
    subroutine prepare_morphology_axes_offsets(nd,connectivity,structure,origin,axes,offsets,status)
        integer, intent(in) :: nd !! rank of the full image
        integer, intent(in) :: connectivity !! connectivity among selected axes
        logical, optional, intent(in) :: structure(..) !! compact structuring element
        integer, optional, intent(in) :: origin(:) !! shifts in axes(:) order
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        integer, allocatable, intent(out) :: offsets(:,:) !! full-rank neighbor offsets
        integer, intent(out) :: status !! validation result
        integer, allocatable :: small(:,:),selected(:),ordered_origin(:)
        integer :: nsel,d,i,j

        if (.not.present(axes)) then
            call prepare_offsets(nd,connectivity,structure,offsets,status,origin=origin)
            return
        end if
        status=morphology_status_invalid_input
        nsel=size(axes)
        if (nsel>nd .or. any(axes<1) .or. any(axes>nd)) return
        do i=1,nsel
            if (count(axes==axes(i))/=1) return
        end do
        if (present(origin)) then
            if (size(origin)/=nsel) return
        end if
        if (nsel==0) then
            if (present(structure)) return
            allocate(offsets(nd,1))
            offsets=0
            status=morphology_status_success
            return
        end if
        allocate(selected(nsel),ordered_origin(nsel))
        ordered_origin=0
        j=0
        do d=1,nd
            do i=1,nsel
                if (axes(i)/=d) cycle
                j=j+1
                selected(j)=d
                if (present(origin)) ordered_origin(j)=origin(i)
                exit
            end do
        end do
        call prepare_offsets(nsel,connectivity,structure,small,status,origin=ordered_origin)
        if (status/=morphology_status_success) return
        allocate(offsets(nd,size(small,2)))
        offsets=0
        do j=1,nsel
            offsets(selected(j),:)=small(j,:)
        end do
    end subroutine prepare_morphology_axes_offsets

    subroutine scatter_logical(buf,y)
        logical, intent(in) :: buf(:) !! column-major source
        logical, intent(out) :: y(..) !! destination of same size
        select rank(y)
        rank(1)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(2)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(3)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(4)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(5)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(6)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(7)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(8)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(9)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(10)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(11)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(12)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(13)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(14)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(15)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank default
        end select
    end subroutine scatter_logical

    subroutine scatter_integer(buf,y)
        integer, intent(in) :: buf(:) !! column-major source
        integer, intent(out) :: y(..) !! destination of same size
        select rank(y)
        rank(1)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(2)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(3)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(4)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(5)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(6)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(7)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(8)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(9)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(10)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(11)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(12)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(13)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(14)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(15)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank default
        end select
    end subroutine scatter_integer

    subroutine morphological_work(x,y,dilate,status,connectivity,iterations,border_value,structure,origin,mask,axes)
        logical, intent(in) :: x(..) !! Boolean input, ranks 1 through 15
        logical, intent(out) :: y(..) !! same-shape morphological result
        logical, intent(in) :: dilate !! select dilation or erosion
        integer, optional, intent(out) :: status !! zero for success
        integer, optional, intent(in) :: connectivity !! 1 (faces) to rank
        logical, optional, intent(in) :: structure(..) !! custom footprint, same rank as image
        integer, optional, intent(in) :: iterations !! count of passes, default one
        logical, optional, intent(in) :: border_value !! outside mask, default false
        logical, allocatable :: flat(:),result(:),editable(:)
        integer, allocatable :: offsets(:,:)
        integer, optional, intent(in) :: origin(:) !! per-axis displacement
        logical, optional, intent(in) :: mask(..) !! only locations allowed to change
        integer, optional, intent(in) :: axes(:) !! selected one-based image axes
        integer :: neighbors,passes,code
        logical :: border

        if (present(status)) status=morphology_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        if (present(mask)) then
            if (rank(mask)/=rank(x)) return
            if (any(shape(mask)/=shape(x))) return
        end if
        neighbors=1
        if (present(connectivity)) neighbors=connectivity
        passes=1
        if (present(iterations)) passes=iterations
        border=.false.
        if (present(border_value)) border=border_value
        call prepare_morphology_axes_offsets(rank(x),neighbors,structure,origin,axes,offsets,code)
        if (code/=morphology_status_success) then
            if (present(status)) status=code
            return
        end if
        if (present(mask)) then
            call pack_logical(mask,editable)
        else
            editable=spread(.true.,1,size(x))
        end if
        call pack_logical(x,flat)
        allocate(result(size(flat)))
        ! Dilation uses the reflected offsets; erosion uses them directly.
        if (dilate) offsets=-offsets
        call morphological_kernel(flat,shape(x),offsets,passes,border,dilate,result,editable,code)
        if (code/=morphology_status_success) then
            if (present(status)) status=code
            return
        end if
        call scatter_logical(result,y)
        if (present(status)) status=morphology_status_success
    end subroutine morphological_work

    subroutine binary_dilation_nd(x,y,status,connectivity,iterations,border_value,structure,origin,mask,axes)
        logical, intent(in) :: x(..) !! foreground pixels
        logical, intent(out) :: y(..) !! dilated foreground
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: connectivity !! maximum Manhattan offset
        logical, optional, intent(in) :: structure(..) !! custom footprint, same rank as image
        integer, optional, intent(in) :: iterations !! number of dilations
        logical, optional, intent(in) :: border_value !! outside value
        integer, optional, intent(in) :: origin(:) !! per-axis displacement
        logical, optional, intent(in) :: mask(..) !! only locations allowed to change
        integer, optional, intent(in) :: axes(:) !! selected one-based image axes
        call morphological_work(x,y,.true.,status,connectivity,iterations,border_value,structure,origin,mask,axes)
    end subroutine binary_dilation_nd

    subroutine binary_erosion_nd(x,y,status,connectivity,iterations,border_value,structure,origin,mask,axes)
        logical, intent(in) :: x(..) !! foreground pixels
        logical, intent(out) :: y(..) !! eroded foreground
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: connectivity !! maximum Manhattan offset
        logical, optional, intent(in) :: structure(..) !! custom footprint, same rank as image
        integer, optional, intent(in) :: iterations !! number of erosions
        logical, optional, intent(in) :: border_value !! outside value
        integer, optional, intent(in) :: origin(:) !! per-axis displacement
        logical, optional, intent(in) :: mask(..) !! only locations allowed to change
        integer, optional, intent(in) :: axes(:) !! selected one-based image axes
        call morphological_work(x,y,.false.,status,connectivity,iterations,border_value,structure,origin,mask,axes)
    end subroutine binary_erosion_nd

    subroutine label_logical_nd(x,y,num_features,status,connectivity,structure,axes)
        logical, intent(in) :: x(..) !! Boolean source mask
        integer, intent(out) :: y(..) !! integer label array, zero for background
        integer, intent(out) :: num_features !! distinct foreground components
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: connectivity !! 1=faces, rank=full neighbors
        logical, optional, intent(in) :: structure(..) !! full or selected-axis footprint
        integer, optional, intent(in) :: axes(:) !! selected one-based axes, all by default
        logical, allocatable :: flat(:)
        integer, allocatable :: labels(:),offsets(:,:)
        integer :: code,neighbors

        num_features=0
        if (present(status)) status=morphology_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        neighbors=1
        if (present(connectivity)) neighbors=connectivity
        call prepare_morphology_axes_offsets(rank(x),neighbors,structure, &
            axes=axes,offsets=offsets,status=code)
        if (code/=morphology_status_success) then
            if (present(status)) status=code
            return
        end if
        if (present(structure)) then
            block
                integer :: j,k
                logical :: paired
                do j=1,size(offsets,2)
                    paired=.false.
                    do k=1,size(offsets,2)
                        if (all(offsets(:,j)==-offsets(:,k))) paired=.true.
                    end do
                    if (.not.paired) return
                end do
            end block
        end if
        call pack_logical(x,flat)
        allocate(labels(size(flat)))
        call label_kernel(flat,shape(x),offsets,labels,num_features)
        call scatter_integer(labels,y)
        if (present(status)) status=morphology_status_success
    end subroutine label_logical_nd
    subroutine morphological_sequence(x,y,closing,status,connectivity,iterations,border_value,structure,origin,mask,axes)
        logical, intent(in) :: x(..) !! input mask, ranks 1 through 15
        logical, intent(out) :: y(..) !! result of opening or closing
        logical, intent(in) :: closing !! close when true, open otherwise
        integer, optional, intent(out) :: status !! status code
        integer, optional, intent(in) :: connectivity !! stencil connectivity
        logical, optional, intent(in) :: structure(..) !! custom footprint, same rank as image
        integer, optional, intent(in) :: iterations !! positive passes in each stage
        logical, optional, intent(in) :: border_value !! out-of-bounds state
        logical, allocatable :: flat(:),tmp(:),result(:),editable(:)
        integer, allocatable :: offsets(:,:)
        integer, optional, intent(in) :: origin(:) !! per-axis displacement
        logical, optional, intent(in) :: mask(..) !! only locations allowed to change
        integer, optional, intent(in) :: axes(:) !! selected one-based image axes
        integer :: neighbors,passes,code
        logical :: border

        if (present(status)) status=morphology_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        if (present(mask)) then
            if (rank(mask)/=rank(x)) return
            if (any(shape(mask)/=shape(x))) return
        end if
        neighbors=1
        if (present(connectivity)) neighbors=connectivity
        passes=1
        if (present(iterations)) passes=iterations
        border=.false.
        if (present(border_value)) border=border_value
        call prepare_morphology_axes_offsets(rank(x),neighbors,structure,origin,axes,offsets,code)
        if (code/=morphology_status_success) then
            if (present(status)) status=code
            return
        end if
        if (present(mask)) then
            call pack_logical(mask,editable)
        else
            editable=spread(.true.,1,size(x))
        end if
        call pack_logical(x,flat)
        allocate(tmp(size(flat)),result(size(flat)))
        if (closing) then
            offsets=-offsets
            call morphological_kernel(flat,shape(x),offsets,passes,border,.true.,tmp,editable,code)
            offsets=-offsets
            if (code/=morphology_status_success) then
                if (present(status)) status=code
                return
            end if
            call morphological_kernel(tmp,shape(x),offsets,passes,border,.false.,result,editable,code)
        else
            call morphological_kernel(flat,shape(x),offsets,passes,border,.false.,tmp,editable,code)
            offsets=-offsets
            if (code/=morphology_status_success) then
                if (present(status)) status=code
                return
            end if
            call morphological_kernel(tmp,shape(x),offsets,passes,border,.true.,result,editable,code)
        end if
        if (code/=morphology_status_success) then
            if (present(status)) status=code
            return
        end if
        call scatter_logical(result,y)
        if (present(status)) status=morphology_status_success
    end subroutine morphological_sequence

    subroutine binary_opening_nd(x,y,status,connectivity,iterations,border_value,structure,origin,mask,axes)
        logical, intent(in) :: x(..) !! input binary image
        logical, intent(out) :: y(..) !! erosion followed by dilation
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: connectivity !! neighborhood radius-one connectivity
        logical, optional, intent(in) :: structure(..) !! custom footprint, same rank as image
        integer, optional, intent(in) :: iterations !! erosion and dilation passes
        logical, optional, intent(in) :: border_value !! constant extension
        integer, optional, intent(in) :: origin(:) !! per-axis displacement
        logical, optional, intent(in) :: mask(..) !! only locations allowed to change
        integer, optional, intent(in) :: axes(:) !! selected one-based image axes
        call morphological_sequence(x,y,.false.,status,connectivity,iterations,border_value,structure,origin,mask,axes)
    end subroutine binary_opening_nd

    subroutine binary_closing_nd(x,y,status,connectivity,iterations,border_value,structure,origin,mask,axes)
        logical, intent(in) :: x(..) !! input binary image
        logical, intent(out) :: y(..) !! dilation followed by erosion
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: connectivity !! neighborhood radius-one connectivity
        logical, optional, intent(in) :: structure(..) !! custom footprint, same rank as image
        integer, optional, intent(in) :: iterations !! dilation and erosion passes
        logical, optional, intent(in) :: border_value !! constant extension
        integer, optional, intent(in) :: origin(:) !! per-axis displacement
        logical, optional, intent(in) :: mask(..) !! only locations allowed to change
        integer, optional, intent(in) :: axes(:) !! selected one-based image axes
        call morphological_sequence(x,y,.true.,status,connectivity,iterations,border_value,structure,origin,mask,axes)
    end subroutine binary_closing_nd

    subroutine flood_fill(seed,mask,dims,offsets,reached)
        logical, intent(in) :: seed(:) !! initial foreground seeds
        logical, intent(in) :: mask(:) !! allowed additional points
        integer, intent(in) :: dims(:) !! grid extents
        integer, intent(in) :: offsets(:,:) !! neighbor offsets
        logical, intent(out) :: reached(:) !! full geodesic propagation
        integer, allocatable :: queue(:)
        integer(int64) :: stride(size(dims)),rem,q,pos
        integer :: head,tail,n,d,k,loc(size(dims))
        logical :: inside

        n=size(seed)
        allocate(queue(n))
        stride(1)=1_int64
        do d=2,size(dims)
            stride(d)=stride(d-1)*int(dims(d-1),int64)
        end do
        reached=seed
        head=1
        tail=0
        do d=1,n
            if (.not.seed(d)) cycle
            tail=tail+1
            queue(tail)=d
        end do
        do while(head<=tail)
            pos=int(queue(head),int64)
            head=head+1
            rem=pos-1_int64
            do d=1,size(dims)
                loc(d)=int(modulo(rem,int(dims(d),int64)))
                rem=rem/int(dims(d),int64)
            end do
            do k=1,size(offsets,2)
                q=pos
                inside=.true.
                do d=1,size(dims)
                    if (loc(d)+offsets(d,k)<0 .or. &
                        loc(d)+offsets(d,k)>=dims(d)) then
                        inside=.false.
                        exit
                    end if
                    q=q+int(offsets(d,k),int64)*stride(d)
                end do
                if (.not.inside) cycle
                if (reached(q) .or. .not.mask(q)) cycle
                reached(q)=.true.
                tail=tail+1
                queue(tail)=int(q)
            end do
        end do
    end subroutine flood_fill

    subroutine binary_propagation_nd(x,mask,y,status,connectivity,structure,axes)
        logical, intent(in) :: x(..) !! foreground seeds
        logical, intent(in) :: mask(..) !! allowed propagation positions
        logical, intent(out) :: y(..) !! propagated foreground
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: connectivity !! neighbor connectivity
        logical, optional, intent(in) :: structure(..) !! full-rank or selected-axis footprint
        integer, optional, intent(in) :: axes(:) !! selected one-based image axes
        logical, allocatable :: seeds(:),allowed(:),reached(:)
        integer, allocatable :: offsets(:,:)
        integer :: neighbors,code

        if (present(status)) status=morphology_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(mask)/=rank(x) .or. rank(y)/=rank(x)) return
        if (any(shape(x)/=shape(mask)) .or. any(shape(x)/=shape(y))) return
        if (any(shape(x)<1) .or. size(x,kind=int64)>int(huge(1),int64)) return
        neighbors=1
        if (present(connectivity)) neighbors=connectivity
        call prepare_morphology_axes_offsets(rank(x),neighbors,structure, &
            axes=axes,offsets=offsets,status=code)
        if (code/=morphology_status_success) then
            if (present(status)) status=code
            return
        end if
        call pack_logical(x,seeds)
        call pack_logical(mask,allowed)
        allocate(reached(size(seeds)))
        call flood_fill(seeds,allowed,shape(x),offsets,reached)
        call scatter_logical(reached,y)
        if (present(status)) status=morphology_status_success
    end subroutine binary_propagation_nd

    subroutine binary_fill_holes_nd(x,y,status,connectivity,structure,axes)
        logical, intent(in) :: x(..) !! binary image with possible internal holes
        logical, intent(out) :: y(..) !! internal holes filled
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: connectivity !! background connectivity
        logical, optional, intent(in) :: structure(..) !! full-rank or selected-axis footprint
        integer, optional, intent(in) :: axes(:) !! spatial dimensions defining cavities
        logical, allocatable :: flat(:),seed(:),background(:),reached(:)
        integer, allocatable :: offsets(:,:)
        integer(int64) :: pos,rem
        integer :: neighbors,code,d,loc,dims(rank(x))
        logical :: edge,selected(rank(x))

        if (present(status)) status=morphology_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        neighbors=1
        if (present(connectivity)) neighbors=connectivity
        call prepare_morphology_axes_offsets(rank(x),neighbors,structure, &
            axes=axes,offsets=offsets,status=code)
        if (code/=morphology_status_success) then
            if (present(status)) status=code
            return
        end if
        selected=.true.
        if (present(axes)) then
            selected=.false.
            if (size(axes)==0) then
                call pack_logical(x,flat)
                call scatter_logical(flat,y)
                if (present(status)) status=morphology_status_success
                return
            end if
            selected(axes)=.true.
        end if
        dims=shape(x)
        call pack_logical(x,flat)
        background=.not.flat
        allocate(seed(size(flat)),reached(size(flat)))
        seed=.false.
        do pos=1_int64,size(flat,kind=int64)
            if (.not.background(pos)) cycle
            rem=pos-1_int64
            edge=.false.
            do d=1,rank(x)
                loc=int(modulo(rem,int(dims(d),int64)))
                if (selected(d)) then
                    if (loc==0 .or. loc==dims(d)-1) edge=.true.
                end if
                rem=rem/int(dims(d),int64)
            end do
            seed(pos)=edge
        end do
        call flood_fill(seed,background,shape(x),offsets,reached)
        call scatter_logical(.not.reached,y)
        if (present(status)) status=morphology_status_success
    end subroutine binary_fill_holes_nd

    subroutine binary_hit_or_miss_nd(x,y,status,structure1,structure2,origin1,origin2)
        logical, intent(in) :: x(..) !! binary image, ranks 1 through 15
        logical, intent(out) :: y(..) !! matches of required hits and misses
        integer, optional, intent(out) :: status !! zero for success
        logical, optional, intent(in) :: structure1(..) !! foreground pattern
        logical, optional, intent(in) :: structure2(..) !! background pattern
        integer, optional, intent(in) :: origin1(:) !! foreground origin
        integer, optional, intent(in) :: origin2(:) !! background origin, defaults to origin1
        integer, allocatable :: hit_offsets(:,:),miss_offsets(:,:),all_offsets(:,:),temporary(:,:)
        logical, allocatable :: flat(:),hits(:),misses(:)
        integer :: code,nd,j,k

        if (present(status)) status=morphology_status_invalid_input
        nd=rank(x)
        if (nd<1 .or. nd>15 .or. rank(y)/=nd) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        call prepare_offsets(nd,1,structure1,hit_offsets,code,origin=origin1)
        if (code/=morphology_status_success) then
            if (present(status)) status=code
            return
        end if
        if (present(structure2)) then
            if (present(origin2)) then
                call prepare_offsets(nd,1,structure2,miss_offsets,code,origin=origin2)
            else
                call prepare_offsets(nd,1,structure2,miss_offsets,code,origin=origin1)
            end if
        else if (present(structure1)) then
            if (present(origin2)) then
                call prepare_offsets(nd,1,structure1,miss_offsets,code,complement=.true.,origin=origin2)
            else
                call prepare_offsets(nd,1,structure1,miss_offsets,code,complement=.true.,origin=origin1)
            end if
        else
            ! Default foreground neighborhood is face-connected; misses
            ! are the remaining entries of its 3-by-...-by-3 box.
            call generate_neighbors(nd,nd,all_offsets,code)
            if (code==morphology_status_success) then
                allocate(temporary(nd,size(all_offsets,2)))
                k=0
                do j=1,size(all_offsets,2)
                    if (count(all_offsets(:,j)/=0)<=1) cycle
                    k=k+1
                    temporary(:,k)=all_offsets(:,j)
                end do
                miss_offsets=temporary(:,:k)
                if (present(origin2)) then
                    if (size(origin2)/=nd) then
                        code=morphology_status_invalid_input
                    else if (any(origin2<(-1)) .or. any(origin2>1)) then
                        code=morphology_status_invalid_input
                    else
                        do j=1,size(miss_offsets,2)
                            miss_offsets(:,j)=miss_offsets(:,j)-origin2
                        end do
                    end if
                else if (present(origin1)) then
                    if (size(origin1)/=nd) then
                        code=morphology_status_invalid_input
                    else
                        do j=1,size(miss_offsets,2)
                            miss_offsets(:,j)=miss_offsets(:,j)-origin1
                        end do
                    end if
                end if
            end if
        end if
        if (code/=morphology_status_success) then
            if (present(status)) status=code
            return
        end if
        call pack_logical(x,flat)
        allocate(hits(size(flat)),misses(size(flat)))
        call morphological_kernel(flat,shape(x),hit_offsets,1,.false.,.false.,hits)
        ! A miss is a pixel hit by a dilation using the background pattern.
        call morphological_kernel(flat,shape(x),-miss_offsets,1,.false.,.true.,misses)
        call scatter_logical(hits.and.(.not.misses),y)
        if (present(status)) status=morphology_status_success
    end subroutine binary_hit_or_miss_nd

    ! Generate the (3,)*rank face/full-connectivity footprint.
    ! Connectivity is clamped to [1,rank], as in scipy.ndimage.
    subroutine generate_binary_structure_nd(connectivity,structure,status)
        integer, intent(in) :: connectivity !! number of coordinate axes allowed to differ
        logical, intent(out) :: structure(..) !! output with three samples in each dimension
        integer, optional, intent(out) :: status !! validation or memory-limit status
        logical, allocatable :: flat(:)
        integer :: nd,degree,j,d,n
        integer(int64) :: rem

        if (present(status)) status=morphology_status_invalid_input
        nd=rank(structure)
        if (nd<1 .or. nd>15) return
        if (any(shape(structure)/=3)) return
        if (size(structure,kind=int64)>20000000_int64) then
            if (present(status)) status=morphology_status_too_expensive
            return
        end if
        degree=max(1,min(nd,connectivity))
        allocate(flat(size(structure)))
        do j=1,size(flat)
            rem=int(j-1,int64)
            n=0
            do d=1,nd
                if (modulo(rem,3_int64)/=1_int64) n=n+1
                rem=rem/3_int64
            end do
            flat(j)=(n<=degree)
        end do
        call scatter_logical(flat,structure)
        if (present(status)) status=morphology_status_success
    end subroutine generate_binary_structure_nd

    ! Repeated Minkowski addition of the set of active footprint coordinates.
    ! This is the support of an iterations-fold binary self-dilation.
    subroutine iterate_structure_nd(structure,iterations,expanded,status,origin,new_origin)
        logical, intent(in) :: structure(..) !! original footprint, ranks 1 through 15
        integer, intent(in) :: iterations !! iteration count, values below two copy the input
        logical, intent(out) :: expanded(..) !! ((shape(structure)-1)*iterations+1) extent for iterations > 1
        integer, optional, intent(out) :: status !! status code
        integer, optional, intent(in) :: origin(:) !! original origin per dimension
        integer, optional, intent(out) :: new_origin(:) !! scaled output origin per dimension
        logical, allocatable :: active(:),current(:),next(:)
        integer, allocatable :: coordinates(:,:),cur_dims(:),target_dims(:)
        integer(int64) :: target_len,remaining,index,steps,work_count,origin_scale
        integer(int64) :: stride(rank(structure)),p,rem
        integer :: nd,nactive,d,j,k,t,now_active,source_active
        integer :: loc(rank(structure))

        if (present(status)) status=morphology_status_invalid_input
        nd=rank(structure)
        if (nd<1 .or. nd>15 .or. rank(expanded)/=nd) return
        if (any(shape(structure)<1)) return
        if (present(origin)) then
            if (size(origin)/=nd) return
        end if
        if (present(new_origin)) then
            if (size(new_origin)/=nd) return
        end if
        if (size(structure,kind=int64)>int(huge(1),int64)) return
        steps=max(1_int64,int(iterations,int64))
        target_len=1_int64
        allocate(target_dims(nd),cur_dims(nd))
        do d=1,nd
            remaining=1_int64+(int(size(structure,d),int64)-1_int64)*steps
            if (remaining>int(huge(1),int64) .or. remaining<1_int64) return
            target_dims(d)=int(remaining)
            if (target_len>20000000_int64/remaining) then
                if (present(status)) status=morphology_status_too_expensive
                return
            end if
            target_len=target_len*remaining
        end do
        if (any(shape(expanded)/=target_dims)) return
        if (steps>10000_int64) then
            if (present(status)) status=morphology_status_too_expensive
            return
        end if
        origin_scale=steps
        if (present(origin)) then
            if (any(abs(int(origin,int64))>int(huge(1),int64)/origin_scale)) return
            if (any(abs(int(origin,int64)*origin_scale)>int(huge(1),int64))) return
        end if
        call pack_logical(structure,active)
        current=active
        cur_dims=shape(structure)
        nactive=count(active)
        allocate(coordinates(nd,nactive))
        j=0
        do k=1,size(active)
            if (.not.active(k)) cycle
            j=j+1
            rem=int(k-1,int64)
            do d=1,nd
                coordinates(d,j)=int(modulo(rem,int(cur_dims(d),int64)))
                rem=rem/int(cur_dims(d),int64)
            end do
        end do
        work_count=0_int64
        do t=2,int(steps)
            do d=1,nd
                target_dims(d)=1+(cur_dims(d)-1)+(size(structure,d)-1)
            end do
            target_len=product(int(target_dims,int64))
            allocate(next(int(target_len)))
            next=.false.
            source_active=count(current)
            if (source_active>0 .and. nactive>0) then
                if (int(source_active,int64)*int(nactive,int64)>30000000_int64-work_count) then
                    if (present(status)) status=morphology_status_too_expensive
                    return
                end if
                work_count=work_count+int(source_active,int64)*int(nactive,int64)
                stride(1)=1_int64
                do d=2,nd
                    stride(d)=stride(d-1)*int(target_dims(d-1),int64)
                end do
                do p=1_int64,size(current,kind=int64)
                    if (.not.current(p)) cycle
                    rem=p-1_int64
                    do d=1,nd
                        loc(d)=int(modulo(rem,int(cur_dims(d),int64)))
                        rem=rem/int(cur_dims(d),int64)
                    end do
                    do j=1,nactive
                        index=1_int64
                        do d=1,nd
                            index=index+int(loc(d)+coordinates(d,j),int64)*stride(d)
                        end do
                        next(index)=.true.
                    end do
                end do
            end if
            call move_alloc(next,current)
            cur_dims=target_dims
        end do
        call scatter_logical(current,expanded)
        if (present(new_origin)) then
            new_origin=0
            if (present(origin)) new_origin=int(int(origin,int64)*origin_scale)
        end if
        if (present(status)) status=morphology_status_success
    end subroutine iterate_structure_nd

    ! Numeric image labeling: any nonzero scalar is foreground.
    subroutine label_integer_nd(x,y,num_features,status,connectivity,structure,axes)
        integer, intent(in) :: x(..) !! integer source image
        integer, intent(out) :: y(..) !! connected-component labels
        integer, intent(out) :: num_features !! component count
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: connectivity !! 1 to rank
        logical, optional, intent(in) :: structure(..) !! compact or full structuring element
        integer, optional, intent(in) :: axes(:) !! one-based connectivity axes
        num_features=0
        if (present(status)) status=morphology_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15) return
        select rank(x)
        rank(1)
            call label_logical_nd(x/=0,y,num_features,status,connectivity,structure,axes)
        rank(2)
            call label_logical_nd(x/=0,y,num_features,status,connectivity,structure,axes)
        rank(3)
            call label_logical_nd(x/=0,y,num_features,status,connectivity,structure,axes)
        rank(4)
            call label_logical_nd(x/=0,y,num_features,status,connectivity,structure,axes)
        rank(5)
            call label_logical_nd(x/=0,y,num_features,status,connectivity,structure,axes)
        rank(6)
            call label_logical_nd(x/=0,y,num_features,status,connectivity,structure,axes)
        rank(7)
            call label_logical_nd(x/=0,y,num_features,status,connectivity,structure,axes)
        rank(8)
            call label_logical_nd(x/=0,y,num_features,status,connectivity,structure,axes)
        rank(9)
            call label_logical_nd(x/=0,y,num_features,status,connectivity,structure,axes)
        rank(10)
            call label_logical_nd(x/=0,y,num_features,status,connectivity,structure,axes)
        rank(11)
            call label_logical_nd(x/=0,y,num_features,status,connectivity,structure,axes)
        rank(12)
            call label_logical_nd(x/=0,y,num_features,status,connectivity,structure,axes)
        rank(13)
            call label_logical_nd(x/=0,y,num_features,status,connectivity,structure,axes)
        rank(14)
            call label_logical_nd(x/=0,y,num_features,status,connectivity,structure,axes)
        rank(15)
            call label_logical_nd(x/=0,y,num_features,status,connectivity,structure,axes)
        rank default
        end select
    end subroutine label_integer_nd

    subroutine label_real_nd(x,y,num_features,status,connectivity,structure,axes)
        real(dp), intent(in) :: x(..) !! real source, nonzero is foreground
        integer, intent(out) :: y(..) !! connected-component labels
        integer, intent(out) :: num_features !! component count
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: connectivity !! 1 to rank
        logical, optional, intent(in) :: structure(..) !! compact or full structuring element
        integer, optional, intent(in) :: axes(:) !! one-based connectivity axes
        num_features=0
        if (present(status)) status=morphology_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15) return
        select rank(x)
        rank(1)
            call label_logical_nd(x/=0.0_dp,y,num_features,status,connectivity,structure,axes)
        rank(2)
            call label_logical_nd(x/=0.0_dp,y,num_features,status,connectivity,structure,axes)
        rank(3)
            call label_logical_nd(x/=0.0_dp,y,num_features,status,connectivity,structure,axes)
        rank(4)
            call label_logical_nd(x/=0.0_dp,y,num_features,status,connectivity,structure,axes)
        rank(5)
            call label_logical_nd(x/=0.0_dp,y,num_features,status,connectivity,structure,axes)
        rank(6)
            call label_logical_nd(x/=0.0_dp,y,num_features,status,connectivity,structure,axes)
        rank(7)
            call label_logical_nd(x/=0.0_dp,y,num_features,status,connectivity,structure,axes)
        rank(8)
            call label_logical_nd(x/=0.0_dp,y,num_features,status,connectivity,structure,axes)
        rank(9)
            call label_logical_nd(x/=0.0_dp,y,num_features,status,connectivity,structure,axes)
        rank(10)
            call label_logical_nd(x/=0.0_dp,y,num_features,status,connectivity,structure,axes)
        rank(11)
            call label_logical_nd(x/=0.0_dp,y,num_features,status,connectivity,structure,axes)
        rank(12)
            call label_logical_nd(x/=0.0_dp,y,num_features,status,connectivity,structure,axes)
        rank(13)
            call label_logical_nd(x/=0.0_dp,y,num_features,status,connectivity,structure,axes)
        rank(14)
            call label_logical_nd(x/=0.0_dp,y,num_features,status,connectivity,structure,axes)
        rank(15)
            call label_logical_nd(x/=0.0_dp,y,num_features,status,connectivity,structure,axes)
        rank default
        end select
    end subroutine label_real_nd

end module scifort_ndimage_morphology
