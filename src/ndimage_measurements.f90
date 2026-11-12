! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Selected-label statistics for real images of ranks 1 through 15.
module scifort_ndimage_measurements
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, &
        ieee_quiet_nan, ieee_positive_inf
    use scifort_kinds, only : dp
    implicit none
    private
    public :: sum_labels_nd, mean_labels_nd, variance_labels_nd, center_of_mass_nd
    public :: standard_deviation_labels_nd, count_labels_nd, histogram_labels_nd
    public :: labeled_comprehension_nd, ndimage_region_reduce
    public :: extrema_labels_nd, minimum_labels_nd, maximum_labels_nd
    public :: minimum_position_nd, maximum_position_nd
    public :: bounding_boxes_nd, find_objects_nd, value_indices_nd
    abstract interface
        function ndimage_region_reduce(values) result(value)
            import dp
            real(dp), intent(in) :: values(:) !! pixels selected by one label, C array order
            real(dp) :: value !! reduced scalar value
        end function ndimage_region_reduce
    end interface
contains
    subroutine pack_real(x,flat)
        real(dp), intent(in) :: x(..) !! scalar input field, arbitrary rank
        real(dp), allocatable, intent(out) :: flat(:) !! flattened field
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
    end subroutine pack_real

    subroutine pack_integer(x,flat)
        integer, intent(in) :: x(..) !! integer labels
        integer, allocatable, intent(out) :: flat(:) !! flattened labels
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

    ! SciPy's value_indices mapping in compressed coordinate-list form:
    ! values(k) owns coordinates(:,starts(k):starts(k+1)-1). Keys are sorted,
    ! locations within a key are in C-order, and coordinates are one-based.
    ! The two index vectors make arbitrary signed integer labels efficient:
    ! sorting costs O(N log N), without allocating an array of max(label).
    subroutine value_indices_nd(x,values,starts,coordinates,status,ignore_value)
        integer, intent(in) :: x(..) !! integer image, ranks 1 through 15
        integer, allocatable, intent(out) :: values(:) !! ascending distinct values
        integer, allocatable, intent(out) :: starts(:) !! 1-based CSR starts, one past last at end
        integer, allocatable, intent(out) :: coordinates(:,:) !! (rank(x), number of included pixels)
        integer, optional, intent(out) :: status !! 0 success, 1 invalid, 2 resource limit
        integer, optional, intent(in) :: ignore_value !! omit this integer key entirely
        integer, allocatable :: flat(:),keys(:),order(:)
        integer :: dims(rank(x)),nd,n,m,d,j,k,groups,idx
        integer(int64) :: rem,linear,step(rank(x))

        if (present(status)) status=1
        nd=rank(x)
        if (nd<1 .or. nd>15) return
        dims=shape(x)
        if (any(dims<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)-1_int64) then
            if (present(status)) status=2
            return
        end if
        n=size(x)
        call pack_integer(x,flat)
        allocate(keys(n),order(n))
        step(1)=1_int64
        do d=2,nd
            step(d)=step(d-1)*int(dims(d-1),int64)
        end do
        m=0
        ! Traverse C-order (last dimension fastest) before sorting values.
        do j=1,n
            rem=int(j-1,int64)
            linear=1_int64
            do d=nd,1,-1
                linear=linear+modulo(rem,int(dims(d),int64))*step(d)
                rem=rem/int(dims(d),int64)
            end do
            idx=flat(linear)
            if (present(ignore_value)) then
                if (idx==ignore_value) cycle
            end if
            m=m+1
            keys(m)=idx
            order(m)=j
        end do
        ! Heap sort by (value, C-order position), with constant stack size.
        do j=m/2,1,-1
            call sift(j,m)
        end do
        do j=m,2,-1
            call swap(1,j)
            call sift(1,j-1)
        end do
        groups=0
        do j=1,m
            if (j==1) then
                groups=groups+1
            else if (keys(j)/=keys(j-1)) then
                groups=groups+1
            end if
        end do
        allocate(values(groups),starts(groups+1),coordinates(nd,m))
        starts(groups+1)=m+1
        k=0
        do j=1,m
            if (j==1) then
                k=k+1
                values(k)=keys(j)
                starts(k)=j
            else if (keys(j)/=keys(j-1)) then
                k=k+1
                values(k)=keys(j)
                starts(k)=j
            end if
            rem=int(order(j)-1,int64)
            do d=nd,1,-1
                coordinates(d,j)=int(modulo(rem,int(dims(d),int64)))+1
                rem=rem/int(dims(d),int64)
            end do
        end do
        if (present(status)) status=0
    contains
        pure logical function greater(i,j)
            integer, intent(in) :: i !! first position to compare
            integer, intent(in) :: j !! second position to compare
            greater=(keys(i)>keys(j))
            if (keys(i)==keys(j)) greater=order(i)>order(j)
        end function greater

        subroutine swap(i,j)
            integer, intent(in) :: i !! first exchange position
            integer, intent(in) :: j !! second exchange position
            integer :: t
            t=keys(i)
            keys(i)=keys(j)
            keys(j)=t
            t=order(i)
            order(i)=order(j)
            order(j)=t
        end subroutine swap

        subroutine sift(first,last)
            integer, intent(in) :: first !! current heap root
            integer, intent(in) :: last !! final heap item
            integer :: parent,child
            parent=first
            do while(parent<=last/2)
                child=2*parent
                if (child<last) then
                    if (greater(child+1,child)) child=child+1
                end if
                if (.not.greater(child,parent)) exit
                call swap(parent,child)
                parent=child
            end do
        end subroutine sift
    end subroutine value_indices_nd

    ! Sort requested labels once, retaining the original output positions.
    ! A binary search per pixel avoids scanning every requested label. Repeated
    ! IDs are accumulated at the first sorted request, then copied afterward.
    subroutine label_lookup_order(index,order)
        integer, intent(in) :: index(:) !! arbitrary signed selected IDs
        integer, allocatable, intent(out) :: order(:) !! permutation sorted by label and original position
        integer :: n,j,last,parent,child,tmp
        n=size(index)
        allocate(order(n))
        order=[(j,j=1,n)]
        do j=n/2,1,-1
            call sift(j,n)
        end do
        do j=n,2,-1
            tmp=order(1)
            order(1)=order(j)
            order(j)=tmp
            call sift(1,j-1)
        end do
    contains
        pure logical function greater(a,b)
            integer, intent(in) :: a !! first selected-index position
            integer, intent(in) :: b !! second selected-index position
            greater=index(a)>index(b)
            if (index(a)==index(b)) greater=a>b
        end function greater

        subroutine sift(start,last)
            integer, intent(in) :: start !! current root
            integer, intent(in) :: last !! last heap position
            parent=start
            do while (parent<=last/2)
                child=2*parent
                if (child<last) then
                    if (greater(order(child+1),order(child))) child=child+1
                end if
                if (.not.greater(order(child),order(parent))) exit
                tmp=order(parent)
                order(parent)=order(child)
                order(child)=tmp
                parent=child
            end do
        end subroutine sift
    end subroutine label_lookup_order

    pure integer function lookup_label(index,order,key) result(j)
        integer, intent(in) :: index(:) !! selected labels
        integer, intent(in) :: order(:) !! permutation sorted by label
        integer, intent(in) :: key !! label in flattened image
        integer :: low,high,mid
        low=1
        high=size(order)+1
        do while (low<high)
            mid=low+(high-low)/2
            if (index(order(mid))<key) then
                low=mid+1
            else
                high=mid
            end if
        end do
        j=0
        if (low<=size(order)) then
            if (index(order(low))==key) j=order(low)
        end if
    end function lookup_label

    subroutine label_statistics(x,labels,index,totals,counts,means,m2,first_moments,status)
        real(dp), intent(in) :: x(..) !! scalar image, ranks 1 through 15
        integer, intent(in) :: labels(..) !! same shape as x
        integer, intent(in) :: index(:) !! label values requested, including zero
        real(dp), intent(out) :: totals(:) !! sums of pixel values per label
        integer, intent(out) :: counts(:) !! number of pixels per label
        real(dp), intent(out) :: means(:) !! mean values per label
        real(dp), intent(out) :: m2(:) !! sum of squared centered deviations
        real(dp), intent(out) :: first_moments(:,:) !! weighted one-based coordinates
        integer, intent(out) :: status !! zero if input is valid
        real(dp), allocatable :: data(:)
        integer, allocatable :: ids(:),order(:)
        integer :: dims(rank(x)),coords(rank(x)),nd,i,k,j,found
        integer(int64) :: pos,rem
        real(dp) :: delta

        status=1
        nd=rank(x)
        if (nd<1 .or. nd>15 .or. rank(labels)/=nd) return
        if (any(shape(x)/=shape(labels)) .or. any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        if (size(index)/=size(totals) .or. size(index)/=size(counts)) return
        if (size(index)/=size(means) .or. size(index)/=size(m2)) return
        if (size(first_moments,1)/=nd .or. size(first_moments,2)/=size(index)) return
        call pack_real(x,data)
        if (.not.all(ieee_is_finite(data))) return
        call pack_integer(labels,ids)
        totals=0.0_dp
        counts=0
        means=0.0_dp
        m2=0.0_dp
        first_moments=0.0_dp
        dims=shape(x)
        call label_lookup_order(index,order)
        do pos=1_int64,size(data,kind=int64)
            found=lookup_label(index,order,ids(pos))
            if (found==0) cycle
            rem=pos-1_int64
            do i=1,nd
                coords(i)=int(modulo(rem,int(dims(i),int64)))+1
                rem=rem/int(dims(i),int64)
            end do
            if (counts(found)==huge(1)) return
            counts(found)=counts(found)+1
            totals(found)=totals(found)+data(pos)
            delta=data(pos)-means(found)
            means(found)=means(found)+delta/real(counts(found),dp)
            m2(found)=m2(found)+delta*(data(pos)-means(found))
            do k=1,nd
                first_moments(k,found)=first_moments(k,found)+data(pos)*real(coords(k),dp)
            end do
        end do
        do j=2,size(order)
            if (index(order(j))/=index(order(j-1))) cycle
            found=order(j-1)
            i=order(j)
            counts(i)=counts(found)
            totals(i)=totals(found)
            means(i)=means(found)
            m2(i)=m2(found)
            first_moments(:,i)=first_moments(:,found)
        end do
        status=0
    end subroutine label_statistics

    subroutine sum_labels_nd(x,labels,index,totals,status)
        real(dp), intent(in) :: x(..) !! image samples
        integer, intent(in) :: labels(..) !! integer labeled regions
        integer, intent(in) :: index(:) !! selected integer label IDs
        real(dp), intent(out) :: totals(:) !! sums for the corresponding IDs
        integer, optional, intent(out) :: status !! success 0, invalid input 1
        integer, allocatable :: counts(:)
        real(dp), allocatable :: means(:),m2(:),moments(:,:)
        integer :: code
        if (present(status)) status=1
        allocate(counts(size(index)),means(size(index)),m2(size(index)))
        allocate(moments(rank(x),size(index)))
        call label_statistics(x,labels,index,totals,counts,means,m2,moments,code)
        if (present(status)) status=code
    end subroutine sum_labels_nd

    subroutine mean_labels_nd(x,labels,index,values,status)
        real(dp), intent(in) :: x(..) !! image samples
        integer, intent(in) :: labels(..) !! integer labeled regions
        integer, intent(in) :: index(:) !! selected label IDs
        real(dp), intent(out) :: values(:) !! arithmetic means
        integer, optional, intent(out) :: status !! success 0, invalid input 1
        integer, allocatable :: counts(:)
        real(dp), allocatable :: totals(:),m2(:),moments(:,:)
        integer :: code,j
        if (present(status)) status=1
        allocate(counts(size(index)),totals(size(index)),m2(size(index)))
        allocate(moments(rank(x),size(index)))
        call label_statistics(x,labels,index,totals,counts,values,m2,moments,code)
        if (code/=0) return
        do j=1,size(index)
            if (counts(j)==0) values(j)=ieee_value(0.0_dp,ieee_quiet_nan)
        end do
        if (present(status)) status=0
    end subroutine mean_labels_nd

    subroutine variance_labels_nd(x,labels,index,values,status)
        real(dp), intent(in) :: x(..) !! image samples
        integer, intent(in) :: labels(..) !! integer labeled regions
        integer, intent(in) :: index(:) !! selected label IDs
        real(dp), intent(out) :: values(:) !! population variance (divisor n)
        integer, optional, intent(out) :: status !! success 0, invalid input 1
        integer, allocatable :: counts(:)
        real(dp), allocatable :: totals(:),means(:),moments(:,:)
        integer :: code,j
        if (present(status)) status=1
        allocate(counts(size(index)),totals(size(index)),means(size(index)))
        allocate(moments(rank(x),size(index)))
        call label_statistics(x,labels,index,totals,counts,means,values,moments,code)
        if (code/=0) return
        do j=1,size(index)
            if (counts(j)==0) then
                values(j)=ieee_value(0.0_dp,ieee_quiet_nan)
            else
                values(j)=max(0.0_dp,values(j)/real(counts(j),dp))
            end if
        end do
        if (present(status)) status=0
    end subroutine variance_labels_nd

    subroutine center_of_mass_nd(x,labels,index,centers,status)
        real(dp), intent(in) :: x(..) !! mass at each image location (may be negative)
        integer, intent(in) :: labels(..) !! integer labeled regions
        integer, intent(in) :: index(:) !! selected label IDs
        real(dp), intent(out) :: centers(:,:) !! (rank, selected IDs), one-based coordinates
        integer, optional, intent(out) :: status !! success 0, invalid input 1
        integer, allocatable :: counts(:)
        real(dp), allocatable :: totals(:),means(:),m2(:)
        integer :: code,j,k
        if (present(status)) status=1
        allocate(counts(size(index)),totals(size(index)),means(size(index)),m2(size(index)))
        call label_statistics(x,labels,index,totals,counts,means,m2,centers,code)
        if (code/=0) return
        do j=1,size(index)
            if (totals(j)/=0.0_dp) then
                centers(:,j)=centers(:,j)/totals(j)
            else
                do k=1,rank(x)
                    if (centers(k,j)==0.0_dp) then
                        centers(k,j)=ieee_value(0.0_dp,ieee_quiet_nan)
                    else
                        centers(k,j)=sign(ieee_value(0.0_dp,ieee_positive_inf),centers(k,j))
                    end if
                end do
            end if
        end do
        if (present(status)) status=0
    end subroutine center_of_mass_nd
    subroutine standard_deviation_labels_nd(x,labels,index,values,status)
        real(dp), intent(in) :: x(..) !! input real image
        integer, intent(in) :: labels(..) !! labels matching x
        integer, intent(in) :: index(:) !! labels to measure
        real(dp), intent(out) :: values(:) !! population standard deviations
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input
        integer :: code

        if (present(status)) status=1
        call variance_labels_nd(x,labels,index,values,code)
        if (code/=0) return
        values=sqrt(values)
        if (present(status)) status=0
    end subroutine standard_deviation_labels_nd

    subroutine count_labels_nd(labels,index,counts,status)
        integer, intent(in) :: labels(..) !! integer region labels, ranks 1..15
        integer, intent(in) :: index(:) !! selected label IDs (duplicates allowed)
        integer, intent(out) :: counts(:) !! number of pixels in each region
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input
        integer, allocatable :: ids(:),order(:)
        integer :: j,k
        integer(int64) :: p

        if (present(status)) status=1
        if (rank(labels)<1 .or. rank(labels)>15) return
        if (any(shape(labels)<1) .or. size(counts)/=size(index)) return
        if (size(labels,kind=int64)>int(huge(1),int64)) return
        call pack_integer(labels,ids)
        call label_lookup_order(index,order)
        counts=0
        do p=1_int64,size(ids,kind=int64)
            j=lookup_label(index,order,ids(p))
            if (j==0) cycle
            if (counts(j)==huge(1)) return
            counts(j)=counts(j)+1
        end do
        do k=2,size(order)
            if (index(order(k))==index(order(k-1))) counts(order(k))=counts(order(k-1))
        end do
        if (present(status)) status=0
    end subroutine count_labels_nd

    subroutine histogram_labels_nd(x,labels,index,low,high,histograms,status)
        real(dp), intent(in) :: x(..) !! finite real-valued image
        integer, intent(in) :: labels(..) !! integer labels matching x
        integer, intent(in) :: index(:) !! IDs for which histograms are returned
        real(dp), intent(in) :: low !! lower edge, inclusive
        real(dp), intent(in) :: high !! upper edge, inclusive in final bin
        integer, intent(out) :: histograms(:,:) !! (number of bins, number of IDs)
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input
        real(dp), allocatable :: data(:)
        integer, allocatable :: ids(:),order(:)
        integer :: bins,j,b,k
        integer(int64) :: p
        real(dp) :: t

        if (present(status)) status=1
        if (rank(x)<1 .or. rank(x)>15 .or. rank(labels)/=rank(x)) return
        if (any(shape(x)<1) .or. any(shape(x)/=shape(labels))) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        if (size(histograms,2)/=size(index)) return
        bins=size(histograms,1)
        if (bins<1 .or. .not.ieee_is_finite(low) .or. .not.ieee_is_finite(high)) return
        if (high<=low .or. .not.ieee_is_finite(high-low)) return
        call pack_real(x,data)
        if (.not.all(ieee_is_finite(data))) return
        call pack_integer(labels,ids)
        histograms=0
        call label_lookup_order(index,order)
        do p=1_int64,size(data,kind=int64)
            if (data(p)<low .or. data(p)>high) cycle
            t=(data(p)-low)/(high-low)
            if (t>=1.0_dp) then
                b=bins
            else
                b=max(1,min(bins,int(t*real(bins,dp))+1))
            end if
            j=lookup_label(index,order,ids(p))
            if (j==0) cycle
            if (histograms(b,j)==huge(1)) return
            histograms(b,j)=histograms(b,j)+1
        end do
        do k=2,size(order)
            if (index(order(k))==index(order(k-1))) histograms(:,order(k))=histograms(:,order(k-1))
        end do
        if (present(status)) status=0
    end subroutine histogram_labels_nd

    ! Region callback reduction. Bucket each selected pixel in one C-order
    ! traversal, retaining SciPy's coordinate order within every group.
    ! Request lookup costs O(log K) per pixel for K selected IDs. Unlike a
    ! repeated full-image scan, this is O(N log K + K log K) plus callbacks.
    ! Duplicate requests invoke the callback separately, as before.
    subroutine labeled_comprehension_nd(x,labels,index,func,default,values,status)
        real(dp), intent(in) :: x(..) !! finite real image, ranks 1 through 15
        integer, intent(in) :: labels(..) !! same-shape integer region labels
        integer, intent(in) :: index(:) !! selected IDs, possibly repeated or missing
        procedure(ndimage_region_reduce) :: func !! callback for each nonempty request
        real(dp), intent(in) :: default !! return this value when an ID is absent
        real(dp), intent(out) :: values(:) !! one callback result per selected ID
        integer, optional, intent(out) :: status !! 0 success, 1 invalid image or output shape
        real(dp), allocatable :: data(:),samples(:)
        integer, allocatable :: ids(:),counts(:),order(:),next_slot(:),starts(:)
        integer :: dims(rank(x)),nd,j,d,coordinate,k,first
        integer(int64) :: pos,fpos,rem,strides(rank(x)),cpos,total

        if (present(status)) status=1
        nd=rank(x)
        if (nd<1 .or. nd>15 .or. rank(labels)/=nd) return
        if (any(shape(x)<1) .or. any(shape(x)/=shape(labels))) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        if (size(values)/=size(index)) return
        call pack_real(x,data)
        if (.not.all(ieee_is_finite(data))) return
        call pack_integer(labels,ids)
        values=default
        if (size(index)==0) then
            if (present(status)) status=0
            return
        end if

        call label_lookup_order(index,order)
        allocate(counts(size(index)),starts(size(index)),next_slot(size(index)))
        counts=0
        starts=0
        do pos=1_int64,size(ids,kind=int64)
            j=lookup_label(index,order,ids(pos))
            if (j/=0) counts(j)=counts(j)+1
        end do
        total=0_int64
        do k=1,size(order)
            j=order(k)
            if (k>1) then
                if (index(j)==index(order(k-1))) cycle
            end if
            starts(j)=int(total)+1
            total=total+int(counts(j),int64)
        end do
        allocate(samples(int(total)))
        next_slot=starts
        dims=shape(x)
        strides(1)=1_int64
        do d=2,nd
            strides(d)=strides(d-1)*int(dims(d-1),int64)
        end do
        do cpos=0_int64,size(ids,kind=int64)-1_int64
            rem=cpos
            fpos=1_int64
            do d=nd,1,-1
                coordinate=int(modulo(rem,int(dims(d),int64)))
                rem=rem/int(dims(d),int64)
                fpos=fpos+int(coordinate,int64)*strides(d)
            end do
            j=lookup_label(index,order,ids(fpos))
            if (j==0) cycle
            samples(next_slot(j))=data(fpos)
            next_slot(j)=next_slot(j)+1
        end do
        do k=1,size(index)
            first=lookup_label(index,order,index(k))
            if (counts(first)==0) cycle
            j=starts(first)
            values(k)=func(samples(j:j+counts(first)-1))
        end do
        if (present(status)) status=0
    end subroutine labeled_comprehension_nd

    ! Position ties are resolved in C-order (axis 1 first), as in SciPy.
    pure logical function earlier_position(candidate,previous) result(earlier)
        integer, intent(in) :: candidate(:) !! candidate one-based position
        integer, intent(in) :: previous(:) !! previous one-based position
        integer :: d
        earlier=.false.
        do d=1,size(candidate)
            if (candidate(d)<previous(d)) then
                earlier=.true.
                return
            else if (candidate(d)>previous(d)) then
                return
            end if
        end do
    end function earlier_position

    subroutine extrema_labels_nd(x,labels,index,minimums,maximums,min_positions,max_positions,status)
        real(dp), intent(in) :: x(..) !! real image, ranks 1 through 15
        integer, intent(in) :: labels(..) !! integer labels with shape(x)
        integer, intent(in) :: index(:) !! requested label IDs
        real(dp), intent(out) :: minimums(:) !! min values; zero if missing
        real(dp), intent(out) :: maximums(:) !! max values; zero if missing
        integer, intent(out) :: min_positions(:,:) !! (rank,number of IDs), one based, zeros if missing
        integer, intent(out) :: max_positions(:,:) !! (rank,number of IDs), one based, zeros if missing
        integer, optional, intent(out) :: status !! 0 success, 1 invalid shape or value
        real(dp), allocatable :: data(:)
        integer, allocatable :: ids(:),order(:)
        logical, allocatable :: found(:)
        integer :: dims(rank(x)),coord(rank(x)),j,d,k
        integer(int64) :: pos,rem

        if (present(status)) status=1
        if (rank(x)<1 .or. rank(x)>15) return
        if (rank(labels)/=rank(x)) return
        if (any(shape(x)<1)) return
        if (any(shape(labels)/=shape(x))) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        if (size(minimums)/=size(index) .or. size(maximums)/=size(index)) return
        if (size(min_positions,1)/=rank(x) .or. size(min_positions,2)/=size(index)) return
        if (size(max_positions,1)/=rank(x) .or. size(max_positions,2)/=size(index)) return
        call pack_real(x,data)
        if (.not.all(ieee_is_finite(data))) return
        call pack_integer(labels,ids)
        allocate(found(size(index)))
        found=.false.
        minimums=0.0_dp
        maximums=0.0_dp
        min_positions=0
        max_positions=0
        dims=shape(x)
        call label_lookup_order(index,order)
        do pos=1_int64,size(data,kind=int64)
            j=lookup_label(index,order,ids(pos))
            if (j==0) cycle
            rem=pos-1_int64
            do d=1,rank(x)
                coord(d)=int(modulo(rem,int(dims(d),int64)))+1
                rem=rem/int(dims(d),int64)
            end do
            if (.not.found(j)) then
                minimums(j)=data(pos)
                maximums(j)=data(pos)
                min_positions(:,j)=coord
                max_positions(:,j)=coord
                found(j)=.true.
            else
                if (data(pos)<minimums(j)) then
                    minimums(j)=data(pos)
                    min_positions(:,j)=coord
                else if (data(pos)==minimums(j)) then
                    if (earlier_position(coord,min_positions(:,j))) min_positions(:,j)=coord
                end if
                if (data(pos)>maximums(j)) then
                    maximums(j)=data(pos)
                    max_positions(:,j)=coord
                else if (data(pos)==maximums(j)) then
                    if (earlier_position(coord,max_positions(:,j))) max_positions(:,j)=coord
                end if
            end if
        end do
        do k=2,size(order)
            if (index(order(k))/=index(order(k-1))) cycle
            j=order(k)
            d=order(k-1)
            minimums(j)=minimums(d)
            maximums(j)=maximums(d)
            min_positions(:,j)=min_positions(:,d)
            max_positions(:,j)=max_positions(:,d)
        end do
        if (present(status)) status=0
    end subroutine extrema_labels_nd

    subroutine minimum_labels_nd(x,labels,index,values,status)
        real(dp), intent(in) :: x(..) !! input real image
        integer, intent(in) :: labels(..) !! region labels
        integer, intent(in) :: index(:) !! selected labels
        real(dp), intent(out) :: values(:) !! minima for selected labels
        integer, optional, intent(out) :: status !! 0 success, 1 invalid
        real(dp) :: maxima(size(index))
        integer :: pmin(rank(x),size(index)),pmax(rank(x),size(index))
        call extrema_labels_nd(x,labels,index,values,maxima,pmin,pmax,status)
    end subroutine minimum_labels_nd

    subroutine maximum_labels_nd(x,labels,index,values,status)
        real(dp), intent(in) :: x(..) !! input real image
        integer, intent(in) :: labels(..) !! region labels
        integer, intent(in) :: index(:) !! selected labels
        real(dp), intent(out) :: values(:) !! maxima for selected labels
        integer, optional, intent(out) :: status !! 0 success, 1 invalid
        real(dp) :: minima(size(index))
        integer :: pmin(rank(x),size(index)),pmax(rank(x),size(index))
        call extrema_labels_nd(x,labels,index,minima,values,pmin,pmax,status)
    end subroutine maximum_labels_nd

    subroutine minimum_position_nd(x,labels,index,positions,status)
        real(dp), intent(in) :: x(..) !! input real image
        integer, intent(in) :: labels(..) !! region labels
        integer, intent(in) :: index(:) !! selected labels
        integer, intent(out) :: positions(:,:) !! (rank, ID count), one-based
        integer, optional, intent(out) :: status !! 0 success, 1 invalid
        real(dp) :: minima(size(index)),maxima(size(index))
        integer :: pmax(rank(x),size(index))
        call extrema_labels_nd(x,labels,index,minima,maxima,positions,pmax,status)
    end subroutine minimum_position_nd

    subroutine maximum_position_nd(x,labels,index,positions,status)
        real(dp), intent(in) :: x(..) !! input real image
        integer, intent(in) :: labels(..) !! region labels
        integer, intent(in) :: index(:) !! selected labels
        integer, intent(out) :: positions(:,:) !! (rank, ID count), one-based
        integer, optional, intent(out) :: status !! 0 success, 1 invalid
        real(dp) :: minima(size(index)),maxima(size(index))
        integer :: pmin(rank(x),size(index))
        call extrema_labels_nd(x,labels,index,minima,maxima,pmin,positions,status)
    end subroutine maximum_position_nd

    subroutine bounding_boxes_nd(labels,index,lower,upper,status)
        integer, intent(in) :: labels(..) !! labeled image of rank 1 through 15
        integer, intent(in) :: index(:) !! requested label values, may include zero
        integer, intent(out) :: lower(:,:) !! inclusive one-based lower bound for each region
        integer, intent(out) :: upper(:,:) !! inclusive one-based upper bound for each region
        integer, optional, intent(out) :: status !! 0 success, 1 invalid shape
        integer, allocatable :: ids(:),order(:)
        integer :: dims(rank(labels)),coord(rank(labels)),j,d,k
        integer(int64) :: pos,rem

        if (present(status)) status=1
        if (rank(labels)<1 .or. rank(labels)>15 .or. any(shape(labels)<1)) return
        if (size(labels,kind=int64)>int(huge(1),int64)) return
        if (size(lower,1)/=rank(labels) .or. size(upper,1)/=rank(labels)) return
        if (size(lower,2)/=size(index) .or. size(upper,2)/=size(index)) return
        lower=0
        upper=0
        call pack_integer(labels,ids)
        dims=shape(labels)
        call label_lookup_order(index,order)
        do pos=1_int64,size(ids,kind=int64)
            j=lookup_label(index,order,ids(pos))
            if (j==0) cycle
            rem=pos-1_int64
            do d=1,rank(labels)
                coord(d)=int(modulo(rem,int(dims(d),int64)))+1
                rem=rem/int(dims(d),int64)
            end do
            if (all(lower(:,j)==0)) then
                lower(:,j)=coord
                upper(:,j)=coord
            else
                lower(:,j)=min(lower(:,j),coord)
                upper(:,j)=max(upper(:,j),coord)
            end if
        end do
        do k=2,size(order)
            if (index(order(k))/=index(order(k-1))) cycle
            j=order(k)
            d=order(k-1)
            lower(:,j)=lower(:,d)
            upper(:,j)=upper(:,d)
        end do
        if (present(status)) status=0
    end subroutine bounding_boxes_nd

    subroutine find_objects_nd(labels,lower,upper,status)
        integer, intent(in) :: labels(..) !! labeled integer image
        integer, intent(out) :: lower(:,:) !! (rank,max_label), lower inclusive, missing zero
        integer, intent(out) :: upper(:,:) !! (rank,max_label), upper inclusive, missing zero
        integer, optional, intent(out) :: status !! zero on success
        integer :: j
        integer, allocatable :: index(:)
        if (present(status)) status=1
        if (size(lower,2)/=size(upper,2)) return
        allocate(index(size(lower,2)))
        do j=1,size(index)
            index(j)=j
        end do
        call bounding_boxes_nd(labels,index,lower,upper,status)
    end subroutine find_objects_nd
end module scifort_ndimage_measurements
