! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Exact separable Euclidean distance transform using lower envelopes of
! sampled one-dimensional squared-distance parabolas.
module scifort_ndimage_distance
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private
    public :: distance_transform_edt_nd, distance_transform_bf_nd
    real(dp), parameter :: no_site=huge(1.0_dp)/32.0_dp
contains
    subroutine pack_binary(x,flat)
        logical, intent(in) :: x(..) !! input mask
        logical, allocatable, intent(out) :: flat(:) !! flattened mask
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
    end subroutine pack_binary

    subroutine scatter_distances(flat,y)
        real(dp), intent(in) :: flat(:) !! flattened distance array
        real(dp), intent(out) :: y(..) !! real image with original shape
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
    end subroutine scatter_distances

    subroutine squared_distance_line(f,spacing2,g,winner)
        real(dp), intent(in) :: f(:) !! squared distances along a grid line
        real(dp), intent(in) :: spacing2 !! squared physical spacing
        real(dp), intent(out) :: g(:) !! minimum squared distance
        integer, optional, intent(out) :: winner(:) !! one-based minimizing site in this line
        integer, allocatable :: sites(:)
        real(dp), allocatable :: cross(:)
        integer :: q,k,current,j,n
        real(dp) :: s

        n=size(f)
        allocate(sites(n),cross(n+1))
        k=0
        do q=0,n-1
            if (f(q+1)==no_site) cycle
            if (k==0) then
                k=1
                sites(k)=q
                cross(1)=-huge(1.0_dp)
                cross(2)=huge(1.0_dp)
                cycle
            end if
            j=sites(k)
            s=((f(q+1)-f(j+1))/spacing2+real(q,dp)**2-real(j,dp)**2) / &
                (2.0_dp*real(q-j,dp))
            do while(k>1)
                if (s>cross(k)) exit
                k=k-1
                j=sites(k)
                s=((f(q+1)-f(j+1))/spacing2+real(q,dp)**2-real(j,dp)**2) / &
                    (2.0_dp*real(q-j,dp))
            end do
            k=k+1
            sites(k)=q
            cross(k)=s
            cross(k+1)=huge(1.0_dp)
        end do
        if (k==0) then
            g=no_site
            if (present(winner)) winner=0
            return
        end if
        current=1
        do q=0,n-1
            do while(current<k)
                if (real(q,dp)<=cross(current+1)) exit
                current=current+1
            end do
            j=sites(current)
            g(q+1)=f(j+1)+spacing2*real(q-j,dp)**2
            if (present(winner)) winner(q+1)=j+1
        end do
    end subroutine squared_distance_line

    subroutine distance_transform_edt_nd(x,y,status,sampling,feature_indices)
        logical, intent(in) :: x(..) !! foreground (nonzero) pixels, ranks 1..15
        real(dp), intent(out) :: y(..) !! distance to nearest zero-valued pixel
        integer, optional, intent(out) :: status !! zero on success, one on invalid input
        real(dp), optional, intent(in) :: sampling(:) !! physical pixel sizes along each axis
        integer, optional, intent(out) :: feature_indices(:,:) !! (rank, total size), one-based nearest background
        integer :: dims(rank(x)),n,d,i
        integer(int64) :: stride,outer,block,o,j,pos,rem
        real(dp), allocatable :: spacing(:),distance2(:),f(:),g(:),result(:)
        integer, allocatable :: features(:,:),line_features(:,:),winning(:)
        logical, allocatable :: flat(:)
        real(dp) :: total

        if (present(status)) status=1
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x)) return
        if (any(shape(y)/=shape(x)) .or. any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        if (present(feature_indices)) then
            if (size(feature_indices,1)/=rank(x)) return
            if (size(feature_indices,2)/=size(x)) return
        end if
        dims=shape(x)
        allocate(spacing(rank(x)))
        spacing=1.0_dp
        if (present(sampling)) then
            if (size(sampling)/=rank(x)) return
            spacing=sampling
        end if
        if (.not.all(ieee_is_finite(spacing))) return
        if (any(spacing<=0.0_dp) .or. any(spacing>1.0e120_dp)) return
        do d=1,size(spacing)
            if (spacing(d)*spacing(d)==0.0_dp) return
        end do
        call pack_binary(x,flat)
        allocate(distance2(size(flat)),result(size(flat)))
        if (present(feature_indices)) then
            allocate(features(rank(x),size(flat)))
            features=0
            do pos=1_int64,size(flat,kind=int64)
                if (flat(pos)) cycle
                rem=pos-1_int64
                do d=1,size(dims)
                    features(d,pos)=int(modulo(rem,int(dims(d),int64)))+1
                    rem=rem/int(dims(d),int64)
                end do
            end do
        end if
        if (all(flat)) then
            ! SciPy's all-foreground convention behaves as if there were a
            ! single virtual background site at (-1,0,...,0).
            if (present(feature_indices)) then
                features=1
                features(1,:)=0
            end if
            do pos=1_int64,size(flat,kind=int64)
                rem=pos-1_int64
                total=0.0_dp
                do d=1,size(dims)
                    i=int(modulo(rem,int(dims(d),int64)))
                    rem=rem/int(dims(d),int64)
                    if (d==1) i=i+1
                    total=total+(real(i,dp)*spacing(d))**2
                end do
                result(pos)=sqrt(total)
            end do
        else
            where(flat)
                distance2=no_site
            elsewhere
                distance2=0.0_dp
            end where
            do d=1,size(dims)
                n=dims(d)
                allocate(f(n),g(n))
                if (present(feature_indices)) allocate(winning(n),line_features(size(dims),n))
                stride=product(int(dims(:d-1),int64))
                outer=product(int(dims(d+1:),int64))
                block=stride*int(n,int64)
                do o=0_int64,outer-1_int64
                    do j=0_int64,stride-1_int64
                        do i=0,n-1
                            pos=o*block+j+int(i,int64)*stride+1_int64
                            f(i+1)=distance2(pos)
                        end do
                        if (present(feature_indices)) then
                            do i=0,n-1
                                pos=o*block+j+int(i,int64)*stride+1_int64
                                line_features(:,i+1)=features(:,pos)
                            end do
                            call squared_distance_line(f,spacing(d)**2,g,winning)
                        else
                            call squared_distance_line(f,spacing(d)**2,g)
                        end if
                        do i=0,n-1
                            pos=o*block+j+int(i,int64)*stride+1_int64
                            distance2(pos)=g(i+1)
                            if (present(feature_indices)) then
                                if (winning(i+1)>0) then
                                    features(:,pos)=line_features(:,winning(i+1))
                                else
                                    features(:,pos)=0
                                end if
                            end if
                        end do
                    end do
                end do
                deallocate(f,g)
                if (present(feature_indices)) deallocate(winning,line_features)
            end do
            result=sqrt(distance2)
        end if
        call scatter_distances(result,y)
        if (present(feature_indices)) feature_indices=features
        if (present(status)) status=0
    end subroutine distance_transform_edt_nd
    ! Deliberately direct N x B search, useful as an independent check of the
    ! faster separable and chamfer transforms (N pixels, B background sites).
    subroutine distance_transform_bf_nd(x,y,status,metric,sampling,feature_indices)
        logical, intent(in) :: x(..) !! foreground (true) / background (false), ranks 1..15
        real(dp), intent(out) :: y(..) !! distances, including integer-valued L1/L-infinity metrics
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 excessive work
        character(len=*), optional, intent(in) :: metric !! euclidean, taxicab, chessboard (aliases allowed)
        real(dp), optional, intent(in) :: sampling(:) !! Euclidean sampling per dimension
        integer, optional, intent(out) :: feature_indices(:,:) !! one-based nearest site (rank, size(x))
        integer(int64), parameter :: max_pairs=30000000_int64
        integer(int64), parameter :: max_uint32=4294967295_int64
        real(dp), allocatable :: distances(:),spacing(:)
        logical, allocatable :: flat(:)
        integer, allocatable :: sites(:,:),coords(:)
        integer :: dims(rank(x)),n,nb,d,k,b,site,naxes
        integer(int64) :: rem
        character(len=:), allocatable :: choice
        real(dp) :: delta,value,best

        if (present(status)) status=1
        naxes=rank(x)
        if (naxes<1 .or. naxes>15 .or. rank(y)/=naxes) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        n=size(x)
        dims=shape(x)
        if (present(feature_indices)) then
            if (size(feature_indices,1)/=naxes) return
            if (size(feature_indices,2)/=n) return
        end if
        choice='euclidean'
        if (present(metric)) choice=trim(metric)
        select case(choice)
        case('euclidean','taxicab','chessboard')
        case('cityblock','manhattan')
            choice='taxicab'
        case default
            return
        end select
        allocate(spacing(naxes))
        spacing=1.0_dp
        if (present(sampling)) then
            if (size(sampling)/=naxes) return
            if (.not.all(ieee_is_finite(sampling))) return
            if (choice=='euclidean') then
                if (any(sampling<=0.0_dp) .or. any(sampling>1.0e120_dp)) return
                spacing=sampling
                if (any(spacing**2==0.0_dp)) return
            end if
        end if
        call pack_binary(x,flat)
        nb=count(.not.flat)
        allocate(distances(n))
        if (present(feature_indices)) feature_indices=0
        if (nb==0) then
            ! SciPy's brute-force routine reports an infinite-distance
            ! sentinel for images with no background, unlike its EDT routine.
            if (choice=='euclidean') then
                distances=sqrt(huge(1.0_dp))
            else
                distances=real(max_uint32,dp)
            end if
        else
            if (int(n,int64)*int(nb,int64)>max_pairs) then
                if (present(status)) status=2
                return
            end if
            allocate(sites(naxes,nb),coords(naxes))
            b=0
            do k=1,n
                if (flat(k)) cycle
                b=b+1
                rem=int(k,int64)-1_int64
                do d=1,naxes
                    sites(d,b)=int(modulo(rem,int(dims(d),int64)))+1
                    rem=rem/int(dims(d),int64)
                end do
            end do
            do k=1,n
                rem=int(k,int64)-1_int64
                do d=1,naxes
                    coords(d)=int(modulo(rem,int(dims(d),int64)))+1
                    rem=rem/int(dims(d),int64)
                end do
                if (.not.flat(k)) then
                    distances(k)=0.0_dp
                    if (present(feature_indices)) feature_indices(:,k)=coords
                    cycle
                end if
                best=huge(1.0_dp)
                site=0
                do b=1,nb
                    value=0.0_dp
                    do d=1,naxes
                        delta=real(abs(coords(d)-sites(d,b)),dp)
                        select case(choice)
                        case('euclidean')
                            value=value+(delta*spacing(d))**2
                        case('taxicab')
                            value=value+delta
                        case('chessboard')
                            value=max(value,delta)
                        end select
                    end do
                    if (value<best) then
                        best=value
                        site=b
                    end if
                end do
                if (choice=='euclidean') best=sqrt(best)
                distances(k)=best
                if (present(feature_indices)) feature_indices(:,k)=sites(:,site)
            end do
        end if
        call scatter_distances(distances,y)
        if (present(status)) status=0
    end subroutine distance_transform_bf_nd
end module scifort_ndimage_distance
