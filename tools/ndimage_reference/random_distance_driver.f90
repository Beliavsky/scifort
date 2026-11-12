! SPDX-License-Identifier: MIT
! Bulk regression driver for randomized distance transforms (rank 1 to 3).
program random_distance_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : distance_transform_edt_nd
    implicit none
    integer :: rd,wr,r,ncase,which,status,dims(3),k
    integer, allocatable :: values(:)
    real(dp) :: spacing(3)
    logical, allocatable :: a1(:),a2(:,:),a3(:,:,:)
    real(dp), allocatable :: y1(:),y2(:,:),y3(:,:,:)

    open(newunit=rd,file='random_distance_input.bin',form='unformatted', &
        access='stream',status='old',action='read')
    open(newunit=wr,file='random_distance_output.bin',form='unformatted', &
        access='stream',status='replace',action='write')
    read(rd) ncase
    do which=1,ncase
        read(rd) r,dims,spacing
        allocate(values(product(dims)))
        read(rd) values
        select case(r)
        case(1)
            allocate(a1(dims(1)),y1(dims(1)))
            a1=values/=0
            call distance_transform_edt_nd(a1,y1,status,spacing(:r))
            if (status/=0) error stop 'random EDT rank1'
            write(wr) y1
            deallocate(a1,y1)
        case(2)
            allocate(a2(dims(1),dims(2)),y2(dims(1),dims(2)))
            a2=reshape(values/=0,shape(a2))
            call distance_transform_edt_nd(a2,y2,status,spacing(:r))
            if (status/=0) error stop 'random EDT rank2'
            write(wr) y2
            deallocate(a2,y2)
        case(3)
            allocate(a3(dims(1),dims(2),dims(3)),y3(dims(1),dims(2),dims(3)))
            a3=reshape(values/=0,shape(a3))
            call distance_transform_edt_nd(a3,y3,status,spacing(:r))
            if (status/=0) error stop 'random EDT rank3'
            write(wr) y3
            deallocate(a3,y3)
        case default
            error stop 'invalid rank'
        end select
        deallocate(values)
    end do
    close(rd)
    close(wr)
end program random_distance_driver
