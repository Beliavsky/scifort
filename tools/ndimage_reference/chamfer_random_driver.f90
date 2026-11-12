! SPDX-License-Identifier: MIT
! Streaming randomized reference inputs for discrete chamfer transforms.
program chamfer_random_driver
    use scifort_ndimage, only : distance_transform_cdt_nd
    implicit none
    integer, parameter :: cases=160
    integer :: inunit,outunit,n,metric,status,encoded2(5,6),encoded3(4,3,2)
    integer :: result2(5,6),result3(4,3,2),feature2(2,30),feature3(3,24)
    logical :: mask2(5,6),mask3(4,3,2)
    character(len=10) :: method

    open(newunit=inunit,file='chamfer_random_inputs.bin',form='unformatted', &
        access='stream',status='old')
    open(newunit=outunit,file='chamfer_random_results.bin',form='unformatted', &
        access='stream',status='replace')
    do n=1,cases
        read(inunit) encoded2,encoded3
        mask2=encoded2/=0
        mask3=encoded3/=0
        do metric=1,2
            if (metric==1) then
                method='taxicab'
            else
                method='chessboard'
            end if
            call distance_transform_cdt_nd(mask2,result2,status,metric=trim(method),feature_indices=feature2)
            if (status/=0) error stop 'random 2D CDT'
            write(outunit) result2,feature2
            call distance_transform_cdt_nd(mask3,result3,status,metric=trim(method),feature_indices=feature3)
            if (status/=0) error stop 'random 3D CDT'
            write(outunit) result3,feature3
        end do
    end do
    close(inunit)
    close(outunit)
end program chamfer_random_driver
