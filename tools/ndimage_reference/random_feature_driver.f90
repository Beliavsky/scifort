! SPDX-License-Identifier: MIT
! Stream random binary images through EDT feature-index comparison.
program random_feature_driver
    use scifort_kinds, only : dp
    use scifort_ndimage_distance, only : distance_transform_edt_nd
    implicit none
    integer, parameter :: ntest=64
    integer :: input2(5,6),input3(4,3,3),indices2(2,30),indices3(3,36)
    logical :: mask2(5,6),mask3(4,3,3)
    real(dp) :: y2(5,6),y3(4,3,3)
    integer :: unit_in,unit_out,status,k
    open(newunit=unit_in,file='feature_random_inputs.bin',form='unformatted',access='stream',status='old')
    open(newunit=unit_out,file='feature_random_results.bin',form='unformatted',access='stream',status='replace')
    do k=1,ntest
        read(unit_in) input2
        mask2=input2/=0
        call distance_transform_edt_nd(mask2,y2,status,sampling=[1.75_dp,0.6_dp], &
            feature_indices=indices2)
        if (status/=0) error stop 'random EDT 2D status'
        write(unit_out) y2
        write(unit_out) indices2
        read(unit_in) input3
        mask3=input3/=0
        call distance_transform_edt_nd(mask3,y3,status,sampling=[0.85_dp,1.45_dp,2.5_dp], &
            feature_indices=indices3)
        if (status/=0) error stop 'random EDT 3D status'
        write(unit_out) y3
        write(unit_out) indices3
    end do
    close(unit_in)
    close(unit_out)
    print '(a)', 'PASS: randomized EDT feature-index Fortran driver'
end program random_feature_driver
