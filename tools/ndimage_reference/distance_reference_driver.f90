! SPDX-License-Identifier: MIT
! Euclidean distance-transform multidimensional reference outputs.
program distance_reference_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : distance_transform_edt_nd
    implicit none
    logical :: m1(7),m2(6,5),m3(4,5,3)
    logical :: m15(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: d1(7),d2(6,5),d3(4,5,3)
    real(dp) :: d15(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: i,j,k,u,s

    do i=1,7
        m1(i)=modulo(3*i+2,5)>1
    end do
    do j=1,5
        do i=1,6
            m2(i,j)=modulo(7*i+3*j+i*j,11)>2
        end do
    end do
    do k=1,3
        do j=1,5
            do i=1,4
                m3(i,j,k)=modulo(i*13+j*5+k*3+i*k,13)>1
            end do
        end do
    end do
    m15=.true.
    m15(1,2,1,1,1,1,1,1,1,1,1,1,1,1,1)=.false.
    open(newunit=u,file='distance_reference_results.bin',form='unformatted', &
         access='stream',status='replace',action='write')
    call distance_transform_edt_nd(m1,d1,s)
    if (s/=0) error stop '1d default'
    write(u) d1
    call distance_transform_edt_nd(m1,d1,s,[2.3_dp])
    if (s/=0) error stop '1d sampling'
    write(u) d1
    call distance_transform_edt_nd(spread(.true.,1,7),d1,s)
    if (s/=0) error stop 'all foreground'
    write(u) d1
    call distance_transform_edt_nd(spread(.false.,1,7),d1,s)
    if (s/=0) error stop 'no foreground'
    write(u) d1
    call distance_transform_edt_nd(m2,d2,s)
    if (s/=0) error stop '2d default'
    write(u) d2
    call distance_transform_edt_nd(m2,d2,s,[1.5_dp,0.6_dp])
    if (s/=0) error stop '2d sampling'
    write(u) d2
    call distance_transform_edt_nd(m3,d3,s)
    if (s/=0) error stop '3d default'
    write(u) d3
    call distance_transform_edt_nd(m3,d3,s,[0.9_dp,1.3_dp,2.2_dp])
    if (s/=0) error stop '3d sampling'
    write(u) d3
    call distance_transform_edt_nd(reshape(spread(.true.,1,30),[6,5]),d2,s, &
        [1.5_dp,0.6_dp])
    if (s/=0) error stop '2d all foreground'
    write(u) d2
    call distance_transform_edt_nd(reshape(spread(.true.,1,60),[4,5,3]),d3,s, &
        [0.9_dp,1.3_dp,2.2_dp])
    if (s/=0) error stop '3d all foreground'
    write(u) d3
    call distance_transform_edt_nd(m15,d15,s)
    if (s/=0) error stop 'rank15 default'
    write(u) d15
    call distance_transform_edt_nd(m15,d15,s, &
        [1.7_dp,0.3_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp, &
         1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp,1.0_dp])
    if (s/=0) error stop 'rank15 sampling'
    write(u) d15
    close(u)
end program distance_reference_driver
