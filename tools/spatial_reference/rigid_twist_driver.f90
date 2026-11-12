! SPDX-License-Identifier: MIT
! Independent SE(3) exponential-coordinate reference driver.
program rigid_twist_driver
    use scifort_kinds, only : dp
    use scifort_spatial_rigid_transform, only : rigid_transform_type, rigid_from_exp_coords
    implicit none
    type(rigid_transform_type) :: tf
    real(dp) :: coords(6)
    integer :: n,i,code
    read(*,*) n
    do i=1,n
        read(*,*) coords
        tf = rigid_from_exp_coords(coords,status=code)
        if (code/=0) error stop 'invalid exponential coordinates'
        write(*,'(*(ES25.16E3,1X))') tf%as_matrix()
    end do
end program rigid_twist_driver
