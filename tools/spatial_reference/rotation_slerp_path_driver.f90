! SPDX-License-Identifier: MIT
program rotation_slerp_path_driver
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_from_quat, &
        rotation_slerp_path
    implicit none
    real(dp), allocatable :: times(:), queries(:)
    type(rotation_type), allocatable :: knots(:), out(:)
    real(dp) :: quat(4)
    integer :: n, nq, i, code
    read(*,*) n, nq
    allocate(times(n),knots(n),queries(nq))
    read(*,*) times
    do i=1,n
        read(*,*) quat
        knots(i) = rotation_from_quat(quat)
    end do
    read(*,*) queries
    call rotation_slerp_path(times,knots,queries,out,code)
    write(*,'(i0)') code
    if (code == 0) then
        do i=1,nq
            write(*,'(*(es25.16e3,1x))') out(i)%as_quat()
        end do
    end if
end program rotation_slerp_path_driver
