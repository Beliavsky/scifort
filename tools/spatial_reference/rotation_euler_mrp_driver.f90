! SPDX-License-Identifier: MIT
program rotation_euler_mrp_driver
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_from_euler, &
        rotation_from_mrp, rotation_from_quat
    implicit none
    type(rotation_type) :: from_euler, from_mrp, base, back
    integer :: n, i, status, s1, s2, s3, s4
    character(len=3) :: seq
    real(dp) :: a(3), q(4), mrp(3), aa(3)
    logical :: lock
    read(*,*) n
    do i=1,n
        read(*,*) seq, a, q, mrp
        from_euler = rotation_from_euler(seq,a,status=s1)
        base = rotation_from_quat(q,status=s2)
        aa = base%as_euler(seq,gimbal_lock=lock,status=s3)
        back = rotation_from_euler(seq,aa,status=s4)
        from_mrp = rotation_from_mrp(mrp,status=status)
        write(*,'(6(i0,1x))') s1,s2,s3,s4,status,merge(1,0,lock)
        write(*,'(24(es26.17e3,1x))') from_euler%as_quat(),aa,back%as_quat(), &
            from_mrp%as_quat(),base%as_mrp(),base%as_quat()
    end do
end program rotation_euler_mrp_driver
