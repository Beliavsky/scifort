! SPDX-License-Identifier: MIT
! Stream-based reference driver for SciPy Rotation comparisons.
program rotation_driver
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_from_quat, &
        rotation_from_rotvec, rotation_from_matrix, rotation_slerp
    implicit none
    type(rotation_type) :: a, b, c, d, e, f
    real(dp) :: qa(4), qb(4), v(3), fraction
    integer :: n, i, status
    read(*,*) n
    do i=1,n
        read(*,*) qa, qb, v, fraction
        a = rotation_from_quat(qa,status=status)
        b = rotation_from_quat(qb)
        c = a%compose(b)
        d = rotation_slerp(a,b,fraction)
        e = rotation_from_matrix(a%as_matrix())
        f = rotation_from_rotvec(a%as_rotvec())
        if (.not. a%valid() .or. .not. b%valid() .or. &
            .not. c%valid() .or. .not. d%valid() .or. &
            .not. e%valid() .or. .not. f%valid()) status = 1
        write(*,'(i0)') status
        write(*,'(*(es25.16e3,1x))') a%as_quat(canonical=.true.), &
            reshape(a%as_matrix(),[9]), a%as_rotvec(), &
            a%apply(v),a%apply(v,inverse=.true.),a%magnitude(), &
            e%as_quat(),f%as_quat(),c%as_quat(),d%as_quat()
    end do
end program rotation_driver
