! SPDX-License-Identifier: MIT
! Driver for independent comparisons with scipy.spatial.transform.RotationSpline.
program rotation_spline_driver
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_from_quat
    use scifort_spatial_rotation_spline, only : rotation_spline_type, rotation_spline
    implicit none
    type(rotation_spline_type) :: spline
    type(rotation_type), allocatable :: keys(:)
    type(rotation_type) :: r
    real(dp), allocatable :: times(:), samples(:)
    real(dp) :: q(4), rate(3), acceleration(3)
    integer :: n, m, i, status, ios
    do
        read(*,*,iostat=ios) n,m
        if (ios /= 0) exit
        allocate(keys(n),times(n),samples(m))
        read(*,*) times
        do i=1,n
            read(*,*) q
            keys(i) = rotation_from_quat(q,status=status)
            if (status /= 0) error stop 'invalid input quaternion'
        end do
        read(*,*) samples
        spline = rotation_spline(times,keys,status=status)
        print '(i0)', status
        if (status == 0) then
            do i=1,m
                r = spline%at(samples(i),status=status)
                rate = spline%angular_rate(samples(i),status=status)
                acceleration = spline%angular_acceleration(samples(i),status=status)
                if (status /= 0) error stop 'invalid sample'
                write(*,'(10(es26.17e3,1x))') r%as_quat(), rate, acceleration
            end do
        end if
        deallocate(keys,times,samples)
    end do
end program rotation_spline_driver
