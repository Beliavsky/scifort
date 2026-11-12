! SPDX-License-Identifier: MIT
! Fortran subprocess driver for independent SciPy comparison.
program rotation_davenport_driver
    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_from_davenport, &
        rotation_from_quat
    implicit none
    type(rotation_type) :: from_angles, from_quat, reconstructed
    real(dp) :: axes(3,3), angles(3), q(4), derived(3), returned(4), qrecon(4)
    integer :: cases, n, i, j, k, code1, code2, code3, code4, lock_int, ios
    logical :: locked
    character(len=16) :: order
    read(input_unit,*,iostat=ios) cases
    if (ios /= 0 .or. cases < 0) error stop 'number of cases'
    do i = 1, cases
        read(input_unit,*,iostat=ios) n, order, ((axes(j,k), k=1,3), j=1,3), angles, q
        if (ios /= 0) error stop 'Davenport input'
        from_angles = rotation_from_davenport(axes(1:n,:),trim(order),angles(1:n),status=code1)
        if (code1 == 0) then
            returned = from_angles%as_quat()
        else
            returned = 0.0_dp
        end if
        derived = 0.0_dp
        qrecon = 0.0_dp
        code2 = 1
        code3 = 1
        code4 = 1
        lock_int = 0
        if (n == 3) then
            from_quat = rotation_from_quat(q,status=code2)
            if (code2 == 0) then
                derived = from_quat%as_davenport(axes,trim(order),gimbal_lock=locked,status=code3)
                if (code3 == 0) then
                    reconstructed = rotation_from_davenport(axes,trim(order),derived,status=code4)
                    if (code4 == 0) then
                        if (.not. reconstructed%valid()) error stop 'bad reconstruction'
                        qrecon = reconstructed%as_quat()
                    end if
                    if (locked) lock_int = 1
                end if
            end if
        end if
        write(output_unit,'(5(i0,1x))') code1, code2, code3, code4, lock_int
        write(output_unit,'(11(es26.17e3,1x))') returned, derived, qrecon
    end do
end program rotation_davenport_driver
