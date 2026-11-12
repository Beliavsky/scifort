! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Binary stream reference driver for SciPy check_COLA/check_NOLA comparisons.
program signal_overlap_checks_reference_driver
    use, intrinsic :: iso_fortran_env, only : int32
    use scifort_kinds, only : dp
    use scifort_signal_istft, only : check_cola, check_nola
    implicit none
    integer(int32) :: dimensions(2), flags(2)
    integer :: input_unit, output_unit, ios
    real(dp) :: tol
    real(dp), allocatable :: w(:)
    character(512) :: input_path, output_path

    if (command_argument_count() /= 2) error stop 'Expected input and output stream file paths'
    call get_command_argument(1, input_path)
    call get_command_argument(2, output_path)
    open(newunit=input_unit, file=trim(input_path), status='old', action='read', &
         access='stream', form='unformatted')
    open(newunit=output_unit, file=trim(output_path), status='replace', action='write', &
         access='stream', form='unformatted')
    do
        read(input_unit, iostat=ios) dimensions
        if (ios /= 0) exit
        read(input_unit) tol
        allocate(w(dimensions(1)))
        read(input_unit) w
        flags(1) = merge(1_int32, 0_int32, check_cola(w,int(dimensions(1)),int(dimensions(2)),tol=tol))
        flags(2) = merge(1_int32, 0_int32, check_nola(w,int(dimensions(1)),int(dimensions(2)),tol=tol))
        write(output_unit) flags
        deallocate(w)
    end do
    close(input_unit)
    close(output_unit)
end program signal_overlap_checks_reference_driver
