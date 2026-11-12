! SPDX-License-Identifier: MIT
! Compare SciFort pinvh with independently generated SciPy reference cases.
! Usage: driver_linalg_pinvh INPUT_TEXT OUTPUT_TEXT
program driver_linalg_pinvh
    use scifort_kinds, only : dp
    use scifort_linalg_pinvh, only : pinvh
    implicit none

    real(dp), allocatable :: a(:, :), inverse(:, :)
    complex(dp), allocatable :: z(:, :), inverse_z(:, :)
    character(len=1024) :: input_path, output_path
    real(dp) :: ar, ai, atol, rtol
    integer :: uin, uout, ios, case_id, n, matrix_type, use_lower, i, j, status, rank

    call get_command_argument(1, input_path)
    call get_command_argument(2, output_path)
    if (len_trim(input_path) == 0 .or. len_trim(output_path) == 0) &
        error stop 'Provide input and output paths'
    open(newunit=uin, file=trim(input_path), status='old', action='read')
    open(newunit=uout, file=trim(output_path), status='replace', action='write')

    do
        read(uin, *, iostat=ios) case_id, n, matrix_type, use_lower, atol, rtol
        if (ios < 0) exit
        if (ios /= 0) error stop 'Cannot read case header'
        if (n < 0 .or. (matrix_type /= 0 .and. matrix_type /= 1)) error stop 'Invalid case header'
        if (matrix_type == 0) then
            allocate(a(n,n),inverse(n,n))
            do j = 1, n
                do i = 1, n
                    read(uin, *) ar, ai
                    a(i,j) = ar
                end do
            end do
            call pinvh(a, inverse, status, rank, atol=atol, rtol=rtol, lower=use_lower==1)
            write(uout,'(3(i0,1x))') case_id, status, rank
            do j = 1, n
                do i = 1, n
                    write(uout,'(2(es27.18e3,1x))') inverse(i,j), 0.0_dp
                end do
            end do
            deallocate(a,inverse)
        else
            allocate(z(n,n),inverse_z(n,n))
            do j = 1, n
                do i = 1, n
                    read(uin, *) ar, ai
                    z(i,j) = cmplx(ar,ai,dp)
                end do
            end do
            call pinvh(z, inverse_z, status, rank, atol=atol, rtol=rtol, lower=use_lower==1)
            write(uout,'(3(i0,1x))') case_id, status, rank
            do j = 1, n
                do i = 1, n
                    write(uout,'(2(es27.18e3,1x))') real(inverse_z(i,j),dp), aimag(inverse_z(i,j))
                end do
            end do
            deallocate(z,inverse_z)
        end if
    end do
    close(uin)
    close(uout)
end program driver_linalg_pinvh
