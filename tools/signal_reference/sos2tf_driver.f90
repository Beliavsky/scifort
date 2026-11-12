! SPDX-License-Identifier: MIT
! Test-only stream driver for SOS to transfer-function conversion.
program sos2tf_driver
    use scifort_kinds, only : dp
    use scifort_signal_design, only : butter_result, sos2tf
    implicit none

    type(butter_result) :: coefficients
    real(dp), allocatable :: sos(:,:)
    integer :: sections, i, ios

    do
        read(*, *, iostat=ios) sections
        if (ios /= 0) exit
        if (sections < 0) exit
        allocate(sos(sections,6))
        do i = 1, sections
            read(*, *) sos(i,:)
        end do
        coefficients = sos2tf(sos)
        write(*, '(i0,1x,i0)') coefficients%status, size(coefficients%b)
        write(*, '(*(es25.16e3,1x))') coefficients%b
        write(*, '(*(es25.16e3,1x))') coefficients%a
        deallocate(sos)
    end do
end program sos2tf_driver
