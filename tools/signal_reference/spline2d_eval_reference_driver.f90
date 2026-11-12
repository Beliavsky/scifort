! SPDX-License-Identifier: MIT
! Text protocol for comparison with scipy.ndimage.map_coordinates.
program spline2d_eval_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_spline2d_eval, only : cspline2d_eval, qspline2d_eval
    implicit none
    integer :: cases, t, m, n, degree, is_complex, count, i, j, status
    real(dp) :: re, im
    real(dp), allocatable :: rows(:), cols(:), input_real(:,:), output_real(:)
    complex(dp), allocatable :: input_complex(:,:), output_complex(:)

    read(*,*) cases
    do t = 1, cases
        read(*,*) m, n, degree, is_complex, count
        allocate(input_real(m,n), input_complex(m,n), rows(count), cols(count))
        do j = 1, n
            do i = 1, m
                read(*,*) re, im
                input_real(i,j) = re
                input_complex(i,j) = cmplx(re, im, dp)
            end do
        end do
        do i = 1, count
            read(*,*) rows(i), cols(i)
        end do
        if (is_complex == 1) then
            if (degree == 3) then
                output_complex = cspline2d_eval(input_complex, rows, cols, status)
            else
                output_complex = qspline2d_eval(input_complex, rows, cols, status)
            end if
        else
            if (degree == 3) then
                output_real = cspline2d_eval(input_real, rows, cols, status)
            else
                output_real = qspline2d_eval(input_real, rows, cols, status)
            end if
            output_complex = cmplx(output_real,0.0_dp,dp)
        end if
        write(*,'(2(i0,1x))') status, size(output_complex)
        if (status == 0) then
            do i = 1, size(output_complex)
                write(*,'(2(es25.17e3,1x))') real(output_complex(i),dp), aimag(output_complex(i))
            end do
        end if
        deallocate(input_real, input_complex, rows, cols)
        if (allocated(output_real)) deallocate(output_real)
        if (allocated(output_complex)) deallocate(output_complex)
    end do
end program spline2d_eval_reference_driver
