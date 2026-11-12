! Standalone comparison-only driver (not an FPM test target).
program place_poles_methods_driver
    use scifort_kinds, only : dp
    use scifort_signal_place_poles, only : place_poles, pole_placement_result
    implicit none
    real(dp), allocatable :: a(:,:), b(:,:)
    complex(dp), allocatable :: poles(:)
    type(pole_placement_result) :: fitted
    character(len=4) :: method
    integer :: n, m, niter, ios, i, j
    real(dp) :: tolerance
    do
        read(*,*,iostat=ios) n, m, method, niter, tolerance
        if (ios /= 0) exit
        allocate(a(n,n),b(n,m),poles(n))
        read(*,*) a
        read(*,*) b
        read(*,*) poles
        fitted = place_poles(a,b,poles,method=trim(method),maxiter=niter,rtol=tolerance)
        write(*,'(i0,1x,i0,2(1x,es24.16e3))',advance='no') fitted%status, &
            fitted%nb_iter, fitted%rtol, fitted%det_transfer
        if (allocated(fitted%gain_matrix)) then
            do i=1,n
                do j=1,m
                    write(*,'(1x,es24.16e3)',advance='no') fitted%gain_matrix(j,i)
                end do
            end do
        end if
        write(*,*)
        deallocate(a,b,poles)
    end do
end program place_poles_methods_driver
