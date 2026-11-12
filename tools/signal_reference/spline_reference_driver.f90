! SPDX-License-Identifier: MIT
! Numerical comparison driver for cubic and quadratic spline routines.
program spline_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_spline, only : cspline1d, qspline1d, cspline1d_eval, qspline1d_eval, &
        spline_status_success
    implicit none
    complex(dp), allocatable :: signal(:), coeff(:), values(:)
    real(dp), allocatable :: xx(:)
    real(dp) :: re, im, lam, dx, x0
    integer :: ncase, icase, n, m, degree, iscomplex, status, i

    read(*,*) ncase
    do icase = 1, ncase
        read(*,*) n, m, degree, iscomplex, lam, dx, x0
        allocate(signal(n),xx(m))
        do i = 1, n
            read(*,*) re, im
            signal(i) = cmplx(re, im, dp)
        end do
        do i = 1, m
            read(*,*) xx(i)
        end do
        if (degree == 3) then
            if (iscomplex == 1) then
                coeff = cspline1d(signal, lamb=lam, status=status)
                if (status == spline_status_success) values = cspline1d_eval(coeff, xx, dx, x0, status)
            else
                coeff = cmplx(cspline1d(real(signal,dp), lamb=lam, status=status),0.0_dp,dp)
                if (status == spline_status_success) then
                    values = cmplx(cspline1d_eval(real(coeff,dp),xx,dx,x0,status),0.0_dp,dp)
                end if
            end if
        else
            if (iscomplex == 1) then
                coeff = qspline1d(signal, lamb=lam, status=status)
                if (status == spline_status_success) values = qspline1d_eval(coeff, xx, dx, x0, status)
            else
                coeff = cmplx(qspline1d(real(signal,dp), lamb=lam, status=status),0.0_dp,dp)
                if (status == spline_status_success) then
                    values = cmplx(qspline1d_eval(real(coeff,dp),xx,dx,x0,status),0.0_dp,dp)
                end if
            end if
        end if
        if (status == spline_status_success) then
            write(*,'(3(i0,1x))') status, size(coeff), size(values)
            do i = 1, n
                write(*,'(2(es25.17e3,1x))') real(coeff(i),dp), aimag(coeff(i))
            end do
            do i = 1, m
                write(*,'(2(es25.17e3,1x))') real(values(i),dp), aimag(values(i))
            end do
        else
            write(*,'(3(i0,1x))') status, 0, 0
        end if
        deallocate(signal,xx)
        if (allocated(coeff)) deallocate(coeff)
        if (allocated(values)) deallocate(values)
    end do
end program spline_reference_driver
