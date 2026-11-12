! SPDX-License-Identifier: MIT
! Stream-based interface used by the independent SciPy comparison script.
program linalg_cdf2rdf_driver
    use scifort_kinds, only : dp
    use scifort_linalg_cdf2rdf, only : cdf2rdf
    implicit none

    complex(dp), allocatable :: w(:), v(:, :)
    real(dp), allocatable :: wr(:, :), vr(:, :), rw(:), rv(:, :)
    real(dp) :: xr, xi
    integer :: cases, idx, variant, n, i, j, status

    read(*,*) cases
    do idx = 1, cases
        read(*,*) variant, n
        allocate(w(n),v(n,n),wr(n,n),vr(n,n))
        if (variant == 1) then
            do i = 1, n
                read(*,*) xr, xi
                w(i) = cmplx(xr,xi,dp)
            end do
            do j = 1, n
                do i = 1, n
                    read(*,*) xr, xi
                    v(i,j) = cmplx(xr,xi,dp)
                end do
            end do
            call cdf2rdf(w,v,wr,vr,status)
        else if (variant == 2) then
            allocate(rw(n),rv(n,n))
            do i = 1, n
                read(*,*) rw(i)
            end do
            do j = 1, n
                do i = 1, n
                    read(*,*) rv(i,j)
                end do
            end do
            call cdf2rdf(rw,rv,wr,vr,status)
            deallocate(rw,rv)
        else
            allocate(rw(n))
            do i = 1, n
                read(*,*) rw(i)
            end do
            do j = 1, n
                do i = 1, n
                    read(*,*) xr, xi
                    v(i,j) = cmplx(xr,xi,dp)
                end do
            end do
            call cdf2rdf(rw,v,wr,vr,status)
            deallocate(rw)
        end if
        write(*,'(a,1x,i0,1x,i0)') 'CASE',idx,status
        do j = 1, n
            do i = 1, n
                write(*,'(es26.17e3)') wr(i,j)
            end do
        end do
        do j = 1, n
            do i = 1, n
                write(*,'(es26.17e3)') vr(i,j)
            end do
        end do
        deallocate(w,v,wr,vr)
    end do
end program linalg_cdf2rdf_driver
