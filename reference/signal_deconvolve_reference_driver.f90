! Reference-only driver, kept outside FPM's automatic test directory.
program signal_deconvolve_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_deconvolve, only : deconvolve
    implicit none
    integer :: cases, i, j, n, m, mode, status
    real(dp) :: r, q
    real(dp), allocatable :: x(:), d(:), qr(:), rr(:)
    complex(dp), allocatable :: z(:), c(:), qc(:), rc(:)
    read(*,*) cases
    do i = 1, cases
        read(*,*) n, m, mode
        allocate(x(n),d(m),z(n),c(m))
        do j=1,n
            read(*,*) r,q
            x(j) = r
            z(j) = cmplx(r,q,dp)
        end do
        do j=1,m
            read(*,*) r,q
            d(j) = r
            c(j) = cmplx(r,q,dp)
        end do
        select case(mode)
        case(0)
            call deconvolve(x,d,qr,rr,status)
            write(*,'(3(i0,1x))') status,size(qr),size(rr)
            do j=1,size(qr)
                write(*,'(2(es27.17,1x))') qr(j),0.0_dp
            end do
            do j=1,size(rr)
                write(*,'(2(es27.17,1x))') rr(j),0.0_dp
            end do
        case(1)
            call deconvolve(x,c,qc,rc,status)
            write(*,'(3(i0,1x))') status,size(qc),size(rc)
            do j=1,size(qc)
                write(*,'(2(es27.17,1x))') real(qc(j),dp),aimag(qc(j))
            end do
            do j=1,size(rc)
                write(*,'(2(es27.17,1x))') real(rc(j),dp),aimag(rc(j))
            end do
        case(2)
            call deconvolve(z,d,qc,rc,status)
            write(*,'(3(i0,1x))') status,size(qc),size(rc)
            do j=1,size(qc)
                write(*,'(2(es27.17,1x))') real(qc(j),dp),aimag(qc(j))
            end do
            do j=1,size(rc)
                write(*,'(2(es27.17,1x))') real(rc(j),dp),aimag(rc(j))
            end do
        case(3)
            call deconvolve(z,c,qc,rc,status)
            write(*,'(3(i0,1x))') status,size(qc),size(rc)
            do j=1,size(qc)
                write(*,'(2(es27.17,1x))') real(qc(j),dp),aimag(qc(j))
            end do
            do j=1,size(rc)
                write(*,'(2(es27.17,1x))') real(rc(j),dp),aimag(rc(j))
            end do
        end select
        deallocate(x,d,z,c)
        if(allocated(qr)) deallocate(qr)
        if(allocated(rr)) deallocate(rr)
        if(allocated(qc)) deallocate(qc)
        if(allocated(rc)) deallocate(rc)
    end do
end program signal_deconvolve_reference_driver
