! SPDX-License-Identifier: MIT
! Stream protocol used by compare_polar_frechet_scipy.py.
program driver_polar_frechet
    use scifort_kinds, only : dp
    use scifort_linalg_polar, only : polar
    use scifort_linalg_expm_frechet, only : expm_frechet
    implicit none
    complex(dp), allocatable :: a(:,:), u(:,:), p(:,:), e(:,:), x(:,:), d(:,:)
    real(dp), allocatable :: ra(:,:), ru(:,:), rp(:,:), re(:,:), rx(:,:), rd(:,:)
    real(dp) :: reval, imval, residual, orthogonality
    integer :: count, which, dtype, m, n, leftflag, i, j, k, status, psize
    character(5) :: side

    read(*,*) count
    do k = 1, count
        read(*,*) which, dtype, m, n, leftflag
        allocate(a(m,n))
        do i = 1,m
            do j = 1,n
                read(*,*) reval, imval
                a(i,j) = cmplx(reval,imval,dp)
            end do
        end do
        if (which == 1) then
            side='right'
            psize=n
            if (leftflag == 1) then
                side='left '
                psize=m
            end if
            allocate(u(m,n),p(psize,psize))
            if (dtype == 1) then
                allocate(ra(m,n),ru(m,n),rp(psize,psize))
                ra = real(a,dp)
                call polar(ra,ru,rp,status,side=trim(side),residual=residual,orthogonality=orthogonality)
                u = cmplx(ru,kind=dp)
                p = cmplx(rp,kind=dp)
                deallocate(ra,ru,rp)
            else
                call polar(a,u,p,status,side=trim(side),residual=residual,orthogonality=orthogonality)
            end if
            write(*,'(i0,1x,es26.17,1x,es26.17)') status, residual, orthogonality
            call write_array(u)
            call write_array(p)
            deallocate(u,p)
        else
            allocate(e(m,n),x(m,n),d(m,n))
            do i=1,m
                do j=1,n
                    read(*,*) reval, imval
                    e(i,j) = cmplx(reval,imval,dp)
                end do
            end do
            if (dtype == 1) then
                allocate(ra(m,n),re(m,n),rx(m,n),rd(m,n))
                ra = real(a,dp)
                re = real(e,dp)
                call expm_frechet(ra,re,rx,rd,status)
                x = cmplx(rx,kind=dp)
                d = cmplx(rd,kind=dp)
                deallocate(ra,re,rx,rd)
            else
                call expm_frechet(a,e,x,d,status)
            end if
            write(*,'(i0)') status
            call write_array(x)
            call write_array(d)
            deallocate(e,x,d)
        end if
        deallocate(a)
    end do
contains
    subroutine write_array(b)
        complex(dp), intent(in) :: b(:,:) !! complex matrix
        integer :: row, col
        do row=1,size(b,1)
            do col=1,size(b,2)
                write(*,'(es26.17,1x,es26.17)') real(b(row,col),dp),aimag(b(row,col))
            end do
        end do
    end subroutine write_array
end program driver_polar_frechet
