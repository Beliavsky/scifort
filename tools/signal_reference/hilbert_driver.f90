! SPDX-License-Identifier: MIT
! Driver for optional Python/SciPy Hilbert cross-checks; not an FPM test.
program hilbert_driver
    use scifort_kinds, only : dp
    use scifort_signal_hilbert, only : hilbert, hilbert2
    implicit none
    integer :: mode, nr, nc, nf, ng, i, j, ios, st
    real(dp), allocatable :: x(:), a(:,:)
    complex(dp), allocatable :: y(:), z(:,:)

    do
        read(*,*,iostat=ios) mode, nr, nc, nf, ng
        if (ios /= 0) exit
        if (mode == 1) then
            allocate(x(nr))
            if (nr > 0) read(*,*) x
            if (nf == -1) then
                y = hilbert(x,status=st)
            else
                y = hilbert(x,n=nf,status=st)
            end if
            print '(2(i0,1x))', st, size(y)
            do i = 1,size(y)
                print '(es27.18,1x,es27.18)', real(y(i),dp), aimag(y(i))
            end do
            deallocate(x)
        else if (mode == 2 .or. mode == 3) then
            allocate(a(nr,nc))
            if (size(a) > 0) read(*,*) a
            if (mode == 3) then
                z = hilbert2(a,n=nf,status=st)
            else if (nf == -1) then
                z = hilbert2(a,status=st)
            else
                z = hilbert2(a,n=[nf,ng],status=st)
            end if
            print '(3(i0,1x))', st, size(z,1), size(z,2)
            do j = 1,size(z,2)
                do i = 1,size(z,1)
                    print '(es27.18,1x,es27.18)', real(z(i,j),dp), aimag(z(i,j))
                end do
            end do
            deallocate(a)
        else
            error stop 'unknown reference mode'
        end if
    end do
end program hilbert_driver
