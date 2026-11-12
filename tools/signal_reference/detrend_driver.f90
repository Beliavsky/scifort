program detrend_driver
    use scifort_signal
    implicit none
    integer :: n, nb, kind, ios, ierr
    real(dp), allocatable :: x(:), y(:)
    integer, allocatable :: bp(:)
    do
        read(*,*,iostat=ios) n,nb,kind
        if (ios /= 0) exit
        allocate(x(n),bp(nb))
        read(*,*) x
        if (nb>0) read(*,*) bp
        if (kind==1) then
            y = detrend(x, 'constant', bp, ierr)
        else
            y = detrend(x, 'linear', bp, ierr)
        end if
        print '(i0,1x,i0)', ierr, size(y)
        print '(*(es23.15e3,1x))',y
        deallocate(x,bp,y)
    end do
end program
