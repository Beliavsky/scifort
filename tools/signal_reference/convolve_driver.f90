program convolve_driver
    use scifort_signal
    implicit none
    integer :: n,m,mode,ios,i,j,status
    complex(dp), allocatable :: a(:),b(:), y(:), z(:)
    real(dp), allocatable :: ar(:),ai(:),br(:),bi(:)
    integer, allocatable :: lags(:)
    character(len=5) :: names(3)
    names=['full ','same ','valid']
    do
        read(*,*,iostat=ios) n,m,mode
        if (ios /= 0) exit
        allocate(ar(n),ai(n),br(m),bi(m),a(n),b(m))
        read(*,*) ar
        read(*,*) ai
        read(*,*) br
        read(*,*) bi
        a=cmplx(ar,ai,dp)
        b=cmplx(br,bi,dp)
        y=convolve(a,b,mode=trim(names(mode)),status=status)
        z=correlate(a,b,mode=trim(names(mode)))
        lags=correlation_lags(n,m,mode=trim(names(mode)))
        print '(i0,1x,i0)', status, size(y)
        do i=1,size(y)
            print '(4(es25.16e3,1x),i0)', real(y(i),dp),aimag(y(i)), &
                real(z(i),dp),aimag(z(i)),lags(i)
        end do
        deallocate(ar,ai,br,bi,a,b,y,z,lags)
    end do
end program
