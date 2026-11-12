program poles_driver
    use scifort_kinds, only : dp
    use scifort_signal_place_poles, only : place_poles, pole_placement_result
    implicit none
    real(dp), allocatable :: a(:,:),b(:,:)
    complex(dp), allocatable :: poles(:)
    type(pole_placement_result) :: result
    integer :: n,m,i,j,ios
    do
        read(*,*,iostat=ios) n,m
        if(ios/=0)exit
        allocate(a(n,n),b(n,m),poles(n))
        read(*,*) a
        read(*,*) b
        read(*,*) poles
        result=place_poles(a,b,poles)
        write(*,'(i0)',advance='no') result%status
        if (allocated(result%gain_matrix)) then
            do i=1,n
                do j=1,m
                    write(*,'(1x,es24.16e3)',advance='no') result%gain_matrix(j,i)
                end do
            end do
        end if
        write(*,*)
        deallocate(a,b,poles)
    end do
end program
