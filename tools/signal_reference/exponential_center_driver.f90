! SPDX-License-Identifier: MIT
program exponential_center_driver
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : get_window
    implicit none
    real(dp) :: center, tau
    real(dp), allocatable :: y(:)
    integer :: n, st, io, j
    character(12) :: kind
    do
        read(*,*,iostat=io) n, center, tau, kind
        if (io<0) exit
        if (io/=0) error stop 'bad reference input'
        y=get_window(trim(kind),n,parameter=tau,parameter2=center,fftbins=.true.,status=st)
        write(*,'(I0,1X,I0)') st,size(y)
        if (st/=0) cycle
        do j=1,size(y)
            write(*,'(ES24.16E3)') y(j)
        end do
    end do
end program exponential_center_driver
