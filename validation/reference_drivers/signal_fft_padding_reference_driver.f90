! SPDX-License-Identifier: MIT
! Binary reference driver for FFT padding across mixed-radix and overlap-add.
program signal_fft_padding_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_fftconvolve, only : fftconvolve_nd, oaconvolve, oaconvolve_nd
    implicit none
    real(dp) :: a(101,251), b(29,59)
    real(dp) :: full(129,309), same(101,251), valid(73,193)
    complex(dp) :: ca(101,251), cb(29,59), cy(129,309)
    real(dp) :: a15(101,1,1,1,1,1,1,1,1,1,1,1,1,1,251)
    real(dp) :: b15(29,1,1,1,1,1,1,1,1,1,1,1,1,1,59)
    real(dp) :: y15(129,1,1,1,1,1,1,1,1,1,1,1,1,1,309)
    real(dp) :: x(9001), h(103), oa1(9103)
    real(dp) :: da(513,257), db(11,7), oa2(523,263)
    character(1024) :: source, dest
    integer :: fu, code

    call get_command_argument(1,source)
    call get_command_argument(2,dest)
    open(newunit=fu,file=trim(source),access='stream',form='unformatted',status='old',action='read')
    read(fu) a,b,ca,cb,a15,b15,x,h,da,db
    close(fu)
    call fftconvolve_nd(a,b,full,status=code)
    if (code/=0) error stop 'large fftconvolve status'
    call fftconvolve_nd(a,b,same,mode='same',status=code)
    if (code/=0) error stop 'large same status'
    call fftconvolve_nd(a,b,valid,mode='valid',status=code)
    if (code/=0) error stop 'large valid status'
    call fftconvolve_nd(ca,cb,cy,status=code)
    if (code/=0) error stop 'complex fftconvolve status'
    call fftconvolve_nd(a15,b15,y15,axes=[1,15],status=code)
    if (code/=0) error stop 'rank15 fftconvolve status'
    oa1 = oaconvolve(x,h,status=code)
    if (code/=0) error stop 'blocked rank1 status'
    call oaconvolve_nd(da,db,oa2,status=code)
    if (code/=0) error stop 'blocked rank2 status'
    open(newunit=fu,file=trim(dest),access='stream',form='unformatted',status='replace',action='write')
    write(fu) full,same,valid,cy,y15,oa1,oa2
    close(fu)
end program signal_fft_padding_reference_driver
