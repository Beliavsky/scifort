! SPDX-License-Identifier: MIT
! Test-only line-stream driver for independent SciPy comparisons.
program irregular_driver
    use scifort_kinds, only : dp
    use scifort_signal_irregular, only : lombscargle, lombscargle_amplitude, &
        vectorstrength, vector_strength_result, irregular_status_success
    implicit none

    real(dp), allocatable :: x(:), y(:), f(:), w(:), p(:), pn(:)
    complex(dp), allocatable :: amp(:)
    type(vector_strength_result) :: vs
    integer :: n, m, i, ierr, ios, icenter, ifloat, iweight

    do
        read(*,*,iostat=ios) n, m, icenter, ifloat, iweight
        if (ios < 0) exit
        if (ios /= 0 .or. n < 1 .or. m < 1) error stop 'bad driver header'
        allocate(x(n),y(n),f(m),w(n))
        read(*,*) x
        read(*,*) y
        read(*,*) f
        read(*,*) w
        if (iweight /= 0) then
            p = lombscargle(x,y,f,precenter=icenter/=0,weights=w, &
                floating_mean=ifloat/=0,status=ierr)
            pn = lombscargle(x,y,f,precenter=icenter/=0,weights=w, &
                floating_mean=ifloat/=0,normalize='normalize')
            amp = lombscargle_amplitude(x,y,f,precenter=icenter/=0,weights=w, &
                floating_mean=ifloat/=0)
        else
            p = lombscargle(x,y,f,precenter=icenter/=0, &
                floating_mean=ifloat/=0,status=ierr)
            pn = lombscargle(x,y,f,precenter=icenter/=0, &
                floating_mean=ifloat/=0,normalize='normalize')
            amp = lombscargle_amplitude(x,y,f,precenter=icenter/=0, &
                floating_mean=ifloat/=0)
        end if
        if (ierr /= irregular_status_success) error stop 'invalid test case'
        vs = vectorstrength(x, 2.0_dp*acos(-1.0_dp)/f)
        if (vs%status /= irregular_status_success) error stop 'bad vector strength'
        write(*,'(i0)') m
        do i = 1,m
            write(*,'(6(es26.17,1x))') p(i),pn(i), &
                real(amp(i),dp),aimag(amp(i)),vs%strength(i),vs%phase(i)
        end do
        deallocate(x,y,f,w)
    end do
end program irregular_driver
