! SPDX-License-Identifier: MIT
! Emit analog prototypes for comparison to scipy.signal.
program analog_prototype_driver
    use scifort_kinds, only : dp
    use scifort_signal_analog_prototypes, only : analog_prototype_result, &
        buttap, cheb1ap, cheb2ap, besselap
    implicit none
    character(32) :: name, arg, normalization
    type(analog_prototype_result) :: result
    integer :: n, i, ios
    real(dp) :: attenuation

    if (command_argument_count() < 2) error stop 'usage: driver family n [rp_or_rs_or_norm]'
    call get_command_argument(1,name)
    call get_command_argument(2,arg)
    read(arg,*,iostat=ios) n
    if (ios /= 0) error stop 'bad order'
    select case (trim(name))
    case ('buttap')
        result = buttap(n)
    case ('cheb1ap','cheb2ap')
        if (command_argument_count() < 3) error stop 'attenuation required'
        call get_command_argument(3,arg)
        read(arg,*,iostat=ios) attenuation
        if (ios /= 0) error stop 'bad attenuation'
        if (name == 'cheb1ap') then
            result = cheb1ap(n,attenuation)
        else
            result = cheb2ap(n,attenuation)
        end if
    case ('besselap')
        normalization = 'phase'
        if (command_argument_count() >= 3) call get_command_argument(3,normalization)
        result = besselap(n,trim(normalization))
    case default
        error stop 'unknown family'
    end select
    print '(3(i0,1x),es25.16e3)', result%status, size_safe(result%zeros), &
        size_safe(result%poles), result%gain
    if (allocated(result%zeros)) then
        do i = 1, size(result%zeros)
            print '(2(es25.16e3,1x))', real(result%zeros(i),dp), aimag(result%zeros(i))
        end do
    end if
    if (allocated(result%poles)) then
        do i = 1, size(result%poles)
            print '(2(es25.16e3,1x))', real(result%poles(i),dp), aimag(result%poles(i))
        end do
    end if
contains
    integer function size_safe(array) result(nitem)
        complex(dp), allocatable, intent(in) :: array(:) !! optionally allocated roots
        nitem = 0
        if (allocated(array)) nitem = size(array)
    end function size_safe
end program analog_prototype_driver
