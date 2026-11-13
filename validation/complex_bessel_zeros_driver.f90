! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program complex_bessel_zeros_driver
    use scifort_kinds, only : dp
    use scifort_bessel_zeros, only : y0_zeros, y1_zeros, y1p_zeros
    implicit none
    integer :: family, nt, status, k
    complex(dp), allocatable :: zeros(:), values(:)
    read(*,*) family, nt
    select case (family)
    case (0)
        call y0_zeros(nt,zeros,values,status)
    case (1)
        call y1_zeros(nt,zeros,values,status)
    case (2)
        call y1p_zeros(nt,zeros,values,status)
    case default
        error stop 'unknown family'
    end select
    print '(i0,1x,i0)', status,size(zeros)
    do k = 1, size(zeros)
        print '(4(es26.17e3,1x))', real(zeros(k),dp),aimag(zeros(k)), &
            real(values(k),dp),aimag(values(k))
    end do
end program complex_bessel_zeros_driver
