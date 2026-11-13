! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program complex_zeta_driver
    use scifort_kinds, only : dp
    use scifort_zeta, only : zeta
    implicit none
    integer :: n, k
    real(dp) :: x, y
    complex(dp) :: value
    read(*,*) n
    do k = 1, n
        read(*,*) x, y
        value = zeta(cmplx(x,y,dp))
        write(*,'(2(es26.17e3,1x))') real(value,dp),aimag(value)
    end do
end program complex_zeta_driver
