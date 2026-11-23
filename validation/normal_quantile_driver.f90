! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program normal_quantile_driver
    use scifort_kinds, only : dp
    use scifort_normal, only : normal_ppf, normal_isf
    implicit none
    real(dp) :: p
    integer :: count, i
    read(*,*) count
    do i=1,count
        read(*,*) p
        write(*,'(2(es26.17e3,1x))') normal_ppf(p),normal_isf(p)
    end do
end program normal_quantile_driver
