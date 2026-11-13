! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program aaa_spectrum_driver
    use scifort_kinds, only : dp
    use scifort_interpolate_aaa_spectrum, only : aaa_arrowhead_eigenvalues
    implicit none
    integer :: m,io
    do
        read(*,*,iostat=io) m
        if (io /= 0) exit
        block
            complex(dp), allocatable :: s(:),u(:),alpha(:),beta(:)
            real(dp) :: sr,si,ur,ui
            integer :: i,status
            allocate(s(m),u(m),alpha(m+1),beta(m+1))
            do i = 1, m
                read(*,*) sr,si,ur,ui
                s(i) = cmplx(sr,si,dp)
                u(i) = cmplx(ur,ui,dp)
            end do
            call aaa_arrowhead_eigenvalues(s,u,alpha,beta,status)
            print '(i0)', status
            if (status /= 0) cycle
            do i = 1, m+1
                print '(4(es26.17e3,1x))', &
                    real(alpha(i),dp),aimag(alpha(i)),real(beta(i),dp),aimag(beta(i))
            end do
        end block
    end do
end program aaa_spectrum_driver
