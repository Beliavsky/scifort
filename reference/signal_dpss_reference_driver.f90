! SPDX-License-Identifier: MIT
! Numerical DPSS driver for test/compare_signal_dpss_scipy.py.
program signal_dpss_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_dpss, only : dpss, dpss_result
    implicit none
    type(dpss_result) :: out
    integer :: m,k,isym,inorm,ios,j
    real(dp) :: nw
    character(16) :: norm
    logical :: symmetric
    do
        read(*,*,iostat=ios) m,nw,k,isym,inorm
        if (ios /= 0) exit
        symmetric = isym /= 0
        norm = 'l2'
        if (inorm == 1) norm = 'approximate'
        if (inorm == 2) norm = 'subsample'
        out = dpss(m,nw,kmax=k,sym=symmetric,norm=trim(norm))
        write(*,'(a,1x,i0,1x,i0,1x,i0)') 'CASE',m,k,out%status
        if (out%status /= 0) cycle
        write(*,'(*(es25.16,1x))') out%ratios
        do j = 1,k
            write(*,'(*(es25.16,1x))') out%tapers(:,j)
        end do
    end do
end program signal_dpss_reference_driver
