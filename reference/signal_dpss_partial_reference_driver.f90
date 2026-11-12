! SPDX-License-Identifier: MIT
! Numerical reference driver for selected-eigenpair DPSS tests.
program signal_dpss_partial_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_dpss, only : dpss, dpss_result
    implicit none
    type(dpss_result) :: out
    integer :: m, k, isym, inorm, isolver, ios, j
    real(dp) :: nw
    character(16) :: mode, algorithm
    logical :: symmetric
    do
        read(*,*,iostat=ios) m, nw, k, isym, inorm, isolver
        if (ios /= 0) exit
        symmetric = isym /= 0
        mode = 'l2'
        if (inorm == 1) mode = 'approximate'
        if (inorm == 2) mode = 'subsample'
        algorithm = 'auto'
        if (isolver == 1) algorithm = 'ql'
        if (isolver == 2) algorithm = 'partial'
        out = dpss(m, nw, kmax=k, sym=symmetric, norm=trim(mode), solver=trim(algorithm))
        write(*,'(a,1x,i0,1x,i0,1x,i0)') 'CASE',m,k,out%status
        if (out%status /= 0) cycle
        write(*,'(*(es25.16,1x))') out%ratios
        do j = 1, k
            write(*,'(*(es25.16,1x))') out%tapers(:,j)
        end do
    end do
end program signal_dpss_partial_reference_driver
