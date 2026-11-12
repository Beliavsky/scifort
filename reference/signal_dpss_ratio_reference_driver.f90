! SPDX-License-Identifier: MIT
! Streaming ratio-only driver for independent SciPy DPSS comparisons.
program signal_dpss_ratio_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_dpss, only : dpss, dpss_result
    implicit none
    type(dpss_result) :: out
    integer :: m, k, isym, inorm, isolver, iratio, ios
    real(dp) :: nw
    character(16) :: mode, eigensolver, ratio_method

    do
        read(*,*,iostat=ios) m, nw, k, isym, inorm, isolver, iratio
        if (ios /= 0) exit
        mode = 'l2'
        if (inorm == 1) mode = 'approximate'
        if (inorm == 2) mode = 'subsample'
        eigensolver = 'auto'
        if (isolver == 1) eigensolver = 'ql'
        if (isolver == 2) eigensolver = 'partial'
        ratio_method = 'auto'
        if (iratio == 1) ratio_method = 'direct'
        if (iratio == 2) ratio_method = 'fft'
        out = dpss(m,nw,kmax=k,sym=isym /= 0,norm=trim(mode), &
            solver=trim(eigensolver),ratio_method=trim(ratio_method))
        write(*,'(a,1x,i0,1x,i0,1x,i0)') 'CASE',m,k,out%status
        if (out%status == 0) write(*,'(*(es25.16,1x))') out%ratios
    end do
end program signal_dpss_ratio_reference_driver
