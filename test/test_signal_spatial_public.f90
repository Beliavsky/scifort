! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_spatial_public
    use scifort_signal, only : dp, order_filter, sepfir2d, &
        spatial_status_success, spatial_status_invalid_input
    implicit none
    integer :: code
    integer, allocatable :: k(:)
    real(dp), allocatable :: r(:,:)
    complex(dp), allocatable :: c(:,:)

    k = order_filter([2,1,3], [1,1,1], 1, status=code)
    if (code /= spatial_status_success) error stop 'public order status'
    if (any(k /= [1,2,1])) error stop 'public order result'
    r = sepfir2d(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp], [2,2]), &
        [1.0_dp], [1.0_dp], status=code)
    if (code /= spatial_status_success) error stop 'public separable status'
    if (maxval(abs(r-reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp],[2,2]))) > 1.0e-12_dp) &
        error stop 'public separable result'
    c = sepfir2d(r, [cmplx(0.0_dp,1.0_dp,dp)], [1.0_dp], code)
    if (code /= spatial_status_success) error stop 'public complex status'
    if (maxval(abs(c-cmplx(0.0_dp,1.0_dp,dp)*r)) > 1.0e-12_dp) &
        error stop 'public mixed result'
    k = order_filter([2,1,3], [1,1,1], 9, code)
    if (code /= spatial_status_invalid_input) error stop 'public invalid rank'
    print *, 'PASS test_signal_spatial_public'
end program test_signal_spatial_public
