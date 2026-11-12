! Generated from scipy.optimize.elementwise, SciPy 1.17.0; do not edit.
! SPDX-License-Identifier: MIT
module optimize_elementwise_rank2_reference
    use scifort_kinds, only : dp
    implicit none
    integer, parameter :: br_root_status_ref(2,3) = &
        reshape([ &
            0, 0, &
            0, 0, &
            0, 0], [2, 3])
    integer, parameter :: br_root_nit_ref(2,3) = &
        reshape([ &
            1, 1, &
            1, 1, &
            1, 1], [2, 3])
    integer, parameter :: br_root_nfev_ref(2,3) = &
        reshape([ &
            4, 4, &
            3, 3, &
            4, 4], [2, 3])
    real(dp), parameter :: br_root_xl_ref(2,3) = &
        reshape([ &
            -3.0000000000000000e+00_dp, -3.0000000000000000e+00_dp, &
            1.0000000000000000e+00_dp, 1.0000000000000000e+00_dp, &
            2.0000000000000000e+00_dp, 2.0000000000000000e+00_dp], [2, 3])
    real(dp), parameter :: br_root_xr_ref(2,3) = &
        reshape([ &
            -2.0000000000000000e+00_dp, -2.0000000000000000e+00_dp, &
            5.5000000000000000e+00_dp, 5.5000000000000000e+00_dp, &
            3.0000000000000000e+00_dp, 3.0000000000000000e+00_dp], [2, 3])
    integer, parameter :: find_root_status_ref(2,3) = &
        reshape([ &
            0, 0, &
            0, 0, &
            -1, 0], [2, 3])
    integer, parameter :: find_root_nit_ref(2,3) = &
        reshape([ &
            6, 6, &
            6, 6, &
            0, 6], [2, 3])
    integer, parameter :: find_root_nfev_ref(2,3) = &
        reshape([ &
            8, 8, &
            8, 8, &
            2, 8], [2, 3])
    real(dp), parameter :: find_root_x_ref(2,3) = &
        reshape([ &
            -2.0000000000000000e+00_dp, 2.0000000000000000e+00_dp, &
            2.0000000000000000e+00_dp, -2.0000000000000000e+00_dp, &
            0.0000000000000000e+00_dp, 2.0000000000000000e+00_dp], [2, 3])
    integer, parameter :: br_min_status_ref(2,3) = &
        reshape([ &
            0, 0, &
            0, 0, &
            0, 0], [2, 3])
    integer, parameter :: br_min_nit_ref(2,3) = &
        reshape([ &
            2, 2, &
            1, 1, &
            0, 0], [2, 3])
    integer, parameter :: br_min_nfev_ref(2,3) = &
        reshape([ &
            5, 5, &
            4, 4, &
            3, 3], [2, 3])
    real(dp), parameter :: br_min_xm_ref(2,3) = &
        reshape([ &
            3.6500000000000004e+00_dp, 3.6500000000000004e+00_dp, &
            8.0000000000000004e-01_dp, 8.0000000000000004e-01_dp, &
            1.5000000000000000e+00_dp, 1.5000000000000000e+00_dp], [2, 3])
    integer, parameter :: find_min_status_ref(2,3) = &
        reshape([ &
            0, 0, &
            0, 0, &
            -1, 0], [2, 3])
    integer, parameter :: find_min_nit_ref(2,3) = &
        reshape([ &
            4, 4, &
            4, 4, &
            0, 4], [2, 3])
    integer, parameter :: find_min_nfev_ref(2,3) = &
        reshape([ &
            7, 7, &
            7, 7, &
            3, 7], [2, 3])
    real(dp), parameter :: find_min_x_ref(2,3) = &
        reshape([ &
            1.5000000000000000e+00_dp, 1.5000000000000000e+00_dp, &
            1.5000000000000000e+00_dp, 1.5000000000000000e+00_dp, &
            0.0000000000000000e+00_dp, 1.5000000000000000e+00_dp], [2, 3])
end module optimize_elementwise_rank2_reference
