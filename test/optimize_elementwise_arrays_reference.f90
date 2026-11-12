! Generated from SciPy 1.17.0; do not edit by hand.
! SPDX-License-Identifier: MIT
module test_optimize_elementwise_arrays_reference
    use scifort_kinds, only : dp
    implicit none
    private
    public :: br_root_status_ref, br_root_nit_ref, br_root_nfev_ref, br_root_xl_ref, br_root_xr_ref
    public :: root_status_ref, root_nit_ref, root_nfev_ref, root_x_valid_ref
    public :: br_min_status_ref, br_min_nit_ref, br_min_nfev_ref, br_min_xm_ref
    public :: min_status_ref, min_nit_ref, min_nfev_ref, min_x_valid_ref
    integer, parameter :: br_root_status_ref(3) = [0, 0, 0]
    integer, parameter :: br_root_nit_ref(3) = [1, 1, 1]
    integer, parameter :: br_root_nfev_ref(3) = [4, 3, 4]
    real(dp), parameter :: br_root_xl_ref(3) = [ &
        -3.0000000000000000e+00_dp, &
        1.0000000000000000e+00_dp, &
        2.0000000000000000e+00_dp]
    real(dp), parameter :: br_root_xr_ref(3) = [ &
        -2.0000000000000000e+00_dp, &
        5.5000000000000000e+00_dp, &
        3.0000000000000000e+00_dp]
    integer, parameter :: root_status_ref(4) = [0, 0, -1, 0]
    integer, parameter :: root_nit_ref(4) = [6, 6, 0, 6]
    integer, parameter :: root_nfev_ref(4) = [8, 8, 2, 8]
    real(dp), parameter :: root_x_valid_ref(3) = [ &
        -2.0000000000000000e+00_dp, &
        2.0000000000000000e+00_dp, &
        2.0000000000000000e+00_dp]
    integer, parameter :: br_min_status_ref(3) = [0, 0, 0]
    integer, parameter :: br_min_nit_ref(3) = [2, 1, 0]
    integer, parameter :: br_min_nfev_ref(3) = [5, 4, 3]
    real(dp), parameter :: br_min_xm_ref(3) = [ &
        3.6500000000000004e+00_dp, &
        8.0000000000000004e-01_dp, &
        1.5000000000000000e+00_dp]
    integer, parameter :: min_status_ref(4) = [0, 0, -1, 0]
    integer, parameter :: min_nit_ref(4) = [4, 4, 0, 4]
    integer, parameter :: min_nfev_ref(4) = [7, 7, 3, 7]
    real(dp), parameter :: min_x_valid_ref(3) = [ &
        1.5000000000000000e+00_dp, &
        1.5000000000000000e+00_dp, &
        1.5000000000000000e+00_dp]
end module test_optimize_elementwise_arrays_reference
