! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_constants
    use scifort_kinds, only : dp
    implicit none
    private

    real(dp), parameter, public :: scifort_pi = &
        3.141592653589793238462643383279502884197_dp
    real(dp), parameter, public :: scifort_sqrt_two = &
        1.414213562373095048801688724209698078570_dp
    real(dp), parameter, public :: scifort_log_two = &
        0.693147180559945309417232121458176568076_dp
    real(dp), parameter, public :: scifort_log_sqrt_two_pi = &
        0.918938533204672741780329736405617639861_dp
    real(dp), parameter, public :: scifort_inv_sqrt_two_pi = &
        0.398942280401432677939946059934381868476_dp
end module scifort_constants
