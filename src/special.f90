! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_special
    use scifort_incomplete_gamma, only : gammainc, gammaincc, gammainccinv, &
        gammaincinv
    use scifort_kinds, only : dp
    implicit none
    private

    public :: dp
    public :: gammainc
    public :: gammaincc
    public :: gammainccinv
    public :: gammaincinv
end module scifort_special
