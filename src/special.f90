! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_special
    use scifort_incomplete_beta, only : betainc, betaincc, betainccinv, &
        betaincinv
    use scifort_incomplete_gamma, only : gammainc, gammaincc, gammainccinv, &
        gammaincinv
    use scifort_kinds, only : dp
    implicit none
    private

    public :: betainc
    public :: betaincc
    public :: betainccinv
    public :: betaincinv
    public :: dp
    public :: gammainc
    public :: gammaincc
    public :: gammainccinv
    public :: gammaincinv
end module scifort_special
