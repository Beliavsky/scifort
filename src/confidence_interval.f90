! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
module scifort_confidence_interval
    !! Shared interval value type, independent of the procedures that compute it.
    use scifort_kinds, only : dp
    implicit none
    private
    public :: confidence_interval

    type :: confidence_interval
        real(dp) :: low !! lower confidence limit
        real(dp) :: high !! upper confidence limit
    end type confidence_interval
end module scifort_confidence_interval
