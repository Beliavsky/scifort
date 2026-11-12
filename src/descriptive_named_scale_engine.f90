! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Named scale conversion for SciPy-compatible descriptive statistics.

module scifort_descriptive_named_scale_engine
    use scifort_kinds, only : dp
    use scifort_math, only : quiet_nan
    implicit none
    private

    real(dp), parameter :: normal_quartile = 0.6744897501960817_dp

    public :: named_iqr_scale
    public :: named_mad_scale

contains

    pure function named_iqr_scale(scale) result(value)
        character(len=*), intent(in) :: scale !! named IQR scale; currently 'normal'
        real(dp) :: value

        if (is_normal_scale(scale)) then
            value = 2.0_dp * normal_quartile
        else
            value = quiet_nan(0.0_dp)
        end if
    end function named_iqr_scale

    pure function named_mad_scale(scale) result(value)
        character(len=*), intent(in) :: scale !! named MAD scale; currently 'normal'
        real(dp) :: value

        if (is_normal_scale(scale)) then
            value = normal_quartile
        else
            value = quiet_nan(0.0_dp)
        end if
    end function named_mad_scale

    pure logical function is_normal_scale(scale) result(is_normal)
        character(len=*), intent(in) :: scale !! scale name to compare case-insensitively
        character(len=len(scale)) :: lower
        integer :: code
        integer :: i

        lower = scale
        do i = 1, len(scale)
            code = iachar(lower(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                lower(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
        end do
        is_normal = trim(lower) == 'normal'
    end function is_normal_scale

end module scifort_descriptive_named_scale_engine
