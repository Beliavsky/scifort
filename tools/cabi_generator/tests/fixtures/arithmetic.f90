! SPDX-License-Identifier: MIT
module example_arithmetic
    use, intrinsic :: iso_c_binding, only : c_double
    implicit none
    private
    public :: affine, total
contains
    pure elemental function affine(t, shift, gain) result(y)
        implicit none
        real(c_double), intent(in) :: t !! evaluation coordinate
        real(c_double), intent(in), optional :: shift !! additive offset, default zero
        real(c_double), intent(in), optional :: gain !! multiplier, default one
        real(c_double) :: y
        y=t
        if (present(gain)) y=y*gain
        if (present(shift)) y=y+shift
    end function affine

    pure function total(t) result(y)
        implicit none
        real(c_double), intent(in) :: t(:) !! coordinates to sum
        real(c_double) :: y
        y=sum(t)
    end function total
end module example_arithmetic
