! SPDX-License-Identifier: MIT
module example_arrays
    use, intrinsic :: iso_c_binding, only : c_double
    implicit none
    private
    public :: row_sums
contains
    pure subroutine row_sums(a, gain, y, status)
        !! Scale the row sums of a rectangular matrix.
        implicit none
        real(c_double), intent(in) :: a(:,:) !! matrix in column-major storage
        real(c_double), intent(in) :: gain !! multiplier for each row sum
        real(c_double), intent(out) :: y(:) !! one scaled sum per row
        integer, intent(out) :: status !! 0 success; 1 incompatible shapes
        if (size(y) /= size(a,1)) then
            status=1
            return
        end if
        y=gain*sum(a,dim=2)
        status=0
    end subroutine row_sums
end module example_arrays
