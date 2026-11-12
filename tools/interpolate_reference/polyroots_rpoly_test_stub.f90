! SPDX-License-Identifier: MIT
! Test double ONLY: verifies the optional adapter's Fortran ABI when upstream
! polyroots-fortran is unavailable. It does not exercise Jenkins--Traub and
! must never be distributed or described as the actual external solver.
module polyroots_module
    use, intrinsic :: iso_fortran_env, only : real64
    use scifort_interpolate_piecewise, only : piecewise_native_real_roots
    implicit none
    integer, parameter :: polyroots_module_rk=real64
contains
    subroutine rpoly(op,degree,zeror,zeroi,istat)
        integer, intent(in) :: degree !! polynomial order
        real(real64), intent(in) :: op(degree+1) !! descending coefficients
        real(real64), intent(out) :: zeror(degree) !! real parts
        real(real64), intent(out) :: zeroi(degree) !! imaginary parts
        integer, intent(out) :: istat !! zero success for mock solver
        real(real64), allocatable :: r(:)
        integer :: n
        r=piecewise_native_real_roots(op,-100.0_real64,100.0_real64)
        n=min(degree,size(r))
        zeror=1.0_real64
        zeroi=1.0_real64
        zeror(:n)=r(:n)
        zeroi(:n)=0.0_real64
        istat=0
    end subroutine rpoly
end module polyroots_module
