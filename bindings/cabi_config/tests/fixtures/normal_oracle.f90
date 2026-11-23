! SPDX-License-Identifier: MIT
module generated_normal_oracle
    use, intrinsic :: iso_c_binding, only : c_double, c_int
    use scifort_kinds, only : dp
    use scifort_normal, only : normal_isf, normal_logcdf, normal_logpdf, normal_logsf, normal_sf
    implicit none
    private
    public :: normal_oracle
contains
    function normal_oracle(method, x, loc, scale) result(y) bind(c)
        implicit none
        integer(c_int), value, intent(in) :: method !! 0 isf, 1 logcdf, 2 logpdf, 3 logsf, 4 sf
        real(c_double), value, intent(in) :: x !! evaluation point or probability for isf
        real(c_double), value, intent(in) :: loc !! mean
        real(c_double), value, intent(in) :: scale !! standard deviation, > 0
        real(c_double) :: y
        select case (method)
        case (0)
            y=real(normal_isf(real(x,dp),real(loc,dp),real(scale,dp)),c_double)
        case (1)
            y=real(normal_logcdf(real(x,dp),real(loc,dp),real(scale,dp)),c_double)
        case (2)
            y=real(normal_logpdf(real(x,dp),real(loc,dp),real(scale,dp)),c_double)
        case (3)
            y=real(normal_logsf(real(x,dp),real(loc,dp),real(scale,dp)),c_double)
        case (4)
            y=real(normal_sf(real(x,dp),real(loc,dp),real(scale,dp)),c_double)
        end select
    end function normal_oracle
end module generated_normal_oracle
