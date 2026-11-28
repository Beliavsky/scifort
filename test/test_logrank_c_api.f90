! SPDX-License-Identifier: MIT
program test_logrank_c_api
    use, intrinsic :: iso_c_binding, only : c_int, c_double
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_empirical_distribution, only : censored_data, logrank, logrank_result
    use scifort_c_api_generated_logrank, only : scifort_logrank_f64
    implicit none
    type(censored_data) :: x, y
    type(logrank_result) :: native
    real(c_double) :: events_x(3), right_x(2), events_y(4), right_y(1), statistic, pvalue
    real(c_double) :: empty(0)
    integer(c_int) :: alternative, status
    character(len=9), parameter :: alternatives(0:2) = ['two-sided','less     ','greater  ']

    events_x = [1.0_c_double,2.0_c_double,2.0_c_double]
    right_x = [2.0_c_double,5.0_c_double]
    events_y = [2.0_c_double,3.0_c_double,4.0_c_double,6.0_c_double]
    right_y = [4.0_c_double]
    x%uncensored = events_x
    x%right = right_x
    y%uncensored = events_y
    y%right = right_y
    do alternative = 0, 2
        native = logrank(x,y,alternative=alternatives(alternative))
        call scifort_logrank_f64(3_c_int,2_c_int,4_c_int,1_c_int,events_x,right_x, &
            events_y,right_y,alternative,statistic,pvalue,status)
        if (status /= 0 .or. statistic /= native%statistic .or. pvalue /= native%pvalue) stop 1
    end do
    call scifort_logrank_f64(3_c_int,0_c_int,4_c_int,0_c_int,events_x,empty, &
        events_y,empty,0_c_int,statistic,pvalue,status)
    native = logrank(events_x,events_y)
    if (status /= 0 .or. statistic /= native%statistic .or. pvalue /= native%pvalue) stop 2
    call scifort_logrank_f64(0_c_int,0_c_int,0_c_int,0_c_int,empty,empty,empty,empty, &
        0_c_int,statistic,pvalue,status)
    if (status /= 0 .or. .not. ieee_is_nan(statistic) .or. .not. ieee_is_nan(pvalue)) stop 3
    call scifort_logrank_f64(0_c_int,0_c_int,0_c_int,0_c_int,empty,empty,empty,empty, &
        3_c_int,statistic,pvalue,status)
    if (status /= 2 .or. .not. ieee_is_nan(statistic) .or. .not. ieee_is_nan(pvalue)) stop 4
    print '(a)', 'Logrank native/C comparison and empty/error checks passed.'
end program test_logrank_c_api
