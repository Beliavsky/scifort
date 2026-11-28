! SPDX-License-Identifier: MIT
program test_confidence_interval_type
    use scifort_kinds, only : dp
    use scifort_confidence_interval, only : shared_interval => confidence_interval
    use scifort_contingency_meta, only : old_interval => confidence_interval, binomtest, binomtest_result
    use scifort_stats, only : facade_interval => confidence_interval
    use scifort_empirical_distribution, only : quantile_test, quantile_test_result
    implicit none
    type(shared_interval) :: direct
    type(old_interval) :: legacy
    type(facade_interval) :: facade
    type(binomtest_result) :: binomial
    type(quantile_test_result) :: quantile

    legacy = old_interval(0.25_dp,0.75_dp)
    direct = legacy
    facade = direct
    call accept_shared(facade)
    call accept_legacy(direct)
    binomial = binomtest(3,10)
    legacy = binomial%proportion_ci()
    direct = legacy
    if (direct%low < 0.0_dp .or. direct%high > 1.0_dp) stop 1
    quantile = quantile_test([1.0_dp,2.0_dp,3.0_dp,4.0_dp],2.5_dp)
    facade = quantile%confidence_interval()
    direct = facade
    if (direct%low > direct%high) stop 2
    print '(a)', 'Shared confidence interval type and existing imports passed.'
contains
    subroutine accept_shared(ci)
        type(shared_interval), intent(in) :: ci !! interval imported through the public facade
        if (ci%low /= 0.25_dp .or. ci%high /= 0.75_dp) stop 3
    end subroutine accept_shared
    subroutine accept_legacy(ci)
        type(old_interval), intent(in) :: ci !! interval imported directly from the shared module
        if (ci%low /= 0.25_dp .or. ci%high /= 0.75_dp) stop 4
    end subroutine accept_legacy
end program test_confidence_interval_type
