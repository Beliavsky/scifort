! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_special
    use scifort_incomplete_beta, only : betainc, betaincc, betainccinv, &
        betaincinv
    use scifort_incomplete_gamma, only : gammainc, gammaincc, gammainccinv, &
        gammaincinv
    use scifort_kinds, only : dp
    use scifort_bessel_i, only : besseli_ratio, log_besseli_scaled
    use scifort_bessel_k, only : log_besselk, log_besselk_scaled, log_besselk_order_derivative, &
        besselk_log_derivative_x
    use scifort_special_elementary, only : betaln, boxcox, boxcox1p, cosm1, &
        digamma, entr, erf, erfc, erfcinv, erfinv, expit, exprel, gammaln, i0e, i1e, &
        inv_boxcox, inv_boxcox1p, &
        kl_div, log_expit, log_ndtr, logit, ndtr, ndtri, psi, rel_entr, &
        xlog1py, xlogy
    use scifort_special_reductions, only : log_softmax, logsumexp, softmax
    use scifort_zeta, only : hurwitz_zeta, hurwitz_zeta_derivative
    implicit none
    private

    public :: besseli_ratio
    public :: betaln
    public :: betainc
    public :: betaincc
    public :: betainccinv
    public :: betaincinv
    public :: boxcox
    public :: boxcox1p
    public :: cosm1
    public :: digamma
    public :: dp
    public :: entr
    public :: erf
    public :: erfc
    public :: erfcinv
    public :: erfinv
    public :: expit
    public :: exprel
    public :: gammaln
    public :: gammainc
    public :: gammaincc
    public :: gammainccinv
    public :: gammaincinv
    public :: i0e
    public :: i1e
    public :: hurwitz_zeta
    public :: hurwitz_zeta_derivative
    public :: inv_boxcox
    public :: inv_boxcox1p
    public :: kl_div
    public :: log_besseli_scaled
    public :: log_besselk
    public :: log_besselk_scaled
    public :: log_besselk_order_derivative
    public :: besselk_log_derivative_x
    public :: log_expit
    public :: log_ndtr
    public :: log_softmax
    public :: logit
    public :: logsumexp
    public :: ndtr
    public :: ndtri
    public :: psi
    public :: rel_entr
    public :: softmax
    public :: xlog1py
    public :: xlogy
end module scifort_special
