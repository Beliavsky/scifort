! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

module scifort_special
    use scifort_incomplete_beta, only : betainc, betaincc, betainccinv, &
        betaincinv
    use scifort_incomplete_gamma, only : gammainc, gammaincc, gammainccinv, &
        gammaincinv
    use scifort_kinds, only : dp
    use scifort_associated_legendre, only : assoc_legendre_p, assoc_legendre_p_all, &
        assoc_legendre_p_derivatives, associated_legendre_status_invalid_argument, &
        associated_legendre_status_success, legendre_p, legendre_p_all, &
        legendre_p_derivatives, lpmv, lqmn, lqn, sph_harm_y, sph_harm_y_all, &
        sph_harm_y_derivatives, sph_legendre_p, sph_legendre_p_all, &
        sph_legendre_p_derivatives
    use scifort_airy, only : ai_zeros, airy, airye, airy_status_invalid_argument, &
        airy_status_numerical_failure, airy_status_success, bi_zeros, itairy
    use scifort_bessel_classic, only : i0, i1, j0, j1, k0, k0e, k1, k1e, y0, y1
    use scifort_bessel_integrals, only : it2i0k0, it2j0y0, iti0k0, itj0y0
    use scifort_bessel_misc, only : bessel_misc_status_invalid_argument, &
        bessel_misc_status_success, besselpoly, lmbda
    use scifort_kelvin, only : bei, bei_zeros, beip, beip_zeros, ber, ber_zeros, berp, &
        berp_zeros, kei, kei_zeros, keip, keip_zeros, kelvin, kelvin_max_zeros, &
        kelvin_status_invalid_argument, kelvin_status_numerical_failure, kelvin_status_success, &
        kelvin_zeros, ker, ker_zeros, kerp, kerp_zeros
    use scifort_mathieu, only : mathieu_a, mathieu_b, mathieu_cem, mathieu_even_coef, &
        mathieu_modcem1, mathieu_modcem2, mathieu_modsem1, mathieu_modsem2, &
        mathieu_odd_coef, mathieu_sem, mathieu_status_invalid_argument, &
        mathieu_status_numerical_failure, mathieu_status_success
    use scifort_carlson_elliptic, only : elliprc, elliprd, elliprf, elliprg, elliprj
    use scifort_ellipsoidal_harmonics, only : ellip_harm, ellip_harm_2, ellip_normal, &
        ellipsoidal_status_invalid_argument, ellipsoidal_status_numerical_failure, &
        ellipsoidal_status_success
    use scifort_complex_digamma, only : digamma_complex, psi_complex
    use scifort_complex_elementary, only : expm1_complex, log1p_complex, &
        log_ndtr_complex, ndtr_complex, xlog1py_complex, xlogy_complex
    use scifort_complex_exponential_integrals, only : exp1_complex, expi_complex
    use scifort_complex_sine_integrals, only : shichi_complex, sici_complex
    use scifort_special_distribution_helpers, only : bdtr, bdtrc, bdtri, bdtrik, bdtrin, &
        btdtria, btdtrib, chdtr, chdtrc, chdtri, chdtriv, fdtr, fdtrc, fdtri, fdtridfd, &
        gdtr, gdtrc, gdtria, gdtrib, gdtrix, nbdtr, nbdtrc, nbdtri, nbdtrik, nbdtrin, &
        nrdtrimn, nrdtrisd, pdtr, pdtrc, pdtri, pdtrik, stdtr, stdtridf, stdtrit
    use scifort_special_noncentral_distribution_helpers, only : chndtr, chndtridf, &
        chndtrinc, chndtrix, ncfdtr, ncfdtri, ncfdtridfd, ncfdtridfn, ncfdtrinc, &
        nctdtr, nctdtridf, nctdtrinc, nctdtrit
    use scifort_error_fresnel, only : dawsn_real => dawsn, erfcx_real => erfcx, &
        erfi_real => erfi, fresnel_real => fresnel
    use scifort_faddeeva, only : dawsn_complex, erf_complex, erfc_complex, &
        erfcx_complex, erfi_complex, fresnel_complex, modfresnelm, modfresnelp, &
        voigt_profile, wofz
    use scifort_erf_zeros, only : erf_zeros, erf_zeros_status_invalid_argument, &
        erf_zeros_status_numerical_failure, erf_zeros_status_success
    use scifort_fresnel_zeros, only : fresnel_zeros, fresnelc_zeros, fresnels_zeros, &
        fresnel_zeros_status_invalid_argument, fresnel_zeros_status_numerical_failure, &
        fresnel_zeros_status_success
    use scifort_exponential_integrals, only : exp1_real => exp1, expi_real => expi, &
        expn, shichi_real => shichi, sici_real => sici
    use scifort_factorials, only : comb, factorial, factorial2, factorialk, perm
    use scifort_gamma_family, only : gamma, gammasgn, loggamma, multigammaln, rgamma
    use scifort_gamma_combinatorics, only : beta, binom, poch, polygamma
    use scifort_integer_sequences, only : bernoulli, euler, &
        integer_sequences_status_invalid_argument, integer_sequences_status_success, stirling2
    use scifort_hypergeometric, only : hyp0f1, hyp1f1, hyp2f1, hyperu
    use scifort_owens_t, only : owens_t
    use scifort_orthogonal_polynomials, only : assoc_laguerre, eval_chebyc, eval_chebys, eval_chebyt, &
        eval_chebyu, eval_genlaguerre, eval_gegenbauer, eval_hermite, eval_hermitenorm, &
        eval_jacobi, eval_laguerre, eval_legendre, eval_sh_chebyt, eval_sh_chebyu, &
        eval_sh_jacobi, eval_sh_legendre
    use scifort_orthogonal_constructors, only : chebyc, chebys, chebyt, chebyu, gegenbauer, &
        genlaguerre, hermite, hermitenorm, jacobi, laguerre, legendre, sh_chebyt, sh_chebyu, &
        sh_jacobi, sh_legendre
    use scifort_orthogonal_roots, only : c_roots, cg_roots, h_roots, he_roots, j_roots, &
        js_roots, l_roots, la_roots, orthogonal_roots_status_invalid_argument, &
        orthogonal_roots_status_numerical_failure, orthogonal_roots_status_success, &
        roots_chebyc, roots_chebys, roots_chebyt, roots_chebyu, roots_genlaguerre, &
        roots_gegenbauer, roots_hermite, roots_hermitenorm, roots_jacobi, roots_laguerre, &
        roots_legendre, roots_sh_chebyt, roots_sh_chebyu, roots_sh_jacobi, roots_sh_legendre, &
        p_roots, ps_roots, s_roots, t_roots, ts_roots, u_roots, us_roots
    use scifort_parabolic_cylinder, only : parabolic_cylinder_status_invalid_argument, &
        parabolic_cylinder_status_success, pbdn_seq, pbdv, pbdv_seq, pbwa, pbvv, pbvv_seq
    use scifort_legendre_elliptic, only : ellipe, ellipeinc, ellipk, ellipkinc, ellipkm1
    use scifort_jacobi_elliptic, only : ellipj
    use scifort_lambert_w, only : lambertw, wrightomega
    use scifort_wright_bessel, only : log_wright_bessel, wright_bessel
    use scifort_bessel_i, only : besseli_ratio, log_besseli_scaled
    use scifort_bessel_jy, only : jn_real => jn, jv_real => jv, jve_real => jve, &
        yv_real => yv, yve_real => yve
    use scifort_bessel_complex, only : iv_complex, ive_complex, ivp_complex, &
        jn_complex, jv_complex, jve_complex, jvp_complex, kv_complex, kve_complex, &
        kvp_complex, yv_complex, yve_complex, yvp_complex
    use scifort_bessel_zeros, only : bessel_zeros_status_invalid_argument, &
        bessel_zeros_status_numerical_failure, bessel_zeros_status_success, &
        jn_zeros, jnjnp_zeros, jnp_zeros, jnyn_zeros, y0_zeros, y1_zeros, y1p_zeros, &
        yn_zeros, ynp_zeros
    use scifort_bessel_derivatives, only : ivp_real => ivp, jvp_real => jvp, &
        kvp_real => kvp, yvp_real => yvp
    use scifort_hankel, only : h1vp, h2vp, hankel1, hankel1e, hankel2, hankel2e
    use scifort_bessel_k, only : log_besselk, log_besselk_scaled, log_besselk_order_derivative, &
        besselk_log_derivative_x
    use scifort_bessel_order, only : iv_real => iv, ive_real => ive, kn, &
        kv_real => kv, kve_real => kve, yn
    use scifort_special_elementary, only : agm, betaln, boxcox, boxcox1p, cbrt, &
        cosdg, cosm1, cotdg, digamma_real => digamma, diric, entr, erf_real => erf, &
        erfc_real => erfc, erfcinv, erfinv, exp10, exp2, expit, expm1_real => expm1, exprel, &
        gammaln, huber, i0e, i1e, inv_boxcox, inv_boxcox1p, kl_div, log1p_real => log1p, &
        log_expit, log_ndtr_real => log_ndtr, logit, ndtr_real => ndtr, ndtri, ndtri_exp, powm1, &
        psi_real => psi, pseudo_huber, radian, rel_entr, round, sinc, sindg, &
        softplus, tandg, xlog1py_real => xlog1py, xlogy_real => xlogy
    use scifort_special_reductions, only : log_softmax, logsumexp, logsumexp_keepdims, &
        logsumexp_signed, logsumexp_signed_keepdims, softmax
    use scifort_special_statistical_legacy, only : kolmogi, kolmogorov, smirnov, smirnovi, tklmbda
    use scifort_spheroidal, only : obl_ang1, obl_ang1_cv, obl_cv, obl_cv_seq, &
        obl_rad1, obl_rad1_cv, obl_rad2, obl_rad2_cv, pro_ang1, pro_ang1_cv, &
        pro_cv, pro_cv_seq, pro_rad1, pro_rad1_cv, pro_rad2, pro_rad2_cv, &
        spheroidal_status_invalid_argument, spheroidal_status_numerical_failure, &
        spheroidal_status_success
    use scifort_spherical_bessel, only : riccati_jn, riccati_yn, &
        spherical_bessel_status_invalid_argument, spherical_bessel_status_success, &
        spherical_in_real => spherical_in, spherical_jn_real => spherical_jn, &
        spherical_kn_real => spherical_kn, spherical_yn_real => spherical_yn
    use scifort_spherical_bessel_complex, only : spherical_in_complex, &
        spherical_jn_complex, spherical_kn_complex, spherical_yn_complex
    use scifort_struve, only : it2struve0, itmodstruve0, itstruve0, modstruve, struve
    use scifort_spence, only : spence
    use scifort_zeta, only : hurwitz_zeta, hurwitz_zeta_derivative, zeta, zetac
    implicit none
    private

    interface dawsn
        module procedure dawsn_real_api
        module procedure dawsn_complex_api
    end interface dawsn

    interface digamma
        module procedure digamma_real_api
        module procedure digamma_complex_api
    end interface digamma

    interface erf
        module procedure erf_real_api
        module procedure erf_complex_api
    end interface erf

    interface erfc
        module procedure erfc_real_api
        module procedure erfc_complex_api
    end interface erfc

    interface erfcx
        module procedure erfcx_real_api
        module procedure erfcx_complex_api
    end interface erfcx

    interface erfi
        module procedure erfi_real_api
        module procedure erfi_complex_api
    end interface erfi

    interface exp1
        module procedure exp1_real_api
        module procedure exp1_complex_api
    end interface exp1

    interface expi
        module procedure expi_real_api
        module procedure expi_complex_api
    end interface expi

    interface expm1
        module procedure expm1_real_api
        module procedure expm1_complex_api
    end interface expm1

    interface log1p
        module procedure log1p_real_api
        module procedure log1p_complex_api
    end interface log1p

    interface log_ndtr
        module procedure log_ndtr_real_api
        module procedure log_ndtr_complex_api
    end interface log_ndtr

    interface ndtr
        module procedure ndtr_real_api
        module procedure ndtr_complex_api
    end interface ndtr

    interface xlog1py
        module procedure xlog1py_real_api
        module procedure xlog1py_complex_api
        module procedure xlog1py_real_complex_api
        module procedure xlog1py_complex_real_api
    end interface xlog1py

    interface xlogy
        module procedure xlogy_real_api
        module procedure xlogy_complex_api
        module procedure xlogy_real_complex_api
        module procedure xlogy_complex_real_api
    end interface xlogy

    interface fresnel
        module procedure fresnel_real_api
        module procedure fresnel_complex_api
    end interface fresnel

    interface shichi
        module procedure shichi_real_api
        module procedure shichi_complex_api
    end interface shichi

    interface sici
        module procedure sici_real_api
        module procedure sici_complex_api
    end interface sici

    interface psi
        module procedure psi_real_api
        module procedure psi_complex_api
    end interface psi

    interface iv
        module procedure iv_real_api
        module procedure iv_complex_api
    end interface iv

    interface ive
        module procedure ive_real_api
        module procedure ive_complex_api
    end interface ive

    interface ivp
        module procedure ivp_real_api
        module procedure ivp_complex_api
    end interface ivp

    interface jn
        module procedure jn_real_api
        module procedure jn_complex_api
    end interface jn

    interface jv
        module procedure jv_real_api
        module procedure jv_complex_api
    end interface jv

    interface jve
        module procedure jve_real_api
        module procedure jve_complex_api
    end interface jve

    interface jvp
        module procedure jvp_real_api
        module procedure jvp_complex_api
    end interface jvp

    interface kv
        module procedure kv_real_api
        module procedure kv_complex_api
    end interface kv

    interface kve
        module procedure kve_real_api
        module procedure kve_complex_api
    end interface kve

    interface kvp
        module procedure kvp_real_api
        module procedure kvp_complex_api
    end interface kvp

    interface yv
        module procedure yv_real_api
        module procedure yv_complex_api
    end interface yv

    interface yve
        module procedure yve_real_api
        module procedure yve_complex_api
    end interface yve

    interface yvp
        module procedure yvp_real_api
        module procedure yvp_complex_api
    end interface yvp

    interface spherical_in
        module procedure spherical_in_real_api
        module procedure spherical_in_complex_api
    end interface spherical_in

    interface spherical_jn
        module procedure spherical_jn_real_api
        module procedure spherical_jn_complex_api
    end interface spherical_jn

    interface spherical_kn
        module procedure spherical_kn_real_api
        module procedure spherical_kn_complex_api
    end interface spherical_kn

    interface spherical_yn
        module procedure spherical_yn_real_api
        module procedure spherical_yn_complex_api
    end interface spherical_yn

    public :: agm
    public :: assoc_laguerre
    public :: ai_zeros
    public :: airy
    public :: airye
    public :: airy_status_invalid_argument
    public :: airy_status_numerical_failure
    public :: airy_status_success
    public :: integer_sequences_status_invalid_argument
    public :: integer_sequences_status_success
    public :: bi_zeros
    public :: binom
    public :: besseli_ratio
    public :: besselpoly
    public :: bessel_misc_status_invalid_argument
    public :: bessel_misc_status_success
    public :: bessel_zeros_status_invalid_argument
    public :: bessel_zeros_status_numerical_failure
    public :: bessel_zeros_status_success
    public :: beta
    public :: bei
    public :: bei_zeros
    public :: beip
    public :: beip_zeros
    public :: ber
    public :: ber_zeros
    public :: bernoulli
    public :: berp
    public :: berp_zeros
    public :: bdtr
    public :: bdtrc
    public :: bdtri
    public :: bdtrik
    public :: bdtrin
    public :: betaln
    public :: betainc
    public :: betaincc
    public :: betainccinv
    public :: betaincinv
    public :: boxcox
    public :: boxcox1p
    public :: cbrt
    public :: btdtria
    public :: btdtrib
    public :: chdtr
    public :: chdtrc
    public :: chdtri
    public :: chdtriv
    public :: chndtr
    public :: chndtridf
    public :: chndtrinc
    public :: chndtrix
    public :: comb
    public :: cosdg
    public :: cotdg
    public :: cosm1
    public :: dawsn
    public :: digamma
    public :: diric
    public :: dp
    public :: entr
    public :: euler
    public :: chebyc
    public :: chebys
    public :: chebyt
    public :: chebyu
    public :: gegenbauer
    public :: genlaguerre
    public :: hermite
    public :: hermitenorm
    public :: jacobi
    public :: laguerre
    public :: legendre
    public :: eval_chebyc
    public :: eval_chebys
    public :: eval_chebyt
    public :: eval_chebyu
    public :: eval_genlaguerre
    public :: eval_gegenbauer
    public :: eval_hermite
    public :: eval_hermitenorm
    public :: eval_jacobi
    public :: eval_laguerre
    public :: eval_legendre
    public :: sh_chebyt
    public :: sh_chebyu
    public :: sh_jacobi
    public :: sh_legendre
    public :: eval_sh_chebyt
    public :: eval_sh_chebyu
    public :: eval_sh_jacobi
    public :: eval_sh_legendre
    public :: orthogonal_roots_status_invalid_argument
    public :: orthogonal_roots_status_numerical_failure
    public :: orthogonal_roots_status_success
    public :: c_roots
    public :: cg_roots
    public :: h_roots
    public :: he_roots
    public :: j_roots
    public :: js_roots
    public :: l_roots
    public :: la_roots
    public :: p_roots
    public :: ps_roots
    public :: s_roots
    public :: t_roots
    public :: ts_roots
    public :: u_roots
    public :: us_roots
    public :: roots_chebyc
    public :: roots_chebys
    public :: roots_chebyt
    public :: roots_chebyu
    public :: roots_genlaguerre
    public :: roots_gegenbauer
    public :: roots_hermite
    public :: roots_hermitenorm
    public :: roots_jacobi
    public :: roots_laguerre
    public :: roots_legendre
    public :: roots_sh_chebyt
    public :: roots_sh_chebyu
    public :: roots_sh_jacobi
    public :: roots_sh_legendre
    public :: exp1
    public :: exp10
    public :: exp2
    public :: expi
    public :: expn
    public :: factorial
    public :: factorial2
    public :: factorialk
    public :: fdtr
    public :: fdtrc
    public :: fdtri
    public :: fdtridfd
    public :: ellip_harm
    public :: ellip_harm_2
    public :: ellip_normal
    public :: ellipsoidal_status_invalid_argument
    public :: ellipsoidal_status_numerical_failure
    public :: ellipsoidal_status_success
    public :: ellipe
    public :: ellipeinc
    public :: ellipk
    public :: ellipkinc
    public :: ellipkm1
    public :: ellipj
    public :: elliprc
    public :: elliprd
    public :: elliprf
    public :: elliprg
    public :: elliprj
    public :: erf
    public :: erfc
    public :: erfcinv
    public :: erfcx
    public :: erfi
    public :: erfinv
    public :: erf_zeros
    public :: erf_zeros_status_invalid_argument
    public :: erf_zeros_status_numerical_failure
    public :: erf_zeros_status_success
    public :: expit
    public :: expm1
    public :: exprel
    public :: fresnel
    public :: fresnel_zeros
    public :: fresnelc_zeros
    public :: fresnels_zeros
    public :: fresnel_zeros_status_invalid_argument
    public :: fresnel_zeros_status_numerical_failure
    public :: fresnel_zeros_status_success
    public :: gammaln
    public :: gammainc
    public :: gammaincc
    public :: gammainccinv
    public :: gammaincinv
    public :: gdtr
    public :: gdtrc
    public :: gdtria
    public :: gdtrib
    public :: gdtrix
    public :: gamma
    public :: gammasgn
    public :: hurwitz_zeta
    public :: hurwitz_zeta_derivative
    public :: hyp0f1
    public :: hyp1f1
    public :: hyp2f1
    public :: hyperu
    public :: h1vp
    public :: h2vp
    public :: hankel1
    public :: hankel1e
    public :: hankel2
    public :: hankel2e
    public :: huber
    public :: i0
    public :: i0e
    public :: i1
    public :: i1e
    public :: inv_boxcox
    public :: inv_boxcox1p
    public :: itairy
    public :: it2i0k0
    public :: it2j0y0
    public :: iti0k0
    public :: itj0y0
    public :: iv
    public :: ivp
    public :: ive
    public :: j0
    public :: j1
    public :: jn
    public :: jn_zeros
    public :: jnjnp_zeros
    public :: jnp_zeros
    public :: jnyn_zeros
    public :: jv
    public :: jvp
    public :: jve
    public :: k0
    public :: kei
    public :: kei_zeros
    public :: keip
    public :: keip_zeros
    public :: kelvin
    public :: kelvin_max_zeros
    public :: kelvin_status_invalid_argument
    public :: kelvin_status_numerical_failure
    public :: kelvin_status_success
    public :: kelvin_zeros
    public :: k0e
    public :: k1
    public :: k1e
    public :: kn
    public :: kv
    public :: kvp
    public :: kve
    public :: kl_div
    public :: lambertw
    public :: ker
    public :: ker_zeros
    public :: kerp
    public :: kerp_zeros
    public :: kolmogi
    public :: kolmogorov
    public :: log_besseli_scaled
    public :: log_besselk
    public :: log_besselk_scaled
    public :: log_besselk_order_derivative
    public :: besselk_log_derivative_x
    public :: log_expit
    public :: log1p
    public :: log_ndtr
    public :: log_softmax
    public :: logit
    public :: lpmv
    public :: logsumexp
    public :: logsumexp_keepdims
    public :: logsumexp_signed
    public :: logsumexp_signed_keepdims
    public :: lmbda
    public :: mathieu_a
    public :: mathieu_b
    public :: mathieu_cem
    public :: mathieu_even_coef
    public :: mathieu_modcem1
    public :: mathieu_modcem2
    public :: mathieu_modsem1
    public :: mathieu_modsem2
    public :: mathieu_odd_coef
    public :: mathieu_sem
    public :: mathieu_status_invalid_argument
    public :: mathieu_status_numerical_failure
    public :: mathieu_status_success
    public :: modfresnelm
    public :: modfresnelp
    public :: log_wright_bessel
    public :: loggamma
    public :: multigammaln
    public :: ndtr
    public :: ndtri
    public :: ndtri_exp
    public :: nrdtrimn
    public :: nrdtrisd
    public :: nbdtr
    public :: nbdtrc
    public :: nbdtri
    public :: nbdtrik
    public :: nbdtrin
    public :: ncfdtr
    public :: ncfdtri
    public :: ncfdtridfd
    public :: ncfdtridfn
    public :: ncfdtrinc
    public :: nctdtr
    public :: nctdtridf
    public :: nctdtrinc
    public :: nctdtrit
    public :: owens_t
    public :: obl_ang1
    public :: obl_ang1_cv
    public :: obl_cv
    public :: obl_cv_seq
    public :: obl_rad1
    public :: obl_rad1_cv
    public :: obl_rad2
    public :: obl_rad2_cv
    public :: perm
    public :: parabolic_cylinder_status_invalid_argument
    public :: parabolic_cylinder_status_success
    public :: pdtr
    public :: pdtrc
    public :: pdtri
    public :: pdtrik
    public :: pbdn_seq
    public :: pbdv
    public :: pbdv_seq
    public :: pbwa
    public :: pbvv
    public :: pbvv_seq
    public :: poch
    public :: pro_ang1
    public :: pro_ang1_cv
    public :: pro_cv
    public :: pro_cv_seq
    public :: pro_rad1
    public :: pro_rad1_cv
    public :: pro_rad2
    public :: pro_rad2_cv
    public :: polygamma
    public :: psi
    public :: powm1
    public :: radian
    public :: pseudo_huber
    public :: rel_entr
    public :: rgamma
    public :: round
    public :: sinc
    public :: sindg
    public :: softplus
    public :: softmax
    public :: smirnov
    public :: smirnovi
    public :: spence
    public :: associated_legendre_status_invalid_argument
    public :: associated_legendre_status_success
    public :: assoc_legendre_p
    public :: assoc_legendre_p_all
    public :: assoc_legendre_p_derivatives
    public :: legendre_p
    public :: legendre_p_all
    public :: legendre_p_derivatives
    public :: lqmn
    public :: lqn
    public :: sph_harm_y
    public :: sph_harm_y_all
    public :: sph_harm_y_derivatives
    public :: sph_legendre_p
    public :: sph_legendre_p_all
    public :: sph_legendre_p_derivatives
    public :: struve
    public :: modstruve
    public :: itstruve0
    public :: it2struve0
    public :: itmodstruve0
    public :: spherical_in
    public :: spherical_jn
    public :: spherical_kn
    public :: spherical_yn
    public :: stdtr
    public :: stdtridf
    public :: stdtrit
    public :: stirling2
    public :: tandg
    public :: tklmbda
    public :: riccati_jn
    public :: riccati_yn
    public :: spherical_bessel_status_invalid_argument
    public :: spherical_bessel_status_success
    public :: spheroidal_status_invalid_argument
    public :: spheroidal_status_numerical_failure
    public :: spheroidal_status_success
    public :: shichi
    public :: sici
    public :: voigt_profile
    public :: wofz
    public :: wright_bessel
    public :: wrightomega
    public :: xlog1py
    public :: xlogy
    public :: y0
    public :: y0_zeros
    public :: y1
    public :: y1_zeros
    public :: y1p_zeros
    public :: yv
    public :: yvp
    public :: yve
    public :: yn
    public :: yn_zeros
    public :: ynp_zeros
    public :: zeta
    public :: zetac

contains

    pure elemental function expm1_real_api(x) result(value)
        real(dp), intent(in) :: x !! real exponent in exp(x)-1
        real(dp) :: value

        value = expm1_real(x)
    end function expm1_real_api

    pure elemental function expm1_complex_api(z) result(value)
        complex(dp), intent(in) :: z !! complex exponent in exp(z)-1
        complex(dp) :: value

        value = expm1_complex(z)
    end function expm1_complex_api

    pure elemental function log1p_real_api(x) result(value)
        real(dp), intent(in) :: x !! real argument in log(1+x)
        real(dp) :: value

        value = log1p_real(x)
    end function log1p_real_api

    pure elemental function log1p_complex_api(z) result(value)
        complex(dp), intent(in) :: z !! complex argument in principal log(1+z)
        complex(dp) :: value

        value = log1p_complex(z)
    end function log1p_complex_api

    pure elemental function log_ndtr_real_api(x) result(value)
        real(dp), intent(in) :: x !! real standard-normal variate
        real(dp) :: value

        value = log_ndtr_real(x)
    end function log_ndtr_real_api

    pure elemental function log_ndtr_complex_api(z) result(value)
        complex(dp), intent(in) :: z !! complex standard-normal variate
        complex(dp) :: value

        value = log_ndtr_complex(z)
    end function log_ndtr_complex_api

    pure elemental function ndtr_real_api(x) result(value)
        real(dp), intent(in) :: x !! real standard-normal variate
        real(dp) :: value

        value = ndtr_real(x)
    end function ndtr_real_api

    pure elemental function ndtr_complex_api(z) result(value)
        complex(dp), intent(in) :: z !! complex standard-normal variate
        complex(dp) :: value

        value = ndtr_complex(z)
    end function ndtr_complex_api

    pure elemental function xlogy_real_api(x, y_arg) result(value)
        real(dp), intent(in) :: x !! real multiplier
        real(dp), intent(in) :: y_arg !! real logarithm argument
        real(dp) :: value

        value = xlogy_real(x, y_arg)
    end function xlogy_real_api

    pure elemental function xlogy_complex_api(x, y_arg) result(value)
        complex(dp), intent(in) :: x !! complex multiplier
        complex(dp), intent(in) :: y_arg !! complex logarithm argument
        complex(dp) :: value

        value = xlogy_complex(x, y_arg)
    end function xlogy_complex_api

    pure elemental function xlogy_real_complex_api(x, y_arg) result(value)
        real(dp), intent(in) :: x !! real multiplier
        complex(dp), intent(in) :: y_arg !! complex logarithm argument
        complex(dp) :: value

        value = xlogy_complex(cmplx(x, 0.0_dp, kind=dp), y_arg)
    end function xlogy_real_complex_api

    pure elemental function xlogy_complex_real_api(x, y_arg) result(value)
        complex(dp), intent(in) :: x !! complex multiplier
        real(dp), intent(in) :: y_arg !! real logarithm argument promoted to complex
        complex(dp) :: value

        value = xlogy_complex(x, cmplx(y_arg, 0.0_dp, kind=dp))
    end function xlogy_complex_real_api

    pure elemental function xlog1py_real_api(x, y_arg) result(value)
        real(dp), intent(in) :: x !! real multiplier
        real(dp), intent(in) :: y_arg !! real argument in log(1+y_arg)
        real(dp) :: value

        value = xlog1py_real(x, y_arg)
    end function xlog1py_real_api

    pure elemental function xlog1py_complex_api(x, y_arg) result(value)
        complex(dp), intent(in) :: x !! complex multiplier
        complex(dp), intent(in) :: y_arg !! complex argument in log(1+y_arg)
        complex(dp) :: value

        value = xlog1py_complex(x, y_arg)
    end function xlog1py_complex_api

    pure elemental function xlog1py_real_complex_api(x, y_arg) result(value)
        real(dp), intent(in) :: x !! real multiplier
        complex(dp), intent(in) :: y_arg !! complex argument in log(1+y_arg)
        complex(dp) :: value

        value = xlog1py_complex(cmplx(x, 0.0_dp, kind=dp), y_arg)
    end function xlog1py_real_complex_api

    pure elemental function xlog1py_complex_real_api(x, y_arg) result(value)
        complex(dp), intent(in) :: x !! complex multiplier
        real(dp), intent(in) :: y_arg !! real argument promoted to complex
        complex(dp) :: value

        value = xlog1py_complex(x, cmplx(y_arg, 0.0_dp, kind=dp))
    end function xlog1py_complex_real_api

    pure elemental function spherical_jn_real_api(n, x, derivative) result(value)
        integer, intent(in) :: n !! nonnegative spherical-Bessel order
        real(dp), intent(in) :: x !! real argument
        logical, intent(in), optional :: derivative !! return d/dx when true; default false
        real(dp) :: value

        value = spherical_jn_real(n, x, derivative)
    end function spherical_jn_real_api

    pure elemental function spherical_jn_complex_api(n, z, derivative) result(value)
        integer, intent(in) :: n !! nonnegative spherical-Bessel order
        complex(dp), intent(in) :: z !! finite complex argument
        logical, intent(in), optional :: derivative !! return d/dz when true; default false
        complex(dp) :: value

        value = spherical_jn_complex(n, z, derivative)
    end function spherical_jn_complex_api

    pure elemental function spherical_yn_real_api(n, x, derivative) result(value)
        integer, intent(in) :: n !! nonnegative spherical-Bessel order
        real(dp), intent(in) :: x !! real argument
        logical, intent(in), optional :: derivative !! return d/dx when true; default false
        real(dp) :: value

        value = spherical_yn_real(n, x, derivative)
    end function spherical_yn_real_api

    pure elemental function spherical_yn_complex_api(n, z, derivative) result(value)
        integer, intent(in) :: n !! nonnegative spherical-Bessel order
        complex(dp), intent(in) :: z !! finite nonzero complex argument
        logical, intent(in), optional :: derivative !! return d/dz when true; default false
        complex(dp) :: value

        value = spherical_yn_complex(n, z, derivative)
    end function spherical_yn_complex_api

    pure elemental function spherical_in_real_api(n, x, derivative) result(value)
        integer, intent(in) :: n !! nonnegative modified spherical-Bessel order
        real(dp), intent(in) :: x !! real argument
        logical, intent(in), optional :: derivative !! return d/dx when true; default false
        real(dp) :: value

        value = spherical_in_real(n, x, derivative)
    end function spherical_in_real_api

    pure elemental function spherical_in_complex_api(n, z, derivative) result(value)
        integer, intent(in) :: n !! nonnegative modified spherical-Bessel order
        complex(dp), intent(in) :: z !! finite complex argument
        logical, intent(in), optional :: derivative !! return d/dz when true; default false
        complex(dp) :: value

        value = spherical_in_complex(n, z, derivative)
    end function spherical_in_complex_api

    pure elemental function spherical_kn_real_api(n, x, derivative) result(value)
        integer, intent(in) :: n !! nonnegative modified spherical-Bessel order
        real(dp), intent(in) :: x !! real argument
        logical, intent(in), optional :: derivative !! return d/dx when true; default false
        real(dp) :: value

        value = spherical_kn_real(n, x, derivative)
    end function spherical_kn_real_api

    pure elemental function spherical_kn_complex_api(n, z, derivative) result(value)
        integer, intent(in) :: n !! nonnegative modified spherical-Bessel order
        complex(dp), intent(in) :: z !! finite nonzero complex argument on the principal branch
        logical, intent(in), optional :: derivative !! return d/dz when true; default false
        complex(dp) :: value

        value = spherical_kn_complex(n, z, derivative)
    end function spherical_kn_complex_api

    pure elemental function jn_real_api(v, x) result(value)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order
        real(dp), intent(in) :: x !! real argument accepted by jn
        real(dp) :: value

        value = jn_real(v, x)
    end function jn_real_api

    pure elemental function jn_complex_api(v, z) result(value)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order
        complex(dp), intent(in) :: z !! finite complex argument
        complex(dp) :: value

        value = jn_complex(v, z)
    end function jn_complex_api

    pure elemental function jv_real_api(v, x) result(value)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order
        real(dp), intent(in) :: x !! real argument accepted by jv
        real(dp) :: value

        value = jv_real(v, x)
    end function jv_real_api

    pure elemental function jv_complex_api(v, z) result(value)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order
        complex(dp), intent(in) :: z !! finite complex argument
        complex(dp) :: value

        value = jv_complex(v, z)
    end function jv_complex_api

    pure elemental function jve_real_api(v, x) result(value)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order
        real(dp), intent(in) :: x !! real argument accepted by jve
        real(dp) :: value

        value = jve_real(v, x)
    end function jve_real_api

    pure elemental function jve_complex_api(v, z) result(value)
        real(dp), intent(in) :: v !! real scaled cylindrical-Bessel-J order
        complex(dp), intent(in) :: z !! finite complex argument
        complex(dp) :: value

        value = jve_complex(v, z)
    end function jve_complex_api

    pure elemental function yv_real_api(v, x) result(value)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order
        real(dp), intent(in) :: x !! real argument accepted by yv
        real(dp) :: value

        value = yv_real(v, x)
    end function yv_real_api

    pure elemental function yv_complex_api(v, z) result(value)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order
        complex(dp), intent(in) :: z !! finite complex argument on the principal branch
        complex(dp) :: value

        value = yv_complex(v, z)
    end function yv_complex_api

    pure elemental function yve_real_api(v, x) result(value)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order
        real(dp), intent(in) :: x !! real argument accepted by yve
        real(dp) :: value

        value = yve_real(v, x)
    end function yve_real_api

    pure elemental function yve_complex_api(v, z) result(value)
        real(dp), intent(in) :: v !! real scaled cylindrical-Bessel-Y order
        complex(dp), intent(in) :: z !! finite nonzero complex argument on the principal branch
        complex(dp) :: value

        value = yve_complex(v, z)
    end function yve_complex_api

    pure elemental function iv_real_api(v, x) result(value)
        real(dp), intent(in) :: v !! real modified-Bessel-I order
        real(dp), intent(in) :: x !! real argument accepted by iv
        real(dp) :: value

        value = iv_real(v, x)
    end function iv_real_api

    pure elemental function iv_complex_api(v, z) result(value)
        real(dp), intent(in) :: v !! real modified-Bessel-I order
        complex(dp), intent(in) :: z !! finite complex argument
        complex(dp) :: value

        value = iv_complex(v, z)
    end function iv_complex_api

    pure elemental function ive_real_api(v, x) result(value)
        real(dp), intent(in) :: v !! real scaled modified-Bessel-I order
        real(dp), intent(in) :: x !! real argument accepted by ive
        real(dp) :: value

        value = ive_real(v, x)
    end function ive_real_api

    pure elemental function ive_complex_api(v, z) result(value)
        real(dp), intent(in) :: v !! real scaled modified-Bessel-I order
        complex(dp), intent(in) :: z !! finite complex argument
        complex(dp) :: value

        value = ive_complex(v, z)
    end function ive_complex_api

    pure elemental function kv_real_api(v, x) result(value)
        real(dp), intent(in) :: v !! real modified-Bessel-K order
        real(dp), intent(in) :: x !! real argument accepted by kv
        real(dp) :: value

        value = kv_real(v, x)
    end function kv_real_api

    pure elemental function kv_complex_api(v, z) result(value)
        real(dp), intent(in) :: v !! real modified-Bessel-K order
        complex(dp), intent(in) :: z !! finite nonzero complex argument on the principal branch
        complex(dp) :: value

        value = kv_complex(v, z)
    end function kv_complex_api

    pure elemental function kve_real_api(v, x) result(value)
        real(dp), intent(in) :: v !! real scaled modified-Bessel-K order
        real(dp), intent(in) :: x !! real argument accepted by kve
        real(dp) :: value

        value = kve_real(v, x)
    end function kve_real_api

    pure elemental function kve_complex_api(v, z) result(value)
        real(dp), intent(in) :: v !! real scaled modified-Bessel-K order
        complex(dp), intent(in) :: z !! finite nonzero complex argument on the principal branch
        complex(dp) :: value

        value = kve_complex(v, z)
    end function kve_complex_api

    pure elemental function jvp_real_api(v, x, n) result(value)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order before differentiation
        real(dp), intent(in) :: x !! real argument accepted by jvp
        integer, intent(in), optional :: n !! derivative order with respect to x, >= 0 (default 1)
        real(dp) :: value

        value = jvp_real(v, x, n)
    end function jvp_real_api

    pure elemental function jvp_complex_api(v, z, n) result(value)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order before differentiation
        complex(dp), intent(in) :: z !! finite complex argument
        integer, intent(in), optional :: n !! derivative order with respect to z, >= 0 (default 1)
        complex(dp) :: value

        value = jvp_complex(v, z, n)
    end function jvp_complex_api

    pure elemental function yvp_real_api(v, x, n) result(value)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order before differentiation
        real(dp), intent(in) :: x !! real argument accepted by yvp
        integer, intent(in), optional :: n !! derivative order with respect to x, >= 0 (default 1)
        real(dp) :: value

        value = yvp_real(v, x, n)
    end function yvp_real_api

    pure elemental function yvp_complex_api(v, z, n) result(value)
        real(dp), intent(in) :: v !! real cylindrical-Bessel order before differentiation
        complex(dp), intent(in) :: z !! finite complex argument on the principal branch
        integer, intent(in), optional :: n !! derivative order with respect to z, >= 0 (default 1)
        complex(dp) :: value

        value = yvp_complex(v, z, n)
    end function yvp_complex_api

    pure elemental function ivp_real_api(v, x, n) result(value)
        real(dp), intent(in) :: v !! real modified-Bessel-I order before differentiation
        real(dp), intent(in) :: x !! real argument accepted by ivp
        integer, intent(in), optional :: n !! derivative order with respect to x, >= 0 (default 1)
        real(dp) :: value

        value = ivp_real(v, x, n)
    end function ivp_real_api

    pure elemental function ivp_complex_api(v, z, n) result(value)
        real(dp), intent(in) :: v !! real modified-Bessel-I order before differentiation
        complex(dp), intent(in) :: z !! finite complex argument
        integer, intent(in), optional :: n !! derivative order with respect to z, >= 0 (default 1)
        complex(dp) :: value

        value = ivp_complex(v, z, n)
    end function ivp_complex_api

    pure elemental function kvp_real_api(v, x, n) result(value)
        real(dp), intent(in) :: v !! real modified-Bessel-K order before differentiation
        real(dp), intent(in) :: x !! real argument accepted by kvp
        integer, intent(in), optional :: n !! derivative order with respect to x, >= 0 (default 1)
        real(dp) :: value

        value = kvp_real(v, x, n)
    end function kvp_real_api

    pure elemental function kvp_complex_api(v, z, n) result(value)
        real(dp), intent(in) :: v !! real modified-Bessel-K order before differentiation
        complex(dp), intent(in) :: z !! finite nonzero complex argument on the principal branch
        integer, intent(in), optional :: n !! derivative order with respect to z, >= 0 (default 1)
        complex(dp) :: value

        value = kvp_complex(v, z, n)
    end function kvp_complex_api

    pure elemental function digamma_real_api(x) result(value)
        real(dp), intent(in) :: x !! real argument of the digamma function
        real(dp) :: value

        value = digamma_real(x)
    end function digamma_real_api

    pure elemental function digamma_complex_api(z) result(value)
        complex(dp), intent(in) :: z !! complex argument of the digamma function
        complex(dp) :: value

        value = digamma_complex(z)
    end function digamma_complex_api

    pure elemental function psi_real_api(x) result(value)
        real(dp), intent(in) :: x !! real argument; alias of digamma
        real(dp) :: value

        value = psi_real(x)
    end function psi_real_api

    pure elemental function psi_complex_api(z) result(value)
        complex(dp), intent(in) :: z !! complex argument; alias of digamma
        complex(dp) :: value

        value = psi_complex(z)
    end function psi_complex_api

    pure elemental function dawsn_real_api(x) result(value)
        real(dp), intent(in) :: x !! real argument of Dawson's integral
        real(dp) :: value

        value = dawsn_real(x)
    end function dawsn_real_api

    pure elemental function dawsn_complex_api(z) result(value)
        complex(dp), intent(in) :: z !! complex argument of Dawson's integral
        complex(dp) :: value

        value = dawsn_complex(z)
    end function dawsn_complex_api

    pure elemental function erf_real_api(x) result(value)
        real(dp), intent(in) :: x !! real argument of the error function
        real(dp) :: value

        value = erf_real(x)
    end function erf_real_api

    pure elemental function erf_complex_api(z) result(value)
        complex(dp), intent(in) :: z !! complex argument of the error function
        complex(dp) :: value

        value = erf_complex(z)
    end function erf_complex_api

    pure elemental function erfc_real_api(x) result(value)
        real(dp), intent(in) :: x !! real argument of the complementary error function
        real(dp) :: value

        value = erfc_real(x)
    end function erfc_real_api

    pure elemental function erfc_complex_api(z) result(value)
        complex(dp), intent(in) :: z !! complex argument of the complementary error function
        complex(dp) :: value

        value = erfc_complex(z)
    end function erfc_complex_api

    pure elemental function erfcx_real_api(x) result(value)
        real(dp), intent(in) :: x !! real argument of exp(x**2)*erfc(x)
        real(dp) :: value

        value = erfcx_real(x)
    end function erfcx_real_api

    pure elemental function erfcx_complex_api(z) result(value)
        complex(dp), intent(in) :: z !! complex argument of exp(z**2)*erfc(z)
        complex(dp) :: value

        value = erfcx_complex(z)
    end function erfcx_complex_api

    pure elemental function erfi_real_api(x) result(value)
        real(dp), intent(in) :: x !! real argument of the imaginary error function
        real(dp) :: value

        value = erfi_real(x)
    end function erfi_real_api

    pure elemental function erfi_complex_api(z) result(value)
        complex(dp), intent(in) :: z !! complex argument of the imaginary error function
        complex(dp) :: value

        value = erfi_complex(z)
    end function erfi_complex_api

    pure elemental function exp1_real_api(x) result(value)
        real(dp), intent(in) :: x !! real argument of the E1 exponential integral
        real(dp) :: value

        value = exp1_real(x)
    end function exp1_real_api

    pure elemental function exp1_complex_api(z) result(value)
        complex(dp), intent(in) :: z !! complex argument of the principal E1 function
        complex(dp) :: value

        value = exp1_complex(z)
    end function exp1_complex_api

    pure elemental function expi_real_api(x) result(value)
        real(dp), intent(in) :: x !! real argument of the principal real Ei function
        real(dp) :: value

        value = expi_real(x)
    end function expi_real_api

    pure elemental function expi_complex_api(z) result(value)
        complex(dp), intent(in) :: z !! complex argument of the principal Ei function
        complex(dp) :: value

        value = expi_complex(z)
    end function expi_complex_api

    pure elemental subroutine shichi_real_api(x, shi, chi)
        real(dp), intent(in) :: x !! real argument of the hyperbolic integrals
        real(dp), intent(out) :: shi !! real hyperbolic sine integral Shi(x)
        real(dp), intent(out) :: chi !! real hyperbolic cosine integral Chi(x)

        call shichi_real(x, shi, chi)
    end subroutine shichi_real_api

    pure elemental subroutine shichi_complex_api(z, shi, chi)
        complex(dp), intent(in) :: z !! complex argument of the hyperbolic integrals
        complex(dp), intent(out) :: shi !! principal hyperbolic sine integral Shi(z)
        complex(dp), intent(out) :: chi !! principal hyperbolic cosine integral Chi(z)

        call shichi_complex(z, shi, chi)
    end subroutine shichi_complex_api

    pure elemental subroutine sici_real_api(x, si, ci)
        real(dp), intent(in) :: x !! real argument of the sine and cosine integrals
        real(dp), intent(out) :: si !! real sine integral Si(x)
        real(dp), intent(out) :: ci !! real cosine-integral branch for x

        call sici_real(x, si, ci)
    end subroutine sici_real_api

    pure elemental subroutine sici_complex_api(z, si, ci)
        complex(dp), intent(in) :: z !! complex argument of the sine and cosine integrals
        complex(dp), intent(out) :: si !! principal sine integral Si(z)
        complex(dp), intent(out) :: ci !! principal cosine integral Ci(z)

        call sici_complex(z, si, ci)
    end subroutine sici_complex_api

    pure elemental subroutine fresnel_real_api(x, s, c)
        real(dp), intent(in) :: x !! real argument of the Fresnel integrals
        real(dp), intent(out) :: s !! Fresnel S integral from zero to x
        real(dp), intent(out) :: c !! Fresnel C integral from zero to x

        call fresnel_real(x, s, c)
    end subroutine fresnel_real_api

    pure elemental subroutine fresnel_complex_api(z, s, c)
        complex(dp), intent(in) :: z !! complex argument of the Fresnel integrals
        complex(dp), intent(out) :: s !! Fresnel S integral from zero to z
        complex(dp), intent(out) :: c !! Fresnel C integral from zero to z

        call fresnel_complex(z, s, c)
    end subroutine fresnel_complex_api

end module scifort_special
