! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Gaussian-quadrature nodes and weights for classical orthogonal polynomials.
!
! The Jacobi-matrix recurrences and Golub-Welsch weight identity follow NIST
! DLMF 3.5(vi) and the classical recurrence coefficients in DLMF Chapter 18.
! This module is an independent SciFort implementation.

module scifort_orthogonal_roots
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_constants, only : scifort_pi
    use scifort_special_elementary, only : gammaln
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, symmetric_eigen_jacobi
    use scifort_math, only : positive_infinity, quiet_nan
    implicit none
    private

    integer, parameter, public :: orthogonal_roots_status_success = 0
    integer, parameter, public :: orthogonal_roots_status_invalid_argument = 1
    integer, parameter, public :: orthogonal_roots_status_numerical_failure = 2

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

    interface p_roots
        module procedure roots_legendre
    end interface p_roots

    interface j_roots
        module procedure roots_jacobi
    end interface j_roots

    interface cg_roots
        module procedure roots_gegenbauer
    end interface cg_roots

    interface l_roots
        module procedure roots_laguerre
    end interface l_roots

    interface la_roots
        module procedure roots_genlaguerre
    end interface la_roots

    interface h_roots
        module procedure roots_hermite
    end interface h_roots

    interface he_roots
        module procedure roots_hermitenorm
    end interface he_roots

    interface t_roots
        module procedure roots_chebyt
    end interface t_roots

    interface u_roots
        module procedure roots_chebyu
    end interface u_roots

    interface c_roots
        module procedure roots_chebyc
    end interface c_roots

    interface s_roots
        module procedure roots_chebys
    end interface s_roots

    interface ps_roots
        module procedure roots_sh_legendre
    end interface ps_roots

    interface js_roots
        module procedure roots_sh_jacobi
    end interface js_roots

    interface ts_roots
        module procedure roots_sh_chebyt
    end interface ts_roots

    interface us_roots
        module procedure roots_sh_chebyu
    end interface us_roots

contains

    pure subroutine roots_legendre(n, nodes, weights, mu, status)
        integer, intent(in) :: n !! quadrature order, >= 1
        real(dp), allocatable, intent(out) :: nodes(:) !! ascending Gauss-Legendre nodes on [-1,1]
        real(dp), allocatable, intent(out) :: weights(:) !! positive quadrature weights
        real(dp), intent(out), optional :: mu !! exact total weight, 2
        integer, intent(out), optional :: status !! zero on success; nonzero otherwise

        real(dp) :: total
        integer :: info

        call roots_jacobi(n, 0.0_dp, 0.0_dp, nodes, weights, total, info)
        if (present(mu)) mu = total
        if (present(status)) status = info
    end subroutine roots_legendre

    pure subroutine roots_jacobi(n, alpha, beta, nodes, weights, mu, status)
        integer, intent(in) :: n !! quadrature order, >= 1
        real(dp), intent(in) :: alpha !! exponent of (1-x), > -1
        real(dp), intent(in) :: beta !! exponent of (1+x), > -1
        real(dp), allocatable, intent(out) :: nodes(:) !! ascending Gauss-Jacobi nodes on [-1,1]
        real(dp), allocatable, intent(out) :: weights(:) !! positive quadrature weights
        real(dp), intent(out), optional :: mu !! integral of the Jacobi weight on [-1,1]
        integer, intent(out), optional :: status !! zero on success; nonzero otherwise

        real(dp), allocatable :: diagonal(:)
        real(dp), allocatable :: offdiagonal(:)
        real(dp) :: ab
        real(dp) :: beta_k
        real(dp) :: denominator
        real(dp) :: log_mu
        real(dp) :: ratio
        real(dp) :: t
        real(dp) :: total
        integer :: info
        integer :: k

        if (n < 1 .or. .not. ieee_is_finite(alpha) .or. .not. ieee_is_finite(beta) .or. &
                alpha <= -1.0_dp .or. beta <= -1.0_dp) then
            call invalid_rule(nodes, weights, mu, status)
            return
        end if

        allocate(diagonal(n), offdiagonal(max(0, n - 1)))
        ab = alpha + beta
        diagonal(1) = (beta - alpha) / (ab + 2.0_dp)
        do k = 1, n - 1
            t = 2.0_dp * real(k, dp) + ab
            diagonal(k + 1) = (beta * beta - alpha * alpha) / (t * (t + 2.0_dp))

            denominator = t - 1.0_dp
            if (real(k, dp) + ab == 0.0_dp .and. denominator == 0.0_dp) then
                ratio = 1.0_dp
            else
                ratio = (real(k, dp) + ab) / denominator
            end if
            beta_k = 4.0_dp * real(k, dp) * (real(k, dp) + alpha) * &
                (real(k, dp) + beta) * ratio / (t * t * (t + 1.0_dp))
            if (.not. ieee_is_finite(beta_k) .or. beta_k <= 0.0_dp) then
                call invalid_rule(nodes, weights, mu, status)
                return
            end if
            offdiagonal(k) = sqrt(beta_k)
        end do

        log_mu = (ab + 1.0_dp) * log(2.0_dp) + gammaln(alpha + 1.0_dp) + &
            gammaln(beta + 1.0_dp) - gammaln(ab + 2.0_dp)
        total = exp_guarded(log_mu)
        call golub_welsch(diagonal, offdiagonal, total, nodes, weights, info)
        if (info == orthogonal_roots_status_success .and. alpha == beta) then
            call symmetrize_rule(nodes, weights)
        end if
        if (present(mu)) mu = total
        if (present(status)) status = info
    end subroutine roots_jacobi

    pure subroutine roots_gegenbauer(n, alpha, nodes, weights, mu, status)
        integer, intent(in) :: n !! quadrature order, >= 1
        real(dp), intent(in) :: alpha !! Gegenbauer parameter, > -1/2
        real(dp), allocatable, intent(out) :: nodes(:) !! ascending nodes on [-1,1]
        real(dp), allocatable, intent(out) :: weights(:) !! weights for (1-x^2)^(alpha-1/2)
        real(dp), intent(out), optional :: mu !! integral of the Gegenbauer weight
        integer, intent(out), optional :: status !! zero on success; nonzero otherwise

        real(dp) :: total
        integer :: info

        if (.not. ieee_is_finite(alpha) .or. alpha <= -0.5_dp) then
            call invalid_rule(nodes, weights, mu, status)
            return
        end if
        call roots_jacobi(n, alpha - 0.5_dp, alpha - 0.5_dp, &
            nodes, weights, total, info)
        if (present(mu)) mu = total
        if (present(status)) status = info
    end subroutine roots_gegenbauer

    pure subroutine roots_laguerre(n, nodes, weights, mu, status)
        integer, intent(in) :: n !! quadrature order, >= 1
        real(dp), allocatable, intent(out) :: nodes(:) !! ascending Gauss-Laguerre nodes on [0,infinity)
        real(dp), allocatable, intent(out) :: weights(:) !! weights for exp(-x)
        real(dp), intent(out), optional :: mu !! exact total weight, 1
        integer, intent(out), optional :: status !! zero on success; nonzero otherwise

        real(dp) :: total
        integer :: info

        call roots_genlaguerre(n, 0.0_dp, nodes, weights, total, info)
        if (present(mu)) mu = total
        if (present(status)) status = info
    end subroutine roots_laguerre

    pure subroutine roots_genlaguerre(n, alpha, nodes, weights, mu, status)
        integer, intent(in) :: n !! quadrature order, >= 1
        real(dp), intent(in) :: alpha !! generalized-Laguerre exponent, > -1
        real(dp), allocatable, intent(out) :: nodes(:) !! ascending nodes on [0,infinity)
        real(dp), allocatable, intent(out) :: weights(:) !! weights for x^alpha*exp(-x)
        real(dp), intent(out), optional :: mu !! integral Gamma(alpha+1) of the weight
        integer, intent(out), optional :: status !! zero on success; nonzero otherwise

        real(dp), allocatable :: diagonal(:)
        real(dp), allocatable :: offdiagonal(:)
        real(dp) :: total
        integer :: info
        integer :: k

        if (n < 1 .or. .not. ieee_is_finite(alpha) .or. alpha <= -1.0_dp) then
            call invalid_rule(nodes, weights, mu, status)
            return
        end if

        allocate(diagonal(n), offdiagonal(max(0, n - 1)))
        do k = 0, n - 1
            diagonal(k + 1) = 2.0_dp * real(k, dp) + alpha + 1.0_dp
        end do
        do k = 1, n - 1
            offdiagonal(k) = sqrt(real(k, dp) * (real(k, dp) + alpha))
        end do
        total = exp_guarded(gammaln(alpha + 1.0_dp))
        call golub_welsch(diagonal, offdiagonal, total, nodes, weights, info)
        if (present(mu)) mu = total
        if (present(status)) status = info
    end subroutine roots_genlaguerre

    pure subroutine roots_hermite(n, nodes, weights, mu, status)
        integer, intent(in) :: n !! quadrature order, >= 1
        real(dp), allocatable, intent(out) :: nodes(:) !! ascending physicists' Hermite nodes
        real(dp), allocatable, intent(out) :: weights(:) !! weights for exp(-x^2)
        real(dp), intent(out), optional :: mu !! exact total weight, sqrt(pi)
        integer, intent(out), optional :: status !! zero on success; nonzero otherwise

        real(dp), allocatable :: diagonal(:)
        real(dp), allocatable :: offdiagonal(:)
        real(dp) :: total
        integer :: info
        integer :: k

        if (n < 1) then
            call invalid_rule(nodes, weights, mu, status)
            return
        end if
        allocate(diagonal(n), offdiagonal(max(0, n - 1)))
        diagonal = 0.0_dp
        do k = 1, n - 1
            offdiagonal(k) = sqrt(0.5_dp * real(k, dp))
        end do
        total = sqrt(scifort_pi)
        call golub_welsch(diagonal, offdiagonal, total, nodes, weights, info)
        if (info == orthogonal_roots_status_success) call symmetrize_rule(nodes, weights)
        if (present(mu)) mu = total
        if (present(status)) status = info
    end subroutine roots_hermite

    pure subroutine roots_hermitenorm(n, nodes, weights, mu, status)
        integer, intent(in) :: n !! quadrature order, >= 1
        real(dp), allocatable, intent(out) :: nodes(:) !! ascending probabilists' Hermite nodes
        real(dp), allocatable, intent(out) :: weights(:) !! weights for exp(-x^2/2)
        real(dp), intent(out), optional :: mu !! exact total weight, sqrt(2*pi)
        integer, intent(out), optional :: status !! zero on success; nonzero otherwise

        real(dp), allocatable :: diagonal(:)
        real(dp), allocatable :: offdiagonal(:)
        real(dp) :: total
        integer :: info
        integer :: k

        if (n < 1) then
            call invalid_rule(nodes, weights, mu, status)
            return
        end if
        allocate(diagonal(n), offdiagonal(max(0, n - 1)))
        diagonal = 0.0_dp
        do k = 1, n - 1
            offdiagonal(k) = sqrt(real(k, dp))
        end do
        total = sqrt(2.0_dp * scifort_pi)
        call golub_welsch(diagonal, offdiagonal, total, nodes, weights, info)
        if (info == orthogonal_roots_status_success) call symmetrize_rule(nodes, weights)
        if (present(mu)) mu = total
        if (present(status)) status = info
    end subroutine roots_hermitenorm

    pure subroutine roots_chebyt(n, nodes, weights, mu, status)
        integer, intent(in) :: n !! quadrature order, >= 1
        real(dp), allocatable, intent(out) :: nodes(:) !! ascending first-kind Chebyshev nodes on [-1,1]
        real(dp), allocatable, intent(out) :: weights(:) !! weights for 1/sqrt(1-x^2)
        real(dp), intent(out), optional :: mu !! exact total weight, pi
        integer, intent(out), optional :: status !! zero on success; nonzero otherwise

        integer :: i
        real(dp) :: theta

        if (n < 1) then
            call invalid_rule(nodes, weights, mu, status)
            return
        end if
        allocate(nodes(n), weights(n))
        do i = 1, n
            theta = real(2 * (n - i) + 1, dp) * scifort_pi / real(2 * n, dp)
            nodes(i) = cos(theta)
        end do
        weights = scifort_pi / real(n, dp)
        call symmetrize_rule(nodes, weights)
        if (present(mu)) mu = scifort_pi
        if (present(status)) status = orthogonal_roots_status_success
    end subroutine roots_chebyt

    pure subroutine roots_chebyu(n, nodes, weights, mu, status)
        integer, intent(in) :: n !! quadrature order, >= 1
        real(dp), allocatable, intent(out) :: nodes(:) !! ascending second-kind Chebyshev nodes on [-1,1]
        real(dp), allocatable, intent(out) :: weights(:) !! weights for sqrt(1-x^2)
        real(dp), intent(out), optional :: mu !! exact total weight, pi/2
        integer, intent(out), optional :: status !! zero on success; nonzero otherwise

        integer :: i
        real(dp) :: theta

        if (n < 1) then
            call invalid_rule(nodes, weights, mu, status)
            return
        end if
        allocate(nodes(n), weights(n))
        do i = 1, n
            theta = real(n - i + 1, dp) * scifort_pi / real(n + 1, dp)
            nodes(i) = cos(theta)
            weights(i) = scifort_pi * sin(theta)**2 / real(n + 1, dp)
        end do
        call symmetrize_rule(nodes, weights)
        if (present(mu)) mu = 0.5_dp * scifort_pi
        if (present(status)) status = orthogonal_roots_status_success
    end subroutine roots_chebyu

    pure subroutine roots_chebyc(n, nodes, weights, mu, status)
        integer, intent(in) :: n !! quadrature order, >= 1
        real(dp), allocatable, intent(out) :: nodes(:) !! ascending C_n nodes on [-2,2]
        real(dp), allocatable, intent(out) :: weights(:) !! weights for 1/sqrt(1-(x/2)^2)
        real(dp), intent(out), optional :: mu !! exact total weight, 2*pi
        integer, intent(out), optional :: status !! zero on success; nonzero otherwise

        real(dp) :: total
        integer :: info

        call roots_chebyt(n, nodes, weights, total, info)
        if (info == orthogonal_roots_status_success) then
            nodes = 2.0_dp * nodes
            weights = 2.0_dp * weights
            total = 2.0_dp * total
        end if
        if (present(mu)) mu = total
        if (present(status)) status = info
    end subroutine roots_chebyc

    pure subroutine roots_chebys(n, nodes, weights, mu, status)
        integer, intent(in) :: n !! quadrature order, >= 1
        real(dp), allocatable, intent(out) :: nodes(:) !! ascending S_n nodes on [-2,2]
        real(dp), allocatable, intent(out) :: weights(:) !! weights for sqrt(1-(x/2)^2)
        real(dp), intent(out), optional :: mu !! exact total weight, pi
        integer, intent(out), optional :: status !! zero on success; nonzero otherwise

        real(dp) :: total
        integer :: info

        call roots_chebyu(n, nodes, weights, total, info)
        if (info == orthogonal_roots_status_success) then
            nodes = 2.0_dp * nodes
            weights = 2.0_dp * weights
            total = 2.0_dp * total
        end if
        if (present(mu)) mu = total
        if (present(status)) status = info
    end subroutine roots_chebys

    pure subroutine roots_sh_legendre(n, nodes, weights, mu, status)
        integer, intent(in) :: n !! quadrature order, >= 1
        real(dp), allocatable, intent(out) :: nodes(:) !! ascending shifted-Legendre nodes on [0,1]
        real(dp), allocatable, intent(out) :: weights(:) !! weights for unit weight on [0,1]
        real(dp), intent(out), optional :: mu !! exact total weight, 1
        integer, intent(out), optional :: status !! zero on success; nonzero otherwise

        real(dp) :: total
        integer :: info

        call roots_legendre(n, nodes, weights, total, info)
        if (info == orthogonal_roots_status_success) then
            nodes = 0.5_dp * (nodes + 1.0_dp)
            weights = 0.5_dp * weights
            total = 0.5_dp * total
        end if
        if (present(mu)) mu = total
        if (present(status)) status = info
    end subroutine roots_sh_legendre

    pure subroutine roots_sh_jacobi(n, p, q, nodes, weights, mu, status)
        integer, intent(in) :: n !! quadrature order, >= 1
        real(dp), intent(in) :: p !! shifted-Jacobi parameter, p-q > -1
        real(dp), intent(in) :: q !! shifted-Jacobi parameter, q > 0
        real(dp), allocatable, intent(out) :: nodes(:) !! ascending shifted-Jacobi nodes on [0,1]
        real(dp), allocatable, intent(out) :: weights(:) !! weights for (1-x)^(p-q)*x^(q-1)
        real(dp), intent(out), optional :: mu !! integral of the shifted-Jacobi weight
        integer, intent(out), optional :: status !! zero on success; nonzero otherwise

        real(dp) :: scale
        real(dp) :: total
        integer :: info

        if (.not. ieee_is_finite(p) .or. .not. ieee_is_finite(q) .or. &
                q <= 0.0_dp .or. p - q <= -1.0_dp) then
            call invalid_rule(nodes, weights, mu, status)
            return
        end if
        call roots_jacobi(n, p - q, q - 1.0_dp, nodes, weights, total, info)
        if (info == orthogonal_roots_status_success) then
            scale = exp_guarded(-p * log(2.0_dp))
            nodes = 0.5_dp * (nodes + 1.0_dp)
            weights = scale * weights
            total = scale * total
        end if
        if (present(mu)) mu = total
        if (present(status)) status = info
    end subroutine roots_sh_jacobi

    pure subroutine roots_sh_chebyt(n, nodes, weights, mu, status)
        integer, intent(in) :: n !! quadrature order, >= 1
        real(dp), allocatable, intent(out) :: nodes(:) !! ascending shifted first-kind nodes on [0,1]
        real(dp), allocatable, intent(out) :: weights(:) !! weights for 1/sqrt(x-x^2)
        real(dp), intent(out), optional :: mu !! exact total weight, pi
        integer, intent(out), optional :: status !! zero on success; nonzero otherwise

        real(dp) :: total
        integer :: info

        call roots_chebyt(n, nodes, weights, total, info)
        if (info == orthogonal_roots_status_success) nodes = 0.5_dp * (nodes + 1.0_dp)
        if (present(mu)) mu = total
        if (present(status)) status = info
    end subroutine roots_sh_chebyt

    pure subroutine roots_sh_chebyu(n, nodes, weights, mu, status)
        integer, intent(in) :: n !! quadrature order, >= 1
        real(dp), allocatable, intent(out) :: nodes(:) !! ascending shifted second-kind nodes on [0,1]
        real(dp), allocatable, intent(out) :: weights(:) !! weights for sqrt(x-x^2)
        real(dp), intent(out), optional :: mu !! exact total weight, pi/8
        integer, intent(out), optional :: status !! zero on success; nonzero otherwise

        real(dp) :: total
        integer :: info

        call roots_chebyu(n, nodes, weights, total, info)
        if (info == orthogonal_roots_status_success) then
            nodes = 0.5_dp * (nodes + 1.0_dp)
            weights = 0.25_dp * weights
            total = 0.25_dp * total
        end if
        if (present(mu)) mu = total
        if (present(status)) status = info
    end subroutine roots_sh_chebyu

    pure subroutine golub_welsch(diagonal, offdiagonal, total, nodes, weights, status)
        real(dp), intent(in) :: diagonal(:) !! Jacobi-matrix diagonal coefficients
        real(dp), intent(in) :: offdiagonal(:) !! positive Jacobi-matrix subdiagonal coefficients
        real(dp), intent(in) :: total !! zeroth moment of the quadrature weight
        real(dp), allocatable, intent(out) :: nodes(:) !! ascending eigenvalues of the Jacobi matrix
        real(dp), allocatable, intent(out) :: weights(:) !! Gaussian quadrature weights
        integer, intent(out) :: status !! zero on success; nonzero on eigensolver failure

        real(dp), allocatable :: eigenvectors(:, :)
        real(dp), allocatable :: matrix(:, :)
        integer :: info
        integer :: k
        integer :: n

        n = size(diagonal)
        allocate(nodes(n), weights(n), eigenvectors(n, n), matrix(n, n))
        nodes = quiet_nan(0.0_dp)
        weights = quiet_nan(0.0_dp)
        matrix = 0.0_dp
        do k = 1, n
            matrix(k, k) = diagonal(k)
        end do
        do k = 1, n - 1
            matrix(k + 1, k) = offdiagonal(k)
        end do
        call symmetric_eigen_jacobi(matrix, nodes, eigenvectors, info)
        if (info /= linalg_status_success) then
            status = orthogonal_roots_status_numerical_failure
            return
        end if
        weights = total * eigenvectors(1, :)**2
        if (ieee_is_finite(total) .and. sum(weights) > 0.0_dp .and. &
                ieee_is_finite(sum(weights))) then
            weights = weights * (total / sum(weights))
        end if
        status = orthogonal_roots_status_success
    end subroutine golub_welsch

    pure subroutine symmetrize_rule(nodes, weights)
        real(dp), intent(inout) :: nodes(:) !! ascending symmetric nodes to average in opposite pairs
        real(dp), intent(inout) :: weights(:) !! corresponding symmetric weights to average in pairs

        integer :: i
        integer :: j
        real(dp) :: x
        real(dp) :: w

        do i = 1, size(nodes) / 2
            j = size(nodes) + 1 - i
            x = 0.5_dp * (nodes(j) - nodes(i))
            w = 0.5_dp * (weights(i) + weights(j))
            nodes(i) = -x
            nodes(j) = x
            weights(i) = w
            weights(j) = w
        end do
        if (mod(size(nodes), 2) == 1) nodes((size(nodes) + 1) / 2) = 0.0_dp
    end subroutine symmetrize_rule

    pure subroutine invalid_rule(nodes, weights, mu, status)
        real(dp), allocatable, intent(out) :: nodes(:) !! empty node array on invalid input
        real(dp), allocatable, intent(out) :: weights(:) !! empty weight array on invalid input
        real(dp), intent(out), optional :: mu !! NaN on invalid input
        integer, intent(out), optional :: status !! invalid-argument status when present

        allocate(nodes(0), weights(0))
        if (present(mu)) mu = quiet_nan(0.0_dp)
        if (present(status)) status = orthogonal_roots_status_invalid_argument
    end subroutine invalid_rule

    pure elemental function exp_guarded(log_value) result(value)
        real(dp), intent(in) :: log_value !! natural logarithm of a nonnegative real quantity
        real(dp) :: value

        if (log_value > log(huge(1.0_dp))) then
            value = positive_infinity(0.0_dp)
        else
            value = exp(log_value)
        end if
    end function exp_guarded

end module scifort_orthogonal_roots
