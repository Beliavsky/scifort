! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_special_orthogonal_aliases
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_special, only : assoc_laguerre, c_roots, cg_roots, dp, eval_genlaguerre, &
        h_roots, he_roots, j_roots, js_roots, l_roots, la_roots, p_roots, ps_roots, &
        roots_chebyc, roots_chebys, roots_chebyt, roots_chebyu, roots_genlaguerre, &
        roots_gegenbauer, roots_hermite, roots_hermitenorm, roots_jacobi, &
        roots_laguerre, roots_legendre, roots_sh_chebyt, roots_sh_chebyu, &
        roots_sh_jacobi, roots_sh_legendre, s_roots, t_roots, ts_roots, u_roots, us_roots
    implicit none

    complex(dp) :: z
    real(dp), allocatable :: x_alias(:)
    real(dp), allocatable :: w_alias(:)
    real(dp), allocatable :: x_root(:)
    real(dp), allocatable :: w_root(:)
    real(dp) :: mu_alias
    real(dp) :: mu_root
    integer :: status_alias
    integer :: status_root

    call check_close('assoc laguerre explicit order', assoc_laguerre(0.7_dp, 4, 1.2_dp), &
        0.6409375000000008_dp, 2.0e-15_dp)
    call check_close('assoc laguerre default order', assoc_laguerre(-2.3_dp, 5), &
        65.59471608333332_dp, 2.0e-13_dp)
    call check_close('assoc laguerre wrapper identity', assoc_laguerre(0.4_dp, 6, 0.75_dp), &
        eval_genlaguerre(6, 0.75_dp, 0.4_dp), 0.0_dp)
    z = assoc_laguerre(cmplx(0.2_dp, 0.3_dp, dp), 3, 0.5_dp)
    call check_complex_close('assoc laguerre complex', z, &
        cmplx(1.2326666666666666_dp, -1.104_dp, dp), 2.0e-15_dp)
    call check_true('assoc laguerre invalid order', ieee_is_nan(assoc_laguerre(0.5_dp, 3, -1.0_dp)))

    call p_roots(5, x_alias, w_alias, mu_alias, status_alias)
    call roots_legendre(5, x_root, w_root, mu_root, status_root)
    call check_rule_alias('p_roots', x_alias, w_alias, mu_alias, status_alias, &
        x_root, w_root, mu_root, status_root)

    call j_roots(5, 0.2_dp, 1.3_dp, x_alias, w_alias, mu_alias, status_alias)
    call roots_jacobi(5, 0.2_dp, 1.3_dp, x_root, w_root, mu_root, status_root)
    call check_rule_alias('j_roots', x_alias, w_alias, mu_alias, status_alias, &
        x_root, w_root, mu_root, status_root)

    call cg_roots(5, 0.7_dp, x_alias, w_alias, mu_alias, status_alias)
    call roots_gegenbauer(5, 0.7_dp, x_root, w_root, mu_root, status_root)
    call check_rule_alias('cg_roots', x_alias, w_alias, mu_alias, status_alias, &
        x_root, w_root, mu_root, status_root)

    call l_roots(5, x_alias, w_alias, mu_alias, status_alias)
    call roots_laguerre(5, x_root, w_root, mu_root, status_root)
    call check_rule_alias('l_roots', x_alias, w_alias, mu_alias, status_alias, &
        x_root, w_root, mu_root, status_root)

    call la_roots(5, 1.2_dp, x_alias, w_alias, mu_alias, status_alias)
    call roots_genlaguerre(5, 1.2_dp, x_root, w_root, mu_root, status_root)
    call check_rule_alias('la_roots', x_alias, w_alias, mu_alias, status_alias, &
        x_root, w_root, mu_root, status_root)

    call h_roots(5, x_alias, w_alias, mu_alias, status_alias)
    call roots_hermite(5, x_root, w_root, mu_root, status_root)
    call check_rule_alias('h_roots', x_alias, w_alias, mu_alias, status_alias, &
        x_root, w_root, mu_root, status_root)

    call he_roots(5, x_alias, w_alias, mu_alias, status_alias)
    call roots_hermitenorm(5, x_root, w_root, mu_root, status_root)
    call check_rule_alias('he_roots', x_alias, w_alias, mu_alias, status_alias, &
        x_root, w_root, mu_root, status_root)

    call t_roots(5, x_alias, w_alias, mu_alias, status_alias)
    call roots_chebyt(5, x_root, w_root, mu_root, status_root)
    call check_rule_alias('t_roots', x_alias, w_alias, mu_alias, status_alias, &
        x_root, w_root, mu_root, status_root)

    call u_roots(5, x_alias, w_alias, mu_alias, status_alias)
    call roots_chebyu(5, x_root, w_root, mu_root, status_root)
    call check_rule_alias('u_roots', x_alias, w_alias, mu_alias, status_alias, &
        x_root, w_root, mu_root, status_root)

    call c_roots(5, x_alias, w_alias, mu_alias, status_alias)
    call roots_chebyc(5, x_root, w_root, mu_root, status_root)
    call check_rule_alias('c_roots', x_alias, w_alias, mu_alias, status_alias, &
        x_root, w_root, mu_root, status_root)

    call s_roots(5, x_alias, w_alias, mu_alias, status_alias)
    call roots_chebys(5, x_root, w_root, mu_root, status_root)
    call check_rule_alias('s_roots', x_alias, w_alias, mu_alias, status_alias, &
        x_root, w_root, mu_root, status_root)

    call ps_roots(5, x_alias, w_alias, mu_alias, status_alias)
    call roots_sh_legendre(5, x_root, w_root, mu_root, status_root)
    call check_rule_alias('ps_roots', x_alias, w_alias, mu_alias, status_alias, &
        x_root, w_root, mu_root, status_root)

    call js_roots(5, 2.5_dp, 1.2_dp, x_alias, w_alias, mu_alias, status_alias)
    call roots_sh_jacobi(5, 2.5_dp, 1.2_dp, x_root, w_root, mu_root, status_root)
    call check_rule_alias('js_roots', x_alias, w_alias, mu_alias, status_alias, &
        x_root, w_root, mu_root, status_root)

    call ts_roots(5, x_alias, w_alias, mu_alias, status_alias)
    call roots_sh_chebyt(5, x_root, w_root, mu_root, status_root)
    call check_rule_alias('ts_roots', x_alias, w_alias, mu_alias, status_alias, &
        x_root, w_root, mu_root, status_root)

    call us_roots(5, x_alias, w_alias, mu_alias, status_alias)
    call roots_sh_chebyu(5, x_root, w_root, mu_root, status_root)
    call check_rule_alias('us_roots', x_alias, w_alias, mu_alias, status_alias, &
        x_root, w_root, mu_root, status_root)

    print *, 'special orthogonal compatibility alias tests passed'

contains

    subroutine check_rule_alias(name, x_alias, w_alias, mu_alias, status_alias, &
            x_root, w_root, mu_root, status_root)
        character(len=*), intent(in) :: name !! alias under test
        real(dp), intent(in) :: x_alias(:) !! nodes returned by compatibility alias
        real(dp), intent(in) :: w_alias(:) !! weights returned by compatibility alias
        real(dp), intent(in) :: mu_alias !! weight sum returned by compatibility alias
        integer, intent(in) :: status_alias !! status returned by compatibility alias
        real(dp), intent(in) :: x_root(:) !! nodes returned by canonical routine
        real(dp), intent(in) :: w_root(:) !! weights returned by canonical routine
        real(dp), intent(in) :: mu_root !! weight sum returned by canonical routine
        integer, intent(in) :: status_root !! status returned by canonical routine

        call check_true(trim(name)//' node size', size(x_alias) == size(x_root))
        call check_true(trim(name)//' weight size', size(w_alias) == size(w_root))
        call check_true(trim(name)//' nodes', all(x_alias == x_root))
        call check_true(trim(name)//' weights', all(w_alias == w_root))
        call check_true(trim(name)//' mu', mu_alias == mu_root)
        call check_true(trim(name)//' status', status_alias == status_root)
    end subroutine check_rule_alias

    subroutine check_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (actual == expected) return
        if (ieee_is_nan(actual) .or. abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_close

    subroutine check_complex_close(name, actual, expected, tolerance)
        character(len=*), intent(in) :: name !! check label printed on failure
        complex(dp), intent(in) :: actual !! computed complex value
        complex(dp), intent(in) :: expected !! reference complex value
        real(dp), intent(in) :: tolerance !! maximum absolute error

        if (actual == expected) return
        if (abs(actual - expected) > tolerance) then
            print *, 'FAIL ', name
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance=', tolerance
            error stop 1
        end if
    end subroutine check_complex_close

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name !! check label printed on failure
        logical, intent(in) :: condition !! condition expected to be true

        if (.not. condition) then
            print *, 'FAIL ', name
            error stop 1
        end if
    end subroutine check_true

end program test_special_orthogonal_aliases
