program test_optimize_root_nonlin
    use scifort_kinds, only : dp
    use scifort_optimize, only : broyden1, newton_krylov, root, vector_root_result
    use test_optimize_root_reference
    implicit none

    integer :: failures
    type(vector_root_result) :: r

    failures = 0

    r = root(nonlinear_system, [0.0_dp, 0.0_dp], method='broyden1')
    call check_result('broyden1', r, broyden1_nonlin_x_ref, broyden1_nonlin_fun_ref, &
        broyden1_nonlin_nfev_ref, broyden1_nonlin_nit_ref, broyden1_nonlin_status_ref, &
        .true., failures)

    r = root(nonlinear_system, [0.0_dp, 0.0_dp], method='broyden2')
    call check_result('broyden2', r, broyden2_nonlin_x_ref, broyden2_nonlin_fun_ref, &
        broyden2_nonlin_nfev_ref, broyden2_nonlin_nit_ref, broyden2_nonlin_status_ref, &
        .true., failures)

    r = root(nonlinear_system, [0.0_dp, 0.0_dp], method='anderson')
    call check_result('anderson', r, anderson_nonlin_x_ref, anderson_nonlin_fun_ref, &
        anderson_nonlin_nfev_ref, anderson_nonlin_nit_ref, anderson_nonlin_status_ref, &
        .true., failures)

    r = root(nonlinear_system, [0.0_dp, 0.0_dp], method='diagbroyden')
    call check_result('diagbroyden', r, diagbroyden_nonlin_x_ref, &
        diagbroyden_nonlin_fun_ref, diagbroyden_nonlin_nfev_ref, &
        diagbroyden_nonlin_nit_ref, diagbroyden_nonlin_status_ref, .true., failures)

    r = root(easy_mixing_system, [0.0_dp, 0.0_dp], method='linearmixing')
    call check_result('linearmixing', r, linearmixing_easy_x_ref, &
        linearmixing_easy_fun_ref, linearmixing_easy_nfev_ref, &
        linearmixing_easy_nit_ref, linearmixing_easy_status_ref, .true., failures)

    r = root(easy_mixing_system, [0.0_dp, 0.0_dp], method='excitingmixing')
    call check_result('excitingmixing', r, excitingmixing_easy_x_ref, &
        excitingmixing_easy_fun_ref, excitingmixing_easy_nfev_ref, &
        excitingmixing_easy_nit_ref, excitingmixing_easy_status_ref, .true., failures)

    r = root(nonlinear_system, [0.0_dp, 0.0_dp], method='broyden1', &
        line_search='wolfe')
    call check_result('broyden1 Wolfe', r, broyden1_wolfe_x_ref, &
        broyden1_wolfe_fun_ref, broyden1_wolfe_nfev_ref, broyden1_wolfe_nit_ref, &
        broyden1_wolfe_status_ref, .true., failures)

    r = root(nonlinear_system, [0.0_dp, 0.0_dp], method='broyden1', &
        max_rank=2, reduction_method='restart')
    call check_result('broyden1 restart rank', r, broyden1_restart_rank_x_ref, &
        broyden1_restart_rank_fun_ref, broyden1_restart_rank_nfev_ref, &
        broyden1_restart_rank_nit_ref, broyden1_restart_rank_status_ref, &
        .true., failures)

    r = root(nonlinear_system, [0.0_dp, 0.0_dp], method='broyden1', &
        max_rank=2, reduction_method='simple')
    call check_result('broyden1 simple rank', r, broyden1_simple_rank_x_ref, &
        broyden1_simple_rank_fun_ref, broyden1_simple_rank_nfev_ref, &
        broyden1_simple_rank_nit_ref, broyden1_simple_rank_status_ref, &
        .true., failures)

    r = root(nonlinear_system, [0.0_dp, 0.0_dp], method='broyden1', &
        max_rank=3, reduction_method='svd')
    call check_result('broyden1 svd rank', r, broyden1_svd_rank_x_ref, &
        broyden1_svd_rank_fun_ref, broyden1_svd_rank_nfev_ref, &
        broyden1_svd_rank_nit_ref, broyden1_svd_rank_status_ref, &
        .true., failures)

    r = root(nonlinear_system, [0.0_dp, 0.0_dp], method='broyden2', &
        max_rank=2, reduction_method='simple')
    call check_result('broyden2 simple rank', r, broyden2_simple_rank_x_ref, &
        broyden2_simple_rank_fun_ref, broyden2_simple_rank_nfev_ref, &
        broyden2_simple_rank_nit_ref, broyden2_simple_rank_status_ref, &
        .false., failures)

    r = broyden1(nonlinear_system, [0.0_dp, 0.0_dp], max_rank=2, &
        reduction_method='restart')
    call check_result('direct broyden1 rank control', r, broyden1_restart_rank_x_ref, &
        broyden1_restart_rank_fun_ref, broyden1_restart_rank_nfev_ref, &
        broyden1_restart_rank_nit_ref, broyden1_restart_rank_status_ref, &
        .true., failures)

    r = root(nonlinear_system, [0.0_dp, 0.0_dp], method='broyden1', &
        line_search='none')
    call check_result('broyden1 no line search', r, broyden1_nolinesearch_x_ref, &
        broyden1_nolinesearch_fun_ref, broyden1_nolinesearch_nfev_ref, &
        broyden1_nolinesearch_nit_ref, broyden1_nolinesearch_status_ref, &
        .true., failures)

    r = root(nonlinear_system, [0.0_dp, 0.0_dp], method='broyden1', tol=1.0e-10_dp)
    call check_result('broyden1 generic tol', r, broyden1_tol_x_ref, &
        broyden1_tol_fun_ref, broyden1_tol_nfev_ref, broyden1_tol_nit_ref, &
        broyden1_tol_status_ref, .true., failures)

    r = root(nonlinear_system, [0.0_dp, 0.0_dp], method='broyden1', maxiter=2)
    call check_result('broyden1 maxiter', r, broyden1_limit_x_ref, &
        broyden1_limit_fun_ref, broyden1_limit_nfev_ref, broyden1_limit_nit_ref, &
        broyden1_limit_status_ref, .false., failures)

    r = broyden1(nonlinear_system, [0.0_dp, 0.0_dp])
    call check_result('direct broyden1', r, broyden1_nonlin_x_ref, &
        broyden1_nonlin_fun_ref, broyden1_nonlin_nfev_ref, broyden1_nonlin_nit_ref, &
        broyden1_nonlin_status_ref, .true., failures)

    r = root(nonlinear_system, [0.0_dp, 0.0_dp], method='df-sane')
    call check_result('df-sane default', r, dfsane_default_x_ref, dfsane_default_fun_ref, &
        dfsane_default_nfev_ref, dfsane_default_nit_ref, dfsane_default_status_ref, &
        .true., failures)

    r = root(nonlinear_system, [0.0_dp, 0.0_dp], method='df-sane', &
        line_search='cheng')
    call check_result('df-sane cheng', r, dfsane_cheng_x_ref, dfsane_cheng_fun_ref, &
        dfsane_cheng_nfev_ref, dfsane_cheng_nit_ref, dfsane_cheng_status_ref, &
        .true., failures)

    r = root(nonlinear_system, [0.0_dp, 0.0_dp], method='df-sane', tol=1.0e-12_dp)
    call check_result('df-sane generic tol', r, dfsane_tol_x_ref, dfsane_tol_fun_ref, &
        dfsane_tol_nfev_ref, dfsane_tol_nit_ref, dfsane_tol_status_ref, .true., failures)

    r = root(nonlinear_system, [0.0_dp, 0.0_dp], method='df-sane', maxfev=2)
    call check_result('df-sane maxfev', r, dfsane_limit_x_ref, dfsane_limit_fun_ref, &
        dfsane_limit_nfev_ref, dfsane_limit_nit_ref, dfsane_limit_status_ref, &
        .false., failures)

    r = root(krylov_linear_system, [0.0_dp, 0.0_dp], method='krylov')
    call check_result('krylov default', r, krylov_default_x_ref, &
        krylov_default_fun_ref, krylov_default_nfev_ref, krylov_default_nit_ref, &
        krylov_default_status_ref, .true., failures)

    r = root(krylov_linear_system, [0.0_dp, 0.0_dp], method='krylov', rdiff=1.0e-6_dp)
    call check_result('krylov rdiff', r, krylov_rdiff_x_ref, krylov_rdiff_fun_ref, &
        krylov_rdiff_nfev_ref, krylov_rdiff_nit_ref, krylov_rdiff_status_ref, &
        .true., failures)

    r = root(diagonal_system, [1.0_dp, 1.0_dp], method='krylov', outer_k=0)
    call check_result('krylov outer_k=0', r, krylov_outer0_x_ref, &
        krylov_outer0_fun_ref, krylov_outer0_nfev_ref, krylov_outer0_nit_ref, &
        krylov_outer0_status_ref, .true., failures)

    r = root(diagonal_system, [1.0_dp, 1.0_dp], method='krylov', inner_method='gmres')
    call check_result('krylov gmres', r, krylov_gmres_x_ref, krylov_gmres_fun_ref, &
        krylov_gmres_nfev_ref, krylov_gmres_nit_ref, krylov_gmres_status_ref, &
        .true., failures)

    r = root(diagonal_system, [1.0_dp, 1.0_dp], method='krylov', &
        line_search='wolfe', outer_k=0)
    call check_result('krylov Wolfe', r, krylov_wolfe_x_ref, krylov_wolfe_fun_ref, &
        krylov_wolfe_nfev_ref, krylov_wolfe_nit_ref, krylov_wolfe_status_ref, &
        .true., failures)

    r = root(diagonal_system, [1.0_dp, 1.0_dp], method='krylov', maxiter=2, outer_k=0)
    call check_result('krylov maxiter', r, krylov_limit_x_ref, krylov_limit_fun_ref, &
        krylov_limit_nfev_ref, krylov_limit_nit_ref, krylov_limit_status_ref, &
        .false., failures)

    r = root(diagonal_system, [1.0_dp, 1.0_dp], method='krylov', &
        fatol=1.0e-10_dp, outer_k=0)
    call check_result('krylov fatol', r, krylov_fatol_x_ref, krylov_fatol_fun_ref, &
        krylov_fatol_nfev_ref, krylov_fatol_nit_ref, krylov_fatol_status_ref, &
        .true., failures)

    r = newton_krylov(krylov_linear_system, [0.0_dp, 0.0_dp])
    call check_result('direct newton_krylov', r, krylov_default_x_ref, &
        krylov_default_fun_ref, krylov_default_nfev_ref, krylov_default_nit_ref, &
        krylov_default_status_ref, .true., failures)

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_root_nonlin: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_optimize_root_nonlin: PASS'

contains

    subroutine nonlinear_system(x, fvec)
        real(dp), intent(in) :: x(:) !! independent variable vector
        real(dp), intent(out) :: fvec(:) !! residual vector

        fvec(1) = x(1) + 0.5_dp * (x(1) - x(2))**3 - 1.0_dp
        fvec(2) = 0.5_dp * (x(2) - x(1))**3 + x(2)
    end subroutine nonlinear_system

    subroutine easy_mixing_system(x, fvec)
        real(dp), intent(in) :: x(:) !! independent variable vector
        real(dp), intent(out) :: fvec(:) !! residual vector

        fvec(1) = 1.0_dp - x(1)
        fvec(2) = 2.0_dp - x(2)
    end subroutine easy_mixing_system

    subroutine krylov_linear_system(x, fvec)
        real(dp), intent(in) :: x(:) !! independent variable vector
        real(dp), intent(out) :: fvec(:) !! residual vector

        fvec(1) = 3.0_dp * x(1) + x(2) - 1.0_dp
        fvec(2) = x(1) + 2.0_dp * x(2) - 2.0_dp
    end subroutine krylov_linear_system

    subroutine diagonal_system(x, fvec)
        real(dp), intent(in) :: x(:) !! independent variable vector
        real(dp), intent(out) :: fvec(:) !! residual vector

        fvec(1) = x(1)**2 - 2.0_dp
        fvec(2) = x(2)**2 - 3.0_dp
    end subroutine diagonal_system

    subroutine check_result(name, actual, expected_x, expected_fun, expected_nfev, &
            expected_nit, expected_status, expected_success, nfail)
        character(len=*), intent(in) :: name !! check label
        type(vector_root_result), intent(in) :: actual !! actual root result
        real(dp), intent(in) :: expected_x(:) !! reference solution
        real(dp), intent(in) :: expected_fun(:) !! reference residual
        integer, intent(in) :: expected_nfev !! reference residual evaluations
        integer, intent(in) :: expected_nit !! reference iteration count
        integer, intent(in) :: expected_status !! reference status
        logical, intent(in) :: expected_success !! reference success flag
        integer, intent(inout) :: nfail !! accumulated failure count

        call check_vector(trim(name)//' x', actual%x, expected_x, 5.0e-12_dp, nfail)
        call check_vector(trim(name)//' fun', actual%fun, expected_fun, 5.0e-12_dp, nfail)
        call check_int(trim(name)//' nfev', actual%nfev, expected_nfev, nfail)
        call check_int(trim(name)//' nit', actual%nit, expected_nit, nfail)
        call check_int(trim(name)//' status', actual%status, expected_status, nfail)
        call check_true(trim(name)//' success', actual%success .eqv. expected_success, nfail)
    end subroutine check_result

    subroutine check_vector(name, actual, expected, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual(:) !! actual vector
        real(dp), intent(in) :: expected(:) !! reference vector
        real(dp), intent(in) :: tol !! absolute tolerance
        integer, intent(inout) :: nfail !! accumulated failure count

        if (size(actual) /= size(expected)) then
            write(*, '(a,2(1x,i0))') trim(name)//' size mismatch:', &
                size(actual), size(expected)
            nfail = nfail + 1
            return
        end if
        if (maxval(abs(actual - expected)) > tol) then
            write(*, '(a)') trim(name)//' mismatch'
            write(*, '(a,*(1x,es24.16))') '  actual:  ', actual
            write(*, '(a,*(1x,es24.16))') '  expected:', expected
            nfail = nfail + 1
        end if
    end subroutine check_vector

    subroutine check_int(name, actual, expected, nfail)
        character(len=*), intent(in) :: name !! check label
        integer, intent(in) :: actual !! actual integer
        integer, intent(in) :: expected !! reference integer
        integer, intent(inout) :: nfail !! accumulated failure count

        if (actual /= expected) then
            write(*, '(a,2(1x,i0))') trim(name)//' mismatch:', actual, expected
            nfail = nfail + 1
        end if
    end subroutine check_int

    subroutine check_true(name, condition, nfail)
        character(len=*), intent(in) :: name !! check label
        logical, intent(in) :: condition !! condition that must hold
        integer, intent(inout) :: nfail !! accumulated failure count

        if (.not. condition) then
            write(*, '(a)') trim(name)//' failed'
            nfail = nfail + 1
        end if
    end subroutine check_true

end program test_optimize_root_nonlin
