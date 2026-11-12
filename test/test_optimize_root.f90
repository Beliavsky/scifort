program test_optimize_root
    use scifort_kinds, only : dp
    use scifort_optimize, only : fsolve, root, vector_root_result
    use test_optimize_root_reference
    implicit none

    integer :: failures
    type(vector_root_result) :: r

    failures = 0

    r = root(system, [1.0_dp, 1.0_dp])
    call check_result('hybr numeric', r, hybr_x_ref, hybr_fun_ref, hybr_nfev_ref, &
        hybr_njev_ref, hybr_status_ref, 5.0e-10_dp, failures)
    call check_true('hybr default method', trim(r%method) == 'hybr', failures)
    call check_true('hybr r allocated', allocated(r%r), failures)

    r = root(system, [1.0_dp, 1.0_dp], method='HyBr', jac=system_jac)
    call check_result('hybr analytic', r, hybr_jac_x_ref, hybr_jac_fun_ref, &
        hybr_jac_nfev_ref, hybr_jac_njev_ref, hybr_jac_status_ref, 5.0e-10_dp, &
        failures)

    r = fsolve(system, [1.0_dp, 1.0_dp])
    call check_result('fsolve', r, fsolve_x_ref, fsolve_fun_ref, fsolve_nfev_ref, &
        0, fsolve_status_ref, 5.0e-10_dp, failures)

    r = root(system, [1.0_dp, 1.0_dp], method='lm')
    call check_result('lm numeric', r, lm_x_ref, lm_fun_ref, lm_nfev_ref, &
        lm_njev_ref, lm_status_ref, 5.0e-10_dp, failures)
    call check_true('lm pivots allocated', allocated(r%ipvt), failures)
    call check_true('lm pivots one based', all(r%ipvt >= 1), failures)

    r = root(system, [1.0_dp, 1.0_dp], method='LM', jac=system_jac)
    call check_result('lm analytic', r, lm_jac_x_ref, lm_jac_fun_ref, &
        lm_jac_nfev_ref, lm_jac_njev_ref, lm_jac_status_ref, 5.0e-10_dp, failures)

    r = root(overdetermined, [1.0_dp, 1.0_dp], method='lm', m=3)
    call check_result('lm overdetermined', r, lm_over_x_ref, lm_over_fun_ref, &
        lm_over_nfev_ref, lm_over_njev_ref, lm_over_status_ref, 5.0e-9_dp, failures)
    call check_int('lm over fjac rows', size(r%fjac, 1), 2, failures)
    call check_int('lm over fjac columns', size(r%fjac, 2), 3, failures)

    r = root(overdetermined, [1.0_dp, 1.0_dp], method='lm', m=3, &
        jac=overdetermined_jac)
    call check_result('lm over analytic', r, lm_over_jac_x_ref, lm_over_jac_fun_ref, &
        lm_over_jac_nfev_ref, lm_over_jac_njev_ref, lm_over_jac_status_ref, &
        5.0e-10_dp, failures)

    r = root(diagonal_system, [1.0_dp, 1.0_dp], method='hybr', band=[0, 0])
    call check_result('hybr banded', r, hybr_band_x_ref, hybr_band_fun_ref, &
        hybr_band_nfev_ref, hybr_band_njev_ref, hybr_band_status_ref, 5.0e-10_dp, &
        failures)

    r = root(system, [1.0_dp, 1.0_dp], method='unknown')
    call check_int('unknown method status', r%status, 0, failures)
    call check_true('unknown method failure', .not. r%success, failures)

    r = root(system, [1.0_dp, 1.0_dp], method='lm', m=1)
    call check_int('lm invalid m status', r%status, 0, failures)
    call check_true('lm invalid m failure', .not. r%success, failures)

    r = root(system, [1.0_dp, 1.0_dp], method='hybr', band=[-1, 0])
    call check_int('invalid band status', r%status, 0, failures)
    call check_true('invalid band failure', .not. r%success, failures)

    if (failures /= 0) then
        write(*, '(a,i0)') 'test_optimize_root: FAIL ', failures
        stop 1
    end if
    write(*, '(a)') 'test_optimize_root: PASS'

contains

    subroutine system(x, fvec)
        real(dp), intent(in) :: x(:) !! independent variable vector
        real(dp), intent(out) :: fvec(:) !! residual vector

        fvec(1) = x(1)*cos(x(2)) - 4.0_dp
        fvec(2) = x(2)*x(1) - x(2) - 5.0_dp
    end subroutine system

    subroutine system_jac(x, fjac)
        real(dp), intent(in) :: x(:) !! independent variable vector
        real(dp), intent(out) :: fjac(:, :) !! residual-by-variable Jacobian

        fjac(1,1) = cos(x(2))
        fjac(1,2) = -x(1)*sin(x(2))
        fjac(2,1) = x(2)
        fjac(2,2) = x(1) - 1.0_dp
    end subroutine system_jac

    subroutine overdetermined(x, fvec)
        real(dp), intent(in) :: x(:) !! independent variable vector
        real(dp), intent(out) :: fvec(:) !! residual vector

        fvec(1) = x(1) + x(2) - 3.0_dp
        fvec(2) = 2.0_dp*x(1) - x(2)
        fvec(3) = x(1) - 2.0_dp*x(2) + 1.0_dp
    end subroutine overdetermined

    subroutine overdetermined_jac(x, fjac)
        real(dp), intent(in) :: x(:) !! independent variable vector
        real(dp), intent(out) :: fjac(:, :) !! residual-by-variable Jacobian

        if (size(x) /= 2) error stop 'overdetermined_jac: unexpected x size'
        fjac = reshape([1.0_dp, 2.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, -2.0_dp], &
            shape(fjac))
    end subroutine overdetermined_jac

    subroutine diagonal_system(x, fvec)
        real(dp), intent(in) :: x(:) !! independent variable vector
        real(dp), intent(out) :: fvec(:) !! residual vector

        fvec(1) = x(1)*x(1) - 2.0_dp
        fvec(2) = x(2)*x(2) - 3.0_dp
    end subroutine diagonal_system

    subroutine check_result(name, actual, expected_x, expected_fun, expected_nfev, &
            expected_njev, expected_status, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        type(vector_root_result), intent(in) :: actual !! actual solver result
        real(dp), intent(in) :: expected_x(:) !! reference solution vector
        real(dp), intent(in) :: expected_fun(:) !! reference residual vector
        integer, intent(in) :: expected_nfev !! reference objective count
        integer, intent(in) :: expected_njev !! reference Jacobian count
        integer, intent(in) :: expected_status !! reference termination status
        real(dp), intent(in) :: tol !! absolute vector tolerance
        integer, intent(inout) :: nfail !! accumulated failure count

        call check_vector(trim(name)//' x', actual%x, expected_x, tol, nfail)
        call check_vector(trim(name)//' fun', actual%fun, expected_fun, tol, nfail)
        call check_int(trim(name)//' nfev', actual%nfev, expected_nfev, nfail)
        call check_int(trim(name)//' njev', actual%njev, expected_njev, nfail)
        call check_int(trim(name)//' status', actual%status, expected_status, nfail)
        call check_true(trim(name)//' success', actual%success, nfail)
    end subroutine check_result

    subroutine check_vector(name, actual, expected, tol, nfail)
        character(len=*), intent(in) :: name !! check label
        real(dp), intent(in) :: actual(:) !! actual vector
        real(dp), intent(in) :: expected(:) !! reference vector
        real(dp), intent(in) :: tol !! absolute tolerance
        integer, intent(inout) :: nfail !! accumulated failure count

        if (size(actual) /= size(expected)) then
            write(*, '(a,2(1x,i0))') trim(name)//' size mismatch:', size(actual), &
                size(expected)
            nfail = nfail + 1
        else if (maxval(abs(actual - expected)) > tol) then
            write(*, '(a,1x,es24.16)') trim(name)//' max error:', &
                maxval(abs(actual - expected))
            nfail = nfail + 1
        end if
    end subroutine check_vector

    subroutine check_int(name, actual, expected, nfail)
        character(len=*), intent(in) :: name !! check label
        integer, intent(in) :: actual !! actual integer
        integer, intent(in) :: expected !! expected integer
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

end program test_optimize_root
