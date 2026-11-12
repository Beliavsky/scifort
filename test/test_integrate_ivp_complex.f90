program test_integrate_ivp_complex
    use scifort_integrate, only : dp, ivp_complex_event_set, &
        ivp_complex_sparse_matrix, ivp_csc_matrix, ivp_status_event, &
        ivp_status_invalid_argument, ivp_status_success, solve_ivp, &
        solve_ivp_complex_result
    use scifort_integrate_bdf_complex, only : bdf_complex_initialize, bdf_complex_state
    implicit none

    integer :: complex_stiff_rhs_calls = 0
    integer :: complex_vectorized_calls = 0
    integer :: complex_vectorized_max_columns = 0
    integer :: complex_sparse_jacobian_calls = 0

    call test_complex_rk23()
    call test_complex_vector_atol()
    call test_complex_vector_rtol()
    call test_complex_constant_jacobian_matrix()
    call test_complex_sparse_jacobian_matrix()
    call test_complex_sparse_jacobian_callback()
    call test_complex_sparse_numerical_storage()
    call test_complex_jacobian_sparsity()
    call test_complex_vectorized_rhs_jacobian()
    call test_complex_invalid_jacobian_options()
    call test_complex_rk45()
    call test_complex_dop853()
    call test_complex_backward_dop853()
    call test_complex_bdf()
    call test_complex_bdf_analytic_jacobian()
    call test_complex_bdf_stiff_finite_difference()
    call test_complex_bdf_coupled_finite_difference()
    call test_complex_backward_bdf()
    call test_complex_bdf_terminal_event()
    call test_complex_terminal_event()
    call test_complex_event_occurrence_limit()
    call test_complex_vector_events()
    call test_complex_invalid_method()

    print *, 'test_integrate_ivp_complex: PASS'

contains

    subroutine test_complex_vector_atol()
        character(len=6), parameter :: methods(4) = [character(len=6) :: &
            'RK23', 'RK45', 'DOP853', 'BDF']
        type(solve_ivp_complex_result) :: relaxed
        type(solve_ivp_complex_result) :: tight
        type(solve_ivp_complex_result) :: invalid
        complex(dp) :: expected
        integer :: i

        expected = exp(cmplx(-2.0_dp, 2.0_dp, dp))
        do i = 1, size(methods)
            relaxed = solve_ivp(complex_multiscale_decay_rhs, [0.0_dp, 2.0_dp], &
                [cmplx(1.0_dp, 0.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp)], &
                method=trim(methods(i)), rtol=1.0e-6_dp, &
                atol=[1.0e-9_dp, 1.0e-3_dp])
            tight = solve_ivp(complex_multiscale_decay_rhs, [0.0_dp, 2.0_dp], &
                [cmplx(1.0_dp, 0.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp)], &
                method=trim(methods(i)), rtol=1.0e-6_dp, atol=1.0e-9_dp)

            call require(relaxed%success, trim(methods(i))//' complex vector atol success')
            call require(tight%success, trim(methods(i))//' complex scalar atol success')
            call check_complex_abs(relaxed%y(1, size(relaxed%t)), expected, 5.0e-6_dp, &
                trim(methods(i))//' complex vector atol first component')
            call require(abs(relaxed%y(2, size(relaxed%t))) <= 1.0e-3_dp, &
                trim(methods(i))//' complex vector atol second component')
            call require(relaxed%nfev < tight%nfev, &
                trim(methods(i))//' complex vector atol reduces RHS work')
        end do

        invalid = solve_ivp(complex_multiscale_decay_rhs, [0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp)], &
            atol=[1.0e-6_dp])
        call require(.not. invalid%success, 'complex short vector atol rejected')
        call require(invalid%status == ivp_status_invalid_argument, &
            'complex short vector atol status')
    end subroutine test_complex_vector_atol

    subroutine test_complex_vector_rtol()
        character(len=6), parameter :: methods(4) = [character(len=6) :: &
            'RK23', 'RK45', 'DOP853', 'BDF']
        type(solve_ivp_complex_result) :: scalar_tol
        type(solve_ivp_complex_result) :: vector_tol
        type(solve_ivp_complex_result) :: invalid
        complex(dp) :: expected_fast
        complex(dp) :: expected_slow
        real(dp) :: rank_two_rtol(1, 2)
        integer :: i

        expected_slow = exp(cmplx(-0.2_dp, 0.2_dp, dp))
        expected_fast = exp(cmplx(-2.0_dp, 0.4_dp, dp))
        do i = 1, size(methods)
            scalar_tol = solve_ivp(complex_multiscale_decay_rhs, [0.0_dp, 0.2_dp], &
                [cmplx(1.0_dp, 0.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp)], &
                method=trim(methods(i)), rtol=1.0e-3_dp, atol=1.0e-12_dp)
            vector_tol = solve_ivp(complex_multiscale_decay_rhs, [0.0_dp, 0.2_dp], &
                [cmplx(1.0_dp, 0.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp)], &
                method=trim(methods(i)), rtol=[1.0e-3_dp, 1.0e-8_dp], atol=1.0e-12_dp)

            call require(scalar_tol%success, trim(methods(i))//' complex scalar rtol success')
            call require(vector_tol%success, trim(methods(i))//' complex vector rtol success')
            call check_complex_abs(vector_tol%y(1, size(vector_tol%t)), expected_slow, &
                5.0e-4_dp, trim(methods(i))//' complex vector rtol slow component')
            call check_complex_abs(vector_tol%y(2, size(vector_tol%t)), expected_fast, &
                2.0e-8_dp, trim(methods(i))//' complex vector rtol fast component')
            call require(vector_tol%nfev > scalar_tol%nfev, &
                trim(methods(i))//' tighter complex vector rtol increases RHS work')
        end do

        invalid = solve_ivp(complex_multiscale_decay_rhs, [0.0_dp, 0.2_dp], &
            [cmplx(1.0_dp, 0.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp)], &
            rtol=[1.0e-3_dp])
        call require(.not. invalid%success, 'complex short vector rtol rejected')
        call require(invalid%status == ivp_status_invalid_argument, &
            'complex short vector rtol status')

        rank_two_rtol = 1.0e-3_dp
        invalid = solve_ivp(complex_multiscale_decay_rhs, [0.0_dp, 0.2_dp], &
            [cmplx(1.0_dp, 0.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp)], &
            rtol=rank_two_rtol)
        call require(.not. invalid%success, 'complex rank-two rtol rejected')
        call require(invalid%status == ivp_status_invalid_argument, &
            'complex rank-two rtol status')
    end subroutine test_complex_vector_rtol

    subroutine test_complex_constant_jacobian_matrix()
        type(solve_ivp_complex_result) :: sol
        complex(dp) :: jac_matrix(2, 2)
        complex(dp) :: y0(2)

        y0 = [cmplx(1.0_dp, 0.5_dp, dp), cmplx(-0.25_dp, 0.75_dp, dp)]
        jac_matrix = cmplx(0.0_dp, 0.0_dp, dp)
        jac_matrix(1, 1) = cmplx(-100.0_dp, 30.0_dp, dp)
        jac_matrix(2, 2) = cmplx(-1.0_dp, 2.0_dp, dp)
        complex_stiff_rhs_calls = 0
        sol = solve_ivp(complex_stiff_rhs, [0.0_dp, 1.0_dp], y0, method='BDF', &
            rtol=1.0e-8_dp, atol=1.0e-10_dp, jac_matrix=jac_matrix)

        call require(sol%success, 'complex BDF constant Jacobian success')
        call require(sol%nfev == 562, 'complex BDF constant Jacobian nfev')
        call require(sol%njev == 0, 'complex BDF constant Jacobian njev')
        call require(sol%nlu == 51, 'complex BDF constant Jacobian nlu')
        call require(sol%nfev == complex_stiff_rhs_calls, &
            'complex BDF constant Jacobian has no hidden RHS probes')
        call check_complex_abs(sol%y(1, size(sol%t)), &
            cmplx(-1.4230211839551403e-14_dp, -1.4265282999130226e-14_dp, dp), &
            5.0e-12_dp, 'complex BDF constant Jacobian first component')
        call check_complex_abs(sol%y(2, size(sol%t)), &
            cmplx(-0.21261093166186007_dp, -0.19844686050431756_dp, dp), &
            5.0e-12_dp, 'complex BDF constant Jacobian second component')
    end subroutine test_complex_constant_jacobian_matrix

    subroutine test_complex_sparse_jacobian_matrix()
        type(ivp_complex_sparse_matrix) :: jac_sparse
        type(solve_ivp_complex_result) :: sol
        complex(dp) :: data(2)
        complex(dp) :: y0(2)

        y0 = [cmplx(1.0_dp, 0.5_dp, dp), cmplx(-0.25_dp, 0.75_dp, dp)]
        data = [cmplx(-100.0_dp, 30.0_dp, dp), cmplx(-1.0_dp, 2.0_dp, dp)]
        jac_sparse = ivp_csc_matrix(2, 2, [1, 2, 3], [1, 2], data)
        call require(jac_sparse%is_valid(), 'complex CSC Jacobian valid')
        call require(jac_sparse%nnz() == 2, 'complex CSC Jacobian nnz')

        complex_stiff_rhs_calls = 0
        sol = solve_ivp(complex_stiff_rhs, [0.0_dp, 1.0_dp], y0, method='BDF', &
            rtol=1.0e-8_dp, atol=1.0e-10_dp, jac_sparse_matrix=jac_sparse)
        call require(sol%success, 'complex BDF sparse matrix Jacobian success')
        call require(sol%nfev == 562, 'complex BDF sparse matrix Jacobian nfev')
        call require(sol%njev == 0, 'complex BDF sparse matrix Jacobian njev')
        call require(sol%nlu == 51, 'complex BDF sparse matrix Jacobian nlu')
        call require(sol%nfev == complex_stiff_rhs_calls, &
            'complex BDF sparse matrix Jacobian has no hidden RHS probes')
        call check_complex_abs(sol%y(1, size(sol%t)), &
            cmplx(-1.4230211839551403e-14_dp, -1.4265282999130226e-14_dp, dp), &
            5.0e-12_dp, 'complex BDF sparse matrix Jacobian first component')
        call check_complex_abs(sol%y(2, size(sol%t)), &
            cmplx(-0.21261093166186007_dp, -0.19844686050431756_dp, dp), &
            5.0e-12_dp, 'complex BDF sparse matrix Jacobian second component')
    end subroutine test_complex_sparse_jacobian_matrix

    subroutine test_complex_sparse_jacobian_callback()
        type(solve_ivp_complex_result) :: sol
        complex(dp) :: y0(2)

        y0 = [cmplx(1.0_dp, 0.5_dp, dp), cmplx(-0.25_dp, 0.75_dp, dp)]
        complex_stiff_rhs_calls = 0
        complex_sparse_jacobian_calls = 0
        sol = solve_ivp(complex_stiff_rhs, [0.0_dp, 1.0_dp], y0, method='BDF', &
            rtol=1.0e-8_dp, atol=1.0e-10_dp, jac_sparse=complex_stiff_sparse_jacobian)
        call require(sol%success, 'complex BDF sparse callback Jacobian success')
        call require(sol%njev == complex_sparse_jacobian_calls, &
            'complex BDF sparse callback Jacobian count')
        call require(sol%njev > 0, 'complex BDF sparse callback Jacobian used')
        call require(sol%nlu > 0, 'complex BDF sparse callback Jacobian LU work')
        call require(sol%nfev == complex_stiff_rhs_calls, &
            'complex BDF sparse callback has no hidden RHS probes')
        call check_complex_abs(sol%y(2, size(sol%t)), &
            cmplx(-0.21261093166186007_dp, -0.19844686050431756_dp, dp), &
            5.0e-12_dp, 'complex BDF sparse callback Jacobian value')
    end subroutine test_complex_sparse_jacobian_callback

    subroutine test_complex_sparse_numerical_storage()
        type(bdf_complex_state) :: bdf
        complex(dp) :: f(2)
        complex(dp) :: y0(2)
        logical :: sparsity(2, 2)
        integer :: njev

        y0 = [cmplx(1.0_dp, 0.5_dp, dp), cmplx(-0.25_dp, 0.75_dp, dp)]
        sparsity = .false.
        sparsity(1, 1) = .true.
        sparsity(2, 2) = .true.
        call complex_stiff_rhs(0.0_dp, y0, f)

        njev = 0
        call bdf_complex_initialize(bdf, complex_stiff_rhs, 0.0_dp, y0, f, 1.0_dp, &
            1.0e-3_dp, 1.0e-8_dp, 1.0e-10_dp, njev, jac_sparsity=sparsity)
        call require(bdf%jacobian_sparse, 'complex BDF sparsity uses sparse Jacobian state')
        call require(.not. allocated(bdf%jac), &
            'complex BDF sparsity avoids dense Jacobian storage')
        call require(size(bdf%lu) == 0, 'complex BDF sparsity avoids dense LU workspace')
        call require(bdf%sparse_jac%is_valid(), 'complex BDF numerical CSC Jacobian valid')
        call require(bdf%sparse_jac%nnz() == 2, 'complex BDF numerical CSC Jacobian nnz')
    end subroutine test_complex_sparse_numerical_storage

    subroutine test_complex_jacobian_sparsity()
        type(solve_ivp_complex_result) :: dense_sol
        type(solve_ivp_complex_result) :: sparse_sol
        complex(dp) :: y0(2)
        logical :: sparsity(2, 2)
        integer :: dense_calls
        integer :: sparse_calls

        y0 = [cmplx(1.0_dp, 0.5_dp, dp), cmplx(-0.25_dp, 0.75_dp, dp)]
        sparsity = .false.
        sparsity(1, 1) = .true.
        sparsity(2, 2) = .true.

        complex_stiff_rhs_calls = 0
        dense_sol = solve_ivp(complex_stiff_rhs, [0.0_dp, 1.0_dp], y0, &
            method='BDF', rtol=1.0e-8_dp, atol=1.0e-10_dp)
        dense_calls = complex_stiff_rhs_calls
        call require(dense_sol%success, 'complex BDF dense Jacobian success')

        complex_stiff_rhs_calls = 0
        sparse_sol = solve_ivp(complex_stiff_rhs, [0.0_dp, 1.0_dp], y0, &
            method='BDF', rtol=1.0e-8_dp, atol=1.0e-10_dp, &
            jac_sparsity=sparsity)
        sparse_calls = complex_stiff_rhs_calls
        call require(sparse_sol%success, 'complex BDF sparse Jacobian success')
        call require(sparse_sol%nfev == dense_sol%nfev, &
            'complex BDF sparsity preserves reported nfev')
        call require(sparse_sol%njev == dense_sol%njev, &
            'complex BDF sparsity preserves Jacobian count')
        call require(sparse_calls < dense_calls, &
            'complex BDF sparsity reduces RHS probes')
        call check_complex_abs(sparse_sol%y(2, size(sparse_sol%t)), &
            dense_sol%y(2, size(dense_sol%t)), 1.0e-13_dp, &
            'complex BDF sparse Jacobian value')
    end subroutine test_complex_jacobian_sparsity

    subroutine test_complex_vectorized_rhs_jacobian()
        type(solve_ivp_complex_result) :: scalar_sol
        type(solve_ivp_complex_result) :: sparse_vec_sol
        type(solve_ivp_complex_result) :: vec_sol
        complex(dp) :: y0(2)
        logical :: sparsity(2, 2)
        integer :: scalar_calls

        y0 = [cmplx(1.0_dp, 0.5_dp, dp), cmplx(-0.25_dp, 0.75_dp, dp)]
        sparsity = .false.
        sparsity(1, 1) = .true.
        sparsity(2, 2) = .true.

        complex_stiff_rhs_calls = 0
        scalar_sol = solve_ivp(complex_stiff_rhs, [0.0_dp, 1.0_dp], y0, &
            method='BDF', rtol=1.0e-8_dp, atol=1.0e-10_dp)
        scalar_calls = complex_stiff_rhs_calls
        call require(scalar_sol%success, 'complex scalar FD success')

        complex_stiff_rhs_calls = 0
        complex_vectorized_calls = 0
        complex_vectorized_max_columns = 0
        vec_sol = solve_ivp(complex_stiff_rhs, [0.0_dp, 1.0_dp], y0, &
            method='BDF', rtol=1.0e-8_dp, atol=1.0e-10_dp, &
            fun_vectorized=complex_stiff_rhs_vectorized)
        call require(vec_sol%success, 'complex vectorized FD success')
        call require(vec_sol%nfev == scalar_sol%nfev, &
            'complex vectorized FD preserves nfev')
        call require(vec_sol%njev == scalar_sol%njev, &
            'complex vectorized FD preserves njev')
        call require(complex_stiff_rhs_calls == vec_sol%nfev, &
            'complex vectorized FD keeps hidden probes batched')
        call require(complex_vectorized_calls > 0, 'complex vectorized FD callback used')
        call require(complex_vectorized_max_columns == size(y0), &
            'complex dense FD batches all perturbed columns')
        call require(complex_stiff_rhs_calls + complex_vectorized_calls < scalar_calls, &
            'complex vectorized FD reduces callback invocations')

        complex_stiff_rhs_calls = 0
        complex_vectorized_calls = 0
        complex_vectorized_max_columns = 0
        sparse_vec_sol = solve_ivp(complex_stiff_rhs, [0.0_dp, 1.0_dp], y0, &
            method='BDF', rtol=1.0e-8_dp, atol=1.0e-10_dp, jac_sparsity=sparsity, &
            fun_vectorized=complex_stiff_rhs_vectorized)
        call require(sparse_vec_sol%success, 'complex vectorized sparse FD success')
        call require(complex_vectorized_max_columns == 1, &
            'complex sparse grouping controls vectorized batch width')
        call check_complex_abs(sparse_vec_sol%y(2, size(sparse_vec_sol%t)), &
            scalar_sol%y(2, size(scalar_sol%t)), 1.0e-13_dp, &
            'complex vectorized sparse FD value')
    end subroutine test_complex_vectorized_rhs_jacobian

    subroutine test_complex_invalid_jacobian_options()
        type(ivp_complex_sparse_matrix) :: bad_sparse
        type(ivp_complex_sparse_matrix) :: sparse_jac
        type(solve_ivp_complex_result) :: sol
        complex(dp) :: bad_matrix(1, 1)
        complex(dp) :: jac_matrix(2, 2)
        logical :: bad_sparsity(1, 1)
        logical :: sparsity(2, 2)

        bad_matrix = cmplx(-1.0_dp, 0.0_dp, dp)
        jac_matrix = cmplx(0.0_dp, 0.0_dp, dp)
        jac_matrix(1, 1) = cmplx(-100.0_dp, 30.0_dp, dp)
        jac_matrix(2, 2) = cmplx(-1.0_dp, 2.0_dp, dp)
        bad_sparsity = .true.
        sparsity = .false.
        sparsity(1, 1) = .true.
        sparsity(2, 2) = .true.
        sparse_jac = ivp_csc_matrix(2, 2, [1, 2, 3], [1, 2], &
            [cmplx(-100.0_dp, 30.0_dp, dp), cmplx(-1.0_dp, 2.0_dp, dp)])
        bad_sparse = ivp_csc_matrix(2, 2, [1, 3, 2], [1, 2], &
            [cmplx(-100.0_dp, 30.0_dp, dp), cmplx(-1.0_dp, 2.0_dp, dp)])

        sol = solve_ivp(complex_stiff_rhs, [0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.5_dp, dp), cmplx(-0.25_dp, 0.75_dp, dp)], &
            method='BDF', jac_matrix=bad_matrix)
        call require(.not. sol%success, 'complex wrong-shape Jacobian rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'complex wrong-shape Jacobian status')

        sol = solve_ivp(complex_stiff_rhs, [0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.5_dp, dp), cmplx(-0.25_dp, 0.75_dp, dp)], &
            method='BDF', jac_sparsity=bad_sparsity)
        call require(.not. sol%success, 'complex wrong-shape sparsity rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'complex wrong-shape sparsity status')

        sol = solve_ivp(complex_stiff_rhs, [0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.5_dp, dp), cmplx(-0.25_dp, 0.75_dp, dp)], &
            method='BDF', jac_sparse_matrix=bad_sparse)
        call require(.not. sol%success, 'complex invalid CSC Jacobian rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'complex invalid CSC Jacobian status')

        sol = solve_ivp(complex_stiff_rhs, [0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.5_dp, dp), cmplx(-0.25_dp, 0.75_dp, dp)], &
            method='RK45', jac_sparse_matrix=sparse_jac)
        call require(.not. sol%success, 'complex explicit-method CSC Jacobian rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'complex explicit-method CSC Jacobian status')

        sol = solve_ivp(complex_growth_rhs, [0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.25_dp, dp)], method='BDF', &
            jac=complex_growth_jacobian, jac_matrix=bad_matrix)
        call require(.not. sol%success, 'complex callback plus matrix rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'complex callback plus matrix status')

        sol = solve_ivp(complex_stiff_rhs, [0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.5_dp, dp), cmplx(-0.25_dp, 0.75_dp, dp)], &
            method='BDF', jac_matrix=jac_matrix, jac_sparse_matrix=sparse_jac)
        call require(.not. sol%success, 'complex dense plus CSC Jacobian rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'complex dense plus CSC Jacobian status')

        sol = solve_ivp(complex_stiff_rhs, [0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.5_dp, dp), cmplx(-0.25_dp, 0.75_dp, dp)], &
            method='RK45', jac_matrix=jac_matrix)
        call require(.not. sol%success, 'complex explicit-method matrix rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'complex explicit-method matrix status')

        sol = solve_ivp(complex_stiff_rhs, [0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.5_dp, dp), cmplx(-0.25_dp, 0.75_dp, dp)], &
            method='RK45', jac_sparsity=sparsity)
        call require(.not. sol%success, 'complex explicit-method sparsity rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'complex explicit-method sparsity status')


        sol = solve_ivp(complex_stiff_rhs, [0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.5_dp, dp), cmplx(-0.25_dp, 0.75_dp, dp)], &
            method='RK45', fun_vectorized=complex_stiff_rhs_vectorized)
        call require(.not. sol%success, 'complex explicit vectorized RHS rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'complex explicit vectorized RHS status')

        sol = solve_ivp(complex_stiff_rhs, [0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.5_dp, dp), cmplx(-0.25_dp, 0.75_dp, dp)], &
            method='RK45', jac_sparse=complex_stiff_sparse_jacobian)
        call require(.not. sol%success, 'complex explicit sparse callback rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'complex explicit sparse callback status')

        sol = solve_ivp(complex_stiff_rhs, [0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.5_dp, dp), cmplx(-0.25_dp, 0.75_dp, dp)], &
            method='BDF', jac=complex_growth_jacobian, &
            jac_sparse=complex_stiff_sparse_jacobian)
        call require(.not. sol%success, 'complex dense plus sparse callback rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'complex dense plus sparse callback status')

        sol = solve_ivp(complex_stiff_rhs, [0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.5_dp, dp), cmplx(-0.25_dp, 0.75_dp, dp)], &
            method='BDF', jac_sparse=invalid_complex_sparse_jacobian)
        call require(.not. sol%success, 'complex invalid sparse callback output rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'complex invalid initial sparse callback status')
    end subroutine test_complex_invalid_jacobian_options

    subroutine test_complex_rk23()
        type(solve_ivp_complex_result) :: sol
        complex(dp), allocatable :: values(:, :)
        real(dp), parameter :: teval(5) = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
        real(dp), parameter :: query(3) = [0.123_dp, 0.321_dp, 0.9_dp]
        complex(dp), parameter :: expected(5) = [ &
            cmplx(1.0_dp, 0.25_dp, dp), &
            cmplx(0.85861623704658852_dp, 0.79186814013494222_dp, dp), &
            cmplx(0.42364434114193883_dp, 1.2539106001269467_dp, dp), &
            cmplx(-0.25991466392814661_dp, 1.4770771235138185_dp, dp), &
            cmplx(-1.0609046591022375_dp, 1.3276504248971872_dp, dp)]
        complex(dp), parameter :: dense_expected(3) = [ &
            cmplx(0.9666718217075273_dp, 0.51682712920068985_dp, dp), &
            cmplx(0.76457263930564212_dp, 0.93813157032494177_dp, dp), &
            cmplx(-0.73814811176260964_dp, 1.4382161145570771_dp, dp)]
        integer :: i
        integer :: status

        sol = solve_ivp(complex_growth_rhs, [0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.25_dp, dp)], method='RK23', t_eval=teval, &
            dense_output=.true., rtol=1.0e-8_dp, atol=1.0e-10_dp)
        call require(sol%success, 'complex RK23 success')
        call require(sol%status == ivp_status_success, 'complex RK23 status')
        call require(sol%nfev == 881, 'complex RK23 SciPy nfev')
        do i = 1, size(teval)
            call check_complex_abs(sol%y(1, i), expected(i), 5.0e-14_dp, &
                'complex RK23 SciPy t_eval value')
        end do
        call sol%sol%evaluate(query, values, status)
        call require(status == ivp_status_success, 'complex RK23 dense status')
        do i = 1, size(query)
            call check_complex_abs(values(1, i), dense_expected(i), 5.0e-14_dp, &
                'complex RK23 SciPy dense value')
        end do
    end subroutine test_complex_rk23

    subroutine test_complex_rk45()
        type(solve_ivp_complex_result) :: sol
        complex(dp), allocatable :: values(:, :)
        real(dp), parameter :: teval(5) = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
        real(dp), parameter :: query(3) = [0.123_dp, 0.321_dp, 0.9_dp]
        complex(dp), parameter :: expected(5) = [ &
            cmplx(1.0_dp, 0.25_dp, dp), &
            cmplx(0.85861624564466799_dp, 0.79186813832332226_dp, dp), &
            cmplx(0.42364436020358398_dp, 1.2539106054744489_dp, dp), &
            cmplx(-0.2599146397541871_dp, 1.4770771468728949_dp, dp), &
            cmplx(-1.0609046439427536_dp, 1.3276504732778893_dp, dp)]
        complex(dp), parameter :: dense_expected(3) = [ &
            cmplx(0.96667182527709494_dp, 0.51682712736630476_dp, dp), &
            cmplx(0.76457265085327142_dp, 0.93813156927392161_dp, dp), &
            cmplx(-0.73814809079446297_dp, 1.4382161526249213_dp, dp)]
        integer :: i
        integer :: status

        sol = solve_ivp(complex_growth_rhs, [0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.25_dp, dp)], method='RK45', t_eval=teval, &
            dense_output=.true., rtol=1.0e-9_dp, atol=1.0e-12_dp)
        call require(sol%success, 'complex RK45 success')
        call require(sol%nfev == 218, 'complex RK45 SciPy nfev')
        do i = 1, size(teval)
            call check_complex_abs(sol%y(1, i), expected(i), 5.0e-13_dp, &
                'complex RK45 SciPy t_eval value')
        end do
        call sol%sol%evaluate(query, values, status)
        call require(status == ivp_status_success, 'complex RK45 dense status')
        do i = 1, size(query)
            call check_complex_abs(values(1, i), dense_expected(i), 5.0e-13_dp, &
                'complex RK45 SciPy dense value')
        end do
    end subroutine test_complex_rk45

    subroutine test_complex_dop853()
        type(solve_ivp_complex_result) :: sol
        complex(dp), allocatable :: values(:, :)
        real(dp), parameter :: teval(5) = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
        real(dp), parameter :: query(3) = [0.123_dp, 0.321_dp, 0.9_dp]
        complex(dp), parameter :: expected(5) = [ &
            cmplx(1.0_dp, 0.25_dp, dp), &
            cmplx(0.85861624558655703_dp, 0.79186813804148559_dp, dp), &
            cmplx(0.42364436045448556_dp, 1.2539106052554709_dp, dp), &
            cmplx(-0.25991463928696856_dp, 1.4770771469121384_dp, dp), &
            cmplx(-1.0609046434000386_dp, 1.3276504737107284_dp, dp)]
        complex(dp), parameter :: dense_expected(3) = [ &
            cmplx(0.96667182526535933_dp, 0.51682712729288549_dp, dp), &
            cmplx(0.76457265097827254_dp, 0.93813156916723517_dp, dp), &
            cmplx(-0.73814809024602746_dp, 1.4382161528960458_dp, dp)]
        integer :: i
        integer :: status

        sol = solve_ivp(complex_growth_rhs, [0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.25_dp, dp)], method='DOP853', t_eval=teval, &
            dense_output=.true., rtol=1.0e-11_dp, atol=1.0e-13_dp)
        call require(sol%success, 'complex DOP853 success')
        call require(sol%nfev == 152, 'complex DOP853 SciPy nfev')
        do i = 1, size(teval)
            call check_complex_abs(sol%y(1, i), expected(i), 5.0e-13_dp, &
                'complex DOP853 SciPy t_eval value')
        end do
        call sol%sol%evaluate(query, values, status)
        call require(status == ivp_status_success, 'complex DOP853 dense status')
        do i = 1, size(query)
            call check_complex_abs(values(1, i), dense_expected(i), 5.0e-13_dp, &
                'complex DOP853 SciPy dense value')
        end do
    end subroutine test_complex_dop853

    subroutine test_complex_backward_dop853()
        type(solve_ivp_complex_result) :: sol
        complex(dp), allocatable :: values(:, :)
        real(dp), parameter :: teval(5) = [1.0_dp, 0.75_dp, 0.5_dp, 0.25_dp, 0.0_dp]
        real(dp), parameter :: query(2) = [0.9_dp, 0.321_dp]
        complex(dp), parameter :: expected(5) = [ &
            cmplx(-1.0609046433999418_dp, 1.327650473712934_dp, dp), &
            cmplx(-0.25991463928759634_dp, 1.4770771469141746_dp, dp), &
            cmplx(0.42364436045553744_dp, 1.2539106052566062_dp, dp), &
            cmplx(0.85861624558968475_dp, 0.79186813803740674_dp, dp), &
            cmplx(1.0000000000011868_dp, 0.24999999999961986_dp, dp)]
        complex(dp), parameter :: dense_expected(2) = [ &
            cmplx(-0.73814809026211248_dp, 1.4382161528968418_dp, dp), &
            cmplx(0.76457265097118587_dp, 0.93813156918188123_dp, dp)]
        complex(dp) :: y1
        integer :: i
        integer :: status

        y1 = cmplx(-1.0609046433999418_dp, 1.327650473712934_dp, dp)
        sol = solve_ivp(complex_growth_rhs, [1.0_dp, 0.0_dp], [y1], method='DOP853', &
            t_eval=teval, dense_output=.true., rtol=1.0e-11_dp, atol=1.0e-13_dp)
        call require(sol%success, 'complex backward DOP853 success')
        call require(sol%nfev == 152, 'complex backward DOP853 SciPy nfev')
        do i = 1, size(teval)
            call check_complex_abs(sol%y(1, i), expected(i), 5.0e-12_dp, &
                'complex backward DOP853 SciPy t_eval value')
        end do
        call sol%sol%evaluate(query, values, status)
        call require(status == ivp_status_success, 'complex backward DOP853 dense status')
        do i = 1, size(query)
            call check_complex_abs(values(1, i), dense_expected(i), 5.0e-12_dp, &
                'complex backward DOP853 SciPy dense value')
        end do
    end subroutine test_complex_backward_dop853

    subroutine test_complex_bdf()
        type(solve_ivp_complex_result) :: sol
        complex(dp), allocatable :: values(:, :)
        real(dp), parameter :: teval(5) = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
        real(dp), parameter :: query(3) = [0.123_dp, 0.321_dp, 0.9_dp]
        complex(dp), parameter :: expected(5) = [ &
            cmplx(1.0_dp, 0.24999999999999997_dp, dp), &
            cmplx(0.858616245498677_dp, 0.7918681400295454_dp, dp), &
            cmplx(0.4236443538875859_dp, 1.253910608835228_dp, dp), &
            cmplx(-0.2599146539591201_dp, 1.4770771456291667_dp, dp), &
            cmplx(-1.0609046626980305_dp, 1.3276504598817922_dp, dp)]
        complex(dp), parameter :: dense_expected(3) = [ &
            cmplx(0.9666718265231606_dp, 0.5168271276844151_dp, dp), &
            cmplx(0.7645726494155605_dp, 0.9381315721018805_dp, dp), &
            cmplx(-0.7381481085501836_dp, 1.4382161449227011_dp, dp)]
        integer :: i
        integer :: status

        sol = solve_ivp(complex_growth_rhs, [0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.25_dp, dp)], method='BDF', t_eval=teval, &
            dense_output=.true., rtol=1.0e-9_dp, atol=1.0e-12_dp)
        call require(sol%success, 'complex BDF success')
        call require(sol%status == ivp_status_success, 'complex BDF status')
        call require(sol%nfev == 174, 'complex BDF SciPy nfev')
        call require(sol%njev == 1, 'complex BDF SciPy njev')
        call require(sol%nlu == 19, 'complex BDF SciPy nlu')
        do i = 1, size(teval)
            call check_complex_abs(sol%y(1, i), expected(i), 5.0e-13_dp, &
                'complex BDF SciPy t_eval value')
        end do
        call sol%sol%evaluate(query, values, status)
        call require(status == ivp_status_success, 'complex BDF dense status')
        do i = 1, size(query)
            call check_complex_abs(values(1, i), dense_expected(i), 5.0e-13_dp, &
                'complex BDF SciPy dense value')
        end do
    end subroutine test_complex_bdf

    subroutine test_complex_bdf_analytic_jacobian()
        type(solve_ivp_complex_result) :: sol

        sol = solve_ivp(complex_growth_rhs, [0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.25_dp, dp)], method='BDF', &
            rtol=1.0e-9_dp, atol=1.0e-12_dp, jac=complex_growth_jacobian)
        call require(sol%success, 'complex BDF analytic Jacobian success')
        call require(sol%nfev == 174, 'complex BDF analytic Jacobian nfev')
        call require(sol%njev == 1, 'complex BDF analytic Jacobian njev')
        call require(sol%nlu == 19, 'complex BDF analytic Jacobian nlu')
        call check_complex_abs(sol%y(1, size(sol%t)), &
            cmplx(-1.0609046626980305_dp, 1.3276504598817922_dp, dp), &
            5.0e-13_dp, 'complex BDF analytic Jacobian final state')
    end subroutine test_complex_bdf_analytic_jacobian

    subroutine test_complex_bdf_stiff_finite_difference()
        type(solve_ivp_complex_result) :: sol
        real(dp), parameter :: teval(5) = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
        complex(dp), parameter :: expected1(5) = [ &
            cmplx(1.0_dp, 0.5_dp, dp), &
            cmplx(1.9715491014134635e-11_dp, 3.792412591411449e-13_dp, dp), &
            cmplx(-1.873335927870132e-12_dp, 1.410491827063643e-12_dp, dp), &
            cmplx(3.7315212480989153e-13_dp, 1.1651951098161332e-13_dp, dp), &
            cmplx(-1.4230211839551403e-14_dp, -1.4265282999130226e-14_dp, dp)]
        complex(dp), parameter :: expected2(5) = [ &
            cmplx(-0.25_dp, 0.75_dp, dp), &
            cmplx(-0.4508982352706193_dp, 0.4192522435980384_dp, dp), &
            cmplx(-0.46471094083010706_dp, 0.11818795444607441_dp, dp), &
            cmplx(-0.36174094030829496_dp, -0.09273539247777629_dp, dp), &
            cmplx(-0.21261093166186007_dp, -0.19844686050431756_dp, dp)]
        integer :: i

        sol = solve_ivp(complex_stiff_rhs, [0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.5_dp, dp), cmplx(-0.25_dp, 0.75_dp, dp)], &
            method='BDF', t_eval=teval, dense_output=.true., &
            rtol=1.0e-8_dp, atol=1.0e-10_dp)
        call require(sol%success, 'complex stiff BDF success')
        call require(sol%nfev == 562, 'complex stiff BDF SciPy nfev')
        call require(sol%njev == 1, 'complex stiff BDF SciPy njev')
        call require(sol%nlu == 51, 'complex stiff BDF SciPy nlu')
        do i = 1, size(teval)
            call check_complex_abs(sol%y(1, i), expected1(i), 5.0e-12_dp, &
                'complex stiff BDF first component')
            call check_complex_abs(sol%y(2, i), expected2(i), 5.0e-12_dp, &
                'complex stiff BDF second component')
        end do
    end subroutine test_complex_bdf_stiff_finite_difference

    subroutine test_complex_bdf_coupled_finite_difference()
        type(solve_ivp_complex_result) :: sol
        complex(dp), allocatable :: values(:, :)
        real(dp), parameter :: teval(5) = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
        real(dp), parameter :: query(3) = [0.123_dp, 0.321_dp, 0.9_dp]
        complex(dp), parameter :: expected1(5) = [ &
            cmplx(1.0_dp, 0.2_dp, dp), &
            cmplx(-0.06852251703140068_dp, 0.08316785458136756_dp, dp), &
            cmplx(-0.05756695508976045_dp, 0.03810001061054903_dp, dp), &
            cmplx(-0.042917217731917465_dp, 0.016090932867193906_dp, dp), &
            cmplx(-0.03019161698322194_dp, 0.0040384926967332105_dp, dp)]
        complex(dp), parameter :: expected2(5) = [ &
            cmplx(-0.3_dp, 0.7_dp, dp), &
            cmplx(-0.46796983734424247_dp, 0.40321860443259727_dp, dp), &
            cmplx(-0.3636012482824385_dp, 0.19032076794640534_dp, dp), &
            cmplx(-0.26381531731015617_dp, 0.06917963731512479_dp, dp), &
            cmplx(-0.1811689970206235_dp, 0.005565420525653306_dp, dp)]
        complex(dp), parameter :: dense1(3) = [ &
            cmplx(0.00407535986090574_dp, 0.1605349513246475_dp, dp), &
            cmplx(-0.0677831762201405_dp, 0.06488045009765897_dp, dp), &
            cmplx(-0.03497859456348309_dp, 0.007942066011746964_dp, dp)]
        complex(dp), parameter :: dense2(3) = [ &
            cmplx(-0.4949707602976541_dp, 0.5658478028938988_dp, dp), &
            cmplx(-0.4396342333695293_dp, 0.33055086532734335_dp, dp), &
            cmplx(-0.2118761316059325_dp, 0.02576443796710964_dp, dp)]
        integer :: i
        integer :: status

        sol = solve_ivp(complex_coupled_rhs, [0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.2_dp, dp), cmplx(-0.3_dp, 0.7_dp, dp)], &
            method='BDF', t_eval=teval, dense_output=.true., &
            rtol=1.0e-9_dp, atol=1.0e-11_dp)
        call require(sol%success, 'complex coupled BDF success')
        call require(sol%nfev == 536, 'complex coupled BDF SciPy nfev')
        call require(sol%njev == 1, 'complex coupled BDF SciPy njev')
        call require(sol%nlu == 49, 'complex coupled BDF SciPy nlu')
        do i = 1, size(teval)
            call check_complex_abs(sol%y(1, i), expected1(i), 5.0e-12_dp, &
                'complex coupled BDF first component')
            call check_complex_abs(sol%y(2, i), expected2(i), 5.0e-12_dp, &
                'complex coupled BDF second component')
        end do
        call sol%sol%evaluate(query, values, status)
        call require(status == ivp_status_success, 'complex coupled BDF dense status')
        do i = 1, size(query)
            call check_complex_abs(values(1, i), dense1(i), 5.0e-12_dp, &
                'complex coupled BDF dense first component')
            call check_complex_abs(values(2, i), dense2(i), 5.0e-12_dp, &
                'complex coupled BDF dense second component')
        end do
    end subroutine test_complex_bdf_coupled_finite_difference

    subroutine test_complex_backward_bdf()
        type(solve_ivp_complex_result) :: sol
        complex(dp), allocatable :: values(:, :)
        real(dp), parameter :: teval(5) = [1.0_dp, 0.75_dp, 0.5_dp, 0.25_dp, 0.0_dp]
        real(dp), parameter :: query(2) = [0.9_dp, 0.321_dp]
        complex(dp), parameter :: expected(5) = [ &
            cmplx(-1.0609046433999418_dp, 1.327650473712934_dp, dp), &
            cmplx(-0.2599146408176843_dp, 1.477077143222394_dp, dp), &
            cmplx(0.4236443520811992_dp, 1.2539106039163095_dp, dp), &
            cmplx(0.8586162346372906_dp, 0.7918681435298038_dp, dp), &
            cmplx(0.9999999921285934_dp, 0.25000001274537564_dp, dp)]
        complex(dp), parameter :: dense_expected(2) = [ &
            cmplx(-0.7381480884436817_dp, 1.4382161503477335_dp, dp), &
            cmplx(0.7645726402206787_dp, 0.9381315724780575_dp, dp)]
        integer :: i
        integer :: status

        sol = solve_ivp(complex_growth_rhs, [1.0_dp, 0.0_dp], &
            [cmplx(-1.0609046433999418_dp, 1.327650473712934_dp, dp)], &
            method='BDF', t_eval=teval, dense_output=.true., &
            rtol=1.0e-9_dp, atol=1.0e-12_dp)
        call require(sol%success, 'complex backward BDF success')
        call require(sol%nfev == 176, 'complex backward BDF SciPy nfev')
        call require(sol%njev == 1, 'complex backward BDF SciPy njev')
        call require(sol%nlu == 20, 'complex backward BDF SciPy nlu')
        do i = 1, size(teval)
            call check_complex_abs(sol%y(1, i), expected(i), 5.0e-12_dp, &
                'complex backward BDF SciPy t_eval value')
        end do
        call sol%sol%evaluate(query, values, status)
        call require(status == ivp_status_success, 'complex backward BDF dense status')
        do i = 1, size(query)
            call check_complex_abs(values(1, i), dense_expected(i), 5.0e-12_dp, &
                'complex backward BDF SciPy dense value')
        end do
    end subroutine test_complex_backward_bdf

    subroutine test_complex_bdf_terminal_event()
        type(solve_ivp_complex_result) :: sol
        real(dp) :: pi

        pi = acos(-1.0_dp)
        sol = solve_ivp(rotation_rhs, [0.0_dp, 4.0_dp], &
            [cmplx(1.0_dp, 0.0_dp, dp)], method='BDF', events=real_zero_event, &
            event_terminal=.true., event_direction=-1.0_dp, &
            rtol=1.0e-10_dp, atol=1.0e-12_dp, dense_output=.true.)
        call require(sol%success, 'complex BDF terminal event success')
        call require(sol%status == ivp_status_event, 'complex BDF terminal event status')
        call require(sol%njev == 1, 'complex BDF terminal event njev')
        call require(sol%nlu == 20, 'complex BDF terminal event nlu')
        call require(size(sol%t_events) == 1, 'complex BDF terminal event count')
        call check_abs(sol%t_events(1), 0.5_dp * pi, 3.0e-10_dp, &
            'complex BDF terminal event time')
        call check_complex_abs(sol%y_events(1, 1), cmplx(0.0_dp, 1.0_dp, dp), &
            3.0e-9_dp, 'complex BDF terminal event state')
    end subroutine test_complex_bdf_terminal_event

    subroutine test_complex_terminal_event()
        type(solve_ivp_complex_result) :: sol
        real(dp) :: pi

        pi = acos(-1.0_dp)
        sol = solve_ivp(rotation_rhs, [0.0_dp, 4.0_dp], &
            [cmplx(1.0_dp, 0.0_dp, dp)], method='RK45', events=real_zero_event, &
            event_terminal=.true., event_direction=-1.0_dp, &
            rtol=1.0e-10_dp, atol=1.0e-12_dp, dense_output=.true.)
        call require(sol%success, 'complex terminal event success')
        call require(sol%status == ivp_status_event, 'complex terminal event status')
        call require(size(sol%t_events) == 1, 'complex terminal event count')
        call check_abs(sol%t_events(1), 0.5_dp * pi, 1.0e-10_dp, &
            'complex terminal event time')
        call check_complex_abs(sol%y_events(1, 1), cmplx(0.0_dp, 1.0_dp, dp), &
            5.0e-10_dp, 'complex terminal event state')
    end subroutine test_complex_terminal_event

    subroutine test_complex_event_occurrence_limit()
        type(solve_ivp_complex_result) :: sol
        real(dp) :: pi

        pi = acos(-1.0_dp)
        sol = solve_ivp(rotation_rhs, [0.0_dp, 5.0_dp], &
            [cmplx(1.0_dp, 0.0_dp, dp)], method='DOP853', events=real_zero_event, &
            event_max_occurrences=2, rtol=1.0e-10_dp, atol=1.0e-12_dp)
        call require(sol%success, 'complex event occurrence success')
        call require(sol%status == ivp_status_event, 'complex event occurrence status')
        call require(size(sol%t_events) == 2, 'complex event occurrence count')
        call require(sol%event_counts(1) == 2, 'complex event occurrence grouped count')
        call check_abs(sol%t_events(1), 0.5_dp * pi, 2.0e-10_dp, &
            'complex first event time')
        call check_abs(sol%t_events(2), 1.5_dp * pi, 2.0e-10_dp, &
            'complex second event time')
    end subroutine test_complex_event_occurrence_limit

    subroutine test_complex_vector_events()
        type(solve_ivp_complex_result) :: sol
        type(ivp_complex_event_set) :: events
        real(dp) :: pi

        pi = acos(-1.0_dp)
        events%n_events = 2
        events%fun => rotation_events
        allocate(events%terminal(2), events%direction(2))
        events%terminal = [.false., .true.]
        events%direction = [1.0_dp, -1.0_dp]
        sol = solve_ivp(rotation_rhs, [0.0_dp, 2.0_dp], &
            [cmplx(1.0_dp, 0.0_dp, dp)], method='RK45', event_set=events, &
            rtol=1.0e-10_dp, atol=1.0e-12_dp)
        call require(sol%success, 'complex vector events success')
        call require(sol%status == ivp_status_event, 'complex vector events status')
        call require(size(sol%t_events) == 2, 'complex vector events count')
        call require(all(sol%event_indices == [1, 2]), 'complex vector event indices')
        call check_abs(sol%t_events(1), pi / 6.0_dp, 2.0e-10_dp, &
            'complex vector first event time')
        call check_abs(sol%t_events(2), pi / 3.0_dp, 2.0e-10_dp, &
            'complex vector second event time')
        call check_complex_abs(sol%y_events(1, 1), &
            cmplx(sqrt(3.0_dp) / 2.0_dp, 0.5_dp, dp), 5.0e-10_dp, &
            'complex vector first event state')
        call check_complex_abs(sol%y_events(1, 2), &
            cmplx(0.5_dp, sqrt(3.0_dp) / 2.0_dp, dp), 5.0e-10_dp, &
            'complex vector second event state')
    end subroutine test_complex_vector_events

    subroutine test_complex_invalid_method()
        type(solve_ivp_complex_result) :: sol

        sol = solve_ivp(rotation_rhs, [0.0_dp, 1.0_dp], &
            [cmplx(1.0_dp, 0.0_dp, dp)], method='Radau')
        call require(.not. sol%success, 'complex Radau rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'complex invalid method status')
    end subroutine test_complex_invalid_method

    subroutine complex_multiscale_decay_rhs(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable, unused by autonomous system
        complex(dp), intent(in) :: y(:) !! two-component complex state
        complex(dp), intent(out) :: dydt(:) !! derivative with separated decay rates

        if (t /= t .or. size(y) /= 2) error stop 'unexpected complex multiscale input'
        dydt(1) = cmplx(-1.0_dp, 1.0_dp, dp) * y(1)
        dydt(2) = cmplx(-10.0_dp, 2.0_dp, dp) * y(2)
    end subroutine complex_multiscale_decay_rhs

    subroutine complex_growth_rhs(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable
        complex(dp), intent(in) :: y(:) !! complex state
        complex(dp), intent(out) :: dydt(:) !! complex state derivative

        dydt = cmplx(0.5_dp, 2.0_dp, dp) * y
    end subroutine complex_growth_rhs

    subroutine complex_growth_jacobian(t, y, jac)
        real(dp), intent(in) :: t !! independent variable
        complex(dp), intent(in) :: y(:) !! complex state at t
        complex(dp), intent(out) :: jac(:, :) !! complex Jacobian d fun / d y

        jac = cmplx(0.0_dp, 0.0_dp, dp)
        jac(1, 1) = cmplx(0.5_dp, 2.0_dp, dp)
    end subroutine complex_growth_jacobian

    subroutine complex_stiff_rhs(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable
        complex(dp), intent(in) :: y(:) !! complex state
        complex(dp), intent(out) :: dydt(:) !! complex state derivative

        complex_stiff_rhs_calls = complex_stiff_rhs_calls + 1
        dydt(1) = cmplx(-100.0_dp, 30.0_dp, dp) * y(1)
        dydt(2) = cmplx(-1.0_dp, 2.0_dp, dp) * y(2)
    end subroutine complex_stiff_rhs

    subroutine complex_stiff_sparse_jacobian(t, y, jac)
        real(dp), intent(in) :: t !! independent variable
        complex(dp), intent(in) :: y(:) !! complex stiff state
        type(ivp_complex_sparse_matrix), intent(out) :: jac !! one-based CSC Jacobian

        complex(dp) :: data(2)

        if (size(y) /= 2) error stop 'unexpected complex sparse Jacobian state'
        complex_sparse_jacobian_calls = complex_sparse_jacobian_calls + 1
        data(1) = cmplx(-100.0_dp, 30.0_dp, dp) + cmplx(0.0_dp * t, 0.0_dp, dp) * y(1)
        data(2) = cmplx(-1.0_dp, 2.0_dp, dp) + cmplx(0.0_dp * t, 0.0_dp, dp) * y(2)
        jac = ivp_csc_matrix(2, 2, [1, 2, 3], [1, 2], data)
    end subroutine complex_stiff_sparse_jacobian

    subroutine invalid_complex_sparse_jacobian(t, y, jac)
        real(dp), intent(in) :: t !! independent variable
        complex(dp), intent(in) :: y(:) !! complex state
        type(ivp_complex_sparse_matrix), intent(out) :: jac !! deliberately malformed CSC matrix

        if (t /= t .or. size(y) < 1) error stop 'unexpected invalid complex sparse input'
        jac = ivp_csc_matrix(2, 2, [1, 3, 2], [1, 2], &
            [cmplx(-100.0_dp, 30.0_dp, dp), cmplx(-1.0_dp, 2.0_dp, dp)])
    end subroutine invalid_complex_sparse_jacobian

    subroutine complex_stiff_rhs_vectorized(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable shared by all state columns
        complex(dp), intent(in) :: y(:, :) !! complex stiff state columns
        complex(dp), intent(out) :: dydt(:, :) !! derivatives for each state column

        integer :: k

        if (t /= t) error stop 'unexpected complex vectorized time'
        complex_vectorized_calls = complex_vectorized_calls + 1
        complex_vectorized_max_columns = max(complex_vectorized_max_columns, size(y, 2))
        do k = 1, size(y, 2)
            dydt(1, k) = cmplx(-100.0_dp, 30.0_dp, dp) * y(1, k)
            dydt(2, k) = cmplx(-1.0_dp, 2.0_dp, dp) * y(2, k)
        end do
    end subroutine complex_stiff_rhs_vectorized

    subroutine complex_coupled_rhs(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable
        complex(dp), intent(in) :: y(:) !! complex state
        complex(dp), intent(out) :: dydt(:) !! complex state derivative

        dydt(1) = cmplx(-20.0_dp, 5.0_dp, dp) * y(1) + &
            cmplx(3.0_dp, -1.0_dp, dp) * y(2)
        dydt(2) = cmplx(-4.0_dp, 2.0_dp, dp) * y(1) + &
            cmplx(-1.0_dp, 0.5_dp, dp) * y(2)
    end subroutine complex_coupled_rhs

    subroutine rotation_rhs(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable
        complex(dp), intent(in) :: y(:) !! complex state
        complex(dp), intent(out) :: dydt(:) !! complex state derivative

        dydt = cmplx(0.0_dp, 1.0_dp, dp) * y
    end subroutine rotation_rhs

    subroutine rotation_events(t, y, values)
        real(dp), intent(in) :: t !! independent variable
        complex(dp), intent(in) :: y(:) !! complex state
        real(dp), intent(out) :: values(:) !! real event-function values

        values(1) = aimag(y(1)) - 0.5_dp
        values(2) = real(y(1), dp) - 0.5_dp
    end subroutine rotation_events

    function real_zero_event(t, y) result(value)
        real(dp), intent(in) :: t !! independent variable
        complex(dp), intent(in) :: y(:) !! complex state
        real(dp) :: value

        value = real(y(1), dp)
    end function real_zero_event

    subroutine require(condition, label)
        logical, intent(in) :: condition !! assertion condition
        character(len=*), intent(in) :: label !! assertion label

        if (.not. condition) then
            print *, 'FAIL:', trim(label)
            error stop 1
        end if
    end subroutine require

    subroutine check_abs(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual !! computed scalar
        real(dp), intent(in) :: expected !! reference scalar
        real(dp), intent(in) :: tolerance !! absolute tolerance
        character(len=*), intent(in) :: label !! assertion label

        if (abs(actual - expected) > tolerance) then
            print *, 'FAIL:', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance= ', tolerance
            error stop 1
        end if
    end subroutine check_abs

    subroutine check_complex_abs(actual, expected, tolerance, label)
        complex(dp), intent(in) :: actual !! computed complex scalar
        complex(dp), intent(in) :: expected !! reference complex scalar
        real(dp), intent(in) :: tolerance !! absolute complex tolerance
        character(len=*), intent(in) :: label !! assertion label

        if (abs(actual - expected) > tolerance) then
            print *, 'FAIL:', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            print *, ' tolerance= ', tolerance
            error stop 1
        end if
    end subroutine check_complex_abs

end program test_integrate_ivp_complex
