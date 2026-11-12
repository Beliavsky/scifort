program test_integrate_ivp
    use scifort_integrate, only : dp, ivp_csc_matrix, ivp_event_set, ivp_sparse_matrix, &
        ivp_status_event, ivp_status_invalid_argument, ivp_status_success, solve_ivp, &
        solve_ivp_result
    use scifort_integrate_bdf, only : bdf_initialize, bdf_state
    use scifort_integrate_lsoda, only : lsoda_advance, lsoda_initialize, lsoda_state, &
        lsoda_status_success
    use scifort_integrate_radau, only : radau_initialize, radau_state
    implicit none

    integer :: diagonal_rhs_calls = 0
    integer :: diagonal_vectorized_calls = 0
    integer :: diagonal_vectorized_max_columns = 0
    integer :: sparse_jacobian_calls = 0

    call test_rk45_exponential()
    call test_vector_atol()
    call test_invalid_vector_atol()
    call test_vector_rtol()
    call test_invalid_vector_rtol()
    call test_constant_jacobian_matrix()
    call test_sparse_jacobian_matrix()
    call test_sparse_jacobian_callback()
    call test_sparse_numerical_storage()
    call test_jacobian_sparsity()
    call test_vectorized_rhs_jacobian()
    call test_invalid_jacobian_options()
    call test_rk45_scipy_dense_output()
    call test_rk23_harmonic()
    call test_dop853_exponential()
    call test_dop853_harmonic_t_eval()
    call test_dop853_vanderpol()
    call test_dop853_native_dense()
    call test_dop853_terminal_event()
    call test_radau_exponential_dense()
    call test_radau_stiff_with_jacobian()
    call test_radau_multicomponent_fd_jacobian()
    call test_radau_backward_dense()
    call test_radau_terminal_event()
    call test_bdf_exponential_dense()
    call test_bdf_stiff_with_jacobian()
    call test_bdf_multicomponent_fd_jacobian()
    call test_bdf_backward_dense()
    call test_bdf_terminal_event()
    call test_lsoda_exponential_dense()
    call test_lsoda_low_order_startup()
    call test_lsoda_startup_convergence_failure()
    call test_lsoda_high_order_adams()
    call test_lsoda_order12_nordsieck()
    call test_lsoda_variable_step_history()
    call test_lsoda_nordsieck_rescale()
    call test_lsoda_backward_nordsieck()
    call test_lsoda_error_failure_recovery()
    call test_lsoda_roundoff_switch_branch()
    call test_lsoda_switch_hysteresis()
    call test_lsoda_stiff_switch()
    call test_lsoda_banded_fd()
    call test_lsoda_banded_jacobian()
    call test_lsoda_tridiagonal_bands()
    call test_lsoda_min_step()
    call test_lsoda_terminal_event()
    call test_multiple_events()
    call test_event_occurrence_limit()
    call test_terminal_event_root_order()
    call test_event_direction_filter()
    call test_initial_event_direction_filter()
    call test_backward_and_t_eval()
    call test_dense_output()
    call test_terminal_event()
    call test_invalid_method()

    print *, 'test_integrate_ivp: PASS'

contains

    subroutine test_vector_atol()
        character(len=6), parameter :: methods(6) = [character(len=6) :: &
            'RK23', 'RK45', 'DOP853', 'Radau', 'BDF', 'LSODA']
        type(solve_ivp_result) :: relaxed
        type(solve_ivp_result) :: tight
        real(dp) :: expected
        integer :: i

        expected = exp(-2.0_dp)
        do i = 1, size(methods)
            relaxed = solve_ivp(multiscale_decay_rhs, [0.0_dp, 2.0_dp], &
                [1.0_dp, 1.0_dp], method=trim(methods(i)), rtol=1.0e-6_dp, &
                atol=[1.0e-9_dp, 1.0e-3_dp])
            tight = solve_ivp(multiscale_decay_rhs, [0.0_dp, 2.0_dp], &
                [1.0_dp, 1.0_dp], method=trim(methods(i)), rtol=1.0e-6_dp, &
                atol=1.0e-9_dp)

            call require(relaxed%success, trim(methods(i))//' vector atol success')
            call require(tight%success, trim(methods(i))//' scalar atol success')
            call check_abs(relaxed%y(1, size(relaxed%t)), expected, 5.0e-6_dp, &
                trim(methods(i))//' vector atol first component')
            call require(abs(relaxed%y(2, size(relaxed%t))) <= 1.0e-3_dp, &
                trim(methods(i))//' vector atol second component')
            call require(relaxed%nfev < tight%nfev, &
                trim(methods(i))//' vector atol reduces RHS work')
        end do
    end subroutine test_vector_atol

    subroutine test_invalid_vector_atol()
        type(solve_ivp_result) :: sol
        real(dp) :: rank_two_atol(1, 2)

        sol = solve_ivp(multiscale_decay_rhs, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], &
            atol=[1.0e-6_dp])
        call require(.not. sol%success, 'short vector atol rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'short vector atol status')

        rank_two_atol = 1.0e-6_dp
        sol = solve_ivp(multiscale_decay_rhs, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], &
            atol=rank_two_atol)
        call require(.not. sol%success, 'rank-two atol rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'rank-two atol status')

        sol = solve_ivp(multiscale_decay_rhs, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], &
            atol=[1.0e-6_dp, -1.0e-6_dp])
        call require(.not. sol%success, 'negative vector atol rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'negative vector atol status')
    end subroutine test_invalid_vector_atol

    subroutine test_vector_rtol()
        character(len=6), parameter :: methods(6) = [character(len=6) :: &
            'RK23', 'RK45', 'DOP853', 'Radau', 'BDF', 'LSODA']
        type(solve_ivp_result) :: scalar_tol
        type(solve_ivp_result) :: vector_tol
        real(dp) :: expected_fast
        real(dp) :: expected_slow
        integer :: i

        expected_slow = exp(-0.2_dp)
        expected_fast = exp(-2.0_dp)
        do i = 1, size(methods)
            scalar_tol = solve_ivp(multiscale_decay_rhs, [0.0_dp, 0.2_dp], &
                [1.0_dp, 1.0_dp], method=trim(methods(i)), rtol=1.0e-3_dp, &
                atol=1.0e-12_dp)
            vector_tol = solve_ivp(multiscale_decay_rhs, [0.0_dp, 0.2_dp], &
                [1.0_dp, 1.0_dp], method=trim(methods(i)), &
                rtol=[1.0e-3_dp, 1.0e-8_dp], atol=1.0e-12_dp)

            call require(scalar_tol%success, trim(methods(i))//' scalar rtol success')
            call require(vector_tol%success, trim(methods(i))//' vector rtol success')
            call check_abs(vector_tol%y(1, size(vector_tol%t)), expected_slow, 5.0e-4_dp, &
                trim(methods(i))//' vector rtol slow component')
            call check_abs(vector_tol%y(2, size(vector_tol%t)), expected_fast, 2.0e-8_dp, &
                trim(methods(i))//' vector rtol fast component')
            call require(vector_tol%nfev > scalar_tol%nfev, &
                trim(methods(i))//' tighter vector rtol increases RHS work')
        end do
    end subroutine test_vector_rtol

    subroutine test_invalid_vector_rtol()
        type(solve_ivp_result) :: sol
        real(dp) :: rank_two_rtol(1, 2)

        sol = solve_ivp(multiscale_decay_rhs, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], &
            rtol=[1.0e-3_dp])
        call require(.not. sol%success, 'short vector rtol rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'short vector rtol status')

        rank_two_rtol = 1.0e-3_dp
        sol = solve_ivp(multiscale_decay_rhs, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], &
            rtol=rank_two_rtol)
        call require(.not. sol%success, 'rank-two rtol rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'rank-two rtol status')

        sol = solve_ivp(multiscale_decay_rhs, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], &
            rtol=[1.0e-3_dp, 0.0_dp])
        call require(.not. sol%success, 'nonpositive vector rtol rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'nonpositive vector rtol status')
    end subroutine test_invalid_vector_rtol

    subroutine test_constant_jacobian_matrix()
        character(len=5), parameter :: methods(2) = [character(len=5) :: 'BDF', 'Radau']
        type(solve_ivp_result) :: sol
        real(dp) :: jac_matrix(8, 8)
        real(dp) :: y0(8)
        integer :: i
        integer :: j

        y0 = 1.0_dp
        jac_matrix = 0.0_dp
        do i = 1, size(y0)
            jac_matrix(i, i) = -(400.0_dp + 100.0_dp * real(i, dp))
        end do

        do j = 1, size(methods)
            diagonal_rhs_calls = 0
            sol = solve_ivp(diagonal_stiff_rhs, [0.0_dp, 0.5_dp], y0, &
                method=trim(methods(j)), jac_matrix=jac_matrix, &
                rtol=2.0e-6_dp, atol=1.0e-9_dp)
            call require(sol%success, trim(methods(j))//' constant Jacobian success')
            call require(sol%njev == 0, trim(methods(j))//' constant Jacobian njev')
            call require(sol%nlu > 0, trim(methods(j))//' constant Jacobian LU work')
            call require(sol%nfev == diagonal_rhs_calls, &
                trim(methods(j))//' constant Jacobian has no hidden RHS probes')
            do i = 1, size(y0)
                call check_abs(sol%y(i, size(sol%t)), cos(0.5_dp), 3.0e-6_dp, &
                    trim(methods(j))//' constant Jacobian value')
            end do
        end do
    end subroutine test_constant_jacobian_matrix

    subroutine test_sparse_jacobian_matrix()
        character(len=5), parameter :: methods(2) = [character(len=5) :: 'BDF', 'Radau']
        type(ivp_sparse_matrix) :: jac_sparse
        type(solve_ivp_result) :: sol
        integer :: indices(8)
        integer :: indptr(9)
        real(dp) :: data(8)
        real(dp) :: y0(8)
        integer :: i
        integer :: j

        y0 = 1.0_dp
        do i = 1, size(y0)
            indptr(i) = i
            indices(i) = i
            data(i) = -(400.0_dp + 100.0_dp * real(i, dp))
        end do
        indptr(size(indptr)) = size(y0) + 1
        jac_sparse = ivp_csc_matrix(size(y0), size(y0), indptr, indices, data)
        call require(jac_sparse%is_valid(), 'real CSC Jacobian valid')
        call require(jac_sparse%nnz() == size(y0), 'real CSC Jacobian nnz')

        do j = 1, size(methods)
            diagonal_rhs_calls = 0
            sol = solve_ivp(diagonal_stiff_rhs, [0.0_dp, 0.5_dp], y0, &
                method=trim(methods(j)), jac_sparse_matrix=jac_sparse, &
                rtol=2.0e-6_dp, atol=1.0e-9_dp)
            call require(sol%success, trim(methods(j))//' sparse matrix Jacobian success')
            call require(sol%njev == 0, trim(methods(j))//' sparse matrix Jacobian njev')
            call require(sol%nlu > 0, trim(methods(j))//' sparse matrix Jacobian LU work')
            call require(sol%nfev == diagonal_rhs_calls, &
                trim(methods(j))//' sparse matrix Jacobian has no hidden RHS probes')
            do i = 1, size(y0)
                call check_abs(sol%y(i, size(sol%t)), cos(0.5_dp), 3.0e-6_dp, &
                    trim(methods(j))//' sparse matrix Jacobian value')
            end do
        end do
    end subroutine test_sparse_jacobian_matrix

    subroutine test_sparse_jacobian_callback()
        character(len=5), parameter :: methods(2) = [character(len=5) :: 'BDF', 'Radau']
        type(solve_ivp_result) :: sol
        real(dp) :: y0(8)
        integer :: i
        integer :: j

        y0 = 1.0_dp
        do j = 1, size(methods)
            diagonal_rhs_calls = 0
            sparse_jacobian_calls = 0
            sol = solve_ivp(diagonal_stiff_rhs, [0.0_dp, 0.5_dp], y0, &
                method=trim(methods(j)), jac_sparse=diagonal_stiff_sparse_jacobian, &
                rtol=2.0e-6_dp, atol=1.0e-9_dp)
            call require(sol%success, trim(methods(j))//' sparse callback Jacobian success')
            call require(sol%njev == sparse_jacobian_calls, &
                trim(methods(j))//' sparse callback Jacobian count')
            call require(sol%njev > 0, trim(methods(j))//' sparse callback Jacobian used')
            call require(sol%nlu > 0, trim(methods(j))//' sparse callback Jacobian LU work')
            call require(sol%nfev == diagonal_rhs_calls, &
                trim(methods(j))//' sparse callback has no hidden RHS probes')
            do i = 1, size(y0)
                call check_abs(sol%y(i, size(sol%t)), cos(0.5_dp), 3.0e-6_dp, &
                    trim(methods(j))//' sparse callback Jacobian value')
            end do
        end do
    end subroutine test_sparse_jacobian_callback

    subroutine test_sparse_numerical_storage()
        type(bdf_state) :: bdf
        type(radau_state) :: radau
        logical :: sparsity(8, 8)
        real(dp) :: f(8)
        real(dp) :: y0(8)
        integer :: i
        integer :: njev

        y0 = 1.0_dp
        sparsity = .false.
        do i = 1, size(y0)
            sparsity(i, i) = .true.
        end do
        call diagonal_stiff_rhs(0.0_dp, y0, f)

        njev = 0
        call bdf_initialize(bdf, diagonal_stiff_rhs, 0.0_dp, y0, f, 1.0_dp, &
            1.0e-3_dp, 2.0e-6_dp, 1.0e-9_dp, njev, jac_sparsity=sparsity)
        call require(bdf%jacobian_sparse, 'BDF sparsity uses sparse Jacobian state')
        call require(.not. allocated(bdf%jac), 'BDF sparsity avoids dense Jacobian storage')
        call require(size(bdf%lu) == 0, 'BDF sparsity avoids dense LU workspace')
        call require(bdf%sparse_jac%is_valid(), 'BDF numerical CSC Jacobian valid')
        call require(bdf%sparse_jac%nnz() == size(y0), 'BDF numerical CSC Jacobian nnz')

        njev = 0
        call radau_initialize(radau, diagonal_stiff_rhs, 0.0_dp, y0, f, &
            1.0e-3_dp, 2.0e-6_dp, 1.0e-9_dp, njev, jac_sparsity=sparsity)
        call require(radau%jacobian_sparse, 'Radau sparsity uses sparse Jacobian state')
        call require(.not. allocated(radau%jac), 'Radau sparsity avoids dense Jacobian storage')
        call require(size(radau%lu_real) == 0 .and. size(radau%lu_complex) == 0, &
            'Radau sparsity avoids dense LU workspaces')
        call require(radau%sparse_jac%is_valid(), 'Radau numerical CSC Jacobian valid')
        call require(radau%sparse_jac%nnz() == size(y0), 'Radau numerical CSC Jacobian nnz')
    end subroutine test_sparse_numerical_storage

    subroutine test_jacobian_sparsity()
        character(len=5), parameter :: methods(2) = [character(len=5) :: 'BDF', 'Radau']
        type(solve_ivp_result) :: dense_sol
        type(solve_ivp_result) :: sparse_sol
        logical :: sparsity(8, 8)
        real(dp) :: y0(8)
        integer :: dense_calls
        integer :: i
        integer :: j
        integer :: sparse_calls

        y0 = 1.0_dp
        sparsity = .false.
        do i = 1, size(y0)
            sparsity(i, i) = .true.
        end do

        do j = 1, size(methods)
            diagonal_rhs_calls = 0
            dense_sol = solve_ivp(diagonal_stiff_rhs, [0.0_dp, 0.5_dp], y0, &
                method=trim(methods(j)), rtol=2.0e-6_dp, atol=1.0e-9_dp)
            dense_calls = diagonal_rhs_calls
            call require(dense_sol%success, trim(methods(j))//' dense Jacobian success')

            diagonal_rhs_calls = 0
            sparse_sol = solve_ivp(diagonal_stiff_rhs, [0.0_dp, 0.5_dp], y0, &
                method=trim(methods(j)), jac_sparsity=sparsity, &
                rtol=2.0e-6_dp, atol=1.0e-9_dp)
            sparse_calls = diagonal_rhs_calls
            call require(sparse_sol%success, trim(methods(j))//' sparse Jacobian success')
            call require(sparse_sol%nfev == dense_sol%nfev, &
                trim(methods(j))//' sparsity preserves reported nfev')
            call require(sparse_sol%njev == dense_sol%njev, &
                trim(methods(j))//' sparsity preserves Jacobian count')
            call require(sparse_calls < dense_calls, &
                trim(methods(j))//' sparsity reduces RHS probes')
            do i = 1, size(y0)
                call check_abs(sparse_sol%y(i, size(sparse_sol%t)), cos(0.5_dp), &
                    3.0e-6_dp, trim(methods(j))//' sparse Jacobian value')
            end do
        end do
    end subroutine test_jacobian_sparsity

    subroutine test_vectorized_rhs_jacobian()
        character(len=5), parameter :: methods(2) = [character(len=5) :: 'BDF', 'Radau']
        type(solve_ivp_result) :: scalar_sol
        type(solve_ivp_result) :: sparse_vec_sol
        type(solve_ivp_result) :: vec_sol
        logical :: sparsity(8, 8)
        real(dp) :: y0(8)
        integer :: i
        integer :: j
        integer :: scalar_calls

        y0 = 1.0_dp
        sparsity = .false.
        do i = 1, size(y0)
            sparsity(i, i) = .true.
        end do

        do j = 1, size(methods)
            diagonal_rhs_calls = 0
            scalar_sol = solve_ivp(diagonal_stiff_rhs, [0.0_dp, 0.5_dp], y0, &
                method=trim(methods(j)), rtol=2.0e-6_dp, atol=1.0e-9_dp)
            scalar_calls = diagonal_rhs_calls
            call require(scalar_sol%success, trim(methods(j))//' scalar FD success')

            diagonal_rhs_calls = 0
            diagonal_vectorized_calls = 0
            diagonal_vectorized_max_columns = 0
            vec_sol = solve_ivp(diagonal_stiff_rhs, [0.0_dp, 0.5_dp], y0, &
                method=trim(methods(j)), rtol=2.0e-6_dp, atol=1.0e-9_dp, &
                fun_vectorized=diagonal_stiff_rhs_vectorized)
            call require(vec_sol%success, trim(methods(j))//' vectorized FD success')
            call require(vec_sol%nfev == scalar_sol%nfev, &
                trim(methods(j))//' vectorized FD preserves nfev')
            call require(vec_sol%njev == scalar_sol%njev, &
                trim(methods(j))//' vectorized FD preserves njev')
            call require(diagonal_rhs_calls == vec_sol%nfev, &
                trim(methods(j))//' vectorized FD keeps hidden probes batched')
            call require(diagonal_vectorized_calls > 0, &
                trim(methods(j))//' vectorized FD callback used')
            call require(diagonal_vectorized_max_columns == size(y0), &
                trim(methods(j))//' dense FD batches all perturbed columns')
            call require(diagonal_rhs_calls + diagonal_vectorized_calls < scalar_calls, &
                trim(methods(j))//' vectorized FD reduces callback invocations')

            diagonal_rhs_calls = 0
            diagonal_vectorized_calls = 0
            diagonal_vectorized_max_columns = 0
            sparse_vec_sol = solve_ivp(diagonal_stiff_rhs, [0.0_dp, 0.5_dp], y0, &
                method=trim(methods(j)), rtol=2.0e-6_dp, atol=1.0e-9_dp, &
                jac_sparsity=sparsity, fun_vectorized=diagonal_stiff_rhs_vectorized)
            call require(sparse_vec_sol%success, &
                trim(methods(j))//' vectorized sparse FD success')
            call require(diagonal_vectorized_max_columns == 1, &
                trim(methods(j))//' sparse grouping controls vectorized batch width')
            do i = 1, size(y0)
                call check_abs(sparse_vec_sol%y(i, size(sparse_vec_sol%t)), cos(0.5_dp), &
                    3.0e-6_dp, trim(methods(j))//' vectorized sparse FD value')
            end do
        end do
    end subroutine test_vectorized_rhs_jacobian

    subroutine test_invalid_jacobian_options()
        type(ivp_sparse_matrix) :: bad_sparse
        type(ivp_sparse_matrix) :: sparse_jac
        type(solve_ivp_result) :: sol
        logical :: bad_sparsity(1, 1)
        logical :: sparsity(2, 2)
        real(dp) :: bad_matrix(1, 1)
        real(dp) :: jac_matrix(2, 2)

        bad_matrix = -1.0_dp
        bad_sparsity = .true.
        jac_matrix = 0.0_dp
        jac_matrix(1, 1) = -1000.0_dp
        jac_matrix(2, 2) = -1.0_dp
        sparsity = .false.
        sparsity(1, 1) = .true.
        sparsity(2, 2) = .true.
        sparse_jac = ivp_csc_matrix(2, 2, [1, 2, 3], [1, 2], &
            [-1000.0_dp, -1.0_dp])
        bad_sparse = ivp_csc_matrix(2, 2, [1, 3, 2], [1, 2], &
            [-1000.0_dp, -1.0_dp])

        sol = solve_ivp(stiff_linear_rhs, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], &
            method='BDF', jac_matrix=bad_matrix)
        call require(.not. sol%success, 'wrong-shape constant Jacobian rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'wrong-shape constant Jacobian status')

        sol = solve_ivp(stiff_linear_rhs, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], &
            method='Radau', jac_sparsity=bad_sparsity)
        call require(.not. sol%success, 'wrong-shape Jacobian sparsity rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'wrong-shape Jacobian sparsity status')

        sol = solve_ivp(stiff_linear_rhs, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], &
            method='BDF', jac_sparse_matrix=bad_sparse)
        call require(.not. sol%success, 'invalid CSC Jacobian rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'invalid CSC Jacobian status')

        sol = solve_ivp(stiff_linear_rhs, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], &
            method='RK45', jac_sparse_matrix=sparse_jac)
        call require(.not. sol%success, 'explicit-method CSC Jacobian rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'explicit-method CSC Jacobian status')

        sol = solve_ivp(stiff_linear_rhs, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], &
            method='LSODA', jac_matrix=jac_matrix)
        call require(.not. sol%success, 'LSODA constant Jacobian rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'LSODA constant Jacobian status')

        sol = solve_ivp(stiff_linear_rhs, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], &
            method='LSODA', jac_sparse_matrix=sparse_jac)
        call require(.not. sol%success, 'LSODA CSC Jacobian rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'LSODA CSC Jacobian status')

        sol = solve_ivp(stiff_linear_rhs, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], &
            method='RK45', jac_sparsity=sparsity)
        call require(.not. sol%success, 'explicit-method Jacobian sparsity rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'explicit-method Jacobian sparsity status')

        sol = solve_ivp(stiff_cosine_rhs, [0.0_dp, 1.0_dp], [1.0_dp], &
            method='BDF', jac=stiff_cosine_jacobian, jac_matrix=bad_matrix)
        call require(.not. sol%success, 'callback plus constant Jacobian rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'callback plus constant Jacobian status')

        sol = solve_ivp(stiff_linear_rhs, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], &
            method='BDF', jac_matrix=jac_matrix, jac_sparse_matrix=sparse_jac)
        call require(.not. sol%success, 'dense plus CSC Jacobian rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'dense plus CSC Jacobian status')

        sol = solve_ivp(stiff_linear_rhs, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], &
            method='RK45', fun_vectorized=diagonal_stiff_rhs_vectorized)
        call require(.not. sol%success, 'explicit-method vectorized RHS rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'explicit-method vectorized RHS status')

        sol = solve_ivp(stiff_linear_rhs, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], &
            method='RK45', jac_sparse=diagonal_stiff_sparse_jacobian)
        call require(.not. sol%success, 'explicit-method sparse callback rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'explicit-method sparse callback status')

        sol = solve_ivp(stiff_linear_rhs, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], &
            method='LSODA', jac_sparse=diagonal_stiff_sparse_jacobian)
        call require(.not. sol%success, 'LSODA sparse callback rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'LSODA sparse callback status')

        sol = solve_ivp(stiff_linear_rhs, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], &
            method='BDF', jac=stiff_cosine_jacobian, &
            jac_sparse=diagonal_stiff_sparse_jacobian)
        call require(.not. sol%success, 'dense plus sparse callback rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'dense plus sparse callback status')

        sol = solve_ivp(stiff_linear_rhs, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], &
            method='BDF', jac_sparse=invalid_sparse_jacobian)
        call require(.not. sol%success, 'invalid sparse callback output rejected')
        call require(sol%status == ivp_status_invalid_argument, &
            'invalid initial sparse callback status')
    end subroutine test_invalid_jacobian_options

    subroutine test_rk45_exponential()
        type(solve_ivp_result) :: sol
        real(dp) :: expected

        sol = solve_ivp(exponential_rhs, [0.0_dp, 2.0_dp], [1.0_dp], method='RK45', &
            rtol=1.0e-9_dp, atol=1.0e-12_dp)
        expected = exp(2.0_dp)
        call require(sol%success, 'RK45 exponential success')
        call require(sol%status == ivp_status_success, 'RK45 exponential status')
        call check_rel(sol%y(1, size(sol%t)), expected, 2.0e-9_dp, 'RK45 exponential value')
        call require(sol%nfev > 0, 'RK45 nfev')
    end subroutine test_rk45_exponential

    subroutine test_rk45_scipy_dense_output()
        type(solve_ivp_result) :: sol
        real(dp), allocatable :: values(:, :)
        real(dp), parameter :: query(5) = [ &
            0.125_dp, 0.25_dp, 0.5_dp, 0.75_dp, 0.875_dp]
        real(dp), parameter :: expected(5) = [ &
            1.1331908920566731_dp, 1.2840961717243726_dp, &
            1.6486478238062923_dp, 2.1168565883910397_dp, &
            2.3988354233066733_dp]
        integer :: status
        integer :: i

        sol = solve_ivp(exponential_rhs, [0.0_dp, 1.0_dp], [1.0_dp], &
            method='RK45', dense_output=.true., first_step=1.0_dp, max_step=1.0_dp, &
            rtol=1.0e-3_dp, atol=1.0e-6_dp)
        call require(sol%success, 'RK45 SciPy dense success')
        call require(sol%nfev == 7, 'RK45 SciPy dense nfev')
        call sol%sol%evaluate(query, values, status)
        call require(status == ivp_status_success, 'RK45 SciPy dense status')
        do i = 1, size(query)
            call check_abs(values(1, i), expected(i), 8.0e-15_dp, &
                'RK45 SciPy dense value')
        end do

        sol = solve_ivp(exponential_rhs, [0.0_dp, 1.0_dp], [1.0_dp], &
            method='RK45', t_eval=query, first_step=1.0_dp, max_step=1.0_dp, &
            rtol=1.0e-3_dp, atol=1.0e-6_dp)
        call require(sol%success, 'RK45 SciPy t_eval success')
        call require(sol%nfev == 7, 'RK45 SciPy t_eval nfev')
        do i = 1, size(query)
            call check_abs(sol%y(1, i), expected(i), 8.0e-15_dp, &
                'RK45 SciPy t_eval value')
        end do

        sol = solve_ivp(exponential_rhs, [0.0_dp, 1.0_dp], [1.0_dp], &
            method='RK45', dense_output=.true., events=event_y_one_point_five, &
            event_terminal=.true., event_direction=1.0_dp, first_step=1.0_dp, &
            max_step=1.0_dp, &
            rtol=1.0e-3_dp, atol=1.0e-6_dp)
        call require(sol%success, 'RK45 SciPy dense event success')
        call require(sol%status == ivp_status_event, 'RK45 SciPy dense event status')
        call check_abs(sol%t_events(1), 0.4054634731951677_dp, 5.0e-14_dp, &
            'RK45 SciPy dense event time')
        call check_abs(sol%y_events(1, 1), 1.5_dp, 2.0e-14_dp, &
            'RK45 SciPy dense event state')
        call sol%sol%evaluate([0.2_dp], values, status)
        call require(status == ivp_status_success, 'RK45 truncated dense status')
        call check_abs(values(1, 1), 1.2214701486310573_dp, 8.0e-15_dp, &
            'RK45 truncated dense value')
    end subroutine test_rk45_scipy_dense_output

    subroutine test_rk23_harmonic()
        type(solve_ivp_result) :: sol
        real(dp) :: y0(2)
        real(dp) :: pi

        pi = acos(-1.0_dp)
        y0 = [1.0_dp, 0.0_dp]
        sol = solve_ivp(harmonic_rhs, [0.0_dp, 0.5_dp * pi], y0, method='RK23', &
            rtol=2.0e-8_dp, atol=1.0e-10_dp)
        call require(sol%success, 'RK23 harmonic success')
        call check_abs(sol%y(1, size(sol%t)), 0.0_dp, 2.0e-7_dp, 'RK23 harmonic position')
        call check_rel(sol%y(2, size(sol%t)), -1.0_dp, 2.0e-7_dp, 'RK23 harmonic velocity')
    end subroutine test_rk23_harmonic


    subroutine test_dop853_exponential()
        type(solve_ivp_result) :: sol
        real(dp) :: expected

        sol = solve_ivp(exponential_rhs, [0.0_dp, 5.0_dp], [1.0_dp], method='DOP853', &
            rtol=1.0e-11_dp, atol=1.0e-13_dp)
        expected = exp(5.0_dp)
        call require(sol%success, 'DOP853 exponential success')
        call require(sol%status == ivp_status_success, 'DOP853 exponential status')
        call check_rel(sol%y(1, size(sol%t)), expected, 5.0e-10_dp, &
            'DOP853 exponential value')
        call require(sol%nfev > 0, 'DOP853 nfev')
    end subroutine test_dop853_exponential

    subroutine test_dop853_harmonic_t_eval()
        type(solve_ivp_result) :: sol
        real(dp) :: teval(9)
        real(dp) :: y0(2)
        real(dp) :: pi
        integer :: i

        pi = acos(-1.0_dp)
        do i = 1, size(teval)
            teval(i) = real(i - 1, dp) * pi / 8.0_dp
        end do
        y0 = [1.0_dp, 0.0_dp]
        sol = solve_ivp(harmonic_rhs, [0.0_dp, pi], y0, method='DOP853', &
            t_eval=teval, rtol=1.0e-10_dp, atol=1.0e-12_dp, max_step=0.05_dp)
        call require(sol%success, 'DOP853 harmonic success')
        do i = 1, size(teval)
            call check_abs(sol%y(1, i), cos(teval(i)), 2.0e-12_dp, &
                'DOP853 harmonic position')
            call check_abs(sol%y(2, i), -sin(teval(i)), 4.0e-12_dp, &
                'DOP853 harmonic velocity')
        end do
    end subroutine test_dop853_harmonic_t_eval

    subroutine test_dop853_vanderpol()
        type(solve_ivp_result) :: sol
        real(dp), parameter :: expected(2) = [ &
            0.6101970163495998_dp, 0.3912935547253552_dp]

        sol = solve_ivp(vanderpol_rhs, [0.0_dp, 20.0_dp], [0.0_dp, 0.1_dp], &
            method='DOP853', rtol=1.0e-9_dp, atol=1.0e-12_dp)
        call require(sol%success, 'DOP853 Van der Pol success')
        call check_abs(sol%y(1, size(sol%t)), expected(1), 2.0e-10_dp, &
            'DOP853 Van der Pol position')
        call check_abs(sol%y(2, size(sol%t)), expected(2), 2.0e-10_dp, &
            'DOP853 Van der Pol velocity')
        call require(sol%nfev == 710, 'DOP853 Van der Pol nfev')
        call require(size(sol%t) == 54, 'DOP853 Van der Pol accepted nodes')
    end subroutine test_dop853_vanderpol


    subroutine test_dop853_native_dense()
        type(solve_ivp_result) :: sol
        real(dp), allocatable :: values(:, :)
        real(dp), parameter :: query(5) = [0.123_dp, 0.777_dp, 1.234_dp, &
            2.345_dp, 4.321_dp]
        integer :: status
        integer :: i

        sol = solve_ivp(exponential_rhs, [0.0_dp, 5.0_dp], [1.0_dp], &
            method='DOP853', dense_output=.true., rtol=1.0e-11_dp, atol=1.0e-13_dp)
        call require(sol%success, 'DOP853 native dense success')
        call sol%sol%evaluate(query, values, status)
        call require(status == ivp_status_success, 'DOP853 native dense status')
        call require(sol%nfev == 332, 'DOP853 native dense nfev')
        call require(size(sol%t) == 23, 'DOP853 native dense accepted nodes')
        do i = 1, size(query)
            call check_rel(values(1, i), exp(query(i)), 5.0e-11_dp, &
                'DOP853 native dense value')
        end do
    end subroutine test_dop853_native_dense

    subroutine test_dop853_terminal_event()
        type(solve_ivp_result) :: sol
        real(dp) :: expected

        expected = log(2.0_dp)
        sol = solve_ivp(exponential_rhs, [0.0_dp, 3.0_dp], [1.0_dp], &
            method='DOP853', events=event_y_two, event_terminal=.true., &
            event_direction=1.0_dp, rtol=1.0e-10_dp, atol=1.0e-12_dp)
        call require(sol%success, 'DOP853 terminal event success')
        call require(sol%status == ivp_status_event, 'DOP853 terminal event status')
        call require(size(sol%t_events) == 1, 'DOP853 terminal event count')
        call check_abs(sol%t_events(1), expected, 1.0e-10_dp, &
            'DOP853 terminal event time')
        call require(sol%event_indices(1) == 1, 'DOP853 terminal event index')
    end subroutine test_dop853_terminal_event

    subroutine test_radau_exponential_dense()
        type(solve_ivp_result) :: sol
        real(dp), allocatable :: values(:, :)
        real(dp), parameter :: teval(5) = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
        real(dp), parameter :: expected(5) = [ &
            1.0_dp, 1.284025415945068_dp, 1.6487212693496391_dp, &
            2.1170000150975525_dp, 2.7182818284709938_dp]
        real(dp), parameter :: query(3) = [0.123_dp, 0.321_dp, 0.9_dp]
        real(dp), parameter :: dense_expected(3) = [ &
            1.130884420270982_dp, 1.3785055814352614_dp, 2.4596031128707283_dp]
        integer :: i
        integer :: status

        sol = solve_ivp(exponential_rhs, [0.0_dp, 1.0_dp], [1.0_dp], method='Radau', &
            t_eval=teval, dense_output=.true., rtol=1.0e-8_dp, atol=1.0e-10_dp)
        call require(sol%success, 'Radau exponential success')
        call require(sol%status == ivp_status_success, 'Radau exponential status')
        call require(sol%nfev == 239, 'Radau exponential nfev')
        call require(sol%njev == 2, 'Radau exponential njev')
        call require(sol%nlu == 6, 'Radau exponential nlu')
        do i = 1, size(teval)
            call check_abs(sol%y(1, i), expected(i), 2.0e-14_dp, &
                'Radau SciPy t_eval value')
        end do
        call sol%sol%evaluate(query, values, status)
        call require(status == ivp_status_success, 'Radau dense status')
        do i = 1, size(query)
            call check_abs(values(1, i), dense_expected(i), 2.0e-14_dp, &
                'Radau SciPy dense value')
        end do
    end subroutine test_radau_exponential_dense

    subroutine test_radau_stiff_with_jacobian()
        type(solve_ivp_result) :: sol
        real(dp), parameter :: teval(5) = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
        real(dp), parameter :: expected(5) = [ &
            1.0_dp, 0.96890585456529_dp, 0.8775833514286183_dp, &
            0.7316990518410265_dp, 0.5403023067018451_dp]
        integer :: i

        sol = solve_ivp(stiff_cosine_rhs, [0.0_dp, 1.0_dp], [1.0_dp], method='Radau', &
            t_eval=teval, jac=stiff_cosine_jacobian, rtol=1.0e-6_dp, atol=1.0e-9_dp)
        call require(sol%success, 'Radau stiff success')
        call require(sol%nfev == 84, 'Radau stiff nfev')
        call require(sol%njev == 1, 'Radau stiff njev')
        call require(sol%nlu == 16, 'Radau stiff nlu')
        do i = 1, size(teval)
            call check_abs(sol%y(1, i), expected(i), 2.0e-14_dp, &
                'Radau stiff SciPy value')
            call check_abs(sol%y(1, i), cos(teval(i)), 1.2e-5_dp, &
                'Radau stiff analytic value')
        end do
    end subroutine test_radau_stiff_with_jacobian

    subroutine test_radau_multicomponent_fd_jacobian()
        type(solve_ivp_result) :: sol
        real(dp), parameter :: teval(5) = [ &
            0.0_dp, 0.001_dp, 0.01_dp, 0.1_dp, 1.0_dp]
        real(dp), parameter :: expected(2, 5) = reshape([ &
            1.0_dp, 1.0_dp, &
            0.367879446156229_dp, 0.999000499833375_dp, &
            4.539990975083933e-5_dp, 0.9900498337491679_dp, &
            -1.4580174607535287e-14_dp, 0.9048374458787757_dp, &
            1.4597911445030359e-28_dp, 0.36787944192248295_dp], [2, 5])
        integer :: i

        sol = solve_ivp(stiff_linear_rhs, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], &
            method='Radau', t_eval=teval, rtol=1.0e-6_dp, atol=1.0e-9_dp)
        call require(sol%success, 'Radau multicomponent finite-difference success')
        call require(sol%nfev == 770, 'Radau multicomponent finite-difference nfev')
        call require(sol%njev == 2, 'Radau multicomponent finite-difference njev')
        call require(sol%nlu == 40, 'Radau multicomponent finite-difference nlu')
        do i = 1, size(teval)
            call check_abs(sol%y(1, i), expected(1, i), 3.0e-15_dp, &
                'Radau multicomponent stiff component')
            call check_abs(sol%y(2, i), expected(2, i), 3.0e-15_dp, &
                'Radau multicomponent slow component')
        end do
    end subroutine test_radau_multicomponent_fd_jacobian

    subroutine test_radau_backward_dense()
        type(solve_ivp_result) :: sol
        real(dp), allocatable :: values(:, :)
        real(dp), parameter :: teval(5) = [1.0_dp, 0.75_dp, 0.5_dp, 0.25_dp, 0.0_dp]
        real(dp), parameter :: expected(5) = [ &
            2.718281828459045_dp, 2.1170000154605906_dp, 1.6487212693831028_dp, &
            1.2840254156615927_dp, 1.0000000000043083_dp]
        real(dp), parameter :: query(2) = [0.9_dp, 0.321_dp]
        real(dp), parameter :: dense_expected(2) = [ &
            2.4596031129655547_dp, 1.3785055812900584_dp]
        integer :: i
        integer :: status

        sol = solve_ivp(exponential_rhs, [1.0_dp, 0.0_dp], [exp(1.0_dp)], &
            method='Radau', t_eval=teval, dense_output=.true., &
            rtol=1.0e-8_dp, atol=1.0e-10_dp)
        call require(sol%success, 'Radau backward success')
        call require(sol%status == ivp_status_success, 'Radau backward status')
        call require(sol%nfev == 239, 'Radau backward nfev')
        call require(sol%njev == 2, 'Radau backward njev')
        call require(sol%nlu == 6, 'Radau backward nlu')
        do i = 1, size(teval)
            call check_abs(sol%y(1, i), expected(i), 2.0e-14_dp, &
                'Radau backward SciPy t_eval value')
        end do
        call sol%sol%evaluate(query, values, status)
        call require(status == ivp_status_success, 'Radau backward dense status')
        do i = 1, size(query)
            call check_abs(values(1, i), dense_expected(i), 2.0e-14_dp, &
                'Radau backward SciPy dense value')
        end do
    end subroutine test_radau_backward_dense

    subroutine test_radau_terminal_event()
        type(solve_ivp_result) :: sol
        real(dp), allocatable :: values(:, :)
        integer :: status

        sol = solve_ivp(exponential_rhs, [0.0_dp, 3.0_dp], [1.0_dp], method='Radau', &
            dense_output=.true., events=event_y_two, event_terminal=.true., &
            event_direction=1.0_dp, rtol=1.0e-9_dp, atol=1.0e-12_dp)
        call require(sol%success, 'Radau terminal event success')
        call require(sol%status == ivp_status_event, 'Radau terminal event status')
        call require(size(sol%t_events) == 1, 'Radau terminal event count')
        call check_abs(sol%t_events(1), 0.6931471805888538_dp, 2.0e-14_dp, &
            'Radau terminal event time')
        call check_abs(sol%y_events(1, 1), 2.0_dp, 2.0e-14_dp, &
            'Radau terminal event state')
        call require(sol%nfev == 282, 'Radau terminal event nfev')
        call require(sol%njev == 1, 'Radau terminal event njev')
        call require(sol%nlu == 4, 'Radau terminal event nlu')
        call sol%sol%evaluate([0.2_dp], values, status)
        call require(status == ivp_status_success, 'Radau terminal dense status')
        call check_abs(values(1, 1), 1.2214027582991944_dp, 2.0e-14_dp, &
            'Radau terminal dense value')
    end subroutine test_radau_terminal_event

    subroutine test_bdf_exponential_dense()
        type(solve_ivp_result) :: sol
        real(dp), allocatable :: values(:, :)
        real(dp), parameter :: teval(5) = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
        real(dp), parameter :: expected(5) = [ &
            1.0_dp, 1.2840255724127794_dp, 1.6487214845261169_dp, &
            2.1170003130263493_dp, 2.7182822404493594_dp]
        real(dp), parameter :: query(2) = [0.123_dp, 0.777_dp]
        real(dp), parameter :: dense_expected(2) = [ &
            1.1308845344193699_dp, 2.1749379623399030_dp]
        integer :: i
        integer :: status

        sol = solve_ivp(exponential_rhs, [0.0_dp, 1.0_dp], [1.0_dp], method='BDF', &
            t_eval=teval, dense_output=.true., rtol=1.0e-8_dp, atol=1.0e-10_dp)
        call require(sol%success, 'BDF exponential success')
        call require(sol%status == ivp_status_success, 'BDF exponential status')
        call require(sol%nfev == 94, 'BDF exponential nfev')
        call require(sol%njev == 1, 'BDF exponential njev')
        call require(sol%nlu == 12, 'BDF exponential nlu')
        do i = 1, size(teval)
            call check_abs(sol%y(1, i), expected(i), 2.0e-14_dp, &
                'BDF SciPy t_eval value')
        end do
        call sol%sol%evaluate(query, values, status)
        call require(status == ivp_status_success, 'BDF dense status')
        do i = 1, size(query)
            call check_abs(values(1, i), dense_expected(i), 2.0e-14_dp, &
                'BDF SciPy dense value')
        end do
    end subroutine test_bdf_exponential_dense

    subroutine test_bdf_stiff_with_jacobian()
        type(solve_ivp_result) :: sol
        real(dp), parameter :: teval(5) = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
        real(dp), parameter :: expected(5) = [ &
            1.0_dp, 0.9689125714160894_dp, 0.8775823408217067_dp, &
            0.7316887958488855_dp, 0.5403023009227299_dp]
        integer :: i

        sol = solve_ivp(stiff_cosine_rhs, [0.0_dp, 1.0_dp], [1.0_dp], method='BDF', &
            t_eval=teval, jac=stiff_cosine_jacobian, rtol=1.0e-6_dp, atol=1.0e-9_dp)
        call require(sol%success, 'BDF stiff success')
        call require(sol%nfev == 54, 'BDF stiff nfev')
        call require(sol%njev == 2, 'BDF stiff njev')
        call require(sol%nlu == 9, 'BDF stiff nlu')
        do i = 1, size(teval)
            call check_abs(sol%y(1, i), expected(i), 2.0e-14_dp, &
                'BDF stiff SciPy value')
            call check_abs(sol%y(1, i), cos(teval(i)), 3.0e-7_dp, &
                'BDF stiff analytic value')
        end do
    end subroutine test_bdf_stiff_with_jacobian

    subroutine test_bdf_multicomponent_fd_jacobian()
        type(solve_ivp_result) :: sol
        real(dp), parameter :: teval(5) = [ &
            0.0_dp, 0.001_dp, 0.01_dp, 0.1_dp, 1.0_dp]
        real(dp), parameter :: expected(2, 5) = reshape([ &
            1.0_dp, 1.0_dp, &
            0.36787925622560708_dp, 0.99900049983501271_dp, &
            4.5402841308535665e-5_dp, 0.99004983375079136_dp, &
            3.2667643322309101e-13_dp, 0.90483522051524323_dp, &
            7.6156448197415878e-19_dp, 0.36787859693589281_dp], [2, 5])
        integer :: i

        sol = solve_ivp(stiff_linear_rhs, [0.0_dp, 1.0_dp], [1.0_dp, 1.0_dp], &
            method='BDF', t_eval=teval, rtol=1.0e-6_dp, atol=1.0e-9_dp)
        call require(sol%success, 'BDF multicomponent finite-difference success')
        call require(sol%nfev == 332, 'BDF multicomponent finite-difference nfev')
        call require(sol%njev == 1, 'BDF multicomponent finite-difference njev')
        call require(sol%nlu == 34, 'BDF multicomponent finite-difference nlu')
        do i = 1, size(teval)
            call check_abs(sol%y(1, i), expected(1, i), 3.0e-15_dp, &
                'BDF multicomponent stiff component')
            call check_abs(sol%y(2, i), expected(2, i), 3.0e-15_dp, &
                'BDF multicomponent slow component')
        end do
    end subroutine test_bdf_multicomponent_fd_jacobian

    subroutine test_bdf_backward_dense()
        type(solve_ivp_result) :: sol
        real(dp), allocatable :: values(:, :)
        real(dp), parameter :: teval(5) = [1.0_dp, 0.75_dp, 0.5_dp, 0.25_dp, 0.0_dp]
        real(dp), parameter :: expected(5) = [ &
            2.7182818284590451_dp, 2.1169999275509164_dp, 1.6487212047838196_dp, &
            1.2840253873538652_dp, 0.9999999912148333_dp]
        real(dp), parameter :: query(2) = [0.9_dp, 0.321_dp]
        real(dp), parameter :: dense_expected(2) = [ &
            2.4596030287907902_dp, 1.3785055433577336_dp]
        integer :: i
        integer :: status

        sol = solve_ivp(exponential_rhs, [1.0_dp, 0.0_dp], [exp(1.0_dp)], method='BDF', &
            t_eval=teval, dense_output=.true., rtol=1.0e-8_dp, atol=1.0e-10_dp)
        call require(sol%success, 'BDF backward success')
        call require(sol%status == ivp_status_success, 'BDF backward status')
        call require(sol%nfev == 92, 'BDF backward nfev')
        call require(sol%njev == 1, 'BDF backward njev')
        call require(sol%nlu == 12, 'BDF backward nlu')
        do i = 1, size(teval)
            call check_abs(sol%y(1, i), expected(i), 2.0e-14_dp, &
                'BDF backward SciPy t_eval value')
        end do
        call sol%sol%evaluate(query, values, status)
        call require(status == ivp_status_success, 'BDF backward dense status')
        do i = 1, size(query)
            call check_abs(values(1, i), dense_expected(i), 2.0e-14_dp, &
                'BDF backward SciPy dense value')
        end do
    end subroutine test_bdf_backward_dense

    subroutine test_bdf_terminal_event()
        type(solve_ivp_result) :: sol
        real(dp), allocatable :: values(:, :)
        integer :: status

        sol = solve_ivp(exponential_rhs, [0.0_dp, 3.0_dp], [1.0_dp], method='BDF', &
            dense_output=.true., events=event_y_two, event_terminal=.true., &
            event_direction=1.0_dp, rtol=1.0e-9_dp, atol=1.0e-12_dp)
        call require(sol%success, 'BDF terminal event success')
        call require(sol%status == ivp_status_event, 'BDF terminal event status')
        call require(size(sol%t_events) == 1, 'BDF terminal event count')
        call check_abs(sol%t_events(1), 0.6931471709305571_dp, 2.0e-14_dp, &
            'BDF terminal event time')
        call check_abs(sol%y_events(1, 1), 2.0_dp, 2.0e-14_dp, &
            'BDF terminal event state')
        call require(sol%nfev == 96, 'BDF terminal event nfev')
        call require(sol%njev == 1, 'BDF terminal event njev')
        call require(sol%nlu == 12, 'BDF terminal event nlu')
        call sol%sol%evaluate([0.2_dp], values, status)
        call require(status == ivp_status_success, 'BDF terminal dense status')
        call check_abs(values(1, 1), 1.2214027659224616_dp, 2.0e-14_dp, &
            'BDF terminal dense value')
    end subroutine test_bdf_terminal_event

    subroutine test_lsoda_exponential_dense()
        type(solve_ivp_result) :: sol
        real(dp), allocatable :: values(:, :)
        real(dp), parameter :: teval(5) = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
        real(dp), parameter :: query(2) = [0.123_dp, 0.777_dp]
        integer :: status
        integer :: i

        sol = solve_ivp(exponential_rhs, [0.0_dp, 1.0_dp], [1.0_dp], method='LSODA', &
            t_eval=teval, dense_output=.true., rtol=1.0e-8_dp, atol=1.0e-10_dp)
        call require(sol%success, 'LSODA exponential success')
        call require(sol%status == ivp_status_success, 'LSODA exponential status')
        call require(sol%nlu == 0, 'LSODA nonstiff branch avoids LU')
        call require(sol%njev == 0, 'LSODA nonstiff branch avoids Jacobian work')
        do i = 1, size(teval)
            call check_rel(sol%y(1, i), exp(teval(i)), 4.0e-8_dp, &
                'LSODA exponential t_eval value')
        end do
        call sol%sol%evaluate(query, values, status)
        call require(status == ivp_status_success, 'LSODA exponential dense status')
        do i = 1, size(query)
            call check_rel(values(1, i), exp(query(i)), 8.0e-8_dp, &
                'LSODA exponential dense value')
        end do
    end subroutine test_lsoda_exponential_dense

    subroutine test_lsoda_low_order_startup()
        type(lsoda_state) :: state
        real(dp) :: y(1)
        real(dp) :: y_new(1)
        real(dp) :: f(1)
        real(dp) :: f_new(1)
        real(dp) :: t_new
        real(dp) :: dense_d(1, 7)
        real(dp) :: dense_h
        real(dp) :: dense_anchor
        integer :: dense_order
        integer :: nfev
        integer :: njev
        integer :: nlu
        integer :: attempts
        integer :: rejected
        integer :: status
        logical :: dense_valid

        y = 0.0_dp
        f = 1.0_dp
        nfev = 1
        njev = 0
        nlu = 0
        call lsoda_initialize(state, y, f, 0.01_dp)
        call lsoda_advance(state, unit_rhs, 0.0_dp, y, f, 1.0_dp, 1.0_dp, &
            1.0e-8_dp, 1.0e-10_dp, 1.0_dp, y_new, f_new, t_new, dense_d, dense_h, &
            dense_order, dense_anchor, dense_valid, nfev, njev, nlu, 100, attempts, &
            rejected, status)
        call require(status == lsoda_status_success, 'LSODA low-order startup status')
        call require(attempts == 1 .and. rejected == 0, &
            'LSODA low-order startup accepts constant RHS')
        call require(nfev == 2, 'LSODA low-order startup converges in one RHS call')
        call require(njev == 0 .and. nlu == 0, 'LSODA low-order startup avoids stiff work')
        call require(state%history_count == 2, 'LSODA low-order startup extends history')
        call check_abs(t_new, 0.01_dp, 8.0_dp * epsilon(1.0_dp), &
            'LSODA low-order startup time')
        call check_abs(y_new(1), 0.01_dp, 8.0_dp * epsilon(1.0_dp), &
            'LSODA low-order startup value')
        call check_abs(f_new(1), 1.0_dp, 8.0_dp * epsilon(1.0_dp), &
            'LSODA low-order startup derivative')
    end subroutine test_lsoda_low_order_startup

    subroutine test_lsoda_startup_convergence_failure()
        type(solve_ivp_result) :: sol

        ! With this deliberately large first step, two fixed Adams corrections give
        ! a deceptively small tolerance-scaled error even though functional iteration
        ! is diverging. The ODEPACK-style convergence test must reject the step before
        ! the user min_step floor permits a retry.
        sol = solve_ivp(rapid_growth_rhs, [0.0_dp, 0.1_dp], [1.0_dp], method='LSODA', &
            first_step=0.1_dp, max_step=0.1_dp, min_step=0.09_dp, &
            rtol=1.0_dp, atol=1.0_dp)
        call require(.not. sol%success, 'LSODA startup corrector detects divergence')
        call require(sol%nfev == 3, 'LSODA startup corrector stops after lost contraction')
        call require(sol%njev == 0 .and. sol%nlu == 0, &
            'LSODA startup corrector failure avoids stiff work')
    end subroutine test_lsoda_startup_convergence_failure

    subroutine test_lsoda_high_order_adams()
        type(lsoda_state) :: state
        real(dp) :: t
        real(dp) :: t_new
        real(dp) :: y(1)
        real(dp) :: y_new(1)
        real(dp) :: f(1)
        real(dp) :: f_new(1)
        real(dp) :: dense_d(1, 7)
        real(dp) :: dense_h
        real(dp) :: dense_anchor
        integer :: dense_order
        integer :: nfev
        integer :: njev
        integer :: nlu
        integer :: attempts
        integer :: rejected
        integer :: status
        integer :: max_history
        integer :: max_order
        integer :: order_before
        integer :: previous_order
        integer :: last_order_change
        integer :: k
        logical :: dense_valid

        t = 0.0_dp
        y = 1.0_dp
        f = 1.0_dp
        nfev = 1
        njev = 0
        nlu = 0
        max_history = 1
        call lsoda_initialize(state, y, f, 0.01_dp)
        max_order = state%adams_method_order
        previous_order = state%adams_method_order
        last_order_change = 0
        do k = 1, 20
            order_before = state%adams_method_order
            call lsoda_advance(state, exponential_rhs, t, y, f, 0.3_dp, 1.0_dp, &
                1.0e-8_dp, 1.0e-10_dp, 0.01_dp, y_new, f_new, t_new, dense_d, &
                dense_h, dense_order, dense_anchor, dense_valid, nfev, njev, nlu, &
                100, attempts, rejected, status)
            call require(status == lsoda_status_success, 'LSODA high-order Adams step status')
            call require(dense_valid .and. dense_order == -order_before, &
                'LSODA high-order Adams returns native Nordsieck dense state')
            if (state%adams_method_order /= previous_order) then
                call require(k - last_order_change >= previous_order + 1, &
                    'LSODA IALTH delays Adams order changes')
                call require(state%adams_ialth == state%adams_method_order + 1, &
                    'LSODA IALTH resets after Adams order change')
                last_order_change = k
                previous_order = state%adams_method_order
            end if
            max_history = max(max_history, state%history_count)
            max_order = max(max_order, state%adams_method_order)
            t = t_new
            y = y_new
            f = f_new
        end do
        call require(max_history == 11, 'LSODA retains full Adams history')
        call require(max_order >= 5, 'LSODA reaches high Adams method order')
        call require(nlu == 0, 'LSODA high-order Adams remains nonstiff')
        call check_rel(y(1), exp(t), 2.0e-8_dp, 'LSODA high-order Adams value')
    end subroutine test_lsoda_high_order_adams

    subroutine test_lsoda_order12_nordsieck()
        type(lsoda_state) :: state
        real(dp), parameter :: h = 1.0e-3_dp
        real(dp) :: dense_d(1, 13)
        real(dp) :: dense_h
        real(dp) :: dense_anchor
        real(dp) :: factorial
        real(dp) :: f(1)
        real(dp) :: f_new(1)
        real(dp) :: t_new
        real(dp) :: y(1)
        real(dp) :: y_new(1)
        integer :: attempts
        integer :: dense_order
        integer :: j
        integer :: nfev
        integer :: njev
        integer :: nlu
        integer :: rejected
        integer :: status
        logical :: dense_valid

        y = 1.0_dp
        f = 1.0_dp
        nfev = 1
        njev = 0
        nlu = 0
        call lsoda_initialize(state, y, f, h)
        state%adams_yh = 0.0_dp
        state%adams_yh(:, 1) = y
        factorial = 1.0_dp
        do j = 1, 12
            factorial = factorial * real(j, dp)
            state%adams_yh(1, j + 1) = h**j / factorial
        end do
        state%adams_method_order = 12
        state%adams_yh_order = 12
        state%adams_yh_h = h
        state%adams_yh_valid = .true.
        state%adams_ialth = 13

        call lsoda_advance(state, exponential_rhs, 0.0_dp, y, f, 1.0_dp, 1.0_dp, &
            1.0e-10_dp, 1.0e-12_dp, h, y_new, f_new, t_new, dense_d, dense_h, &
            dense_order, dense_anchor, dense_valid, nfev, njev, nlu, 100, attempts, &
            rejected, status)
        call require(status == lsoda_status_success, 'LSODA order-12 Adams step status')
        call require(attempts == 1 .and. rejected == 0, &
            'LSODA order-12 Adams accepts exact Nordsieck state')
        call require(state%adams_method_order == 12, 'LSODA retains Adams order 12')
        call require(dense_valid .and. dense_order == -12, &
            'LSODA returns order-12 Adams dense state')
        call check_rel(y_new(1), exp(h), 2.0e-14_dp, 'LSODA order-12 Adams value')
        call check_rel(f_new(1), exp(h), 2.0e-14_dp, 'LSODA order-12 Adams derivative')
        call require(njev == 0 .and. nlu == 0, 'LSODA order-12 Adams remains nonstiff')
    end subroutine test_lsoda_order12_nordsieck

    subroutine test_lsoda_variable_step_history()
        type(lsoda_state) :: state
        real(dp) :: t
        real(dp) :: t_new
        real(dp) :: y(1)
        real(dp) :: y_new(1)
        real(dp) :: f(1)
        real(dp) :: f_new(1)
        real(dp) :: dense_d(1, 7)
        real(dp) :: dense_h
        real(dp) :: dense_anchor
        real(dp) :: previous_h
        real(dp) :: accepted_h
        integer :: dense_order
        integer :: nfev
        integer :: njev
        integer :: nlu
        integer :: attempts
        integer :: rejected
        integer :: status
        integer :: max_order
        integer :: min_history_after_six
        integer :: k
        logical :: dense_valid
        logical :: saw_step_change
        logical :: reached_six

        t = 0.0_dp
        y = 1.0_dp
        f = 1.0_dp
        nfev = 1
        njev = 0
        nlu = 0
        max_order = 0
        min_history_after_six = huge(1)
        previous_h = -1.0_dp
        saw_step_change = .false.
        reached_six = .false.
        call lsoda_initialize(state, y, f, 0.001_dp)
        do k = 1, 1000
            if (t >= 1.0_dp - 8.0_dp * epsilon(1.0_dp)) exit
            call lsoda_advance(state, exponential_rhs, t, y, f, 1.0_dp, 1.0_dp, &
                1.0e-7_dp, 1.0e-9_dp, 0.2_dp, y_new, f_new, t_new, dense_d, &
                dense_h, dense_order, dense_anchor, dense_valid, nfev, njev, nlu, &
                100, attempts, rejected, status)
            call require(status == lsoda_status_success, &
                'LSODA variable-step Adams step status')
            accepted_h = t_new - t
            if (previous_h > 0.0_dp .and. state%history_count >= 4 .and. &
                    abs(accepted_h - previous_h) > 1.0e-12_dp) then
                saw_step_change = .true.
            end if
            previous_h = accepted_h
            if (state%history_count >= 6) reached_six = .true.
            if (reached_six) then
                min_history_after_six = min(min_history_after_six, state%history_count)
            end if
            max_order = max(max_order, state%adams_method_order)
            t = t_new
            y = y_new
            f = f_new
        end do
        call require(k <= 1000, 'LSODA variable-step Adams terminates')
        call require(saw_step_change, 'LSODA variable-step Adams changes step size')
        call require(min_history_after_six >= 6, &
            'LSODA variable-step Adams keeps nonuniform history')
        call require(max_order >= 6, 'LSODA variable-step Adams raises order')
        call require(njev == 0, 'LSODA variable-step Adams avoids nonstiff Jacobians')
        call require(nlu == 0, 'LSODA variable-step Adams remains nonstiff')
        call check_rel(y(1), exp(1.0_dp), 4.0e-7_dp, &
            'LSODA variable-step Adams value')
    end subroutine test_lsoda_variable_step_history

    subroutine test_lsoda_nordsieck_rescale()
        type(lsoda_state) :: state
        real(dp) :: t
        real(dp) :: t_new
        real(dp) :: y(1)
        real(dp) :: y_new(1)
        real(dp) :: f(1)
        real(dp) :: f_new(1)
        real(dp) :: dense_d(1, 7)
        real(dp) :: dense_h
        real(dp) :: dense_anchor
        real(dp) :: forced_step
        real(dp) :: old_reference_h
        integer :: dense_order
        integer :: nfev
        integer :: njev
        integer :: nlu
        integer :: attempts
        integer :: rejected
        integer :: status
        integer :: k
        logical :: dense_valid

        t = 0.0_dp
        y = 0.0_dp
        f = 0.0_dp
        nfev = 1
        njev = 0
        nlu = 0
        call lsoda_initialize(state, y, f, 0.01_dp)
        do k = 1, 8
            call lsoda_advance(state, quadratic_rhs, t, y, f, 0.2_dp, 1.0_dp, &
                1.0e-10_dp, 1.0e-12_dp, 0.01_dp, y_new, f_new, t_new, dense_d, &
                dense_h, dense_order, dense_anchor, dense_valid, nfev, njev, nlu, &
                100, attempts, rejected, status)
            call require(status == lsoda_status_success, 'LSODA Nordsieck buildup status')
            t = t_new
            y = y_new
            f = f_new
        end do

        call require(state%adams_yh_valid, 'LSODA builds Adams Nordsieck state')
        call require(state%adams_yh_order == state%adams_method_order, &
            'LSODA Nordsieck state tracks Adams method order')
        call check_abs(state%adams_yh(1, 1), y(1), 2.0e-14_dp, &
            'LSODA Nordsieck state value column')
        call check_abs(state%adams_yh(1, 2), state%adams_yh_h * f(1), 2.0e-14_dp, &
            'LSODA Nordsieck first derivative column')
        call check_abs(state%adams_yh(1, 3), state%adams_yh_h**2, 2.0e-13_dp, &
            'LSODA Nordsieck second derivative column')

        old_reference_h = abs(state%adams_yh_h)
        forced_step = old_reference_h / 4.0_dp
        call lsoda_advance(state, quadratic_rhs, t, y, f, 0.2_dp, 1.0_dp, &
            1.0e-10_dp, 1.0e-12_dp, forced_step, y_new, f_new, t_new, dense_d, &
            dense_h, dense_order, dense_anchor, dense_valid, nfev, njev, nlu, &
            100, attempts, rejected, status)
        call require(status == lsoda_status_success, 'LSODA Nordsieck rescale status')
        call require(abs((t_new - t) - forced_step) <= 32.0_dp * epsilon(1.0_dp), &
            'LSODA Nordsieck rescale uses forced smaller step')
        call require(state%adams_yh_valid, 'LSODA keeps Nordsieck state after rescale')
        call require(abs(state%adams_yh_h - forced_step) <= 32.0_dp * epsilon(1.0_dp), &
            'LSODA Nordsieck state adopts new step scale')
        call check_abs(state%adams_yh(1, 1), y_new(1), 2.0e-14_dp, &
            'LSODA rescaled Nordsieck value column')
        call check_abs(state%adams_yh(1, 2), state%adams_yh_h * f_new(1), 2.0e-14_dp, &
            'LSODA rescaled Nordsieck first derivative column')
        call check_abs(state%adams_yh(1, 3), state%adams_yh_h**2, 2.0e-13_dp, &
            'LSODA rescaled Nordsieck second derivative column')
        call check_rel(y_new(1), t_new**2, 2.0e-12_dp, &
            'LSODA Nordsieck rescale quadratic value')
    end subroutine test_lsoda_nordsieck_rescale

    subroutine test_lsoda_backward_nordsieck()
        type(solve_ivp_result) :: sol
        real(dp), parameter :: teval(5) = [1.0_dp, 0.75_dp, 0.5_dp, 0.25_dp, 0.0_dp]
        integer :: i

        sol = solve_ivp(exponential_rhs, [1.0_dp, 0.0_dp], [exp(1.0_dp)], &
            method='LSODA', t_eval=teval, rtol=1.0e-8_dp, atol=1.0e-10_dp)
        call require(sol%success, 'LSODA backward Nordsieck success')
        call require(sol%status == ivp_status_success, 'LSODA backward Nordsieck status')
        call require(sol%njev == 0, 'LSODA backward Nordsieck avoids Jacobian work')
        call require(sol%nlu == 0, 'LSODA backward Nordsieck avoids LU work')
        do i = 1, size(teval)
            call check_rel(sol%y(1, i), exp(teval(i)), 6.0e-8_dp, &
                'LSODA backward Nordsieck value')
        end do
    end subroutine test_lsoda_backward_nordsieck

    subroutine test_lsoda_error_failure_recovery()
        type(lsoda_state) :: state
        real(dp) :: t
        real(dp) :: t_new
        real(dp) :: y(1)
        real(dp) :: y_new(1)
        real(dp) :: f(1)
        real(dp) :: f_new(1)
        real(dp) :: dense_d(1, 7)
        real(dp) :: dense_h
        real(dp) :: dense_anchor
        integer :: dense_order
        integer :: nfev
        integer :: njev
        integer :: nlu
        integer :: attempts
        integer :: rejected
        integer :: status
        integer :: k
        logical :: dense_valid

        t = 0.0_dp
        y = 1.0_dp
        f = 1.0_dp
        nfev = 1
        njev = 0
        nlu = 0
        call lsoda_initialize(state, y, f, 0.01_dp)
        do k = 1, 8
            call lsoda_advance(state, exponential_rhs, t, y, f, 10.0_dp, 1.0_dp, &
                1.0e-8_dp, 1.0e-10_dp, 0.01_dp, y_new, f_new, t_new, dense_d, &
                dense_h, dense_order, dense_anchor, dense_valid, nfev, njev, nlu, &
                100, attempts, rejected, status)
            call require(status == lsoda_status_success, &
                'LSODA failure-recovery buildup status')
            t = t_new
            y = y_new
            f = f_new
        end do

        state%h_abs = 0.49_dp
        call lsoda_advance(state, exponential_rhs, t, y, f, 10.0_dp, 1.0_dp, &
            1.0e-14_dp, 1.0e-16_dp, 0.49_dp, y_new, f_new, t_new, dense_d, &
            dense_h, dense_order, dense_anchor, dense_valid, nfev, njev, nlu, &
            100, attempts, rejected, status)
        call require(status == lsoda_status_success, 'LSODA error-failure recovery status')
        call require(rejected >= 4, 'LSODA error-failure recovery exercises repeated failure')
        call require(state%history_count == 2, &
            'LSODA error-failure recovery restarts low-order history')
        call require(njev == 0, 'LSODA rejected Adams trials avoid Jacobian probes')
        call require(nlu == 0, 'LSODA rejected Adams trials avoid BDF work')
        call check_rel(y_new(1), exp(t_new), 2.0e-8_dp, &
            'LSODA error-failure recovery value')
    end subroutine test_lsoda_error_failure_recovery

    subroutine test_lsoda_roundoff_switch_branch()
        type(lsoda_state) :: state
        real(dp) :: t
        real(dp) :: t_new
        real(dp) :: y(1)
        real(dp) :: y_new(1)
        real(dp) :: f(1)
        real(dp) :: f_new(1)
        real(dp) :: dense_d(1, 7)
        real(dp) :: dense_h
        real(dp) :: dense_anchor
        integer :: dense_order
        integer :: nfev
        integer :: njev
        integer :: nlu
        integer :: attempts
        integer :: rejected
        integer :: status
        integer :: k
        logical :: dense_valid

        t = 0.0_dp
        y = 0.0_dp
        f = 1.0_dp
        nfev = 1
        njev = 0
        nlu = 0
        call lsoda_initialize(state, y, f, 0.01_dp)
        do k = 1, 3
            call lsoda_advance(state, unit_rhs, t, y, f, 1.0_dp, 1.0_dp, &
                1.0e-8_dp, 1.0e-10_dp, 0.01_dp, y_new, f_new, t_new, dense_d, &
                dense_h, dense_order, dense_anchor, dense_valid, nfev, njev, nlu, &
                100, attempts, rejected, status, unit_jacobian)
            call require(status == lsoda_status_success, &
                'LSODA roundoff-switch buildup status')
            t = t_new
            y = y_new
            f = f_new
        end do

        ! Model an ODEPACK PDLAST estimate from a preceding correction.  The constant
        ! RHS makes the present local-error estimate exactly roundoff-polluted.
        state%switch_countdown = -1
        state%adams_pdlast = 1000.0_dp
        state%h_abs = 0.01_dp
        call lsoda_advance(state, unit_rhs, t, y, f, 1.0_dp, 1.0_dp, &
            1.0e-8_dp, 1.0e-10_dp, 0.01_dp, y_new, f_new, t_new, dense_d, &
            dense_h, dense_order, dense_anchor, dense_valid, nfev, njev, nlu, &
            100, attempts, rejected, status, unit_jacobian)
        call require(status == lsoda_status_success, 'LSODA roundoff-switch status')
        call require(rejected == 0, 'LSODA roundoff-switch accepts Adams step')
        call require(njev == 0, 'LSODA roundoff-switch defers Jacobian until BDF trial')
        call require(nlu == 0, 'LSODA roundoff-switch defers BDF factorization')
        call check_abs(state%adams_pdlast, 0.0_dp, 0.0_dp, &
            'LSODA roundoff-switch resets Adams Lipschitz history')
        call check_abs(y_new(1), t_new, 16.0_dp * epsilon(1.0_dp), &
            'LSODA roundoff-switch Adams value')

        t = t_new
        y = y_new
        f = f_new
        call lsoda_advance(state, unit_rhs, t, y, f, 1.0_dp, 1.0_dp, &
            1.0e-8_dp, 1.0e-10_dp, 0.01_dp, y_new, f_new, t_new, dense_d, &
            dense_h, dense_order, dense_anchor, dense_valid, nfev, njev, nlu, &
            100, attempts, rejected, status, unit_jacobian)
        call require(status == lsoda_status_success, 'LSODA roundoff-switch BDF status')
        call require(nlu > 0, 'LSODA roundoff-switch enters BDF branch')
        call check_abs(y_new(1), t_new, 64.0_dp * epsilon(1.0_dp), &
            'LSODA roundoff-switch BDF value')
    end subroutine test_lsoda_roundoff_switch_branch

    subroutine test_lsoda_switch_hysteresis()
        type(lsoda_state) :: state
        real(dp) :: t
        real(dp) :: t_new
        real(dp) :: y(1)
        real(dp) :: y_new(1)
        real(dp) :: f(1)
        real(dp) :: f_new(1)
        real(dp) :: f_reference(1)
        real(dp) :: dense_d(1, 7)
        real(dp) :: dense_h
        real(dp) :: dense_anchor
        integer :: dense_order
        integer :: nfev
        integer :: njev
        integer :: nlu
        integer :: nfev_before
        integer :: attempts
        integer :: rejected
        integer :: status
        integer :: k
        logical :: dense_valid

        t = 0.0_dp
        y = 1.0_dp
        call stiff_cosine_rhs(t, y, f)
        nfev = 1
        njev = 0
        nlu = 0
        call lsoda_initialize(state, y, f, 1.0e-4_dp)
        do k = 1, 20
            call lsoda_advance(state, stiff_cosine_rhs, t, y, f, 1.0_dp, 1.0_dp, &
                1.0e-6_dp, 1.0e-9_dp, huge(1.0_dp), y_new, f_new, t_new, dense_d, &
                dense_h, dense_order, dense_anchor, dense_valid, nfev, njev, nlu, &
                100, attempts, rejected, status, stiff_cosine_jacobian)
            call require(status == lsoda_status_success, &
                'LSODA switch hysteresis Adams status')
            call require(njev == 0, 'LSODA switch hysteresis delays Jacobian work')
            call require(nlu == 0, 'LSODA switch hysteresis delays BDF work')
            if (k >= 4) then
                call require(state%h_abs <= 5.76e-4_dp, &
                    'LSODA Adams stability boundary limits stiff step')
            end if
            t = t_new
            y = y_new
            f = f_new
        end do
        call lsoda_advance(state, stiff_cosine_rhs, t, y, f, 1.0_dp, 1.0_dp, &
            1.0e-6_dp, 1.0e-9_dp, huge(1.0_dp), y_new, f_new, t_new, dense_d, &
            dense_h, dense_order, dense_anchor, dense_valid, nfev, njev, nlu, &
            100, attempts, rejected, status, stiff_cosine_jacobian)
        call require(status == lsoda_status_success, 'LSODA switch hysteresis switch status')
        call require(njev == 0, 'LSODA switch hysteresis defers Jacobian until BDF trial')
        call require(nlu == 0, 'LSODA switch hysteresis switches after accepted Adams step')

        t = t_new
        y = y_new
        f = f_new
        nfev_before = nfev
        call lsoda_advance(state, stiff_cosine_rhs, t, y, f, 1.0_dp, 1.0_dp, &
            1.0e-6_dp, 1.0e-9_dp, huge(1.0_dp), y_new, f_new, t_new, dense_d, &
            dense_h, dense_order, dense_anchor, dense_valid, nfev, njev, nlu, &
            100, attempts, rejected, status, stiff_cosine_jacobian)
        call require(status == lsoda_status_success, 'LSODA switch hysteresis BDF status')
        call require(nlu > 0, 'LSODA switch hysteresis enters BDF factorization')
        call require(nfev - nfev_before == 1, &
            'LSODA BDF first chord correction can converge from predictor RHS')
        call stiff_cosine_rhs(t_new, y_new, f_reference)
        call check_abs(f_new(1), f_reference(1), 2.0e-12_dp, &
            'LSODA BDF implied endpoint derivative')
    end subroutine test_lsoda_switch_hysteresis

    subroutine test_lsoda_stiff_switch()
        type(solve_ivp_result) :: sol
        real(dp), parameter :: teval(5) = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
        integer :: i

        sol = solve_ivp(stiff_cosine_rhs, [0.0_dp, 1.0_dp], [1.0_dp], method='LSODA', &
            t_eval=teval, jac=stiff_cosine_jacobian, rtol=1.0e-6_dp, atol=1.0e-9_dp)
        call require(sol%success, 'LSODA stiff success')
        call require(sol%status == ivp_status_success, 'LSODA stiff status')
        call require(sol%njev > 0, 'LSODA stiff Jacobian evaluations')
        call require(sol%nlu > 0, 'LSODA switched to BDF branch')
        do i = 1, size(teval)
            call check_abs(sol%y(1, i), cos(teval(i)), 7.0e-7_dp, &
                'LSODA stiff analytic value')
        end do
    end subroutine test_lsoda_stiff_switch

    subroutine test_lsoda_banded_fd()
        type(solve_ivp_result) :: dense_sol
        type(solve_ivp_result) :: band_sol
        real(dp) :: y0(8)
        integer :: dense_calls
        integer :: band_calls
        integer :: i

        y0 = 1.0_dp
        diagonal_rhs_calls = 0
        dense_sol = solve_ivp(diagonal_stiff_rhs, [0.0_dp, 0.5_dp], y0, &
            method='LSODA', rtol=2.0e-6_dp, atol=1.0e-9_dp)
        dense_calls = diagonal_rhs_calls
        call require(dense_sol%success, 'LSODA dense finite-difference success')
        call require(dense_sol%nlu > 0, 'LSODA dense finite-difference switched to BDF')
        call require(dense_sol%nfev == dense_calls, &
            'LSODA dense finite-difference nfev includes Jacobian probes')

        diagonal_rhs_calls = 0
        band_sol = solve_ivp(diagonal_stiff_rhs, [0.0_dp, 0.5_dp], y0, &
            method='LSODA', lband=0, uband=0, rtol=2.0e-6_dp, atol=1.0e-9_dp)
        band_calls = diagonal_rhs_calls
        call require(band_sol%success, 'LSODA banded finite-difference success')
        call require(band_sol%nlu > 0, 'LSODA banded finite-difference switched to BDF')
        call require(band_sol%nfev == band_calls, &
            'LSODA banded finite-difference nfev includes Jacobian probes')
        call require(band_calls < dense_calls, 'LSODA banded finite differences reduce RHS calls')
        do i = 1, size(y0)
            call check_abs(band_sol%y(i, size(band_sol%t)), cos(0.5_dp), 3.0e-6_dp, &
                'LSODA banded finite-difference value')
        end do
    end subroutine test_lsoda_banded_fd

    subroutine test_lsoda_banded_jacobian()
        type(solve_ivp_result) :: sol
        real(dp) :: y0(6)
        integer :: i

        y0 = 1.0_dp
        sol = solve_ivp(diagonal_stiff_rhs, [0.0_dp, 0.5_dp], y0, method='LSODA', &
            jac=diagonal_stiff_jacobian_packed, lband=0, uband=0, &
            rtol=2.0e-6_dp, atol=1.0e-9_dp)
        call require(sol%success, 'LSODA packed banded Jacobian success')
        call require(sol%njev > 0, 'LSODA packed banded Jacobian evaluations')
        call require(sol%nlu > 0, 'LSODA packed banded Jacobian switched to BDF')
        do i = 1, size(y0)
            call check_abs(sol%y(i, size(sol%t)), cos(0.5_dp), 3.0e-6_dp, &
                'LSODA packed banded Jacobian value')
        end do
    end subroutine test_lsoda_banded_jacobian

    subroutine test_lsoda_tridiagonal_bands()
        type(solve_ivp_result) :: sol
        real(dp) :: y0(7)
        integer :: i

        y0 = 1.0_dp
        sol = solve_ivp(tridiagonal_stiff_rhs, [0.0_dp, 0.5_dp], y0, method='LSODA', &
            lband=1, uband=1, rtol=2.0e-6_dp, atol=1.0e-9_dp)
        call require(sol%success, 'LSODA tridiagonal finite-difference success')
        call require(sol%nlu > 0, 'LSODA tridiagonal finite-difference switched to BDF')
        do i = 1, size(y0)
            call check_abs(sol%y(i, size(sol%t)), cos(0.5_dp), 4.0e-6_dp, &
                'LSODA tridiagonal finite-difference value')
        end do

        sol = solve_ivp(tridiagonal_stiff_rhs, [0.0_dp, 0.5_dp], y0, method='LSODA', &
            jac=tridiagonal_stiff_jacobian_packed, lband=1, uband=1, &
            rtol=2.0e-6_dp, atol=1.0e-9_dp)
        call require(sol%success, 'LSODA tridiagonal packed Jacobian success')
        call require(sol%njev > 0, 'LSODA tridiagonal packed Jacobian evaluations')
        call require(sol%nlu > 0, 'LSODA tridiagonal packed Jacobian switched to BDF')
        do i = 1, size(y0)
            call check_abs(sol%y(i, size(sol%t)), cos(0.5_dp), 4.0e-6_dp, &
                'LSODA tridiagonal packed Jacobian value')
        end do
    end subroutine test_lsoda_tridiagonal_bands

    subroutine test_lsoda_min_step()
        type(solve_ivp_result) :: sol

        sol = solve_ivp(exponential_rhs, [0.0_dp, 1.0_dp], [1.0_dp], method='LSODA', &
            first_step=1.0_dp, max_step=1.0_dp, min_step=0.9_dp, &
            rtol=1.0e-12_dp, atol=1.0e-14_dp)
        call require(.not. sol%success, 'LSODA min_step blocks required retry')

        sol = solve_ivp(exponential_rhs, [0.0_dp, 1.0_dp], [1.0_dp], method='LSODA', &
            min_step=-1.0_dp)
        call require(.not. sol%success, 'LSODA negative min_step failure')
        call require(sol%status == ivp_status_invalid_argument, &
            'LSODA negative min_step status')

        sol = solve_ivp(exponential_rhs, [0.0_dp, 1.0_dp], [1.0_dp], method='LSODA', &
            lband=1, uband=0)
        call require(.not. sol%success, 'LSODA invalid bandwidth failure')
        call require(sol%status == ivp_status_invalid_argument, &
            'LSODA invalid bandwidth status')
    end subroutine test_lsoda_min_step

    subroutine test_lsoda_terminal_event()
        type(solve_ivp_result) :: sol
        real(dp), allocatable :: values(:, :)
        integer :: status

        sol = solve_ivp(exponential_rhs, [0.0_dp, 3.0_dp], [1.0_dp], method='LSODA', &
            dense_output=.true., events=event_y_two, event_terminal=.true., &
            event_direction=1.0_dp, rtol=1.0e-8_dp, atol=1.0e-10_dp)
        call require(sol%success, 'LSODA terminal event success')
        call require(sol%status == ivp_status_event, 'LSODA terminal event status')
        call require(size(sol%t_events) == 1, 'LSODA terminal event count')
        call check_abs(sol%t_events(1), log(2.0_dp), 8.0e-8_dp, &
            'LSODA terminal event time')
        call check_abs(sol%y_events(1, 1), 2.0_dp, 2.0e-12_dp, &
            'LSODA terminal event state')
        call sol%sol%evaluate([0.2_dp], values, status)
        call require(status == ivp_status_success, 'LSODA terminal dense status')
        call check_rel(values(1, 1), exp(0.2_dp), 8.0e-8_dp, &
            'LSODA terminal dense value')
    end subroutine test_lsoda_terminal_event

    subroutine test_multiple_events()
        type(solve_ivp_result) :: sol
        type(ivp_event_set) :: event_set
        real(dp), allocatable :: event_times(:)
        real(dp), allocatable :: event_states(:, :)

        event_set%n_events = 2
        event_set%fun => two_events
        event_set%terminal = [.false., .true.]
        event_set%direction = [1.0_dp, 1.0_dp]
        sol = solve_ivp(unit_rhs, [0.0_dp, 3.0_dp], [0.0_dp], method='DOP853', &
            first_step=3.0_dp, event_set=event_set, rtol=1.0e-10_dp, atol=1.0e-12_dp)
        call require(sol%success, 'multiple events success')
        call require(sol%status == ivp_status_event, 'multiple events terminal status')
        call require(size(sol%t_events) == 2, 'multiple events count')
        call check_abs(sol%t_events(1), 1.0_dp, 1.0e-12_dp, 'first vector event time')
        call check_abs(sol%t_events(2), 2.0_dp, 1.0e-12_dp, 'second vector event time')
        call require(all(sol%event_indices == [1, 2]), 'multiple event indices')
        call require(sol%n_event_functions == 2, 'multiple event function count')
        call require(all(sol%event_counts == [1, 1]), 'multiple per-event counts')
        event_times = sol%event_times_for(1)
        event_states = sol%event_states_for(2)
        call require(size(event_times) == 1, 'multiple grouped event time count')
        call check_abs(event_times(1), 1.0_dp, 1.0e-12_dp, &
            'multiple grouped event time')
        call require(size(event_states, 1) == 1 .and. size(event_states, 2) == 1, &
            'multiple grouped event state shape')
        call check_abs(event_states(1, 1), 2.0_dp, 1.0e-12_dp, &
            'multiple grouped event state')
        call check_abs(sol%t(size(sol%t)), 2.0_dp, 1.0e-12_dp, &
            'multiple events terminal endpoint')
    end subroutine test_multiple_events

    subroutine test_event_occurrence_limit()
        type(solve_ivp_result) :: sol
        real(dp), allocatable :: event_times(:)
        real(dp), allocatable :: event_states(:, :)

        sol = solve_ivp(unit_rhs, [0.0_dp, 3.0_dp], [0.0_dp], method='RK45', &
            events=periodic_event, event_terminal=.true., event_max_occurrences=2, &
            max_step=0.2_dp, &
            rtol=1.0e-10_dp, atol=1.0e-12_dp)
        call require(sol%success, 'event occurrence limit success')
        call require(sol%status == ivp_status_event, 'event occurrence limit status')
        call require(size(sol%t_events) == 2, 'event occurrence limit root count')
        call require(sol%n_event_functions == 1, 'event occurrence limit function count')
        call require(all(sol%event_counts == [2]), 'event occurrence limit count')
        call check_abs(sol%t_events(1), 0.25_dp, 5.0e-13_dp, &
            'event occurrence first root')
        call check_abs(sol%t_events(2), 1.25_dp, 5.0e-13_dp, &
            'event occurrence second root')
        call check_abs(sol%t(size(sol%t)), 1.25_dp, 5.0e-13_dp, &
            'event occurrence terminal endpoint')
        event_times = sol%event_times_for(1)
        event_states = sol%event_states_for(1)
        call require(size(event_times) == 2, 'event occurrence grouped times')
        call require(size(event_states, 1) == 1 .and. size(event_states, 2) == 2, &
            'event occurrence grouped states')
        call check_abs(event_states(1, 1), 0.25_dp, 5.0e-13_dp, &
            'event occurrence first state')
        call check_abs(event_states(1, 2), 1.25_dp, 5.0e-13_dp, &
            'event occurrence second state')
    end subroutine test_event_occurrence_limit

    subroutine test_terminal_event_root_order()
        type(solve_ivp_result) :: sol
        type(ivp_event_set) :: event_set
        real(dp), allocatable :: later_times(:)

        event_set%n_events = 2
        event_set%fun => close_events
        event_set%direction = [1.0_dp, 1.0_dp]
        event_set%max_occurrences = [1, 0]
        sol = solve_ivp(unit_rhs, [0.0_dp, 1.5_dp], [0.0_dp], method='RK45', &
            first_step=1.0_dp, max_step=1.0_dp, event_set=event_set, &
            rtol=1.0e-10_dp, atol=1.0e-12_dp)
        call require(sol%success, 'terminal event root order success')
        call require(sol%status == ivp_status_event, 'terminal event root order status')
        call require(size(sol%t_events) == 1, 'terminal event root order stored roots')
        call require(sol%event_indices(1) == 1, 'terminal event root order index')
        call require(all(sol%event_counts == [1, 0]), 'terminal event root order counts')
        call check_abs(sol%t_events(1), 0.6_dp, 5.0e-13_dp, &
            'terminal event root order time')
        later_times = sol%event_times_for(2)
        call require(size(later_times) == 0, 'terminal event root order truncation')
    end subroutine test_terminal_event_root_order

    subroutine test_event_direction_filter()
        type(solve_ivp_result) :: sol
        type(ivp_event_set) :: event_set

        event_set%n_events = 2
        event_set%fun => two_events
        event_set%terminal = [.false., .true.]
        event_set%direction = [1.0_dp, -1.0_dp]
        sol = solve_ivp(unit_rhs, [0.0_dp, 3.0_dp], [0.0_dp], method='RK45', &
            first_step=3.0_dp, event_set=event_set, rtol=1.0e-10_dp, atol=1.0e-12_dp)
        call require(sol%success, 'event direction filter success')
        call require(sol%status == ivp_status_success, 'event direction filter status')
        call require(size(sol%t_events) == 1, 'event direction filter count')
        call require(sol%event_indices(1) == 1, 'event direction filter index')
        call check_abs(sol%t(size(sol%t)), 3.0_dp, 1.0e-14_dp, &
            'event direction filter endpoint')
    end subroutine test_event_direction_filter

    subroutine test_initial_event_direction_filter()
        type(solve_ivp_result) :: sol

        sol = solve_ivp(unit_rhs, [0.0_dp, 1.0_dp], [0.0_dp], method='RK45', &
            events=event_y_zero, event_max_occurrences=1, event_direction=-1.0_dp, &
            first_step=0.5_dp, rtol=1.0e-10_dp, atol=1.0e-12_dp)
        call require(sol%success, 'initial event direction success')
        call require(sol%status == ivp_status_success, 'initial event direction status')
        call require(size(sol%t_events) == 0, 'initial event direction root count')
        call require(all(sol%event_counts == [0]), 'initial event direction count')
    end subroutine test_initial_event_direction_filter

    subroutine test_backward_and_t_eval()
        type(solve_ivp_result) :: sol
        real(dp) :: teval(5)
        integer :: i

        teval = [2.0_dp, 1.5_dp, 1.0_dp, 0.5_dp, 0.0_dp]
        sol = solve_ivp(exponential_rhs, [2.0_dp, 0.0_dp], [exp(2.0_dp)], &
            method='RK45', t_eval=teval, rtol=1.0e-9_dp, atol=1.0e-12_dp)
        call require(sol%success, 'backward integration success')
        call require(size(sol%t) == size(teval), 'backward t_eval size')
        do i = 1, size(teval)
            call check_rel(sol%y(1, i), exp(teval(i)), 5.0e-8_dp, 'backward t_eval value')
        end do
    end subroutine test_backward_and_t_eval

    subroutine test_dense_output()
        type(solve_ivp_result) :: sol
        real(dp), allocatable :: values(:, :)
        real(dp) :: query(3)
        integer :: status
        integer :: i

        query = [0.15_dp, 0.73_dp, 1.61_dp]
        sol = solve_ivp(exponential_rhs, [0.0_dp, 2.0_dp], [1.0_dp], method='RK45', &
            dense_output=.true., rtol=1.0e-8_dp, atol=1.0e-11_dp)
        call sol%sol%evaluate(query, values, status)
        call require(status == ivp_status_success, 'dense output status')
        do i = 1, size(query)
            call check_rel(values(1, i), exp(query(i)), 2.0e-6_dp, 'dense output value')
        end do
    end subroutine test_dense_output

    subroutine test_terminal_event()
        type(solve_ivp_result) :: sol
        real(dp) :: expected

        expected = log(2.0_dp)
        sol = solve_ivp(exponential_rhs, [0.0_dp, 3.0_dp], [1.0_dp], method='RK45', &
            events=event_y_two, event_terminal=.true., event_direction=1.0_dp, &
            rtol=1.0e-9_dp, atol=1.0e-12_dp)
        call require(sol%success, 'terminal event success')
        call require(sol%status == ivp_status_event, 'terminal event status')
        call require(size(sol%t_events) == 1, 'terminal event count')
        call check_abs(sol%t_events(1), expected, 5.0e-8_dp, 'terminal event time')
        call check_abs(sol%y_events(1, 1), 2.0_dp, 2.0e-8_dp, 'terminal event state')
        call check_abs(sol%t(size(sol%t)), sol%t_events(1), 1.0e-14_dp, 'terminal endpoint')
    end subroutine test_terminal_event

    subroutine test_invalid_method()
        type(solve_ivp_result) :: sol

        sol = solve_ivp(exponential_rhs, [0.0_dp, 1.0_dp], [1.0_dp], method='bogus')
        call require(.not. sol%success, 'invalid method failure')
        call require(sol%status == ivp_status_invalid_argument, 'invalid method status')

        sol = solve_ivp(exponential_rhs, [0.0_dp, 1.0_dp], [1.0_dp], &
            events=event_y_two, event_max_occurrences=-1)
        call require(.not. sol%success, 'invalid event occurrence limit failure')
        call require(sol%status == ivp_status_invalid_argument, &
            'invalid event occurrence limit status')
    end subroutine test_invalid_method

    subroutine exponential_rhs(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable, unused by autonomous system
        real(dp), intent(in) :: y(:) !! current state
        real(dp), intent(out) :: dydt(:) !! derivative of current state

        if (t /= t) error stop 'unexpected NaN time'
        dydt = y
    end subroutine exponential_rhs

    subroutine rapid_growth_rhs(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable, unused by autonomous system
        real(dp), intent(in) :: y(:) !! current state
        real(dp), intent(out) :: dydt(:) !! rapidly expanding derivative

        if (t /= t) error stop 'unexpected NaN time'
        dydt = 100.0_dp * y
    end subroutine rapid_growth_rhs

    subroutine quadratic_rhs(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable in y=t**2
        real(dp), intent(in) :: y(:) !! current scalar state, unused analytically
        real(dp), intent(out) :: dydt(:) !! exact derivative 2*t

        dydt(1) = 2.0_dp * t + 0.0_dp * y(1)
    end subroutine quadratic_rhs

    subroutine harmonic_rhs(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable, unused by autonomous system
        real(dp), intent(in) :: y(:) !! position and velocity
        real(dp), intent(out) :: dydt(:) !! velocity and acceleration

        if (t /= t) error stop 'unexpected NaN time'
        dydt(1) = y(2)
        dydt(2) = -y(1)
    end subroutine harmonic_rhs

    subroutine vanderpol_rhs(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable, unused by autonomous system
        real(dp), intent(in) :: y(:) !! position and velocity
        real(dp), intent(out) :: dydt(:) !! velocity and acceleration

        if (t /= t) error stop 'unexpected NaN time'
        dydt(1) = y(2)
        dydt(2) = 0.2_dp * (1.0_dp - y(1)**2) * y(2) - y(1)
    end subroutine vanderpol_rhs

    function event_y_one_point_five(t, y) result(value)
        real(dp), intent(in) :: t !! independent variable, unused by event definition
        real(dp), intent(in) :: y(:) !! current state
        real(dp) :: value

        if (t /= t) error stop 'unexpected NaN time'
        value = y(1) - 1.5_dp
    end function event_y_one_point_five

    function event_y_two(t, y) result(value)
        real(dp), intent(in) :: t !! independent variable, unused by event definition
        real(dp), intent(in) :: y(:) !! current state
        real(dp) :: value

        if (t /= t) error stop 'unexpected NaN time'
        value = y(1) - 2.0_dp
    end function event_y_two

    function event_y_zero(t, y) result(value)
        real(dp), intent(in) :: t !! independent variable, unused by event definition
        real(dp), intent(in) :: y(:) !! current state
        real(dp) :: value

        if (t /= t) error stop 'unexpected NaN time'
        value = y(1)
    end function event_y_zero

    function periodic_event(t, y) result(value)
        real(dp), intent(in) :: t !! independent variable
        real(dp), intent(in) :: y(:) !! current state, used only for shape validation
        real(dp) :: value

        if (size(y) /= 1) error stop 'unexpected periodic event state'
        value = sin(acos(-1.0_dp) * (t - 0.25_dp))
    end function periodic_event


    subroutine stiff_linear_rhs(t, y, dydt)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydt(:)

        if (t /= t .or. size(y) /= 2) error stop 'unexpected stiff linear input'
        dydt(1) = -1000.0_dp * y(1)
        dydt(2) = -y(2)
    end subroutine stiff_linear_rhs

    subroutine diagonal_stiff_rhs(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable
        real(dp), intent(in) :: y(:) !! diagonal stiff state
        real(dp), intent(out) :: dydt(:) !! derivative with exact solution cos(t)

        integer :: i
        real(dp) :: lambda

        diagonal_rhs_calls = diagonal_rhs_calls + 1
        do i = 1, size(y)
            lambda = 400.0_dp + 100.0_dp * real(i, dp)
            dydt(i) = -lambda * (y(i) - cos(t)) - sin(t)
        end do
    end subroutine diagonal_stiff_rhs

    subroutine diagonal_stiff_rhs_vectorized(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable shared by all state columns
        real(dp), intent(in) :: y(:, :) !! diagonal stiff state columns
        real(dp), intent(out) :: dydt(:, :) !! derivatives for each state column

        integer :: i
        integer :: k
        real(dp) :: lambda

        diagonal_vectorized_calls = diagonal_vectorized_calls + 1
        diagonal_vectorized_max_columns = max(diagonal_vectorized_max_columns, size(y, 2))
        do k = 1, size(y, 2)
            do i = 1, size(y, 1)
                lambda = 400.0_dp + 100.0_dp * real(i, dp)
                dydt(i, k) = -lambda * (y(i, k) - cos(t)) - sin(t)
            end do
        end do
    end subroutine diagonal_stiff_rhs_vectorized

    subroutine diagonal_stiff_sparse_jacobian(t, y, jac)
        real(dp), intent(in) :: t !! independent variable
        real(dp), intent(in) :: y(:) !! current diagonal stiff state
        type(ivp_sparse_matrix), intent(out) :: jac !! one-based diagonal CSC Jacobian

        integer :: indices(size(y))
        integer :: indptr(size(y) + 1)
        real(dp) :: data(size(y))
        integer :: i

        sparse_jacobian_calls = sparse_jacobian_calls + 1
        do i = 1, size(y)
            indptr(i) = i
            indices(i) = i
            data(i) = -(400.0_dp + 100.0_dp * real(i, dp)) + 0.0_dp * t + 0.0_dp * y(i)
        end do
        indptr(size(y) + 1) = size(y) + 1
        jac = ivp_csc_matrix(size(y), size(y), indptr, indices, data)
    end subroutine diagonal_stiff_sparse_jacobian

    subroutine invalid_sparse_jacobian(t, y, jac)
        real(dp), intent(in) :: t !! independent variable
        real(dp), intent(in) :: y(:) !! current state
        type(ivp_sparse_matrix), intent(out) :: jac !! deliberately malformed CSC Jacobian

        if (t /= t .or. size(y) < 1) error stop 'unexpected invalid sparse Jacobian input'
        jac = ivp_csc_matrix(2, 2, [1, 3, 2], [1, 2], [-1000.0_dp, -1.0_dp])
    end subroutine invalid_sparse_jacobian

    subroutine diagonal_stiff_jacobian_packed(t, y, jac)
        real(dp), intent(in) :: t !! independent variable, unused by constant-in-y Jacobian
        real(dp), intent(in) :: y(:) !! current diagonal stiff state
        real(dp), intent(out) :: jac(:, :) !! packed diagonal Jacobian, shape (1,n)

        integer :: i

        if (t /= t .or. size(jac, 1) /= 1 .or. size(jac, 2) /= size(y)) then
            jac = huge(1.0_dp)
            return
        end if
        do i = 1, size(y)
            jac(1, i) = -(400.0_dp + 100.0_dp * real(i, dp))
        end do
    end subroutine diagonal_stiff_jacobian_packed

    subroutine tridiagonal_stiff_rhs(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable
        real(dp), intent(in) :: y(:) !! tridiagonal stiff state
        real(dp), intent(out) :: dydt(:) !! derivative with exact solution cos(t)

        real(dp) :: shifted(size(y))
        integer :: i

        shifted = y - cos(t)
        dydt = -1000.0_dp * shifted - sin(t)
        if (size(y) > 1) then
            dydt(1) = dydt(1) + 100.0_dp * shifted(2)
            do i = 2, size(y) - 1
                dydt(i) = dydt(i) + 100.0_dp * (shifted(i - 1) + shifted(i + 1))
            end do
            dydt(size(y)) = dydt(size(y)) + 100.0_dp * shifted(size(y) - 1)
        end if
    end subroutine tridiagonal_stiff_rhs

    subroutine tridiagonal_stiff_jacobian_packed(t, y, jac)
        real(dp), intent(in) :: t !! independent variable, unused by constant Jacobian
        real(dp), intent(in) :: y(:) !! current tridiagonal stiff state
        real(dp), intent(out) :: jac(:, :) !! packed Jacobian, shape (3,n)

        integer :: j

        if (t /= t .or. size(jac, 1) /= 3 .or. size(jac, 2) /= size(y)) then
            jac = huge(1.0_dp)
            return
        end if
        jac = 0.0_dp
        do j = 1, size(y)
            jac(2, j) = -1000.0_dp
            if (j > 1) jac(1, j) = 100.0_dp
            if (j < size(y)) jac(3, j) = 100.0_dp
        end do
    end subroutine tridiagonal_stiff_jacobian_packed

    subroutine multiscale_decay_rhs(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable, unused by autonomous system
        real(dp), intent(in) :: y(:) !! two-component state
        real(dp), intent(out) :: dydt(:) !! derivative with separated decay rates

        if (t /= t .or. size(y) /= 2) error stop 'unexpected multiscale decay input'
        dydt(1) = -y(1)
        dydt(2) = -10.0_dp * y(2)
    end subroutine multiscale_decay_rhs

    subroutine stiff_cosine_rhs(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable
        real(dp), intent(in) :: y(:) !! current scalar state
        real(dp), intent(out) :: dydt(:) !! stiff derivative with exact solution cos(t)

        dydt(1) = -1000.0_dp * (y(1) - cos(t)) - sin(t)
    end subroutine stiff_cosine_rhs

    subroutine stiff_cosine_jacobian(t, y, jac)
        real(dp), intent(in) :: t !! independent variable, unused by constant Jacobian
        real(dp), intent(in) :: y(:) !! current state, unused by constant Jacobian
        real(dp), intent(out) :: jac(:, :) !! Jacobian of stiff_cosine_rhs with respect to y

        if (t /= t .or. size(y) /= 1) error stop 'unexpected stiff Jacobian input'
        jac = 0.0_dp
        jac(1, 1) = -1000.0_dp
    end subroutine stiff_cosine_jacobian

    subroutine unit_rhs(t, y, dydt)
        real(dp), intent(in) :: t !! independent variable, unused by constant system
        real(dp), intent(in) :: y(:) !! current state, unused by constant system
        real(dp), intent(out) :: dydt(:) !! constant derivative

        if (t /= t .or. size(y) <= 0) error stop 'unexpected unit_rhs input'
        dydt = 1.0_dp
    end subroutine unit_rhs


    subroutine unit_jacobian(t, y, jac)
        real(dp), intent(in) :: t !! independent variable
        real(dp), intent(in) :: y(:) !! scalar state
        real(dp), intent(out) :: jac(:, :) !! zero Jacobian, shape (1,1)

        if (t /= t .or. size(y) /= 1 .or. size(jac, 1) /= 1 .or. size(jac, 2) /= 1) then
            error stop 'unexpected unit Jacobian input'
        end if
        jac = 0.0_dp
    end subroutine unit_jacobian

    subroutine two_events(t, y, values)
        real(dp), intent(in) :: t !! independent variable, unused by event definitions
        real(dp), intent(in) :: y(:) !! current scalar state
        real(dp), intent(out) :: values(:) !! event values y-1 and y-2

        if (t /= t) error stop 'unexpected NaN time'
        values(1) = y(1) - 1.0_dp
        values(2) = y(1) - 2.0_dp
    end subroutine two_events

    subroutine close_events(t, y, values)
        real(dp), intent(in) :: t !! independent variable, unused by event definitions
        real(dp), intent(in) :: y(:) !! current scalar state
        real(dp), intent(out) :: values(:) !! event values y-0.6 and y-0.8

        if (t /= t) error stop 'unexpected NaN time'
        values(1) = y(1) - 0.6_dp
        values(2) = y(1) - 0.8_dp
    end subroutine close_events

    subroutine require(condition, label)
        logical, intent(in) :: condition !! condition that must be true
        character(len=*), intent(in) :: label !! failure label

        if (.not. condition) then
            print *, 'FAIL: ', trim(label)
            error stop 1
        end if
    end subroutine require

    subroutine check_abs(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tolerance !! maximum absolute error
        character(len=*), intent(in) :: label !! failure label

        if (abs(actual - expected) > tolerance) then
            print *, 'FAIL: ', trim(label)
            print *, 'actual  = ', actual
            print *, 'expected= ', expected
            error stop 1
        end if
    end subroutine check_abs

    subroutine check_rel(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! nonzero reference value
        real(dp), intent(in) :: tolerance !! maximum relative error
        character(len=*), intent(in) :: label !! failure label

        if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
            print *, 'FAIL: ', trim(label)
            print *, 'actual  = ', actual
            print *, 'expected= ', expected
            error stop 1
        end if
    end subroutine check_rel

end program test_integrate_ivp
