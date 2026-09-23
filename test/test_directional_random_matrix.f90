! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors

program test_directional_random_matrix
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_bessel_i, only : besseli_ratio, log_besseli_scaled
    use scifort_kinds, only : dp
    use scifort_linalg, only : symmetric_eigen_jacobi, linalg_status_success
    use scifort_random, only : rng_seed, rng_state
    use scifort_random_matrix_helpers, only : real_determinant
    use scifort_stats
    use test_checks, only : check_close, check_true
    use test_directional_random_matrix_reference
    implicit none

    integer :: failures

    failures = 0
    call test_bessel_i(failures)
    call test_vonmises_fisher(failures)
    call test_uniform_direction(failures)
    call test_orthogonal_groups(failures)
    call test_unitary_group(failures)
    call test_random_correlation(failures)
    call test_random_table(failures)

    if (failures /= 0) then
        print '(a,i0)', 'test_directional_random_matrix: FAIL, failures=', failures
        error stop 1
    end if
    print '(a)', 'test_directional_random_matrix: PASS'

contains

    subroutine test_bessel_i(failures)
        integer, intent(inout) :: failures !! accumulated test failure count

        call check_close('besseli scaled log', log_besseli_scaled(1.5_dp, 7.0_dp), &
            besseli_log_ref, 2.0e-14_dp, 2.0e-14_dp, failures)
        call check_close('besseli ratio', besseli_ratio(1.5_dp, 7.0_dp), &
            besseli_ratio_ref, 2.0e-14_dp, 2.0e-14_dp, failures)
        call check_close('besseli large scaled log', log_besseli_scaled(50.0_dp, 1000.0_dp), &
            besseli_large_log_ref, 3.0e-11_dp, 3.0e-11_dp, failures)
        call check_close('besseli large ratio', besseli_ratio(50.0_dp, 1000.0_dp), &
            besseli_large_ratio_ref, 8.0e-13_dp, 8.0e-13_dp, failures)
        call check_close('besseli ratio zero', besseli_ratio(2.0_dp, 0.0_dp), &
            0.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_true('besseli invalid order nan', &
            ieee_is_nan(log_besseli_scaled(-1.0_dp, 2.0_dp)), failures)
    end subroutine test_bessel_i

    subroutine test_vonmises_fisher(failures)
        integer, intent(inout) :: failures !! accumulated test failure count

        integer :: info
        integer :: j
        real(dp) :: fit_kappa
        real(dp) :: fit_mu(3)
        real(dp) :: invalid_mu(3)
        real(dp) :: sample1(5)
        real(dp) :: sample2d(2)
        real(dp) :: sample2(5)
        real(dp) :: samples3(3,20)
        real(dp) :: samples5(5,20)
        type(rng_state) :: state1
        type(rng_state) :: state2

        call check_close('vmf 2d logpdf', vonmises_fisher_logpdf(vmf_x2_ref, vmf_mu2_ref, vmf_kappa2_ref), &
            vmf_logpdf2_ref, 5.0e-13_dp, 5.0e-13_dp, failures)
        call check_close('vmf 2d pdf', vonmises_fisher_pdf(vmf_x2_ref, vmf_mu2_ref, vmf_kappa2_ref), &
            vmf_pdf2_ref, 5.0e-13_dp, 5.0e-13_dp, failures)
        call check_close('vmf 2d entropy', vonmises_fisher_entropy(vmf_mu2_ref, vmf_kappa2_ref), &
            vmf_entropy2_ref, 5.0e-13_dp, 5.0e-13_dp, failures)
        call check_close('vmf 3d logpdf', vonmises_fisher_logpdf(vmf_x3_ref, vmf_mu3_ref, vmf_kappa3_ref), &
            vmf_logpdf3_ref, 5.0e-13_dp, 5.0e-13_dp, failures)
        call check_close('vmf 3d pdf', vonmises_fisher_pdf(vmf_x3_ref, vmf_mu3_ref, vmf_kappa3_ref), &
            vmf_pdf3_ref, 5.0e-13_dp, 5.0e-13_dp, failures)
        call check_close('vmf 3d entropy', vonmises_fisher_entropy(vmf_mu3_ref, vmf_kappa3_ref), &
            vmf_entropy3_ref, 5.0e-13_dp, 5.0e-13_dp, failures)
        call check_close('vmf 5d logpdf', vonmises_fisher_logpdf(vmf_x5_ref, vmf_mu5_ref, vmf_kappa5_ref), &
            vmf_logpdf5_ref, 2.0e-12_dp, 2.0e-12_dp, failures)
        call check_close('vmf 5d pdf', vonmises_fisher_pdf(vmf_x5_ref, vmf_mu5_ref, vmf_kappa5_ref), &
            vmf_pdf5_ref, 2.0e-12_dp, 2.0e-12_dp, failures)
        call check_close('vmf 5d entropy', vonmises_fisher_entropy(vmf_mu5_ref, vmf_kappa5_ref), &
            vmf_entropy5_ref, 2.0e-12_dp, 2.0e-12_dp, failures)

        call vonmises_fisher_fit(vmf_fit_data_ref, fit_mu, fit_kappa, info)
        call check_true('vmf fit status', info == vonmises_fisher_status_success, failures)
        do j = 1, 3
            call check_close('vmf fit mu', fit_mu(j), vmf_fit_mu_ref(j), &
                2.0e-12_dp, 2.0e-12_dp, failures)
        end do
        call check_close('vmf fit kappa', fit_kappa, vmf_fit_kappa_ref, &
            2.0e-10_dp, 2.0e-10_dp, failures)

        invalid_mu = [1.0_dp, 1.0_dp, 0.0_dp]
        call check_true('vmf invalid mu nan', &
            ieee_is_nan(vonmises_fisher_pdf(vmf_x3_ref, invalid_mu, 2.0_dp)), failures)
        call check_true('vmf zero kappa nan', &
            ieee_is_nan(vonmises_fisher_entropy(vmf_mu3_ref, 0.0_dp)), failures)

        call rng_seed(state1, 44001)
        call rng_seed(state2, 44001)
        call vonmises_fisher_rvs(state1, vmf_mu5_ref, vmf_kappa5_ref, sample1, info)
        call check_true('vmf rvs status', info == vonmises_fisher_status_success, failures)
        call vonmises_fisher_rvs(state2, vmf_mu5_ref, vmf_kappa5_ref, sample2, info)
        call check_true('vmf rvs deterministic', all(sample1 == sample2), failures)
        call check_close('vmf rvs norm', sqrt(dot_product(sample1,sample1)), &
            1.0_dp, 2.0e-14_dp, 2.0e-14_dp, failures)

        call rng_seed(state1, 44004)
        call vonmises_fisher_rvs(state1, vmf_mu2_ref, vmf_kappa2_ref, sample2d, info)
        call check_true('vmf 2d rvs status', info == vonmises_fisher_status_success, failures)
        call check_close('vmf 2d rvs norm', sqrt(dot_product(sample2d,sample2d)), &
            1.0_dp, 3.0e-14_dp, 3.0e-14_dp, failures)

        call rng_seed(state1, 44002)
        call vonmises_fisher_rvs_array(state1, vmf_mu3_ref, vmf_kappa3_ref, samples3, info)
        call check_true('vmf 3d array status', info == vonmises_fisher_status_success, failures)
        do j = 1, size(samples3,2)
            call check_close('vmf 3d array norm', sqrt(dot_product(samples3(:,j),samples3(:,j))), &
                1.0_dp, 3.0e-14_dp, 3.0e-14_dp, failures)
        end do
        call rng_seed(state1, 44003)
        call vonmises_fisher_rvs_array(state1, vmf_mu5_ref, vmf_kappa5_ref, samples5, info)
        call check_true('vmf 5d array status', info == vonmises_fisher_status_success, failures)
        do j = 1, size(samples5,2)
            call check_close('vmf 5d array norm', sqrt(dot_product(samples5(:,j),samples5(:,j))), &
                1.0_dp, 3.0e-14_dp, 3.0e-14_dp, failures)
        end do
    end subroutine test_vonmises_fisher

    subroutine test_uniform_direction(failures)
        integer, intent(inout) :: failures !! accumulated test failure count

        integer :: info
        integer :: j
        real(dp) :: sample1(4)
        real(dp) :: sample2(4)
        real(dp) :: samples(4,10)
        type(rng_state) :: state1
        type(rng_state) :: state2

        call rng_seed(state1, 45001)
        call rng_seed(state2, 45001)
        call uniform_direction_rvs(state1, sample1, info)
        call check_true('uniform direction status', info == uniform_direction_status_success, failures)
        call uniform_direction_rvs(state2, sample2, info)
        call check_true('uniform direction deterministic', all(sample1 == sample2), failures)
        call check_close('uniform direction norm', sqrt(dot_product(sample1,sample1)), &
            1.0_dp, 2.0e-14_dp, 2.0e-14_dp, failures)
        call uniform_direction_rvs_array(state1, samples, info)
        call check_true('uniform direction array status', info == uniform_direction_status_success, failures)
        do j = 1, size(samples,2)
            call check_close('uniform direction array norm', sqrt(dot_product(samples(:,j),samples(:,j))), &
                1.0_dp, 2.0e-14_dp, 2.0e-14_dp, failures)
        end do
    end subroutine test_uniform_direction

    subroutine test_orthogonal_groups(failures)
        integer, intent(inout) :: failures !! accumulated test failure count

        integer :: info
        integer :: k
        real(dp) :: identity(4,4)
        real(dp) :: q1(4,4)
        real(dp) :: q2(4,4)
        real(dp) :: qs(4,4,3)
        real(dp) :: so1(4,4)
        real(dp) :: so2(4,4)
        real(dp) :: sos(4,4,3)
        type(rng_state) :: state1
        type(rng_state) :: state2

        identity = real_identity(4)
        call rng_seed(state1, 46001)
        call rng_seed(state2, 46001)
        call ortho_group_rvs(state1, q1, info)
        call check_true('ortho status', info == ortho_group_status_success, failures)
        call ortho_group_rvs(state2, q2, info)
        call check_true('ortho deterministic', all(q1 == q2), failures)
        call check_true('ortho orthogonal', &
            maxval(abs(matmul(transpose(q1),q1)-identity)) < 2.0e-13_dp, failures)
        call check_close('ortho determinant magnitude', abs(real_determinant(q1)), &
            1.0_dp, 2.0e-13_dp, 2.0e-13_dp, failures)
        call ortho_group_rvs_array(state1, qs, info)
        call check_true('ortho array status', info == ortho_group_status_success, failures)
        do k = 1, size(qs,3)
            call check_true('ortho array orthogonal', &
                maxval(abs(matmul(transpose(qs(:,:,k)),qs(:,:,k))-identity)) < 3.0e-13_dp, failures)
        end do

        call rng_seed(state1, 46002)
        call rng_seed(state2, 46002)
        call special_ortho_group_rvs(state1, so1, info)
        call check_true('special ortho status', info == special_ortho_group_status_success, failures)
        call special_ortho_group_rvs(state2, so2, info)
        call check_true('special ortho deterministic', all(so1 == so2), failures)
        call check_true('special ortho orthogonal', &
            maxval(abs(matmul(transpose(so1),so1)-identity)) < 2.0e-13_dp, failures)
        call check_close('special ortho determinant', real_determinant(so1), &
            1.0_dp, 2.0e-13_dp, 2.0e-13_dp, failures)
        call special_ortho_group_rvs_array(state1, sos, info)
        call check_true('special ortho array status', info == special_ortho_group_status_success, failures)
        do k = 1, size(sos,3)
            call check_close('special ortho array determinant', real_determinant(sos(:,:,k)), &
                1.0_dp, 3.0e-13_dp, 3.0e-13_dp, failures)
        end do
    end subroutine test_orthogonal_groups

    subroutine test_unitary_group(failures)
        integer, intent(inout) :: failures !! accumulated test failure count

        integer :: info
        integer :: k
        complex(dp) :: identity(3,3)
        complex(dp) :: q1(3,3)
        complex(dp) :: q2(3,3)
        complex(dp) :: qs(3,3,3)
        type(rng_state) :: state1
        type(rng_state) :: state2

        identity = complex_identity(3)
        call rng_seed(state1, 47001)
        call rng_seed(state2, 47001)
        call unitary_group_rvs(state1, q1, info)
        call check_true('unitary status', info == unitary_group_status_success, failures)
        call unitary_group_rvs(state2, q2, info)
        call check_true('unitary deterministic', all(q1 == q2), failures)
        call check_true('unitary orthonormal', &
            maxval(abs(matmul(transpose(conjg(q1)),q1)-identity)) < 2.0e-13_dp, failures)
        call unitary_group_rvs_array(state1, qs, info)
        call check_true('unitary array status', info == unitary_group_status_success, failures)
        do k = 1, size(qs,3)
            call check_true('unitary array orthonormal', &
                maxval(abs(matmul(transpose(conjg(qs(:,:,k))),qs(:,:,k))-identity)) < 3.0e-13_dp, failures)
        end do
    end subroutine test_unitary_group

    subroutine test_random_correlation(failures)
        integer, intent(inout) :: failures !! accumulated test failure count

        integer :: eig_info
        integer :: i
        integer :: info
        real(dp) :: eigenvalues(4)
        real(dp) :: eigenvectors(4,4)
        real(dp) :: eigs(4)
        real(dp) :: r1(4,4)
        real(dp) :: r2(4,4)
        real(dp) :: rank_eigs(3)
        real(dp) :: rank_values(3)
        real(dp) :: rank_vectors(3,3)
        real(dp) :: rank_sample(3,3)
        type(rng_state) :: state1
        type(rng_state) :: state2

        eigs = [0.5_dp, 0.8_dp, 1.2_dp, 1.5_dp]
        call rng_seed(state1, 48001)
        call rng_seed(state2, 48001)
        call random_correlation_rvs(state1, eigs, r1, info)
        call check_true('random correlation status', info == random_correlation_status_success, failures)
        call random_correlation_rvs(state2, eigs, r2, info)
        call check_true('random correlation deterministic', all(r1 == r2), failures)
        call check_true('random correlation symmetric', maxval(abs(r1-transpose(r1))) < 2.0e-14_dp, failures)
        do i = 1, 4
            call check_close('random correlation diagonal', r1(i,i), 1.0_dp, &
                0.0_dp, 0.0_dp, failures)
        end do
        call symmetric_eigen_jacobi(r1, eigenvalues, eigenvectors, eig_info)
        call check_true('random correlation eig status', eig_info == linalg_status_success, failures)
        do i = 1, 4
            call check_close('random correlation eigenvalue', eigenvalues(i), eigs(i), &
                2.0e-11_dp, 2.0e-11_dp, failures)
        end do

        rank_eigs = [0.0_dp, 1.0_dp, 2.0_dp]
        call random_correlation_rvs(state1, rank_eigs, rank_sample, info)
        call check_true('random correlation singular status', &
            info == random_correlation_status_success, failures)
        call symmetric_eigen_jacobi(rank_sample, rank_values, rank_vectors, eig_info)
        do i = 1, 3
            call check_close('random correlation singular eigenvalue', rank_values(i), rank_eigs(i), &
                3.0e-11_dp, 3.0e-11_dp, failures)
        end do
        call random_correlation_rvs(state1, [0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp], r1, info)
        call check_true('random correlation invalid trace', &
            info == random_correlation_status_invalid_parameter, failures)
    end subroutine test_random_correlation

    subroutine test_random_table(failures)
        integer, intent(inout) :: failures !! accumulated test failure count

        integer :: i
        integer :: info
        integer :: j
        integer :: sample1(2,3)
        integer :: sample2(2,3)
        integer :: samples(2,3,5)
        integer :: wrong(2,3)
        real(dp) :: mu(2,3)
        type(rng_state) :: state1
        type(rng_state) :: state2

        call check_close('random table logpmf', &
            random_table_logpmf(random_table_x_ref, random_table_row_ref, random_table_col_ref), &
            random_table_logpmf_ref, 2.0e-14_dp, 2.0e-14_dp, failures)
        call check_close('random table pmf', &
            random_table_pmf(random_table_x_ref, random_table_row_ref, random_table_col_ref), &
            random_table_pmf_ref, 2.0e-14_dp, 2.0e-14_dp, failures)
        mu = random_table_mean(random_table_row_ref, random_table_col_ref)
        do j = 1, size(mu,2)
            do i = 1, size(mu,1)
                call check_close('random table mean', mu(i,j), random_table_mean_ref(i,j), &
                    2.0e-14_dp, 2.0e-14_dp, failures)
            end do
        end do
        wrong = random_table_x_ref
        wrong(1,1) = wrong(1,1) + 1
        call check_close('random table wrong margins pmf', &
            random_table_pmf(wrong, random_table_row_ref, random_table_col_ref), &
            0.0_dp, 0.0_dp, 0.0_dp, failures)
        call check_true('random table invalid margins nan', &
            ieee_is_nan(random_table_pmf(random_table_x_ref, [7,5], random_table_col_ref)), failures)
        call check_true('random table zero-total mean nan', &
            all(ieee_is_nan(random_table_mean([0], [0]))), failures)

        call rng_seed(state1, 49001)
        call rng_seed(state2, 49001)
        call random_table_rvs(state1, random_table_row_ref, random_table_col_ref, sample1, info)
        call check_true('random table rvs status', info == random_table_status_success, failures)
        call random_table_rvs(state2, random_table_row_ref, random_table_col_ref, sample2, info)
        call check_true('random table deterministic', all(sample1 == sample2), failures)
        call check_true('random table nonnegative', all(sample1 >= 0), failures)
        do i = 1, size(sample1,1)
            call check_true('random table row margins', &
                sum(sample1(i,:)) == random_table_row_ref(i), failures)
        end do
        do j = 1, size(sample1,2)
            call check_true('random table col margins', &
                sum(sample1(:,j)) == random_table_col_ref(j), failures)
        end do
        call random_table_rvs_array(state1, random_table_row_ref, random_table_col_ref, samples, info)
        call check_true('random table array status', info == random_table_status_success, failures)
        do j = 1, size(samples,3)
            call check_true('random table array total', &
                sum(samples(:,:,j)) == sum(random_table_row_ref), failures)
        end do
    end subroutine test_random_table

    pure function real_identity(n) result(a)
        integer, intent(in) :: n !! identity matrix dimension
        real(dp) :: a(n,n)

        integer :: i

        a = 0.0_dp
        do i = 1, n
            a(i,i) = 1.0_dp
        end do
    end function real_identity

    pure function complex_identity(n) result(a)
        integer, intent(in) :: n !! identity matrix dimension
        complex(dp) :: a(n,n)

        integer :: i

        a = cmplx(0.0_dp, 0.0_dp, kind=dp)
        do i = 1, n
            a(i,i) = cmplx(1.0_dp, 0.0_dp, kind=dp)
        end do
    end function complex_identity

end program test_directional_random_matrix
