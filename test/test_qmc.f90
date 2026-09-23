program test_qmc
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_random, only : rng_seed, rng_state
    use scifort_stats, only : qmc_discrepancy, qmc_geometric_discrepancy, &
        qmc_halton, qmc_halton_init, qmc_latin_hypercube, qmc_latin_hypercube_init, &
        qmc_multinomial, qmc_multinomial_init, qmc_multivariate_normal, &
        qmc_multivariate_normal_init, qmc_poisson_disk, qmc_poisson_disk_init, &
        qmc_scale, qmc_sobol, qmc_sobol_init, qmc_status_invalid_parameter, &
        qmc_status_out_of_bounds, qmc_status_success, qmc_to_integers, &
        qmc_update_discrepancy, qmc_van_der_corput
    use qmc_reference
    implicit none

    integer :: failures
    integer :: info
    integer(int64) :: integer_sample(3,2)
    integer :: multinomial_sample(4,3)
    integer :: multinomial_zero(3,3)
    integer :: poisson_drawn1
    integer :: poisson_drawn2
    real(dp) :: disc0
    real(dp) :: halton_first(8,3)
    real(dp) :: halton_later(4,3)
    real(dp) :: halton_replay(8,3)
    real(dp) :: lhs1(7,3)
    real(dp) :: lhs2(7,3)
    real(dp) :: oa(9,3)
    real(dp) :: scaled(3,2)
    real(dp) :: reversed(3,2)
    real(dp) :: scale_input(3,2)
    real(dp) :: scrambled1(6,3)
    real(dp) :: scrambled2(6,3)
    real(dp) :: vdc(8)
    real(dp) :: sobol_first(8,3)
    real(dp) :: sobol_later(4,3)
    real(dp) :: sobol_scrambled1(8,3)
    real(dp) :: sobol_scrambled2(8,3)
    real(dp) :: sobol_deep(6,7)
    real(dp) :: sobol_first8_2(8,2)
    real(dp) :: sobol_second8(8,2)
    real(dp) :: sobol_exhaust(1,2)
    real(dp), allocatable :: sobol_high(:, :)
    real(dp) :: mvn_sample(8,2)
    real(dp) :: mvn_singular(4,2)
    real(dp) :: mvn_box1(8,3)
    real(dp) :: mvn_box2(8,3)
    real(dp) :: poisson1(60,2)
    real(dp) :: poisson2(60,2)
    real(dp) :: poisson_negative(40,2)
    real(dp) :: poisson_surface(40,2)
    real(dp) :: poisson_after_fill(2,2)
    real(dp), allocatable :: poisson_fill(:, :)
    type(qmc_halton) :: halton
    type(qmc_halton) :: halton2
    type(qmc_latin_hypercube) :: lhs
    type(qmc_sobol) :: sobol
    type(qmc_sobol) :: sobol2
    type(qmc_sobol) :: sobol_hi
    type(qmc_multinomial) :: multinomial_qmc
    type(qmc_multivariate_normal) :: mvn_qmc
    type(qmc_multivariate_normal) :: mvn_qmc2
    type(qmc_poisson_disk) :: poisson_engine1
    type(qmc_poisson_disk) :: poisson_engine2
    type(rng_state) :: state1
    type(rng_state) :: state2

    failures = 0
    call check_close('CD', qmc_discrepancy(qmc_ref_space, 'CD', status=info), &
        qmc_ref_discrepancy(1), 2.0e-14_dp, failures)
    call check_int('CD status', info, qmc_status_success, failures)
    call check_close('WD', qmc_discrepancy(qmc_ref_space, 'WD'), &
        qmc_ref_discrepancy(2), 2.0e-14_dp, failures)
    call check_close('MD', qmc_discrepancy(qmc_ref_space, 'MD'), &
        qmc_ref_discrepancy(3), 3.0e-14_dp, failures)
    call check_close('L2-star', qmc_discrepancy(qmc_ref_space, 'L2-star'), &
        qmc_ref_discrepancy(4), 2.0e-14_dp, failures)
    disc0 = qmc_discrepancy(qmc_ref_space(1:5,:), 'CD', iterative=.true.)
    call check_close('iterative CD', disc0, qmc_ref_iterative, 2.0e-14_dp, failures)
    call check_close('updated CD', qmc_update_discrepancy(qmc_ref_space(6,:), &
        qmc_ref_space(1:5,:), disc0), qmc_ref_update, 2.0e-14_dp, failures)

    call check_close('geometric mindist', qmc_geometric_discrepancy(qmc_ref_space), &
        qmc_ref_geo(1), 2.0e-14_dp, failures)
    call check_close('geometric mst', qmc_geometric_discrepancy(qmc_ref_space, 'mst'), &
        qmc_ref_geo(2), 2.0e-14_dp, failures)
    call check_close('geometric cityblock', qmc_geometric_discrepancy(qmc_ref_space, &
        'mindist', 'cityblock'), qmc_ref_geo(3), 2.0e-14_dp, failures)

    scale_input(1,:) = [0.5_dp, 0.75_dp]
    scale_input(2,:) = [0.5_dp, 0.5_dp]
    scale_input(3,:) = [0.75_dp, 0.25_dp]
    call qmc_scale(scale_input, [-2.0_dp,0.0_dp], [6.0_dp,5.0_dp], scaled, info)
    call check_int('scale status', info, qmc_status_success, failures)
    call check_mat('scale', scaled, qmc_ref_scale, 0.0_dp, failures)
    call qmc_scale(scaled, [-2.0_dp,0.0_dp], [6.0_dp,5.0_dp], reversed, info, reverse=.true.)
    call check_mat('reverse scale', reversed, scale_input, 0.0_dp, failures)
    call qmc_to_integers(scale_input, [0_int64,2_int64], [10_int64,5_int64], &
        integer_sample, info)
    call check_int('integer map status', info, qmc_status_success, failures)
    call check_int64_mat('integer map', integer_sample, reshape([5_int64,5_int64,7_int64, &
        4_int64,3_int64,2_int64],[3,2]), failures)

    call qmc_van_der_corput(8, 3, vdc, info, start_index=2)
    call check_int('vdc status', info, qmc_status_success, failures)
    call check_vec('vdc', vdc, qmc_ref_vdc, 2.0e-16_dp, failures)

    call qmc_halton_init(halton, 3, info, scramble=.false.)
    call check_int('halton init', info, qmc_status_success, failures)
    call halton%random(8, halton_first, info)
    call check_mat('halton first', halton_first, qmc_ref_halton_first, 2.0e-16_dp, failures)
    call halton%fast_forward(5, info)
    call halton%random(4, halton_later, info)
    call check_mat('halton later', halton_later, qmc_ref_halton_later, 2.0e-16_dp, failures)
    call halton%reset()
    call halton%random(8, halton_replay, info)
    call check_mat('halton reset', halton_replay, qmc_ref_halton_first, 2.0e-16_dp, failures)

    call rng_seed(state1, 24680)
    call rng_seed(state2, 24680)
    call qmc_halton_init(halton, 3, info, scramble=.true., state=state1)
    call qmc_halton_init(halton2, 3, info, scramble=.true., state=state2)
    call halton%random(6, scrambled1, info)
    call halton2%random(6, scrambled2, info)
    call check_mat('scrambled halton replay', scrambled1, scrambled2, 0.0_dp, failures)
    call check_unit_cube('scrambled halton bounds', scrambled1, failures)

    call qmc_sobol_init(sobol, 3, info, scramble=.false., bits=30)
    call check_int('sobol init', info, qmc_status_success, failures)
    call sobol%random_base2(3, sobol_first, info)
    call check_int('sobol base2 status', info, qmc_status_success, failures)
    call check_mat('sobol first', sobol_first, qmc_ref_sobol_first, 0.0_dp, failures)
    call sobol%reset()
    call sobol%fast_forward(4, info)
    call sobol%random_base2(2, sobol_later, info)
    call check_mat('sobol later', sobol_later, qmc_ref_sobol_later, 0.0_dp, failures)

    call rng_seed(state1, 97531)
    call rng_seed(state2, 97531)
    call qmc_sobol_init(sobol, 3, info, scramble=.true., state=state1)
    call qmc_sobol_init(sobol2, 3, info, scramble=.true., state=state2)
    call sobol%random(8, sobol_scrambled1, info)
    call sobol2%random(8, sobol_scrambled2, info)
    call check_mat('scrambled sobol replay', sobol_scrambled1, sobol_scrambled2, 0.0_dp, failures)
    call check_unit_cube('scrambled sobol bounds', sobol_scrambled1, failures)

    allocate(sobol_high(16,21201))
    call qmc_sobol_init(sobol_hi, 21201, info, scramble=.false., bits=64)
    call check_int('sobol maximum dimension init', info, qmc_status_success, failures)
    call sobol_hi%random(16, sobol_high, info)
    call check_vec('sobol maximum dimension data', &
        [sobol_high(16,1), sobol_high(16,100), sobol_high(16,1000), sobol_high(16,21201)], &
        qmc_ref_sobol_hi, 0.0_dp, failures)
    deallocate(sobol_high)

    call qmc_sobol_init(sobol, 7, info, scramble=.false., bits=20)
    call sobol%fast_forward(12345, info)
    call check_int('sobol deep fast-forward', info, qmc_status_success, failures)
    call sobol%random(6, sobol_deep, info)
    call check_mat('sobol deep sequence', sobol_deep, qmc_ref_sobol_deep, 0.0_dp, failures)

    call qmc_sobol_init(sobol, 2, info, scramble=.false., bits=4)
    call sobol%random_base2(3, sobol_first8_2, info)
    call check_int('sobol bits first block', info, qmc_status_success, failures)
    call sobol%random_base2(3, sobol_second8, info)
    call check_int('sobol bits second block', info, qmc_status_success, failures)
    call check_mat('sobol continued base2', sobol_second8, qmc_ref_sobol_second8, 0.0_dp, failures)
    call sobol%random(1, sobol_exhaust, info)
    call check_int('sobol bit exhaustion', info, qmc_status_invalid_parameter, failures)

    call qmc_multinomial_init(multinomial_qmc, [0.2_dp,0.4_dp,0.4_dp], 8, &
        status=info, scramble=.false.)
    call check_int('multinomial QMC init', info, qmc_status_success, failures)
    call multinomial_qmc%random(4, multinomial_sample, info)
    call check_int_mat('multinomial QMC sample', multinomial_sample, qmc_ref_multinomial, failures)
    call qmc_multinomial_init(multinomial_qmc, [0.2_dp,0.3_dp,0.5_dp], 0, &
        status=info, scramble=.false.)
    call multinomial_qmc%random(3, multinomial_zero, info)
    call check_int('multinomial zero-trial status', info, qmc_status_success, failures)
    if (any(multinomial_zero /= 0)) then
        write(*,'(a)') 'multinomial zero-trial mismatch'
        failures = failures + 1
    end if

    call qmc_multivariate_normal_init(mvn_qmc, [1.0_dp,-2.0_dp], info, &
        cov=reshape([2.0_dp,0.5_dp,0.5_dp,1.0_dp],[2,2]), scramble=.false.)
    call check_int('multivariate normal QMC init', info, qmc_status_success, failures)
    call mvn_qmc%random(8, mvn_sample, info)
    call check_mat('multivariate normal QMC', mvn_sample, qmc_ref_mvn, 3.0e-14_dp, failures)
    call qmc_multivariate_normal_init(mvn_qmc, [0.5_dp,-1.0_dp], info, &
        cov=reshape([1.0_dp,1.0_dp,1.0_dp,1.0_dp],[2,2]), scramble=.false.)
    call check_int('singular multivariate normal QMC init', info, qmc_status_success, failures)
    call mvn_qmc%random(4, mvn_singular, info)
    call check_mat('singular multivariate normal QMC', mvn_singular, &
        qmc_ref_mvn_singular, 5.0e-14_dp, failures)

    call rng_seed(state1, 86420)
    call rng_seed(state2, 86420)
    call qmc_multivariate_normal_init(mvn_qmc, [0.0_dp,1.0_dp,-1.0_dp], info, &
        state=state1, inv_transform=.false., scramble=.true.)
    call qmc_multivariate_normal_init(mvn_qmc2, [0.0_dp,1.0_dp,-1.0_dp], info, &
        state=state2, inv_transform=.false., scramble=.true.)
    call mvn_qmc%random(8, mvn_box1, info)
    call mvn_qmc2%random(8, mvn_box2, info)
    call check_mat('Box-Muller multivariate normal QMC replay', mvn_box1, mvn_box2, 0.0_dp, failures)
    if (.not. all(ieee_is_finite(mvn_box1))) then
        write(*,'(a)') 'Box-Muller multivariate normal QMC finite failure'
        failures = failures + 1
    end if

    call qmc_poisson_disk_init(poisson_engine1, 2, info, radius=0.1_dp)
    call qmc_poisson_disk_init(poisson_engine2, 2, info, radius=0.1_dp)
    call rng_seed(state1, 24681357)
    call rng_seed(state2, 24681357)
    call poisson_engine1%random(state1, 60, poisson1, info, poisson_drawn1)
    call check_int('poisson disk status', info, qmc_status_success, failures)
    call poisson_engine2%random(state2, 60, poisson2, info, poisson_drawn2)
    call check_int('poisson disk replay count', poisson_drawn2, poisson_drawn1, failures)
    call check_mat('poisson disk replay', poisson2(1:poisson_drawn2,:), &
        poisson1(1:poisson_drawn1,:), 0.0_dp, failures)
    call check_poisson_disk('poisson disk geometry', poisson1(1:poisson_drawn1,:), &
        0.1_dp, [0.0_dp,0.0_dp], [1.0_dp,1.0_dp], failures)

    call qmc_poisson_disk_init(poisson_engine1, 2, info, radius=0.18_dp, &
        lower=[-2.0_dp,-3.0_dp], upper=[-1.0_dp,-1.0_dp])
    call rng_seed(state1, 34567)
    call poisson_engine1%random(state1, 40, poisson_negative, info, poisson_drawn1)
    call check_int('negative-bound poisson status', info, qmc_status_success, failures)
    call check_poisson_disk('negative-bound poisson geometry', &
        poisson_negative(1:poisson_drawn1,:), 0.18_dp, &
        [-2.0_dp,-3.0_dp], [-1.0_dp,-1.0_dp], failures)

    call qmc_poisson_disk_init(poisson_engine1, 2, info, radius=0.12_dp, hypersphere='surface')
    call rng_seed(state1, 45678)
    call poisson_engine1%random(state1, 40, poisson_surface, info, poisson_drawn1)
    call check_int('surface poisson status', info, qmc_status_success, failures)
    call check_poisson_disk('surface poisson geometry', poisson_surface(1:poisson_drawn1,:), &
        0.12_dp, [0.0_dp,0.0_dp], [1.0_dp,1.0_dp], failures)

    call qmc_poisson_disk_init(poisson_engine1, 2, info, radius=0.28_dp)
    call rng_seed(state1, 56789)
    call poisson_engine1%fill_space(state1, poisson_fill, info)
    call check_int('poisson fill-space status', info, qmc_status_success, failures)
    call check_poisson_disk('poisson fill-space geometry', poisson_fill, 0.28_dp, &
        [0.0_dp,0.0_dp], [1.0_dp,1.0_dp], failures)
    call poisson_engine1%random(state1, 2, poisson_after_fill, info, poisson_drawn1)
    call check_int('poisson exhausted count', poisson_drawn1, 0, failures)

    call qmc_latin_hypercube_init(lhs, 3, info, scramble=.true., strength=1)
    call check_int('lhs init', info, qmc_status_success, failures)
    call rng_seed(state1, 13579)
    call rng_seed(state2, 13579)
    call lhs%random(state1, 7, lhs1, info)
    call check_int('lhs random', info, qmc_status_success, failures)
    call check_latin('lhs strata', lhs1, failures)
    call lhs%reset()
    call lhs%random(state2, 7, lhs2, info)
    call check_mat('lhs replay', lhs1, lhs2, 0.0_dp, failures)

    call qmc_latin_hypercube_init(lhs, 3, info, scramble=.false., strength=2)
    call rng_seed(state1, 112233)
    call lhs%random(state1, 9, oa, info)
    call check_int('OA-LHS status', info, qmc_status_success, failures)
    call check_latin('OA-LHS strata', oa, failures)
    call check_oa_strength_two('OA-LHS strength', oa, 3, failures)

    call invalid_cases(failures)

    if (failures /= 0) then
        write(*,'(a,i0)') 'test_qmc: FAIL ', failures
        error stop 1
    end if
    write(*,'(a)') 'test_qmc: PASS'

contains

    subroutine invalid_cases(failures)
        integer, intent(inout) :: failures !! failure counter

        integer :: info
        real(dp) :: bad(2,2)
        real(dp) :: out(2,2)
        real(dp) :: oa_bad(9,5)
        type(qmc_halton) :: bad_halton
        type(qmc_latin_hypercube) :: bad_lhs
        type(qmc_sobol) :: bad_sobol
        type(qmc_multinomial) :: bad_multinomial
        type(qmc_multivariate_normal) :: bad_mvn
        type(qmc_poisson_disk) :: bad_poisson
        type(rng_state) :: state
        real(dp) :: sobol_three(3,2)
        real(dp) :: sobol_two(2,2)

        bad = reshape([-0.1_dp,0.2_dp,0.3_dp,0.4_dp],[2,2])
        call qmc_scale(bad, [0.0_dp,0.0_dp], [1.0_dp,1.0_dp], out, info)
        call check_int('scale out of bounds', info, qmc_status_out_of_bounds, failures)
        call qmc_halton_init(bad_halton, 2, info, scramble=.true.)
        call check_int('scramble needs state', info, qmc_status_invalid_parameter, failures)
        call qmc_latin_hypercube_init(bad_lhs, 5, info, strength=2)
        call rng_seed(state, 5)
        call bad_lhs%random(state, 9, oa_bad, info)
        call check_int('OA dimension restriction', info, qmc_status_invalid_parameter, failures)

        call qmc_sobol_init(bad_sobol, 21202, info, scramble=.false.)
        call check_int('sobol dimension restriction', info, qmc_status_invalid_parameter, failures)
        call qmc_sobol_init(bad_sobol, 2, info, scramble=.false., bits=65)
        call check_int('sobol bit restriction', info, qmc_status_invalid_parameter, failures)
        call qmc_sobol_init(bad_sobol, 2, info, scramble=.false., bits=8)
        call bad_sobol%random(3, sobol_three, info)
        call bad_sobol%random_base2(1, sobol_two, info)
        call check_int('sobol base2 balance restriction', info, qmc_status_invalid_parameter, failures)

        call qmc_multinomial_init(bad_multinomial, [0.2_dp,0.2_dp], 3, &
            status=info, scramble=.false.)
        call check_int('multinomial probability restriction', info, qmc_status_invalid_parameter, failures)
        call qmc_multivariate_normal_init(bad_mvn, [0.0_dp,0.0_dp], info, &
            cov=reshape([1.0_dp,2.0_dp,2.0_dp,1.0_dp],[2,2]), scramble=.false.)
        call check_int('QMC normal PSD restriction', info, qmc_status_invalid_parameter, failures)
        call qmc_poisson_disk_init(bad_poisson, 2, info, radius=0.1_dp, &
            lower=[0.0_dp,1.0_dp], upper=[1.0_dp,0.0_dp])
        call check_int('poisson bounds restriction', info, qmc_status_invalid_parameter, failures)
    end subroutine invalid_cases


    subroutine check_latin(name, sample, failures)
        character(len=*), intent(in) :: name !! test label
        real(dp), intent(in) :: sample(:, :) !! candidate Latin-hypercube sample
        integer, intent(inout) :: failures !! failure counter

        integer :: i
        integer :: j
        integer :: n
        integer, allocatable :: counts(:)

        n = size(sample,1)
        allocate(counts(0:n-1))
        do j = 1, size(sample,2)
            counts = 0
            do i = 1, n
                if (sample(i,j) < 0.0_dp .or. sample(i,j) >= 1.0_dp) then
                    failures = failures + 1
                    write(*,'(a)') trim(name)//' out of bounds'
                    return
                end if
                counts(min(n-1,int(floor(real(n,dp)*sample(i,j))))) = &
                    counts(min(n-1,int(floor(real(n,dp)*sample(i,j))))) + 1
            end do
            if (any(counts /= 1)) then
                failures = failures + 1
                write(*,'(a)') trim(name)//' missing/repeated stratum'
                return
            end if
        end do
    end subroutine check_latin

    subroutine check_oa_strength_two(name, sample, p, failures)
        character(len=*), intent(in) :: name !! test label
        real(dp), intent(in) :: sample(:, :) !! candidate strength-two OA-LHS sample
        integer, intent(in) :: p !! prime number of coarse symbols per coordinate
        integer, intent(inout) :: failures !! failure counter

        integer :: a
        integer :: b
        integer :: i
        integer :: j
        integer :: k
        integer :: counts(0:p-1,0:p-1)

        do j = 1, size(sample,2) - 1
            do k = j + 1, size(sample,2)
                counts = 0
                do i = 1, size(sample,1)
                    a = min(p-1,int(floor(real(p,dp)*sample(i,j))))
                    b = min(p-1,int(floor(real(p,dp)*sample(i,k))))
                    counts(a,b) = counts(a,b) + 1
                end do
                if (any(counts /= 1)) then
                    failures = failures + 1
                    write(*,'(a,2i4)') trim(name)//' failed columns ', j, k
                    return
                end if
            end do
        end do
    end subroutine check_oa_strength_two

    subroutine check_unit_cube(name, sample, failures)
        character(len=*), intent(in) :: name !! test label
        real(dp), intent(in) :: sample(:, :) !! sample expected in the unit cube
        integer, intent(inout) :: failures !! failure counter

        if (any(sample < 0.0_dp) .or. any(sample >= 1.0_dp)) then
            failures = failures + 1
            write(*,'(a)') trim(name)//' failed'
        end if
    end subroutine check_unit_cube

    subroutine check_close(name, got, expected, tol, failures)
        character(len=*), intent(in) :: name !! test label
        real(dp), intent(in) :: got !! observed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: tol !! absolute tolerance
        integer, intent(inout) :: failures !! failure counter

        if (abs(got - expected) > tol) then
            write(*,'(a,2es24.15,es12.3)') trim(name)//' mismatch: ', got, expected, abs(got-expected)
            failures = failures + 1
        end if
    end subroutine check_close

    subroutine check_vec(name, got, expected, tol, failures)
        character(len=*), intent(in) :: name !! test label
        real(dp), intent(in) :: got(:) !! observed vector
        real(dp), intent(in) :: expected(:) !! reference vector
        real(dp), intent(in) :: tol !! maximum absolute tolerance
        integer, intent(inout) :: failures !! failure counter

        if (size(got) /= size(expected) .or. maxval(abs(got-expected)) > tol) then
            write(*,'(a,es12.3)') trim(name)//' max error: ', maxval(abs(got-expected))
            failures = failures + 1
        end if
    end subroutine check_vec

    subroutine check_mat(name, got, expected, tol, failures)
        character(len=*), intent(in) :: name !! test label
        real(dp), intent(in) :: got(:, :) !! observed matrix
        real(dp), intent(in) :: expected(:, :) !! reference matrix
        real(dp), intent(in) :: tol !! maximum absolute tolerance
        integer, intent(inout) :: failures !! failure counter

        if (any(shape(got) /= shape(expected)) .or. maxval(abs(got-expected)) > tol) then
            write(*,'(a,es12.3)') trim(name)//' max error: ', maxval(abs(got-expected))
            failures = failures + 1
        end if
    end subroutine check_mat

    subroutine check_int(name, got, expected, failures)
        character(len=*), intent(in) :: name !! test label
        integer, intent(in) :: got !! observed integer
        integer, intent(in) :: expected !! expected integer
        integer, intent(inout) :: failures !! failure counter

        if (got /= expected) then
            write(*,'(a,2i8)') trim(name)//' mismatch: ', got, expected
            failures = failures + 1
        end if
    end subroutine check_int

    subroutine check_int_mat(name, got, expected, failures)
        character(len=*), intent(in) :: name !! test label
        integer, intent(in) :: got(:, :) !! observed integer matrix
        integer, intent(in) :: expected(:, :) !! expected integer matrix
        integer, intent(inout) :: failures !! failure counter

        if (any(shape(got) /= shape(expected)) .or. any(got /= expected)) then
            write(*,'(a)') trim(name)//' mismatch'
            failures = failures + 1
        end if
    end subroutine check_int_mat

    subroutine check_poisson_disk(name, sample, radius, lower, upper, failures)
        character(len=*), intent(in) :: name !! test label
        real(dp), intent(in) :: sample(:, :) !! candidate Poisson-disk sample
        real(dp), intent(in) :: radius !! required minimum Euclidean distance
        real(dp), intent(in) :: lower(:) !! lower coordinate bounds
        real(dp), intent(in) :: upper(:) !! upper coordinate bounds
        integer, intent(inout) :: failures !! failure counter

        integer :: i
        integer :: j

        do i = 1, size(sample,1)
            if (any(sample(i,:) < lower) .or. any(sample(i,:) > upper)) then
                write(*,'(a)') trim(name)//' bounds failure'
                failures = failures + 1
                return
            end if
            do j = 1, i - 1
                if (sqrt(sum((sample(i,:) - sample(j,:)) ** 2)) < radius - 1.0e-13_dp) then
                    write(*,'(a)') trim(name)//' separation failure'
                    failures = failures + 1
                    return
                end if
            end do
        end do
    end subroutine check_poisson_disk

    subroutine check_int64_mat(name, got, expected, failures)
        character(len=*), intent(in) :: name !! test label
        integer(int64), intent(in) :: got(:, :) !! observed integer matrix
        integer(int64), intent(in) :: expected(:, :) !! expected integer matrix
        integer, intent(inout) :: failures !! failure counter

        if (any(shape(got) /= shape(expected)) .or. any(got /= expected)) then
            write(*,'(a)') trim(name)//' mismatch'
            failures = failures + 1
        end if
    end subroutine check_int64_mat

end program test_qmc
