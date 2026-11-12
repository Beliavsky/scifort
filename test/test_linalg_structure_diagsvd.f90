! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_linalg_structure_diagsvd
    use, intrinsic :: iso_fortran_env, only : int32, int64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    use scifort_linalg_structure, only : bandwidth, issymmetric, ishermitian
    use scifort_linalg_diagsvd, only : diagsvd
    implicit none
    integer :: lower, upper, status, i
    integer :: lower_batch(2), upper_batch(2), lower_bad(1)
    logical :: is_symmetric, is_hermitian, flags(2), bad_flags(1)
    real(dp) :: a(5,5), empty(0,0), nonsquare(2,3), nanval, infval
    real(dp) :: batch(2,3,3), rsv(3), sigma(3,5), zerosigma(0,3)
    complex(dp) :: c(2,2), cs(2), cds(2,3), cd_batch(2,2), cd_out(2,2,3)
    integer(int64) :: big(2,2), bigsingular(2), bigdiag(2,3)
    integer(int32) :: ids(2), diag32(3,2)
    real(dp) :: large(2,2)

    a = 0.0_dp
    do i = 1, 5
        a(i,i) = real(i, dp)
    end do
    a(4,1) = 8.0_dp
    a(5,2) = 9.0_dp
    a(3,4) = 1.0_dp
    call bandwidth(a, lower, upper, status)
    call check(status == 0 .and. lower == 3 .and. upper == 1, 'real bandwidth')

    call bandwidth(empty, lower, upper, status)
    call check(status == 0 .and. lower == 0 .and. upper == 0, 'empty bandwidth')
    call issymmetric(empty, is_symmetric, status)
    call check(status == 0 .and. is_symmetric, 'empty symmetric')
    call ishermitian(empty, is_hermitian, status)
    call check(status == 0 .and. is_hermitian, 'empty Hermitian')

    nonsquare = 0.0_dp
    call issymmetric(nonsquare, is_symmetric, status)
    call check(status == linalg_status_invalid_shape .and. .not. is_symmetric, 'nonsquare symmetry')
    call ishermitian(nonsquare, is_hermitian, status)
    call check(status == linalg_status_invalid_shape .and. .not. is_hermitian, 'nonsquare Hermitian')

    nanval = ieee_value(0.0_dp, ieee_quiet_nan)
    infval = ieee_value(0.0_dp, ieee_positive_inf)
    large = reshape([nanval, 1.0_dp, 1.0_dp, 2.0_dp], [2,2])
    call issymmetric(large, is_symmetric, status)
    call check(status == 0 .and. is_symmetric, 'ignore diagonal NaN in exact symmetry')
    call issymmetric(large, is_symmetric, status, atol=1.0e-6_dp)
    call check(status == 0 .and. .not. is_symmetric, 'compare NaN diagonal with tolerance')
    large = reshape([1.0_dp, infval, infval, 2.0_dp], [2,2])
    call ishermitian(large, is_hermitian, status)
    call check(status == 0 .and. is_hermitian, 'matching infinite offdiagonals')
    large(1,2) = nanval
    call bandwidth(large, lower, upper, status)
    call check(status == 0 .and. lower == 1 .and. upper == 1, 'NaN is a nonzero band element')
    call issymmetric(large, is_symmetric, status)
    call check(status == 0 .and. .not. is_symmetric, 'offdiagonal NaN differs')

    c = cmplx(0.0_dp, 0.0_dp, dp)
    c(1,1) = cmplx(1.0_dp, 2.0_dp, dp)
    c(2,2) = cmplx(3.0_dp, 0.0_dp, dp)
    c(1,2) = cmplx(0.0_dp, 5.0_dp, dp)
    c(2,1) = c(1,2)
    call issymmetric(c, is_symmetric, status)
    call check(status == 0 .and. is_symmetric, 'complex symmetry is not Hermitian')
    call ishermitian(c, is_hermitian, status)
    call check(status == 0 .and. .not. is_hermitian, 'complex diagonal imaginary part')
    c(1,1) = cmplx(1.0_dp, 0.0_dp, dp)
    c(2,1) = conjg(c(1,2))
    call ishermitian(c, is_hermitian, status)
    call check(status == 0 .and. is_hermitian, 'complex Hermitian')
    call issymmetric(c, is_symmetric, status)
    call check(status == 0 .and. .not. is_symmetric, 'Hermitian not necessarily symmetric')
    c(2,1) = c(2,1) + cmplx(1.0e-9_dp, 0.0_dp, dp)
    call ishermitian(c, is_hermitian, status, atol=1.0e-8_dp)
    call check(status == 0 .and. is_hermitian, 'Hermitian within tolerance')
    call ishermitian(c, is_hermitian, status, atol=1.0e-10_dp)
    call check(status == 0 .and. .not. is_hermitian, 'Hermitian outside tolerance')
    call ishermitian(c, is_hermitian, status, atol=-1.0_dp)
    call check(status == linalg_status_invalid_shape, 'reject negative tolerance')

    big = 0_int64
    big(1,1) = 9007199254740993_int64
    big(2,2) = 17_int64
    big(1,2) = 9007199254740993_int64
    big(2,1) = 9007199254740994_int64
    call issymmetric(big, is_symmetric, status)
    call check(status == 0 .and. .not. is_symmetric, 'exact int64 above 2**53')
    call issymmetric(big, is_symmetric, status, atol=0.1_dp)
    call check(status == 0 .and. .not. is_symmetric, 'int64 tolerant comparison retains low bits')
    call issymmetric(big, is_symmetric, status, atol=1.0_dp)
    call check(status == 0 .and. .not. is_symmetric, 'integer tolerance ignored as in SciPy')
    call bandwidth(big, lower, upper, status)
    call check(status == 0 .and. lower == 1 .and. upper == 1, 'int64 bandwidth')

    batch = 0.0_dp
    batch(1,3,1) = 2.0_dp
    batch(2,1,3) = 3.0_dp
    call bandwidth(batch, lower_batch, upper_batch, status)
    call check(status == 0, 'batch bandwidth status')
    call check(all(lower_batch == [2,0]) .and. all(upper_batch == [0,2]), 'batch bandwidth values')
    call bandwidth(batch, lower_bad, upper_batch, status)
    call check(status == linalg_status_invalid_shape, 'batch bandwidth shape check')
    call issymmetric(batch, flags, status)
    call check(status == 0 .and. all(.not. flags), 'batch symmetry values')
    batch(1,1,3) = 2.0_dp
    batch(2,3,1) = 3.0_dp
    call ishermitian(batch, flags, status)
    call check(status == 0 .and. all(flags), 'batched Hermitian values')
    call issymmetric(batch, bad_flags, status)
    call check(status == linalg_status_invalid_shape, 'batched symmetry shape check')

    rsv = [3.0_dp, 2.0_dp, 1.0_dp]
    call diagsvd(rsv, 3, 5, sigma, status)
    call check(status == 0, 'real diagsvd shape')
    call check(all([sigma(1,1),sigma(2,2),sigma(3,3)] == rsv), 'real diagsvd diagonal')
    call check(all(sigma(:,4:5) == 0.0_dp), 'real diagsvd padding')
    call diagsvd(rsv, 0, 3, zerosigma, status)
    call check(status == linalg_status_invalid_shape, 'singular value length mismatch')

    ids = [3_int32, -5_int32]
    diag32 = 0_int32
    call diagsvd(ids, 3, 2, diag32, status)
    call check(status == 0, 'int32 diagsvd')
    call check(diag32(1,1) == 3_int32 .and. diag32(2,2) == -5_int32, 'int32 diagonals')
    bigsingular = [9007199254740993_int64, -9007199254740995_int64]
    call diagsvd(bigsingular, 2, 3, bigdiag, status)
    call check(status == 0 .and. bigdiag(1,1) == bigsingular(1), 'exact int64 singular values')
    call check(bigdiag(2,2) == bigsingular(2), 'exact int64 second singular value')

    cs = [cmplx(3.0_dp,1.0_dp,dp), cmplx(-4.0_dp,2.0_dp,dp)]
    call diagsvd(cs, 2, 3, cds, status)
    call check(status == 0 .and. cds(1,1) == cs(1) .and. cds(2,2) == cs(2), &
        'complex diagsvd values')
    cd_batch(1,:) = cs
    cd_batch(2,:) = 2.0_dp*cs
    call diagsvd(cd_batch, 2, 3, cd_out, status)
    call check(status == 0 .and. cd_out(2,1,1) == 2.0_dp*cs(1), 'batched complex diagsvd')

    print *, 'SciFort structure and diagsvd native tests passed'
contains
    subroutine check(ok, label)
        logical, intent(in) :: ok !! asserted condition
        character(*), intent(in) :: label !! assertion label
        if (.not. ok) then
            print *, 'FAILED: ', label
            error stop 1
        end if
    end subroutine check
end program test_linalg_structure_diagsvd
