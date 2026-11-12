! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Additional scipy.sparse-style constructors.  Kept separate from scifort_sparse
! so deterministic random generation can reuse scifort_random without making
! the base CSR/CSC module depend on the random-number subsystem.
module scifort_sparse_constructors
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_random, only : rng_state, rng_uniform
    use scifort_sparse, only : csr_matrix, csr_diags, csr_from_coo, sparse_success, sparse_invalid_input
    implicit none
    private

    public :: sparse_diagonal, csr_diags_array, csr_random, csr_random_array, csr_rand

    interface csr_random_array
        module procedure csr_random
    end interface csr_random_array

    ! A ragged list of unpadded diagonal values, compatible with
    ! sparse_diagonal([1.0_dp, 2.0_dp, ...]) structure constructors.
    type :: sparse_diagonal
        real(dp), allocatable :: values(:)
    end type sparse_diagonal

contains

    function csr_diags_array(diagonals, offsets, n_rows, n_cols, status) result(a)
        type(sparse_diagonal), intent(in) :: diagonals(:) !! Ragged, unpadded diagonal arrays.
        integer, intent(in) :: offsets(:) !! Column-minus-row offsets; must be distinct.
        integer, intent(in), optional :: n_rows !! Explicit number of rows; requires n_cols.
        integer, intent(in), optional :: n_cols !! Explicit number of columns; requires n_rows.
        integer, intent(out), optional :: status !! sparse_success or sparse_invalid_input.
        type(csr_matrix) :: a
        integer, allocatable :: pointers(:)
        real(dp), allocatable :: packed(:)
        integer(int64) :: first, last, count64, total64, inferred
        integer :: m, n, t, cursor, count_t, ierr
        logical :: has_shape

        if (present(status)) status = sparse_invalid_input
        if (size(diagonals) /= size(offsets) .or. size(diagonals) == 0) return
        if (present(n_rows) .neqv. present(n_cols)) return
        has_shape = present(n_rows)
        if (has_shape) then
            m = n_rows
            n = n_cols
        else
            if (.not. allocated(diagonals(1)%values)) return
            inferred = int(size(diagonals(1)%values),int64) + abs(int(offsets(1),int64))
            if (inferred >= int(huge(0),int64)) return
            m = int(inferred)
            n = m
        end if
        if (m < 0 .or. m == huge(0) .or. n < 0) return
        total64 = 0_int64
        do t = 1, size(offsets)
            if (.not. allocated(diagonals(t)%values)) return
            if (.not. all(ieee_is_finite(diagonals(t)%values))) return
            if (any(offsets(:t-1) == offsets(t))) return
            ! For nonempty axes, match SciPy's offsets in [-m, n].
            if (m > 0 .and. n > 0) then
                if (int(offsets(t),int64) < -int(m,int64) .or. &
                    int(offsets(t),int64) > int(n,int64)) return
            end if
            first = max(1_int64,1_int64-int(offsets(t),int64))
            last = min(int(m,int64),int(n,int64)-int(offsets(t),int64))
            count64 = max(0_int64,last-first+1_int64)
            count_t = size(diagonals(t)%values)
            if (count64 > 0_int64) then
                if (int(count_t,int64) < count64) then
                    if (.not. (has_shape .and. count_t == 1)) return
                end if
            end if
            total64 = total64 + count64
            if (total64 >= int(huge(0),int64)) return
        end do

        allocate(pointers(size(offsets)+1),packed(int(total64)))
        cursor = 1
        do t = 1, size(offsets)
            pointers(t) = cursor
            first = max(1_int64,1_int64-int(offsets(t),int64))
            last = min(int(m,int64),int(n,int64)-int(offsets(t),int64))
            count_t = int(max(0_int64,last-first+1_int64))
            if (count_t == 0) cycle
            if (size(diagonals(t)%values) == 1 .and. has_shape) then
                packed(cursor:cursor+count_t-1) = diagonals(t)%values(1)
            else
                packed(cursor:cursor+count_t-1) = diagonals(t)%values(:count_t)
            end if
            cursor = cursor+count_t
        end do
        pointers(size(offsets)+1) = cursor
        a = csr_diags(m,n,offsets,pointers,packed,ierr)
        if (present(status)) status = ierr
    end function csr_diags_array

    ! Samples exactly NINT(density * n_rows * n_cols) DISTINCT coordinates,
    ! with no O(n_rows*n_cols) work array.  A copy of the generator is advanced
    ! first, so invalid arguments or construction failures do not alter state.
    ! The PRNG is SciFort's MRG32k3a; it does not duplicate NumPy's bit stream.
    function csr_random(n_rows, n_cols, density, state, status, distribution, data_values) result(a)
        integer, intent(in) :: n_rows !! Matrix rows, nonnegative.
        integer, intent(in) :: n_cols !! Matrix columns, nonnegative.
        real(dp), intent(in) :: density !! Fraction of stored positions in [0, 1].
        type(rng_state), intent(inout) :: state !! Explicit RNG, committed only on success.
        integer, intent(out), optional :: status !! sparse_success or sparse_invalid_input.
        character(len=*), intent(in), optional :: distribution !! 'uniform' (default) or 'normal'.
        real(dp), intent(in), optional :: data_values(:) !! Optional explicit stored values, size nnz.
        type(csr_matrix) :: a
        type(rng_state) :: work
        integer(int64), allocatable :: choices(:), table(:)
        integer, allocatable :: rows(:), columns(:)
        real(dp), allocatable :: values(:)
        integer(int64) :: total, target, j, draw, chosen, capacity, whole, used
        real(dp) :: exact, fractional, u, v
        character(len=:), allocatable :: mode
        integer :: i, k, ierr

        if (present(status)) status = sparse_invalid_input
        if (n_rows < 0 .or. n_rows == huge(0) .or. n_cols < 0) return
        if (.not. ieee_is_finite(density)) return
        if (density < 0.0_dp .or. density > 1.0_dp) return
        mode = 'uniform'
        if (present(distribution)) mode = trim(adjustl(distribution))
        if (mode /= 'uniform' .and. mode /= 'normal') return
        total = int(n_rows,int64)*int(n_cols,int64)
        if (density >= 1.0_dp) then
            target = total
        else
            exact = density * real(total,dp)
            if (.not. ieee_is_finite(exact)) return
            if (exact >= real(huge(0),dp)) return
            whole = int(exact,int64)
            fractional = exact-real(whole,dp)
            target = whole
            ! Python round() used by scipy.sparse.random: ties to even.
            if (fractional >= 0.5_dp) then
                if (fractional > 0.5_dp .or. mod(whole,2_int64) /= 0) target = target+1_int64
            end if
        end if
        if (target >= int(huge(0),int64)) return
        if (present(data_values)) then
            if (size(data_values,kind=int64) /= target) return
            if (.not. all(ieee_is_finite(data_values))) return
        end if
        if (target /= total) then
            ! Open-addressed set for Floyd sampling, O(nnz) expected memory.
            if (target > int((huge(0)-3)/2,int64)) return
        end if

        allocate(choices(int(target)),rows(int(target)),columns(int(target)),values(int(target)))
        work = state
        if (target == total) then
            do i = 1, int(target)
                choices(i) = int(i,int64)-1_int64
            end do
        else if (target > 0_int64) then
            capacity = 2_int64*target+1_int64
            allocate(table(int(capacity)))
            table = -1_int64
            used = 0_int64
            do j = total-target, total-1_int64
                ! rng_uniform is in (0,1); guard rounding to j+1.
                draw = sparse_random_below(work,j+1_int64)
                if (sparse_hash_contains(table,draw)) then
                    chosen = j
                else
                    chosen = draw
                end if
                call sparse_hash_insert(table,chosen)
                used = used+1_int64
                choices(int(used)) = chosen
            end do
        end if

        do k = 1, int(target)
            rows(k) = int(choices(k)/int(n_cols,int64))+1
            columns(k) = int(mod(choices(k),int(n_cols,int64)))+1
            if (present(data_values)) then
                values(k) = data_values(k)
            else if (mode == 'normal') then
                u = rng_uniform(work)
                v = rng_uniform(work)
                values(k) = sqrt(-2.0_dp*log(u))*cos(2.0_dp*acos(-1.0_dp)*v)
            else
                values(k) = rng_uniform(work)
            end if
        end do
        a = csr_from_coo(n_rows,n_cols,rows,columns,values,ierr)
        if (ierr /= sparse_success) return
        state = work
        if (present(status)) status = sparse_success
    end function csr_random

    function csr_rand(n_rows, n_cols, density, state, status) result(a)
        integer, intent(in) :: n_rows !! Matrix rows, nonnegative.
        integer, intent(in) :: n_cols !! Matrix columns, nonnegative.
        real(dp), intent(in) :: density !! Fraction of stored positions in [0, 1].
        type(rng_state), intent(inout) :: state !! Explicit SciFort RNG.
        integer, intent(out), optional :: status !! sparse_success or sparse_invalid_input.
        type(csr_matrix) :: a
        integer :: ierr

        a = csr_random(n_rows,n_cols,density,state,ierr)
        if (present(status)) status = ierr
    end function csr_rand

    ! A bounded integer draw uses two independent 31-bit chunks. This avoids
    ! the 52-bit address-range limitation of converting a single real uniform
    ! directly to a coordinate when the logical matrix is exceptionally large.
    ! Every product of two default integer matrix dimensions is below 2**62.
    function sparse_random_below(state, bound) result(value)
        type(rng_state), intent(inout) :: state !! RNG advanced by two or more uniforms.
        integer(int64), intent(in) :: bound !! Positive exclusive upper bound, <= 2**62.
        integer(int64) :: value
        integer(int64), parameter :: base = 2147483648_int64
        integer(int64), parameter :: span = base*base
        integer(int64) :: high, low, bits, accepted

        accepted = (span/bound)*bound
        do
            high = min(base-1_int64,int(rng_uniform(state)*real(base,dp),int64))
            low = min(base-1_int64,int(rng_uniform(state)*real(base,dp),int64))
            bits = high*base+low
            if (bits < accepted) exit
        end do
        value = mod(bits,bound)
    end function sparse_random_below

    pure logical function sparse_hash_contains(table, key) result(found)
        integer(int64), intent(in) :: table(:) !! Open-addressed set; -1 is empty.
        integer(int64), intent(in) :: key !! Nonnegative linear coordinate.
        integer :: slot, capacity

        capacity = size(table)
        slot = int(modulo(key,int(capacity,int64)))+1
        do
            if (table(slot) == -1_int64) then
                found = .false.
                return
            end if
            if (table(slot) == key) then
                found = .true.
                return
            end if
            slot = slot+1
            if (slot > capacity) slot = 1
        end do
    end function sparse_hash_contains

    subroutine sparse_hash_insert(table, key)
        integer(int64), intent(inout) :: table(:) !! Hash set to extend.
        integer(int64), intent(in) :: key !! A distinct nonnegative key.
        integer :: slot, capacity

        capacity = size(table)
        slot = int(modulo(key,int(capacity,int64)))+1
        do
            if (table(slot) == -1_int64) then
                table(slot) = key
                return
            end if
            if (table(slot) == key) return
            slot = slot+1
            if (slot > capacity) slot = 1
        end do
    end subroutine sparse_hash_insert
end module scifort_sparse_constructors
