! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Minimal compressed-sparse-column containers for solve_ivp Jacobians.
! Constant, callback, and sparsity-structured numerical CSC Jacobians are factored
! by scifort_integrate_sparse_lu without assembling a dense Newton matrix.

module scifort_integrate_sparse_jacobian
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    type, public :: ivp_sparse_matrix
        integer :: nrows = 0
        integer :: ncols = 0
        integer, allocatable :: indptr(:)
        integer, allocatable :: indices(:)
        real(dp), allocatable :: data(:)
    contains
        procedure :: is_valid => real_sparse_is_valid
        procedure :: nnz => real_sparse_nnz
    end type ivp_sparse_matrix

    type, public :: ivp_complex_sparse_matrix
        integer :: nrows = 0
        integer :: ncols = 0
        integer, allocatable :: indptr(:)
        integer, allocatable :: indices(:)
        complex(dp), allocatable :: data(:)
    contains
        procedure :: is_valid => complex_sparse_is_valid
        procedure :: nnz => complex_sparse_nnz
    end type ivp_complex_sparse_matrix

    interface ivp_csc_matrix
        module procedure make_real_sparse_matrix
        module procedure make_complex_sparse_matrix
    end interface ivp_csc_matrix

    interface ivp_csc_from_pattern
        module procedure initialize_real_sparse_pattern
        module procedure initialize_complex_sparse_pattern
    end interface ivp_csc_from_pattern

    public :: ivp_csc_matrix
    public :: ivp_csc_from_pattern

contains

    function make_real_sparse_matrix(nrows, ncols, indptr, indices, data) result(matrix)
        integer, intent(in) :: nrows !! number of matrix rows, nonnegative
        integer, intent(in) :: ncols !! number of matrix columns, nonnegative
        integer, intent(in) :: indptr(:) !! one-based CSC column pointers, length ncols+1
        integer, intent(in) :: indices(:) !! one-based CSC row indices
        real(dp), intent(in) :: data(:) !! nonzero values corresponding to indices
        type(ivp_sparse_matrix) :: matrix

        matrix%nrows = nrows
        matrix%ncols = ncols
        allocate(matrix%indptr(size(indptr)))
        allocate(matrix%indices(size(indices)))
        allocate(matrix%data(size(data)))
        matrix%indptr = indptr
        matrix%indices = indices
        matrix%data = data
    end function make_real_sparse_matrix

    function make_complex_sparse_matrix(nrows, ncols, indptr, indices, data) result(matrix)
        integer, intent(in) :: nrows !! number of matrix rows, nonnegative
        integer, intent(in) :: ncols !! number of matrix columns, nonnegative
        integer, intent(in) :: indptr(:) !! one-based CSC column pointers, length ncols+1
        integer, intent(in) :: indices(:) !! one-based CSC row indices
        complex(dp), intent(in) :: data(:) !! nonzero values corresponding to indices
        type(ivp_complex_sparse_matrix) :: matrix

        matrix%nrows = nrows
        matrix%ncols = ncols
        allocate(matrix%indptr(size(indptr)))
        allocate(matrix%indices(size(indices)))
        allocate(matrix%data(size(data)))
        matrix%indptr = indptr
        matrix%indices = indices
        matrix%data = data
    end function make_complex_sparse_matrix


    subroutine initialize_real_sparse_pattern(sparsity, matrix)
        logical, intent(in) :: sparsity(:, :) !! structural nonzero pattern, shape (nrows,ncols)
        type(ivp_sparse_matrix), intent(out) :: matrix !! zero-valued one-based CSC matrix

        integer :: j
        integer :: k
        integer :: nnz

        matrix%nrows = size(sparsity, 1)
        matrix%ncols = size(sparsity, 2)
        nnz = count(sparsity)
        allocate(matrix%indptr(matrix%ncols + 1), matrix%indices(nnz), matrix%data(nnz))
        matrix%indptr(1) = 1
        k = 1
        do j = 1, matrix%ncols
            call append_pattern_column(sparsity(:, j), matrix%indices, k)
            matrix%indptr(j + 1) = k
        end do
        matrix%data = 0.0_dp
    end subroutine initialize_real_sparse_pattern

    subroutine initialize_complex_sparse_pattern(sparsity, matrix)
        logical, intent(in) :: sparsity(:, :) !! structural nonzero pattern, shape (nrows,ncols)
        type(ivp_complex_sparse_matrix), intent(out) :: matrix !! zero-valued one-based CSC matrix

        integer :: j
        integer :: k
        integer :: nnz

        matrix%nrows = size(sparsity, 1)
        matrix%ncols = size(sparsity, 2)
        nnz = count(sparsity)
        allocate(matrix%indptr(matrix%ncols + 1), matrix%indices(nnz), matrix%data(nnz))
        matrix%indptr(1) = 1
        k = 1
        do j = 1, matrix%ncols
            call append_pattern_column(sparsity(:, j), matrix%indices, k)
            matrix%indptr(j + 1) = k
        end do
        matrix%data = cmplx(0.0_dp, 0.0_dp, kind=dp)
    end subroutine initialize_complex_sparse_pattern

    pure subroutine append_pattern_column(column, indices, k)
        logical, intent(in) :: column(:) !! structural nonzeros for one matrix column
        integer, intent(inout) :: indices(:) !! CSC row-index storage being filled
        integer, intent(inout) :: k !! next one-based position in indices

        integer :: i

        do i = 1, size(column)
            if (.not. column(i)) cycle
            indices(k) = i
            k = k + 1
        end do
    end subroutine append_pattern_column

    logical function real_sparse_is_valid(self) result(valid)
        class(ivp_sparse_matrix), intent(in) :: self !! CSC matrix to validate

        valid = sparse_structure_is_valid(self%nrows, self%ncols, self%indptr, self%indices)
        if (.not. valid) return
        valid = allocated(self%data)
        if (.not. valid) return
        valid = size(self%data) == size(self%indices)
        if (.not. valid) return
        valid = all(ieee_is_finite(self%data))
    end function real_sparse_is_valid

    logical function complex_sparse_is_valid(self) result(valid)
        class(ivp_complex_sparse_matrix), intent(in) :: self !! complex CSC matrix to validate

        valid = sparse_structure_is_valid(self%nrows, self%ncols, self%indptr, self%indices)
        if (.not. valid) return
        valid = allocated(self%data)
        if (.not. valid) return
        valid = size(self%data) == size(self%indices)
        if (.not. valid) return
        valid = all(ieee_is_finite(real(self%data, dp))) .and. &
            all(ieee_is_finite(aimag(self%data)))
    end function complex_sparse_is_valid

    integer function real_sparse_nnz(self) result(value)
        class(ivp_sparse_matrix), intent(in) :: self !! CSC matrix whose stored count is requested

        if (allocated(self%data)) then
            value = size(self%data)
        else
            value = 0
        end if
    end function real_sparse_nnz

    integer function complex_sparse_nnz(self) result(value)
        class(ivp_complex_sparse_matrix), intent(in) :: &
            self !! complex CSC matrix whose stored count is requested

        if (allocated(self%data)) then
            value = size(self%data)
        else
            value = 0
        end if
    end function complex_sparse_nnz

    logical function sparse_structure_is_valid(nrows, ncols, indptr, indices) result(valid)
        integer, intent(in) :: nrows !! declared row count, nonnegative
        integer, intent(in) :: ncols !! declared column count, nonnegative
        integer, allocatable, intent(in) :: indptr(:) !! one-based CSC column pointers
        integer, allocatable, intent(in) :: indices(:) !! one-based CSC row indices

        integer :: j
        integer :: k
        integer :: nnz

        valid = .false.
        if (nrows < 0 .or. ncols < 0) return
        if (.not. allocated(indptr) .or. .not. allocated(indices)) return
        if (size(indptr) /= ncols + 1) return
        nnz = size(indices)
        if (indptr(1) /= 1 .or. indptr(ncols + 1) /= nnz + 1) return
        if (any(indptr < 1) .or. any(indptr > nnz + 1)) return
        do j = 1, ncols
            if (indptr(j + 1) < indptr(j)) return
            do k = indptr(j), indptr(j + 1) - 1
                if (indices(k) < 1 .or. indices(k) > nrows) return
                if (k > indptr(j)) then
                    if (indices(k) <= indices(k - 1)) return
                end if
            end do
        end do
        valid = .true.
    end function sparse_structure_is_valid

end module scifort_integrate_sparse_jacobian
