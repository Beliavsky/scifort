! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Grid Laplacian matching scipy.sparse.linalg.LaplacianNd for the
! mathematical (matrix-free) operator. Grid coordinates use C ordering:
! the final dimension varies fastest in a one-based Fortran vector.
! All boundary conditions use the standard second difference, with
! negative-semidefinite Neumann and periodic operators.
module scifort_sparse_laplacian_nd
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_coo, sparse_success
    implicit none
    private

    integer, parameter, public :: laplacian_nd_success = 0
    integer, parameter, public :: laplacian_nd_invalid_input = 1
    integer, parameter :: bc_dirichlet = 1, bc_neumann = 2, bc_periodic = 3

    type, public :: laplacian_nd_operator
        private
        integer :: n = 0
        integer :: boundary = 0
        integer, allocatable :: grid_shape(:)
        integer, allocatable :: strides(:)
    end type laplacian_nd_operator

    public :: laplacian_nd, valid_laplacian_nd, laplacian_nd_shape
    public :: laplacian_nd_matvec, laplacian_nd_matmat, laplacian_nd_tocsr
    public :: laplacian_nd_eigenvalues, laplacian_nd_eigenvectors

contains

    function laplacian_nd(grid_shape, boundary_conditions, status) result(op)
        integer, intent(in) :: grid_shape(:) !! Positive grid dimensions; last axis varies fastest.
        character(len=*), intent(in), optional :: boundary_conditions !! 'neumann' (default), 'dirichlet', or 'periodic'.
        integer, intent(out), optional :: status !! Zero on success, invalid-input otherwise.
        type(laplacian_nd_operator) :: op
        integer(int64) :: total
        integer :: d, j
        character(len=:), allocatable :: bc

        if (present(status)) status = laplacian_nd_invalid_input
        d = size(grid_shape)
        if (d == 0 .or. any(grid_shape < 1)) return
        bc = 'neumann'
        if (present(boundary_conditions)) bc = trim(boundary_conditions)
        select case (bc)
        case ('dirichlet')
            op%boundary = bc_dirichlet
        case ('neumann')
            op%boundary = bc_neumann
        case ('periodic')
            op%boundary = bc_periodic
        case default
            return
        end select
        total = 1_int64
        do j = 1, d
            if (total > int(huge(0)-1, int64)/int(grid_shape(j), int64)) return
            total = total*int(grid_shape(j), int64)
        end do
        op%n = int(total)
        op%grid_shape = grid_shape
        allocate(op%strides(d))
        op%strides(d) = 1
        do j = d-1, 1, -1
            op%strides(j) = op%strides(j+1)*grid_shape(j+1)
        end do
        if (present(status)) status = laplacian_nd_success
    end function laplacian_nd

    pure logical function valid_laplacian_nd(op) result(ok)
        type(laplacian_nd_operator), intent(in) :: op !! Candidate operator.
        integer :: j, d
        integer(int64) :: total

        ok = .false.
        if (.not. allocated(op%grid_shape) .or. .not. allocated(op%strides)) return
        d = size(op%grid_shape)
        if (d == 0 .or. size(op%strides) /= d) return
        if (op%boundary < bc_dirichlet .or. op%boundary > bc_periodic) return
        if (any(op%grid_shape < 1) .or. op%n < 1) return
        total = 1_int64
        do j = 1, d
            if (total > int(huge(0)-1, int64)/int(op%grid_shape(j), int64)) return
            total = total*int(op%grid_shape(j), int64)
        end do
        if (total /= int(op%n, int64)) return
        if (op%strides(d) /= 1) return
        do j = d-1, 1, -1
            if (int(op%strides(j),int64) /= &
                int(op%strides(j+1),int64)*int(op%grid_shape(j+1),int64)) return
        end do
        ok = .true.
    end function valid_laplacian_nd

    subroutine laplacian_nd_shape(op, shape, size_total, status)
        type(laplacian_nd_operator), intent(in) :: op !! Initialized grid Laplacian.
        integer, allocatable, intent(out) :: shape(:) !! Grid shape, in C order.
        integer, intent(out) :: size_total !! Number of grid points.
        integer, intent(out), optional :: status !! Zero on success.

        allocate(shape(0))
        size_total = 0
        if (present(status)) status = laplacian_nd_invalid_input
        if (.not. valid_laplacian_nd(op)) return
        shape = op%grid_shape
        size_total = op%n
        if (present(status)) status = laplacian_nd_success
    end subroutine laplacian_nd_shape

    subroutine laplacian_nd_matvec(op, x, y, status)
        type(laplacian_nd_operator), intent(in) :: op !! Initialized Laplacian.
        real(dp), intent(in) :: x(:) !! Finite vector of length product(grid_shape).
        real(dp), allocatable, intent(out) :: y(:) !! Matrix-free Laplacian times x.
        integer, intent(out), optional :: status !! Zero on success.
        integer :: i, j, n_axis, stride, coordinate

        allocate(y(0))
        if (present(status)) status = laplacian_nd_invalid_input
        if (.not. valid_laplacian_nd(op)) return
        if (size(x) /= op%n .or. .not. all(ieee_is_finite(x))) return
        deallocate(y)
        allocate(y(op%n))
        y = -2.0_dp*real(size(op%grid_shape),dp)*x
        do i = 1, op%n
            do j = 1, size(op%grid_shape)
                n_axis = op%grid_shape(j)
                stride = op%strides(j)
                coordinate = mod((i-1)/stride, n_axis)
                if (coordinate > 0) then
                    y(i) = y(i)+x(i-stride)
                else if (op%boundary == bc_periodic) then
                    y(i) = y(i)+x(i+(n_axis-1)*stride)
                else if (op%boundary == bc_neumann) then
                    y(i) = y(i)+x(i)
                end if
                if (coordinate < n_axis-1) then
                    y(i) = y(i)+x(i+stride)
                else if (op%boundary == bc_periodic) then
                    y(i) = y(i)+x(i-(n_axis-1)*stride)
                else if (op%boundary == bc_neumann) then
                    y(i) = y(i)+x(i)
                end if
            end do
        end do
        if (present(status)) status = laplacian_nd_success
    end subroutine laplacian_nd_matvec

    subroutine laplacian_nd_matmat(op, x, y, status)
        type(laplacian_nd_operator), intent(in) :: op !! Initialized Laplacian.
        real(dp), intent(in) :: x(:,:) !! Finite matrix with product(grid_shape) rows.
        real(dp), allocatable, intent(out) :: y(:,:) !! Operator applied independently to each column.
        integer, intent(out), optional :: status !! Zero on success.
        real(dp), allocatable :: column(:)
        integer :: k

        allocate(y(0,0))
        if (present(status)) status = laplacian_nd_invalid_input
        if (.not. valid_laplacian_nd(op)) return
        if (size(x,1) /= op%n .or. .not. all(ieee_is_finite(x))) return
        deallocate(y)
        allocate(y(op%n,size(x,2)))
        do k = 1, size(x,2)
            call laplacian_nd_matvec(op, x(:,k), column)
            y(:,k) = column
        end do
        if (present(status)) status = laplacian_nd_success
    end subroutine laplacian_nd_matmat

    function laplacian_nd_tocsr(op, status) result(a)
        type(laplacian_nd_operator), intent(in) :: op !! Initialized grid Laplacian.
        integer, intent(out), optional :: status !! Zero on success.
        type(csr_matrix) :: a
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: values(:)
        integer :: i, j, n_axis, stride, coordinate, n_entry, istat, d
        integer(int64) :: max_entry
        real(dp) :: diag

        if (present(status)) status = laplacian_nd_invalid_input
        if (.not. valid_laplacian_nd(op)) return
        d = size(op%grid_shape)
        if (int(d,int64)*2_int64+1_int64 > int(huge(0),int64)) return
        max_entry = int(op%n,int64)*(2_int64*int(d,int64)+1_int64)
        if (max_entry > int(huge(0)-1,int64)) return
        allocate(rows(int(max_entry)), cols(int(max_entry)), values(int(max_entry)))
        n_entry = 0
        do i = 1, op%n
            diag = -2.0_dp*real(d,dp)
            do j = 1, d
                n_axis = op%grid_shape(j)
                stride = op%strides(j)
                coordinate = mod((i-1)/stride, n_axis)
                if (coordinate > 0) then
                    call add(i,i-stride,1.0_dp)
                else if (op%boundary == bc_periodic) then
                    call add(i,i+(n_axis-1)*stride,1.0_dp)
                else if (op%boundary == bc_neumann) then
                    diag = diag+1.0_dp
                end if
                if (coordinate < n_axis-1) then
                    call add(i,i+stride,1.0_dp)
                else if (op%boundary == bc_periodic) then
                    call add(i,i-(n_axis-1)*stride,1.0_dp)
                else if (op%boundary == bc_neumann) then
                    diag = diag+1.0_dp
                end if
            end do
            call add(i,i,diag)
        end do
        a = csr_from_coo(op%n,op%n,rows(:n_entry),cols(:n_entry),values(:n_entry),istat)
        if (present(status)) then
            if (istat == sparse_success) status = laplacian_nd_success
        end if
    contains
        subroutine add(row,col,value)
            integer, intent(in) :: row !! One-based row.
            integer, intent(in) :: col !! One-based column.
            real(dp), intent(in) :: value !! Entry to accumulate.
            n_entry = n_entry+1
            rows(n_entry) = row
            cols(n_entry) = col
            values(n_entry) = value
        end subroutine add
    end function laplacian_nd_tocsr

    subroutine laplacian_nd_eigenvalues(op, eigenvalues, status, m)
        type(laplacian_nd_operator), intent(in) :: op !! Initialized grid Laplacian.
        real(dp), allocatable, intent(out) :: eigenvalues(:) !! Ascending; m selects m largest eigenvalues, like SciPy.
        integer, intent(out), optional :: status !! Zero on success.
        integer, intent(in), optional :: m !! Number of largest eigenvalues; default all.
        real(dp), allocatable :: spectrum(:)
        integer, allocatable :: indices(:)
        integer :: nret

        allocate(eigenvalues(0))
        if (present(status)) status = laplacian_nd_invalid_input
        if (.not. valid_laplacian_nd(op)) return
        nret = op%n
        if (present(m)) nret = m
        if (nret < 1 .or. nret > op%n) return
        call laplacian_spectrum(op,spectrum,indices)
        eigenvalues = spectrum(op%n-nret+1:)
        if (present(status)) status = laplacian_nd_success
    end subroutine laplacian_nd_eigenvalues

    subroutine laplacian_nd_eigenvectors(op, eigenvectors, status, m)
        type(laplacian_nd_operator), intent(in) :: op !! Initialized grid Laplacian.
        real(dp), allocatable, intent(out) :: eigenvectors(:,:) !! Orthonormal columns for ascending largest eigenvalues.
        integer, intent(out), optional :: status !! Zero on success.
        integer, intent(in), optional :: m !! Number of largest eigenpairs; default all.
        real(dp), allocatable :: spectrum(:)
        integer, allocatable :: indices(:)
        integer :: nret, i, k, j, mode, coordinate
        real(dp) :: phi

        allocate(eigenvectors(0,0))
        if (present(status)) status = laplacian_nd_invalid_input
        if (.not. valid_laplacian_nd(op)) return
        nret = op%n
        if (present(m)) nret = m
        if (nret < 1 .or. nret > op%n) return
        call laplacian_spectrum(op,spectrum,indices)
        deallocate(eigenvectors)
        allocate(eigenvectors(op%n,nret))
        do k = 1, nret
            do i = 1, op%n
                phi = 1.0_dp
                do j = 1, size(op%grid_shape)
                    mode = mod((indices(op%n-nret+k)-1)/op%strides(j),op%grid_shape(j))
                    coordinate = mod((i-1)/op%strides(j),op%grid_shape(j))
                    phi = phi*laplacian_1d_eigenvector(op%boundary,op%grid_shape(j),mode,coordinate)
                end do
                if (abs(phi) < epsilon(1.0_dp)) phi = 0.0_dp
                eigenvectors(i,k) = phi
            end do
        end do
        if (present(status)) status = laplacian_nd_success
    end subroutine laplacian_nd_eigenvectors

    pure real(dp) function laplacian_1d_eigenvector(bc,n,mode,coordinate) result(phi)
        integer, intent(in) :: bc !! Boundary-condition identifier.
        integer, intent(in) :: n !! Dimension length.
        integer, intent(in) :: mode !! Zero-based spectral mode.
        integer, intent(in) :: coordinate !! Zero-based spatial coordinate.
        real(dp) :: pi

        pi = acos(-1.0_dp)
        select case (bc)
        case (bc_dirichlet)
            phi = sqrt(2.0_dp/real(n+1,dp))* &
                sin(pi*real(coordinate+1,dp)*real(mode+1,dp)/real(n+1,dp))
        case (bc_neumann)
            phi = sqrt(merge(1.0_dp,2.0_dp,mode == 0)/real(n,dp))* &
                cos(pi*(real(coordinate,dp)+0.5_dp)*real(mode,dp)/real(n,dp))
        case default
            if (mode == 0) then
                phi = sqrt(1.0_dp/real(n,dp))
            else if (mod(n,2) == 0 .and. mode == n-1) then
                phi = sqrt(1.0_dp/real(n,dp))*merge(1.0_dp,-1.0_dp,mod(coordinate,2) == 0)
            else if (mod(mode+1,2) == 0) then
                phi = sqrt(2.0_dp/real(n,dp))* &
                    sin(pi*(real(coordinate,dp)+0.5_dp)*real(mode+1,dp)/real(n,dp))
            else
                phi = sqrt(2.0_dp/real(n,dp))* &
                    cos(pi*(real(coordinate,dp)+0.5_dp)*real(mode,dp)/real(n,dp))
            end if
        end select
    end function laplacian_1d_eigenvector

    subroutine laplacian_spectrum(op,values,indices)
        type(laplacian_nd_operator), intent(in) :: op !! Valid grid Laplacian.
        real(dp), allocatable, intent(out) :: values(:) !! Sorted eigenvalues.
        integer, allocatable, intent(out) :: indices(:) !! Original C-order spectral indices.
        real(dp), allocatable :: scratch(:)
        integer, allocatable :: idx_scratch(:)
        real(dp) :: pi, angle
        integer :: i, j, mode, n_axis, n, width, first, center, last, left, right, out

        n = op%n
        pi = acos(-1.0_dp)
        allocate(values(n), indices(n), scratch(n), idx_scratch(n))
        do i = 1, n
            values(i) = 0.0_dp
            indices(i) = i
            do j = 1, size(op%grid_shape)
                n_axis = op%grid_shape(j)
                mode = mod((i-1)/op%strides(j),n_axis)
                select case (op%boundary)
                case (bc_dirichlet)
                    angle = pi*real(mode+1,dp)/(2.0_dp*real(n_axis+1,dp))
                case (bc_neumann)
                    angle = pi*real(mode,dp)/(2.0_dp*real(n_axis,dp))
                case default
                    angle = pi*real((mode+1)/2,dp)/real(n_axis,dp)
                end select
                values(i) = values(i)-4.0_dp*sin(angle)**2
            end do
        end do
        ! Stable bottom-up merge sort. It also retains the spectral mode for
        ! each eigenvalue so that eigenvectors can be generated analytically.
        width = 1
        do while (width < n)
            first = 1
            do while (first <= n)
                center = first+min(width,n-first+1)-1
                last = center+min(width,n-center)
                left = first
                right = center+1
                do out = first, last
                    if (left > center) then
                        scratch(out) = values(right)
                        idx_scratch(out) = indices(right)
                        right = right+1
                    else if (right > last) then
                        scratch(out) = values(left)
                        idx_scratch(out) = indices(left)
                        left = left+1
                    else if (values(left) <= values(right)) then
                        scratch(out) = values(left)
                        idx_scratch(out) = indices(left)
                        left = left+1
                    else
                        scratch(out) = values(right)
                        idx_scratch(out) = indices(right)
                        right = right+1
                    end if
                end do
                if (last == n) exit
                first = last+1
            end do
            values = scratch
            indices = idx_scratch
            if (width > n/2) exit
            width = 2*width
        end do
    end subroutine laplacian_spectrum

end module scifort_sparse_laplacian_nd
