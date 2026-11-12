! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Real state-feedback pole assignment for controllable SISO and MIMO systems.
! The partially actuated MIMO path supports KNV0 rank-one and YT
! determinant-ascent updates, with an alternative deterministic sweep order
! and transfer-vector initialization to SciPy's reference implementation.
module scifort_signal_place_poles
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_place_poles_mimo, only : place_poles_mimo_nullspace
    implicit none
    private
    integer, parameter, public :: pole_status_success = 0
    integer, parameter, public :: pole_status_invalid_input = 1
    integer, parameter, public :: pole_status_uncontrollable = 2
    integer, parameter, public :: pole_status_unsupported = 3
    integer, parameter, public :: pole_status_numerical_failure = 4
    type, public :: pole_placement_result
        real(dp), allocatable :: gain_matrix(:,:) !! feedback K such that A-B*K has requested poles
        complex(dp), allocatable :: requested_poles(:) !! input eigenvalues, not computed poles
        integer :: status = pole_status_invalid_input !! success, invalid or unsupported
        integer :: nb_iter = 0 !! determinant-ascent sweeps (MIMO only)
        real(dp) :: rtol = 0.0_dp !! achieved relative determinant change
        real(dp) :: det_transfer = 0.0_dp !! abs(det(real transfer basis)); MIMO only
    end type pole_placement_result
    interface place_poles
        module procedure place_poles_real
        module procedure place_poles_complex
    end interface
    public :: place_poles
contains
    function place_poles_real(a, b, poles, method, rtol, maxiter) result(out)
        real(dp), intent(in) :: a(:,:) !! real state matrix, shape n by n
        real(dp), intent(in) :: b(:,:) !! real control matrix, shape n by m
        real(dp), intent(in) :: poles(:) !! desired real closed-loop poles
        character(len=*), optional, intent(in) :: method !! "YT" (default) or "KNV0"
        real(dp), optional, intent(in) :: rtol !! determinant convergence target (1e-3)
        integer, optional, intent(in) :: maxiter !! maximum number of sweeps (30)
        type(pole_placement_result) :: out
        out = place_poles_complex(a,b,cmplx(poles,0.0_dp,dp),method,rtol,maxiter)
    end function place_poles_real

    function place_poles_complex(a, b, poles, method, rtol, maxiter) result(out)
        real(dp), intent(in) :: a(:,:) !! real state matrix, shape n by n
        real(dp), intent(in) :: b(:,:) !! real control matrix, shape n by m
        complex(dp), intent(in) :: poles(:) !! desired real poles and conjugate pairs
        character(len=*), optional, intent(in) :: method !! "YT" (default) or "KNV0"
        real(dp), optional, intent(in) :: rtol !! convergence threshold
        integer, optional, intent(in) :: maxiter !! iteration limit
        type(pole_placement_result) :: out
        real(dp), allocatable :: ctrl(:,:), rhs(:,:), v(:,:), phi(:,:), id(:,:), &
            gram(:,:), target(:,:)
        complex(dp), allocatable :: coeff(:), next_coeff(:)
        logical, allocatable :: used(:)
        logical :: ok
        real(dp) :: tol
        integer :: n, m, i, j, mate, max_iterations
        real(dp) :: target_rtol
        character(len=4) :: algorithm

        algorithm = 'YT'
        if (present(method)) then
            select case (trim(adjustl(method)))
            case ('YT','yt','Yt','yT')
                algorithm = 'YT'
            case ('KNV0','knv0','Knv0')
                algorithm = 'KNV0'
            case default
                return
            end select
        end if
        target_rtol = 1.0e-3_dp
        if (present(rtol)) target_rtol = rtol
        max_iterations = 30
        if (present(maxiter)) max_iterations = maxiter
        if (.not. ieee_is_finite(target_rtol) .or. target_rtol > 1.0_dp .or. &
            max_iterations < 1) return
        if (algorithm == 'KNV0') then
            if (any(abs(aimag(poles)) > 128.0_dp*epsilon(1.0_dp))) return
        end if
        n = size(a,1)
        m = size(b,2)
        if (n < 1 .or. m < 1 .or. size(a,2) /= n) return
        if (size(b,1) /= n .or. size(poles) /= n) return
        if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(b))) return
        if (.not. all(ieee_is_finite(real(poles,dp))) .or. &
                .not. all(ieee_is_finite(aimag(poles)))) return
        tol = 256.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(poles)))
        allocate(used(n))
        used = .false.
        ! A real gain matrix requires complex roots in conjugate pairs.
        do i=1,n
            if (used(i)) cycle
            used(i) = .true.
            if (abs(aimag(poles(i))) <= tol) cycle
            mate = 0
            do j=i+1,n
                if (used(j)) cycle
                if (abs(poles(j)-conjg(poles(i))) <= tol) then
                    mate = j
                    exit
                end if
            end do
            if (mate == 0) return
            used(mate) = .true.
        end do
        out%requested_poles = poles
        allocate(id(n,n))
        id = 0.0_dp
        do i=1,n
            id(i,i) = 1.0_dp
        end do

        if (m == 1) then
            ! SciPy rejects pole multiplicity larger than rank(B)=1.
            do i=1,n
                do j=i+1,n
                    if (abs(poles(i)-poles(j)) <= tol) return
                end do
            end do
            allocate(ctrl(n,n))
            ctrl(:,1) = b(:,1)
            do i=2,n
                ctrl(:,i) = matmul(a,ctrl(:,i-1))
            end do
            allocate(rhs(n,1))
            rhs = 0.0_dp
            rhs(n,1) = 1.0_dp
            call solve_real(transpose(ctrl),rhs,v,ok)
            if (.not. ok) then
                out%status = pole_status_uncontrollable
                return
            end if
            allocate(coeff(1))
            coeff(1) = (1.0_dp,0.0_dp)
            do i=1,n
                allocate(next_coeff(i+1))
                next_coeff = (0.0_dp,0.0_dp)
                do j=1,i
                    next_coeff(j) = next_coeff(j) + coeff(j)
                    next_coeff(j+1) = next_coeff(j+1)-poles(i)*coeff(j)
                end do
                deallocate(coeff)
                call move_alloc(next_coeff,coeff)
            end do
            if (any(abs(aimag(coeff)) > 10.0_dp*tol)) then
                out%status = pole_status_numerical_failure
                return
            end if
            phi = id
            do i=2,n+1
                phi = matmul(phi,a) + real(coeff(i),dp)*id
            end do
            allocate(out%gain_matrix(1,n))
            out%gain_matrix(1,:) = matmul(v(:,1),phi)
        else if (m >= n) then
            ! Fully actuated B: solve B*K = A-T using a minimum-norm
            ! real right inverse. No iterative eigenvector optimization.
            allocate(target(n,n))
            target = 0.0_dp
            used = .false.
            do i=1,n
                if (used(i)) cycle
                used(i) = .true.
                target(i,i) = real(poles(i),dp)
                if (abs(aimag(poles(i))) <= tol) cycle
                mate = 0
                do j=i+1,n
                    if (used(j)) cycle
                    if (abs(poles(j)-conjg(poles(i))) <= tol) then
                        mate = j
                        exit
                    end if
                end do
                if (mate == 0) return
                used(mate) = .true.
                target(mate,mate) = real(poles(i),dp)
                target(i,mate) = -aimag(poles(i))
                target(mate,i) = aimag(poles(i))
            end do
            gram = matmul(b,transpose(b))
            rhs = a-target
            call solve_real(gram,rhs,v,ok)
            if (.not. ok) then
                ! Even when m>=n, B may be rank deficient while (A,B)
                ! is controllable. Try the general MIMO assignment path.
                call place_poles_mimo_nullspace(a,b,poles,out%gain_matrix,out%status, &
                    algorithm,target_rtol,max_iterations,out%nb_iter,out%rtol,out%det_transfer)
                return
            end if
            out%gain_matrix = matmul(transpose(b),v)
        else
            call place_poles_mimo_nullspace(a,b,poles,out%gain_matrix,out%status, &
                algorithm,target_rtol,max_iterations,out%nb_iter,out%rtol,out%det_transfer)
            return
        end if
        if (.not. all(ieee_is_finite(out%gain_matrix))) then
            deallocate(out%gain_matrix)
            out%status = pole_status_numerical_failure
            return
        end if
        out%status = pole_status_success
    end function place_poles_complex

    subroutine solve_real(a, rhs, x, ok)
        real(dp), intent(in) :: a(:,:) !! n by n coefficient matrix
        real(dp), intent(in) :: rhs(:,:) !! n by p right-hand sides
        real(dp), allocatable, intent(out) :: x(:,:) !! solution matrix
        logical, intent(out) :: ok !! true if numerically nonsingular
        real(dp), allocatable :: mat(:,:), work(:,:), swap(:)
        real(dp) :: t, tol
        integer :: n, j, k, pivot
        ok = .false.
        n = size(a,1)
        if (size(a,2) /= n .or. size(rhs,1) /= n) return
        allocate(mat(n,n),work(n,size(rhs,2)),swap(max(n,size(rhs,2))))
        mat = a
        work = rhs
        tol = 128.0_dp*epsilon(1.0_dp)*real(n,dp)*maxval(abs(a))
        if (tol <= 0.0_dp) return
        do k=1,n
            pivot = k-1+maxloc(abs(mat(k:n,k)),dim=1)
            if (abs(mat(pivot,k)) <= tol) return
            if (pivot /= k) then
                swap(:n) = mat(k,:)
                mat(k,:) = mat(pivot,:)
                mat(pivot,:) = swap(:n)
                j = size(rhs,2)
                swap(:j) = work(k,:)
                work(k,:) = work(pivot,:)
                work(pivot,:) = swap(:j)
            end if
            t = mat(k,k)
            mat(k,:) = mat(k,:)/t
            work(k,:) = work(k,:)/t
            do j=1,n
                if (j == k) cycle
                t = mat(j,k)
                mat(j,:) = mat(j,:) - t*mat(k,:)
                work(j,:) = work(j,:) - t*work(k,:)
            end do
        end do
        if (.not. all(ieee_is_finite(work))) return
        x = work
        ok = .true.
    end subroutine solve_real
end module scifort_signal_place_poles
