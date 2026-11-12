! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Greedy approximate multiplicity grouping for real and complex roots.
! Each group compares distances to its first unused member, not to its mean.
module scifort_signal_unique_roots
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: unique_roots_status_success = 0
    integer, parameter, public :: unique_roots_status_invalid_input = 1
    interface unique_roots
        module procedure unique_roots_real
        module procedure unique_roots_complex
    end interface
    public :: unique_roots

contains

    pure logical function less_complex(x, y) result(less)
        complex(dp), intent(in) :: x !! first value
        complex(dp), intent(in) :: y !! second value
        less = real(x,dp) < real(y,dp)
        if (real(x,dp) == real(y,dp)) less = aimag(x) < aimag(y)
    end function less_complex

    subroutine unique_roots_complex(p, roots, multiplicities, tol, rtype, status)
        complex(dp), intent(in) :: p(:) !! original roots in greedy iteration order
        complex(dp), allocatable, intent(out) :: roots(:) !! representative of each cluster
        integer, allocatable, intent(out) :: multiplicities(:) !! number of roots per cluster
        real(dp), optional, intent(in) :: tol !! absolute Euclidean distance threshold
        character(*), optional, intent(in) :: rtype !! 'min', 'max' or 'avg'
        integer, optional, intent(out) :: status !! success or invalid input
        complex(dp), allocatable :: representatives(:)
        integer, allocatable :: counts(:)
        logical, allocatable :: used(:)
        real(dp) :: distance
        character(:), allocatable :: how
        complex(dp) :: representative
        integer :: n, i, j, k, count

        allocate(roots(0), multiplicities(0))
        if (present(status)) status = unique_roots_status_invalid_input
        distance = 0.001_dp
        if (present(tol)) distance = tol
        how = 'min'
        if (present(rtype)) how = trim(rtype)
        if (.not. ieee_is_finite(distance) .or. distance < 0.0_dp) return
        if (.not. all(ieee_is_finite(real(p,dp))) .or. &
            .not. all(ieee_is_finite(aimag(p)))) return
        select case (how)
        case ('min', 'minimum', 'max', 'maximum', 'avg', 'mean')
        case default
            return
        end select
        n = size(p)
        allocate(representatives(n), counts(n), used(n))
        used = .false.
        k = 0
        do i = 1, n
            if (used(i)) cycle
            representative = p(i)
            count = 0
            do j = i,n
                if (used(j)) cycle
                if (abs(p(j)-p(i)) > distance) cycle
                used(j) = .true.
                count = count+1
                select case (how)
                case ('min', 'minimum')
                    if (less_complex(p(j),representative)) representative = p(j)
                case ('max', 'maximum')
                    if (less_complex(representative,p(j))) representative = p(j)
                case ('avg', 'mean')
                    if (count == 1) then
                        representative = p(j)
                    else
                        ! Online average prevents overflowing a sum of large roots.
                        representative = representative + (p(j)-representative)/real(count,dp)
                    end if
                end select
            end do
            k = k+1
            representatives(k) = representative
            counts(k) = count
        end do
        roots = representatives(:k)
        multiplicities = counts(:k)
        if (present(status)) status = unique_roots_status_success
    end subroutine unique_roots_complex

    subroutine unique_roots_real(p, roots, multiplicities, tol, rtype, status)
        real(dp), intent(in) :: p(:) !! original real roots
        real(dp), allocatable, intent(out) :: roots(:) !! representative real roots
        integer, allocatable, intent(out) :: multiplicities(:) !! multiplicity per group
        real(dp), optional, intent(in) :: tol !! distance threshold
        character(*), optional, intent(in) :: rtype !! min/max/avg or synonyms
        integer, optional, intent(out) :: status !! zero on success
        complex(dp), allocatable :: representatives(:)

        call unique_roots_complex(cmplx(p,0.0_dp,dp),representatives,multiplicities,tol,rtype,status)
        roots = real(representatives,dp)
    end subroutine unique_roots_real
end module scifort_signal_unique_roots
