! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Finite proper-rotation symmetry groups and symmetry-aware reduction in SO(3).
! Quaternion convention: (x,y,z,w). Group ordering is intentionally not
! guaranteed to be identical to SciPy's ordering; indices are one-based.
module scifort_spatial_rotation_group
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, &
        rotation_from_quat, rotation_from_rotvec, spatial_rotation_success, &
        spatial_rotation_invalid_input
    use scifort_spatial_rotation_stack, only : rotation_stack_type, rotation_stack_from_rotations
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: max_cyclic_order = 10000

    public :: rotation_create_group, rotation_reduce

    interface rotation_reduce
        module procedure reduce_rotation
        module procedure reduce_rotation_stack
    end interface rotation_reduce

contains

    pure subroutine mark(status, code)
        integer, intent(out), optional :: status !! zero for success, one for invalid input
        integer, intent(in) :: code !! spatial rotation status code
        if (present(status)) status = code
    end subroutine mark

    function invalid_rotation() result(r)
        type(rotation_type) :: r
        r = rotation_from_quat([0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp])
    end function invalid_rotation

    ! Add a quaternion unless it is equal up to the q/-q equivalence. The
    ! tetrahedral, octahedral and icosahedral groups contain at most 60 elements.
    subroutine add_unique(quats, n, q)
        real(dp), intent(inout) :: quats(4,60) !! accepted scalar-last unit quaternions
        integer, intent(inout) :: n !! number of accepted rotations, updated in place
        real(dp), intent(in) :: q(4) !! new quaternion in scalar-last convention
        real(dp) :: v(4), mag
        integer :: i
        mag = sqrt(sum(q*q))
        if (mag <= 0.0_dp) return
        v = q/mag
        do i = 1, n
            if (abs(dot_product(v,quats(:,i))) >= 1.0_dp-1.0e-13_dp) return
        end do
        if (n >= size(quats,2)) return
        n = n+1
        quats(:,n) = v
    end subroutine add_unique

    ! Unit quaternions for the binary tetrahedral group (24 elements, 12
    ! proper rotations), extended to the binary octahedral and binary
    ! icosahedral groups. The latter uses the 120 vertices of the 600-cell:
    ! coordinate permutations are EVEN and the three nonzero signs independent.
    subroutine polyhedral_quaternions(group, quats, n)
        character(len=1), intent(in) :: group !! one of T, O, I
        real(dp), intent(out) :: quats(4,60) !! distinct elements, scalar-last
        integer, intent(out) :: n !! group order
        real(dp) :: q(4), base(4), root_half, phi
        integer :: i, j, k, mask, axis, other_axis, sign, perms(4), inversion, idx
        integer :: p1, p2, p3, p4
        quats = 0.0_dp
        n = 0
        ! Coordinate unit quaternions give the identity and three half-turns.
        do i = 1, 4
            q = 0.0_dp
            q(i) = 1.0_dp
            call add_unique(quats, n, q)
        end do
        ! Eight distinct one-third-turn rotations and their inverses.
        do mask = 0, 15
            q = 0.5_dp
            do i = 1, 4
                if (btest(mask,i-1)) q(i) = -q(i)
            end do
            call add_unique(quats, n, q)
        end do
        if (group == 'T') return
        if (group == 'O') then
            root_half = sqrt(0.5_dp)
            ! Quarter-turns about the three coordinate axes.
            do axis = 1, 3
                do sign = -1, 1, 2
                    q = 0.0_dp
                    q(axis) = real(sign,dp)*root_half
                    q(4) = root_half
                    call add_unique(quats, n, q)
                end do
            end do
            ! Half-turns about the six edge-midpoint axes of a cube.
            do axis = 1, 2
                do other_axis = axis+1, 3
                    do sign = -1, 1, 2
                        q = 0.0_dp
                        q(axis) = root_half
                        q(other_axis) = real(sign,dp)*root_half
                        call add_unique(quats, n, q)
                    end do
                end do
            end do
            return
        end if
        phi = (1.0_dp+sqrt(5.0_dp))/2.0_dp
        base = [0.0_dp, 0.5_dp, 0.5_dp*phi, 0.5_dp/phi]
        do p1 = 1, 4
            do p2 = 1, 4
                if (p2 == p1) cycle
                do p3 = 1, 4
                    if (p3 == p1 .or. p3 == p2) cycle
                    do p4 = 1, 4
                        if (p4 == p1 .or. p4 == p2 .or. p4 == p3) cycle
                        perms = [p1, p2, p3, p4]
                        inversion = 0
                        do i = 1, 3
                            do j = i+1, 4
                                if (perms(i) > perms(j)) inversion = inversion+1
                            end do
                        end do
                        if (mod(inversion,2) /= 0) cycle
                        do mask = 0, 7
                            k = 0
                            q = base(perms)
                            do idx = 1, 4
                                if (q(idx) == 0.0_dp) cycle
                                if (btest(mask,k)) q(idx) = -q(idx)
                                k = k+1
                            end do
                            call add_unique(quats, n, q)
                        end do
                    end do
                end do
            end do
        end do
    end subroutine polyhedral_quaternions

    ! Match SciPy Rotation.create_group('Cn'|'Dn'|'T'|'O'|'I'). Group
    ! elements may have a different enumeration order than in SciPy.
    function rotation_create_group(group, axis, status) result(out)
        character(len=*), intent(in) :: group !! 'Cn', 'Dn', 'T', 'O', or 'I', case-sensitive
        character(len=*), intent(in), optional :: axis !! X/Y/Z for Cn/Dn (default Z); ignored for T/O/I
        integer, intent(out), optional :: status !! zero on success, invalid input otherwise
        type(rotation_stack_type) :: out
        type(rotation_type), allocatable :: rotations(:)
        type(rotation_type) :: perpendicular
        real(dp) :: quats(4,60), v(3), theta
        integer :: n, total, i, ios, code, axis_id, cross_id
        character(len=1) :: group_type, axis_char
        call mark(status, spatial_rotation_invalid_input)
        if (group == 'T' .or. group == 'O' .or. group == 'I') then
            call polyhedral_quaternions(group, quats, n)
            if ((group == 'T' .and. n /= 12) .or. &
                (group == 'O' .and. n /= 24) .or. &
                (group == 'I' .and. n /= 60)) return
            allocate(rotations(n))
            do i = 1, n
                rotations(i) = rotation_from_quat(quats(:,i), status=code)
                if (code /= spatial_rotation_success) return
            end do
        else
            if (len_trim(group) < 2) return
            group_type = group(1:1)
            if (group_type /= 'C' .and. group_type /= 'D') return
            do i = 2, len_trim(group)
                if (group(i:i) < '0' .or. group(i:i) > '9') return
            end do
            read(group(2:len_trim(group)),*,iostat=ios) n
            if (ios /= 0) return
            if (n < 1 .or. n > max_cyclic_order) return
            axis_char = 'Z'
            if (present(axis)) then
                if (len_trim(axis) /= 1) return
                axis_char = axis(1:1)
            end if
            select case (axis_char)
            case ('X','x')
                axis_id = 1
                cross_id = 2
            case ('Y','y')
                axis_id = 2
                cross_id = 3
            case ('Z','z')
                axis_id = 3
                cross_id = 1
            case default
                return
            end select
            total = n
            if (group_type == 'D') total = 2*n
            allocate(rotations(total))
            do i = 1, n
                v = 0.0_dp
                theta = 2.0_dp*pi*real(i-1,dp)/real(n,dp)
                v(axis_id) = theta
                rotations(i) = rotation_from_rotvec(v, status=code)
                if (code /= spatial_rotation_success) return
            end do
            if (group_type == 'D') then
                v = 0.0_dp
                v(cross_id) = pi
                perpendicular = rotation_from_rotvec(v, status=code)
                if (code /= spatial_rotation_success) return
                do i = 1, n
                    rotations(n+i) = rotations(i)%compose(perpendicular)
                end do
            end if
        end if
        out = rotation_stack_from_rotations(rotations, code)
        if (code /= spatial_rotation_success) return
        call mark(status, spatial_rotation_success)
    end function rotation_create_group

    ! Return the minimum-angle representative of l*p*r over the supplied
    ! symmetries. Abs(q_w) is monotonic in the principal rotation angle;
    ! maximizing it avoids repeated acos and is stable near zero rotation.
    function reduce_rotation(p, left, right, index_left, index_right, status) result(out)
        type(rotation_type), intent(in) :: p !! input rotation
        type(rotation_stack_type), intent(in), optional :: left !! left symmetry group
        type(rotation_stack_type), intent(in), optional :: right !! right symmetry group
        integer, intent(out), optional :: index_left !! 1-based group index, or 0 if omitted
        integer, intent(out), optional :: index_right !! 1-based group index, or 0 if omitted
        integer, intent(out), optional :: status !! zero on success, invalid input otherwise
        type(rotation_type) :: out, lr, trial, symmetry
        integer :: nl, nr, i, j, best_i, best_j
        real(dp) :: score, best, quat(4)
        out = invalid_rotation()
        call mark(status, spatial_rotation_invalid_input)
        if (present(index_left)) index_left = 0
        if (present(index_right)) index_right = 0
        if (.not. p%valid()) return
        nl = 1
        nr = 1
        if (present(left)) then
            if (.not. left%valid()) return
            nl = left%size()
            if (nl == 0) return
        end if
        if (present(right)) then
            if (.not. right%valid()) return
            nr = right%size()
            if (nr == 0) return
        end if
        best = -1.0_dp
        best_i = 0
        best_j = 0
        do i = 1, nl
            lr = p
            if (present(left)) then
                symmetry = left%at(i)
                lr = symmetry%compose(p)
            end if
            do j = 1, nr
                trial = lr
                if (present(right)) then
                    symmetry = right%at(j)
                    trial = lr%compose(symmetry)
                end if
                quat = trial%as_quat()
                score = abs(quat(4))
                if (score > best) then
                    best = score
                    out = trial
                    best_i = i
                    best_j = j
                end if
            end do
        end do
        if (present(index_left)) then
            if (present(left)) index_left = best_i
        end if
        if (present(index_right)) then
            if (present(right)) index_right = best_j
        end if
        call mark(status, spatial_rotation_success)
    end function reduce_rotation

    ! Batch version corresponding to Rotation.reduce on a Rotation stack.
    ! Each input rotation is reduced independently using the same groups.
    function reduce_rotation_stack(p, left, right, index_left, index_right, status) result(out)
        type(rotation_stack_type), intent(in) :: p !! input stack
        type(rotation_stack_type), intent(in), optional :: left !! left symmetry group
        type(rotation_stack_type), intent(in), optional :: right !! right symmetry group
        integer, allocatable, intent(out), optional :: index_left(:) !! 1-based group indices; zeros if omitted
        integer, allocatable, intent(out), optional :: index_right(:) !! 1-based group indices; zeros if omitted
        integer, intent(out), optional :: status !! zero on success, invalid input otherwise
        type(rotation_stack_type) :: out
        type(rotation_type), allocatable :: result(:)
        integer :: i, il, ir, code
        call mark(status, spatial_rotation_invalid_input)
        if (.not. p%valid() .or. p%size() < 1) return
        if (present(left)) then
            if (.not. left%valid() .or. left%size() == 0) return
        end if
        if (present(right)) then
            if (.not. right%valid() .or. right%size() == 0) return
        end if
        allocate(result(p%size()))
        if (present(index_left)) allocate(index_left(p%size()))
        if (present(index_right)) allocate(index_right(p%size()))
        do i = 1, p%size()
            result(i) = reduce_rotation(p%at(i), left, right, il, ir, code)
            if (code /= spatial_rotation_success) return
            if (present(index_left)) index_left(i) = il
            if (present(index_right)) index_right(i) = ir
        end do
        out = rotation_stack_from_rotations(result, code)
        if (code /= spatial_rotation_success) return
        call mark(status, spatial_rotation_success)
    end function reduce_rotation_stack

end module scifort_spatial_rotation_group
