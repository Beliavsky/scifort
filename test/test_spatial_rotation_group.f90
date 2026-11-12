! SPDX-License-Identifier: MIT
program test_spatial_rotation_group
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_identity, &
        rotation_from_rotvec, rotation_from_quat
    use scifort_spatial_rotation_stack, only : rotation_stack_type, rotation_stack_from_quat
    use scifort_spatial_rotation_group, only : rotation_create_group, rotation_reduce
    implicit none
    type(rotation_stack_type) :: group, left, right, stack, reduced, bad_group, empty
    type(rotation_type) :: p, q, ri, rj, prod, best, test, broken
    real(dp) :: angle
    real(dp), allocatable :: quats(:,:)
    integer :: status, i, j, k, il, ir, jl, jr
    integer, allocatable :: left_idx(:), right_idx(:)
    character(len=6), parameter :: names(9) = [character(len=6) :: &
        'C1', 'C2', 'C7', 'D1', 'D2', 'D7', 'T', 'O', 'I']
    integer, parameter :: sizes(9) = [1,2,7,2,4,14,12,24,60]

    do k=1,size(names)
        group = rotation_create_group(trim(names(k)),status=status)
        if (status /= 0 .or. .not. group%valid()) error stop 'group creation'
        if (group%size() /= sizes(k)) error stop 'group size'
        quats = group%as_quat()
        do i=1,group%size()
            ri = group%at(i)
            if (abs(sqrt(sum(quats(i,:)**2))-1.0_dp) > 2.0e-14_dp) error stop 'quaternion norm'
            prod = ri%inverse()
            call require_member(group,prod,4.0e-12_dp)
            do j=i+1,group%size()
                if (abs(dot_product(quats(i,:),quats(j,:))) > 1.0_dp-1.0e-12_dp) &
                    error stop 'duplicate rotation'
            end do
            do j=1,group%size()
                rj = group%at(j)
                prod = ri%compose(rj)
                call require_member(group,prod,4.0e-12_dp)
            end do
        end do
    end do

    do k=1,3
        select case(k)
        case(1)
            group = rotation_create_group('D5',axis='X',status=status)
        case(2)
            group = rotation_create_group('D5',axis='y',status=status)
        case(3)
            group = rotation_create_group('D5',axis='Z',status=status)
        end select
        if (status /= 0 .or. group%size() /= 10) error stop 'dihedral axes'
        do i=1,group%size()
            ri=group%at(i)
            do j=1,group%size()
                rj=group%at(j)
                call require_member(group,ri%compose(rj),4.0e-12_dp)
            end do
        end do
    end do

    group = rotation_create_group('O',axis='invalid',status=status)
    if (status /= 0 .or. group%size() /= 24) error stop 'polyhedral ignores axis'
    group = rotation_create_group('C4',axis='X',status=status)
    if (status /= 0 .or. group%size() /= 4) error stop 'cyclic axis'
    ri=group%at(2)
    call assert_close(ri%apply([0.0_dp,1.0_dp,0.0_dp]),[0.0_dp,0.0_dp,1.0_dp],3.0e-14_dp)

    ! Group reduction must agree with brute-force enumeration and report the
    ! indices of the actual left/right symmetries that yield its result.
    left = rotation_create_group('O')
    right = rotation_create_group('D7')
    p = rotation_from_rotvec([0.3_dp,-1.2_dp,1.7_dp])
    best = rotation_reduce(p,left,right,il,ir,status)
    if (status /= 0 .or. il<1 .or. il>left%size() .or. ir<1 .or. ir>right%size()) &
        error stop 'reduce indices'
    ri=left%at(il)
    rj=right%at(ir)
    q=ri%compose(p)
    prod=q%compose(rj)
    call assert_rotation(best,prod,3.0e-14_dp)
    angle=best%magnitude()
    do i=1,left%size()
        ri=left%at(i)
        q=ri%compose(p)
        do j=1,right%size()
            rj=right%at(j)
            prod=q%compose(rj)
            if (angle > prod%magnitude()+1.0e-12_dp) error stop 'nonoptimal reduction'
        end do
    end do

    q=rotation_reduce(p,index_left=il,index_right=ir,status=status)
    if (status /= 0 .or. il/=0 .or. ir/=0) error stop 'no groups'
    call assert_rotation(p,q,1.0e-14_dp)
    q=rotation_reduce(p,left=left,index_left=il,index_right=ir,status=status)
    if (status /= 0 .or. il<1 .or. ir/=0) error stop 'left only'
    q=rotation_reduce(p,right=right,index_left=il,index_right=ir,status=status)
    if (status /= 0 .or. ir<1 .or. il/=0) error stop 'right only'

    stack = rotation_stack_from_quat(reshape([0.1_dp,0.3_dp,0.4_dp,1.0_dp, &
        0.6_dp,0.2_dp,-0.1_dp,1.0_dp, -0.3_dp,0.1_dp,0.7_dp,1.0_dp], &
        [3,4],order=[2,1]))
    if (.not. stack%valid()) error stop 'make stack'
    reduced = rotation_reduce(stack,left,right,left_idx,right_idx,status)
    if (status /= 0 .or. reduced%size() /= stack%size()) error stop 'batch reduction'
    if (size(left_idx)/=stack%size() .or. size(right_idx)/=stack%size()) error stop 'batch indices'
    do i=1,stack%size()
        p=stack%at(i)
        test=rotation_reduce(p,left,right,jl,jr)
        q=reduced%at(i)
        call assert_rotation(q,test,3.0e-14_dp)
        if (left_idx(i)/=jl .or. right_idx(i)/=jr) error stop 'batch witnesses'
    end do
    reduced=rotation_reduce(stack,status=status)
    if (status/=0 .or. .not. all(stack%approx_equal(reduced))) error stop 'stack default'

    bad_group = rotation_create_group('C0',status=status)
    if (status == 0 .or. bad_group%valid()) error stop 'C0 accepted'
    bad_group = rotation_create_group('D-2',status=status)
    if (status == 0 .or. bad_group%valid()) error stop 'negative group accepted'
    bad_group = rotation_create_group('D10001',status=status)
    if (status == 0 .or. bad_group%valid()) error stop 'excess order accepted'
    bad_group = rotation_create_group('c3',status=status)
    if (status == 0 .or. bad_group%valid()) error stop 'lowercase group accepted'
    bad_group = rotation_create_group('C3',axis='Q',status=status)
    if (status == 0 .or. bad_group%valid()) error stop 'invalid axis accepted'
    bad_group = rotation_create_group('Cbad',status=status)
    if (status == 0 .or. bad_group%valid()) error stop 'invalid order accepted'
    bad_group = rotation_create_group('TQ',status=status)
    if (status == 0 .or. bad_group%valid()) error stop 'invalid polyhedral group accepted'
    broken = rotation_from_quat([0.0_dp,0.0_dp,0.0_dp,0.0_dp])
    q=rotation_reduce(broken,left,right,status=status)
    if (status == 0 .or. q%valid()) error stop 'invalid rotation accepted'
    q=rotation_reduce(rotation_identity(),left=empty,status=status)
    if (status == 0 .or. q%valid()) error stop 'invalid symmetry accepted'
    reduced=rotation_reduce(stack,left=empty,status=status)
    if (status == 0 .or. reduced%valid()) error stop 'invalid stack symmetry accepted'
    print '(A)', 'PASS: finite rotation groups, closure, reduction, indices and invalid inputs'
contains
    subroutine require_member(g,r,tol)
        type(rotation_stack_type), intent(in) :: g !! group to check
        type(rotation_type), intent(in) :: r !! candidate rotation
        real(dp), intent(in) :: tol !! maximum quaternion distance
        type(rotation_type) :: test_rotation
        real(dp) :: qa(4), qb(4), error
        integer :: j
        qa=r%as_quat()
        do j=1,g%size()
            test_rotation=g%at(j)
            qb=test_rotation%as_quat()
            error=min(maxval(abs(qa-qb)),maxval(abs(qa+qb)))
            if (error < tol) return
        end do
        error stop 'group closure/inverse'
    end subroutine require_member

    subroutine assert_close(a,b,tol)
        real(dp), intent(in) :: a(:) !! computed vector
        real(dp), intent(in) :: b(:) !! expected vector
        real(dp), intent(in) :: tol !! maximum absolute component discrepancy
        if (size(a)/=size(b)) error stop 'vector dimensions'
        if (maxval(abs(a-b)) > tol) error stop 'vector mismatch'
    end subroutine assert_close

    subroutine assert_rotation(a,b,tol)
        type(rotation_type), intent(in) :: a !! computed rotation
        type(rotation_type), intent(in) :: b !! expected rotation
        real(dp), intent(in) :: tol !! maximum quaternion discrepancy
        real(dp) :: qa(4),qb(4)
        qa=a%as_quat()
        qb=b%as_quat()
        if (min(maxval(abs(qa-qb)),maxval(abs(qa+qb)))>tol) error stop 'rotation mismatch'
    end subroutine assert_rotation
end program test_spatial_rotation_group
