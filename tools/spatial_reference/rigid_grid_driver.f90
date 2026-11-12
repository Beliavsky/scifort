! SPDX-License-Identifier: MIT
! Independent text-I/O driver for comparing SE(3) grid operations with SciPy.
program rigid_grid_driver
    use scifort_kinds, only : dp
    use scifort_spatial_rigid_transform, only : rigid_transform_type, rigid_transform_stack_type
    use scifort_spatial_rigid_transform_grid, only : rigid_transform_grid_type, &
        rigid_grid_from_exp_coords, rigid_grid_from_dual_quat, rigid_grid_from_matrix
    implicit none
    type(rigid_transform_grid_type) :: a, b, composed, inverse, halfway, reconstructed
    type(rigid_transform_stack_type) :: m1, m2
    type(rigid_transform_type) :: mean_tf
    real(dp), allocatable :: x(:,:,:), y(:,:,:), v(:,:,:), w(:,:), mats(:,:,:,:), xyz(:,:,:)
    integer :: n,m,p,q,c,d,i,j,status

    read(*,*) n,m,p,q,c,d
    if (min(n,m,p,q,c,d)<1) error stop 'nonpositive dimensions'
    allocate(x(n,m,6),y(p,q,6),v(c,d,3),w(n,m))
    do j=1,m
        do i=1,n
            read(*,*) x(i,j,:)
        end do
    end do
    do j=1,q
        do i=1,p
            read(*,*) y(i,j,:)
        end do
    end do
    do j=1,d
        do i=1,c
            read(*,*) v(i,j,:)
        end do
    end do
    do j=1,m
        do i=1,n
            read(*,*) w(i,j)
        end do
    end do

    a = rigid_grid_from_exp_coords(x,status)
    if (status/=0 .or. .not. a%valid()) error stop 'construct a'
    b = rigid_grid_from_exp_coords(y,status)
    if (status/=0 .or. .not. b%valid()) error stop 'construct b'
    composed = a%compose(b)
    if (.not. composed%valid()) error stop 'compose'
    inverse = a%inverse()
    if (.not. inverse%valid()) error stop 'inverse'
    halfway = a%power(0.5_dp,status)
    if (status/=0) error stop 'power'
    xyz = a%apply(v)
    if (any(shape(xyz) < 1)) error stop 'apply'

    mats = a%as_matrix()
    call print_matrices(mats)
    call print_matrices(composed%as_matrix())
    call print_vectors(xyz)
    call print_vectors(a%apply(v,inverse=.true.))
    call print_matrices(inverse%as_matrix())
    reconstructed = rigid_grid_from_dual_quat(a%as_dual_quat(scalar_first=.true.), &
        scalar_first=.true., status=status)
    if (status/=0) error stop 'dual roundtrip'
    call print_matrices(reconstructed%as_matrix())
    reconstructed = rigid_grid_from_matrix(mats,status)
    if (status/=0) error stop 'matrix roundtrip'
    call print_matrices(reconstructed%as_matrix())
    reconstructed = rigid_grid_from_exp_coords(a%as_exp_coords(),status)
    if (status/=0) error stop 'exp roundtrip'
    call print_matrices(reconstructed%as_matrix())
    mean_tf = a%mean(w,status)
    if (status/=0 .or. .not. mean_tf%valid()) error stop 'mean'
    call print_matrix(mean_tf%as_matrix())
    m1 = a%mean_axis(1,w,status)
    if (status/=0 .or. .not. m1%valid()) error stop 'mean axis 1'
    call print_stack(m1)
    m2 = a%mean_axis(2,w,status)
    if (status/=0 .or. .not. m2%valid()) error stop 'mean axis 2'
    call print_stack(m2)
    call print_matrices(halfway%as_matrix())
contains
    subroutine print_matrix(mat)
        real(dp), intent(in) :: mat(4,4) !! homogeneous transform
        integer :: k
        do k=1,4
            write(*,'(4(es25.16e3,1x))') mat(k,:)
        end do
    end subroutine print_matrix
    subroutine print_matrices(arr)
        real(dp), intent(in) :: arr(:,:,:,:) !! matrices (n,m,4,4)
        integer :: i,j
        do j=1,size(arr,2)
            do i=1,size(arr,1)
                call print_matrix(arr(i,j,:,:))
            end do
        end do
    end subroutine print_matrices
    subroutine print_vectors(arr)
        real(dp), intent(in) :: arr(:,:,:) !! vectors (n,m,3)
        integer :: i,j
        do j=1,size(arr,2)
            do i=1,size(arr,1)
                write(*,'(3(es25.16e3,1x))') arr(i,j,:)
            end do
        end do
    end subroutine print_vectors
    subroutine print_stack(stack)
        type(rigid_transform_stack_type), intent(in) :: stack !! stack of rigid transforms
        type(rigid_transform_type) :: tf
        integer :: i
        do i=1,stack%size()
            tf=stack%at(i)
            call print_matrix(tf%as_matrix())
        end do
    end subroutine print_stack
end program rigid_grid_driver
