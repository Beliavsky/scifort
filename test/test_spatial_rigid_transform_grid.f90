! SPDX-License-Identifier: MIT
program test_spatial_rigid_transform_grid
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_identity, rotation_from_rotvec
    use scifort_spatial_rotation_grid, only : rotation_grid_type, rotation_grid_from_rotvec
    use scifort_spatial_rigid_transform, only : rigid_transform_type, rigid_transform_stack_type, rigid_identity, &
        rigid_from_translation, rigid_stack_from_exp_coords
    use scifort_spatial_rigid_transform_grid
    implicit none
    type(rigid_transform_grid_type) :: grid, other, combined, back
    type(rigid_transform_stack_type) :: row, col, mean_row, mean_col, singleton
    type(rotation_type) :: rotations(1,3)
    type(rotation_grid_type) :: rg
    type(rigid_transform_type) :: tf, avg
    real(dp) :: t(2,1,3), x(1,3,3), mat(4,4), expected(3)
    real(dp) :: weights(2,1), rv(1,3,3)
    real(dp), allocatable :: mats(:,:,:,:), e(:,:,:), q(:,:,:), a(:,:,:)
    integer :: i, j, status

    t(1,1,:) = [1.0_dp,0.0_dp,2.0_dp]
    t(2,1,:) = [3.0_dp,0.0_dp,4.0_dp]
    do j=1,3
        rv(1,j,:) = [0.0_dp,0.0_dp,0.2_dp*real(j-1,dp)]
        rotations(1,j) = rotation_from_rotvec(rv(1,j,:))
    end do
    grid = rigid_grid_from_components(t,rotations,status)
    if (status/=0 .or. .not. grid%valid()) error stop 'components broadcast'
    if (any(grid%shape()/=[2,3])) error stop 'grid shape'
    rg = rotation_grid_from_rotvec(rv)
    other = rigid_grid_from_components(t,rg,status)
    if (status/=0 .or. maxval(abs(other%as_matrix()-grid%as_matrix()))>1.0e-14_dp) &
        error stop 'rotation_grid components'
    do j=1,3
        do i=1,2
            tf = grid%at(i,j)
            if (.not. tf%valid()) error stop 'at valid'
            expected = t(i,1,:) + rotations(1,j)%apply([1.0_dp,0.0_dp,0.0_dp])
            call check3(tf%apply([1.0_dp,0.0_dp,0.0_dp]),expected,1.0e-13_dp)
        end do
    end do
    tf = grid%at(3,1)
    if (tf%valid()) error stop 'invalid index'
    mats = grid%as_matrix()
    back = rigid_grid_from_matrix(mats,status)
    if (status/=0) error stop 'matrix roundtrip'
    call check4(back%as_matrix(),mats,1.0e-13_dp)
    e = grid%as_exp_coords()
    back = rigid_grid_from_exp_coords(e,status)
    if (status/=0) error stop 'exp roundtrip'
    call check4(back%as_matrix(),mats,1.0e-12_dp)
    q = grid%as_dual_quat(scalar_first=.true.)
    back = rigid_grid_from_dual_quat(q,scalar_first=.true.,status=status)
    if (status/=0) error stop 'dual quaternion roundtrip'
    call check4(back%as_matrix(),mats,1.0e-13_dp)
    other = rigid_grid_from_rotation(rg,status)
    if (status/=0 .or. any(other%shape()/=[1,3])) error stop 'rotation constructor'
    back = rigid_grid_from_translation(t,status)
    if (status/=0 .or. any(back%shape()/=[2,1])) error stop 'translation constructor'
    tf = rigid_identity()
    mat = tf%as_matrix()
    back = rigid_grid_identity([1,1],status)
    if (status/=0 .or. .not. back%valid()) error stop 'identity'
    combined = grid%compose(back)
    call check4(combined%as_matrix(),mats,1.0e-13_dp)
    combined = back%compose(grid)
    call check4(combined%as_matrix(),mats,1.0e-13_dp)
    combined = grid%compose(rigid_from_translation([0.0_dp,0.0_dp,0.0_dp]))
    call check4(combined%as_matrix(),mats,1.0e-13_dp)
    combined = grid%compose(grid%inverse())
    do j=1,3
        do i=1,2
            tf = combined%at(i,j)
            call checkmat(tf%as_matrix(),mat,2.0e-13_dp)
        end do
    end do
    combined = grid%power(0.5_dp,status)
    if (status/=0) error stop 'power'
    back = combined%compose(combined)
    call check4(back%as_matrix(),mats,1.0e-12_dp)

    x(1,1,:) = [1.0_dp,2.0_dp,3.0_dp]
    x(1,2,:) = [-2.0_dp,1.0_dp,0.0_dp]
    x(1,3,:) = [0.5_dp,-0.5_dp,1.0_dp]
    a = grid%apply(x)
    if (any(shape(a)/=[2,3,3])) error stop 'apply grid shape'
    do j=1,3
        do i=1,2
            tf = grid%at(i,j)
            call check3(a(i,j,:),tf%apply(x(1,j,:)),1.0e-13_dp)
            call check3(tf%apply(a(i,j,:),inverse=.true.),x(1,j,:),1.0e-13_dp)
        end do
    end do
    a = grid%apply(x,inverse=.true.)
    do j=1,3
        do i=1,2
            tf = grid%at(i,j)
            call check3(a(i,j,:),tf%apply(x(1,j,:),inverse=.true.),1.0e-13_dp)
        end do
    end do
    a = grid%apply([2.0_dp,-1.0_dp,0.3_dp])
    if (any(shape(a)/=[2,3,3])) error stop 'apply vector shape'
    do j=1,3
        do i=1,2
            tf = grid%at(i,j)
            call check3(a(i,j,:),tf%apply([2.0_dp,-1.0_dp,0.3_dp]),1.0e-13_dp)
        end do
    end do
    a = grid%apply(x(1,:,:))
    if (any(shape(a)/=[2,3,3])) error stop 'apply rank-two points shape'

    row = grid%row(0)
    if (row%valid()) error stop 'row invalid index'
    col = grid%column(4)
    if (col%valid()) error stop 'column invalid index'
    row = grid%row(1)
    col = grid%column(2)
    if (.not. row%valid() .or. row%size()/=3) error stop 'row'
    if (.not. col%valid() .or. col%size()/=2) error stop 'column'
    tf = row%at(3)
    call checkmat(tf%as_matrix(),mats(1,3,:,:),1.0e-13_dp)
    tf = col%at(2)
    call checkmat(tf%as_matrix(),mats(2,2,:,:),1.0e-13_dp)
    singleton = rigid_stack_from_exp_coords(e(1,:,:))
    combined = grid%compose(singleton)
    if (any(combined%shape()/=[2,3])) error stop 'compose stack broadcast'
    do j=1,3
        do i=1,2
            tf = combined%at(i,j)
            call checkmat(tf%as_matrix(),matmul(mats(i,j,:,:),mats(1,j,:,:)),2.0e-13_dp)
        end do
    end do

    weights(:,1) = [1.0_dp,3.0_dp]
    avg = grid%mean(weights=weights,status=status)
    if (status/=0) error stop 'weighted global mean'
    call check3(avg%translation(),[2.5_dp,0.0_dp,3.5_dp],1.0e-13_dp)
    mean_row = grid%mean_axis(1,weights=weights,status=status)
    if (status/=0 .or. .not. mean_row%valid() .or. mean_row%size()/=3) error stop 'mean axis 1'
    mean_col = grid%mean_axis(2,status=status)
    if (status/=0 .or. .not. mean_col%valid() .or. mean_col%size()/=2) error stop 'mean axis 2'
    do j=1,3
        tf = mean_row%at(j)
        call check3(tf%translation(),[2.5_dp,0.0_dp,3.5_dp],1.0e-13_dp)
    end do
    tf = mean_col%at(1)
    call check3(tf%translation(),[1.0_dp,0.0_dp,2.0_dp],1.0e-13_dp)

    other = rigid_grid_from_matrix(reshape([real(dp)::],[0,0,4,4]),status)
    if (status==0 .or. other%valid()) error stop 'zero-sized invalid'
    other = rigid_grid_identity([0,2],status)
    if (status==0 .or. other%valid()) error stop 'identity zero size'
    other = rigid_grid_from_components(t,spread(rotations(1,1:2),dim=1,ncopies=3),status)
    if (status==0 .or. other%valid()) error stop 'bad constructor dimensions'
    other = rigid_grid_from_exp_coords(e(:,:,1:5),status)
    if (status==0 .or. other%valid()) error stop 'bad exp shape'
    back = rigid_grid_from_matrix(mats)
    other = grid%compose(rigid_grid_identity([3,1]))
    if (other%valid()) error stop 'bad compose'
    a = grid%apply(reshape([(real(i,dp),i=1,12)], [4,3]))
    if (any(shape(a)/=[0,0,3])) error stop 'invalid apply shape'
    avg = grid%mean(weights=reshape([0.0_dp],[1,1]),status=status)
    if (status==0 .or. avg%valid()) error stop 'zero weights'
    avg = grid%mean(weights=reshape([ieee_value(0.0_dp,ieee_quiet_nan)],[1,1]),status=status)
    if (status==0 .or. avg%valid()) error stop 'NaN weights'
    mean_row = grid%mean_axis(1,weights=reshape([0.0_dp,1.0_dp,0.0_dp], [1,3]),status=status)
    if (status==0 .or. mean_row%valid()) error stop 'empty weighted slice'
    mean_row = grid%mean_axis(3,status=status)
    if (status==0 .or. mean_row%valid()) error stop 'bad axis'
    other = grid%power(ieee_value(0.0_dp,ieee_quiet_nan),status=status)
    if (status==0 .or. other%valid()) error stop 'bad exponent'
    print '(A)', 'PASS: 2D rigid grids, broadcasting, SE(3) conversions, reduction, validation'
contains
    subroutine check3(v,w,tol)
        real(dp), intent(in) :: v(:) !! actual vector
        real(dp), intent(in) :: w(:) !! reference vector
        real(dp), intent(in) :: tol !! tolerance
        if (size(v)/=3 .or. size(w)/=3) error stop 'vector dimensions'
        if (.not. all(abs(v-w)<tol)) then
            print *, 'actual ',v, 'expected ',w
            error stop 'vector mismatch'
        end if
    end subroutine check3
    subroutine checkmat(v,w,tol)
        real(dp), intent(in) :: v(:,:) !! actual matrix
        real(dp), intent(in) :: w(:,:) !! reference matrix
        real(dp), intent(in) :: tol !! tolerance
        if (any(shape(v)/=[4,4]) .or. any(shape(w)/=[4,4])) error stop 'matrix dimensions'
        if (.not. all(abs(v-w)<tol)) then
            print *, 'max matrix discrepancy ', maxval(abs(v-w))
            error stop 'matrix mismatch'
        end if
    end subroutine checkmat
    subroutine check4(v,w,tol)
        real(dp), intent(in) :: v(:,:,:,:) !! actual matrix grid
        real(dp), intent(in) :: w(:,:,:,:) !! reference matrix grid
        real(dp), intent(in) :: tol !! tolerance
        if (any(shape(v)/=shape(w))) error stop 'matrix grid dimensions'
        if (.not. all(abs(v-w)<tol)) then
            print *, 'max matrix grid discrepancy ', maxval(abs(v-w))
            error stop 'matrix grid mismatch'
        end if
    end subroutine check4
end program test_spatial_rigid_transform_grid
