! SPDX-License-Identifier: MIT
program test_spatial_transform_nd
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_from_quat
    use scifort_spatial_rigid_transform, only : rigid_transform_type
    use scifort_spatial_transform_nd
    implicit none
    type(rotation_nd_type) :: r1, r2, rr, back, ident, bad
    type(rigid_nd_type) :: t1, t2, tt, trev, tback, tp, invalid_t
    type(rotation_type) :: rot, ra, rb, checkrot
    type(rigid_transform_type) :: tf, tf2
    real(dp) :: qa(2,1,3,4), qb(1,4,1,4), trans(2,1,3,3)
    real(dp) :: twists(1,4,1,6), points(1,4,1,3), qbad(3,1,3,4)
    real(dp) :: mats3(2,4,3,3,3), mats4(2,4,3,4,4)
    real(dp) :: nan
    real(dp) :: q15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,4)
    real(dp) :: twists15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,6)
    real(dp), allocatable :: out(:,:), dq(:,:), vectors(:,:)
    integer, allocatable :: dims(:)
    integer :: i,j,k,code
    real(dp), parameter :: tol=2.0e-12_dp

    do k=1,3
        do i=1,2
            rot = rotation_from_quat([0.1_dp*real(i,dp),0.02_dp,0.08_dp*real(k,dp),1.0_dp])
            qa(i,1,k,:) = rot%as_quat()
            trans(i,1,k,:) = [real(i,dp),0.4_dp,real(k,dp)*0.5_dp]
        end do
    end do
    do j=1,4
        rot = rotation_from_quat([0.0_dp,0.04_dp*real(j,dp),0.08_dp,1.0_dp])
        qb(1,j,1,:) = rot%as_quat()
        twists(1,j,1,:) = [0.0_dp,0.01_dp*real(j,dp),0.02_dp,1.0_dp,-0.25_dp,0.3_dp]
        points(1,j,1,:) = [0.5_dp,real(j,dp),-0.1_dp]
    end do

    r1=rotation_nd_from_quat(qa,status=code)
    if (code/=0 .or. .not. r1%valid() .or. any(r1%shape()/=[2,1,3])) error stop 'rotation constructor'
    r2=rotation_nd_from_quat(qb)
    rr=r1%compose(r2,status=code)
    if (code/=0 .or. .not. rr%valid() .or. any(rr%shape()/=[2,4,3])) error stop 'rotation broadcasting'
    if (rr%size()/=24) error stop 'rotation count'
    out=rr%apply(points,result_shape=dims,status=code)
    if (code/=0 .or. any(dims/=[2,4,3]) .or. any(shape(out)/=[24,3])) error stop 'rotation apply shape'
    kloop: do k=1,3
        do j=1,4
            do i=1,2
                ra=r1%at([i,1,k]); rb=r2%at([1,j,1]); rot=rr%at([i,j,k])
                checkrot=ra%compose(rb)
                if (maxval(abs(rot%as_matrix()-checkrot%as_matrix()))>tol) error stop 'rotation composition'
                if (maxval(abs(out(i+2*(j-1)+8*(k-1),:)-rot%apply(points(1,j,1,:))))>tol) &
                    error stop 'rotation application'
                mats3(i,j,k,:,:)=rot%as_matrix()
            end do
        end do
    end do kloop
    back=rotation_nd_from_matrix(mats3,status=code)
    if (code/=0 .or. .not. back%valid()) error stop 'rotation matrix roundtrip'
    if (maxval(abs(back%as_matrix()-rr%as_matrix()))>tol) error stop 'rotation matrix values'
    back=rr%compose(rr%inverse())
    if (maxval(abs(back%as_rotvec()))>tol) error stop 'rotation inverse'
    out=rr%apply([1.0_dp,0.0_dp,0.0_dp],status=code)
    if (code/=0 .or. any(shape(out)/=[24,3])) error stop 'rotation scalar point'
    rot=rr%at([3,4,2])
    if (rot%valid()) error stop 'rotation invalid index'
    qbad = 0.0_dp
    qbad(:,:,:,4) = 1.0_dp
    bad=rotation_nd_from_quat(qbad)
    back=r1%compose(bad,status=code)
    if (back%valid() .or. code/=spatial_nd_invalid_input) error stop 'incompatible rotations'
    nan=ieee_value(0.0_dp,ieee_quiet_nan)
    qbad(1,1,1,1)=nan
    bad=rotation_nd_from_quat(qbad,status=code)
    if (bad%valid() .or. code/=spatial_nd_invalid_input) error stop 'non-finite quaternion'
    bad=rotation_nd_identity([2,3,4])
    if (.not. bad%valid() .or. bad%size()/=24) error stop 'rotation identity'
    back=rotation_nd_from_rotvec(reshape(rr%as_rotvec(),[2,4,3,3]))
    if (.not. back%valid() .or. maxval(abs(back%as_matrix()-rr%as_matrix()))>tol) &
        error stop 'rotation-vector roundtrip'

    t1=rigid_nd_from_components(trans,r1,status=code)
    if (code/=0 .or. .not. t1%valid() .or. any(t1%shape()/=[2,1,3])) error stop 'rigid components'
    t2=rigid_nd_from_exp_coords(twists,status=code)
    if (code/=0 .or. .not. t2%valid() .or. any(t2%shape()/=[1,4,1])) error stop 'rigid exp coords'
    tt=t1%compose(t2,status=code)
    if (code/=0 .or. .not. tt%valid() .or. any(tt%shape()/=[2,4,3])) error stop 'rigid broadcasting'
    vectors=tt%apply(points,result_shape=dims,status=code)
    if (code/=0 .or. any(dims/=[2,4,3])) error stop 'rigid apply shape'
    do k=1,3
        do j=1,4
            do i=1,2
                tf=t1%at([i,1,k]); tf2=t2%at([1,j,1])
                tf=tf%compose(tf2)
                tf2=tt%at([i,j,k])
                if (maxval(abs(tf%as_matrix()-tf2%as_matrix()))>tol) error stop 'rigid composition'
                if (maxval(abs(vectors(i+2*(j-1)+8*(k-1),:)-tf%apply(points(1,j,1,:))))>tol) &
                    error stop 'rigid point application'
                mats4(i,j,k,:,:)=tf%as_matrix()
            end do
        end do
    end do
    tback=rigid_nd_from_matrix(mats4,status=code)
    if (code/=0 .or. .not. tback%valid()) error stop 'matrix input'
    if (maxval(abs(tback%as_matrix()-tt%as_matrix()))>tol) error stop 'rigid matrix roundtrip'
    dq=tt%as_dual_quat()
    tback=rigid_nd_from_dual_quat(reshape(dq,[2,4,3,8]),status=code)
    if (code/=0 .or. maxval(abs(tback%as_matrix()-tt%as_matrix()))>tol) error stop 'dual quaternion roundtrip'
    tback=rigid_nd_from_exp_coords(reshape(tt%as_exp_coords(),[2,4,3,6]),status=code)
    if (code/=0 .or. maxval(abs(tback%as_matrix()-tt%as_matrix()))>tol) error stop 'SE3 exp coords roundtrip'
    trev=tt%compose(tt%inverse())
    do i=1,trev%size()
        tf=trev%at([mod(i-1,2)+1,mod((i-1)/2,4)+1,(i-1)/8+1])
        if (maxval(abs(tf%translation()))>tol) error stop 'rigid inverse translation'
        rot=tf%rotation()
        if (maxval(abs(rot%as_rotvec()))>tol) error stop 'rigid inverse rotation'
    end do
    tp=tt%power(0.5_dp,status=code)
    if (code/=0 .or. .not. tp%valid()) error stop 'rigid power'
    tback=tp%compose(tp)
    if (maxval(abs(tback%as_matrix()-tt%as_matrix()))>1.0e-11_dp) error stop 'rigid power reconstruction'
    vectors=tt%apply([1.0_dp,2.0_dp,3.0_dp],status=code)
    if (code/=0 .or. any(shape(vectors)/=[24,3])) error stop 'rigid scalar point'
    invalid_t=rigid_nd_identity([1,2])
    invalid_t=t1%compose(invalid_t,status=code)
    if (invalid_t%valid() .or. code/=spatial_nd_invalid_input) error stop 'rigid incompatible shapes'
    tf=tt%at([0,1,1])
    if (tf%valid()) error stop 'rigid invalid index'
    invalid_t=rigid_nd_from_exp_coords([nan,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp],status=code)
    if (invalid_t%valid() .or. code/=spatial_nd_invalid_input) error stop 'rigid nonfinite input'
    invalid_t=rigid_nd_from_matrix(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp],[2,2]),status=code)
    if (invalid_t%valid() .or. code/=spatial_nd_invalid_input) error stop 'invalid matrix dimensions'
    ident=rotation_nd_from_quat([0.0_dp,0.0_dp,0.0_dp,1.0_dp])
    if (.not. ident%valid() .or. size(ident%shape())/=0) error stop 'scalar quaternion'
    out=ident%apply([1.0_dp,0.0_dp,0.0_dp],result_shape=dims,status=code)
    if (code/=0 .or. size(dims)/=0 .or. maxval(abs(out(1,:)-[1.0_dp,0.0_dp,0.0_dp]))>tol) &
        error stop 'scalar rotation application'
    q15=0.0_dp
    q15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,4)=1.0_dp
    ident=rotation_nd_from_quat(q15,status=code)
    if (code/=0 .or. size(ident%shape())/=14 .or. ident%size()/=1) error stop 'rank15 quaternion'
    twists15=0.0_dp
    invalid_t=rigid_nd_from_exp_coords(twists15,status=code)
    if (code/=0 .or. size(invalid_t%shape())/=14 .or. invalid_t%size()/=1) error stop 'rank15 twist'
    print *, 'transform ND PASS' 
end program test_spatial_transform_nd
