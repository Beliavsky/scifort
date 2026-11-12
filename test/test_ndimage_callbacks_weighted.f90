! SPDX-License-Identifier: MIT
! Native regression tests for generic kernels and weighted distance transforms.
program test_ndimage_callbacks_weighted
    use scifort_kinds, only : dp
    use scifort_ndimage, only : correlate1d_nd, convolve1d_nd, generic_filter_nd, &
        labeled_comprehension_nd, distance_transform_cdt_nd, distance_transform_weighted_cdt_nd, &
        ndimage_status_invalid_input
    implicit none
    real(dp) :: a(5),b(5),aa(2,2),bb(2,2),result(3),labelimage(5)
    complex(dp) :: z(5),zr(5),c(1)
    integer :: status,labels(5),dist(5),features(1,5),weights1(3),missing(3),j
    integer :: weights2(3,3),dist2(3,4),feat2(2,12),weighted(2)
    logical :: mask(5),mask2(3,4),foot(2,2)
    logical :: mask15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: x15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: y15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: d15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1),cost15(15)

    a=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
    call correlate1d_nd(a,[1.0_dp,0.0_dp],1,b,status,mode='nearest')
    if (status/=0 .or. any(b/=[1.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp])) &
        error stop 'correlate1d even kernel'
    call convolve1d_nd(a,[1.0_dp,0.0_dp],1,b,status,mode='nearest')
    if (status/=0 .or. any(b/=[2.0_dp,3.0_dp,4.0_dp,5.0_dp,5.0_dp])) &
        error stop 'convolve1d even kernel'
    call correlate1d_nd(a,[1.0_dp],2,b,status)
    if (status/=ndimage_status_invalid_input) error stop 'invalid axis'
    call correlate1d_nd(a,[1.0_dp,2.0_dp],1,b,status,origin=1)
    if (status/=ndimage_status_invalid_input) error stop 'invalid origin'
    call convolve1d_nd(a,[1.0_dp],1,b,status,mode='bad')
    if (status/=ndimage_status_invalid_input) error stop 'invalid extension'
    z=cmplx(a,0.0_dp,dp)
    c=[cmplx(2.0_dp,3.0_dp,dp)]
    call correlate1d_nd(z,c,1,zr,status)
    if (status/=0 .or. any(abs(zr-conjg(c(1))*z)>1.0e-14_dp)) &
        error stop 'complex coefficient conjugation'
    call convolve1d_nd(z,c,1,zr,status)
    if (status/=0 .or. any(abs(zr-c(1)*z)>1.0e-14_dp)) &
        error stop 'complex convolution coefficient'

    aa=reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp],[2,2])
    call generic_filter_nd(aa,weighted_sum,bb,status,sizes=[2,2],origin=[-1,-1], &
        mode='constant',cval=0.0_dp)
    if (status/=0 .or. abs(bb(1,1)-29.0_dp)>1.0e-13_dp) &
        error stop 'generic callback must receive C-order pixels'
    foot=reshape([.true.,.false.,.true.,.false.],[2,2])
    call generic_filter_nd(aa,weighted_sum,bb,status,footprint=foot,origin=[-1,0])
    if (status/=0) error stop 'sparse generic footprint'
    foot=.false.
    call generic_filter_nd(aa,weighted_sum,bb,status,footprint=foot)
    if (status/=ndimage_status_invalid_input) error stop 'empty generic footprint'
    call generic_filter_nd(aa,weighted_sum,bb,status)
    if (status/=ndimage_status_invalid_input) error stop 'missing generic footprint'

    labelimage=[1.0_dp,2.0_dp,7.0_dp,4.0_dp,5.0_dp]
    labels=[1,1,2,1,2]
    call labeled_comprehension_nd(labelimage,labels,[1,2,9],weighted_sum, &
        -17.0_dp,result,status)
    if (status/=0 .or. any(abs(result-[17.0_dp,17.0_dp,-17.0_dp])>1.0e-14_dp)) &
        error stop 'labeled callback values'
    call labeled_comprehension_nd(labelimage,labels,[1,2],weighted_sum, &
        -17.0_dp,result,status)
    if (status/=ndimage_status_invalid_input) error stop 'invalid labeled output'
    call labeled_comprehension_nd(labelimage,labels,[integer ::],weighted_sum, &
        -17.0_dp,result(:0),status)
    if (status/=0) error stop 'empty labeled selection'

    mask=[.true.,.true.,.false.,.true.,.true.]
    call distance_transform_cdt_nd(mask,dist,status,axis_weights=[3],feature_indices=features)
    if (status/=0 .or. any(dist/=[6,3,0,3,6])) error stop 'weighted L1 distances'
    if (any(features(1,:)/=3)) error stop 'weighted L1 feature sites'
    call distance_transform_cdt_nd(mask,dist,status,metric='chessboard',axis_weights=[3])
    if (status/=ndimage_status_invalid_input) error stop 'weighted chessboard invalid'
    call distance_transform_cdt_nd(mask,dist,status,axis_weights=[0])
    if (status/=ndimage_status_invalid_input) error stop 'nonpositive axis cost'
    weights1=[3,0,5]
    call distance_transform_weighted_cdt_nd(mask,weights1,dist,status,features)
    if (status/=0 .or. any(dist/=[10,5,0,3,6])) error stop 'directed chamfer distance'
    if (any(features(1,:)/=3)) error stop 'directed chamfer witness'
    weights1=[0,1,1]
    call distance_transform_weighted_cdt_nd(mask,weights1,dist,status)
    if (status/=ndimage_status_invalid_input) error stop 'nonzero stencil center'
    weights1=[-1,0,1]
    call distance_transform_weighted_cdt_nd(mask,weights1,dist,status)
    if (status/=ndimage_status_invalid_input) error stop 'negative stencil step'
    missing=[0,0,0]
    call distance_transform_weighted_cdt_nd(mask,missing,dist,status,features)
    if (status/=0 .or. any(dist/=[-1,-1,0,-1,-1])) error stop 'disconnected stencil'
    if (any(features(:,[1,2,4,5])/=0)) error stop 'unreachable stencil sites'
    mask2=.true.
    mask2(2,2)=.false.
    weights2=4
    weights2(2,2)=0
    weights2(1,2)=3
    weights2(3,2)=3
    weights2(2,1)=3
    weights2(2,3)=3
    call distance_transform_weighted_cdt_nd(mask2,weights2,dist2,status,feat2)
    if (status/=0 .or. dist2(3,3)/=4 .or. dist2(3,4)/=7) &
        error stop '3-4 chamfer 2D'
    if (any(feat2(:,12)/=[2,2])) error stop '3-4 chamfer feature'
    weighted=[2,3]
    call distance_transform_cdt_nd(mask2,dist2,status,axis_weights=weighted)
    if (status/=0 .or. dist2(3,4)/=8) error stop 'anisotropic weighted L1'

    mask15=.true.
    mask15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=.false.
    cost15=1
    cost15(1)=3
    call distance_transform_cdt_nd(mask15,d15,status,axis_weights=cost15)
    if (status/=0 .or. d15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/=3) &
        error stop 'rank15 weighted distances'
    x15=4.0_dp
    call generic_filter_nd(x15,weighted_sum,y15,status,sizes=spread(1,1,15))
    if (status/=0 .or. any(abs(reshape(y15,[size(y15)])-4.0_dp)>1.0e-14_dp)) &
        error stop 'rank15 generic callback'
    print '(a)', 'PASS: generic filtering, 1D kernels, weighted chamfer, labeled callbacks'
contains
    real(dp) function weighted_sum(values) result(y)
        real(dp), intent(in) :: values(:)
        integer :: i
        y=0.0_dp
        do i=1,size(values)
            y=y+real(i,dp)*values(i)
        end do
    end function weighted_sum
end program test_ndimage_callbacks_weighted
