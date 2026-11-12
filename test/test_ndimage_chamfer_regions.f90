! SPDX-License-Identifier: MIT
! Focused regression checks for chamfer distances, extrema footprints, and histograms.
program test_ndimage_chamfer_regions
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_kinds, only : dp
    use scifort_ndimage, only : minimum_filter_nd, maximum_filter_nd, &
        distance_transform_cdt_nd, count_labels_nd, histogram_labels_nd, &
        standard_deviation_labels_nd, ndimage_status_invalid_input
    implicit none
    logical :: mask(5),foot(3),empty(3),mask2(3,4)
    real(dp) :: a(5),b(5),image(2,3),sds(3),expected
    integer :: d(5),features(1,5),d2(3,4),indices(2,12),status,i,j
    integer :: labels(2,3),counts(3),hist(3,3),ids(3)
    logical :: mask15(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: distance15(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1),features15(15,4)

    a=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
    foot=[.true.,.false.,.false.]
    call minimum_filter_nd(a,y=b,status=status,mode='nearest',footprint=foot)
    if (status/=0 .or. any(b/=[1.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp])) &
        error stop 'minimum asymmetric footprint'
    call maximum_filter_nd(a,y=b,status=status,mode='nearest',footprint=foot)
    if (status/=0 .or. any(b/=[1.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp])) &
        error stop 'maximum footprint orientation'
    call maximum_filter_nd(a,2,b,status,mode='nearest',origin=[-1])
    if (status/=0 .or. any(b/=[2.0_dp,3.0_dp,4.0_dp,5.0_dp,5.0_dp])) &
        error stop 'even maximum origin'
    empty=.false.
    call minimum_filter_nd(a,y=b,status=status,footprint=empty)
    if (status/=ndimage_status_invalid_input) error stop 'empty stencil'
    call maximum_filter_nd(a,y=b,status=status,footprint=foot,origin=[2])
    if (status/=ndimage_status_invalid_input) error stop 'invalid origin'

    mask=[.true.,.true.,.false.,.true.,.true.]
    call distance_transform_cdt_nd(mask,d,status,feature_indices=features)
    if (status/=0 .or. any(d/=[2,1,0,1,2])) error stop 'taxicab 1D'
    if (any(features(1,:)/=3)) error stop 'taxicab feature index'
    call distance_transform_cdt_nd(mask,d,status,metric='chessboard')
    if (status/=0 .or. any(d/=[2,1,0,1,2])) error stop 'chessboard 1D'
    mask2=.true.
    mask2(1,1)=.false.
    call distance_transform_cdt_nd(mask2,d2,status,metric='taxicab',feature_indices=indices)
    if (status/=0) error stop 'taxicab 2D status'
    do j=1,4
        do i=1,3
            if (d2(i,j)/=i+j-2) error stop 'taxicab 2D distance'
            if (any(indices(:,i+3*(j-1))/=[1,1])) error stop 'taxicab 2D index'
        end do
    end do
    call distance_transform_cdt_nd(mask2,d2,status,metric='chessboard',feature_indices=indices)
    if (status/=0) error stop 'chessboard 2D status'
    do j=1,4
        do i=1,3
            if (d2(i,j)/=max(i-1,j-1)) error stop 'chessboard 2D distance'
            if (any(indices(:,i+3*(j-1))/=[1,1])) error stop 'chessboard 2D index'
        end do
    end do
    mask2=.true.
    call distance_transform_cdt_nd(mask2,d2,status)
    if (status/=0 .or. any(d2/=-1)) error stop 'all foreground chamfer'
    call distance_transform_cdt_nd(mask2,d2,status,metric='euclidean')
    if (status/=ndimage_status_invalid_input) error stop 'invalid chamfer metric'

    image=reshape([0.0_dp,2.0_dp,4.0_dp,6.0_dp,3.0_dp,5.0_dp],[2,3])
    labels=reshape([1,1,2,2,1,2],[2,3])
    ids=[1,2,99]
    call count_labels_nd(labels,ids,counts,status)
    if (status/=0 .or. any(counts/=[3,3,0])) error stop 'count labels'
    call histogram_labels_nd(image,labels,ids,0.0_dp,6.0_dp,hist,status)
    if (status/=0) error stop 'histogram status'
    if (any(hist(:,1)/=[1,2,0])) error stop 'first label histogram'
    if (any(hist(:,2)/=[0,0,3])) error stop 'second label histogram'
    if (any(hist(:,3)/=0)) error stop 'missing histogram zeros'
    call standard_deviation_labels_nd(image,labels,ids,sds,status)
    if (status/=0) error stop 'standard deviation status'
    expected=sqrt(14.0_dp/9.0_dp)
    if (abs(sds(1)-expected)>1.0e-13_dp) error stop 'standard deviation value'
    if (.not.ieee_is_nan(sds(3))) error stop 'missing standard deviation'
    call histogram_labels_nd(image,labels,ids,3.0_dp,3.0_dp,hist,status)
    if (status/=ndimage_status_invalid_input) error stop 'invalid histogram range'
    call count_labels_nd(labels,ids,counts(:2),status)
    if (status/=ndimage_status_invalid_input) error stop 'invalid count output'

    mask15=.true.
    mask15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=.false.
    call distance_transform_cdt_nd(mask15,distance15,status,metric='chessboard',feature_indices=features15)
    if (status/=0) error stop 'rank15 chamfer status'
    if (distance15(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)/=1) error stop 'rank15 chebyshev'
    if (any(features15(:,4)/=1)) error stop 'rank15 chamfer features'
    print '(a)', 'PASS: extrema footprints, chamfer and selected-region histograms'
end program test_ndimage_chamfer_regions
