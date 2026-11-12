! SPDX-License-Identifier: MIT
! New labeled-region, masked morphology and order-statistic regression tests.
program test_ndimage_regions_masks
    use scifort_kinds, only : dp
    use scifort_ndimage, only : extrema_labels_nd, minimum_labels_nd, &
        maximum_labels_nd, minimum_position_nd, maximum_position_nd, &
        bounding_boxes_nd, find_objects_nd, binary_dilation_nd, binary_erosion_nd, &
        binary_opening_nd, binary_closing_nd, median_filter_nd, rank_filter_nd, &
        percentile_filter_nd, ndimage_status_invalid_input
    implicit none
    real(dp) :: data(2,3),mins(3),maxs(3),ranked(2,3),med(2,3),percent(2,3)
    integer :: ids(3),labs(2,3),pmin(2,3),pmax(2,3),lo(2,3),hi(2,3)
    integer :: objects_lo(2,4),objects_hi(2,4),status
    logical :: input(5),editable(5),out(5),e(5),d(5),opened(5)
    logical :: wrong_mask(5,1),short_mask(4),foot(2,3),empty_foot(2,3)
    real(dp) :: high(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: high_out(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: base(2,2),base_out(2,2)
    logical :: high_foot(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: high_size(15),high_origin(15)
    integer :: high_labels(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: high_minpos(15,1),high_maxpos(15,1),high_lower(15,2),high_upper(15,2)
    real(dp) :: high_min(1),high_max(1)

    data=reshape([6.0_dp,0.0_dp,0.0_dp,3.0_dp,2.0_dp,2.0_dp],[2,3])
    labs=1
    labs(2,2)=3
    ids=[1,3,9]
    call extrema_labels_nd(data,labs,ids,mins,maxs,pmin,pmax,status)
    if (status/=0) error stop 'extrema status'
    if (any(mins/=[0.0_dp,3.0_dp,0.0_dp])) error stop 'minima values'
    if (any(maxs/=[6.0_dp,3.0_dp,0.0_dp])) error stop 'maxima values'
    if (any(pmin(:,1)/=[1,2])) error stop 'C-order minimum tie position'
    if (any(pmax(:,1)/=[1,1])) error stop 'maximum position'
    if (any(pmin(:,3)/=0) .or. any(pmax(:,3)/=0)) error stop 'missing extrema positions'
    call minimum_labels_nd(data,labs,ids,mins,status)
    if (status/=0 .or. any(mins/=[0.0_dp,3.0_dp,0.0_dp])) error stop 'min wrapper'
    call maximum_labels_nd(data,labs,ids,maxs,status)
    if (status/=0 .or. any(maxs/=[6.0_dp,3.0_dp,0.0_dp])) error stop 'max wrapper'
    call minimum_position_nd(data,labs,ids,pmin,status)
    if (status/=0 .or. any(pmin(:,1)/=[1,2])) error stop 'min pos wrapper'
    call maximum_position_nd(data,labs,ids,pmax,status)
    if (status/=0 .or. any(pmax(:,1)/=[1,1])) error stop 'max pos wrapper'
    call bounding_boxes_nd(labs,ids,lo,hi,status)
    if (status/=0 .or. any(lo(:,1)/=[1,1]) .or. any(hi(:,1)/=[2,3])) &
        error stop 'region bounding box'
    if (any(lo(:,3)/=0) .or. any(hi(:,3)/=0)) error stop 'missing region box'
    call find_objects_nd(labs,objects_lo,objects_hi,status)
    if (status/=0) error stop 'find objects status'
    if (any(objects_lo(:,1)/=[1,1]) .or. any(objects_hi(:,1)/=[2,3])) &
        error stop 'find objects first label'
    if (any(objects_lo(:,2)/=0) .or. any(objects_hi(:,2)/=0)) &
        error stop 'find objects missing label'
    if (any(objects_lo(:,3)/=[2,2]) .or. any(objects_hi(:,3)/=[2,2])) &
        error stop 'find objects third label'

    call extrema_labels_nd(data,reshape(labs,[6]),ids,mins,maxs,pmin,pmax,status)
    if (status/=ndimage_status_invalid_input) error stop 'extrema label rank mismatch'
    call bounding_boxes_nd(labs,ids,lo(1:1,:),hi,status)
    if (status/=ndimage_status_invalid_input) error stop 'invalid boxes output shape'
    call find_objects_nd(labs,objects_lo(1:1,:),objects_hi,status)
    if (status/=ndimage_status_invalid_input) error stop 'invalid objects output shape'

    input=[.false.,.false.,.true.,.false.,.false.]
    editable=[.true.,.false.,.true.,.false.,.true.]
    call binary_dilation_nd(input,out,status,mask=editable,iterations=3)
    if (status/=0 .or. any(out.neqv.input)) error stop 'frozen neighbors'
    call binary_dilation_nd(input,out,status,iterations=0)
    if (status/=0 .or. .not.all(out)) error stop 'dilation until stable'
    call binary_erosion_nd(.not.input,out,status,iterations=-1)
    if (status/=0 .or. any(out)) error stop 'erosion until stable'
    editable=[.true.,.false.,.true.,.true.,.true.]
    call binary_erosion_nd(input,e,status,mask=editable)
    if (status/=0) error stop 'masked erosion'
    call binary_dilation_nd(e,d,status,mask=editable)
    if (status/=0) error stop 'masked dilation'
    call binary_opening_nd(input,opened,status,mask=editable)
    if (status/=0 .or. any(opened.neqv.d)) error stop 'opening compositions under mask'
    call binary_closing_nd(input,out,status,mask=editable)
    if (status/=0) error stop 'masked closing'
    wrong_mask=.true.
    call binary_dilation_nd(input,out,status,mask=wrong_mask)
    if (status/=ndimage_status_invalid_input) error stop 'wrong mask rank'
    short_mask=.true.
    call binary_erosion_nd(input,out,status,mask=short_mask)
    if (status/=ndimage_status_invalid_input) error stop 'wrong mask shape'

    foot=reshape([.true.,.false.,.true.,.true.,.true.,.false.],[2,3])
    call median_filter_nd(data,[2,3],med,status,mode='nearest',footprint=foot,origin=[0,1])
    if (status/=0) error stop 'asymmetric median footprint'
    call median_filter_nd(data,y=percent,status=status,mode='nearest', &
        footprint=foot,origin=[0,1])
    if (status/=0 .or. any(percent/=med)) error stop 'footprint-only median'
    call rank_filter_nd(data,-2,y=percent,status=status,mode='nearest', &
        footprint=foot,origin=[0,1])
    if (status/=0) error stop 'footprint-only rank'
    call rank_filter_nd(data,-2,[2,3],ranked,status,mode='nearest',footprint=foot,origin=[0,1])
    if (status/=0) error stop 'asymmetric rank footprint'
    call percentile_filter_nd(data,50.0_dp,y=percent,status=status,mode='nearest', &
        footprint=foot,origin=[0,1])
    if (status/=0 .or. any(percent/=med)) error stop 'footprint-only percentile'
    call percentile_filter_nd(data,50.0_dp,[2,3],percent,status, &
        mode='nearest',footprint=foot,origin=[0,1])
    if (status/=0 .or. any(percent/=med)) error stop 'median percentile identity'
    empty_foot=.false.
    call median_filter_nd(data,[2,3],med,status,footprint=empty_foot)
    if (status/=ndimage_status_invalid_input) error stop 'empty footprint'
    call median_filter_nd(data,[2,3],med,status,footprint=foot,origin=[1,0])
    if (status/=ndimage_status_invalid_input) error stop 'invalid even-width origin'
    call rank_filter_nd(data,5,[2,3],ranked,status,footprint=foot)
    if (status/=ndimage_status_invalid_input) error stop 'rank exceeds active footprint'

    high=0.0_dp
    high(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=4.0_dp
    high(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=9.0_dp
    high(1,2,1,1,1,1,1,1,1,1,1,1,1,1,1)=1.0_dp
    high(2,2,1,1,1,1,1,1,1,1,1,1,1,1,1)=2.0_dp
    base=reshape([4.0_dp,9.0_dp,1.0_dp,2.0_dp],[2,2])
    high_labels=1
    call extrema_labels_nd(high,high_labels,[1],high_min,high_max, &
        high_minpos,high_maxpos,status)
    if (status/=0 .or. high_min(1)/=1.0_dp .or. high_max(1)/=9.0_dp) &
        error stop 'rank15 extrema values'
    if (any(high_minpos(1:2,1)/=[1,2]) .or. any(high_maxpos(1:2,1)/=[2,1])) &
        error stop 'rank15 extrema positions'
    call bounding_boxes_nd(high_labels,[1,8],high_lower,high_upper,status)
    if (status/=0) error stop 'rank15 bounding boxes'
    if (any(high_lower(1:2,1)/=[1,1]) .or. any(high_upper(1:2,1)/=[2,2])) &
        error stop 'rank15 nonempty bounds'
    if (any(high_lower(:,2)/=0) .or. any(high_upper(:,2)/=0)) &
        error stop 'rank15 missing bounds'
    high_foot=.true.
    high_size=[2,2,1,1,1,1,1,1,1,1,1,1,1,1,1]
    high_origin=0
    call rank_filter_nd(high,2,high_size,high_out,status,mode='nearest', &
        footprint=high_foot,origin=high_origin)
    if (status/=0) error stop 'rank15 rank filtering'
    call rank_filter_nd(base,2,[2,2],base_out,status,mode='nearest')
    if (status/=0 .or. any(reshape(high_out,[4])/=reshape(base_out,[4]))) &
        error stop 'rank15 filter singleton axes'
    print '(a)', 'PASS: extrema, bounds, masked morphology, rank15 order statistics'
end program test_ndimage_regions_masks
