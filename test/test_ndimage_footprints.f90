! SPDX-License-Identifier: MIT
! Tests for custom structuring elements, grayscale morphology, EDT sites.
program test_ndimage_footprints
    use scifort_kinds, only : dp
    use scifort_ndimage, only : binary_dilation_nd, binary_erosion_nd, &
        binary_opening_nd, binary_hit_or_miss_nd, binary_propagation_nd, &
        label_nd, grey_dilation_nd, grey_erosion_nd, grey_opening_nd, &
        grey_closing_nd, morphological_gradient_nd, morphological_laplace_nd, &
        white_tophat_nd, black_tophat_nd, distance_transform_edt_nd, &
        ndimage_status_invalid_input, ndimage_status_too_expensive
    implicit none
    logical :: x(5),y(5),foot(3),too_big(226,226),allzero(3)
    logical :: a(3,4),b(3,4),symmetric(3,3),asymmetric(3,3)
    logical :: high(3,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    logical :: high_y(3,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    logical :: high_fp(3,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: stat,nfeat,labels(3,4),indices(1,5),features15(15,6)
    real(dp) :: signal(5),gd(5),ge(5),work(5),r(5),d(5)
    real(dp) :: high_d(3,2,1,1,1,1,1,1,1,1,1,1,1,1,1)

    x=[.false.,.false.,.true.,.false.,.false.]
    foot=[.true.,.false.,.false.]
    call binary_dilation_nd(x,y,stat,structure=foot)
    if (stat/=0 .or. any(y.neqv.[.false.,.true.,.false.,.false.,.false.])) &
        error stop 'asymmetric dilation 1D'
    call binary_erosion_nd(.not.x,y,stat,structure=foot)
    if (stat/=0 .or. any(y.neqv.[.false.,.true.,.true.,.false.,.true.])) &
        error stop 'asymmetric erosion 1D'
    call binary_opening_nd(x,y,stat,structure=foot)
    if (stat/=0) error stop 'asymmetric opening'
    call binary_hit_or_miss_nd(x,y,stat,structure1=[.false.,.true.,.false.])
    if (stat/=0 .or. any(y.neqv.x)) error stop 'hit or miss isolated point'
    call binary_hit_or_miss_nd(x,y,stat,structure1=[.true.,.true.,.true.])
    if (stat/=0 .or. any(y)) error stop 'hit or miss full foreground stencil'
    call binary_propagation_nd(x,.true..eqv.x,y,stat,structure=foot)
    if (stat/=0) error stop 'custom propagation'
    call binary_dilation_nd(x,y,stat,structure=reshape([.true.,.false.,.true.,.false.],[2,2]))
    if (stat/=ndimage_status_invalid_input) error stop 'invalid footprint rank'
    too_big=.false.
    call binary_dilation_nd(x,y,stat,structure=too_big)
    if (stat/=ndimage_status_invalid_input) error stop 'invalid footprint rank precedence'
    a=.false.
    a(2,1)=.true.
    a(3,2)=.true.
    symmetric=.false.
    symmetric(1,1)=.true.
    symmetric(2,2)=.true.
    symmetric(3,3)=.true.
    call label_nd(a,labels,nfeat,stat,structure=symmetric)
    if (stat/=0 .or. nfeat/=1) error stop 'label custom diagonal connectivity'
    asymmetric=.false.
    asymmetric(1,1)=.true.
    call label_nd(a,labels,nfeat,stat,structure=asymmetric)
    if (stat/=ndimage_status_invalid_input) error stop 'invalid asymmetric labeling'
    call binary_dilation_nd(a,b,stat,structure=too_big)
    if (stat/=ndimage_status_too_expensive) error stop 'oversized binary footprint'

    signal=[2.0_dp,-1.0_dp,4.0_dp,1.0_dp,3.0_dp]
    call grey_dilation_nd(signal,gd,stat,sizes=[3])
    if (stat/=0 .or. any(gd/=[2.0_dp,4.0_dp,4.0_dp,4.0_dp,3.0_dp])) &
        error stop 'grey dilation values'
    call grey_erosion_nd(signal,ge,stat,sizes=[3])
    if (stat/=0 .or. any(ge/=[-1.0_dp,-1.0_dp,-1.0_dp,1.0_dp,1.0_dp])) &
        error stop 'grey erosion values'
    call morphological_gradient_nd(signal,r,stat,sizes=[3])
    if (stat/=0 .or. any(r/=(gd-ge))) error stop 'morphological gradient identity'
    call morphological_laplace_nd(signal,r,stat,sizes=[3])
    if (stat/=0 .or. any(r/=(gd+ge-2.0_dp*signal))) error stop 'morphological Laplace identity'
    call grey_opening_nd(signal,work,stat,sizes=[3])
    if (stat/=0) error stop 'grey opening'
    call white_tophat_nd(signal,r,stat,sizes=[3])
    if (stat/=0 .or. any(r/=(signal-work))) error stop 'white top-hat identity'
    call grey_closing_nd(signal,work,stat,sizes=[3])
    if (stat/=0) error stop 'grey closing'
    call black_tophat_nd(signal,r,stat,sizes=[3])
    if (stat/=0 .or. any(r/=(work-signal))) error stop 'black top-hat identity'
    allzero=.false.
    call grey_erosion_nd(signal,r,stat,footprint=allzero)
    if (stat/=ndimage_status_invalid_input) error stop 'empty grayscale footprint'
    call grey_dilation_nd(signal,r,stat,sizes=[huge(1)])
    if (stat/=ndimage_status_too_expensive) error stop 'oversized grayscale footprint'
    call grey_dilation_nd(signal,r,stat,sizes=[3],mode='invalid')
    if (stat/=ndimage_status_invalid_input) error stop 'invalid grayscale mode'
    call grey_erosion_nd(signal,r,stat,sizes=[3,3])
    if (stat/=ndimage_status_invalid_input) error stop 'invalid grayscale size rank'

    x=[.true.,.true.,.false.,.true.,.true.]
    call distance_transform_edt_nd(x,d,stat,feature_indices=indices)
    if (stat/=0 .or. any(indices(1,:)/=[3,3,3,3,3])) &
        error stop 'EDT feature indices one-based'
    call distance_transform_edt_nd(x,d,stat,feature_indices=features15)
    if (stat/=ndimage_status_invalid_input) error stop 'invalid EDT feature array shape'
    high=.true.
    high(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=.false.
    call distance_transform_edt_nd(high,high_d,stat,feature_indices=features15)
    if (stat/=0) error stop 'rank-15 EDT site features'
    if (features15(1,1)/=2 .or. any(features15(2:15,1)/=1)) &
        error stop 'rank-15 EDT nearest site'
    high=.false.
    high(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=.true.
    high_fp=.false.
    high_fp(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=.true.
    call binary_dilation_nd(high,high_y,stat,structure=high_fp)
    if (stat/=0 .or. .not.high_y(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)) &
        error stop 'rank-15 custom structuring element'
    print '(a)', 'PASS: custom footprints, grayscale morphology, feature indices and rank15'
end program test_ndimage_footprints
