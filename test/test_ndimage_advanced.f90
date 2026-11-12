! SPDX-License-Identifier: MIT
! Regression tests for shifted stencils, non-flat morphology, measurements, and IFT.
program test_ndimage_advanced
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_ndimage, only : grey_dilation_nd, grey_erosion_nd, &
        grey_opening_nd, binary_dilation_nd, binary_hit_or_miss_nd, &
        sum_labels_nd, mean_labels_nd, variance_labels_nd, center_of_mass_nd, &
        watershed_ift_nd, ndimage_status_invalid_input, ndimage_status_too_expensive
    implicit none
    real(dp) :: x(5),y(5),nonflat(3),sums(3),means(3),vars(3),centers(1,3)
    real(dp) :: two(2),cent1(1,1)
    logical :: mask(5),result(5),struct(3)
    integer :: labels(5),region(5),ift(7),intensity(7),seeds(7)
    integer :: high(3,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: high_mark(3,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: high_labels(3,2,1,1,1,1,1,1,1,1,1,1,1,1,1)
    logical :: bad_struct(3,3,3,3,3,3,3,3,3,3)
    integer :: small10(1,1,1,1,1,1,1,1,1,1)
    integer :: markers10(1,1,1,1,1,1,1,1,1,1)
    integer :: result10(1,1,1,1,1,1,1,1,1,1)
    integer :: st

    x=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
    nonflat=[0.0_dp,3.0_dp,-1.0_dp]
    call grey_dilation_nd(x,y,st,structure=nonflat,origin=[1],mode='nearest')
    if (st/=0 .or. maxval(abs(y-[5.0_dp,6.0_dp,7.0_dp,8.0_dp,8.0_dp]))>1e-12_dp) &
        error stop 'non-flat dilation origin'
    call grey_erosion_nd(x,y,st,structure=nonflat,origin=[-1],mode='nearest')
    if (st/=0) error stop 'non-flat erosion shifted origin'
    call grey_opening_nd(x,y,st,structure=nonflat,origin=[1])
    if (st/=0) error stop 'non-flat opening shifted origin'
    call grey_dilation_nd(x,y,st,sizes=[3],origin=[2])
    if (st/=ndimage_status_invalid_input) error stop 'origin outside footprint'
    call grey_dilation_nd(x,y,st,structure=reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp],[2,2]))
    if (st/=ndimage_status_invalid_input) error stop 'grayscale structure rank'
    call grey_dilation_nd(x,y,st,structure=[1.0_dp,2.0_dp],footprint=[.true.,.true.,.true.])
    if (st/=ndimage_status_invalid_input) error stop 'mismatched gray shapes'
    nonflat(1)=ieee_value(0.0_dp,ieee_quiet_nan)
    call grey_dilation_nd(x,y,st,structure=nonflat)
    if (st/=ndimage_status_invalid_input) error stop 'nonfinite structure'

    mask=.false.
    mask(4)=.true.
    struct=[.true.,.false.,.false.]
    call binary_dilation_nd(mask,result,st,structure=struct,origin=[1])
    if (st/=0 .or. any(result.neqv.[.false.,.true.,.false.,.false.,.false.])) &
        error stop 'binary dilation shifted origin'
    call binary_dilation_nd(mask,result,st,structure=struct,origin=[2])
    if (st/=ndimage_status_invalid_input) error stop 'invalid binary origin'
    call binary_hit_or_miss_nd(mask,result,st,structure1=struct,origin1=[1])
    if (st/=0) error stop 'hit-or-miss origin'

    labels=[1,1,2,2,0]
    region=[1,1,2,2,0]
    call sum_labels_nd(x,region,[1,2,9],sums,st)
    if (st/=0 .or. any(sums/=[3.0_dp,7.0_dp,0.0_dp])) error stop 'label sums'
    call mean_labels_nd(x,region,[1,2,9],means,st)
    if (st/=0 .or. maxval(abs(means(:2)-[1.5_dp,3.5_dp]))>1e-12_dp) &
        error stop 'label means'
    if (.not.ieee_is_nan(means(3))) error stop 'missing mean NaN'
    call variance_labels_nd(x,region,[1,2,9],vars,st)
    if (st/=0 .or. maxval(abs(vars(:2)-0.25_dp))>1e-12_dp) &
        error stop 'label variance'
    if (.not.ieee_is_nan(vars(3))) error stop 'missing variance NaN'
    call center_of_mass_nd(x,region,[1,2,9],centers,st)
    if (st/=0 .or. abs(centers(1,1)-5.0_dp/3.0_dp)>1e-12_dp) &
        error stop 'weighted one-based COM'
    if (.not.ieee_is_nan(centers(1,3))) error stop 'missing center NaN'
    two=[-1.0_dp,1.0_dp]
    call center_of_mass_nd(two,[1,1],[1],cent1,st)
    if (st/=0 .or. cent1(1,1)<huge(1.0_dp)) error stop 'zero-mass center infinity'
    call sum_labels_nd(x,region,[1,2,9],sums(:2),st)
    if (st/=ndimage_status_invalid_input) error stop 'wrong total output length'

    intensity=[0,1,2,3,4,5,6]
    seeds=[1,0,0,0,0,0,2]
    call watershed_ift_nd(intensity,seeds,ift,st)
    if (st/=0 .or. any(ift/=[1,1,1,1,1,2,2])) error stop 'watershed 1D'
    call watershed_ift_nd([-1,0,1,2,3,4,5],seeds,ift,st)
    if (st/=ndimage_status_invalid_input) error stop 'negative watershed intensity'
    call watershed_ift_nd([65536,0,1,2,3,4,5],seeds,ift,st)
    if (st/=ndimage_status_invalid_input) error stop 'oversized watershed intensity'
    high=7
    high_mark=0
    high_mark(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=17
    call watershed_ift_nd(high,high_mark,high_labels,st)
    if (st/=0 .or. any(high_labels/=17)) error stop 'rank15 watershed'
    bad_struct=.true.
    small10=1
    markers10=0
    call watershed_ift_nd(small10,markers10,result10,st,structure=bad_struct)
    if (st/=ndimage_status_too_expensive) error stop 'large watershed stencil guard'
    print '(a)', 'PASS: advanced ndimage grayscale, binary, measurements, watershed rank15'
end program test_ndimage_advanced
