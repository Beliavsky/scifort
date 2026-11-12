! SPDX-License-Identifier: MIT
! Default and full-connectivity binary morphology and labeling.
program test_ndimage_morphology
    use scifort_ndimage, only : binary_dilation_nd, binary_erosion_nd, &
        label_nd, morphology_status_success, morphology_status_invalid_input
    implicit none
    logical :: a(5,6),b(5,6),c(5,6),v(5,6)
    logical :: mask15(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)
    logical :: dilated15(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: labels(5,6), labels15(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: label_strided(10,12)
    logical :: mask_strided(10,12)
    integer :: status,n

    a=.false.
    a(2,2)=.true.
    a(2,3)=.true.
    a(5,6)=.true.
    call label_nd(a,labels,n,status)
    if (status/=morphology_status_success .or. n/=2) error stop 'face labels'
    if (labels(2,2)/=labels(2,3)) error stop 'same component'
    if (labels(5,6)==labels(2,2)) error stop 'distinct components'
    if (labels(1,1)/=0) error stop 'background label'
    call binary_dilation_nd(a,b,status)
    if (status/=morphology_status_success) error stop 'binary dilation status'
    if (count(b)<=count(a)) error stop 'dilation cardinality'
    call binary_dilation_nd(a,mask_strided(1:10:2,1:12:2),status)
    if (status/=morphology_status_success) error stop 'strided dilation status'
    if (any(b.neqv.mask_strided(1:10:2,1:12:2))) &
        error stop 'strided dilation values'
    call label_nd(a,label_strided(1:10:2,1:12:2),n,status)
    if (status/=morphology_status_success .or. &
        any(labels/=label_strided(1:10:2,1:12:2))) &
        error stop 'strided labels'
    call binary_erosion_nd(b,c,status)
    if (status/=morphology_status_success) error stop 'binary erosion status'
    call binary_dilation_nd(a,v,status,connectivity=2)
    if (status/=morphology_status_success) error stop 'full dilation'
    if (any(b.and..not.v)) error stop 'cross subset of full dilation'

    mask15=.false.
    mask15(1,2,1,1,1,1,1,1,1,1,1,1,1,1,1)=.true.
    mask15(2,3,1,1,1,1,1,1,1,1,1,1,1,1,1)=.true.
    call label_nd(mask15,labels15,n,status)
    if (status/=morphology_status_success .or. n/=2) error stop 'rank15 labels'
    call label_nd(mask15,labels15,n,status,connectivity=15)
    if (status/=2) error stop 'rank15 enormous neighborhood guard'
    call binary_dilation_nd(mask15,dilated15,status)
    if (status/=morphology_status_success) error stop 'rank15 dilation'
    if (count(dilated15)<2) error stop 'rank15 results'
    call label_nd(a,labels,n,status,connectivity=0)
    if (status/=morphology_status_invalid_input) error stop 'invalid connectivity'
    call binary_dilation_nd(a,b,status,iterations=0)
    if (status/=morphology_status_success .or. .not.all(b)) &
        error stop 'dilation convergence for zero iterations'
    print '(a)', 'PASS: ndimage binary morphology and components'
end program test_ndimage_morphology
