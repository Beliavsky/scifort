! SPDX-License-Identifier: MIT
! Generate/iterate structure, numeric labels and selected derivative axes.
program test_ndimage_structure_axes
    use scifort_kinds, only : dp
    use scifort_ndimage, only : generate_binary_structure_nd, iterate_structure_nd, &
        label_nd, generic_gradient_magnitude_nd, generic_laplace_nd, &
        correlate1d_nd, morphology_status_invalid_input, &
        morphology_status_too_expensive, ndimage_status_invalid_input
    implicit none
    logical :: cross(3,3),full(3,3),diamond(5,5),sphere(5,5,5)
    logical :: cross3(3,3,3),edge(2),expanded(4),wrong(4,5)
    logical :: stencil15(3,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    logical :: extended15(5,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: ints(4,5),labels(4,5),labels_copy(4,5),n,s,i,j
    integer :: int15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    integer :: label15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    logical :: singleton(1),same_singleton(1)
    integer :: origins(2)
    real(dp) :: a(4,5),b(4,5),term(4,5),term2(4,5),reals(4,5)
    real(dp) :: x15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: y15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

    call generate_binary_structure_nd(1,cross,status=s)
    if (s/=0 .or. count(cross)/=5) error stop 'cross generation'
    if (.not.cross(2,2) .or. cross(1,1)) error stop 'cross stencil'
    call generate_binary_structure_nd(2,full,status=s)
    if (s/=0 .or. .not.all(full)) error stop 'full stencil'
    call generate_binary_structure_nd(0,full,status=s)
    if (s/=0 .or. any(full.neqv.cross)) error stop 'SciPy clamped connectivity'
    call generate_binary_structure_nd(2,cross3,status=s)
    if (s/=0 .or. count(cross3)/=19) error stop '3d structure'
    call generate_binary_structure_nd(1,wrong,status=s)
    if (s/=morphology_status_invalid_input) error stop 'wrong structure shape'
    call iterate_structure_nd(cross,2,diamond,status=s)
    if (s/=0 .or. count(diamond)/=13) error stop '2d diamond iteration'
    if (diamond(1,1) .or. .not.diamond(3,1)) error stop 'diamond coordinates'
    call iterate_structure_nd(cross3,2,sphere,status=s)
    if (s/=0 .or. count(sphere)/=93) error stop '3d diagonal iteration'
    edge=[.false.,.true.]
    call iterate_structure_nd(edge,3,expanded,status=s)
    if (s/=0 .or. count(expanded)/=1 .or. .not.expanded(4)) &
        error stop 'even asymmetric iteration'
    call iterate_structure_nd(cross,2,diamond,status=s,origin=[1,-1],new_origin=origins)
    if (s/=0 .or. any(origins/=[2,-2])) error stop 'origin scaling'
    call iterate_structure_nd(cross,2,wrong,status=s)
    if (s/=morphology_status_invalid_input) error stop 'wrong expanded shape'
    stencil15=.false.
    stencil15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=.true.
    stencil15(3,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=.true.
    call iterate_structure_nd(stencil15,2,extended15,status=s)
    if (s/=0 .or. count(extended15)/=3) error stop 'rank15 iteration'
    call iterate_structure_nd(stencil15,10001,stencil15,status=s)
    if (s/=morphology_status_invalid_input) error stop 'output shape validates early'
    do j=1,5
        do i=1,4
            ints(i,j)=0
            if ((i==2 .and. j==2) .or. (i==2 .and. j==3)) ints(i,j)=-2
            if (i==4 .and. j==5) ints(i,j)=5
            a(i,j)=sin(0.4_dp*i+0.3_dp*j)+real(i*j,dp)/10.0_dp
        end do
    end do
    call label_nd(ints,labels,n,status=s)
    if (s/=0 .or. n/=2) error stop 'integer labeling'
    reals=real(ints,dp)*0.25_dp
    call label_nd(reals,labels_copy,n,status=s)
    if (s/=0 .or. n/=2 .or. any(labels/=labels_copy)) error stop 'real labeling'
    call label_nd(ints,labels,n,status=s,connectivity=0)
    if (s/=morphology_status_invalid_input) error stop 'invalid integer connectivity'
    int15=0
    int15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=9
    call label_nd(int15,label15,n,status=s)
    if (s/=0 .or. n/=1) error stop 'rank15 integer labels'
    if (label15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)/=0) &
        error stop 'rank15 integer background'
    singleton=.true.
    call iterate_structure_nd(singleton,10001,same_singleton,status=s)
    if (s/=morphology_status_too_expensive) error stop 'iteration limit status'

    call generic_gradient_magnitude_nd(a,first_derivative,b,status=s,axes=[2])
    call correlate1d_nd(a,[-0.5_dp,0.0_dp,0.5_dp],2,term,status=n)
    if (s/=0 .or. n/=0 .or. maxval(abs(b-abs(term)))>1.e-14_dp) &
        error stop 'gradient axis 2'
    call generic_gradient_magnitude_nd(a,first_derivative,b,status=s,axes=[2,1])
    call correlate1d_nd(a,[-0.5_dp,0.0_dp,0.5_dp],1,term2,status=n)
    if (s/=0 .or. maxval(abs(b-sqrt(term**2+term2**2)))>1.e-14_dp) &
        error stop 'gradient permuted axes'
    call generic_laplace_nd(a,second_derivative,b,status=s,axes=[1], &
        modes=['nearest         ','constant        '],cval=-2.0_dp)
    call correlate1d_nd(a,[1.0_dp,-2.0_dp,1.0_dp],1,term,status=n,mode='nearest')
    if (s/=0 .or. n/=0 .or. maxval(abs(b-term))>1.e-14_dp) &
        error stop 'selected Laplacian mode'
    call generic_laplace_nd(a,second_derivative,b,status=s,axes=[1,1])
    if (s/=ndimage_status_invalid_input) error stop 'duplicate axes validation'
    call generic_laplace_nd(a,second_derivative,b,status=s,axes=[0])
    if (s/=ndimage_status_invalid_input) error stop 'zero-based axis validation'
    call generic_laplace_nd(a,second_derivative,b,status=s,axes=[3])
    if (s/=ndimage_status_invalid_input) error stop 'out-of-range axis validation'
    call generic_laplace_nd(a,second_derivative,b,status=s,axes=[integer::])
    if (s/=ndimage_status_invalid_input) error stop 'empty axes validation'
    x15=0.0_dp
    x15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=3.0_dp
    call generic_gradient_magnitude_nd(x15,first_derivative,y15,status=s,axes=[1])
    if (s/=0 .or. any(abs(y15-1.5_dp)>1.e-14_dp)) error stop 'rank15 selected axis'
    print '(a)', 'PASS: binary structures, numeric labeling, selected derivative axes'
contains
    subroutine first_derivative(input,axis,output,mode,cval,status)
        real(dp), intent(in) :: input(..)
        integer, intent(in) :: axis
        real(dp), intent(out) :: output(..)
        character(len=*), intent(in) :: mode
        real(dp), intent(in) :: cval
        integer, intent(out) :: status
        call correlate1d_nd(input,[-0.5_dp,0.0_dp,0.5_dp],axis,output, &
            status=status,mode=mode,cval=cval)
    end subroutine first_derivative
    subroutine second_derivative(input,axis,output,mode,cval,status)
        real(dp), intent(in) :: input(..)
        integer, intent(in) :: axis
        real(dp), intent(out) :: output(..)
        character(len=*), intent(in) :: mode
        real(dp), intent(in) :: cval
        integer, intent(out) :: status
        call correlate1d_nd(input,[1.0_dp,-2.0_dp,1.0_dp],axis,output, &
            status=status,mode=mode,cval=cval)
    end subroutine second_derivative
end program test_ndimage_structure_axes
