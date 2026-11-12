! SPDX-License-Identifier: MIT
! Complex line callbacks and axis-selective binary morphology regression.
program test_ndimage_complex_line_binary_axes
    use scifort_kinds, only : dp
    use scifort_ndimage, only : generic_filter1d_nd, binary_dilation_nd, &
        binary_erosion_nd, binary_opening_nd, binary_closing_nd, &
        ndimage_status_invalid_input, ndimage_status_too_expensive
    implicit none
    complex(dp) :: c(4,3),z(4,3),z_copy(4,3),cval
    complex(dp) :: strided(8,3)
    complex(dp) :: c15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    complex(dp) :: z15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    logical :: a(4,3,2),b(4,3,2),expected(4,3,2),one(4)
    logical :: mask(4,3,2),s3(3),s23(2,3)
    logical :: binary15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    logical :: result15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    logical :: wider(8,3,2)
    integer :: i,j,k,s

    do j=1,3
        do i=1,4
            c(i,j)=cmplx(0.1_dp*i*j,real(i-2*j,dp),dp)
        end do
    end do
    call generic_filter1d_nd(c,identity_line,1,2,z,status=s)
    if (s/=0 .or. any(z/=c)) error stop 'complex line identity'
    strided=(0.0_dp,0.0_dp)
    call generic_filter1d_nd(c,identity_line,1,2,strided(::2,:),status=s)
    if (s/=0 .or. any(strided(::2,:)/=c)) error stop 'complex strided output'
    cval=cmplx(2.0_dp,-3.0_dp,dp)
    call generic_filter1d_nd(c,first_sample,2,1,z,status=s, &
        mode='constant',cval=cval)
    if (s/=0) error stop 'complex constant line status'
    if (any(z(1,:)/=cval) .or. any(z(2:,:)/=c(:3,:))) &
        error stop 'complex line constant exterior'
    call generic_filter1d_nd(c,complex_weighted,4,2,z,status=s, &
        mode='grid-wrap',origin=-1)
    if (s/=0) error stop 'complex weighted line status'
    z_copy=z
    call generic_filter1d_nd(c,complex_weighted,4,2,z,status=s, &
        mode='grid-wrap',origin=-1)
    if (s/=0 .or. any(z/=z_copy)) error stop 'complex line deterministic result'
    call generic_filter1d_nd(c,identity_line,1,3,z,status=s)
    if (s/=ndimage_status_invalid_input) error stop 'complex invalid axis'
    call generic_filter1d_nd(c,identity_line,0,1,z,status=s)
    if (s/=ndimage_status_invalid_input) error stop 'complex invalid width'
    call generic_filter1d_nd(c,identity_line,2,1,z,status=s,origin=1)
    if (s/=ndimage_status_invalid_input) error stop 'complex invalid origin'
    call generic_filter1d_nd(c,identity_line,2,1,z,status=s,mode='bad')
    if (s/=ndimage_status_invalid_input) error stop 'complex invalid boundary'
    call generic_filter1d_nd(c,identity_line,100001,1,z,status=s)
    if (s/=ndimage_status_too_expensive) error stop 'complex size guard'
    call generic_filter1d_nd(c,identity_line,1,1,z(:2,:),status=s)
    if (s/=ndimage_status_invalid_input) error stop 'complex output shape'
    call generic_filter1d_nd(c,identity_line,1,1,z(::2,:),status=s)
    if (s/=ndimage_status_invalid_input) error stop 'complex output strided shape'
    c15=cmplx(3.0_dp,-5.0_dp,dp)
    call generic_filter1d_nd(c15,identity_line,1,15,z15,status=s)
    if (s/=0 .or. any(z15/=c15)) error stop 'rank15 complex callback'

    do k=1,2
        do j=1,3
            do i=1,4
                a(i,j,k)=mod(i+3*j+5*k,7)<3
                mask(i,j,k)=mod(2*i+j+k,5)/=0
            end do
        end do
    end do
    s3=[.true.,.true.,.true.]
    s23=reshape([.true.,.false.,.true.,.true.,.true.,.false.],shape(s23))
    call binary_dilation_nd(a,b,status=s,axes=[integer::])
    if (s/=0 .or. any(b.neqv.a)) error stop 'empty axes identity'
    call binary_erosion_nd(a,b,status=s,axes=[integer::])
    if (s/=0 .or. any(b.neqv.a)) error stop 'empty axes erosion identity'
    call binary_dilation_nd(a,b,status=s,axes=[1],structure=s3)
    if (s/=0) error stop 'axis1 morphology status'
    do k=1,2
        do j=1,3
            call binary_dilation_nd(a(:,j,k),one,status=s,structure=s3)
            if (s/=0) error stop '1d binary slice status'
            expected(:,j,k)=one
        end do
    end do
    if (any(expected.neqv.b)) error stop 'binary selected-axis slices'
    wider=.false.
    call binary_dilation_nd(a,wider(::2,:,:),status=s,axes=[1],structure=s3)
    if (s/=0 .or. any(wider(::2,:,:).neqv.expected)) &
        error stop 'selected-axis strided binary output'
    binary15=.false.
    binary15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=.true.
    call binary_dilation_nd(binary15,result15,status=s,axes=[1])
    if (s/=0 .or. count(result15)/=2) error stop 'rank15 selected binary dilation'
    call binary_closing_nd(binary15,result15,status=s,axes=[integer::])
    if (s/=0 .or. any(result15.neqv.binary15)) error stop 'rank15 empty axes'
    call binary_opening_nd(a,b,status=s,axes=[3,1],structure=s23,origin=[1,-1],mask=mask)
    if (s/=0) error stop 'masked selected-axis opening'
    call binary_closing_nd(a,b,status=s,axes=[1,3],structure=s23,origin=[-1,1],mask=mask)
    if (s/=0) error stop 'masked selected-axis closing'
    call binary_dilation_nd(a,b,status=s,axes=[1,1])
    if (s/=ndimage_status_invalid_input) error stop 'duplicate axes'
    call binary_dilation_nd(a,b,status=s,axes=[0])
    if (s/=ndimage_status_invalid_input) error stop 'zero axis'
    call binary_dilation_nd(a,b,status=s,axes=[4])
    if (s/=ndimage_status_invalid_input) error stop 'axis too large'
    call binary_dilation_nd(a,b,status=s,axes=[1,2],structure=s3)
    if (s/=ndimage_status_invalid_input) error stop 'footprint rank vs axes'
    call binary_erosion_nd(a,b,status=s,axes=[1,2],origin=[0])
    if (s/=ndimage_status_invalid_input) error stop 'origin shape vs axes'
    call binary_opening_nd(a,b,status=s,axes=[1,2],connectivity=3)
    if (s/=ndimage_status_invalid_input) error stop 'selected connectivity bounds'
    call binary_closing_nd(a,b,status=s,axes=[integer::],structure=s3)
    if (s/=ndimage_status_invalid_input) error stop 'empty axes footprint'

    print '(a)', 'PASS: complex generic_filter1d and selected-axis binary morphology'
contains
    subroutine identity_line(input_line,output_line)
        complex(dp), intent(in) :: input_line(:)
        complex(dp), intent(out) :: output_line(:)
        output_line=input_line(1:size(output_line))
    end subroutine identity_line

    subroutine first_sample(input_line,output_line)
        complex(dp), intent(in) :: input_line(:)
        complex(dp), intent(out) :: output_line(:)
        output_line=input_line(1:size(output_line))
    end subroutine first_sample

    subroutine complex_weighted(input_line,output_line)
        complex(dp), intent(in) :: input_line(:)
        complex(dp), intent(out) :: output_line(:)
        integer :: j,t,width
        width=size(input_line)-size(output_line)+1
        output_line=(0.0_dp,0.0_dp)
        do j=1,size(output_line)
            do t=1,width
                output_line(j)=output_line(j)+cmplx(0.2_dp*t,0.03_dp*(-1)**t*t,dp)*input_line(j+t-1)
            end do
        end do
    end subroutine complex_weighted
end program test_ndimage_complex_line_binary_axes
