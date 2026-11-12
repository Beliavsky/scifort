! SPDX-License-Identifier: MIT
! Produce binary reference arrays for footprint, numeric label and axes comparisons.
program structure_axes_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : generate_binary_structure_nd, iterate_structure_nd, &
        label_nd, generic_gradient_magnitude_nd, generic_laplace_nd, correlate1d_nd
    implicit none
    logical :: g1(3),g2(3,3),g3(3,3,3),g4(3,3,3,3)
    logical :: two(2,2),iter2(5,5),iter3(7,7),iter_even(4,4)
    logical :: asymmetric(2,3,3),asymmetric2(3,5,5)
    logical :: span(3),span2(5),span3(7)
    integer :: ints(6,5),labels2(6,5),features,s,i,j,k,u
    real(dp) :: real3(4,3,2),a2(5,6),a3(4,3,2),b2(5,6),b3(4,3,2)
    integer :: labels3(4,3,2)

    open(newunit=u,file='structure_axes_results.bin',form='unformatted', &
        access='stream',status='replace',action='write')
    call generate_binary_structure_nd(1,g1,status=s)
    if (s/=0) error stop 'generate rank1'
    write(u) merge(1,0,g1)
    call generate_binary_structure_nd(1,g2,status=s)
    if (s/=0) error stop 'generate rank2 conn1'
    write(u) merge(1,0,g2)
    call generate_binary_structure_nd(2,g2,status=s)
    if (s/=0) error stop 'generate rank2 conn2'
    write(u) merge(1,0,g2)
    call generate_binary_structure_nd(2,g3,status=s)
    if (s/=0) error stop 'generate rank3 conn2'
    write(u) merge(1,0,g3)
    call generate_binary_structure_nd(3,g4,status=s)
    if (s/=0) error stop 'generate rank4 conn3'
    write(u) merge(1,0,g4)

    call generate_binary_structure_nd(1,g2,status=s)
    call iterate_structure_nd(g2,2,iter2,status=s)
    if (s/=0) error stop 'iterate rank2 twice'
    write(u) merge(1,0,iter2)
    call iterate_structure_nd(g2,3,iter3,status=s)
    if (s/=0) error stop 'iterate rank2 thrice'
    write(u) merge(1,0,iter3)
    span=[.true.,.false.,.true.]
    call iterate_structure_nd(span,2,span2,status=s)
    if (s/=0) error stop 'iterate 1d 2'
    write(u) merge(1,0,span2)
    call iterate_structure_nd(span,3,span3,status=s)
    if (s/=0) error stop 'iterate 1d 3'
    write(u) merge(1,0,span3)
    two=reshape([.true.,.false.,.true.,.true.],shape(two))
    call iterate_structure_nd(two,3,iter_even,status=s)
    if (s/=0) error stop 'iterate even rank2'
    write(u) merge(1,0,iter_even)
    do k=1,3
        do j=1,3
            do i=1,2
                asymmetric(i,j,k)=mod(2*i+3*j+5*k+i*j,7)>3
            end do
        end do
    end do
    call iterate_structure_nd(asymmetric,2,asymmetric2,status=s)
    if (s/=0) error stop 'iterate 3d asymmetric'
    write(u) merge(1,0,asymmetric2)

    do j=1,5
        do i=1,6
            ints(i,j)=mod(i*i+3*j+2*i*j,7)-4
        end do
    end do
    call label_nd(ints,labels2,features,status=s)
    if (s/=0) error stop 'label 2d integer'
    write(u) features,labels2
    call label_nd(ints,labels2,features,status=s,connectivity=2)
    if (s/=0) error stop 'label 2d integer full connectivity'
    write(u) features,labels2
    do k=1,2
        do j=1,3
            do i=1,4
                real3(i,j,k)=real(mod(2*i+3*j+5*k+i*j,5)-2,dp)/3.0_dp
            end do
        end do
    end do
    call label_nd(real3,labels3,features,status=s)
    if (s/=0) error stop 'label 3d real'
    write(u) features,labels3
    call label_nd(real3,labels3,features,status=s,connectivity=3)
    if (s/=0) error stop 'label 3d real full connectivity'
    write(u) features,labels3

    do j=1,6
        do i=1,5
            a2(i,j)=sin(.32_dp*i+.41_dp*j)+real(i*j,dp)/13.0_dp
        end do
    end do
    do k=1,2
        do j=1,3
            do i=1,4
                a3(i,j,k)=cos(.23_dp*i+.14_dp*j+.77_dp*k)+real(i*j*k,dp)/19.0_dp
            end do
        end do
    end do
    call generic_gradient_magnitude_nd(a2,first_derivative,b2,status=s,mode='nearest',axes=[1])
    if (s/=0) error stop 'axes 2d first'
    write(u) b2
    call generic_laplace_nd(a2,second_derivative,b2,status=s,mode='reflect',axes=[2])
    if (s/=0) error stop 'axes 2d second'
    write(u) b2
    call generic_gradient_magnitude_nd(a2,first_derivative,b2,status=s,mode='wrap',axes=[2,1])
    if (s/=0) error stop 'axes 2d both'
    write(u) b2
    call generic_gradient_magnitude_nd(a3,first_derivative,b3,status=s,axes=[1,3], &
        modes=['nearest         ','constant        '],cval=-0.7_dp)
    if (s/=0) error stop 'axes 3d gradient compact'
    write(u) b3
    call generic_laplace_nd(a3,second_derivative,b3,status=s,axes=[3], &
        modes=['grid-wrap       '],cval=0.25_dp)
    if (s/=0) error stop 'axes 3d laplace compact'
    write(u) b3
    call generic_laplace_nd(a3,second_derivative,b3,status=s,axes=[2,1], &
        modes=['wrap            ','constant        ','nearest         '],cval=1.3_dp)
    if (s/=0) error stop 'axes 3d laplace full'
    write(u) b3
    call generic_gradient_magnitude_nd(a3,first_derivative,b3,status=s,axes=[2,3,1], &
        modes=['reflect         ','grid-wrap       ','nearest         '],cval=0.5_dp)
    if (s/=0) error stop 'axes 3d gradient permuted'
    write(u) b3
    close(u)
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
end program structure_axes_driver
