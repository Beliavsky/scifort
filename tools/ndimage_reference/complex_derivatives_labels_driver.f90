! SPDX-License-Identifier: MIT
! Binary reference output for complex derivatives and indexed label reductions.
program complex_derivatives_labels_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : generic_gradient_magnitude_nd, generic_laplace_nd, &
        correlate1d_nd, sum_labels_nd, mean_labels_nd, variance_labels_nd, &
        standard_deviation_labels_nd, count_labels_nd, histogram_labels_nd, &
        center_of_mass_nd, extrema_labels_nd, bounding_boxes_nd
    implicit none
    complex(dp) :: a2(6,5),b2(6,5),a3(4,3,2),b3(4,3,2)
    real(dp) :: x(5,4),sums(8),means(8),vars(8),stds(8),centers(2,8),mins(8),maxs(8)
    integer :: labels(5,4),req(8),counts(8),hist(7,8),mpos(2,8),xpos(2,8),lower(2,8),upper(2,8)
    integer :: i,j,k,p,s,u
    character(len=16), parameter :: modes(8)=[character(len=16) :: &
        'reflect','grid-mirror','mirror','nearest','wrap','grid-wrap','constant','grid-constant']
    character(len=16) :: selected_modes(2)
    complex(dp), parameter :: fill=cmplx(0.37_dp,-0.6_dp,dp)

    do j=1,5
        do i=1,6
            a2(i,j)=cmplx(sin(.25_dp*i+.37_dp*j),cos(.43_dp*i-.16_dp*j),dp)
        end do
    end do
    do k=1,2
        do j=1,3
            do i=1,4
                a3(i,j,k)=cmplx(sin(.29_dp*i+.41_dp*j+.14_dp*k), &
                    cos(.11_dp*i-.17_dp*j+.31_dp*k),dp)
            end do
        end do
    end do
    open(newunit=u,file='complex_derivatives_labels.bin',status='replace', &
        action='write',form='unformatted',access='stream')
    do p=1,8
        call generic_gradient_magnitude_nd(a2,first_derivative,b2,s,mode=trim(modes(p)),cval=fill)
        call check(s)
        write(u) b2
        call generic_laplace_nd(a2,second_derivative,b2,s,mode=trim(modes(p)),cval=fill)
        call check(s)
        write(u) b2
    end do
    do p=1,8
        selected_modes=[modes(p),modes(modulo(p+2,8)+1)]
        call generic_gradient_magnitude_nd(a3,first_derivative,b3,s, &
            modes=selected_modes,axes=[3,1],cval=fill)
        call check(s)
        write(u) b3
        call generic_laplace_nd(a3,second_derivative,b3,s, &
            modes=selected_modes,axes=[3,1],cval=fill)
        call check(s)
        write(u) b3
    end do
    do j=1,4
        do i=1,5
            x(i,j)=real(3*i-j,dp)/5.0_dp+cos(real(i*j,dp))
            select case(modulo(2*i+3*j,6))
            case(0)
                labels(i,j)=-7
            case(1)
                labels(i,j)=0
            case(2)
                labels(i,j)=500000001
            case default
                labels(i,j)=3
            end select
        end do
    end do
    req=[3,-7,0,3,500000001,-999,-7,huge(1)]
    call sum_labels_nd(x,labels,req,sums,s)
    call check(s)
    write(u) sums
    call mean_labels_nd(x,labels,req,means,s)
    call check(s)
    write(u) means
    call variance_labels_nd(x,labels,req,vars,s)
    call check(s)
    write(u) vars
    call standard_deviation_labels_nd(x,labels,req,stds,s)
    call check(s)
    write(u) stds
    call center_of_mass_nd(x,labels,req,centers,s)
    call check(s)
    write(u) centers
    call count_labels_nd(labels,req,counts,s)
    call check(s)
    write(u) counts
    call histogram_labels_nd(x,labels,req,-5.0_dp,5.0_dp,hist,s)
    call check(s)
    write(u) hist
    call extrema_labels_nd(x,labels,req,mins,maxs,mpos,xpos,s)
    call check(s)
    write(u) mins,maxs,mpos,xpos
    call bounding_boxes_nd(labels,req,lower,upper,s)
    call check(s)
    write(u) lower,upper
    close(u)
    print '(a)', 'PASS: complex and labeled region reference driver'
contains
    subroutine check(status)
        integer, intent(in) :: status !! operation status
        if (status/=0) error stop 'failed complex derivative or label reference operation'
    end subroutine check
    subroutine first_derivative(input,axis,output,mode,cval,status)
        complex(dp), intent(in) :: input(..) !! sampled image
        integer, intent(in) :: axis !! axis index
        complex(dp), intent(out) :: output(..) !! derivative
        character(len=*), intent(in) :: mode !! extension
        complex(dp), intent(in) :: cval !! extension value
        integer, intent(out) :: status !! return status
        call correlate1d_nd(input,cmplx([-0.5_dp,0.0_dp,0.5_dp],0.0_dp,dp), &
            axis,output,status=status,mode=mode,cval=cval)
    end subroutine first_derivative
    subroutine second_derivative(input,axis,output,mode,cval,status)
        complex(dp), intent(in) :: input(..) !! sampled image
        integer, intent(in) :: axis !! axis index
        complex(dp), intent(out) :: output(..) !! derivative
        character(len=*), intent(in) :: mode !! extension
        complex(dp), intent(in) :: cval !! extension value
        integer, intent(out) :: status !! return status
        call correlate1d_nd(input,cmplx([1.0_dp,-2.0_dp,1.0_dp],0.0_dp,dp), &
            axis,output,status=status,mode=mode,cval=cval)
    end subroutine second_derivative
end program complex_derivatives_labels_driver
