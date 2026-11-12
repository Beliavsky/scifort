! SPDX-License-Identifier: MIT
! Binary reference stream for brute-force transforms and derivative callbacks.
program bruteforce_generic_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : distance_transform_bf_nd, &
        generic_gradient_magnitude_nd, generic_laplace_nd, correlate1d_nd
    implicit none
    logical :: m1(8),m2(7,6),m3(4,5,3)
    logical :: m15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: d1(8),d2(7,6),d3(4,5,3)
    real(dp) :: d15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: a2(5,6),b2(5,6),a3(4,3,2),b3(4,3,2)
    integer :: idx1(1,8),idx2(2,42),idx3(3,60),idx15(15,2)
    character(len=16), parameter :: metrics(3)=[character(len=16) :: &
        'euclidean','taxicab','chessboard']
    character(len=16), parameter :: modes(8)=[character(len=16) :: &
        'reflect','grid-mirror','mirror','nearest', &
        'wrap','grid-wrap','constant','grid-constant']
    character(len=16) :: mode3(3)
    integer :: u,s,i,j,k,p

    open(newunit=u,file='bruteforce_generic_results.bin',status='replace', &
        action='write',form='unformatted',access='stream')
    do i=1,8
        m1(i)=modulo(7*i+2,9)>2
    end do
    do j=1,6
        do i=1,7
            m2(i,j)=modulo(3*i+5*j+i*j,11)>3
        end do
    end do
    do k=1,3
        do j=1,5
            do i=1,4
                m3(i,j,k)=modulo(13*i+5*j+3*k+i*k,13)>2
            end do
        end do
    end do
    do p=1,3
        call distance_transform_bf_nd(m1,d1,status=s,metric=metrics(p),feature_indices=idx1)
        call check(s,'one-dimensional brute force')
        write(u) d1,idx1
    end do
    do p=1,3
        call distance_transform_bf_nd(m2,d2,status=s,metric=metrics(p),feature_indices=idx2)
        call check(s,'two-dimensional brute force')
        write(u) d2,idx2
    end do
    do p=1,3
        call distance_transform_bf_nd(m3,d3,status=s,metric=metrics(p),feature_indices=idx3)
        call check(s,'three-dimensional brute force')
        write(u) d3,idx3
    end do
    call distance_transform_bf_nd(m2,d2,status=s,sampling=[1.2_dp,2.7_dp],feature_indices=idx2)
    call check(s,'two-dimensional sampled Euclidean')
    write(u) d2,idx2
    call distance_transform_bf_nd(m3,d3,status=s,sampling=[0.7_dp,1.8_dp,2.3_dp],feature_indices=idx3)
    call check(s,'three-dimensional sampled Euclidean')
    write(u) d3,idx3
    do p=1,3
        call distance_transform_bf_nd(spread(.true.,1,8),d1,status=s, &
            metric=metrics(p),feature_indices=idx1)
        call check(s,'all-foreground brute force')
        write(u) d1,idx1
    end do
    do p=1,3
        call distance_transform_bf_nd(spread(.false.,1,8),d1,status=s, &
            metric=metrics(p),feature_indices=idx1)
        call check(s,'all-background brute force')
        write(u) d1,idx1
    end do
    m15=.true.
    m15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=.false.
    do p=1,3
        call distance_transform_bf_nd(m15,d15,status=s,metric=metrics(p),feature_indices=idx15)
        call check(s,'rank15 brute force')
        write(u) d15,idx15
    end do
    do j=1,6
        do i=1,5
            a2(i,j)=sin(.27_dp*i+.41_dp*j)+real(i*j-j*j,dp)/13.0_dp
        end do
    end do
    do k=1,2
        do j=1,3
            do i=1,4
                a3(i,j,k)=sin(.2_dp*i+.5_dp*j+.7_dp*k)+real(i*j*k,dp)/11.0_dp
            end do
        end do
    end do
    do p=1,8
        call generic_gradient_magnitude_nd(a2,first_derivative,b2,status=s, &
            mode=modes(p),cval=-0.9_dp)
        call check(s,'two-dimensional gradient callback')
        write(u) b2
        call generic_laplace_nd(a2,second_derivative,b2,status=s, &
            mode=modes(p),cval=-0.9_dp)
        call check(s,'two-dimensional Laplace callback')
        write(u) b2
        mode3=[modes(p),modes(modulo(p,8)+1),modes(modulo(p+3,8)+1)]
        call generic_gradient_magnitude_nd(a3,first_derivative,b3,status=s, &
            modes=mode3,cval=0.75_dp)
        call check(s,'three-dimensional mixed gradient callback')
        write(u) b3
        call generic_laplace_nd(a3,second_derivative,b3,status=s, &
            modes=mode3,cval=0.75_dp)
        call check(s,'three-dimensional mixed Laplace callback')
        write(u) b3
    end do
    close(u)
    print '(a)', 'PASS: brute-force and generic derivative reference output'
contains
    subroutine check(status,name)
        integer, intent(in) :: status
        character(len=*), intent(in) :: name
        if (status/=0) then
            print '(3a,i0)', 'FAIL ',name,': status=',status
            error stop 1
        end if
    end subroutine check
    subroutine first_derivative(x,axis,y,mode,cval,status)
        real(dp), intent(in) :: x(..)
        integer, intent(in) :: axis
        real(dp), intent(out) :: y(..)
        character(len=*), intent(in) :: mode
        real(dp), intent(in) :: cval
        integer, intent(out) :: status
        call correlate1d_nd(x,[-0.5_dp,0.0_dp,0.5_dp],axis,y,status=status,mode=mode,cval=cval)
    end subroutine first_derivative
    subroutine second_derivative(x,axis,y,mode,cval,status)
        real(dp), intent(in) :: x(..)
        integer, intent(in) :: axis
        real(dp), intent(out) :: y(..)
        character(len=*), intent(in) :: mode
        real(dp), intent(in) :: cval
        integer, intent(out) :: status
        call correlate1d_nd(x,[1.0_dp,-2.0_dp,1.0_dp],axis,y,status=status,mode=mode,cval=cval)
    end subroutine second_derivative
end program bruteforce_generic_driver
