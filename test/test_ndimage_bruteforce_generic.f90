! SPDX-License-Identifier: MIT
! Native brute-force distances and generic derivative callback regression.
program test_ndimage_bruteforce_generic
    use scifort_kinds, only : dp
    use scifort_ndimage, only : distance_transform_bf_nd, &
        generic_gradient_magnitude_nd, generic_laplace_nd, &
        correlate1d_nd, laplace_nd, ndimage_status_invalid_input, &
        ndimage_status_too_expensive
    implicit none
    logical :: mask(3,4),large_mask(300,300)
    logical :: mask15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: d(3,4),a(3,4),b(3,4),c(3,4),t1(3,4),t2(3,4)
    real(dp) :: d15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: x15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: y15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    real(dp) :: large_d(300,300)
    integer :: indices(2,12),features15(15,2),s,i,j

    mask=.true.
    mask(2,3)=.false.
    call distance_transform_bf_nd(mask,d,status=s,feature_indices=indices)
    if (s/=0 .or. abs(d(1,1)-sqrt(5.0_dp))>1.e-14_dp) &
        error stop 'Euclidean brute-force distance'
    if (any(indices(:,1)/=[2,3])) error stop 'one-based brute-force witness'
    call distance_transform_bf_nd(mask,d,status=s,sampling=[2.0_dp,3.0_dp])
    if (s/=0 .or. abs(d(1,1)-sqrt(40.0_dp))>1.e-14_dp) &
        error stop 'anisotropic brute-force Euclidean distance'
    call distance_transform_bf_nd(mask,d,status=s,metric='cityblock')
    if (s/=0 .or. d(1,1)/=3.0_dp) error stop 'taxicab alias distance'
    call distance_transform_bf_nd(mask,d,status=s,metric='chessboard')
    if (s/=0 .or. d(1,1)/=2.0_dp) error stop 'chessboard distance'
    call distance_transform_bf_nd(mask,d,status=s,metric='invalid')
    if (s/=ndimage_status_invalid_input) error stop 'invalid brute-force metric'
    call distance_transform_bf_nd(mask,d,status=s,sampling=[1.0_dp])
    if (s/=ndimage_status_invalid_input) error stop 'invalid sampling shape'
    call distance_transform_bf_nd(mask,d,status=s,sampling=[-1.0_dp,1.0_dp])
    if (s/=ndimage_status_invalid_input) error stop 'invalid Euclidean sampling'
    call distance_transform_bf_nd(mask,d,status=s,feature_indices=features15)
    if (s/=ndimage_status_invalid_input) error stop 'invalid feature shape'
    mask=.false.
    call distance_transform_bf_nd(mask,d,status=s,feature_indices=indices)
    if (s/=0 .or. any(d/=0.0_dp)) error stop 'all-background distance'
    mask=.true.
    call distance_transform_bf_nd(mask,d,status=s,feature_indices=indices)
    if (s/=0 .or. any(d/=sqrt(huge(1.0_dp))) .or. any(indices/=0)) &
        error stop 'all-foreground brute-force sentinel'
    call distance_transform_bf_nd(mask,d,status=s,metric='taxicab')
    if (s/=0 .or. any(d/=4294967295.0_dp)) error stop 'no-site taxicab sentinel'
    large_mask=.true.
    large_mask(:,1:150)=.false.
    call distance_transform_bf_nd(large_mask,large_d,status=s)
    if (s/=ndimage_status_too_expensive) error stop 'brute-force pair guard'
    mask15=.true.
    mask15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=.false.
    call distance_transform_bf_nd(mask15,d15,status=s,feature_indices=features15)
    if (s/=0) error stop 'brute-force rank15'
    if (abs(d15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)-1.0_dp)>1.e-14_dp) &
        error stop 'rank15 background distance'
    if (features15(1,1)/=2 .or. any(features15(2:,1)/=1)) &
        error stop 'rank15 background coordinate'

    do j=1,4
        do i=1,3
            a(i,j)=sin(0.3_dp*i+0.7_dp*j)+real(i*j,dp)/11.0_dp
        end do
    end do
    call generic_gradient_magnitude_nd(a,first_derivative,b,status=s)
    if (s/=0) error stop 'generic gradient callback'
    call correlate1d_nd(a,[-0.5_dp,0.0_dp,0.5_dp],1,t1,status=s)
    if (s/=0) error stop 'axis derivative reference 1'
    call correlate1d_nd(a,[-0.5_dp,0.0_dp,0.5_dp],2,t2,status=s)
    if (s/=0) error stop 'axis derivative reference 2'
    c=sqrt(t1**2+t2**2)
    if (maxval(abs(b-c))>1.e-14_dp) error stop 'generic gradient value'
    call generic_laplace_nd(a,second_derivative,b,status=s)
    if (s/=0) error stop 'generic Laplacian callback'
    call laplace_nd(a,c,status=s)
    if (s/=0 .or. maxval(abs(b-c))>1.e-14_dp) error stop 'generic Laplacian value'
    call generic_laplace_nd(a,second_derivative,b,status=s, &
        modes=['constant        ','nearest         '],cval=-2.0_dp)
    if (s/=0) error stop 'generic mixed derivative'
    call correlate1d_nd(a,[1.0_dp,-2.0_dp,1.0_dp],1,t1,status=s,mode='constant',cval=-2.0_dp)
    if (s/=0) error stop 'mixed derivative reference 1'
    call correlate1d_nd(a,[1.0_dp,-2.0_dp,1.0_dp],2,t2,status=s,mode='nearest')
    if (s/=0 .or. maxval(abs(b-(t1+t2)))>1.e-14_dp) &
        error stop 'generic mixed derivative value'
    call generic_laplace_nd(a,second_derivative,b,status=s,mode='nearest', &
        modes=['nearest         ','reflect         '])
    if (s/=ndimage_status_invalid_input) error stop 'conflicting callback modes'
    call generic_gradient_magnitude_nd(a,first_derivative,b,status=s, &
        modes=['nearest         '])
    if (s/=ndimage_status_invalid_input) error stop 'invalid callback mode count'
    call generic_laplace_nd(a,failing_derivative,b,status=s)
    if (s/=ndimage_status_invalid_input) error stop 'callback error propagation'
    x15=2.0_dp
    x15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)=5.0_dp
    call generic_gradient_magnitude_nd(x15,first_derivative,y15,status=s)
    if (s/=0 .or. any(abs(y15-1.5_dp)>1.e-14_dp)) &
        error stop 'rank15 generic gradient'
    call generic_laplace_nd(x15,second_derivative,y15,status=s)
    if (s/=0) error stop 'rank15 generic Laplacian status'
    if (abs(y15(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)-3.0_dp)>1.e-14_dp) &
        error stop 'rank15 generic Laplacian first pixel'
    if (abs(y15(2,1,1,1,1,1,1,1,1,1,1,1,1,1,1)+3.0_dp)>1.e-14_dp) &
        error stop 'rank15 generic Laplacian second pixel'
    print '(a)', 'PASS: brute-force distances and generic gradient/Laplacian callbacks'
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
    subroutine failing_derivative(input,axis,output,mode,cval,status)
        real(dp), intent(in) :: input(..)
        integer, intent(in) :: axis
        real(dp), intent(out) :: output(..)
        character(len=*), intent(in) :: mode
        real(dp), intent(in) :: cval
        integer, intent(out) :: status
        status=ndimage_status_invalid_input
    end subroutine failing_derivative
end program test_ndimage_bruteforce_generic
