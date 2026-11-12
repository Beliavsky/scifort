! SPDX-License-Identifier: MIT
program test_spatial_slerp
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_spatial, only : geometric_slerp, spatial_distance_success, &
        spatial_distance_invalid_input
    implicit none
    real(dp), parameter :: pi=acos(-1.0_dp)
    real(dp) :: a(3), b(3), fractions(7), reference(3), t
    real(dp), allocatable :: points(:,:), point(:), reversed(:,:), empty(:,:)
    integer :: code, i

    a = [1.0_dp,0.0_dp,0.0_dp]
    b = [0.0_dp,1.0_dp,0.0_dp]
    fractions = [0.0_dp,0.25_dp,0.5_dp,0.75_dp,1.0_dp,-1.0_dp,2.0_dp]
    points = geometric_slerp(a,b,fractions,status=code)
    if (code/=spatial_distance_success) error stop 'slerp expected success'
    if (any(shape(points)/=[7,3])) error stop 'slerp output shape'
    do i=1,size(fractions)
        t = fractions(i)
        reference = [cos(t*pi/2.0_dp),sin(t*pi/2.0_dp),0.0_dp]
        call check(points(i,:),reference,2.0e-15_dp)
        if (abs(norm2(points(i,:))-1.0_dp)>2.0e-15_dp) error stop 'slerp sphere norm'
    end do
    point = geometric_slerp(a,b,0.5_dp,status=code)
    if (code/=0 .or. size(point)/=3) error stop 'scalar slerp shape/status'
    call check(point,points(3,:),1.0e-15_dp)
    reversed = geometric_slerp(b,a,1.0_dp-fractions,status=code)
    if (code/=0) error stop 'reverse slerp status'
    do i=1,size(fractions)
        call check(points(i,:),reversed(i,:),2.0e-15_dp)
    end do
    points = geometric_slerp(a,a,fractions,status=code)
    if (code/=0) error stop 'coincident slerp status'
    do i=1,size(fractions)
        call check(points(i,:),a,0.0_dp)
    end do
    empty = geometric_slerp(a,b,[real(dp)::],status=code)
    if (code/=0 .or. any(shape(empty)/=[0,3])) error stop 'empty fraction shape'

    ! A tiny but unambiguous arc must still interpolate at high precision.
    point = geometric_slerp([1.0_dp,0.0_dp], &
        [sqrt(1.0_dp-1.0e-20_dp),1.0e-10_dp],0.5_dp,status=code)
    if (code/=0) error stop 'near coincident path rejected'
    if (abs(point(2)-5.0e-11_dp)>1.0e-14_dp) &
        error stop 'near coincident path inaccuracy'
    ! An almost-antipodal path is ambiguous at the default tolerance, but
    ! may be requested when the caller intentionally specifies a tighter one.
    point = geometric_slerp(a,[-sqrt(1.0_dp-1.0e-16_dp),1.0e-8_dp,0.0_dp], &
                            0.5_dp,status=code)
    call invalid_vector(point,code)
    point = geometric_slerp(a,[-sqrt(1.0_dp-1.0e-16_dp),1.0e-8_dp,0.0_dp], &
                            0.5_dp,tol=1.0e-10_dp,status=code)
    if (code/=0 .or. size(point)/=3) error stop 'tolerance override rejected'
    if (abs(norm2(point)-1.0_dp)>1.0e-12_dp) error stop 'near-antipodal unit norm'

    point = geometric_slerp(a,-a,0.5_dp,status=code)
    call invalid_vector(point,code)
    points = geometric_slerp(a,-a,fractions,status=code)
    if (code/=spatial_distance_invalid_input .or. any(shape(points)/=[0,0])) &
        error stop 'antipodal path accepted'
    point = geometric_slerp(a,[1.1_dp,0.0_dp,0.0_dp],0.5_dp,status=code)
    call invalid_vector(point,code)
    point = geometric_slerp(a,[0.0_dp,1.0_dp],0.5_dp,status=code)
    call invalid_vector(point,code)
    point = geometric_slerp(a,b,ieee_value(0.0_dp,ieee_quiet_nan),status=code)
    call invalid_vector(point,code)
    point = geometric_slerp(a,b,0.5_dp,tol=-1.0_dp,status=code)
    call invalid_vector(point,code)
    point = geometric_slerp([1.0_dp],[-1.0_dp],0.5_dp,status=code)
    call invalid_vector(point,code)
    points = geometric_slerp(a,b,[0.0_dp,ieee_value(0.0_dp,ieee_quiet_nan)],status=code)
    if (code/=1 .or. any(shape(points)/=[0,0])) error stop 'nonfinite fractions accepted'
    print '(a)', 'test_spatial_slerp PASS'

contains
    subroutine check(actual,expected,tol)
        real(dp), intent(in) :: actual(:) !! computed spherical point
        real(dp), intent(in) :: expected(:) !! reference unit-sphere point
        real(dp), intent(in) :: tol !! absolute coordinate tolerance
        if (size(actual)/=size(expected)) error stop 'slerp length mismatch'
        if (.not. all(abs(actual-expected)<=tol)) then
            print *,actual,expected
            error stop 'slerp disagreement'
        end if
    end subroutine check

    subroutine invalid_vector(v,status)
        real(dp), intent(in) :: v(:) !! returned empty vector on error
        integer, intent(in) :: status !! invalid-input status code
        if (size(v)/=0 .or. status/=spatial_distance_invalid_input) &
            error stop 'invalid scalar slerp accepted'
    end subroutine invalid_vector
end program test_spatial_slerp
