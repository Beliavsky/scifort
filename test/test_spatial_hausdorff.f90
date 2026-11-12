! SPDX-License-Identifier: MIT
program test_spatial_hausdorff
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, &
        ieee_positive_inf, ieee_value
    use scifort_kinds, only : dp
    use scifort_spatial, only : directed_hausdorff, hausdorff_distance, &
        spatial_distance_success, spatial_distance_invalid_input
    implicit none
    real(dp) :: u(2,2), v(2,2), x(2,2), z(2,2), d, r
    integer :: i, j, code
    logical :: from_u

    u(1,:) = [0.0_dp,0.0_dp]
    u(2,:) = [3.0_dp,0.0_dp]
    v(1,:) = [0.0_dp,0.0_dp]
    v(2,:) = [1.0_dp,0.0_dp]
    d = directed_hausdorff(u,v,index_u=i,index_v=j,status=code)
    call check('directed',d,2.0_dp,1.0e-14_dp)
    if (code/=spatial_distance_success .or. i/=2 .or. j/=2) &
        error stop 'directed Hausdorff witness'
    r = directed_hausdorff(v,u)
    call check('reverse',r,1.0_dp,1.0e-14_dp)
    d = hausdorff_distance(u,v,index_u=i,index_v=j,from_u=from_u,status=code)
    call check('symmetric',d,2.0_dp,1.0e-14_dp)
    if (code/=0 .or. .not. from_u .or. i/=2 .or. j/=2) &
        error stop 'symmetric Hausdorff witness'
    d = hausdorff_distance(v,u,index_u=i,index_v=j,from_u=from_u,status=code)
    call check('symmetric reversal',d,2.0_dp,1.0e-14_dp)
    if (code/=0 .or. from_u .or. i/=2 .or. j/=2) &
        error stop 'symmetric reversed witness'

    ! The two implementations must agree even when early abandonment occurs.
    x(1,:) = [0.0_dp,0.0_dp]
    x(2,:) = [0.0_dp,0.0_dp]
    z(1,:) = [1.0_dp,0.0_dp]
    z(2,:) = [1.0_dp,0.0_dp]
    d = directed_hausdorff(x,z,index_u=i,index_v=j)
    r = directed_hausdorff(x,z,early_break=.false.)
    call check('ties',d,r,0.0_dp)
    if (i/=1 .or. j/=1) error stop 'tie is not first witness'
    d = directed_hausdorff(u,u,index_u=i,index_v=j)
    call check('identical clouds',d,0.0_dp,0.0_dp)
    if (i/=1 .or. j/=1) error stop 'identical witness'

    d = directed_hausdorff(reshape([1.0e200_dp,0.0_dp],[1,2]), &
                           reshape([0.0_dp,0.0_dp],[1,2]))
    call check('large magnitude',d/1.0e200_dp,1.0_dp,1.0e-14_dp)
    d = directed_hausdorff(reshape([1.0e-200_dp,1.0e-200_dp],[1,2]), &
                           reshape([0.0_dp,0.0_dp],[1,2]))
    call check('small magnitude',d/1.0e-200_dp,sqrt(2.0_dp),1.0e-14_dp)
    d = directed_hausdorff(reshape([0.0_dp,0.0_dp],[1,2]), &
                           reshape([0.0_dp,0.0_dp],[1,2]))
    call check('one-point',d,0.0_dp,0.0_dp)

    d = directed_hausdorff(u(1:0,:),v,index_u=i,index_v=j,status=code)
    call invalid(d,i,j,code)
    d = directed_hausdorff(u,v(1:0,:),index_u=i,index_v=j,status=code)
    call invalid(d,i,j,code)
    d = directed_hausdorff(u(:,1:0),v(:,1:0),index_u=i,index_v=j,status=code)
    call invalid(d,i,j,code)
    d = directed_hausdorff(u,v(:,1:1),index_u=i,index_v=j,status=code)
    call invalid(d,i,j,code)
    x = u
    x(2,2) = ieee_value(0.0_dp,ieee_quiet_nan)
    d = directed_hausdorff(x,v,index_u=i,index_v=j,status=code)
    call invalid(d,i,j,code)
    x(2,2) = ieee_value(0.0_dp,ieee_positive_inf)
    d = directed_hausdorff(x,v,index_u=i,index_v=j,status=code)
    call invalid(d,i,j,code)
    d = hausdorff_distance(x,v,index_u=i,index_v=j,from_u=from_u,status=code)
    call invalid(d,i,j,code)
    if (from_u) error stop 'invalid symmetric Hausdorff direction' 
    print '(a)', 'test_spatial_hausdorff PASS'

contains
    subroutine check(label,actual,expected,tol)
        character(len=*), intent(in) :: label !! diagnostic description
        real(dp), intent(in) :: actual !! computed distance
        real(dp), intent(in) :: expected !! independently expected distance
        real(dp), intent(in) :: tol !! allowed absolute error
        if (.not. abs(actual-expected)<=tol) then
            print *,label,actual,expected
            error stop 'Hausdorff disagreement'
        end if
    end subroutine check

    subroutine invalid(actual,iu,iv,status)
        real(dp), intent(in) :: actual !! invalid distance must be NaN
        integer, intent(in) :: iu !! invalid source witness must be zero
        integer, intent(in) :: iv !! invalid target witness must be zero
        integer, intent(in) :: status !! invalid status code
        if (status/=spatial_distance_invalid_input .or. .not.ieee_is_nan(actual) &
            .or. iu/=0 .or. iv/=0) error stop 'invalid Hausdorff input accepted'
    end subroutine invalid
end program test_spatial_hausdorff
