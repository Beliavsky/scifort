! SPDX-License-Identifier: MIT
program test_cluster_hierarchy_shortcuts
    use scifort_kinds, only : dp
    use scifort_cluster_hierarchy, only : linkage, single, complete, average, &
        weighted, centroid, median, ward, cluster_success
    implicit none
    real(dp) :: x(5,2), y(10)
    real(dp), allocatable :: z(:,:), reference(:,:)
    integer :: status, i

    x(:,1) = [0.0_dp, 1.0_dp, 1.5_dp, 4.0_dp, 6.0_dp]
    x(:,2) = [0.0_dp, 2.0_dp, 1.0_dp, 5.0_dp, 3.0_dp]
    i = 0
    block
        integer :: a, b
        do a = 1, 4
            do b = a+1, 5
                i = i+1
                y(i) = sqrt(sum((x(a,:)-x(b,:))**2))
            end do
        end do
    end block

    call compare('single',single(y,status),single(x,status=status),y,x)
    call compare('complete',complete(y,status),complete(x,status=status),y,x)
    call compare('average',average(y,status),average(x,status=status),y,x)
    call compare('weighted',weighted(y,status),weighted(x,status=status),y,x)
    call compare('centroid',centroid(y,status),centroid(x,status=status),y,x)
    call compare('median',median(y,status),median(x,status=status),y,x)
    call compare('ward',ward(y,status),ward(x,status=status),y,x)

    z = single(x,metric='not-a-metric',status=status)
    if (status == cluster_success .or. size(z,1) /= 0) error stop 'invalid metric accepted'
    print *, 'hierarchy named-linkage shortcuts PASS'
contains
    subroutine compare(rule, actual_condensed, actual_obs, condensed, observations)
        character(len=*), intent(in) :: rule !! Linkage name.
        real(dp), intent(in) :: actual_condensed(:,:) !! Named condensed result.
        real(dp), intent(in) :: actual_obs(:,:) !! Named observation result.
        real(dp), intent(in) :: condensed(:) !! Input condensed distances.
        real(dp), intent(in) :: observations(:,:) !! Input observations.
        real(dp), allocatable :: expected(:,:)

        expected = linkage(condensed,method=rule,status=status)
        if (status /= cluster_success) error stop 'condensed method error'
        if (any(shape(actual_condensed) /= shape(expected))) error stop 'condensed shortcut shape'
        if (any(abs(actual_condensed-expected) > 1.0e-12_dp)) error stop 'condensed shortcut'
        expected = linkage(observations,method=rule,status=status)
        if (status /= cluster_success) error stop 'observation method error'
        if (any(shape(actual_obs) /= shape(expected))) error stop 'observation shortcut shape'
        if (any(abs(actual_obs-expected) > 1.0e-12_dp)) error stop 'observation shortcut'
    end subroutine compare
end program test_cluster_hierarchy_shortcuts
