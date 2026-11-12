! SPDX-License-Identifier: MIT
! For randomized SciPy reference testing, one case per process.
program cluster_hierarchy_flatdata_driver
    use scifort_kinds, only : dp
    use scifort_spatial_distance, only : pdist
    use scifort_cluster_hierarchy, only : fcluster, fclusterdata, linkage, &
        inconsistent, cophenetic_correlation, cluster_success
    implicit none

    integer :: n, nfeatures, depth, i, status, link_status, corr_status, r_status
    real(dp) :: threshold, corr
    real(dp), allocatable :: x(:,:), z(:,:), y(:), r(:,:)
    integer, allocatable :: labels(:), manual(:)
    character(len=60) :: method, metric, criterion

    read(*,*) n, nfeatures
    read(*,'(a)') method
    read(*,'(a)') metric
    read(*,'(a)') criterion
    read(*,*) threshold, depth
    allocate(x(n,nfeatures))
    do i = 1,n
        read(*,*) x(i,:)
    end do

    labels = fclusterdata(x,threshold,criterion=trim(criterion), &
        method=trim(method),metric=trim(metric),depth=depth,status=status)
    z = linkage(x,method=trim(method),metric=trim(metric),status=link_status)
    corr_status = link_status
    r_status = link_status
    corr = 0.0_dp
    allocate(manual(0))
    if (link_status == cluster_success) then
        y = pdist(x,metric=trim(metric))
        corr = cophenetic_correlation(z,y,corr_status)
        r = inconsistent(z,depth=depth,status=r_status)
        if (r_status == cluster_success) then
            manual = fcluster(z,threshold,criterion=trim(criterion),r=r,depth=depth)
        end if
    end if
    print '(4(i0,1x),es25.17)',status,size(labels),corr_status,size(manual),corr
    if (size(labels)>0) print '(*(i0,1x))',labels
    if (size(manual)>0) print '(*(i0,1x))',manual
end program cluster_hierarchy_flatdata_driver
