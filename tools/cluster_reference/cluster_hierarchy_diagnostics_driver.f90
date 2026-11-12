! SPDX-License-Identifier: MIT
! Independent SciPy cross-check driver for extended hierarchy diagnostics.
program cluster_hierarchy_diagnostics_driver
    use scifort_kinds, only : dp
    use scifort_cluster_hierarchy, only : cluster_success, inconsistent, maxdists, &
        maxrstat, maxinconsts, fcluster_inconsistent, fcluster_monocrit, &
        fcluster_maxclust_monocrit, cut_tree, is_monotonic
    implicit none
    real(dp), allocatable :: z(:,:), r(:,:), d(:), stat(:), hc(:)
    integer, allocatable :: labels(:), cuts(:,:)
    real(dp) :: incon_cut, mono_cut
    integer :: n, depth, k, status, i, j
    logical :: monotone

    read(*,*) n, depth, incon_cut, mono_cut, k
    if (n < 2) error stop 'n < 2'
    allocate(z(n-1,4),hc(4))
    do i = 1,n-1
        read(*,*) z(i,:)
    end do
    read(*,*) hc
    r = inconsistent(z,depth,status)
    if (status /= cluster_success) error stop 'inconsistent'
    do i = 1,n-1
        write(*,'(*(es26.17e3,1x))') r(i,:)
    end do
    d = maxdists(z,status)
    if (status /= cluster_success) error stop 'maxdists'
    write(*,'(*(es26.17e3,1x))') d
    do j = 1,4
        stat = maxrstat(z,r,j,status)
        if (status /= cluster_success) error stop 'maxRstat'
        write(*,'(*(es26.17e3,1x))') stat
    end do
    stat = maxinconsts(z,r,status)
    if (status /= cluster_success) error stop 'maxinconsts'
    write(*,'(*(es26.17e3,1x))') stat
    labels = fcluster_inconsistent(z,incon_cut,status,depth)
    if (status /= cluster_success) error stop 'fcluster_inconsistent'
    write(*,'(*(i0,1x))') labels
    labels = fcluster_inconsistent(z,incon_cut,status,r=r)
    if (status /= cluster_success) error stop 'fcluster_inconsistent_R'
    write(*,'(*(i0,1x))') labels
    labels = fcluster_monocrit(z,d,mono_cut,status)
    if (status /= cluster_success) error stop 'fcluster_monocrit'
    write(*,'(*(i0,1x))') labels
    labels = fcluster_maxclust_monocrit(z,d,k,status)
    if (status /= cluster_success) error stop 'fcluster_maxclust_monocrit'
    write(*,'(*(i0,1x))') labels
    cuts = cut_tree(z,status=status)
    if (status /= cluster_success) error stop 'cut_tree_full'
    do i = 1,n
        write(*,'(*(i0,1x))') cuts(i,:)
    end do
    cuts = cut_tree(z,n_clusters=[1,k,n,k],status=status)
    if (status /= cluster_success) error stop 'cut_tree_counts'
    do i = 1,n
        write(*,'(*(i0,1x))') cuts(i,:)
    end do
    monotone = is_monotonic(z,status)
    if (status /= cluster_success) error stop 'is_monotonic'
    if (.not. monotone) then
        write(*,'(a)') 'NONMONOTONE'
    else
        cuts = cut_tree(z,height=hc,status=status)
        if (status /= cluster_success) error stop 'cut_tree_height'
        do i = 1,n
            write(*,'(*(i0,1x))') cuts(i,:)
        end do
    end if
end program cluster_hierarchy_diagnostics_driver
