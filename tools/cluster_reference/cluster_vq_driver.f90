! SPDX-License-Identifier: MIT
! Standalone protocol for independent SciPy numerical comparisons.
program cluster_vq_driver
    use scifort_kinds, only : dp
    use scifort_cluster_vq, only : whiten, vq, kmeans2, kmeans
    implicit none

    character(len=16) :: operation
    integer :: n, p, k, steps, status, i
    integer, allocatable :: labels(:)
    real(dp) :: distortion
    real(dp), allocatable :: observations(:,:), guess(:,:), centers(:,:), distances(:), whitened(:,:)

    read(*,*) operation, n, p, k, steps
    allocate(observations(n,p),guess(k,p))
    do i=1,n
        read(*,*) observations(i,:)
    end do
    if (operation/='whiten') then
        do i=1,k
            read(*,*) guess(i,:)
        end do
    end if

    select case(operation)
    case('whiten')
        whitened = whiten(observations,status)
        print '(i0)', status
        do i=1,size(whitened,1)
            print '(*(es26.17,1x))', whitened(i,:)
        end do
    case('vq')
        call vq(observations,guess,labels,distances,status)
        print '(i0)', status
        print '(*(i0,1x))', labels
        print '(*(es26.17,1x))', distances
    case('kmeans2')
        call kmeans2(observations,guess,centers,labels,status,iterations=steps,missing='raise')
        print '(i0)', status
        if (status/=0) stop
        do i=1,size(centers,1)
            print '(*(es26.17,1x))', centers(i,:)
        end do
        print '(*(i0,1x))', labels
    case('kmeans')
        call kmeans(observations,guess,centers,distortion,status,threshold=1.0e-10_dp)
        print '(i0)', status
        if (status/=0) stop
        print '(i0)', size(centers,1)
        do i=1,size(centers,1)
            print '(*(es26.17,1x))', centers(i,:)
        end do
        print '(es26.17)', distortion
    end select
end program cluster_vq_driver
