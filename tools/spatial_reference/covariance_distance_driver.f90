! SPDX-License-Identifier: MIT
! Independent driver for SciPy comparisons of three spatial distance metrics.
program covariance_distance_driver
    use scifort_kinds, only : dp
    use scifort_spatial, only : cdist, pdist
    implicit none
    real(dp), allocatable :: xa(:,:), xb(:,:), variance(:), inv_cov(:,:), d(:,:), v(:)
    character(len=32) :: name
    real(dp) :: base
    integer :: ncase, k, m, n, ndim, i, j, parameter_mode, cd, pd

    read(*,*) ncase
    do k=1,ncase
        read(*,*) n,m,ndim,name,parameter_mode,base
        allocate(xa(n,ndim),xb(m,ndim),variance(ndim),inv_cov(ndim,ndim))
        if (name=='seuclidean' .and. parameter_mode==1) read(*,*) variance
        if (name=='mahalanobis' .and. parameter_mode==1) then
            do i=1,ndim
                read(*,*) inv_cov(i,:)
            end do
        end if
        do i=1,n
            read(*,*) xa(i,:)
        end do
        do i=1,m
            read(*,*) xb(i,:)
        end do
        select case (name)
        case ('seuclidean')
            if (parameter_mode==1) then
                d=cdist(xa,xb,name,variance=variance,status=cd)
                v=pdist(xa,name,variance=variance,status=pd)
            else
                d=cdist(xa,xb,name,status=cd)
                v=pdist(xa,name,status=pd)
            end if
        case ('mahalanobis')
            if (parameter_mode==1) then
                d=cdist(xa,xb,name,inv_cov=inv_cov,status=cd)
                v=pdist(xa,name,inv_cov=inv_cov,status=pd)
            else
                d=cdist(xa,xb,name,status=cd)
                v=pdist(xa,name,status=pd)
            end if
        case ('jensenshannon')
            if (base>0.0_dp) then
                d=cdist(xa,xb,name,base=base,status=cd)
                v=pdist(xa,name,base=base,status=pd)
            else
                d=cdist(xa,xb,name,status=cd)
                v=pdist(xa,name,status=pd)
            end if
        case default
            error stop 'unsupported driver metric'
        end select
        write(*,'(i0,1x,i0)') cd,pd
        if (cd==0) then
            do i=1,n
                do j=1,m
                    write(*,'(es27.17e3)') d(i,j)
                end do
            end do
        end if
        if (pd==0) then
            do i=1,size(v)
                write(*,'(es27.17e3)') v(i)
            end do
        end if
        deallocate(xa,xb,variance,inv_cov,d,v)
    end do
end program covariance_distance_driver
