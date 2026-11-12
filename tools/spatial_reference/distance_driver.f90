! SPDX-License-Identifier: MIT
! Numerical comparison driver; one input stream with many cdist/pdist cases.
program spatial_distance_driver
    use scifort_kinds, only : dp
    use scifort_spatial, only : cdist, pdist
    implicit none
    real(dp), allocatable :: xa(:,:),xb(:,:),w(:),d(:,:),v(:)
    character(len=32) :: name
    real(dp) :: p
    integer :: ncase, k, m, n, ndim, i, j, weighted, cd, pd

    read(*,*) ncase
    do k=1,ncase
        read(*,*) n,m,ndim,name,weighted,p
        allocate(xa(n,ndim),xb(m,ndim),w(ndim))
        if (weighted == 1) read(*,*) w
        do i=1,n
            read(*,*) xa(i,:)
        end do
        do i=1,m
            read(*,*) xb(i,:)
        end do
        if (name=='minkowski') then
            if (weighted==1) then
                d=cdist(xa,xb,name,p=p,w=w,status=cd)
                v=pdist(xa,name,p=p,w=w,status=pd)
            else
                d=cdist(xa,xb,name,p=p,status=cd)
                v=pdist(xa,name,p=p,status=pd)
            end if
        else
            if (weighted==1) then
                d=cdist(xa,xb,name,w=w,status=cd)
                v=pdist(xa,name,w=w,status=pd)
            else
                d=cdist(xa,xb,name,status=cd)
                v=pdist(xa,name,status=pd)
            end if
        end if
        write(*,'(i0,1x,i0)') cd,pd
        do i=1,n
            do j=1,m
                write(*,'(es27.17e3)') d(i,j)
            end do
        end do
        do i=1,size(v)
            write(*,'(es27.17e3)') v(i)
        end do
        deallocate(xa,xb,w,d,v)
    end do
end program spatial_distance_driver
