! SPDX-License-Identifier: MIT
program rbf_reference_driver
    use scifort_kinds, only:dp
    use scifort_interpolate_rbf, only:rbf_interpolator,complex_rbf_interpolator_type
    implicit none
    type(complex_rbf_interpolator_type) :: p
    integer :: n,nd,nch,nq,degree,neighbors,i,j,code
    real(dp) :: eps
    character(32) :: kernel
    real(dp), allocatable :: y(:,:),queries(:,:),smoothing(:)
    complex(dp), allocatable :: vals(:,:),out(:,:)
    real(dp) :: r,im
    read(*,*) n,nd,nch,nq,degree,neighbors
    read(*,*) kernel,eps
    allocate(y(n,nd),vals(n,nch),queries(nq,nd),smoothing(n))
    do i=1,n
        read(*,*) smoothing(i),y(i,:)
    end do
    ! Values are deliberately read separately to make malformed input obvious.
    do i=1,n
        do j=1,nch
            read(*,*) r,im
            vals(i,j)=cmplx(r,im,dp)
        end do
    end do
    do i=1,nq
        read(*,*) queries(i,:)
    end do
    if (neighbors==0) then
        p=rbf_interpolator(y,vals,kernel=trim(kernel),epsilon=eps,degree=degree, &
            smoothing_values=smoothing,status=code)
    else
        p=rbf_interpolator(y,vals,kernel=trim(kernel),epsilon=eps,degree=degree, &
            smoothing_values=smoothing,neighbors=neighbors,status=code)
    end if
    write(*,'(i0)') code
    if (code/=0) stop
    out=p%evaluate(queries,status=code)
    write(*,'(i0)') code
    do i=1,nq
        do j=1,nch
            write(*,'(es25.16,1x,es25.16)') real(out(i,j),dp),aimag(out(i,j))
        end do
    end do
end program rbf_reference_driver
