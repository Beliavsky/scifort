! SPDX-License-Identifier: MIT
program legacy_reference_driver
    use scifort_kinds, only:dp
    use scifort_interpolate_legacy, only:interp1d,complex_interp1d_type
    implicit none
    type(complex_interp1d_type) :: p
    integer :: n,nch,nq,extr,i,j,code
    character(16) :: kind
    real(dp), allocatable :: x(:),q(:)
    complex(dp), allocatable :: vals(:,:),out(:,:)
    real(dp) :: r,im,low_r,low_i,high_r,high_i
    read(*,*) n,nch,nq,extr
    read(*,*) kind,low_r,low_i,high_r,high_i
    allocate(x(n),vals(n,nch),q(nq))
    read(*,*) x
    do i=1,n
        do j=1,nch
            read(*,*) r,im
            vals(i,j)=cmplx(r,im,dp)
        end do
    end do
    read(*,*) q
    if (extr==1) then
        p=interp1d(x,vals,kind=kind,extrapolate=.true.,status=code)
    else
        p=interp1d(x,vals,kind=kind,bounds_error=.false., &
            fill_below=cmplx(low_r,low_i,dp),fill_above=cmplx(high_r,high_i,dp),status=code)
    end if
    write(*,'(i0)') code
    if (code/=0) stop
    out=p%evaluate(q,status=code)
    write(*,'(i0)') code
    do i=1,nq
        do j=1,nch
            write(*,'(es25.16,1x,es25.16)') real(out(i,j),dp),aimag(out(i,j))
        end do
    end do
end program legacy_reference_driver
