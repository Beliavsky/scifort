! SPDX-License-Identifier: MIT
! Streaming test interface for real/complex barycentric interpolation.
program barycentric_driver
    use scifort_kinds, only : dp
    use scifort_interpolate, only : barycentric_interpolator, &
        barycentric_interpolator_type, complex_barycentric_interpolator_type
    implicit none
    type(barycentric_interpolator_type) :: p
    type(complex_barycentric_interpolator_type) :: cp
    real(dp), allocatable :: x(:), q(:), nx(:), yr(:,:), yi(:,:), nyr(:,:), nyi(:,:)
    real(dp), allocatable :: result_real(:,:)
    complex(dp), allocatable :: result_complex(:,:)
    integer :: cases, idx, kind, n, nch, nq, na, phase, nu, status
    read(*,*) cases
    do idx=1,cases
        read(*,*) kind,n,nch,nq,na
        allocate(x(n),q(nq),nx(na),yr(n,nch),yi(n,nch),nyr(na,nch),nyi(na,nch))
        read(*,*) x
        read(*,*) yr
        read(*,*) yi
        read(*,*) q
        if (na>0) then
            read(*,*) nx
            read(*,*) nyr
            read(*,*) nyi
        end if
        if (kind==1) then
            p=barycentric_interpolator(x,yr,status=status)
        else
            cp=barycentric_interpolator(x,cmplx(yr,yi,dp),status=status)
        end if
        if (status/=0) error stop 'invalid barycentric data'
        write(*,'(a,1x,i0)') 'CASE',idx
        do phase=0,1
            if (phase==1 .and. na>0) then
                if (kind==1) then
                    call p%add_nodes(nx,nyr,status)
                else
                    call cp%add_nodes(nx,cmplx(nyr,nyi,dp),status)
                end if
                if (status/=0) error stop 'invalid barycentric insertion'
            end if
            do nu=0,3
                write(*,'(a,1x,i0,1x,i0)') 'D',phase,nu
                if (kind==1) then
                    result_real=p%values(q,nu)
                    write(*,'(*(1x,es26.17e3))') result_real
                    write(*,'(*(1x,es26.17e3))') 0.0_dp*result_real
                else
                    result_complex=cp%values(q,nu)
                    write(*,'(*(1x,es26.17e3))') real(result_complex,dp)
                    write(*,'(*(1x,es26.17e3))') aimag(result_complex)
                end if
            end do
        end do
        deallocate(x,q,nx,yr,yi,nyr,nyi)
    end do
end program barycentric_driver
