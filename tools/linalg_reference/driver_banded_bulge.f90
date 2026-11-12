! SPDX-License-Identifier: MIT
! Streaming driver for opt-in full-spectrum banded inverse iteration.
program driver_banded_bulge
    use scifort_kinds, only : dp
    use scifort_linalg_banded_eigen, only : eig_banded
    implicit none
    integer :: count, item, typ, n, kd, found, status, i, j
    real(dp), allocatable :: ab(:, :), w(:), zreal(:, :)
    complex(dp), allocatable :: cab(:, :), zcomplex(:, :)
    real(dp) :: re, im
    logical :: low

    read(*,*) count
    do item=1,count
        read(*,*) typ,n,kd
        allocate(cab(kd+1,n),w(n),zcomplex(n,n))
        do j=1,n
            do i=1,kd+1
                read(*,*) re, im
                cab(i,j)=cmplx(re,im,dp)
            end do
        end do
        low=typ == 1 .or. typ == 3
        if (typ <= 2) then
            allocate(ab(kd+1,n),zreal(n,n))
            ab=real(cab,dp)
            call eig_banded(ab,w,status,zreal,lower=low, &
                num_selected=found,method='bulge')
            zcomplex=cmplx(zreal,0.0_dp,dp)
            deallocate(ab,zreal)
        else
            call eig_banded(cab,w,status,zcomplex,lower=low, &
                num_selected=found,method='bulge')
        end if
        write(*,'(i0,1x,i0)') status,found
        if (status == 0) then
            do j=1,n
                write(*,'(es25.16e3)') w(j)
            end do
            do j=1,n
                do i=1,n
                    write(*,'(es25.16e3,1x,es25.16e3)') &
                        real(zcomplex(i,j),dp),aimag(zcomplex(i,j))
                end do
            end do
        end if
        deallocate(cab,w,zcomplex)
    end do
end program driver_banded_bulge
