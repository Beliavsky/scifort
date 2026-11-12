! SPDX-License-Identifier: MIT
! Streaming reference driver for banded full-spectrum auto dispatch.
program driver_banded_auto
    use scifort_kinds, only : dp
    use scifort_linalg_banded_eigen, only : eigvals_banded, eig_banded
    implicit none
    integer :: count, item, typ, n, kd, i, j
    integer :: sa, sb, sd, se
    real(dp), allocatable :: ab(:, :), wa(:), wb(:), wd(:), we(:)
    complex(dp), allocatable :: cab(:, :)
    real(dp) :: re, im
    logical :: low

    read(*, *) count
    do item = 1, count
        read(*, *) typ, n, kd
        allocate(cab(kd+1,n),wa(n),wb(n),wd(n),we(n))
        do j = 1,n
            do i = 1,kd+1
                read(*, *) re, im
                cab(i,j) = cmplx(re,im,dp)
            end do
        end do
        low = typ == 1 .or. typ == 3
        if (typ <= 2) then
            allocate(ab(kd+1,n))
            ab = real(cab,dp)
            call eigvals_banded(ab,wa,sa,lower=low,method='auto')
            call eigvals_banded(ab,wb,sb,lower=low,method='bulge')
            call eigvals_banded(ab,wd,sd,lower=low,method='dense')
            call eig_banded(ab,we,se,lower=low,method='auto')
            deallocate(ab)
        else
            call eigvals_banded(cab,wa,sa,lower=low,method='auto')
            call eigvals_banded(cab,wb,sb,lower=low,method='bulge')
            call eigvals_banded(cab,wd,sd,lower=low,method='dense')
            call eig_banded(cab,we,se,lower=low,method='auto')
        end if
        write(*,'(4(i0,1x))') sa,sb,sd,se
        do j = 1,n
            write(*,'(4(es25.16e3,1x))') wa(j),wb(j),wd(j),we(j)
        end do
        deallocate(cab,wa,wb,wd,we)
    end do
end program driver_banded_auto
