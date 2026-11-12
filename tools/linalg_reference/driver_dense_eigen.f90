! SPDX-License-Identifier: MIT
! Protocol-based driver for independent SciPy eigh comparisons.
program driver_dense_eigen
    use scifort_kinds, only : dp
    use scifort_linalg_dense_eigen, only : eigh, eigvalsh
    implicit none
    integer :: count, case_id, n, kind_id, lower_id, first, last, want_z
    integer :: i, j, m, status
    real(dp) :: xr, xi
    real(dp), allocatable :: a(:, :), w(:), z(:, :)
    complex(dp), allocatable :: ca(:, :), cz(:, :)
    logical :: low

    read(*,*) count
    do case_id = 1, count
        read(*,*) kind_id,n,lower_id,first,last,want_z
        low = lower_id == 1
        m = n
        if (first /= 0) m = last-first+1
        allocate(a(n,n),ca(n,n),w(m))
        do j=1,n
            do i=1,n
                read(*,*) xr, xi
                a(i,j) = xr
                ca(i,j) = cmplx(xr,xi,dp)
            end do
        end do
        if (kind_id == 1) then
            if (want_z == 1) then
                allocate(z(n,m))
                if (first == 0) then
                    call eigh(a,w,status,z,lower=low)
                else
                    call eigh(a,w,status,z,lower=low,subset_by_index=[first,last])
                end if
            else
                if (first == 0) then
                    call eigvalsh(a,w,status,lower=low)
                else
                    call eigvalsh(a,w,status,lower=low,subset_by_index=[first,last])
                end if
            end if
        else
            if (want_z == 1) then
                allocate(cz(n,m))
                if (first == 0) then
                    call eigh(ca,w,status,cz,lower=low)
                else
                    call eigh(ca,w,status,cz,lower=low,subset_by_index=[first,last])
                end if
            else
                if (first == 0) then
                    call eigvalsh(ca,w,status,lower=low)
                else
                    call eigvalsh(ca,w,status,lower=low,subset_by_index=[first,last])
                end if
            end if
        end if
        print '(i0)',status
        if (status == 0) then
            do i=1,m
                print '(es26.17e3)',w(i)
            end do
            if (want_z == 1) then
                do j=1,m
                    do i=1,n
                        if (kind_id == 1) then
                            print '(2(es26.17e3,1x))',z(i,j),0.0_dp
                        else
                            print '(2(es26.17e3,1x))',real(cz(i,j),dp),aimag(cz(i,j))
                        end if
                    end do
                end do
            end if
        end if
        if (allocated(z)) deallocate(z)
        if (allocated(cz)) deallocate(cz)
        deallocate(a,ca,w)
    end do
end program driver_dense_eigen
