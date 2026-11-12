! SPDX-License-Identifier: MIT
! Bulk I/O for independent SciPy comparisons of nonsymmetric eigenproblems.
program driver_general_eigen
    use scifort_kinds, only : dp
    use scifort_linalg_general_eigen, only : eig, eigvals
    implicit none
    real(dp), allocatable :: a(:, :)
    complex(dp), allocatable :: ca(:, :), w(:), w2(:), vr(:, :), vl(:, :)
    real(dp) :: xr, xi
    integer :: count, idx, kind_id, n, i, j, status, status2

    read(*, *) count
    do idx = 1, count
        read(*, *) kind_id, n
        allocate(a(n,n), ca(n,n), w(n), w2(n), vr(n,n), vl(n,n))
        do j = 1, n
            do i = 1, n
                read(*, *) xr, xi
                a(i,j) = xr
                ca(i,j) = cmplx(xr,xi,dp)
            end do
        end do
        if (kind_id == 1) then
            call eig(a,w,status,vr=vr,vl=vl)
            call eigvals(a,w2,status2)
        else
            call eig(ca,w,status,vr=vr,vl=vl)
            call eigvals(ca,w2,status2)
        end if
        print '(2(i0,1x))', status, status2
        if (status == 0 .and. status2 == 0) then
            do i = 1, n
                print '(4(es26.17e3,1x))', real(w(i),dp), aimag(w(i)), real(w2(i),dp), aimag(w2(i))
            end do
            do j = 1, n
                do i = 1, n
                    print '(4(es26.17e3,1x))', real(vr(i,j),dp), aimag(vr(i,j)), &
                        real(vl(i,j),dp), aimag(vl(i,j))
                end do
            end do
        end if
        deallocate(a, ca, w, w2, vr, vl)
    end do
end program driver_general_eigen
