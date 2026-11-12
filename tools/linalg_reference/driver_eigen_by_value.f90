! SPDX-License-Identifier: MIT
! Batch protocol for half-open eigenvalue-interval SciPy comparisons.
program driver_eigen_by_value
    use scifort_kinds, only : dp
    use scifort_linalg_dense_eigen, only : eigh, eigvalsh
    implicit none
    integer :: ncases, icase, kind_id, n, lower_id, problem_type, want_vectors
    integer :: cap_w, cap_z, i, j, status, nselected
    real(dp) :: lo, hi, xr, xi
    real(dp), allocatable :: a(:, :), b(:, :), w(:), z(:, :)
    complex(dp), allocatable :: ca(:, :), cb(:, :), cz(:, :)
    logical :: lower

    read(*, *) ncases
    do icase = 1, ncases
        read(*, *) kind_id, n, lower_id, problem_type, want_vectors, cap_w, cap_z
        read(*, *) lo, hi
        lower = lower_id == 1
        allocate(a(n,n), b(n,n), ca(n,n), cb(n,n), w(cap_w))
        do j = 1, n
            do i = 1, n
                read(*, *) xr, xi
                a(i,j) = xr
                ca(i,j) = cmplx(xr,xi,dp)
            end do
        end do
        do j = 1, n
            do i = 1, n
                read(*, *) xr, xi
                b(i,j) = xr
                cb(i,j) = cmplx(xr,xi,dp)
            end do
        end do
        if (kind_id == 1) then
            if (want_vectors == 1) then
                allocate(z(n,cap_z))
                if (problem_type == 0) then
                    call eigh(a,w,status,z,lower=lower, &
                        subset_by_value=[lo,hi],num_selected=nselected)
                else
                    call eigh(a,w,status,z,lower=lower,b=b,problem_type=problem_type, &
                        subset_by_value=[lo,hi],num_selected=nselected)
                end if
            else
                if (problem_type == 0) then
                    call eigvalsh(a,w,status,lower=lower, &
                        subset_by_value=[lo,hi],num_selected=nselected)
                else
                    call eigvalsh(a,w,status,lower=lower,b=b,problem_type=problem_type, &
                        subset_by_value=[lo,hi],num_selected=nselected)
                end if
            end if
        else
            if (want_vectors == 1) then
                allocate(cz(n,cap_z))
                if (problem_type == 0) then
                    call eigh(ca,w,status,cz,lower=lower, &
                        subset_by_value=[lo,hi],num_selected=nselected)
                else
                    call eigh(ca,w,status,cz,lower=lower,b=cb,problem_type=problem_type, &
                        subset_by_value=[lo,hi],num_selected=nselected)
                end if
            else
                if (problem_type == 0) then
                    call eigvalsh(ca,w,status,lower=lower, &
                        subset_by_value=[lo,hi],num_selected=nselected)
                else
                    call eigvalsh(ca,w,status,lower=lower,b=cb,problem_type=problem_type, &
                        subset_by_value=[lo,hi],num_selected=nselected)
                end if
            end if
        end if
        print '(i0,1x,i0)', status,nselected
        if (status == 0) then
            do i = 1, nselected
                print '(es26.17e3)', w(i)
            end do
            if (want_vectors == 1) then
                do j = 1, nselected
                    do i = 1, n
                        if (kind_id == 1) then
                            print '(2(es26.17e3,1x))', z(i,j), 0.0_dp
                        else
                            print '(2(es26.17e3,1x))', real(cz(i,j),dp),aimag(cz(i,j))
                        end if
                    end do
                end do
            end if
        end if
        if (allocated(z)) deallocate(z)
        if (allocated(cz)) deallocate(cz)
        deallocate(a,b,ca,cb,w)
    end do
end program driver_eigen_by_value
