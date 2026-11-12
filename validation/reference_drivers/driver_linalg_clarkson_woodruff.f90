! SPDX-License-Identifier: MIT
! Line-oriented protocol for compare_linalg_clarkson_woodruff_scipy.py.
program driver_linalg_clarkson_woodruff
    use scifort_kinds, only : dp
    use scifort_linalg_clarkson_woodruff, only : countsketch_plan, &
        clarkson_woodruff_transform, clarkson_woodruff_csr
    implicit none
    type(countsketch_plan) :: plan
    real(dp), allocatable :: a(:, :), res(:, :), ab(:, :, :), rb(:, :, :), vals(:)
    complex(dp), allocatable :: ac(:, :), rc(:, :), acb(:, :, :), rcb(:, :, :), cvals(:)
    integer, allocatable :: ptr(:), cols(:)
    real(dp) :: re, im
    integer :: test_count, testcase, kind, n, d, k, batch_size, nnz, status, b, i, j

    read(*,*) test_count
    do testcase = 1, test_count
        read(*,*) kind, n, d, k, batch_size, nnz
        plan%sketch_rows = k
        allocate(plan%bucket(n), plan%sign(n))
        do i = 1, n
            read(*,*) plan%bucket(i), plan%sign(i)
        end do
        if (kind == 2 .or. kind == 3) then
            allocate(ptr(n+1), cols(nnz), vals(nnz), cvals(nnz))
            do i = 1, n+1
                read(*,*) ptr(i)
            end do
            do i = 1, nnz
                read(*,*) cols(i)
            end do
            do i = 1, nnz
                read(*,*) re, im
                vals(i) = re
                cvals(i) = cmplx(re, im, dp)
            end do
            if (kind == 2) then
                allocate(res(k,d))
                call clarkson_woodruff_csr(ptr,cols,vals,d,k,res,status,plan)
            else
                allocate(rc(k,d))
                call clarkson_woodruff_csr(ptr,cols,cvals,d,k,rc,status,plan)
            end if
        else if (kind == 4 .or. kind == 5) then
            allocate(ab(batch_size,n,d),acb(batch_size,n,d))
            do b = 1, batch_size
                do i = 1, n
                    do j = 1, d
                        read(*,*) re, im
                        ab(b,i,j) = re
                        acb(b,i,j) = cmplx(re, im, dp)
                    end do
                end do
            end do
            if (kind == 4) then
                allocate(rb(batch_size,k,d))
                call clarkson_woodruff_transform(ab,k,rb,status,plan)
            else
                allocate(rcb(batch_size,k,d))
                call clarkson_woodruff_transform(acb,k,rcb,status,plan)
            end if
        else
            allocate(a(n,d),ac(n,d))
            do i = 1, n
                do j = 1, d
                    read(*,*) re, im
                    a(i,j) = re
                    ac(i,j) = cmplx(re, im, dp)
                end do
            end do
            if (kind == 0) then
                allocate(res(k,d))
                call clarkson_woodruff_transform(a,k,res,status,plan)
            else
                allocate(rc(k,d))
                call clarkson_woodruff_transform(ac,k,rc,status,plan)
            end if
        end if
        write(*,'(a,1x,i0,1x,i0)') 'CASE',testcase,status
        if (status == 0) then
            do b = 1, batch_size
                do i = 1, k
                    do j = 1, d
                        if (kind == 4) then
                            re = rb(b,i,j)
                            im = 0.0_dp
                        else if (kind == 5) then
                            re = real(rcb(b,i,j),dp)
                            im = aimag(rcb(b,i,j))
                        else if (kind == 0 .or. kind == 2) then
                            re = res(i,j)
                            im = 0.0_dp
                        else
                            re = real(rc(i,j),dp)
                            im = aimag(rc(i,j))
                        end if
                        write(*,'(2(es26.17e3,1x))') re, im
                    end do
                end do
            end do
        end if
        if (allocated(a)) deallocate(a)
        if (allocated(ac)) deallocate(ac)
        if (allocated(res)) deallocate(res)
        if (allocated(rc)) deallocate(rc)
        if (allocated(ab)) deallocate(ab)
        if (allocated(acb)) deallocate(acb)
        if (allocated(rb)) deallocate(rb)
        if (allocated(rcb)) deallocate(rcb)
        if (allocated(vals)) deallocate(vals)
        if (allocated(cvals)) deallocate(cvals)
        if (allocated(ptr)) deallocate(ptr)
        if (allocated(cols)) deallocate(cols)
        deallocate(plan%bucket, plan%sign)
    end do
end program driver_linalg_clarkson_woodruff
