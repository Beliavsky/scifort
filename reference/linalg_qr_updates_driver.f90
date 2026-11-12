! SPDX-License-Identifier: MIT
! Text protocol for independent Python/SciPy validation of QR modifications.
program linalg_qr_updates_driver
    use scifort_kinds, only : dp
    use scifort_linalg_qr_updates, only : qr_update, qr_insert, qr_delete
    implicit none
    complex(dp), allocatable :: q(:, :), r(:, :), u(:, :), v(:, :), q1(:, :), r1(:, :)
    real(dp), allocatable :: rq(:, :), rr(:, :), ru(:, :), rv(:, :), rq1(:, :), rr1(:, :)
    real(dp) :: re, im
    integer :: number_cases, t, dtype, op, m, n, h, p, k, mnew, nnew, hnew, s, i, j
    character(len=3) :: direction

    read(*,*) number_cases
    do t = 1, number_cases
        read(*,*) dtype, op, m, n, h, p, k
        allocate(q(m,h),r(h,n))
        do j = 1, h
            do i = 1, m
                read(*,*) re, im
                q(i,j) = cmplx(re,im,dp)
            end do
        end do
        do j = 1, n
            do i = 1, h
                read(*,*) re, im
                r(i,j) = cmplx(re,im,dp)
            end do
        end do
        select case(op)
        case(1,2)
            allocate(u(m,p),v(n,p))
        case(3,4)
            allocate(u(p,n),v(0,0))
        case(5,6)
            allocate(u(m,p),v(0,0))
        case default
            allocate(u(0,0),v(0,0))
        end select
        do j = 1, size(u,2)
            do i = 1, size(u,1)
                read(*,*) re, im
                u(i,j) = cmplx(re,im,dp)
            end do
        end do
        do j = 1, size(v,2)
            do i = 1, size(v,1)
                read(*,*) re, im
                v(i,j) = cmplx(re,im,dp)
            end do
        end do
        mnew = m
        nnew = n
        if (op == 3 .or. op == 4) mnew = m + p
        if (op == 5 .or. op == 6) nnew = n + p
        if (op == 7) mnew = m - p
        if (op == 8) nnew = n - p
        hnew = min(mnew,nnew)
        if (h == m) hnew = mnew
        allocate(q1(mnew,hnew),r1(hnew,nnew))
        if (dtype == 0) then
            rq = real(q,dp)
            rr = real(r,dp)
            ru = real(u,dp)
            rv = real(v,dp)
            allocate(rq1(mnew,hnew),rr1(hnew,nnew))
            select case(op)
            case(1)
                call qr_update(rq,rr,ru(:,1),rv(:,1),rq1,rr1,s)
            case(2)
                call qr_update(rq,rr,ru,rv,rq1,rr1,s)
            case(3)
                call qr_insert(rq,rr,ru(1,:),k,rq1,rr1,s,which='row')
            case(4)
                call qr_insert(rq,rr,ru,k,rq1,rr1,s,which='row')
            case(5)
                call qr_insert(rq,rr,ru(:,1),k,rq1,rr1,s,which='col')
            case(6)
                call qr_insert(rq,rr,ru,k,rq1,rr1,s,which='col')
            case(7)
                call qr_delete(rq,rr,k,rq1,rr1,s,p=p,which='row')
            case(8)
                call qr_delete(rq,rr,k,rq1,rr1,s,p=p,which='col')
            end select
            q1 = cmplx(rq1,0.0_dp,dp)
            r1 = cmplx(rr1,0.0_dp,dp)
            deallocate(rq1,rr1)
        else
            select case(op)
            case(1)
                call qr_update(q,r,u(:,1),v(:,1),q1,r1,s)
            case(2)
                call qr_update(q,r,u,v,q1,r1,s)
            case(3)
                call qr_insert(q,r,u(1,:),k,q1,r1,s,which='row')
            case(4)
                call qr_insert(q,r,u,k,q1,r1,s,which='row')
            case(5)
                call qr_insert(q,r,u(:,1),k,q1,r1,s,which='col')
            case(6)
                call qr_insert(q,r,u,k,q1,r1,s,which='col')
            case(7)
                call qr_delete(q,r,k,q1,r1,s,p=p,which='row')
            case(8)
                call qr_delete(q,r,k,q1,r1,s,p=p,which='col')
            end select
        end if
        write(*,'(4(I0,1X))') s,mnew,nnew,hnew
        do j = 1, size(q1,2)
            do i = 1, size(q1,1)
                write(*,'(2(ES26.17E3,1X))') real(q1(i,j),dp),aimag(q1(i,j))
            end do
        end do
        do j = 1, size(r1,2)
            do i = 1, size(r1,1)
                write(*,'(2(ES26.17E3,1X))') real(r1(i,j),dp),aimag(r1(i,j))
            end do
        end do
        deallocate(q,r,u,v,q1,r1)
    end do
end program linalg_qr_updates_driver
