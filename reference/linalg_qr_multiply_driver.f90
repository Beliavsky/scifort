! SPDX-License-Identifier: MIT
! Text protocol for independent Python/SciPy QR multiplication comparisons.
program linalg_qr_multiply_driver
    use scifort_kinds, only : dp
    use scifort_linalg_qr_multiply, only : qr_multiply
    implicit none
    real(dp), allocatable :: ar(:, :), cr(:, :), yr(:, :), rr(:, :)
    complex(dp), allocatable :: a(:, :), c(:, :), y(:, :), r(:, :)
    integer, allocatable :: piv(:)
    integer :: cases, t, m, n, k, p, da, dc, side, cp, pivoted, vec, status, i, j
    integer :: crow, ccol, yrow, ycol
    real(dp) :: re, im
    character(len=5) :: mode

    read(*,*) cases
    do t = 1, cases
        read(*,*) da, dc, side, cp, pivoted, vec, m, n, p
        k = min(m,n)
        mode = 'left'
        if (side == 1) mode = 'right'
        if (side == 0) then
            crow = k
            ccol = p
            yrow = m
            ycol = p
        else
            crow = p
            ccol = m
            yrow = p
            ycol = k
        end if
        allocate(a(m,n),c(crow,ccol),y(yrow,ycol),r(k,n),piv(n))
        do j = 1, n
            do i = 1, m
                read(*,*) re, im
                a(i,j) = cmplx(re,im,dp)
            end do
        end do
        do j = 1, ccol
            do i = 1, crow
                read(*,*) re, im
                c(i,j) = cmplx(re,im,dp)
            end do
        end do
        if (da == 0) allocate(ar(m,n))
        if (dc == 0) allocate(cr(crow,ccol))
        if (da == 0) ar = real(a,dp)
        if (dc == 0) cr = real(c,dp)
        if (da == 0 .and. dc == 0) then
            allocate(yr(yrow,ycol),rr(k,n))
            if (vec == 0) then
                if (pivoted == 1) then
                    call qr_multiply(ar,cr,yr,rr,status,mode,piv,cp==1)
                else
                    call qr_multiply(ar,cr,yr,rr,status,mode,conjugate=cp==1)
                end if
            else
                if (side == 0) then
                    if (pivoted == 1) then
                        call qr_multiply(ar,cr(:,1),yr(:,1),rr,status,mode,piv,cp==1)
                    else
                        call qr_multiply(ar,cr(:,1),yr(:,1),rr,status,mode,conjugate=cp==1)
                    end if
                else
                    if (pivoted == 1) then
                        call qr_multiply(ar,cr(1,:),yr(1,:),rr,status,mode,piv,cp==1)
                    else
                        call qr_multiply(ar,cr(1,:),yr(1,:),rr,status,mode,conjugate=cp==1)
                    end if
                end if
            end if
            r = cmplx(rr,0.0_dp,dp)
            y = cmplx(yr,0.0_dp,dp)
        else if (da == 0) then
            if (vec == 0) then
                if (pivoted == 1) then
                    call qr_multiply(ar,c,y,r,status,mode,piv,cp==1)
                else
                    call qr_multiply(ar,c,y,r,status,mode,conjugate=cp==1)
                end if
            else
                if (side == 0) then
                    if (pivoted == 1) then
                        call qr_multiply(ar,c(:,1),y(:,1),r,status,mode,piv,cp==1)
                    else
                        call qr_multiply(ar,c(:,1),y(:,1),r,status,mode,conjugate=cp==1)
                    end if
                else
                    if (pivoted == 1) then
                        call qr_multiply(ar,c(1,:),y(1,:),r,status,mode,piv,cp==1)
                    else
                        call qr_multiply(ar,c(1,:),y(1,:),r,status,mode,conjugate=cp==1)
                    end if
                end if
            end if
        else if (dc == 0) then
            if (vec == 0) then
                if (pivoted == 1) then
                    call qr_multiply(a,cr,y,r,status,mode,piv,cp==1)
                else
                    call qr_multiply(a,cr,y,r,status,mode,conjugate=cp==1)
                end if
            else
                if (side == 0) then
                    if (pivoted == 1) then
                        call qr_multiply(a,cr(:,1),y(:,1),r,status,mode,piv,cp==1)
                    else
                        call qr_multiply(a,cr(:,1),y(:,1),r,status,mode,conjugate=cp==1)
                    end if
                else
                    if (pivoted == 1) then
                        call qr_multiply(a,cr(1,:),y(1,:),r,status,mode,piv,cp==1)
                    else
                        call qr_multiply(a,cr(1,:),y(1,:),r,status,mode,conjugate=cp==1)
                    end if
                end if
            end if
        else
            if (vec == 0) then
                if (pivoted == 1) then
                    call qr_multiply(a,c,y,r,status,mode,piv,cp==1)
                else
                    call qr_multiply(a,c,y,r,status,mode,conjugate=cp==1)
                end if
            else
                if (side == 0) then
                    if (pivoted == 1) then
                        call qr_multiply(a,c(:,1),y(:,1),r,status,mode,piv,cp==1)
                    else
                        call qr_multiply(a,c(:,1),y(:,1),r,status,mode,conjugate=cp==1)
                    end if
                else
                    if (pivoted == 1) then
                        call qr_multiply(a,c(1,:),y(1,:),r,status,mode,piv,cp==1)
                    else
                        call qr_multiply(a,c(1,:),y(1,:),r,status,mode,conjugate=cp==1)
                    end if
                end if
            end if
        end if
        if (pivoted == 0) piv = [(i,i=1,n)]
        write(*,'(i0,*(1x,i0))') status, piv
        do j = 1, n
            do i = 1, k
                write(*,'(es25.17,1x,es25.17)') real(r(i,j),dp), aimag(r(i,j))
            end do
        end do
        do j = 1, ycol
            do i = 1, yrow
                write(*,'(es25.17,1x,es25.17)') real(y(i,j),dp), aimag(y(i,j))
            end do
        end do
        if (allocated(ar)) deallocate(ar)
        if (allocated(cr)) deallocate(cr)
        if (allocated(yr)) deallocate(yr)
        if (allocated(rr)) deallocate(rr)
        deallocate(a,c,y,r,piv)
    end do
end program linalg_qr_multiply_driver
