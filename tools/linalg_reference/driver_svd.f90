! SPDX-License-Identifier: MIT
! Streaming real/complex SVD, pinv, and rank-deficient lstsq reference driver.
program driver_svd
    use scifort_kinds, only : dp
    use scifort_linalg_svd, only : svd, pinv
    use scifort_linalg_qr, only : lstsq
    implicit none
    complex(dp), allocatable :: a(:, :), b(:, :), u(:, :), vh(:, :)
    complex(dp), allocatable :: ap(:, :), x(:, :), xauto(:, :)
    real(dp), allocatable :: s(:)
    real(dp) :: re, im, rcond
    integer :: trials, case_number, type_id, m, n, nrhs, i, j, k
    integer :: st_svd, st_pinv, st_lstsq, st_auto, rank_pinv, rank_lstsq, rank_auto

    read(*,*) trials
    do case_number = 1, trials
        read(*,*) type_id, m, n, nrhs, rcond
        k = min(m,n)
        allocate(a(m,n),b(m,nrhs),u(m,k),s(k),vh(k,n),ap(n,m),x(n,nrhs),xauto(n,nrhs))
        do j = 1, n
            do i = 1, m
                read(*,*) re, im
                a(i,j) = cmplx(re,im,dp)
            end do
        end do
        do j = 1, nrhs
            do i = 1, m
                read(*,*) re, im
                b(i,j) = cmplx(re,im,dp)
            end do
        end do
        if (type_id == 1) then
            block
                real(dp), allocatable :: ar(:, :), br(:, :), ur(:, :), vhr(:, :)
                real(dp), allocatable :: apr(:, :), xr(:, :), xautor(:, :)
                allocate(ar(m,n),br(m,nrhs),ur(m,k),vhr(k,n))
                allocate(apr(n,m),xr(n,nrhs),xautor(n,nrhs))
                ar = real(a,dp)
                br = real(b,dp)
                call svd(ar,ur,s,vhr,st_svd)
                call pinv(ar,apr,st_pinv,rank_pinv,rcond)
                call lstsq(ar,br,xr,st_lstsq,rank_lstsq,rcond,method='svd')
                call lstsq(ar,br,xautor,st_auto,rank_auto,rcond)
                u = cmplx(ur,0.0_dp,dp)
                vh = cmplx(vhr,0.0_dp,dp)
                ap = cmplx(apr,0.0_dp,dp)
                x = cmplx(xr,0.0_dp,dp)
                xauto = cmplx(xautor,0.0_dp,dp)
            end block
        else
            call svd(a,u,s,vh,st_svd)
            call pinv(a,ap,st_pinv,rank_pinv,rcond)
            call lstsq(a,b,x,st_lstsq,rank_lstsq,rcond,method='svd')
            call lstsq(a,b,xauto,st_auto,rank_auto,rcond)
        end if
        print '(7(I0,1X))', st_svd,st_pinv,st_lstsq,st_auto,rank_pinv,rank_lstsq,rank_auto
        do i = 1, k
            print '(ES26.17E3)', s(i)
        end do
        do j = 1, k
            do i = 1, m
                call print_complex(u(i,j))
            end do
        end do
        do j = 1, n
            do i = 1, k
                call print_complex(vh(i,j))
            end do
        end do
        do j = 1, m
            do i = 1, n
                call print_complex(ap(i,j))
            end do
        end do
        do j = 1, nrhs
            do i = 1, n
                call print_complex(x(i,j))
            end do
        end do
        do j = 1, nrhs
            do i = 1, n
                call print_complex(xauto(i,j))
            end do
        end do
        deallocate(a,b,u,s,vh,ap,x,xauto)
    end do
contains
    subroutine print_complex(z)
        complex(dp), intent(in) :: z !! finite complex value
        print '(ES26.17E3,1X,ES26.17E3)', real(z,dp),aimag(z)
    end subroutine print_complex
end program driver_svd
