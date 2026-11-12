! SPDX-License-Identifier: MIT
! Machine-readable driver for independent SciPy comparisons of CSD.
program driver_cossin
    use scifort_kinds, only : dp
    use scifort_linalg_cossin, only : cossin, cossin_separate
    use scifort_linalg, only : linalg_status_success
    implicit none
    integer :: nt, it, m, p, q, flip, kind_case, status, parts_status, i, j, r
    real(dp) :: re, im, err, ue, ve
    complex(dp), allocatable :: xc(:, :), uc(:, :), vhc(:, :)
    complex(dp), allocatable :: u1(:, :),u2(:, :),v1h(:, :),v2h(:, :)
    real(dp), allocatable :: xr(:, :), ur(:, :), vtr(:, :), cs(:, :), th(:)
    complex(dp), allocatable :: gram(:, :)

    read(*,*) nt
    do it = 1, nt
        read(*,*) kind_case,m,p,q,flip
        allocate(xc(m,m),uc(m,m),vhc(m,m),xr(m,m),ur(m,m),vtr(m,m),cs(m,m),gram(m,m))
        do j = 1, m
            do i = 1, m
                read(*,*) re,im
                xc(i,j) = cmplx(re,im,dp)
            end do
        end do
        xr = real(xc,dp)
        if (kind_case == 0) then
            call cossin(xr,p,q,ur,cs,vtr,status,swap_sign=(flip==1))
        else
            call cossin(xc,p,q,uc,cs,vhc,status,swap_sign=(flip==1))
        end if
        err = 0.0_dp
        ue = 0.0_dp
        ve = 0.0_dp
        if (status == linalg_status_success) then
            if (kind_case == 0) then
                err = maxval(abs(matmul(matmul(ur,cs),vtr)-xr))
                uc = cmplx(ur,0.0_dp,dp)
                vhc = cmplx(vtr,0.0_dp,dp)
            else
                err = maxval(abs(matmul(matmul(uc,cmplx(cs,0.0_dp,dp)),vhc)-xc))
            end if
            gram = matmul(transpose(conjg(uc)),uc)
            do i=1,m
                gram(i,i) = gram(i,i)-1.0_dp
            end do
            ue = maxval(abs(gram))
            gram = matmul(transpose(conjg(vhc)),vhc)
            do i=1,m
                gram(i,i) = gram(i,i)-1.0_dp
            end do
            ve = maxval(abs(gram))
        end if
        write(*,'(i0,3(1x,es24.16))') status,err,ue,ve
        if (status == linalg_status_success) then
            r = min(p,m-p,q,m-q)
            allocate(u1(p,p),u2(m-p,m-p),v1h(q,q),v2h(m-q,m-q),th(r))
            call cossin_separate(xc,p,q,u1,u2,th,v1h,v2h,parts_status,swap_sign=(flip==1))
            if (parts_status /= linalg_status_success) error stop 'separate factor mismatch'
            write(*,'(*(es24.16,1x))') th
            write(*,'(*(es24.16,1x))') reshape(cs,[m*m])
            deallocate(u1,u2,v1h,v2h,th)
        end if
        deallocate(xc,uc,vhc,xr,ur,vtr,cs,gram)
    end do
end program driver_cossin
