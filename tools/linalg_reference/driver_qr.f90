! SPDX-License-Identifier: MIT
! Text transport for independent real/complex QR and least-squares tests.
program driver_qr
    use scifort_kinds, only : dp
    use scifort_linalg_qr, only : qr, qr_factor, qr_unpack, qr_solve, lstsq
    implicit none
    real(dp), allocatable :: a(:, :), q(:, :), r(:, :), b(:, :), x(:, :)
    real(dp), allocatable :: af(:, :), tau(:), xr(:, :)
    complex(dp), allocatable :: ac(:, :), qc(:, :), rc(:, :), bc(:, :), xc(:, :)
    complex(dp), allocatable :: afc(:, :), tauc(:), xrc(:, :)
    integer, allocatable :: piv(:)
    integer :: ncase, icase, variant, m, n, nrhs, pivoted, complete
    integer :: kmax, qcols, i, j, status, qr_status, ls_status, rank
    real(dp) :: ar, ai

    read(*,*) ncase
    do icase = 1, ncase
        read(*,*) variant, m, n, nrhs, pivoted, complete
        kmax = min(m,n)
        qcols = kmax
        if (complete /= 0) qcols = m
        if (allocated(a)) deallocate(a,q,r,b,x,af,tau,xr,ac,qc,rc,bc,xc,afc,tauc,xrc,piv)
        allocate(a(m,n),q(m,qcols),r(qcols,n),b(m,nrhs),x(n,nrhs),xr(n,nrhs),af(m,n),tau(kmax))
        allocate(ac(m,n),qc(m,qcols),rc(qcols,n),bc(m,nrhs),xc(n,nrhs),xrc(n,nrhs),afc(m,n),tauc(kmax),piv(n))
        do j = 1, n
            do i = 1, m
                read(*,*) ar, ai
                a(i,j) = ar
                ac(i,j) = cmplx(ar,ai,dp)
            end do
        end do
        do j = 1, nrhs
            do i = 1, m
                read(*,*) ar, ai
                b(i,j) = ar
                bc(i,j) = cmplx(ar,ai,dp)
            end do
        end do
        do i = 1, n
            piv(i) = i
        end do
        if (variant == 1) then
            if (pivoted /= 0) then
                call qr_factor(a,af,tau,status,piv)
            else
                call qr_factor(a,af,tau,status)
            end if
            if (status == 0) call qr_unpack(af,tau,q,r,status)
            qr_status = status
            call lstsq(a,b,x,ls_status,rank)
            if (ls_status == 0 .and. m >= n .and. pivoted /= 0) then
                call qr_solve(af,tau,piv,b,xr,status)
                if (status /= 0) error stop 'reused real QR solve'
                if (maxval(abs(xr-x)) > 2.0e-7_dp*max(1.0_dp,maxval(abs(x)))) &
                    error stop 'reused real QR mismatch'
            end if
        else
            if (pivoted /= 0) then
                call qr_factor(ac,afc,tauc,status,piv)
            else
                call qr_factor(ac,afc,tauc,status)
            end if
            if (status == 0) call qr_unpack(afc,tauc,qc,rc,status)
            qr_status = status
            call lstsq(ac,bc,xc,ls_status,rank)
            if (ls_status == 0 .and. m >= n .and. pivoted /= 0) then
                call qr_solve(afc,tauc,piv,bc,xrc,status)
                if (status /= 0) error stop 'reused complex QR solve'
                if (maxval(abs(xrc-xc)) > 2.0e-7_dp*max(1.0_dp,maxval(abs(xc)))) &
                    error stop 'reused complex QR mismatch'
            end if
        end if
        write(*,'(I0,1X,I0,1X,I0)') qr_status, ls_status, rank
        if (qr_status /= 0) cycle
        write(*,'(*(I0,1X))') piv
        do j = 1, qcols
            do i = 1, m
                if (variant == 1) then
                    write(*,'(ES28.18,1X,ES28.18)') q(i,j),0.0_dp
                else
                    write(*,'(ES28.18,1X,ES28.18)') real(qc(i,j),dp),aimag(qc(i,j))
                end if
            end do
        end do
        do j = 1, n
            do i = 1, qcols
                if (variant == 1) then
                    write(*,'(ES28.18,1X,ES28.18)') r(i,j),0.0_dp
                else
                    write(*,'(ES28.18,1X,ES28.18)') real(rc(i,j),dp),aimag(rc(i,j))
                end if
            end do
        end do
        if (ls_status /= 0) cycle
        do j = 1, nrhs
            do i = 1, n
                if (variant == 1) then
                    write(*,'(ES28.18,1X,ES28.18)') x(i,j),0.0_dp
                else
                    write(*,'(ES28.18,1X,ES28.18)') real(xc(i,j),dp),aimag(xc(i,j))
                end if
            end do
        end do
    end do
end program driver_qr
