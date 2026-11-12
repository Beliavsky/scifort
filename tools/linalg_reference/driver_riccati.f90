! SPDX-License-Identifier: MIT
! Streamed reference driver: one-based Fortran arrays, both real/complex APIs.
program driver_riccati
    use scifort_kinds, only : dp
    use scifort_linalg_riccati, only : solve_continuous_are, solve_discrete_are
    implicit none
    character(len=1) :: dtype
    integer :: equation, n, m, i, j, ios, status, steps
    real(dp) :: residual
    real(dp), allocatable :: ar(:,:), br(:,:), qr(:,:), rr(:,:), er(:,:), sr(:,:), xr(:,:)
    complex(dp), allocatable :: ac(:,:), bc(:,:), qc(:,:), rc(:,:), ec(:,:), sc(:,:), xc(:,:)
    do
        read(*,*,iostat=ios) equation,dtype,n,m
        if (ios/=0) exit
        if (n<1 .or. m<1 .or. n>20 .or. m>20) error stop 'invalid driver case'
        if (dtype=='R') then
            allocate(ar(n,n),br(n,m),qr(n,n),rr(m,m),er(n,n),sr(n,m),xr(n,n))
            read(*,*) ar
            read(*,*) br
            read(*,*) qr
            read(*,*) rr
            read(*,*) er
            read(*,*) sr
            if (equation==1) then
                call solve_continuous_are(ar,br,qr,rr,xr,status,e=er,s=sr,iterations=steps,residual=residual)
            else
                call solve_discrete_are(ar,br,qr,rr,xr,status,e=er,s=sr,iterations=steps,residual=residual)
            end if
            print '(i0,1x,i0,1x,es23.14)', status,steps,residual
            do j=1,n
                do i=1,n
                    print '(2es26.16)', xr(i,j), 0.0_dp
                end do
            end do
            deallocate(ar,br,qr,rr,er,sr,xr)
        else
            allocate(ac(n,n),bc(n,m),qc(n,n),rc(m,m),ec(n,n),sc(n,m),xc(n,n))
            read(*,*) ac
            read(*,*) bc
            read(*,*) qc
            read(*,*) rc
            read(*,*) ec
            read(*,*) sc
            if (equation==1) then
                call solve_continuous_are(ac,bc,qc,rc,xc,status,e=ec,s=sc,iterations=steps,residual=residual)
            else
                call solve_discrete_are(ac,bc,qc,rc,xc,status,e=ec,s=sc,iterations=steps,residual=residual)
            end if
            print '(i0,1x,i0,1x,es23.14)', status,steps,residual
            do j=1,n
                do i=1,n
                    print '(2es26.16)', real(xc(i,j),dp), aimag(xc(i,j))
                end do
            end do
            deallocate(ac,bc,qc,rc,ec,sc,xc)
        end if
    end do
end program driver_riccati
