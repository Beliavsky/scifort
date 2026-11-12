! SPDX-License-Identifier: MIT
! Batch numerical reference driver for SciFort's range, kernel and angle APIs.
program driver_subspaces
    use scifort_kinds, only : dp
    use scifort_linalg_subspaces, only : orth, null_space, subspace_angles
    implicit none
    complex(dp), allocatable :: a(:, :), b(:, :), q(:, :), z(:, :)
    real(dp), allocatable :: qr(:, :), zr(:, :), angles(:)
    real(dp) :: xr, xi, rcond
    integer :: ncases, testno, kind, m, n, p, i, j, so, sz, sa

    read(*,*) ncases
    do testno = 1, ncases
        read(*,*) kind, m, n, p, rcond
        allocate(a(m,n),b(m,p))
        do j = 1, n
            do i = 1, m
                read(*,*) xr, xi
                a(i,j) = cmplx(xr,xi,dp)
            end do
        end do
        do j = 1, p
            do i = 1, m
                read(*,*) xr, xi
                b(i,j) = cmplx(xr,xi,dp)
            end do
        end do

        if (kind == 1 .or. kind == 3) then
            call orth(real(a,dp),qr,so,rcond)
            call null_space(real(a,dp),zr,sz,rcond)
            q = cmplx(qr,0.0_dp,dp)
            z = cmplx(zr,0.0_dp,dp)
            deallocate(qr,zr)
        else
            call orth(a,q,so,rcond)
            call null_space(a,z,sz,rcond)
        end if
        select case (kind)
        case (1)
            call subspace_angles(real(a,dp),real(b,dp),angles,sa,rcond)
        case (2)
            call subspace_angles(a,b,angles,sa,rcond)
        case (3)
            call subspace_angles(real(a,dp),b,angles,sa,rcond)
        case (4)
            call subspace_angles(a,real(b,dp),angles,sa,rcond)
        case default
            error stop 'invalid input kind'
        end select
        write(*,'(6(I0,1X))') so,sz,sa,size(q,2),size(z,2),size(angles)
        do j = 1, size(q,2)
            do i = 1, size(q,1)
                write(*,'(2(ES27.18E3,1X))') real(q(i,j),dp),aimag(q(i,j))
            end do
        end do
        do j = 1, size(z,2)
            do i = 1, size(z,1)
                write(*,'(2(ES27.18E3,1X))') real(z(i,j),dp),aimag(z(i,j))
            end do
        end do
        do i = 1, size(angles)
            write(*,'(ES27.18E3)') angles(i)
        end do
        deallocate(a,b,q,z,angles)
    end do
end program driver_subspaces
