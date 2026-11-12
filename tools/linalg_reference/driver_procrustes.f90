! SPDX-License-Identifier: MIT
! Stream real, complex, and mixed orthogonal/unitary Procrustes inputs.
program driver_procrustes
    use scifort_kinds, only : dp
    use scifort_linalg_procrustes, only : orthogonal_procrustes
    implicit none
    integer :: cases, kind_case, m, n, c, i, j, status
    real(dp) :: ar, ai, br, bi, factor
    real(dp), allocatable :: a(:, :), b(:, :), r(:, :)
    complex(dp), allocatable :: ac(:, :), bc(:, :), rc(:, :)

    read(*, *) cases
    do c = 1, cases
        read(*, *) kind_case, m, n
        allocate(a(m,n), b(m,n), r(n,n), ac(m,n), bc(m,n), rc(n,n))
        do j = 1, n
            do i = 1, m
                read(*, *) ar, ai, br, bi
                a(i,j) = ar
                b(i,j) = br
                ac(i,j) = cmplx(ar,ai,dp)
                bc(i,j) = cmplx(br,bi,dp)
            end do
        end do
        select case (kind_case)
        case (1)
            call orthogonal_procrustes(a,b,r,factor,status)
            rc = cmplx(r,0.0_dp,dp)
        case (2)
            call orthogonal_procrustes(ac,bc,rc,factor,status)
        case (3)
            call orthogonal_procrustes(a,bc,rc,factor,status)
        case (4)
            call orthogonal_procrustes(ac,b,rc,factor,status)
        case default
            error stop 'unknown input case'
        end select
        write(*,'(I0,1X,ES25.16E3)') status,factor
        do j = 1, n
            do i = 1, n
                write(*,'(ES25.16E3,1X,ES25.16E3)') real(rc(i,j),dp),aimag(rc(i,j))
            end do
        end do
        deallocate(a,b,r,ac,bc,rc)
    end do
end program driver_procrustes
