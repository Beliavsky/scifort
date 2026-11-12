! SPDX-License-Identifier: MIT
! Stream real and complex matrices into the Hessenberg reduction.
program driver_hessenberg
    use scifort_kinds, only : dp
    use scifort_linalg_hessenberg, only : hessenberg
    implicit none
    integer :: cases, kind_case, n, c, i, j, status
    real(dp) :: ar, ai
    real(dp), allocatable :: a(:, :), h(:, :), q(:, :)
    complex(dp), allocatable :: ac(:, :), hc(:, :), qc(:, :)

    read(*, *) cases
    do c = 1, cases
        read(*, *) kind_case, n
        allocate(a(n,n),h(n,n),q(n,n),ac(n,n),hc(n,n),qc(n,n))
        do j = 1, n
            do i = 1, n
                read(*, *) ar, ai
                a(i,j) = ar
                ac(i,j) = cmplx(ar,ai,dp)
            end do
        end do
        if (kind_case == 1) then
            call hessenberg(a,h,status,q)
            hc = cmplx(h,0.0_dp,dp)
            qc = cmplx(q,0.0_dp,dp)
        else if (kind_case == 2) then
            call hessenberg(ac,hc,status,qc)
        else
            error stop 'unknown Hessenberg kind'
        end if
        write(*,'(I0)') status
        do j = 1, n
            do i = 1, n
                write(*,'(ES25.16E3,1X,ES25.16E3)') real(hc(i,j),dp),aimag(hc(i,j))
            end do
        end do
        do j = 1, n
            do i = 1, n
                write(*,'(ES25.16E3,1X,ES25.16E3)') real(qc(i,j),dp),aimag(qc(i,j))
            end do
        end do
        deallocate(a,h,q,ac,hc,qc)
    end do
end program driver_hessenberg
