! SPDX-License-Identifier: MIT
! I/O driver for independent SciPy comparisons of complex Schur form.
program driver_schur
    use scifort_kinds, only : dp
    use scifort_linalg_schur, only : schur
    implicit none

    real(dp), allocatable :: r(:, :)
    complex(dp), allocatable :: a(:, :), t(:, :), z(:, :)
    real(dp) :: re, im
    integer :: k, cases, kind, n, rule, i, j, status, selected
    character(len=4) :: sortcode

    read (*,*) cases
    do k=1,cases
        read (*,*) kind,n,rule
        allocate(a(n,n),t(n,n),z(n,n))
        do j=1,n
            do i=1,n
                read (*,*) re,im
                a(i,j) = cmplx(re,im,dp)
            end do
        end do
        select case (rule)
        case (1)
            sortcode = 'lhp'
        case (2)
            sortcode = 'rhp'
        case (3)
            sortcode = 'iuc'
        case (4)
            sortcode = 'ouc'
        case default
            sortcode = 'none'
        end select
        if (kind == 1) then
            allocate(r(n,n))
            r = real(a,dp)
            if (rule == 5) then
                call schur(r,t,z,status,selector=positive_imag,sdim=selected)
            else
                call schur(r,t,z,status,sort=trim(sortcode),sdim=selected)
            end if
            deallocate(r)
        else
            if (rule == 5) then
                call schur(a,t,z,status,selector=positive_imag,sdim=selected)
            else
                call schur(a,t,z,status,sort=trim(sortcode),sdim=selected)
            end if
        end if
        print '(2(i0,1x))', status,selected
        do j=1,n
            do i=1,n
                print '(4(es25.16e3,1x))', real(t(i,j),dp),aimag(t(i,j)), &
                    real(z(i,j),dp),aimag(z(i,j))
            end do
        end do
        deallocate(a,t,z)
    end do
contains
    logical function positive_imag(w)
        complex(dp), intent(in) :: w !! eigenvalue to classify
        positive_imag = aimag(w) > 0.0_dp
    end function positive_imag
end program driver_schur
