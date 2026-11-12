! SPDX-License-Identifier: MIT
! Batch probe used by compare_signal_partial_fraction_scipy.py.
program probe_partial_fraction
    use scifort_signal, only : dp, residue, residuez, invres, invresz, &
        partial_fraction_result, partial_fraction_tf_result, partial_fraction_success
    implicit none
    complex(dp), allocatable :: b(:), a(:)
    type(partial_fraction_result) :: pf
    type(partial_fraction_tf_result) :: tf
    real(dp) :: re, im
    integer :: cases, run, kind, nb, na, i

    read(*,*) cases
    do run = 1, cases
        read(*,*) kind, nb, na
        allocate(b(nb),a(na))
        do i = 1, nb
            read(*,*) re, im
            b(i) = cmplx(re,im,dp)
        end do
        do i = 1, na
            read(*,*) re, im
            a(i) = cmplx(re,im,dp)
        end do
        if (kind == 1) then
            pf = residue(b,a)
        else
            pf = residuez(b,a)
        end if
        if (pf%status /= partial_fraction_success) then
            print '(4(i0,1x))', pf%status,0,0,0
        else
            print '(4(i0,1x))', pf%status,size(pf%residues),size(pf%direct),size(pf%poles)
            do i = 1, size(pf%residues)
                call emit(pf%residues(i))
            end do
            do i = 1, size(pf%poles)
                call emit(pf%poles(i))
            end do
            do i = 1, size(pf%direct)
                call emit(pf%direct(i))
            end do
            if (kind == 1) then
                tf = invres(pf%residues,pf%poles,pf%direct)
            else
                tf = invresz(pf%residues,pf%poles,pf%direct)
            end if
            print '(3(i0,1x))', tf%status,size(tf%b),size(tf%a)
            do i = 1, size(tf%b)
                call emit(tf%b(i))
            end do
            do i = 1, size(tf%a)
                call emit(tf%a(i))
            end do
        end if
        deallocate(b,a)
    end do
contains
    subroutine emit(z)
        complex(dp), intent(in) :: z !! complex value to print
        print '(es25.16e3,1x,es25.16e3)', real(z,dp),aimag(z)
    end subroutine emit
end program probe_partial_fraction
