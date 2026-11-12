! SPDX-License-Identifier: MIT
! Frequency-domain filters: deterministic complex reference arrays.
program fourier_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : fourier_gaussian_nd, fourier_uniform_nd, fourier_shift_nd
    implicit none
    complex(dp) :: a(5,6),b(5,6),c(4,3,2),d(4,3,2)
    complex(dp) :: r(5,4),u(5,4),t(2,3,5),v(2,3,5)
    integer :: i,j,k,unit,s,op
    open(newunit=unit,file='fourier_results.bin',access='stream',form='unformatted',status='replace')
    do j=1,6
        do i=1,5
            a(i,j)=cmplx(sin(0.4_dp*i+0.23_dp*j),cos(0.27_dp*i-0.49_dp*j),dp)
        end do
    end do
    do op=1,3
        select case(op)
        case(1)
            call fourier_gaussian_nd(a,[1.2_dp,0.8_dp],b,s)
        case(2)
            call fourier_uniform_nd(a,[2.1_dp,3.3_dp],b,s)
        case(3)
            call fourier_shift_nd(a,[0.75_dp,-1.35_dp],b,s)
        end select
        call check(s)
        write(unit) b
        select case(op)
        case(1)
            call fourier_gaussian_nd(a,1.1_dp,b,s)
        case(2)
            call fourier_uniform_nd(a,2.4_dp,b,s)
        case(3)
            call fourier_shift_nd(a,0.3_dp,b,s)
        end select
        call check(s)
        write(unit) b
    end do
    do k=1,2
        do j=1,3
            do i=1,4
                c(i,j,k)=cmplx(sin(.4_dp*i+.7_dp*j-.2_dp*k),sin(.31_dp*i*j*k),dp)
            end do
        end do
    end do
    call fourier_gaussian_nd(c,[0.8_dp,1.9_dp,0.35_dp],d,s)
    call check(s)
    write(unit) d
    call fourier_uniform_nd(c,[3.2_dp,1.6_dp,2.7_dp],d,s)
    call check(s)
    write(unit) d
    call fourier_shift_nd(c,[0.4_dp,-0.75_dp,1.9_dp],d,s)
    call check(s)
    write(unit) d
    do j=1,4
        do i=1,5
            r(i,j)=cmplx(sin(.3_dp*i+j),cos(.44_dp*i-j),dp)
        end do
    end do
    call fourier_gaussian_nd(r,[1.2_dp,0.7_dp],u,s,n=8,axis=1)
    call check(s)
    write(unit) u
    call fourier_uniform_nd(r,[2.3_dp,3.1_dp],u,s,n=8,axis=1)
    call check(s)
    write(unit) u
    call fourier_shift_nd(r,[0.9_dp,-1.4_dp],u,s,n=8,axis=1)
    call check(s)
    write(unit) u
    do k=1,5
        do j=1,3
            do i=1,2
                t(i,j,k)=cmplx(sin(.5_dp*i+.6_dp*j+.3_dp*k),cos(.41_dp*i-.53_dp*j*k),dp)
            end do
        end do
    end do
    call fourier_gaussian_nd(t,[1.6_dp,0.4_dp,2.1_dp],v,s,n=9,axis=-1)
    call check(s)
    write(unit) v
    call fourier_uniform_nd(t,[2.1_dp,3.1_dp,0.8_dp],v,s,n=9,axis=-1)
    call check(s)
    write(unit) v
    call fourier_shift_nd(t,[1.25_dp,-0.4_dp,2.2_dp],v,s,n=9,axis=-1)
    call check(s)
    write(unit) v
    close(unit)
    print *, 'PASS: Fourier operation driver'
contains
    subroutine check(status)
        integer, intent(in) :: status
        if (status/=0) then
            print *, 'FAIL Fourier ',status
            stop 1
        end if
    end subroutine check
end program fourier_driver
