! SPDX-License-Identifier: MIT
! Machine-readable driver for SciPy analog frequency-transform comparisons.
program analog_frequency_transform_driver
    use scifort_kinds, only : dp
    use scifort_signal_frequency_transforms, only : lp2lp, lp2hp, lp2bp, lp2bs, &
        lp2lp_zpk, lp2hp_zpk, lp2bp_zpk, lp2bs_zpk, &
        analog_tf_transform_result, analog_zpk_transform_result
    implicit none
    integer :: cases, q, form, kind, n1, n2, j
    real(dp) :: wo, bw, re, im
    complex(dp) :: k
    complex(dp), allocatable :: x(:), y(:)
    type(analog_tf_transform_result) :: ba
    type(analog_zpk_transform_result) :: zpk

    read(*,*) cases
    do q = 1, cases
        read(*,*) form,kind,n1,n2,wo,bw
        allocate(x(n1),y(n2))
        do j = 1, n1
            read(*,*) re, im
            x(j) = cmplx(re,im,dp)
        end do
        do j = 1, n2
            read(*,*) re, im
            y(j) = cmplx(re,im,dp)
        end do
        if (form == 1) then
            select case(kind)
            case(1)
                ba = lp2lp(x,y,wo)
            case(2)
                ba = lp2hp(x,y,wo)
            case(3)
                ba = lp2bp(x,y,wo,bw)
            case(4)
                ba = lp2bs(x,y,wo,bw)
            end select
            if (ba%status == 0) then
                write(*, '(a,4(1x,i0))') 'CASE',q,ba%status,size(ba%b),size(ba%a)
                do j = 1, size(ba%b)
                    call show_complex(ba%b(j))
                end do
                do j = 1, size(ba%a)
                    call show_complex(ba%a(j))
                end do
            else
                write(*, '(a,4(1x,i0))') 'CASE',q,ba%status,0,0
            end if
        else
            read(*,*) re, im
            k = cmplx(re,im,dp)
            select case(kind)
            case(1)
                zpk = lp2lp_zpk(x,y,k,wo)
            case(2)
                zpk = lp2hp_zpk(x,y,k,wo)
            case(3)
                zpk = lp2bp_zpk(x,y,k,wo,bw)
            case(4)
                zpk = lp2bs_zpk(x,y,k,wo,bw)
            end select
            if (zpk%status == 0) then
                write(*, '(a,4(1x,i0))') 'CASE',q,zpk%status,size(zpk%zeros),size(zpk%poles)
                do j = 1, size(zpk%zeros)
                    call show_complex(zpk%zeros(j))
                end do
                do j = 1, size(zpk%poles)
                    call show_complex(zpk%poles(j))
                end do
                call show_complex(zpk%gain)
            else
                write(*, '(a,4(1x,i0))') 'CASE',q,zpk%status,0,0
            end if
        end if
        deallocate(x,y)
    end do
contains
    subroutine show_complex(value)
        complex(dp), intent(in) :: value !! complex quantity to print
        write(*, '(es26.17e3,1x,es26.17e3)') real(value,dp),aimag(value)
    end subroutine show_complex
end program analog_frequency_transform_driver
