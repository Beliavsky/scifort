! SPDX-License-Identifier: MIT
! Standard-input reference driver for bilinear and elliptic prototypes.
program driver_signal_bilinear_ellipap
    use scifort_kinds, only : dp
    use scifort_signal_bilinear, only : bilinear, bilinear_zpk
    use scifort_signal_ellipap, only : ellipap
    use scifort_signal_analog_prototypes, only : analog_prototype_result, besselap
    use scifort_signal_frequency_transforms, only : analog_tf_transform_result, analog_zpk_transform_result
    implicit none
    type(analog_tf_transform_result) :: tf
    type(analog_zpk_transform_result) :: zpk
    type(analog_prototype_result) :: proto
    complex(dp), allocatable :: b(:), a(:), z(:), p(:)
    complex(dp) :: gain
    real(dp) :: rate, rp, rs, x, y
    integer :: nb, na, nz, np, i, status
    character(1) :: op
    character(5) :: norm

    do
        read(*,*,iostat=status) op
        if (status /= 0) exit
        select case (op)
        case ('B')
            read(*,*) nb,na,rate
            allocate(b(nb),a(na))
            do i=1,nb
                read(*,*) x,y
                b(i)=cmplx(x,y,dp)
            end do
            do i=1,na
                read(*,*) x,y
                a(i)=cmplx(x,y,dp)
            end do
            tf=bilinear(b,a,rate)
            if (tf%status == 0) then
                write(*,'(a,3(1x,i0))') 'B',tf%status,size(tf%b),size(tf%a)
                call print_array(tf%b)
                call print_array(tf%a)
            else
                write(*,'(a,3(1x,i0))') 'B',tf%status,0,0
            end if
            deallocate(b,a)
        case ('Z')
            read(*,*) nz,np,rate,x,y
            gain=cmplx(x,y,dp)
            allocate(z(nz),p(np))
            do i=1,nz
                read(*,*) x,y
                z(i)=cmplx(x,y,dp)
            end do
            do i=1,np
                read(*,*) x,y
                p(i)=cmplx(x,y,dp)
            end do
            zpk=bilinear_zpk(z,p,gain,rate)
            if (zpk%status == 0) then
                write(*,'(a,3(1x,i0))') 'Z',zpk%status,size(zpk%zeros),size(zpk%poles)
                write(*,'(2(es26.17,1x))') real(zpk%gain,dp),aimag(zpk%gain)
                call print_array(zpk%zeros)
                call print_array(zpk%poles)
            else
                write(*,'(a,3(1x,i0))') 'Z',zpk%status,0,0
            end if
            deallocate(z,p)
        case ('E')
            read(*,*) nb,rp,rs
            proto=ellipap(nb,rp,rs)
            call print_proto('E',proto)
        case ('S')
            read(*,*) nb,norm
            proto=besselap(nb,trim(norm))
            call print_proto('S',proto)
        end select
    end do
contains
    subroutine print_array(arr)
        complex(dp), intent(in) :: arr(:) !! complex output data
        integer :: j
        do j=1,size(arr)
            write(*,'(2(es26.17,1x))') real(arr(j),dp),aimag(arr(j))
        end do
    end subroutine print_array

    subroutine print_proto(label, res)
        character(1), intent(in) :: label !! operation label
        type(analog_prototype_result), intent(in) :: res !! resulting prototype
        if (res%status == 0) then
            write(*,'(a,3(1x,i0))') label,res%status,size(res%zeros),size(res%poles)
            write(*,'(es26.17)') res%gain
            call print_array(res%zeros)
            call print_array(res%poles)
        else
            write(*,'(a,3(1x,i0))') label,res%status,0,0
        end if
    end subroutine print_proto
end program driver_signal_bilinear_ellipap
