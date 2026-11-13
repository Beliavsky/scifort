! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program fftlog_reference_driver
    use scifort_fft
    implicit none
    real(dp), allocatable :: a(:), b(:), x(:,:,:), y(:,:,:)
    real(dp) :: dln, mu, shift, bias, value
    integer :: op, n, dims(3), axis, code, ios, i
    do
        read(*,*,iostat=ios) op,n,dln,mu,shift,bias,dims,axis
        if (ios /= 0) exit
        allocate(a(n))
        if (n > 0) read(*,*) a
        select case(op)
        case(1)
            b = fht(a,dln,mu,shift,bias,status=code)
        case(2)
            b = ifht(a,dln,mu,shift,bias,status=code)
        case(3)
            value = fhtoffset(dln,mu,shift,bias,status=code)
            b = [value]
        case(4,5)
            x = reshape(a,dims)
            allocate(y(dims(1),dims(2),dims(3)))
            if (op == 4) then
                call fht_nd(x,y,dln,mu,shift,bias,axis,status=code)
            else
                call ifht_nd(x,y,dln,mu,shift,bias,axis,status=code)
            end if
            b = reshape(y,[size(y)])
            deallocate(x,y)
        end select
        write(*,'(2(i0,1x))') code,size(b)
        do i = 1, size(b)
            write(*,'(es26.17e3)') b(i)
        end do
        deallocate(a,b)
    end do
end program fftlog_reference_driver
