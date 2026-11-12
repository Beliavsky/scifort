! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! CSV-friendly deterministic batches for comparison with scipy.signal.czt
! and scipy.signal.zoom_fft. Input rows: kind,n,m,complex_flag,variant.
program signal_czt_nd_reference_driver
    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, iostat_end
    use scifort_kinds, only : dp
    use scifort_signal_czt_nd, only : czt_into, zoom_fft_into
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), allocatable :: x(:,:,:)
    complex(dp), allocatable :: z(:,:,:), y(:,:,:)
    complex(dp) :: w,a
    integer :: kind, n, m, is_complex, variant, stat, ios
    integer :: i, j, k

    do
        read(input_unit,*,iostat=ios) kind,n,m,is_complex,variant
        if (ios == iostat_end) exit
        if (ios /= 0) error stop 'malformed reference case'
        if (n < 1 .or. m < 1 .or. n > 100 .or. m > 150) error stop 'reference case bounds'
        allocate(x(2,n,3),z(2,n,3),y(2,m,3))
        do k=1,3
            do j=1,n
                do i=1,2
                    x(i,j,k) = cos(0.21_dp*real(i+3*j+7*k,dp)) + &
                        0.11_dp*real(i-j+2*k,dp)
                    z(i,j,k) = cmplx(x(i,j,k), &
                        0.35_dp*sin(0.11_dp*real(2*i+j+5*k,dp)),dp)
                end do
            end do
        end do
        w = cmplx(cos(-2.0_dp*pi/real(m+4,dp)), &
            sin(-2.0_dp*pi/real(m+4,dp)),dp)
        a = cmplx(cos(0.17_dp),sin(0.17_dp),dp)
        select case (kind)
        case (0)
            select case (variant)
            case (0)
                if (is_complex == 1) then
                    call czt_into(z,y,m=m,axis=2,status=stat)
                else
                    call czt_into(x,y,m=m,axis=2,status=stat)
                end if
            case (1)
                if (is_complex == 1) then
                    call czt_into(z,y,m=m,w=w,a=a,axis=-2,status=stat)
                else
                    call czt_into(x,y,m=m,w=w,a=a,axis=-2,status=stat)
                end if
            case (2)
                w = cmplx(0.985_dp,-0.009_dp,dp)
                a = cmplx(1.016_dp,0.035_dp,dp)
                if (is_complex == 1) then
                    call czt_into(z,y,m=m,w=w,a=a,axis=2,status=stat)
                else
                    call czt_into(x,y,m=m,w=w,a=a,axis=2,status=stat)
                end if
            case default
                error stop 'unknown CZT variant'
            end select
        case (1)
            select case (variant)
            case (0)
                if (is_complex == 1) then
                    call zoom_fft_into(z,y,[-1.2_dp,3.1_dp],m=m,fs=12.0_dp,axis=2,status=stat)
                else
                    call zoom_fft_into(x,y,[-1.2_dp,3.1_dp],m=m,fs=12.0_dp,axis=2,status=stat)
                end if
            case (1)
                if (is_complex == 1) then
                    call zoom_fft_into(z,y,[-0.8_dp,1.75_dp],m=m,fs=7.0_dp, &
                        endpoint=.true.,axis=2,status=stat)
                else
                    call zoom_fft_into(x,y,[-0.8_dp,1.75_dp],m=m,fs=7.0_dp, &
                        endpoint=.true.,axis=2,status=stat)
                end if
            case (2)
                if (is_complex == 1) then
                    call zoom_fft_into(z,y,2.8_dp,m=m,fs=14.0_dp,axis=2,status=stat)
                else
                    call zoom_fft_into(x,y,2.8_dp,m=m,fs=14.0_dp,axis=2,status=stat)
                end if
            case default
                error stop 'unknown zoom variant'
            end select
        case default
            error stop 'unknown transform'
        end select
        write(output_unit,'(a,6(1x,i0))') 'CASE',kind,n,m,is_complex,variant,stat
        if (stat /= 0) error stop 'reference transform failed'
        do k=1,3
            do j=1,m
                do i=1,2
                    write(output_unit,'(2(es26.17,1x))') real(y(i,j,k),dp),aimag(y(i,j,k))
                end do
            end do
        end do
        deallocate(x,z,y)
    end do
end program signal_czt_nd_reference_driver
