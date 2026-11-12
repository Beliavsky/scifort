! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Deterministic complex FIR resampling and boundary extension references.
program multirate_complex_driver
    use scifort_kinds, only : dp
    use scifort_signal_resample, only : upfirdn, resample_poly
    implicit none
    character(16), parameter :: modes(*) = [character(16) :: &
        'constant', 'mean', 'median', 'minimum', 'maximum', 'antireflect', &
        'antisymmetric', 'line', 'smooth', 'wrap', 'reflect', 'symmetric', 'edge']
    character(16), parameter :: short_modes(*) = [character(16) :: &
        'antireflect', 'antisymmetric', 'wrap', 'line']
    real(dp) :: r1(7), r2(7,4), r3(2,7,3), rh(5)
    complex(dp) :: c1(7), c2(7,4), c3(2,7,3), h(5)
    complex(dp), parameter :: outside=cmplx(0.2_dp,-0.12_dp,dp)
    integer :: a,i,j,k,m
    interface emit
        procedure :: emit_r1, emit_c1, emit_r2, emit_c2, emit_r3, emit_c3
    end interface

    rh=[0.08_dp,0.17_dp,0.51_dp,0.17_dp,0.08_dp]
    h=cmplx(rh,[0.04_dp,-0.08_dp,0.13_dp,-0.03_dp,0.06_dp],dp)
    do i=1,7
        r1(i)=sin(0.31_dp*real(i,dp))+0.07_dp*real(i,dp)
        c1(i)=cmplx(r1(i),cos(0.41_dp*real(i,dp)),dp)
        do j=1,4
            r2(i,j)=r1(i)+0.09_dp*real(j,dp)+0.03_dp*real(i*j,dp)
            c2(i,j)=cmplx(r2(i,j),cos(0.19_dp*real(i+3*j,dp)),dp)
        end do
    end do
    do k=1,3
        do j=1,7
            do i=1,2
                r3(i,j,k)=r1(j)+0.11_dp*real(k,dp)+0.04_dp*real(i,dp)
                c3(i,j,k)=cmplx(r3(i,j,k),cos(0.23_dp*real(i+j+2*k,dp)),dp)
            end do
        end do
    end do

    do m=1,size(modes)
        call emit('ctaps_r1_a1_'//trim(modes(m)),resample_poly(r1,3,2,taps=h,padtype=trim(modes(m))))
        call emit('ctaps_c1_a1_'//trim(modes(m)),resample_poly(c1,3,2,taps=h,padtype=trim(modes(m))))
        do a=1,2
            call emit('ctaps_r2_a'//digit(a)//'_'//trim(modes(m)), &
                resample_poly(r2,3,2,taps=h,padtype=trim(modes(m)),axis=a))
            call emit('ctaps_c2_a'//digit(a)//'_'//trim(modes(m)), &
                resample_poly(c2,3,2,taps=h,padtype=trim(modes(m)),axis=a))
        end do
        do a=1,3
            call emit('ctaps_r3_a'//digit(a)//'_'//trim(modes(m)), &
                resample_poly(r3,3,2,taps=h,padtype=trim(modes(m)),axis=a))
            call emit('ctaps_c3_a'//digit(a)//'_'//trim(modes(m)), &
                resample_poly(c3,3,2,taps=h,padtype=trim(modes(m)),axis=a))
        end do
    end do
    do m=2,5
        call emit('realtaps_r1_a1_'//trim(modes(m)),resample_poly(r1,3,2,taps=rh,padtype=trim(modes(m))))
        call emit('realtaps_c1_a1_'//trim(modes(m)),resample_poly(c1,3,2,taps=rh,padtype=trim(modes(m))))
    end do
    do a=1,2
        call emit('outside_r2_a'//digit(a)//'_constant', &
            resample_poly(r2,3,2,taps=h,padtype='constant',cval=outside,axis=a))
        call emit('outside_c2_a'//digit(a)//'_constant', &
            resample_poly(c2,3,2,taps=h,padtype='constant',cval=outside,axis=a))
    end do
    do m=1,size(short_modes)
        call emit('up_rr1_a1_'//trim(short_modes(m)),upfirdn(rh,r1,up=2,down=3,mode=trim(short_modes(m))))
        call emit('up_rc1_a1_'//trim(short_modes(m)),upfirdn(rh,c1,up=2,down=3,mode=trim(short_modes(m))))
        call emit('up_cr1_a1_'//trim(short_modes(m)),upfirdn(h,r1,up=2,down=3,mode=trim(short_modes(m))))
        call emit('up_cc1_a1_'//trim(short_modes(m)),upfirdn(h,c1,up=2,down=3,mode=trim(short_modes(m))))
    end do
contains
    function digit(n) result(s)
        integer, intent(in) :: n !! number in [1,3]
        character(1) :: s
        write(s,'(i1)') n
    end function digit

    subroutine emit_r1(label,x)
        character(*), intent(in) :: label !! test label
        real(dp), intent(in) :: x(:) !! real samples
        integer :: i
        print '(a,1x,a,1x,i0,1x,i0)', 'CASE',label,1,size(x)
        do i=1,size(x)
            print '(es25.16e3,1x,es25.16e3)',x(i),0.0_dp
        end do
    end subroutine emit_r1

    subroutine emit_c1(label,x)
        character(*), intent(in) :: label !! test label
        complex(dp), intent(in) :: x(:) !! complex samples
        integer :: i
        print '(a,1x,a,1x,i0,1x,i0)', 'CASE',label,1,size(x)
        do i=1,size(x)
            print '(es25.16e3,1x,es25.16e3)',real(x(i),dp),aimag(x(i))
        end do
    end subroutine emit_c1

    subroutine emit_r2(label,x)
        character(*), intent(in) :: label !! test label
        real(dp), intent(in) :: x(:,:) !! real sample planes
        integer :: i,j
        print '(a,1x,a,1x,i0,1x,i0,1x,i0)', 'CASE',label,2,shape(x)
        do j=1,size(x,2)
            do i=1,size(x,1)
                print '(es25.16e3,1x,es25.16e3)',x(i,j),0.0_dp
            end do
        end do
    end subroutine emit_r2

    subroutine emit_c2(label,x)
        character(*), intent(in) :: label !! test label
        complex(dp), intent(in) :: x(:,:) !! complex sample planes
        integer :: i,j
        print '(a,1x,a,1x,i0,1x,i0,1x,i0)', 'CASE',label,2,shape(x)
        do j=1,size(x,2)
            do i=1,size(x,1)
                print '(es25.16e3,1x,es25.16e3)',real(x(i,j),dp),aimag(x(i,j))
            end do
        end do
    end subroutine emit_c2

    subroutine emit_r3(label,x)
        character(*), intent(in) :: label !! test label
        real(dp), intent(in) :: x(:,:,:) !! rank-3 real samples
        integer :: i,j,k
        print '(a,1x,a,1x,i0,1x,i0,1x,i0,1x,i0)', 'CASE',label,3,shape(x)
        do k=1,size(x,3)
            do j=1,size(x,2)
                do i=1,size(x,1)
                    print '(es25.16e3,1x,es25.16e3)',x(i,j,k),0.0_dp
                end do
            end do
        end do
    end subroutine emit_r3

    subroutine emit_c3(label,x)
        character(*), intent(in) :: label !! test label
        complex(dp), intent(in) :: x(:,:,:) !! rank-3 complex samples
        integer :: i,j,k
        print '(a,1x,a,1x,i0,1x,i0,1x,i0,1x,i0)', 'CASE',label,3,shape(x)
        do k=1,size(x,3)
            do j=1,size(x,2)
                do i=1,size(x,1)
                    print '(es25.16e3,1x,es25.16e3)',real(x(i,j,k),dp),aimag(x(i,j,k))
                end do
            end do
        end do
    end subroutine emit_c3
end program multirate_complex_driver
