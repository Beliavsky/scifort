! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Emit deterministic multidimensional multirate samples for SciPy reference comparisons.
program multirate_nd_driver
    use scifort_kinds, only : dp
    use scifort_signal_resample, only : upfirdn, resample_poly, decimate
    use scifort_signal_fourier_resample, only : resample
    implicit none
    real(dp) :: r2(12,4), r3(5,12,3), taps(3)
    complex(dp) :: c2(12,4), c3(5,12,3), ctaps(3)
    integer :: a, i, j, k
    interface emit
        procedure :: emit_real2, emit_complex2, emit_real3, emit_complex3
    end interface

    taps=[0.2_dp,0.6_dp,0.2_dp]
    ctaps=cmplx(taps,[0.1_dp,-0.2_dp,0.3_dp],dp)
    do j=1,4
        do i=1,12
            r2(i,j)=sin(0.17_dp*real(i*j,dp))+0.03_dp*real(i+j,dp)
            c2(i,j)=cmplx(r2(i,j),cos(0.2_dp*real(i+2*j,dp)),dp)
        end do
    end do
    do k=1,3
        do j=1,12
            do i=1,5
                r3(i,j,k)=sin(real(i+j+k,dp)*0.12_dp)
                c3(i,j,k)=cmplx(r3(i,j,k),cos(real(i+2*j+3*k,dp)*0.13_dp),dp)
            end do
        end do
    end do

    do a=1,2
        call emit('fourier_r2_a'//num(a),resample(r2,9,axis=a))
        call emit('fourier_c2_a'//num(a),resample(c2,7,window='hann',axis=a))
        call emit('poly_r2_a'//num(a),resample_poly(r2,2,3,taps=taps,axis=a))
        call emit('poly_c2_a'//num(a),resample_poly(c2,3,2,taps=taps,padtype='mean',axis=a))
        call emit('up_rr2_a'//num(a),upfirdn(taps,r2,up=2,down=3,axis=a,mode='edge'))
        call emit('up_rc2_a'//num(a),upfirdn(taps,c2,up=2,down=3,axis=a,mode='symmetric'))
        call emit('up_cr2_a'//num(a),upfirdn(ctaps,r2,up=2,down=3,axis=a,mode='reflect'))
        call emit('up_cc2_a'//num(a),upfirdn(ctaps,c2,up=2,down=3,axis=a))
        call emit('dec_r2_a'//num(a),decimate(r2,2,n=12,ftype='fir',axis=a))
        call emit('dec_c2_a'//num(a),decimate(c2,2,n=12,ftype='fir',zero_phase=.false.,axis=a))
    end do
    call emit('poly_r2_edge',resample_poly(r2,2,3,taps=taps,padtype='edge',axis=2))
    call emit('poly_c2_const',resample_poly(c2,2,3,taps=taps,padtype='constant', &
        cval=cmplx(0.75_dp,-0.25_dp,dp),axis=2))
    call emit('dec_r2_iir',decimate(r2,2,n=4,ftype='iir',zero_phase=.false.,axis=1))
    do a=1,3
        call emit('fourier_r3_a'//num(a),resample(r3,7,window='hann',axis=a))
        call emit('fourier_c3_a'//num(a),resample(c3,6,axis=a))
        call emit('poly_r3_a'//num(a),resample_poly(r3,2,3,axis=a))
        call emit('poly_c3_a'//num(a),resample_poly(c3,3,2,taps=taps,axis=a))
        call emit('up_rr3_a'//num(a),upfirdn(taps,r3,up=2,down=3,axis=a,mode='edge'))
        call emit('up_rc3_a'//num(a),upfirdn(taps,c3,up=2,down=3,axis=a,mode='symmetric'))
        call emit('up_cr3_a'//num(a),upfirdn(ctaps,r3,up=2,down=3,axis=a,mode='reflect'))
        call emit('up_cc3_a'//num(a),upfirdn(ctaps,c3,up=2,down=3,axis=a))
        call emit('dec_r3_a'//num(a),decimate(r3,2,n=12,ftype='fir',axis=a))
        call emit('dec_c3_a'//num(a),decimate(c3,2,n=12,ftype='fir',zero_phase=.false.,axis=a))
    end do
contains
    function num(v) result(s)
        integer, intent(in) :: v !! axis number
        character(1) :: s
        write(s,'(i1)') v
    end function num

    subroutine emit_real2(label, x)
        character(*), intent(in) :: label !! case id
        real(dp), intent(in) :: x(:,:) !! values
        integer :: i, j
        print '(a,1x,a,1x,i0,1x,i0,1x,i0)', 'CASE',label,2,size(x,1),size(x,2)
        do j=1,size(x,2)
            do i=1,size(x,1)
                print '(es25.16e3,1x,es25.16e3)', x(i,j), 0.0_dp
            end do
        end do
    end subroutine emit_real2

    subroutine emit_complex2(label, x)
        character(*), intent(in) :: label !! case id
        complex(dp), intent(in) :: x(:,:) !! values
        integer :: i, j
        print '(a,1x,a,1x,i0,1x,i0,1x,i0)', 'CASE',label,2,size(x,1),size(x,2)
        do j=1,size(x,2)
            do i=1,size(x,1)
                print '(es25.16e3,1x,es25.16e3)', real(x(i,j),dp), aimag(x(i,j))
            end do
        end do
    end subroutine emit_complex2

    subroutine emit_real3(label, x)
        character(*), intent(in) :: label !! case id
        real(dp), intent(in) :: x(:,:,:) !! values
        integer :: i, j, k
        print '(a,1x,a,1x,i0,1x,i0,1x,i0,1x,i0)', 'CASE',label,3,size(x,1),size(x,2),size(x,3)
        do k=1,size(x,3)
            do j=1,size(x,2)
                do i=1,size(x,1)
                    print '(es25.16e3,1x,es25.16e3)', x(i,j,k), 0.0_dp
                end do
            end do
        end do
    end subroutine emit_real3

    subroutine emit_complex3(label, x)
        character(*), intent(in) :: label !! case id
        complex(dp), intent(in) :: x(:,:,:) !! values
        integer :: i, j, k
        print '(a,1x,a,1x,i0,1x,i0,1x,i0,1x,i0)', 'CASE',label,3,size(x,1),size(x,2),size(x,3)
        do k=1,size(x,3)
            do j=1,size(x,2)
                do i=1,size(x,1)
                    print '(es25.16e3,1x,es25.16e3)', real(x(i,j,k),dp), aimag(x(i,j,k))
                end do
            end do
        end do
    end subroutine emit_complex3
end program multirate_nd_driver
