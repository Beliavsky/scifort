! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Array order in this reference driver is Fortran/column-major.
program signal_envelope_nd_driver
    use scifort_kinds, only : dp
    use scifort_signal, only : envelope_into
    implicit none
    real(dp), allocatable :: x(:,:,:,:), y(:,:,:,:), r(:,:,:,:), yf(:), rf(:)
    complex(dp), allocatable :: z(:,:,:,:), cr(:,:,:,:), cf(:)
    integer :: cases, case_number, sh(4), sy(4), ax, nout, low, high
    integer :: is_complex, is_squared, res_mode, st, i, j, k, m, q
    real(dp) :: re, im
    character(len=7) :: mode

    read(*,*) cases
    do case_number=1,cases
        read(*,*) sh, ax, nout, low, high, is_complex, is_squared, res_mode
        sy=sh
        sy(ax)=nout
        allocate(x(sh(1),sh(2),sh(3),sh(4)), z(sh(1),sh(2),sh(3),sh(4)))
        allocate(y(sy(1),sy(2),sy(3),sy(4)), r(sy(1),sy(2),sy(3),sy(4)))
        allocate(cr(sy(1),sy(2),sy(3),sy(4)))
        do m=1,sh(4)
            do k=1,sh(3)
                do j=1,sh(2)
                    do i=1,sh(1)
                        read(*,*) re,im
                        x(i,j,k,m)=re
                        z(i,j,k,m)=cmplx(re,im,dp)
                    end do
                end do
            end do
        end do
        select case(res_mode)
        case(0)
            mode='none'
        case(1)
            mode='lowpass'
        case default
            mode='all'
        end select
        y=-99.0_dp
        r=-99.0_dp
        cr=cmplx(-99.0_dp,-99.0_dp,dp)
        if (is_complex == 0) then
            call envelope_into(x,y,bp_low=low,bp_high=high,n_out=nout, &
                squared=is_squared==1,residual=trim(mode),axis=ax,remainder=r,status=st)
        else
            call envelope_into(z,y,bp_low=low,bp_high=high,n_out=nout, &
                squared=is_squared==1,residual=trim(mode),axis=ax,remainder=cr,status=st)
        end if
        write(*,'(i0,1x,4(i0,1x))') st, sy
        yf=reshape(y,[size(y)])
        rf=reshape(r,[size(r)])
        cf=reshape(cr,[size(cr)])
        do q=1,size(yf)
            if (res_mode == 0) then
                re=0.0_dp
                im=0.0_dp
            else if (is_complex == 0) then
                re=rf(q)
                im=0.0_dp
            else
                re=real(cf(q),dp)
                im=aimag(cf(q))
            end if
            write(*,'(3(es27.18,1x))') yf(q),re,im
        end do
        deallocate(x,y,z,r,cr,yf,rf,cf)
    end do
end program signal_envelope_nd_driver
