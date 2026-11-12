! SPDX-License-Identifier: MIT
! Emit deterministic test cases for independent SciPy filtfilt_nd comparisons.
program filtfilt_nd_driver
    use, intrinsic :: iso_c_binding, only : c_loc, c_f_pointer
    use scifort_kinds, only : dp
    use scifort_signal_filtfilt_nd, only : filtfilt_nd, sosfiltfilt_nd
    implicit none
    real(dp) :: b(3), a(3), sos(2,6)
    complex(dp) :: cb(3), ca(3), csos(2,6)
    real(dp) :: xr2(12,23), yr2(12,23), xr5(4,2,23,3,2), yr5(4,2,23,3,2)
    real(dp) :: xr15(24,1,2,1,1,1,1,1,1,1,1,1,1,1,2)
    real(dp) :: yr15(24,1,2,1,1,1,1,1,1,1,1,1,1,1,2)
    complex(dp) :: xc2(12,23), yc2(12,23), xc5(4,2,23,3,2), yc5(4,2,23,3,2)
    complex(dp) :: xc15(24,1,2,1,1,1,1,1,1,1,1,1,1,1,2)
    complex(dp) :: yc15(24,1,2,1,1,1,1,1,1,1,1,1,1,1,2)
    integer :: cid, rk, typ, meth, mode, variant, plen, ax, ios, status
    logical :: default_pad
    character(len=8) :: pad
    character(len=8), parameter :: modes(4) = [character(len=8) :: &
        'odd', 'even', 'constant', 'none']
    b = [0.20_dp,0.12_dp,0.08_dp]
    a = [1.0_dp,-0.3_dp,0.05_dp]
    sos(1,:) = [0.30_dp,0.20_dp,0.10_dp,1.0_dp,-0.27_dp,0.04_dp]
    sos(2,:) = [0.90_dp,-0.12_dp,0.08_dp,1.0_dp,-0.18_dp,0.09_dp]
    cb = cmplx(b,0.0_dp,dp)
    ca = cmplx(a,0.0_dp,dp)
    cb(2) = cb(2)+(0.0_dp,0.06_dp)
    ca(2) = ca(2)-(0.0_dp,0.03_dp)
    csos = cmplx(sos,0.0_dp,dp)
    csos(1,2) = csos(1,2)+(0.0_dp,0.04_dp)
    csos(2,5) = csos(2,5)-(0.0_dp,0.04_dp)
    call initialize2()
    call initialize5()
    call initialize15()
    do
        read(*,*,iostat=ios) cid
        if (ios /= 0) exit
        default_pad = cid >= 768
        if (default_pad) cid = cid-768
        rk = mod(cid,3)
        typ = mod(cid/3,8)
        meth = mod(cid/24,2)
        mode = mod(cid/48,4)
        variant = mod(cid/192,2)
        plen = mod(cid/384,2)
        if (default_pad) then
            mode = mod(cid/48,3)
            variant = 0
        end if
        pad = modes(mode+1)
        select case (rk)
        case (0)
            ax = 2
            if (variant == 1) ax = 1
            if (default_pad) then
                call run_case(xr2,xc2,yr2,yc2,ax,typ,meth,pad,stat=status)
            else
                call run_case(xr2,xc2,yr2,yc2,ax,typ,meth,pad,plen,status)
            end if
        case (1)
            ax = 3
            if (variant == 1) ax = 1
            if (default_pad) then
                call run_case(xr5,xc5,yr5,yc5,ax,typ,meth,pad,stat=status)
            else
                call run_case(xr5,xc5,yr5,yc5,ax,typ,meth,pad,plen,status)
            end if
        case (2)
            ax = 1
            if (variant == 1) ax = 15
            if (default_pad) then
                call run_case(xr15,xc15,yr15,yc15,ax,typ,meth,pad,stat=status)
            else
                call run_case(xr15,xc15,yr15,yc15,ax,typ,meth,pad,plen,status)
            end if
        end select
        if (status /= 0) error stop 'unexpected SciFort status'
    end do
contains
    subroutine initialize2()
        integer :: k
        xr2 = reshape([(sin(0.117_dp*k)+0.0003_dp*k,k=1,size(xr2))],shape(xr2))
        xc2 = cmplx(xr2,reshape([(0.16_dp*cos(0.092_dp*k)-0.004_dp, &
            k=1,size(xr2))],shape(xr2)),dp)
    end subroutine initialize2
    subroutine initialize5()
        integer :: k
        xr5 = reshape([(sin(0.117_dp*k)+0.0003_dp*k,k=1,size(xr5))],shape(xr5))
        xc5 = cmplx(xr5,reshape([(0.16_dp*cos(0.092_dp*k)-0.004_dp, &
            k=1,size(xr5))],shape(xr5)),dp)
    end subroutine initialize5
    subroutine initialize15()
        integer :: k
        xr15 = reshape([(sin(0.117_dp*k)+0.0003_dp*k,k=1,size(xr15))],shape(xr15))
        xc15 = cmplx(xr15,reshape([(0.16_dp*cos(0.092_dp*k)-0.004_dp, &
            k=1,size(xr15))],shape(xr15)),dp)
    end subroutine initialize15
    subroutine run_case(xr,xc,yr,yc,ax,typ,meth,pad,plen,stat)
        real(dp), contiguous, target, intent(in) :: xr(..) !! real samples
        complex(dp), contiguous, target, intent(in) :: xc(..) !! complex samples
        real(dp), contiguous, target, intent(out) :: yr(..) !! real result
        complex(dp), contiguous, target, intent(out) :: yc(..) !! complex result
        integer, intent(in) :: ax !! axis index
        integer, intent(in) :: typ !! real/complex combination
        integer, intent(in) :: meth !! 0 for TF, 1 for SOS
        character(len=*), intent(in) :: pad !! pad mode
        integer, optional, intent(in) :: plen !! pad length, absent for SciPy default
        integer, intent(out) :: stat !! status code
        real(dp), pointer :: re(:)
        complex(dp), pointer :: cm(:)
        integer :: k
        logical :: real_result
        real_result = typ == 0
        if (meth == 0) then
            select case (typ)
            case (0)
                call filtfilt_nd(b,a,xr,yr,padtype=trim(pad),padlen=plen,status=stat,axis=ax)
            case (1)
                call filtfilt_nd(cb,ca,xc,yc,padtype=trim(pad),padlen=plen,status=stat,axis=ax)
            case (2)
                call filtfilt_nd(b,a,xc,yc,padtype=trim(pad),padlen=plen,status=stat,axis=ax)
            case (3)
                call filtfilt_nd(b,ca,xr,yc,padtype=trim(pad),padlen=plen,status=stat,axis=ax)
            case (4)
                call filtfilt_nd(b,ca,xc,yc,padtype=trim(pad),padlen=plen,status=stat,axis=ax)
            case (5)
                call filtfilt_nd(cb,a,xr,yc,padtype=trim(pad),padlen=plen,status=stat,axis=ax)
            case (6)
                call filtfilt_nd(cb,a,xc,yc,padtype=trim(pad),padlen=plen,status=stat,axis=ax)
            case (7)
                call filtfilt_nd(cb,ca,xr,yc,padtype=trim(pad),padlen=plen,status=stat,axis=ax)
            end select
        else
            select case (typ)
            case (0)
                call sosfiltfilt_nd(sos,xr,yr,padtype=trim(pad),padlen=plen,status=stat,axis=ax)
            case (1)
                call sosfiltfilt_nd(csos,xc,yc,padtype=trim(pad),padlen=plen,status=stat,axis=ax)
            case (2)
                call sosfiltfilt_nd(sos,xc,yc,padtype=trim(pad),padlen=plen,status=stat,axis=ax)
            case (3)
                call sosfiltfilt_nd(csos,xr,yc,padtype=trim(pad),padlen=plen,status=stat,axis=ax)
            end select
        end if
        write(*,'(a,1x,i0,1x,i0)') 'CASE', size(xr), stat
        if (real_result) then
            call c_f_pointer(c_loc(yr),re,[size(yr)])
            do k = 1,size(xr)
                write(*,'(es25.16,1x,es25.16)') re(k),0.0_dp
            end do
        else
            call c_f_pointer(c_loc(yc),cm,[size(yc)])
            do k = 1,size(xr)
                write(*,'(es25.16,1x,es25.16)') real(cm(k),dp),aimag(cm(k))
            end do
        end if
    end subroutine run_case
end program filtfilt_nd_driver
