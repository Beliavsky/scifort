! SPDX-License-Identifier: MIT
! Deterministic rank-1/2/5/15 reference driver for SciPy Gustafsson comparisons.
program driver_signal_filtfilt_nd_gust
    use, intrinsic :: iso_c_binding, only : c_loc, c_f_pointer
    use scifort_kinds, only : dp
    use scifort_signal_filtfilt_nd, only : filtfilt_nd
    implicit none
    integer, parameter :: lengths(5) = [1,2,5,17,47]
    real(dp), allocatable, target :: xr1(:), yr1(:), xr2(:,:), yr2(:,:)
    real(dp), allocatable, target :: xr5(:,:,:,:,:), yr5(:,:,:,:,:)
    real(dp), allocatable, target :: xr15(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:), &
        yr15(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
    complex(dp), allocatable, target :: xc1(:), yc1(:), xc2(:,:), yc2(:,:)
    complex(dp), allocatable, target :: xc5(:,:,:,:,:), yc5(:,:,:,:,:)
    complex(dp), allocatable, target :: xc15(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:), &
        yc15(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
    integer :: cid, ios, rk, typ, coeff, ircase, n, axis, status, irlen
    real(dp) :: b(4), a(4)
    complex(dp) :: cb(4), ca(4)
    integer :: nb, na

    do
        read(*,*,iostat=ios) cid
        if (ios /= 0) exit
        rk = mod(cid,6)
        typ = mod(cid/6,8)
        coeff = mod(cid/48,4)
        ircase = mod(cid/192,5)
        n = lengths(mod(cid/6 + cid/48 + 7*(cid/192),5)+1)
        call set_coefficients(coeff, b, a, cb, ca, nb, na)
        irlen = 10000
        select case (ircase)
        case (1)
            irlen = 1
        case (2)
            irlen = 2
        case (3)
            irlen = 7
        case (4)
            irlen = 28
        end select
        if (allocated(xr1)) deallocate(xr1,yr1,xc1,yc1)
        if (allocated(xr2)) deallocate(xr2,yr2,xc2,yc2)
        if (allocated(xr5)) deallocate(xr5,yr5,xc5,yc5)
        if (allocated(xr15)) deallocate(xr15,yr15,xc15,yc15)
        select case (rk)
        case (0)
            axis = 1
            allocate(xr1(n),yr1(n),xc1(n),yc1(n))
            call fill_input(xr1,xc1,cid)
            call run_case(xr1,xc1,yr1,yc1,axis,typ,irlen,status)
        case (1)
            axis = 2
            allocate(xr2(2,n),yr2(2,n),xc2(2,n),yc2(2,n))
            call fill_input(xr2,xc2,cid)
            call run_case(xr2,xc2,yr2,yc2,axis,typ,irlen,status)
        case (2)
            axis = 1
            allocate(xr2(n,2),yr2(n,2),xc2(n,2),yc2(n,2))
            call fill_input(xr2,xc2,cid)
            call run_case(xr2,xc2,yr2,yc2,axis,typ,irlen,status)
        case (3)
            axis = 3
            allocate(xr5(2,2,n,2,1),yr5(2,2,n,2,1), &
                xc5(2,2,n,2,1),yc5(2,2,n,2,1))
            call fill_input(xr5,xc5,cid)
            call run_case(xr5,xc5,yr5,yc5,axis,typ,irlen,status)
        case (4)
            axis = 1
            allocate(xr15(n,1,2,1,1,1,1,1,1,1,1,1,1,1,2), &
                yr15(n,1,2,1,1,1,1,1,1,1,1,1,1,1,2), &
                xc15(n,1,2,1,1,1,1,1,1,1,1,1,1,1,2), &
                yc15(n,1,2,1,1,1,1,1,1,1,1,1,1,1,2))
            call fill_input(xr15,xc15,cid)
            call run_case(xr15,xc15,yr15,yc15,axis,typ,irlen,status)
        case (5)
            axis = 15
            allocate(xr15(2,1,2,1,1,1,1,1,1,1,1,1,1,1,n), &
                yr15(2,1,2,1,1,1,1,1,1,1,1,1,1,1,n), &
                xc15(2,1,2,1,1,1,1,1,1,1,1,1,1,1,n), &
                yc15(2,1,2,1,1,1,1,1,1,1,1,1,1,1,n))
            call fill_input(xr15,xc15,cid)
            call run_case(xr15,xc15,yr15,yc15,axis,typ,irlen,status)
        end select
        if (status /= 0) error stop 'unexpected filtfilt_nd Gustafsson failure'
    end do
contains
    subroutine set_coefficients(choice,br,ar,bc,ac,nbr,nar)
        integer, intent(in) :: choice !! filter coefficient family
        real(dp), intent(out) :: br(4),ar(4) !! real coefficient tables
        complex(dp), intent(out) :: bc(4),ac(4) !! complex coefficients
        integer, intent(out) :: nbr,nar !! active lengths
        br = 0.0_dp
        ar = 0.0_dp
        select case (choice)
        case (0)
            nbr=1
            nar=1
            br(1)=0.37_dp
            ar(1)=1.2_dp
        case (1)
            nbr=3
            nar=3
            br(1:3) = [0.31_dp,0.12_dp,0.03_dp]
            ar(1:3) = [1.0_dp,-0.23_dp,0.07_dp]
        case (2)
            nbr=2
            nar=4
            br(1:2) = [0.25_dp,-0.11_dp]
            ar(1:4) = [1.1_dp,-0.29_dp,0.07_dp,-0.014_dp]
        case (3)
            nbr=4
            nar=2
            br(1:4) = [0.17_dp,0.07_dp,0.05_dp,-0.008_dp]
            ar(1:2) = [0.9_dp,-0.21_dp]
        end select
        bc = cmplx(br,0.0_dp,dp)
        ac = cmplx(ar,0.0_dp,dp)
        if (nbr > 1) bc(2) = bc(2)+(0.0_dp,0.032_dp)
        if (nar > 1) ac(2) = ac(2)-(0.0_dp,0.019_dp)
        if (nar == 1) ac(1) = ac(1)+(0.0_dp,0.07_dp)
    end subroutine set_coefficients

    subroutine fill_input(xr,xc,case_id)
        real(dp), contiguous, target, intent(out) :: xr(..) !! samples
        complex(dp), contiguous, target, intent(out) :: xc(..) !! complex samples
        integer, intent(in) :: case_id !! case identifier
        real(dp), pointer :: rr(:)
        complex(dp), pointer :: cc(:)
        integer :: k
        call c_f_pointer(c_loc(xr),rr,[size(xr)])
        call c_f_pointer(c_loc(xc),cc,[size(xc)])
        do k=1,size(xr)
            rr(k) = sin(0.087_dp*k+0.03_dp*case_id)+0.003_dp*k
            cc(k) = cmplx(rr(k),0.1_dp*cos(0.12_dp*k+0.013_dp*case_id),dp)
        end do
    end subroutine fill_input

    subroutine run_case(xr,xc,yr,yc,axis,typ,irlen,code)
        real(dp), contiguous, target, intent(in) :: xr(..) !! samples
        complex(dp), contiguous, target, intent(in) :: xc(..) !! samples
        real(dp), contiguous, target, intent(out) :: yr(..) !! real output
        complex(dp), contiguous, target, intent(out) :: yc(..) !! complex output
        integer, intent(in) :: axis,typ,irlen !! axis, type, impulse response length
        integer, intent(out) :: code !! filter result
        logical :: real_output
        real(dp), pointer :: rr(:)
        complex(dp), pointer :: cc(:)
        integer :: k
        real_output = typ == 0
        select case (typ)
        case (0)
            call filtfilt_nd(b(:nb),a(:na),xr,yr,axis=axis,method='gust',irlen=irlen,status=code)
        case (1)
            call filtfilt_nd(cb(:nb),ca(:na),xc,yc,axis=axis,method='gust',irlen=irlen,status=code)
        case (2)
            call filtfilt_nd(b(:nb),a(:na),xc,yc,axis=axis,method='gust',irlen=irlen,status=code)
        case (3)
            call filtfilt_nd(b(:nb),ca(:na),xr,yc,axis=axis,method='gust',irlen=irlen,status=code)
        case (4)
            call filtfilt_nd(b(:nb),ca(:na),xc,yc,axis=axis,method='gust',irlen=irlen,status=code)
        case (5)
            call filtfilt_nd(cb(:nb),a(:na),xr,yc,axis=axis,method='gust',irlen=irlen,status=code)
        case (6)
            call filtfilt_nd(cb(:nb),a(:na),xc,yc,axis=axis,method='gust',irlen=irlen,status=code)
        case (7)
            call filtfilt_nd(cb(:nb),ca(:na),xr,yc,axis=axis,method='gust',irlen=irlen,status=code)
        end select
        write(*,'(a,1x,i0,1x,i0)') 'CASE',size(xr),code
        if (real_output) then
            call c_f_pointer(c_loc(yr),rr,[size(yr)])
            do k=1,size(rr)
                write(*,'(2(es28.18e3,1x))') rr(k),0.0_dp
            end do
        else
            call c_f_pointer(c_loc(yc),cc,[size(yc)])
            do k=1,size(cc)
                write(*,'(2(es28.18e3,1x))') real(cc(k),dp),aimag(cc(k))
            end do
        end if
    end subroutine run_case
end program driver_signal_filtfilt_nd_gust
