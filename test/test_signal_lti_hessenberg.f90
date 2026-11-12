! SPDX-License-Identifier: MIT
! Compare Hessenberg MIMO frequency evaluation with the direct shifted solve.
program test_signal_lti_hessenberg
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_signal, only : dp, freqresp_mimo, dfreqresp_mimo, &
        bode_mimo, dbode_mimo, mimo_frequency_response, mimo_bode_result, &
        design_status_success, design_status_invalid_input, design_status_singular
    implicit none
    integer, parameter :: ns = 16, ni = 3, no = 2, nw = 48
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: a(ns,ns), b(ns,ni), c(no,ns), d(no,ni), w(nw)
    real(dp) :: at(2,2), bt(2,1), ct(1,2), dt(1,1)
    complex(dp) :: ac(ns,ns), bc(ns,ni), cc(no,ns), dc(no,ni)
    real(dp) :: difference
    integer :: i, j, k
    type(mimo_frequency_response) :: direct, reduced, automatic
    type(mimo_bode_result) :: plot1, plot2

    do j = 1,ns
        do i = 1,ns
            a(i,j) = 0.025_dp*sin(real(i+3*j,dp))
            if (i == j) a(i,j) = a(i,j)-real(i,dp)*0.4_dp-0.4_dp
            ac(i,j) = cmplx(a(i,j),0.03_dp*cos(real(2*i+j,dp)),dp)
        end do
    end do
    do j = 1,ni
        do i = 1,ns
            b(i,j) = sin(real(3*i+j,dp))
            bc(i,j) = cmplx(b(i,j),0.2_dp*cos(real(i+2*j,dp)),dp)
        end do
    end do
    do j = 1,ns
        do i = 1,no
            c(i,j) = cos(real(i+2*j,dp))
            cc(i,j) = cmplx(c(i,j),0.2_dp*sin(real(i+j,dp)),dp)
        end do
    end do
    d = 0.2_dp
    dc = cmplx(d,0.3_dp,dp)
    do k = 1,nw
        w(k) = real(k-1,dp)*pi/real(nw-1,dp)
    end do

    direct = freqresp_mimo(a,b,c,d,w=w,solver='direct')
    reduced = freqresp_mimo(a,b,c,d,w=w,solver='hessenberg')
    automatic = freqresp_mimo(a,b,c,d,w=w)
    call check(direct%status == design_status_success,'real analog direct')
    call check(reduced%status == design_status_success,'real analog Hessenberg')
    call check(automatic%status == design_status_success,'real analog auto')
    difference = maxval(abs(reduced%response-direct%response)/(1.0_dp+abs(direct%response)))
    call check(difference < 2.0e-12_dp,'real analog agreement')
    call check(maxval(abs(automatic%response-reduced%response)) == 0.0_dp,'real auto chooses Hessenberg')

    direct = dfreqresp_mimo(a,b,c,d,w=w,solver='direct')
    reduced = dfreqresp_mimo(a,b,c,d,w=w,solver='hessenberg')
    difference = maxval(abs(reduced%response-direct%response)/(1.0_dp+abs(direct%response)))
    call check(difference < 2.0e-12_dp,'real digital agreement')

    direct = freqresp_mimo(ac,bc,cc,dc,w=w,solver='direct')
    reduced = freqresp_mimo(ac,bc,cc,dc,w=w,solver='hessenberg')
    difference = maxval(abs(reduced%response-direct%response)/(1.0_dp+abs(direct%response)))
    call check(difference < 2.0e-12_dp,'complex analog agreement')
    plot1 = bode_mimo(ac,bc,cc,dc,w=w,solver='direct')
    plot2 = bode_mimo(ac,bc,cc,dc,w=w,solver='hessenberg')
    call check(plot2%status == design_status_success,'complex analog Bode')
    call check(maxval(abs(plot1%magnitude-plot2%magnitude)) < 1.0e-10_dp,'Bode magnitude')
    call check(maxval(abs(plot1%phase-plot2%phase)) < 1.0e-10_dp,'Bode phase')

    direct = dfreqresp_mimo(ac,bc,cc,dc,w=w,solver='direct')
    reduced = dfreqresp_mimo(ac,bc,cc,dc,w=w,solver='hessenberg')
    difference = maxval(abs(reduced%response-direct%response)/(1.0_dp+abs(direct%response)))
    call check(difference < 2.0e-12_dp,'complex digital agreement')
    plot1 = dbode_mimo(ac,bc,cc,dc,w=w,dt=0.1_dp,solver='direct')
    plot2 = dbode_mimo(ac,bc,cc,dc,w=w,dt=0.1_dp,solver='hessenberg')
    call check(maxval(abs(plot1%magnitude-plot2%magnitude)) < 1.0e-10_dp,'digital Bode magnitude')
    call check(maxval(abs(plot1%phase-plot2%phase)) < 1.0e-10_dp,'digital Bode phase')
    call check(maxval(abs(plot2%frequencies-w/0.1_dp)) < 2.0e-14_dp,'digital frequency scaling')

    reduced = freqresp_mimo(a,b,c,d,w=w,solver='unknown')
    call check(reduced%status == design_status_invalid_input,'invalid solver name')
    call check(.not. allocated(reduced%response),'invalid solver no response allocation')
    call check(.not. allocated(reduced%frequencies),'invalid solver no grid allocation')
    reduced = freqresp_mimo(a,b,c,d,w=w(1:3),solver='hessenberg')
    direct = freqresp_mimo(a,b,c,d,w=w(1:3),solver='direct')
    call check(maxval(abs(reduced%response-direct%response)) < 1.0e-11_dp,'forced short grid')

    ! A frequency singularity must invalidate only its own response slice.
    at = 0.0_dp
    bt = 1.0_dp
    ct = 1.0_dp
    dt = 0.0_dp
    reduced = freqresp_mimo(at,bt,ct,dt,w=[0.0_dp,1.0_dp],solver='hessenberg')
    call check(reduced%status == design_status_singular,'singular status')
    call check(ieee_is_nan(real(reduced%response(1,1,1),dp)),'singular value')
    call check(abs(reduced%response(1,1,2)-cmplx(0.0_dp,-2.0_dp,dp)) < 1.0e-12_dp, &
        'nonsingular response retained')

    write(*,'(a)') 'PASS: Hessenberg MIMO reduction and shifted solves'
contains
    subroutine check(pass, description)
        logical, intent(in) :: pass !! test predicate
        character(len=*), intent(in) :: description !! test name
        if (.not. pass) then
            print '(a)', 'FAIL: '//description
            error stop 1
        end if
    end subroutine check
end program test_signal_lti_hessenberg
