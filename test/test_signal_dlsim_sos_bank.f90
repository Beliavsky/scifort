! SPDX-License-Identifier: MIT
program test_signal_dlsim_sos_bank
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_discrete_simulation, only : dlsim, dlsim_sos, dlsim_sos_bank, &
        discrete_simulation_result, discrete_sos_bank_result, &
        discrete_sim_success, discrete_sim_invalid_input, discrete_sim_numerical_failure
    implicit none
    real(dp), parameter :: tol = 2.0e-12_dp
    integer, parameter :: ns = 53, nc = 3, sec = 2
    real(dp) :: sos(sec,6), rsos_bank(sec,6,nc), u(ns,nc), zi(sec,2,nc), t(ns)
    real(dp) :: ustrided(2*ns,2*nc), sosstrided(sec,12,2*nc)
    real(dp) :: badsos(sec,6,nc), badu(ns,nc), badzi(sec,2,nc)
    complex(dp) :: csos(sec,6,nc), cu(ns,nc), czi(sec,2,nc)
    complex(dp) :: csos_shared(sec,6)
    type(discrete_sos_bank_result) :: r, r2, first, last
    type(discrete_simulation_result) :: one
    integer :: i, j

    sos(1,:) = [0.18_dp,0.13_dp,0.05_dp,1.0_dp,-0.27_dp,0.14_dp]
    sos(2,:) = [0.79_dp,0.10_dp,0.03_dp,1.0_dp,0.08_dp,0.03_dp]
    do j = 1, nc
        rsos_bank(:,:,j) = sos
        rsos_bank(1,1,j) = sos(1,1) + 0.01_dp*real(j,dp)
        zi(:,:,j) = reshape([(0.01_dp*real(i-j,dp),i=1,2*sec)], [sec,2])
        do i = 1, ns
            u(i,j) = sin(0.17_dp*i+0.35_dp*j)+0.1_dp*cos(0.9_dp*i*j)
            cu(i,j) = cmplx(u(i,j),0.2_dp*sin(0.31_dp*i-0.2_dp*j),dp)
        end do
    end do
    do i = 1, ns
        t(i) = 0.06_dp*real(i-1,dp)*(0.85_dp+0.003_dp*real(i,dp))
    end do
    csos = cmplx(rsos_bank,0.0_dp,dp)
    csos(1,2,:) = csos(1,2,:) + cmplx(0.0_dp,0.03_dp,dp)
    csos(2,5,:) = csos(2,5,:) - cmplx(0.0_dp,0.01_dp,dp)
    csos_shared = csos(:,:,1)
    czi = cmplx(zi,0.015_dp,dp)

    ! Shared SOS: every channel is an independent SISO result.
    r = dlsim_sos_bank(sos,u,0.05_dp,zi=zi)
    call check(r%status == discrete_sim_success,'shared RR status')
    call check(all(shape(r%y)==[ns,nc]),'bank sample shape')
    call check(all(shape(r%x)==[ns,2*sec,nc]),'bank state history shape')
    call check(all(shape(r%zf)==[sec,2,nc]),'bank final state shape')
    do j = 1, nc
        one = dlsim_sos(sos,u(:,j),0.05_dp,zi=zi(:,:,j))
        call match_channel(r,one,j,'shared RR vs SISO')
    end do
    r2 = dlsim(sos,u,0.05_dp,zi=zi)
    call check(maxval(abs(r%y-r2%y))<tol,'dlsim generic bank alias')
    r2 = dlsim_sos(sos,u,0.05_dp,zi=zi)
    call check(maxval(abs(r%zf-r2%zf))<tol,'dlsim_sos generic bank alias')
    r = dlsim_sos_bank(rsos_bank,u,0.05_dp,zi=zi)
    do j = 1, nc
        one = dlsim_sos(rsos_bank(:,:,j),u(:,j),0.05_dp,zi=zi(:,:,j))
        call match_channel(r,one,j,'independent RR vs SISO')
    end do
    r = dlsim_sos_bank(rsos_bank,cu,0.05_dp,zi=czi)
    do j = 1, nc
        one = dlsim_sos(rsos_bank(:,:,j),cu(:,j),0.05_dp,zi=czi(:,:,j))
        call match_channel(r,one,j,'independent RC vs SISO')
    end do
    r = dlsim_sos_bank(csos,u,0.05_dp,zi=czi)
    do j = 1, nc
        one = dlsim_sos(csos(:,:,j),u(:,j),0.05_dp,zi=czi(:,:,j))
        call match_channel(r,one,j,'independent CR vs SISO')
    end do
    r = dlsim_sos_bank(csos,cu,0.05_dp,zi=czi)
    do j = 1, nc
        one = dlsim_sos(csos(:,:,j),cu(:,j),0.05_dp,zi=czi(:,:,j))
        call match_channel(r,one,j,'independent CC vs SISO')
    end do
    r = dlsim_sos_bank(sos,cu,0.05_dp,zi=czi)
    do j = 1, nc
        one = dlsim_sos(sos,cu(:,j),0.05_dp,zi=czi(:,:,j))
        call match_channel(r,one,j,'shared RC vs SISO')
    end do
    r = dlsim_sos_bank(csos_shared,u,0.05_dp,zi=czi)
    do j = 1, nc
        one = dlsim_sos(csos_shared,u(:,j),0.05_dp,zi=czi(:,:,j))
        call match_channel(r,one,j,'shared CR vs SISO')
    end do
    r = dlsim_sos_bank(csos_shared,cu,0.05_dp,zi=czi)
    do j = 1, nc
        one = dlsim_sos(csos_shared,cu(:,j),0.05_dp,zi=czi(:,:,j))
        call match_channel(r,one,j,'shared CC vs SISO')
    end do

    ! Irregular timestamps and initially zero filter states.
    r = dlsim_sos_bank(rsos_bank,cu,0.04_dp,t=t)
    call check(r%status==discrete_sim_success,'timestamp simulation')
    call check(size(r%t)==int(floor(t(ns)/0.04_dp))+1,'timestamp output size')
    call check(all(abs(r%x(1,:,:)) < tol),'default zero states')
    do j = 1, nc
        one = dlsim_sos(rsos_bank(:,:,j),cu(:,j),0.04_dp,t=t)
        call match_channel(r,one,j,'timestamp interpolation vs SISO')
    end do
    r = dlsim_sos_bank(csos,cu,0.05_dp,zi=czi)
    first = dlsim_sos_bank(csos,cu(:24,:),0.05_dp,zi=czi)
    last = dlsim_sos_bank(csos,cu(25:,:),0.05_dp,zi=first%zf)
    call check(first%status==discrete_sim_success .and. last%status==discrete_sim_success, &
        'streaming status')
    call check(maxval(abs(r%y(:24,:)-first%y))<tol,'streaming first block')
    call check(maxval(abs(r%y(25:,:)-last%y))<tol,'streaming second block')
    call check(maxval(abs(r%zf-last%zf))<tol,'streaming final state')

    ! Both arrays may be noncontiguous views.
    ustrided = -999.0_dp
    ustrided(1::2,2::2) = u
    sosstrided = -999.0_dp
    sosstrided(:,1::2,2::2) = rsos_bank
    r = dlsim_sos_bank(sosstrided(:,1::2,2::2),ustrided(1::2,2::2),0.05_dp)
    r2 = dlsim_sos_bank(rsos_bank,u,0.05_dp)
    call check(r%status==discrete_sim_success,'strided arrays status')
    call check(maxval(abs(r%y-r2%y))<tol,'strided array output')

    r = dlsim_sos_bank(rsos_bank,u,-0.1_dp)
    call invalid(r,'negative timestep')
    r = dlsim_sos_bank(rsos_bank,u,ieee_value(0.0_dp,ieee_quiet_nan))
    call invalid(r,'nonfinite timestep')
    r = dlsim_sos_bank(rsos_bank,u(:0,:),0.1_dp)
    call invalid(r,'empty time axis')
    r = dlsim_sos_bank(rsos_bank,u(:,:0),0.1_dp)
    call invalid(r,'empty channel axis')
    r = dlsim_sos_bank(rsos_bank(:0,:,:),u,0.1_dp)
    call invalid(r,'empty section axis')
    r = dlsim_sos_bank(rsos_bank(:,:5,:),u,0.1_dp)
    call invalid(r,'incorrect SOS columns')
    r = dlsim_sos_bank(rsos_bank(:,:,:2),u,0.1_dp)
    call invalid(r,'unequal bank channel count')
    r = dlsim_sos_bank(rsos_bank,u,0.1_dp,zi=zi(:1,:,:))
    call invalid(r,'state section count')
    r = dlsim_sos_bank(sos,u,0.1_dp,zi=zi(:,:,:2))
    call invalid(r,'shared bank state channel count')
    r = dlsim_sos_bank(rsos_bank,u,0.1_dp,t=t(:ns-1))
    call invalid(r,'timestamp shape')
    r = dlsim_sos_bank(rsos_bank,u,0.1_dp,t=t(ns:1:-1))
    call invalid(r,'timestamp monotonicity')
    badsos = rsos_bank
    badsos(1,4,2) = 2.0_dp
    r = dlsim_sos_bank(badsos,u,0.1_dp)
    call invalid(r,'non-unit SOS a0')
    badsos = rsos_bank
    badsos(1,3,3) = ieee_value(0.0_dp,ieee_quiet_nan)
    r = dlsim_sos_bank(badsos,u,0.1_dp)
    call invalid(r,'nonfinite SOS')
    badu = u
    badu(7,2) = ieee_value(0.0_dp,ieee_quiet_nan)
    r = dlsim_sos_bank(rsos_bank,badu,0.1_dp)
    call invalid(r,'nonfinite channel input')
    badzi = zi
    badzi(1,1,3) = ieee_value(0.0_dp,ieee_quiet_nan)
    r = dlsim_sos_bank(rsos_bank,u,0.1_dp,zi=badzi)
    call invalid(r,'nonfinite initial channel state')
    r = dlsim_sos_bank(sos,u(:2,:),0.1_dp,t=[0.0_dp,0.04_dp])
    call check(r%status==discrete_sim_success,'single output timestamp')
    call check(size(r%t)==1,'single output time length')
    call check(maxval(abs(r%y(1,:)-cmplx(sos(1,1)*sos(2,1)*u(1,:),0.0_dp,dp))) < tol, &
        'single output input sample')
    badsos = rsos_bank
    badsos(1,1,2) = huge(1.0_dp)
    badsos(2,1,2) = huge(1.0_dp)
    r = dlsim_sos_bank(badsos,u,0.1_dp)
    call check(r%status==discrete_sim_numerical_failure,'overflow status')

    print '(a)', 'PASS multichannel SOS simulation'
contains
    subroutine check(ok, message)
        logical, intent(in) :: ok !! assertion outcome
        character(*), intent(in) :: message !! assertion description
        if (.not. ok) then
            print '(2a)', 'FAIL: ',message
            error stop 1
        end if
    end subroutine check

    subroutine match_channel(bank, scalar, channel, message)
        type(discrete_sos_bank_result), intent(in) :: bank !! bank result
        type(discrete_simulation_result), intent(in) :: scalar !! SISO reference
        integer, intent(in) :: channel !! one-based channel index
        character(*), intent(in) :: message !! assertion description
        call check(bank%status==discrete_sim_success,message//' status')
        call check(scalar%status==discrete_sim_success,message//' scalar status')
        call check(maxval(abs(bank%t-scalar%t))<tol,message//' times')
        call check(maxval(abs(bank%y(:,channel)-scalar%y(:,1)))<tol,message//' output')
        call check(maxval(abs(bank%x(:,:,channel)-scalar%x))<tol,message//' states')
        call check(maxval(abs(bank%zf(:,:,channel)-scalar%zf))<tol,message//' final states')
    end subroutine match_channel

    subroutine invalid(bank, message)
        type(discrete_sos_bank_result), intent(in) :: bank !! invalid result
        character(*), intent(in) :: message !! assertion description
        call check(bank%status==discrete_sim_invalid_input,message)
        call check(.not. allocated(bank%y),'invalid output not allocated')
        call check(.not. allocated(bank%zf),'invalid states not allocated')
    end subroutine invalid
end program test_signal_dlsim_sos_bank
