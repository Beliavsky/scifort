! SPDX-License-Identifier: MIT
! Verify mixed real/complex zero-pole-gain overloads against canonical casts.
program test_signal_continuous_mixed_zpk
    use scifort_signal, only : dp, step_zpk, impulse_zpk, lsim_zpk, &
        continuous_simulation_result, continuous_response_result, continuous_sim_success
    implicit none
    real(dp), parameter :: tol=2.0e-12_dp
    real(dp) :: zr(1), pr(2), gain_real, ur(11), t(11)
    complex(dp) :: zc(1), pc(2), gain_complex, uc(11), x0(2)
    type(continuous_response_result) :: got, ref
    type(continuous_simulation_result) :: sim, expected
    integer :: i

    zr = [-0.2_dp]
    pr = [-1.0_dp,-1.8_dp]
    zc(1) = cmplx(-0.2_dp,0.4_dp,dp)
    pc = [cmplx(-1.0_dp,0.5_dp,dp),cmplx(-1.8_dp,-0.1_dp,dp)]
    gain_real = 1.5_dp
    gain_complex = cmplx(1.0_dp,0.3_dp,dp)
    x0 = [cmplx(0.2_dp,-0.3_dp,dp),cmplx(-0.2_dp,0.1_dp,dp)]
    do i=1,11
        t(i) = real(i-1,dp)/10.0_dp
    end do
    ur = sin(t)
    uc = cmplx(ur,0.4_dp*ur,dp)

    got = step_zpk(zr,pr,gain_complex,n=47)
    ref = step_zpk(cmplx(zr,0.0_dp,dp),cmplx(pr,0.0_dp,dp),gain_complex,n=47)
    call check_response(got,ref,'real roots with complex gain: step')
    got = impulse_zpk(zr,pr,gain_complex,t,x0)
    ref = impulse_zpk(cmplx(zr,0.0_dp,dp),cmplx(pr,0.0_dp,dp),gain_complex,t,x0)
    call check_response(got,ref,'real roots with complex gain: impulse')

    got = step_zpk(zc,pc,gain_real,n=47)
    ref = step_zpk(zc,pc,cmplx(gain_real,0.0_dp,dp),n=47)
    call check_response(got,ref,'complex roots with real gain: step')
    got = impulse_zpk(zc,pc,gain_real,t,x0)
    ref = impulse_zpk(zc,pc,cmplx(gain_real,0.0_dp,dp),t,x0)
    call check_response(got,ref,'complex roots with real gain: impulse')

    sim = lsim_zpk(zr,pr,gain_real,uc,t,x0)
    expected = lsim_zpk(cmplx(zr,0.0_dp,dp),cmplx(pr,0.0_dp,dp), &
        cmplx(gain_real,0.0_dp,dp),uc,t,x0)
    call check_sim(sim,expected,'real roots/gain, complex input')
    sim = lsim_zpk(zc,pc,gain_complex,ur,t,x0)
    expected = lsim_zpk(zc,pc,gain_complex,cmplx(ur,0.0_dp,dp),t,x0)
    call check_sim(sim,expected,'complex roots/gain, real input')

    sim = lsim_zpk(zr,pr,gain_complex,uc,t,x0)
    expected = lsim_zpk(cmplx(zr,0.0_dp,dp),cmplx(pr,0.0_dp,dp),gain_complex,uc,t,x0)
    call check_sim(sim,expected,'real roots, complex gain/input')
    sim = lsim_zpk(zr,pr,gain_complex,ur,t,x0)
    expected = lsim_zpk(cmplx(zr,0.0_dp,dp),cmplx(pr,0.0_dp,dp), &
        gain_complex,cmplx(ur,0.0_dp,dp),t,x0)
    call check_sim(sim,expected,'real roots, complex gain, real input')
    sim = lsim_zpk(zc,pc,gain_real,uc,t,x0)
    expected = lsim_zpk(zc,pc,cmplx(gain_real,0.0_dp,dp),uc,t,x0)
    call check_sim(sim,expected,'complex roots, real gain, complex input')
    sim = lsim_zpk(zc,pc,gain_real,ur,t,x0)
    expected = lsim_zpk(zc,pc,cmplx(gain_real,0.0_dp,dp),cmplx(ur,0.0_dp,dp),t,x0)
    call check_sim(sim,expected,'complex roots, real gain/input')
    print '(a)', 'PASS: mixed real/complex ZPK continuous responses'
contains
    subroutine check_response(x,y,description)
        type(continuous_response_result), intent(in) :: x !! actual result
        type(continuous_response_result), intent(in) :: y !! reference result
        character(*), intent(in) :: description !! test description
        if (x%status/=continuous_sim_success .or. y%status/=continuous_sim_success) then
            print '(a)', 'FAIL: '//description//' status'
            error stop 1
        end if
        if (size(x%t)/=size(y%t) .or. any(x%t/=y%t) .or. &
            maxval(abs(x%y-y%y))>tol) then
            print '(a)', 'FAIL: '//description//' values'
            error stop 1
        end if
    end subroutine check_response
    subroutine check_sim(x,y,description)
        type(continuous_simulation_result), intent(in) :: x !! actual simulation
        type(continuous_simulation_result), intent(in) :: y !! reference simulation
        character(*), intent(in) :: description !! test description
        if (x%status/=continuous_sim_success .or. y%status/=continuous_sim_success) then
            print '(a)', 'FAIL: '//description//' status'
            error stop 1
        end if
        if (maxval(abs(x%y-y%y))>tol .or. maxval(abs(x%x-y%x))>tol) then
            print '(a)', 'FAIL: '//description//' values'
            error stop 1
        end if
    end subroutine check_sim
end program test_signal_continuous_mixed_zpk
