! SPDX-License-Identifier: MIT
program test_signal_cont2discrete
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_signal, only : dp, cont2discrete, cont2discrete_zpk, &
        discrete_state_space_result, discrete_tf_result, discrete_zpk_result, &
        cont2discrete_success, cont2discrete_invalid_input, cont2discrete_numerical_failure
    implicit none
    real(dp), parameter :: tol = 3.0e-12_dp
    real(dp), parameter :: dt = 0.1_dp
    real(dp) :: a(2,2), b(2,1), c(1,2), d(1,1), nanvalue
    real(dp) :: aa(0,0), bb(0,1), cc(1,0), dd(1,1)
    complex(dp) :: az(1,1), bz(1,1), cz(1,1), dz(1,1)
    type(discrete_state_space_result) :: sd
    type(discrete_tf_result) :: tf
    type(discrete_zpk_result) :: zp

    a = reshape([0.0_dp,-10.0_dp,1.0_dp,-3.0_dp],[2,2])
    b(:,1) = [0.0_dp,10.0_dp]
    c(1,:) = [1.0_dp,0.0_dp]
    d = 0.0_dp

    sd = cont2discrete(a,b,c,d,dt)
    call check(sd%status==cont2discrete_success,'default ZOH status')
    call check(abs(sd%dt-dt)<tol,'sample interval')
    call check(abs(sd%a(1,1)-0.955015412674261_dp)<tol,'ZOH A11')
    call check(abs(sd%a(1,2)-0.084963349921582_dp)<tol,'ZOH A12')
    call check(abs(sd%a(2,1)+0.849633499215824_dp)<tol,'ZOH A21')
    call check(abs(sd%a(2,2)-0.700125362909513_dp)<tol,'ZOH A22')
    call check(abs(sd%b(1,1)-0.044984587325739_dp)<tol,'ZOH B1')
    call check(abs(sd%b(2,1)-0.849633499215824_dp)<tol,'ZOH B2')
    call check(maxval(abs(sd%c-cmplx(c,0.0_dp,dp)))<tol,'ZOH C')
    call check(maxval(abs(sd%d))<tol,'ZOH D')

    sd = cont2discrete(a,b,c,d,dt,method='foh')
    call check(sd%status==cont2discrete_success,'FOH status')
    call check(abs(sd%b(1,1)-0.082511663971297_dp)<tol,'FOH B1')
    call check(abs(sd%b(2,1)-0.701640952021056_dp)<tol,'FOH B2')
    call check(abs(sd%d(1,1)-0.015412738806957_dp)<tol,'FOH D')

    sd = cont2discrete(a,b,c,d,dt,method='gbt',alpha=0.35_dp)
    call check(sd%status==cont2discrete_success,'GBT status')
    call check(abs(sd%a(1,1)-0.968673081226225_dp)<tol,'GBT A')
    call check(abs(sd%b(2,1)-0.895054822107854_dp)<tol,'GBT B')
    call check(abs(sd%c(1,1)-0.989035578429179_dp)<tol,'GBT C')
    call check(abs(sd%d(1,1)-0.010964421570821_dp)<tol,'GBT D')

    sd = cont2discrete(a,b,c,d,dt,method='bilinear')
    call check(sd%status==cont2discrete_success,'bilinear status')
    call check(abs(sd%a(1,1)-0.957446808510638_dp)<tol,'bilinear A')
    call check(abs(sd%d(1,1)-0.021276595744681_dp)<tol,'bilinear D')
    sd = cont2discrete(a,b,c,d,dt,method='TUSTIN')
    call check(abs(sd%a(1,1)-0.957446808510638_dp)<tol,'Tustin alias')

    sd = cont2discrete(a,b,c,d,dt,method='euler')
    call check(sd%status==cont2discrete_success,'Euler status')
    call check(abs(sd%a(1,2)-0.1_dp)<tol,'Euler A')
    call check(abs(sd%b(2,1)-1.0_dp)<tol,'Euler B')
    sd = cont2discrete(a,b,c,d,dt,method='forward_diff')
    call check(abs(sd%b(2,1)-1.0_dp)<tol,'forward_diff alias')

    sd = cont2discrete(a,b,c,d,dt,method='backward_diff')
    call check(sd%status==cont2discrete_success,'backward_diff status')
    call check(abs(sd%a(1,1)-0.928571428571429_dp)<tol,'backward A')
    call check(abs(sd%d(1,1)-0.071428571428571_dp)<tol,'backward D')

    sd = cont2discrete(a,b,c,d,dt,method='impulse')
    call check(sd%status==cont2discrete_success,'impulse status')
    call check(abs(sd%b(1,1)-0.084963349921582_dp)<tol,'impulse B')
    call check(maxval(abs(sd%d))<tol,'impulse D')

    tf = cont2discrete([1.0_dp],[1.0_dp,2.0_dp,1.0_dp],0.2_dp,method='bilinear')
    call check(tf%status==cont2discrete_success,'TF bilinear status')
    call check(size(tf%b,1)==1 .and. size(tf%a)==3,'TF shapes')
    call check(maxval(abs(tf%b(1,:)-cmplx([1.0_dp,2.0_dp,1.0_dp]/121.0_dp,0.0_dp,dp)))<tol, &
        'TF bilinear numerator')
    call check(maxval(abs(tf%a-cmplx([1.0_dp,-18.0_dp/11.0_dp,81.0_dp/121.0_dp],0.0_dp,dp)))<tol, &
        'TF bilinear denominator')

    tf = cont2discrete(reshape([1.0_dp,2.0_dp,0.5_dp,1.5_dp],[2,2]), &
        [1.0_dp,2.0_dp,1.0_dp],dt,method='zoh')
    call check(tf%status==cont2discrete_success,'SIMO transfer function status')
    call check(size(tf%b,1)==2 .and. size(tf%a)==3,'SIMO shapes')

    zp = cont2discrete_zpk([real(dp) :: ],[-1.0_dp,-2.0_dp,-3.0_dp,-4.0_dp], &
        1.0_dp,dt,method='bilinear')
    call check(zp%status==cont2discrete_success,'exact repeated bilinear zeros')
    call check(size(zp%zeros)==4,'four bilinear zeros')
    call check(maxval(abs(zp%zeros+1.0_dp))<tol,'bilinear zeros are exactly -1')
    zp = cont2discrete_zpk([real(dp) :: ],[-1.0_dp,-2.0_dp,-3.0_dp,-4.0_dp], &
        1.0_dp,dt,method='gbt',alpha=0.35_dp)
    call check(zp%status==cont2discrete_success,'exact repeated GBT zeros')
    call check(maxval(abs(zp%zeros+0.65_dp/0.35_dp))<tol,'GBT artificial zeros')

    zp = cont2discrete_zpk([-0.5_dp],[-1.0_dp,-2.0_dp], &
        cmplx(2.0_dp,0.3_dp,dp),dt,method='bilinear')
    call check(zp%status==cont2discrete_success,'complex gain with real roots')
    zp = cont2discrete_zpk([cmplx(-0.5_dp,0.2_dp,dp)], &
        [cmplx(-1.0_dp,0.1_dp,dp)],2.0_dp,dt,method='bilinear')
    call check(zp%status==cont2discrete_success,'complex roots with real gain')

    zp = cont2discrete_zpk([-0.5_dp],[-1.0_dp,-2.0_dp],2.0_dp,dt)
    call check(zp%status==cont2discrete_success,'ZPK status')
    call check(size(zp%poles)==2 .and. size(zp%zeros)==1,'ZPK shapes')
    call check(minval(abs(zp%poles-exp(-dt)))<tol,'ZPK pole 1')
    call check(minval(abs(zp%poles-exp(-2.0_dp*dt)))<tol,'ZPK pole 2')

    az(1,1) = cmplx(-0.5_dp,0.4_dp,dp)
    bz(1,1) = cmplx(2.0_dp,-0.2_dp,dp)
    cz(1,1) = cmplx(1.0_dp,0.6_dp,dp)
    dz = (0.0_dp,0.0_dp)
    sd = cont2discrete(az,bz,cz,dz,0.3_dp,method='zoh')
    call check(sd%status==cont2discrete_success,'complex state-space')
    call check(abs(sd%a(1,1)-exp(az(1,1)*0.3_dp))<tol,'complex exponential')

    dd = 4.0_dp
    sd = cont2discrete(aa,bb,cc,dd,dt,method='foh')
    call check(sd%status==cont2discrete_success,'empty state')
    call check(size(sd%a,1)==0 .and. abs(sd%d(1,1)-4.0_dp)<tol,'empty state matrices')

    sd = cont2discrete(a,b,c,d,0.0_dp)
    call check(sd%status==cont2discrete_invalid_input,'zero dt rejection')
    sd = cont2discrete(a,b,c,d,-0.1_dp)
    call check(sd%status==cont2discrete_invalid_input,'negative dt rejection')
    nanvalue = ieee_value(0.0_dp,ieee_quiet_nan)
    sd = cont2discrete(a,b,c,d,nanvalue)
    call check(sd%status==cont2discrete_invalid_input,'nan dt rejection')
    sd = cont2discrete(a,b,c,d,dt,method='gbt')
    call check(sd%status==cont2discrete_invalid_input,'missing alpha rejection')
    sd = cont2discrete(a,b,c,d,dt,method='gbt',alpha=1.2_dp)
    call check(sd%status==cont2discrete_invalid_input,'invalid alpha rejection')
    sd = cont2discrete(a,b,c,d,dt,method='unknown')
    call check(sd%status==cont2discrete_invalid_input,'unknown method rejection')
    dd = 4.0_dp
    sd = cont2discrete(a,b,c,dd,dt,method='impulse')
    call check(sd%status==cont2discrete_invalid_input,'impulse D rejection')
    sd = cont2discrete(a(1:1,:),b,c,d,dt)
    call check(sd%status==cont2discrete_invalid_input,'bad dimensions rejection')
    sd = cont2discrete(reshape([20.0_dp],[1,1]),reshape([1.0_dp],[1,1]), &
        reshape([1.0_dp],[1,1]),reshape([0.0_dp],[1,1]),dt,method='bilinear')
    call check(sd%status==cont2discrete_numerical_failure,'singular bilinear matrix')

    print '(a)', 'PASS: cont2discrete'
contains
    subroutine check(ok,label)
        logical, intent(in) :: ok !! assertion condition
        character(*), intent(in) :: label !! diagnostic name
        if (.not. ok) then
            print '(a)', 'FAIL: '//label
            error stop 1
        end if
    end subroutine check
end program test_signal_cont2discrete
