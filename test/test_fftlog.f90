! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_fftlog
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_nan
    use scifort_fft, only : dp, fht, ifht, fhtoffset, fht_nd, ifht_nd, &
        fft_status_singular, fft_status_range_error
    implicit none
    real(dp), allocatable :: a(:), b(:), c(:), r(:), k(:), expected(:), saved(:), gap(:)
    real(dp), parameter :: orders(8) = [-4.0_dp,-3.0_dp,-1.2_dp,-0.5_dp,0.0_dp,0.5_dp,2.0_dp,50.0_dp]
    real(dp), parameter :: biases(3) = [-0.3_dp,0.0_dp,0.3_dp]
    real(dp) :: batch(3,7,2), out(3,7,2), back(3,7,2), strided(6,7,2)
    real(dp) :: high(2,1,1,1,1,1,1,1,1,1,1,1,1,1,3)
    real(dp) :: hout(2,1,1,1,1,1,1,1,1,1,1,1,1,1,3)
    real(dp) :: hback(2,1,1,1,1,1,1,1,1,1,1,1,1,1,3)
    real(dp) :: empty(0,3), empty_out(0,3)
    real(dp) :: dln, mu, q, shift, value, nan, inf, center, factor, tolerance
    integer :: n, i, j, m, code, checks, failures
    checks = 0
    failures = 0
    nan = ieee_value(0.0_dp,ieee_quiet_nan)
    inf = ieee_value(0.0_dp,ieee_positive_inf)
    dln = 0.2_dp
    do n = 1, 65
        a = [(sin(0.7_dp*real(i,dp))+0.3_dp*cos(0.31_dp*real(i,dp)),i=1,n)]
        do j = 1, size(orders)
            mu = orders(j)
            do m = 1, size(biases)
                q = biases(m)
                shift = fhtoffset(dln,mu,initial=0.3_dp,bias=q,status=code)
                call check(code == 0,'offset status')
                call check(abs(shift-0.3_dp) <= 0.5_dp*dln+epsilon(dln),'nearest low-ringing offset')
                b = fht(a,dln,mu,shift,q,status=code)
                call check(code == 0 .and. size(b) == n,'forward status and size')
                if (code /= 0) cycle
                c = ifht(b,dln,mu,shift,q,status=code)
                call check(code == 0,'inverse status')
                tolerance = 3e-12_dp*exp(2*abs(q)*dln*real(n-1,dp))
                call close(c,a,tolerance,'round trip')
                if (q == 0) then
                    call check(abs(sum(b*b)-sum(a*a)) < 3e-11_dp*max(1.0_dp,sum(a*a)), &
                        'low-ringing unbiased norm preservation')
                    c = fht(b,dln,mu,shift,q)
                    call close(c,a,3e-11_dp,'unbiased self inverse')
                end if
            end do
        end do
    end do
    do n = 1, 17
        q = 0.5_dp
        mu = 2.0_dp
        shift = 0.1_dp
        center = 0.5_dp*real(n-1,dp)
        a = [(exp(q*(real(i-1,dp)-center)*dln),i=1,n)]
        factor = 2.0_dp**q*gamma(0.5_dp*(mu+1+q))/gamma(0.5_dp*(mu+1-q))
        expected = [(factor*exp(-q*((real(i-1,dp)-center)*dln+shift)),i=1,n)]
        b = fht(a,dln,mu,shift,q)
        call close(b,expected,3e-12_dp,'analytic Mellin power law')
    end do
    n = 512
    dln = 10*log(10.0_dp)/real(n-1,dp)
    r = [(exp(-8*log(10.0_dp)+real(i-1,dp)*dln),i=1,n)]
    shift = fhtoffset(dln,0.0_dp,initial=-6*log(10.0_dp))
    k = exp(shift)/r(n:1:-1)
    a = r*exp(-0.5_dp*r*r)
    expected = k*exp(-0.5_dp*k*k)
    b = fht(a,dln,0.0_dp,shift)
    call check(maxval(abs(b-expected),mask=k >= 0.01_dp .and. k <= 5.0_dp) < 2e-7_dp, &
        'continuous Gaussian Hankel identity on resolved interior')
    dln = -0.2_dp
    a = [(sin(real(i,dp)),i=1,9)]
    shift = fhtoffset(dln,-0.5_dp)
    b = fht(a,dln,-0.5_dp,shift,bias=0.2_dp)
    c = ifht(b,dln,-0.5_dp,shift,bias=0.2_dp)
    call close(c,a,3e-12_dp,'descending logarithmic grid')
    gap = [(real(i,dp),i=1,18)]
    saved = gap
    b = fht(gap(1:18:2),0.2_dp,0.0_dp)
    c = fht(saved(1:18:2),0.2_dp,0.0_dp)
    call close(b,c,0.0_dp,'strided rank-one input')
    call check(all(gap == saved),'input preservation')
    batch = reshape([(sin(real(i,dp)),i=1,42)],shape(batch))
    shift = fhtoffset(0.2_dp,0.0_dp)
    call fht_nd(batch,out,0.2_dp,0.0_dp,shift,axis=2,status=code)
    call check(code == 0,'batched status')
    do j = 1, 2
        do i = 1, 3
            b = fht(batch(i,:,j),0.2_dp,0.0_dp,shift)
            call close(out(i,:,j),b,0.0_dp,'batch/scalar agreement')
        end do
    end do
    call ifht_nd(out,back,0.2_dp,0.0_dp,shift,axis=2,status=code)
    call check(code == 0 .and. maxval(abs(back-batch)) < 3e-12_dp,'batched round trip')
    call fht_nd(batch,out,0.2_dp,0.0_dp,status=code)
    call check(code == 0,'default last axis')
    strided(1:6:2,:,:) = batch
    call fht_nd(strided(1:6:2,:,:),strided(2:6:2,:,:),0.2_dp,0.0_dp,status=code)
    call check(code == 0 .and. maxval(abs(strided(2:6:2,:,:)-out)) == 0, &
        'strided batched input and output')
    high = reshape([(real(i,dp),i=1,6)],shape(high))
    call fht_nd(high,hout,0.2_dp,0.0_dp,status=code)
    call ifht_nd(hout,hback,0.2_dp,0.0_dp,status=code)
    call check(code == 0 .and. maxval(abs(hback-high)) < 3e-12_dp,'rank fifteen')
    call fht_nd(empty,empty_out,0.2_dp,0.0_dp,status=code)
    call check(code == 0,'empty batch')
    call fht_nd(batch,out,0.2_dp,0.0_dp,axis=0,status=code)
    call check(code == 1 .and. all(ieee_is_nan(out)),'invalid batch axis')
    call fht_nd(batch,empty_out,0.2_dp,0.0_dp,status=code)
    call check(code == 1,'wrong rank output')
    a = [3.0_dp]
    b = fht(a,0.2_dp,-4.0_dp,bias=1.0_dp,status=code)
    call check(code == 0 .and. abs(b(1)+12.0_dp) < 1e-13_dp,'simultaneous gamma poles')
    c = ifht(b,0.2_dp,-4.0_dp,bias=1.0_dp)
    call close(c,a,1e-14_dp,'removable pole inverse')
    b = fht(a,0.2_dp,0.0_dp,bias=-1.0_dp,status=code)
    call check(code == fft_status_singular .and. size(b) == 0,'singular numerator pole')
    b = ifht(a,0.2_dp,0.0_dp,bias=1.0_dp,status=code)
    call check(code == fft_status_singular .and. size(b) == 0,'singular denominator pole')
    b = fht(a,0.2_dp,0.0_dp,bias=1.0_dp,status=code)
    call check(code == 0 .and. all(b == 0),'well-defined forward zero DC coefficient')
    b = ifht(a,0.2_dp,0.0_dp,bias=-1.0_dp,status=code)
    call check(code == 0 .and. all(b == 0),'well-defined inverse zero DC coefficient')
    shift = fhtoffset(0.2_dp,0.0_dp)+0.1_dp
    b = fht([1.0_dp,2.0_dp],0.2_dp,0.0_dp,shift,status=code)
    call check(code == fft_status_singular,'singular even Nyquist projection')
    b = fht([real(dp) ::],0.2_dp,0.0_dp,status=code)
    call check(code == 1 .and. size(b) == 0,'empty transform')
    b = fht(a,0.0_dp,0.0_dp,status=code)
    call check(code == 1,'zero spacing')
    b = fht([nan],0.2_dp,0.0_dp,status=code)
    call check(code == 1,'NaN input')
    b = ifht([inf],0.2_dp,0.0_dp,status=code)
    call check(code == 1,'infinite input')
    value = fhtoffset(0.2_dp,nan,status=code)
    call check(code == 1 .and. ieee_is_nan(value),'invalid offset parameters')
    b = fht(a,1e-300_dp,0.0_dp,status=code)
    call check(code == fft_status_range_error,'unrepresentable frequency range guarded')
    b = fht(a,0.2_dp,real(huge(0),dp),status=code)
    call check(code == fft_status_range_error,'log-gamma recurrence integer range guarded')
    b = fht(a,0.2_dp,0.0_dp,offset=1e-308_dp,status=code)
    call check(code == 0,'subnormal offset with zero bias')
    b = fht(a,0.2_dp,0.0_dp,offset=1000.0_dp,bias=2.0_dp,status=code)
    call check(code == fft_status_range_error,'excessive logarithmic bias span guarded')
    if (failures > 0) error stop 'FFTLog tests failed'
    print '(a,i0,a)', 'FFTLog: ',checks,' checks passed'
contains
    subroutine check(ok,label)
        implicit none
        logical, intent(in) :: ok !! success condition
        character(*), intent(in) :: label !! failure description
        checks = checks+1
        if (ok) return
        failures = failures+1
        print '(a)', 'FAIL: '//label
    end subroutine check
    subroutine close(actual,wanted,tolerance,label)
        implicit none
        real(dp), intent(in) :: actual(:) !! computed values
        real(dp), intent(in) :: wanted(:) !! reference values
        real(dp), intent(in) :: tolerance !! relative maximum-norm bound
        character(*), intent(in) :: label !! failure description
        call check(size(actual) == size(wanted),label//' size')
        if (size(actual) /= size(wanted)) return
        call check(maxval(abs(actual-wanted)) <= tolerance*max(1.0_dp,maxval(abs(wanted))),label)
    end subroutine close
end program test_fftlog
