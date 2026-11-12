! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_group_delay_sos
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_design, only : group_delay, group_delay_sos, group_delay_result, &
        design_status_success, design_status_invalid_input, design_status_singular
    implicit none

    type(group_delay_result) :: gd, reference
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: w(5) = [0.1_dp, 0.5_dp, 1.0_dp, 2.1_dp, 3.0_dp]
    real(dp) :: sos(2,6), bad(2,6)
    complex(dp) :: csos(2,6)

    sos = 0.0_dp
    sos(:,4) = 1.0_dp
    sos(1,2) = 1.0_dp
    sos(2,3) = 1.0_dp
    gd = group_delay_sos(sos,w=w)
    call check(gd%status == design_status_success, 'pure delay status')
    call close(gd%frequencies,w, 'pure delay grid')
    call close(gd%delay, [3.0_dp,3.0_dp,3.0_dp,3.0_dp,3.0_dp], 'pure delay values')

    ! The overloaded name group_delay accepts SOS rank-2 arrays as well.
    gd = group_delay(sos,whole=.true.,worn=24)
    call check(gd%status == design_status_success,'generic SOS status')
    call check(size(gd%delay) == 24, 'generic SOS grid size')
    call check(all(abs(gd%delay-3.0_dp) < 1.0e-12_dp),'generic SOS delay')
    call check(abs(gd%frequencies(24)-2.0_dp*pi*23.0_dp/24.0_dp) < 1.0e-12_dp, &
        'generic SOS whole grid')

    sos = 0.0_dp
    sos(1,:) = [0.8_dp,-0.2_dp,0.1_dp,1.0_dp,-0.5_dp,0.2_dp]
    sos(2,:) = [1.1_dp,0.3_dp,-0.05_dp,1.0_dp,0.2_dp,0.04_dp]
    gd = group_delay_sos(sos,w=w)
    reference = group_delay( &
        [sos(1,1)*sos(2,1),sos(1,1)*sos(2,2)+sos(1,2)*sos(2,1), &
         sos(1,1)*sos(2,3)+sos(1,2)*sos(2,2)+sos(1,3)*sos(2,1), &
         sos(1,2)*sos(2,3)+sos(1,3)*sos(2,2),sos(1,3)*sos(2,3)], &
        [sos(1,4)*sos(2,4),sos(1,4)*sos(2,5)+sos(1,5)*sos(2,4), &
         sos(1,4)*sos(2,6)+sos(1,5)*sos(2,5)+sos(1,6)*sos(2,4), &
         sos(1,5)*sos(2,6)+sos(1,6)*sos(2,5),sos(1,6)*sos(2,6)],w=w)
    call check(gd%status == design_status_success, 'two-section status')
    call check(reference%status == design_status_success,'transfer-function reference')
    call close(gd%delay,reference%delay,'two-section TF equivalence')

    csos = cmplx(sos,0.0_dp,dp)
    csos(1,2) = cmplx(-0.2_dp,0.15_dp,dp)
    csos(2,6) = cmplx(0.04_dp,-0.1_dp,dp)
    gd = group_delay_sos(csos,w=w)
    call check(gd%status == design_status_success, 'complex SOS status')
    reference = group_delay(csos(1,1:3),csos(1,4:6),w=w)
    ! Full cascade delay equals the sum of its sections' delays.
    gd = group_delay(csos,w=w)
    call check(gd%status == design_status_success, 'complex generic status')
    call close(gd%delay(1:5),reference%delay(1:5)+ &
        group_delay_single(csos(2,1:3),csos(2,4:6)), 'complex section sum')

    gd = group_delay_sos(sos,fs=200.0_dp,w=[2.0_dp,12.0_dp,50.0_dp,90.0_dp])
    call check(gd%status == design_status_success, 'sampling-rate status')
    call check(abs(gd%frequencies(4)-90.0_dp)<1.0e-12_dp,'sampling-rate grid')

    bad = sos
    bad(1,1:3) = 0.0_dp
    gd = group_delay_sos(bad)
    call check(gd%status == design_status_invalid_input, 'zero numerator')
    call check(.not. allocated(gd%delay), 'zero numerator result')
    bad = sos
    bad(2,4) = 2.0_dp
    gd = group_delay_sos(bad)
    call check(gd%status == design_status_invalid_input, 'a0 must be one')
    gd = group_delay_sos(sos(:,1:5))
    call check(gd%status == design_status_invalid_input, 'invalid SOS width')
    gd = group_delay_sos(sos(1:0,:))
    call check(gd%status == design_status_invalid_input, 'empty SOS')
    gd = group_delay_sos(sos,worn=0)
    call check(gd%status == design_status_invalid_input,'invalid frequency count')
    gd = group_delay_sos(sos,fs=0.0_dp)
    call check(gd%status == design_status_invalid_input,'invalid sampling rate')
    gd = group_delay_sos(sos,w=[ieee_value(0.0_dp,ieee_quiet_nan)])
    call check(gd%status == design_status_invalid_input,'invalid frequencies')
    bad = sos
    bad(1,2) = ieee_value(0.0_dp,ieee_quiet_nan)
    gd = group_delay_sos(bad)
    call check(gd%status == design_status_invalid_input,'NaN SOS')

    bad = 0.0_dp
    bad(:,4) = 1.0_dp
    bad(:,1) = 1.0_dp
    bad(1,2) = 1.0_dp
    gd = group_delay_sos(bad,w=[pi,0.3_dp])
    call check(gd%status == design_status_singular,'unit-circle zero status')
    call check(ieee_is_nan(gd%delay(1)),'unit-circle zero NaN')
    call check(.not. ieee_is_nan(gd%delay(2)),'unit-circle other bin')
    bad(1,2) = 0.0_dp
    bad(1,5) = -1.0_dp
    gd = group_delay_sos(bad,w=[0.0_dp,0.3_dp])
    call check(gd%status == design_status_singular,'unit-circle pole status')
    call check(ieee_is_nan(gd%delay(1)),'unit-circle pole NaN')
    call check(.not. ieee_is_nan(gd%delay(2)),'unit-circle pole other bin')

    print '(a)', 'test_signal_group_delay_sos passed'
contains
    subroutine check(ok,label)
        logical, intent(in) :: ok !! validated condition
        character(*), intent(in) :: label !! name of the failed condition
        if (.not. ok) then
            print '(a)', 'failure: '//label
            error stop 'group_delay_sos check failed'
        end if
    end subroutine check

    subroutine close(actual,expected,label)
        real(dp), intent(in) :: actual(:) !! computed group delay samples
        real(dp), intent(in) :: expected(:) !! comparison samples
        character(*), intent(in) :: label !! failing comparison name
        call check(size(actual) == size(expected),label//' shape')
        call check(all(abs(actual-expected)<2.0e-11_dp),label//' values')
    end subroutine close

    function group_delay_single(b,a) result(delay)
        complex(dp), intent(in) :: b(:) !! numerator of one section
        complex(dp), intent(in) :: a(:) !! denominator of one section
        real(dp), allocatable :: delay(:) !! group delay along w
        type(group_delay_result) :: values
        values = group_delay(b,a,w=w)
        call check(values%status == design_status_success,'single complex section')
        delay = values%delay
    end function group_delay_single
end program test_signal_group_delay_sos
