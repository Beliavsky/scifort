! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_signal_resample
    use scifort_signal, only : dp, upfirdn, resample_poly, decimate, &
        resample_status_success, resample_status_invalid_input
    implicit none
    real(dp), allocatable :: y(:), h(:), d(:), part(:)
    complex(dp), allocatable :: z(:)
    integer :: st, k
    real(dp), parameter :: tol=2.0e-12_dp

    y=upfirdn([1.0_dp,2.0_dp,1.0_dp], [1.0_dp,3.0_dp,2.0_dp,5.0_dp], &
        up=2, down=3, status=st)
    if (st/=resample_status_success) error stop 'upfirdn status'
    call check(y,[1.0_dp,6.0_dp,7.0_dp],tol)
    y=upfirdn([0.5_dp,0.3_dp,0.2_dp],[2.0_dp,1.0_dp,4.0_dp], &
        up=3,down=2,mode='edge',status=st)
    call check(y,[1.0_dp,0.4_dp,0.3_dp,2.0_dp,0.8_dp],tol)
    if (st/=resample_status_success) error stop 'upfirdn edge status'

    z=upfirdn([cmplx(1.0_dp,2.0_dp,dp),cmplx(0.25_dp,-0.5_dp,dp)], &
        [cmplx(1.0_dp,2.0_dp,dp),cmplx(3.0_dp,-1.0_dp,dp)],up=2,status=st)
    if (st/=resample_status_success) error stop 'complex upfirdn status'
    if (maxval(abs(z-[cmplx(-3.0_dp,4.0_dp,dp),cmplx(1.25_dp,0.0_dp,dp), &
        cmplx(5.0_dp,5.0_dp,dp),cmplx(0.25_dp,-1.75_dp,dp)]))>tol) &
        error stop 'complex upfirdn values'
    z=upfirdn([cmplx(1.0_dp,1.0_dp,dp)], [2.0_dp,3.0_dp],status=st)
    if (st/=resample_status_success) error stop 'complex kernel, real signal status'
    if (maxval(abs(z-[cmplx(2.0_dp,2.0_dp,dp),cmplx(3.0_dp,3.0_dp,dp)]))>tol) &
        error stop 'complex kernel, real signal values'

    y=resample_poly([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp],2,3, &
        taps=[0.25_dp,0.5_dp,0.25_dp],status=st)
    if (st/=resample_status_success) error stop 'resample_poly status'
    call check(y,[1.0_dp,2.5_dp,4.0_dp,2.5_dp],tol)
    y=resample_poly([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp],2,3, &
        taps=[0.25_dp,0.5_dp,0.25_dp],padtype='mean',status=st)
    call check(y,[1.0_dp,2.5_dp,4.0_dp,4.0_dp],tol)
    y=resample_poly([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp],2,3, &
        taps=[0.25_dp,0.5_dp,0.25_dp],padtype='constant',cval=4.0_dp,status=st)
    call check(y,[1.0_dp,2.5_dp,4.0_dp,4.5_dp],tol)

    ! Reduced ratios yield exactly the same resampling result.
    h=[0.125_dp,0.375_dp,0.375_dp,0.125_dp]
    y=resample_poly([2.0_dp,3.0_dp,-1.0_dp,4.0_dp],6,4,taps=h,status=st)
    part=resample_poly([2.0_dp,3.0_dp,-1.0_dp,4.0_dp],3,2,taps=h)
    call check(y,part,tol)
    y=resample_poly([1.0_dp,2.0_dp,3.0_dp],4,4,status=st)
    call check(y,[1.0_dp,2.0_dp,3.0_dp],tol)

    ! FIR-only decimation: both zero-phase and causal filter paths.
    y=decimate([(real(k,dp),k=1,16)],3,n=12,ftype='fir',zero_phase=.true.,status=st)
    if (st/=resample_status_success) error stop 'decimate zero-phase status'
    call check(y,[1.0107558210645804_dp,3.9922273885324255_dp, &
        7.0_dp,10.0_dp,13.501431312558658_dp,10.30959926124143_dp],tol)
    y=decimate([(real(k,dp),k=1,16)],3,n=12,ftype='fir',zero_phase=.false.,status=st)
    call check(y,[0.0_dp,-0.06585010571358445_dp,1.0107558210645804_dp, &
        3.9922273885324255_dp,7.0_dp,10.0_dp],tol)

    ! Invalid requests return zero-length outputs and a defined status.
    y=upfirdn([1.0_dp],[1.0_dp],up=0,status=st)
    if (st/=resample_status_invalid_input .or. size(y)/=0) error stop 'bad up factor'
    y=upfirdn([1.0_dp],[1.0_dp],mode='unknown',status=st)
    if (st/=resample_status_invalid_input .or. size(y)/=0) error stop 'bad extension mode'
    y=resample_poly([1.0_dp],1,0,status=st)
    if (st/=resample_status_invalid_input .or. size(y)/=0) error stop 'bad resampling factor'
    y=resample_poly([1.0_dp],2,3,window='unknown',status=st)
    if (st/=resample_status_invalid_input .or. size(y)/=0) error stop 'bad design window'
    y=resample_poly([1.0_dp],2,3,padtype='edge',cval=1.0_dp,status=st)
    if (st/=resample_status_invalid_input .or. size(y)/=0) error stop 'bad cval for padtype'
    y=decimate([1.0_dp,2.0_dp],2,ftype='unknown',status=st)
    if (st/=resample_status_invalid_input .or. size(y)/=0) error stop 'unsupported filter type'

    print '(a)', 'test_signal_resample passed'
contains
    subroutine check(actual,expected,eps)
        real(dp), intent(in) :: actual(:) !! observed samples
        real(dp), intent(in) :: expected(:) !! SciPy reference samples
        real(dp), intent(in) :: eps !! absolute tolerance
        if (size(actual)/=size(expected)) error stop 'wrong resampling length'
        if (size(actual)>0) then
            if (maxval(abs(actual-expected))>eps) error stop 'resampling reference mismatch'
        end if
    end subroutine check
end program test_signal_resample
