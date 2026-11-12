! SPDX-License-Identifier: MIT
! Explicit reusable Bluestein FFT plans and public API regression.
program test_signal_fft_plans
    use scifort_kinds, only : dp
    use scifort_signal, only : fft_plan, forward_fft, forward_fft_planned, spectral_status_success, &
        spectral_status_invalid_input
    implicit none
    integer, parameter :: lengths(*) = [1,2,3,5,11,16,17,19,31,33,63,64,65,97, &
        101,257,509,1009,2053,4093]
    type(fft_plan) :: plan
    complex(dp), allocatable :: x(:), y(:), z(:), scratch(:), initial(:)
    real(dp) :: tolerance
    integer :: i, j, k, n, status, cases

    cases = 0
    do i = 1, size(lengths)
        n = lengths(i)
        x = [(cmplx(sin(real(j,dp)*0.21_dp)+0.1_dp*cos(real(j*j,dp)*0.004_dp), &
            cos(real(j,dp)*0.11_dp)-0.2_dp*sin(real(j,dp)*0.14_dp),dp),j=1,n)]
        tolerance = 2.0e-12_dp*max(1.0_dp,real(n,dp))
        call plan%initialize(n,status)
        if (status /= spectral_status_success) error stop 'plan initialize'
        if (plan%length() /= n) error stop 'plan length'
        initial = x
        ! Multiple successive transformations on the same plan.
        do k = 1, 4
            y = x*cmplx(1.0_dp+real(k,dp)*0.01_dp,0.04_dp*real(k,dp),dp)
            z = y
            call forward_fft(z)
            call forward_fft_planned(y,plan,status=status)
            if (status /= spectral_status_success) error stop 'planned forward status'
            if (maxval(abs(y-z)) > tolerance) error stop 'planned/unplanned mismatch'
            call plan%inverse(y,status)
            if (status /= spectral_status_success) error stop 'plan inverse status'
            if (maxval(abs(y-x*cmplx(1.0_dp+real(k,dp)*0.01_dp, &
                0.04_dp*real(k,dp),dp))) > tolerance) error stop 'planned inverse mismatch'
            cases = cases+1
        end do
        ! A rejected operation must not change the input or the plan.
        if (n /= 3) then
            scratch = [(cmplx(real(j,dp),-real(j,dp),dp),j=1,3)]
            initial = scratch
            call forward_fft_planned(scratch,plan,status)
            if (status /= spectral_status_invalid_input) error stop 'rejected mismatch status'
            if (any(scratch /= initial)) error stop 'rejected mismatch altered input'
        end if
        ! An invalid attempt to reinitialize must keep the existing good plan.
        call plan%initialize(0,status)
        if (status /= spectral_status_invalid_input) error stop 'invalid size accepted'
        if (plan%length() /= n) error stop 'invalid init modified plan'
        ! Strided arrays and unchanged gaps.
        if (allocated(scratch)) deallocate(scratch)
        allocate(scratch(2*n))
        scratch = (-123.0_dp, 7.0_dp)
        scratch(1:2*n:2) = x
        z = x
        call forward_fft(z)
        call plan%execute(scratch(1:2*n:2),status)
        if (status /= spectral_status_success) error stop 'strided planned FFT'
        if (maxval(abs(scratch(1:2*n:2)-z)) > tolerance) error stop 'strided mismatch'
        if (any(scratch(2:2*n:2) /= cmplx(-123.0_dp,7.0_dp,dp))) error stop 'strided gap'
        deallocate(scratch)
    end do
    call plan%clear()
    if (plan%length() /= 0) error stop 'plan clear length'
    x = [(cmplx(real(i,dp),-real(i,dp),dp),i=1,2)]
    y = x
    call plan%execute(x,status)
    if (status /= spectral_status_invalid_input .or. any(x /= y)) error stop 'cleared plan'
    call plan%inverse(x,status)
    if (status /= spectral_status_invalid_input .or. any(x /= y)) error stop 'cleared inverse'
    print '(a,i0,a)', 'test_signal_fft_plans passed (',cases,' repeated transforms)'
end program test_signal_fft_plans
