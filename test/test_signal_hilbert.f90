! SPDX-License-Identifier: MIT
! Tests for SciPy-compatible Hilbert analytic signals in one and two dimensions.
program test_signal_hilbert
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_signal_hilbert, only : hilbert, hilbert2, &
        hilbert_status_success, hilbert_status_invalid_input
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: tol = 1.0e-12_dp
    real(dp) :: x(4), a(3,4), angle, nan_value
    real(dp), allocatable :: tone(:), plane(:,:)
    complex(dp), allocatable :: y(:), z(:,:)
    integer :: st, k, j, nr, nc

    x = [0.2_dp, 1.0_dp, -0.5_dp, 2.0_dp]
    y = hilbert(x,status=st)
    if (st /= hilbert_status_success) error stop 'hilbert default status'
    call check_real(real(y,dp), x)
    call check_real(aimag(y),[0.5_dp,0.35_dp,-0.5_dp,-0.35_dp])

    y = hilbert(x,n=5,status=st)
    if (st /= hilbert_status_success) error stop 'hilbert padded status'
    if (size(y) /= 5) error stop 'hilbert padded shape'
    call check_real(real(y,dp),[x,0.0_dp])
    if (abs(aimag(y(1))+0.9788079714377309_dp) > tol) error stop 'hilbert padded reference'

    y = hilbert(x,n=3,status=st)
    if (st /= hilbert_status_success) error stop 'hilbert truncated status'
    call check_real(real(y,dp), x(:3))
    if (abs(aimag(y(1))+0.866025403784439_dp) > tol) error stop 'hilbert truncated reference'

    do nr = 1, 34
        allocate(tone(nr))
        do k = 1, nr
            tone(k) = 2.1_dp
            if (nr >= 3) tone(k) = tone(k) + 0.8_dp * cos(2.0_dp*pi*real(k-1,dp)/real(nr,dp))
        end do
        y = hilbert(tone,status=st)
        if (st /= hilbert_status_success) error stop 'hilbert periodic signal status'
        do k = 1,nr
            angle = 0.0_dp
            if (nr >= 3) angle = 0.8_dp * sin(2.0_dp*pi*real(k-1,dp)/real(nr,dp))
            if (abs(aimag(y(k))-angle) > tol) error stop 'hilbert periodic imaginary part'
        end do
        call check_real(real(y,dp),tone)
        deallocate(tone)
    end do

    a(1,:) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    a(2,:) = [0.0_dp,-1.0_dp, 2.0_dp, 0.0_dp]
    a(3,:) = [1.0_dp, 1.0_dp, 0.0_dp,-2.0_dp]
    z = hilbert2(a,status=st)
    if (st /= hilbert_status_success) error stop 'hilbert2 default status'
    if (any(shape(z) /= [3,4])) error stop 'hilbert2 default shape'
    call check_complex_scalar(z(1,1),cmplx(2.6547005383792515_dp,1.7216878364870323_dp,dp))
    call check_complex_scalar(z(2,3),cmplx(2.693375672974064_dp,1.809401076758503_dp,dp))
    if (abs(sum(z)-cmplx(11.0_dp,0.0_dp,dp)) > tol) error stop 'hilbert2 DC sum'

    z = hilbert2(a,n=[5,6],status=st)
    if (st /= hilbert_status_success) error stop 'hilbert2 padded status'
    if (any(shape(z) /= [5,6])) error stop 'hilbert2 padded shape'
    call check_complex_scalar(z(1,1),cmplx(1.7726075218913768_dp,-0.7500598475943444_dp,dp))
    call check_complex_scalar(z(2,3),cmplx(3.2769014186686123_dp,1.679617658072226_dp,dp))

    z = hilbert2(a,n=[2,3],status=st)
    if (st /= hilbert_status_success) error stop 'hilbert2 truncated status'
    if (any(shape(z) /= [2,3])) error stop 'hilbert2 truncated shape'
    call check_complex_scalar(z(1,1),cmplx(0.5_dp,1.1547005383792515_dp,dp))
    call check_complex_scalar(z(2,3),cmplx(2.5_dp,0.0_dp,dp))

    z = hilbert2(a,n=4,status=st)
    if (st /= hilbert_status_success) error stop 'hilbert2 square status'
    if (any(shape(z) /= [4,4])) error stop 'hilbert2 square shape'
    call check_complex_scalar(z(1,1),cmplx(1.0625_dp,1.625_dp,dp))

    ! Single-orthant analytic signal of cos(row frequency) cos(column frequency).
    do nr = 3,6
        do nc = 3,7
            allocate(plane(nr,nc))
            do j = 1,nc
                do k = 1,nr
                    plane(k,j) = cos(2.0_dp*pi*real(k-1,dp)/real(nr,dp)) * &
                        cos(2.0_dp*pi*real(j-1,dp)/real(nc,dp))
                end do
            end do
            z = hilbert2(plane,status=st)
            if (st /= hilbert_status_success) error stop 'hilbert2 separable status'
            do j = 1,nc
                do k = 1,nr
                    angle = 2.0_dp*pi*(real(k-1,dp)/real(nr,dp)+real(j-1,dp)/real(nc,dp))
                    call check_complex_scalar(z(k,j),cmplx(cos(angle),sin(angle),dp))
                end do
            end do
            deallocate(plane)
        end do
    end do

    y = hilbert(x,n=0,status=st)
    if (st /= hilbert_status_invalid_input .or. size(y) /= 0) error stop 'hilbert invalid length'
    z = hilbert2(a,n=[3,0],status=st)
    if (st /= hilbert_status_invalid_input .or. size(z) /= 0) error stop 'hilbert2 invalid lengths'
    z = hilbert2(a,n=[4],status=st)
    if (st /= hilbert_status_invalid_input .or. size(z) /= 0) error stop 'hilbert2 invalid shape'
    y = hilbert([real(dp)::],n=3,status=st)
    if (st /= hilbert_status_success) error stop 'hilbert empty padded status'
    call check_real(real(y,dp),[0.0_dp,0.0_dp,0.0_dp])
    nan_value = ieee_value(0.0_dp,ieee_quiet_nan)
    y = hilbert([1.0_dp,nan_value],status=st)
    if (st /= hilbert_status_invalid_input .or. size(y) /= 0) error stop 'hilbert nonfinite'
    a(1,1) = nan_value
    z = hilbert2(a,status=st)
    if (st /= hilbert_status_invalid_input .or. size(z) /= 0) error stop 'hilbert2 nonfinite'

    print '(a)', 'test_signal_hilbert passed'
contains
    subroutine check_real(got,expected)
        real(dp), intent(in) :: got(:) !! returned samples
        real(dp), intent(in) :: expected(:) !! target samples
        if (size(got) /= size(expected)) error stop 'hilbert real shape mismatch'
        if (size(got) > 0) then
            if (maxval(abs(got-expected)) > tol) error stop 'hilbert real mismatch'
        end if
    end subroutine check_real

    subroutine check_complex_scalar(got,expected)
        complex(dp), intent(in) :: got !! calculated value
        complex(dp), intent(in) :: expected !! expected value
        if (abs(got-expected) > tol) then
            print *, 'hilbert2 mismatch', got, expected
            error stop 'hilbert2 complex mismatch'
        end if
    end subroutine check_complex_scalar
end program test_signal_hilbert
