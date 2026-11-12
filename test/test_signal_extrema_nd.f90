! SPDX-License-Identifier: MIT
! Multidimensional extrema tests through the public SciFort signal interface.
program test_signal_extrema_nd
    use scifort_signal, only : dp, argrelmax_mask_into, argrelmin_mask_into, &
        argrelextrema_mask_into, extrema_status_success, extrema_status_invalid_input
    implicit none
    real(dp) :: x(5), x2(2, 5), x15(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5)
    real(dp) :: matrix(5, 8)
    logical :: mask(5), m2(2, 5), m15(2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5)
    logical :: output(5, 8), bad_shape(4), mzero(0, 3)
    real(dp) :: zero(0, 3)
    integer :: code, i

    x = [0.0_dp, 3.0_dp, 1.0_dp, 2.0_dp, 0.0_dp]
    mask = .true.
    call argrelmax_mask_into(x, mask, status=code)
    call require(code == extrema_status_success, '1-D status')
    call require(all(mask .eqv. [.false.,.true.,.false.,.true.,.false.]), '1-D maxima')
    call argrelmin_mask_into(x, mask, order=2, mode='wrap', status=code)
    call require(code == extrema_status_success, '1-D wrapped status')
    call require(all(mask .eqv. [.false.,.false.,.false.,.false.,.false.]), '1-D minima')

    ! A non-strict comparator deliberately admits flat plateaus.
    x = [1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    call argrelextrema_mask_into(x, mask, greater_equal, status=code)
    call require(code == extrema_status_success, 'custom comparator status')
    call require(all(mask .eqv. [.true.,.true.,.false.,.true.,.true.]), 'plateaus')

    x2(1,:) = [0.0_dp, 4.0_dp, 1.0_dp, 3.0_dp, 0.0_dp]
    x2(2,:) = [0.0_dp, 1.0_dp, 4.0_dp, 1.0_dp, 0.0_dp]
    call argrelmax_mask_into(x2, m2, axis=-1, status=code)
    call require(code == extrema_status_success, '2-D status')
    call require(all(m2(1,:) .eqv. [.false.,.true.,.false.,.true.,.false.]), '2-D row 1')
    call require(all(m2(2,:) .eqv. [.false.,.false.,.true.,.false.,.false.]), '2-D row 2')
    call argrelmax_mask_into(x2, m2, axis=1, status=code)
    call require(code == extrema_status_success, '2-D default axis')
    call require(.not. any(m2), '2-D axis 1 endpoint clipping')

    x15 = reshape([0.0_dp, 0.0_dp, 4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, &
        2.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], shape(x15))
    call argrelmax_mask_into(x15, m15, axis=15, status=code)
    call require(code == extrema_status_success, 'rank-15 status')
    call require(count(m15) == 3, 'rank-15 axis 15 extrema count')

    do i = 1, 8
        matrix(:,i) = real(i,dp) * x
    end do
    output = .true.
    call argrelmax_mask_into(matrix(1:5:2,1:8:2), output(1:5:2,1:8:2), &
        axis=2, status=code)
    call require(code == extrema_status_success, 'strided status')
    call require(all(output(2:4:2,:)), 'strided unrelated rows unchanged')
    call require(all(output(:,2:8:2)), 'strided unrelated columns unchanged')
    call require(.not. any(output(1:5:2,1:8:2)), 'strided axis mask')

    mask = .true.
    call argrelmax_mask_into(x, mask, order=0, status=code)
    call require(code == extrema_status_invalid_input, 'invalid order')
    call require(all(mask), 'unchanged mask after invalid order')
    call argrelmax_mask_into(x, mask, axis=2, status=code)
    call require(code == extrema_status_invalid_input, 'invalid axis')
    call require(all(mask), 'unchanged mask after invalid axis')
    call argrelmax_mask_into(x, mask, mode='reflect', status=code)
    call require(code == extrema_status_invalid_input, 'invalid mode')
    call require(all(mask), 'unchanged mask after invalid mode')
    bad_shape = .true.
    call argrelmax_mask_into(x, bad_shape, status=code)
    call require(code == extrema_status_invalid_input, 'invalid output shape')
    call require(all(bad_shape), 'unchanged mask after invalid shape')

    call argrelmin_mask_into(zero, mzero, axis=2, status=code)
    call require(code == extrema_status_success, 'empty input succeeds')
    call require(size(mzero) == 0, 'empty shape preserved')
    print '(a)', 'PASS multidimensional relative-extrema tests'
contains
    pure logical function greater_equal(a, b)
        real(dp), intent(in) :: a !! candidate amplitude
        real(dp), intent(in) :: b !! comparison amplitude
        greater_equal = a >= b
    end function greater_equal

    subroutine require(ok, why)
        logical, intent(in) :: ok !! assertion condition
        character(*), intent(in) :: why !! failure description
        if (.not. ok) then
            print '(a)', 'FAILED: '//why
            error stop 1
        end if
    end subroutine require
end program test_signal_extrema_nd
