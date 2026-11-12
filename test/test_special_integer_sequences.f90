program test_special_integer_sequences
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_special, only : bernoulli, dp, euler, &
        integer_sequences_status_invalid_argument, integer_sequences_status_success, stirling2
    implicit none

    real(dp), allocatable :: b(:)
    real(dp), allocatable :: e(:)
    real(dp) :: array_values(3)
    integer :: status

    call bernoulli(20, b, status)
    call require(status == integer_sequences_status_success, 'bernoulli status')
    call require(lbound(b, 1) == 0 .and. ubound(b, 1) == 20, 'bernoulli bounds')
    call check_close(b(0), 1.0_dp, 0.0_dp, 'B0')
    call check_close(b(1), -0.5_dp, 0.0_dp, 'B1')
    call check_close(b(2), 1.0_dp / 6.0_dp, 8.0e-15_dp, 'B2')
    call check_close(b(4), -1.0_dp / 30.0_dp, 8.0e-15_dp, 'B4')
    call check_close(b(10), 5.0_dp / 66.0_dp, 1.0e-14_dp, 'B10')
    call check_close(b(20), -529.124242424242424242424242424242_dp, 2.0e-14_dp, 'B20')
    call require(maxval(abs(b(3:19:2))) <= tiny(1.0_dp), 'odd Bernoulli numbers')

    call bernoulli(100, b, status)
    call check_close(b(50), 7.5008667460769643668557200757576e24_dp, 3.0e-14_dp, 'B50')
    call check_close(b(100), -2.8382249570693706959264156336482e78_dp, 8.0e-14_dp, 'B100')

    call bernoulli(-1, b, status)
    call require(status == integer_sequences_status_invalid_argument, 'bernoulli invalid status')
    call require(size(b) == 0, 'bernoulli invalid size')

    call euler(22, e, status)
    call require(status == integer_sequences_status_success, 'euler status')
    call require(lbound(e, 1) == 0 .and. ubound(e, 1) == 22, 'euler bounds')
    call check_close(e(0), 1.0_dp, 0.0_dp, 'E0')
    call check_close(e(2), -1.0_dp, 0.0_dp, 'E2')
    call check_close(e(4), 5.0_dp, 0.0_dp, 'E4')
    call check_close(e(10), -50521.0_dp, 0.0_dp, 'E10')
    call check_close(e(20), 370371188237525.0_dp, 0.0_dp, 'E20')
    call check_close(e(22), -69348874393137901.0_dp, 2.0e-15_dp, 'E22')
    call require(maxval(abs(e(1:21:2))) <= tiny(1.0_dp), 'odd Euler numbers')

    call euler(100, e, status)
    call check_close(e(50), -6.0532852481886218963143837851116e54_dp, 2.0e-14_dp, 'E50')
    call check_close(e(100), 2.9035283466610974970546038347644e138_dp, 3.0e-14_dp, 'E100')

    call euler(188, e, status)
    call require(.not. ieee_is_finite(e(188)) .and. e(188) > 0.0_dp, 'E188 overflow sign')

    call euler(-1, e, status)
    call require(status == integer_sequences_status_invalid_argument, 'euler invalid status')
    call require(size(e) == 0, 'euler invalid size')

    call check_close(stirling2(0, 0), 1.0_dp, 0.0_dp, 'S00')
    call check_close(stirling2(10, 3), 9330.0_dp, 0.0_dp, 'S10,3')
    call check_close(stirling2(50, 25), 7.4538021532732000833796262348376e42_dp, &
        3.0e-15_dp, 'S50,25')
    call check_close(stirling2(100, 50), 4.3098323700936634042151430154726e101_dp, &
        6.0e-15_dp, 'S100,50')
    call require(abs(stirling2(-1, 0)) <= tiny(1.0_dp), 'negative n')
    call require(abs(stirling2(10, -1)) <= tiny(1.0_dp), 'negative k')
    call require(abs(stirling2(5, 6)) <= tiny(1.0_dp), 'k greater than n')

    array_values = stirling2([10, 9, 8], [3, 3, 2])
    call require(maxval(abs(array_values - [9330.0_dp, 3025.0_dp, 127.0_dp])) <= tiny(1.0_dp), &
        'stirling elemental')

    print *, 'test_special_integer_sequences: PASS'

contains

    subroutine require(condition, label)
        logical, intent(in) :: condition !! condition that must hold
        character(len=*), intent(in) :: label !! failure label

        if (.not. condition) then
            print *, 'FAIL: ', trim(label)
            error stop 1
        end if
    end subroutine require

    subroutine check_close(actual, expected, relative_tolerance, label)
        real(dp), intent(in) :: actual !! computed value
        real(dp), intent(in) :: expected !! reference value
        real(dp), intent(in) :: relative_tolerance !! allowed relative error
        character(len=*), intent(in) :: label !! failure label
        real(dp) :: scale

        scale = max(1.0_dp, abs(expected))
        if (abs(actual - expected) > relative_tolerance * scale) then
            print *, 'FAIL: ', trim(label)
            print *, ' actual  = ', actual
            print *, ' expected= ', expected
            error stop 1
        end if
    end subroutine check_close

end program test_special_integer_sequences
