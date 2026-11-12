! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Two-dimensional convolution and correlation reference and edge cases.
program test_signal_convolve2d
    use scifort_signal, only : dp, convolve2d, correlate2d, &
        convolution_status_success, convolution_status_invalid_input
    implicit none

    real(dp) :: a(2,2), b(2,2)
    real(dp), allocatable :: y(:,:)
    complex(dp) :: c(2,2), d(2,2)
    complex(dp), allocatable :: z(:,:)
    integer :: status, i, j

    a(:,1) = [1.0_dp, 3.0_dp]
    a(:,2) = [2.0_dp, 4.0_dp]
    b(:,1) = [5.0_dp, 7.0_dp]
    b(:,2) = [6.0_dp, 8.0_dp]

    ! SciPy convolve2d reference for 2-by-2 matrices.
    y = convolve2d(a,b)
    call expect_real(y,reshape([5.0_dp,22.0_dp,21.0_dp,16.0_dp,60.0_dp,52.0_dp, &
        12.0_dp,40.0_dp,32.0_dp],[3,3]))
    y = convolve2d(a,b,mode='same')
    call expect_real(y,reshape([5.0_dp,22.0_dp,16.0_dp,60.0_dp],[2,2]))
    y = convolve2d(a,b,mode='valid')
    call expect_real(y,reshape([60.0_dp],[1,1]))

    y = correlate2d(a,b,mode='full')
    call expect_real(y,reshape([8.0_dp,30.0_dp,18.0_dp,23.0_dp,70.0_dp,39.0_dp, &
        14.0_dp,38.0_dp,20.0_dp],[3,3]))
    y = correlate2d(a,b,mode='same')
    call expect_real(y,reshape([70.0_dp,39.0_dp,38.0_dp,20.0_dp],[2,2]))

    y = convolve2d(a,b,boundary='wrap')
    call expect_real(y,reshape([70.0_dp,62.0_dp,70.0_dp,68.0_dp,60.0_dp,68.0_dp, &
        70.0_dp,62.0_dp,70.0_dp],[3,3]))
    y = correlate2d(a,b,boundary='symm',mode='same')
    call expect_real(y,reshape([70.0_dp,92.0_dp,82.0_dp,104.0_dp],[2,2]))

    c = cmplx(a, a/3.0_dp, dp)
    d = cmplx(b, -b/5.0_dp, dp)
    z = correlate2d(c,d)
    call expect_complex(z, convolve2d(c,conjg(d(2:1:-1,2:1:-1))))
    z = convolve2d(a,d)
    call expect_complex(z,convolve2d(cmplx(a,kind=dp),d))
    z = correlate2d(c,b)
    call expect_complex(z,correlate2d(c,cmplx(b,kind=dp)))

    ! Correlation with a larger kernel ignores nonzero fill and symm/wrap
    ! extensions, and must preserve the direction of the correlation.
    block
        real(dp) :: big(3,3)
        big(:,1) = [5.0_dp,8.0_dp,11.0_dp]
        big(:,2) = [6.0_dp,9.0_dp,12.0_dp]
        big(:,3) = [7.0_dp,10.0_dp,13.0_dp]
        y = correlate2d(a,big,mode='valid',boundary='wrap',fillvalue=17.0_dp)
        call expect_real(y,reshape([117.0_dp,87.0_dp,107.0_dp,77.0_dp],[2,2]))
        y = correlate2d(a,big,mode='valid',boundary='symm',fillvalue=-4.0_dp)
        call expect_real(y,reshape([117.0_dp,87.0_dp,107.0_dp,77.0_dp],[2,2]))
        z = correlate2d(c,cmplx(big,kind=dp),mode='valid', &
            boundary='wrap',fillvalue=(1.0_dp,2.0_dp),status=status)
        if (status /= convolution_status_success) error stop 'complex valid status'
        ! For a larger second operand, SciPy conjugates the result of
        ! the mathematical full-correlation slice in valid mode.
        block
            complex(dp), allocatable :: full(:,:)
            full = correlate2d(c,cmplx(big,kind=dp),mode='full')
            call expect_complex(z,conjg(full(2:3,2:3)))
        end block
    end block

    ! For valid output, boundary extension and fill values have no effect.
    y = convolve2d(a,b,mode='valid',boundary='symm',fillvalue=4.0_dp)
    call expect_real(y,reshape([60.0_dp],[1,1]))
    y = convolve2d(a,b,mode='none',status=status)
    if (status /= convolution_status_invalid_input .or. size(y) /= 0) error stop 'invalid mode'
    y = correlate2d(a,b,boundary='bad',status=status)
    if (status /= convolution_status_invalid_input .or. size(y) /= 0) error stop 'invalid boundary'
    y = convolve2d(a,reshape([real(dp) :: ],[0,1]),status=status)
    if (status /= convolution_status_invalid_input .or. size(y) /= 0) error stop 'empty kernel'
    y = convolve2d(reshape([1.0_dp,2.0_dp],[2,1]),reshape([1.0_dp,2.0_dp],[1,2]), &
        mode='valid',status=status)
    if (status /= convolution_status_invalid_input .or. size(y) /= 0) error stop 'mixed valid shapes'

    ! Exhaustive one-by-one kernels are useful for checking all boundaries.
    do i = 1, 5
        do j = 1, 5
            block
                real(dp) :: x(i,j)
                x = 1.0_dp
                y = convolve2d(x,reshape([2.0_dp],[1,1]),mode='same')
                if (any(shape(y) /= [i,j]) .or. any(abs(y - 2.0_dp) > 1.e-12_dp)) &
                    error stop 'unit kernel mismatch'
            end block
        end do
    end do
    print '(a)', 'test_signal_convolve2d passed'

contains

    subroutine expect_real(actual, expected)
        real(dp), intent(in) :: actual(:,:) !! computed real matrix
        real(dp), intent(in) :: expected(:,:) !! reference real matrix
        if (any(shape(actual) /= shape(expected))) error stop 'real shape mismatch'
        if (any(abs(actual - expected) > 1.e-12_dp)) then
            print *, 'actual:', actual, 'expected:', expected
            error stop 'real value mismatch'
        end if
    end subroutine expect_real

    subroutine expect_complex(actual, expected)
        complex(dp), intent(in) :: actual(:,:) !! computed complex matrix
        complex(dp), intent(in) :: expected(:,:) !! reference complex matrix
        if (any(shape(actual) /= shape(expected))) error stop 'complex shape mismatch'
        if (any(abs(actual - expected) > 1.e-12_dp)) error stop 'complex value mismatch'
    end subroutine expect_complex
end program test_signal_convolve2d
