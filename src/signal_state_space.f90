! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Continuous- or discrete-time LTI transfer-function/state-space conversions.
! The algebra is identical for s and z: polynomial coefficients are in
! descending powers. This module does not choose a sampling interval.
module scifort_signal_state_space
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_zpk_conversion, only : zpk2tf, tf2zpk, &
        filter_tf_result, filter_zpk_result, zpk_conversion_success, &
        zpk_conversion_invalid_input
    implicit none
    private

    integer, parameter, public :: state_space_success = 0
    integer, parameter, public :: state_space_invalid_input = 1
    integer, parameter, public :: state_space_numerical_failure = 2

    type, public :: state_space_result
        complex(dp), allocatable :: a(:,:) !! state matrix, (n,n)
        complex(dp), allocatable :: b(:,:) !! input matrix, (n,nin)
        complex(dp), allocatable :: c(:,:) !! output matrix, (nout,n)
        complex(dp), allocatable :: d(:,:) !! direct feedthrough, (nout,nin)
        integer :: status = state_space_invalid_input !! conversion status
    end type state_space_result

    type, public :: state_space_tf_result
        complex(dp), allocatable :: b(:,:) !! one polynomial per output, descending powers
        complex(dp), allocatable :: a(:) !! common denominator, descending powers
        integer :: status = state_space_invalid_input !! conversion status
    end type state_space_tf_result

    interface tf2ss
        module procedure tf2ss_real, tf2ss_complex
        module procedure tf2ss_real_multi, tf2ss_complex_multi
    end interface tf2ss
    interface ss2tf
        module procedure ss2tf_real, ss2tf_complex
    end interface ss2tf
    interface zpk2ss
        module procedure zpk2ss_real, zpk2ss_real_complex_gain
        module procedure zpk2ss_complex_real_gain, zpk2ss_complex
    end interface zpk2ss
    interface ss2zpk
        module procedure ss2zpk_real, ss2zpk_complex
    end interface ss2zpk

    public :: tf2ss, ss2tf, zpk2ss, ss2zpk

contains

    pure logical function finite_vector(v) result(ok)
        complex(dp), intent(in) :: v(:) !! complex values
        ok = all(ieee_is_finite(real(v,dp))) .and. all(ieee_is_finite(aimag(v)))
    end function finite_vector

    pure logical function finite_matrix(v) result(ok)
        complex(dp), intent(in) :: v(:,:) !! complex matrix
        ok = all(ieee_is_finite(real(v,dp))) .and. all(ieee_is_finite(aimag(v)))
    end function finite_matrix

    function tf2ss_real(b, a) result(out)
        real(dp), intent(in) :: b(:) !! numerator, descending powers
        real(dp), intent(in) :: a(:) !! denominator, descending powers
        type(state_space_result) :: out
        out = tf2ss_complex(cmplx(b,0.0_dp,dp), cmplx(a,0.0_dp,dp))
    end function tf2ss_real

    function tf2ss_complex(b, a) result(out)
        complex(dp), intent(in) :: b(:) !! numerator, descending powers
        complex(dp), intent(in) :: a(:) !! denominator, descending powers
        type(state_space_result) :: out
        complex(dp) :: numer(1,size(b))
        numer(1,:) = b
        out = tf2ss_complex_multi(numer,a)
    end function tf2ss_complex

    function tf2ss_real_multi(b, a) result(out)
        real(dp), intent(in) :: b(:,:) !! numerator rows for each output
        real(dp), intent(in) :: a(:) !! common denominator
        type(state_space_result) :: out
        out = tf2ss_complex_multi(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp))
    end function tf2ss_real_multi

    function tf2ss_complex_multi(b, a) result(out)
        complex(dp), intent(in) :: b(:,:) !! numerator rows for each output
        complex(dp), intent(in) :: a(:) !! denominator polynomial
        type(state_space_result) :: out
        complex(dp), allocatable :: den(:), numer(:,:)
        integer :: i, j, first_a, first_b, n, nout, nb

        nout = size(b,1)
        if (nout == 0 .or. size(b,2) == 0 .or. size(a) == 0) return
        if (.not. finite_matrix(b) .or. .not. finite_vector(a)) return
        first_a = 1
        do while (first_a <= size(a))
            if (abs(a(first_a)) > 0.0_dp) exit
            first_a = first_a+1
        end do
        if (first_a > size(a)) return
        den = a(first_a:)/a(first_a)
        if (.not. finite_vector(den)) then
            out%status = state_space_numerical_failure
            return
        end if
        first_b = 1
        do while (first_b < size(b,2))
            if (any(abs(b(:,first_b)) > 0.0_dp)) exit
            first_b = first_b+1
        end do
        nb = size(b,2)-first_b+1
        if (nb > size(den)) return  ! improper rational transfer function
        allocate(numer(nout,size(den)))
        numer = (0.0_dp,0.0_dp)
        numer(:,size(den)-nb+1:) = b(:,first_b:)/a(first_a)
        if (.not. finite_matrix(numer)) then
            out%status = state_space_numerical_failure
            return
        end if

        n = size(den)-1
        if (n == 0) then
            ! Match SciPy's static-gain representation: one dummy zero state.
            allocate(out%a(1,1),out%b(1,1),out%c(nout,1),out%d(nout,1))
            out%a = (0.0_dp,0.0_dp)
            out%b = (0.0_dp,0.0_dp)
            out%c = (0.0_dp,0.0_dp)
            out%d(:,1) = numer(:,1)
        else
            allocate(out%a(n,n),out%b(n,1),out%c(nout,n),out%d(nout,1))
            out%a = (0.0_dp,0.0_dp)
            out%b = (0.0_dp,0.0_dp)
            out%a(1,:) = -den(2:)
            do i = 2, n
                out%a(i,i-1) = (1.0_dp,0.0_dp)
            end do
            out%b(1,1) = (1.0_dp,0.0_dp)
            out%d(:,1) = numer(:,1)
            do j = 1, n
                out%c(:,j) = numer(:,j+1) - numer(:,1)*den(j+1)
            end do
        end if
        if (.not. finite_matrix(out%a) .or. .not. finite_matrix(out%c)) then
            out%status = state_space_numerical_failure
            return
        end if
        out%status = state_space_success
    end function tf2ss_complex_multi

    function ss2tf_real(a, b, c, d, input) result(out)
        real(dp), intent(in) :: a(:,:) !! state matrix
        real(dp), intent(in) :: b(:,:) !! input matrix
        real(dp), intent(in) :: c(:,:) !! output matrix
        real(dp), intent(in) :: d(:,:) !! feedthrough matrix
        integer, intent(in), optional :: input !! one-based selected input, default 1
        type(state_space_tf_result) :: out
        out = ss2tf_complex(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp), &
            cmplx(c,0.0_dp,dp),cmplx(d,0.0_dp,dp),input)
    end function ss2tf_real

    function ss2tf_complex(a, b, c, d, input) result(out)
        ! Faddeev-LeVerrier recursion produces characteristic coefficients
        ! and adjugate factors together; no matrix inversion is required.
        complex(dp), intent(in) :: a(:,:) !! state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! feedthrough matrix
        integer, intent(in), optional :: input !! one-based selected input, default 1
        type(state_space_tf_result) :: out
        complex(dp), allocatable :: factors(:,:), product(:,:), bi(:)
        complex(dp) :: coefficient
        integer :: n, nout, nin, selected, k, j

        n = size(a,1)
        nout = size(c,1)
        nin = size(b,2)
        selected = 1
        if (present(input)) selected = input
        if (size(a,2) /= n .or. size(b,1) /= n .or. size(c,2) /= n) return
        if (nout == 0 .or. nin == 0 .or. selected < 1 .or. selected > nin) return
        if (size(d,1) /= nout .or. size(d,2) /= nin) return
        if (.not. finite_matrix(a) .or. .not. finite_matrix(b) .or. &
            .not. finite_matrix(c) .or. .not. finite_matrix(d)) return

        allocate(out%a(n+1),out%b(nout,n+1))
        out%a(1) = (1.0_dp,0.0_dp)
        out%b(:,1) = d(:,selected)
        if (n > 0) then
            allocate(factors(n,n),product(n,n),bi(n))
            factors = (0.0_dp,0.0_dp)
            do j = 1, n
                factors(j,j) = (1.0_dp,0.0_dp)
            end do
            bi = b(:,selected)
            do k = 1, n
                product = matmul(a,factors)
                coefficient = (0.0_dp,0.0_dp)
                do j = 1, n
                    coefficient = coefficient-product(j,j)
                end do
                coefficient = coefficient/real(k,dp)
                out%a(k+1) = coefficient
                out%b(:,k+1) = matmul(c,matmul(factors,bi)) + &
                    d(:,selected)*coefficient
                factors = product
                do j = 1, n
                    factors(j,j) = factors(j,j)+coefficient
                end do
            end do
        end if
        if (.not. finite_vector(out%a) .or. .not. finite_matrix(out%b)) then
            out%status = state_space_numerical_failure
            return
        end if
        out%status = state_space_success
    end function ss2tf_complex

    function zpk2ss_real(z, p, k) result(out)
        real(dp), intent(in) :: z(:) !! real zeros
        real(dp), intent(in) :: p(:) !! real poles
        real(dp), intent(in) :: k !! real gain
        type(state_space_result) :: out
        out = zpk2ss_complex(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),cmplx(k,0.0_dp,dp))
    end function zpk2ss_real

    function zpk2ss_real_complex_gain(z, p, k) result(out)
        real(dp), intent(in) :: z(:) !! real zeros
        real(dp), intent(in) :: p(:) !! real poles
        complex(dp), intent(in) :: k !! complex gain
        type(state_space_result) :: out
        out = zpk2ss_complex(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),k)
    end function zpk2ss_real_complex_gain

    function zpk2ss_complex_real_gain(z, p, k) result(out)
        complex(dp), intent(in) :: z(:) !! complex zeros
        complex(dp), intent(in) :: p(:) !! complex poles
        real(dp), intent(in) :: k !! real gain
        type(state_space_result) :: out
        out = zpk2ss_complex(z,p,cmplx(k,0.0_dp,dp))
    end function zpk2ss_complex_real_gain

    function zpk2ss_complex(z, p, k) result(out)
        complex(dp), intent(in) :: z(:) !! zeros
        complex(dp), intent(in) :: p(:) !! poles
        complex(dp), intent(in) :: k !! gain
        type(state_space_result) :: out
        type(filter_tf_result) :: tf
        tf = zpk2tf(z,p,k)
        if (tf%status /= zpk_conversion_success) then
            if (tf%status /= zpk_conversion_invalid_input) &
                out%status = state_space_numerical_failure
            return
        end if
        out = tf2ss_complex(tf%b,tf%a)
    end function zpk2ss_complex

    function ss2zpk_real(a, b, c, d, input, output) result(out)
        real(dp), intent(in) :: a(:,:) !! state matrix
        real(dp), intent(in) :: b(:,:) !! input matrix
        real(dp), intent(in) :: c(:,:) !! output matrix
        real(dp), intent(in) :: d(:,:) !! feedthrough matrix
        integer, intent(in), optional :: input !! one-based input, default 1
        integer, intent(in), optional :: output !! one-based output, default 1
        type(filter_zpk_result) :: out
        out = ss2zpk_complex(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp), &
            cmplx(c,0.0_dp,dp),cmplx(d,0.0_dp,dp),input,output)
    end function ss2zpk_real

    function ss2zpk_complex(a, b, c, d, input, output) result(out)
        complex(dp), intent(in) :: a(:,:) !! state matrix
        complex(dp), intent(in) :: b(:,:) !! input matrix
        complex(dp), intent(in) :: c(:,:) !! output matrix
        complex(dp), intent(in) :: d(:,:) !! feedthrough matrix
        integer, intent(in), optional :: input !! one-based input, default 1
        integer, intent(in), optional :: output !! one-based output, default 1
        type(filter_zpk_result) :: out
        type(state_space_tf_result) :: tf
        integer :: selected

        selected = 1
        if (present(output)) selected = output
        if (selected < 1 .or. selected > size(c,1)) return
        tf = ss2tf_complex(a,b,c,d,input)
        if (tf%status /= state_space_success) return
        out = tf2zpk(tf%b(selected,:),tf%a)
    end function ss2zpk_complex

end module scifort_signal_state_space
