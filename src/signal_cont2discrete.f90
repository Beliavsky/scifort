! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Continuous-to-discrete LTI conversion (state-space, TF, and ZPK).
! ZOH/FOH use an augmented matrix exponential with a scaled Pade(13)
! approximant. Generalized bilinear methods use pivoted linear solves.
module scifort_signal_cont2discrete
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_state_space, only : tf2ss, ss2tf, zpk2ss, ss2zpk, &
        state_space_result, state_space_tf_result, state_space_success
    use scifort_signal_zpk_conversion, only : filter_zpk_result, zpk_conversion_success
    implicit none
    private

    integer, parameter, public :: cont2discrete_success = 0
    integer, parameter, public :: cont2discrete_invalid_input = 1
    integer, parameter, public :: cont2discrete_numerical_failure = 2

    type, public :: discrete_state_space_result
        complex(dp), allocatable :: a(:,:) !! discrete state transition
        complex(dp), allocatable :: b(:,:) !! discrete input matrix
        complex(dp), allocatable :: c(:,:) !! discrete output matrix
        complex(dp), allocatable :: d(:,:) !! discrete feedthrough matrix
        real(dp) :: dt = 0.0_dp !! sampling interval
        integer :: status = cont2discrete_invalid_input !! conversion status
    end type discrete_state_space_result

    type, public :: discrete_tf_result
        complex(dp), allocatable :: b(:,:) !! one numerator per output, descending powers of z
        complex(dp), allocatable :: a(:) !! common denominator, descending powers of z
        real(dp) :: dt = 0.0_dp !! sampling interval
        integer :: status = cont2discrete_invalid_input !! conversion status
    end type discrete_tf_result

    type, public :: discrete_zpk_result
        complex(dp), allocatable :: zeros(:) !! discrete zeros
        complex(dp), allocatable :: poles(:) !! discrete poles
        complex(dp) :: gain = (0.0_dp,0.0_dp) !! discrete gain
        real(dp) :: dt = 0.0_dp !! sampling interval
        integer :: status = cont2discrete_invalid_input !! conversion status
    end type discrete_zpk_result

    interface cont2discrete
        module procedure cont2discrete_ss_real, cont2discrete_ss_complex
        module procedure cont2discrete_tf_real, cont2discrete_tf_complex
        module procedure cont2discrete_tf_real_multi, cont2discrete_tf_complex_multi
    end interface cont2discrete
    interface cont2discrete_zpk
        module procedure cont2discrete_zpk_real, cont2discrete_zpk_complex
        module procedure cont2discrete_zpk_real_complex_gain
        module procedure cont2discrete_zpk_complex_real_gain
    end interface cont2discrete_zpk
    ! Shared internal matrix exponential for continuous-time propagation.
    public :: cont2discrete, cont2discrete_zpk, matrix_exp_pade

contains

    pure logical function finite_matrix(a) result(ok)
        complex(dp), intent(in) :: a(:,:) !! matrix to check
        ok = all(ieee_is_finite(real(a,dp))) .and. all(ieee_is_finite(aimag(a)))
    end function finite_matrix

    pure function lower_case(value) result(out)
        character(*), intent(in) :: value !! case-insensitive method name
        character(len(value)) :: out
        integer :: j, code
        out = value
        do j = 1, len(out)
            code = iachar(out(j:j))
            if (code >= iachar('A') .and. code <= iachar('Z')) out(j:j) = achar(code+32)
        end do
    end function lower_case

    function cont2discrete_ss_real(a,b,c,d,dt,method,alpha) result(out)
        real(dp), intent(in) :: a(:,:) !! continuous state matrix
        real(dp), intent(in) :: b(:,:) !! continuous input matrix
        real(dp), intent(in) :: c(:,:) !! continuous output matrix
        real(dp), intent(in) :: d(:,:) !! continuous feedthrough matrix
        real(dp), intent(in) :: dt !! positive sampling interval
        character(*), intent(in), optional :: method !! zoh (default), foh, gbt, bilinear, euler, backward_diff, impulse
        real(dp), intent(in), optional :: alpha !! gbt weight, in [0,1]
        type(discrete_state_space_result) :: out
        out = cont2discrete_ss_complex(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp), &
            cmplx(c,0.0_dp,dp),cmplx(d,0.0_dp,dp),dt,method,alpha)
    end function cont2discrete_ss_real

    function cont2discrete_ss_complex(a,b,c,d,dt,method,alpha) result(out)
        complex(dp), intent(in) :: a(:,:) !! continuous state matrix, (n,n)
        complex(dp), intent(in) :: b(:,:) !! input matrix, (n,nin)
        complex(dp), intent(in) :: c(:,:) !! output matrix, (nout,n)
        complex(dp), intent(in) :: d(:,:) !! feedthrough matrix, (nout,nin)
        real(dp), intent(in) :: dt !! positive sampling interval
        character(*), intent(in), optional :: method !! conversion method
        real(dp), intent(in), optional :: alpha !! gbt weight
        type(discrete_state_space_result) :: out
        complex(dp), allocatable :: e(:,:), ex(:,:), w(:,:), rhs(:,:), solution(:,:), cd_trans(:,:)
        real(dp) :: weight
        character(:), allocatable :: selected
        integer :: n, nin, nout, j, info

        out%dt = dt
        n = size(a,1)
        nin = size(b,2)
        nout = size(c,1)
        if (.not. ieee_is_finite(dt)) return
        if (dt <= 0.0_dp) return
        if (size(a,2) /= n .or. size(b,1) /= n .or. size(c,2) /= n) return
        if (size(d,1) /= nout .or. size(d,2) /= nin) return
        if (nin < 1 .or. nout < 1) return
        if (.not. finite_matrix(a) .or. .not. finite_matrix(b) .or. &
            .not. finite_matrix(c) .or. .not. finite_matrix(d)) return
        selected = 'zoh'
        if (present(method)) selected = trim(adjustl(lower_case(method)))
        weight = 0.0_dp
        select case(selected)
        case('gbt')
            if (.not. present(alpha)) return
            weight = alpha
            if (.not. ieee_is_finite(weight)) return
            if (weight < 0.0_dp .or. weight > 1.0_dp) return
        case('bilinear','tustin')
            weight = 0.5_dp
        case('euler','forward_diff')
            weight = 0.0_dp
        case('backward_diff')
            weight = 1.0_dp
        case('zoh','foh','impulse')
            continue
        case default
            return
        end select

        allocate(out%a(n,n),out%b(n,nin),out%c(nout,n),out%d(nout,nin))
        if (selected == 'gbt' .or. selected == 'bilinear' .or. selected == 'tustin' .or. &
            selected == 'euler' .or. selected == 'forward_diff' .or. selected == 'backward_diff') then
            ! (I - alpha*dt*A)^{-1} times (I + (1-alpha)*dt*A, dt*B).
            allocate(w(n,n),rhs(n,n+nin))
            w = -weight*dt*a
            rhs(:,1:n) = (1.0_dp-weight)*dt*a
            rhs(:,n+1:) = dt*b
            do j = 1, n
                w(j,j) = w(j,j)+1.0_dp
                rhs(j,j) = rhs(j,j)+1.0_dp
            end do
            call solve_matrix(w,rhs,solution,info)
            if (info /= 0) then
                out%status = cont2discrete_numerical_failure
                return
            end if
            out%a = solution(:,1:n)
            out%b = solution(:,n+1:)
            call solve_matrix(transpose(w),transpose(c),cd_trans,info)
            if (info /= 0) then
                out%status = cont2discrete_numerical_failure
                return
            end if
            out%c = transpose(cd_trans)
            out%d = d+weight*matmul(c,out%b)
        else if (selected == 'zoh') then
            allocate(e(n+nin,n+nin))
            e = (0.0_dp,0.0_dp)
            e(1:n,1:n) = dt*a
            e(1:n,n+1:) = dt*b
            call matrix_exp_pade(e,ex,info)
            if (info /= 0) then
                out%status = cont2discrete_numerical_failure
                return
            end if
            out%a = ex(1:n,1:n)
            out%b = ex(1:n,n+1:)
            out%c = c
            out%d = d
        else if (selected == 'foh') then
            ! Augmented exponential also integrates the input ramp exactly.
            allocate(e(n+2*nin,n+2*nin))
            e = (0.0_dp,0.0_dp)
            e(1:n,1:n) = dt*a
            e(1:n,n+1:n+nin) = dt*b
            do j = 1, nin
                e(n+j,n+nin+j) = (1.0_dp,0.0_dp)
            end do
            call matrix_exp_pade(e,ex,info)
            if (info /= 0) then
                out%status = cont2discrete_numerical_failure
                return
            end if
            out%a = ex(1:n,1:n)
            out%b = ex(1:n,n+1:n+nin)-ex(1:n,n+nin+1:)+ &
                matmul(out%a,ex(1:n,n+nin+1:))
            out%c = c
            out%d = d+matmul(c,ex(1:n,n+nin+1:))
        else
            ! Impulse-invariant discretization is defined only for a strictly
            ! proper system (all direct-feedthrough coefficients are zero).
            if (any(abs(d)>0.0_dp)) return
            call matrix_exp_pade(dt*a,ex,info)
            if (info /= 0) then
                out%status = cont2discrete_numerical_failure
                return
            end if
            out%a = ex
            out%b = matmul(ex,b)*dt
            out%c = c
            out%d = matmul(c,b)*dt
        end if
        if (.not. finite_matrix(out%a) .or. .not. finite_matrix(out%b) .or. &
            .not. finite_matrix(out%c) .or. .not. finite_matrix(out%d)) then
            out%status = cont2discrete_numerical_failure
            return
        end if
        out%status = cont2discrete_success
    end function cont2discrete_ss_complex

    function cont2discrete_tf_real(b,a,dt,method,alpha) result(out)
        real(dp), intent(in) :: b(:) !! numerator, descending powers of s
        real(dp), intent(in) :: a(:) !! denominator, descending powers of s
        real(dp), intent(in) :: dt !! positive sampling interval
        character(*), intent(in), optional :: method !! conversion method
        real(dp), intent(in), optional :: alpha !! gbt weight
        type(discrete_tf_result) :: out
        out = cont2discrete_tf_complex(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),dt,method,alpha)
    end function cont2discrete_tf_real

    function cont2discrete_tf_complex(b,a,dt,method,alpha) result(out)
        complex(dp), intent(in) :: b(:) !! numerator, descending powers of s
        complex(dp), intent(in) :: a(:) !! denominator, descending powers of s
        real(dp), intent(in) :: dt !! positive sampling interval
        character(*), intent(in), optional :: method !! conversion method
        real(dp), intent(in), optional :: alpha !! gbt weight
        type(discrete_tf_result) :: out
        type(state_space_result) :: ss
        type(discrete_state_space_result) :: sd
        type(state_space_tf_result) :: tf
        out%dt = dt
        ss = tf2ss(b,a)
        if (ss%status /= state_space_success) return
        sd = cont2discrete_ss_complex(ss%a,ss%b,ss%c,ss%d,dt,method,alpha)
        if (sd%status /= cont2discrete_success) then
            out%status = sd%status
            return
        end if
        tf = ss2tf(sd%a,sd%b,sd%c,sd%d)
        if (tf%status /= state_space_success) then
            out%status = cont2discrete_numerical_failure
            return
        end if
        out%b = tf%b
        out%a = tf%a
        out%status = cont2discrete_success
    end function cont2discrete_tf_complex

    function cont2discrete_tf_real_multi(b,a,dt,method,alpha) result(out)
        real(dp), intent(in) :: b(:,:) !! numerator rows for single-input/multiple-output system
        real(dp), intent(in) :: a(:) !! common denominator, descending powers of s
        real(dp), intent(in) :: dt !! positive sampling interval
        character(*), intent(in), optional :: method !! conversion method
        real(dp), intent(in), optional :: alpha !! gbt weight
        type(discrete_tf_result) :: out
        out = cont2discrete_tf_complex_multi(cmplx(b,0.0_dp,dp),cmplx(a,0.0_dp,dp),dt,method,alpha)
    end function cont2discrete_tf_real_multi

    function cont2discrete_tf_complex_multi(b,a,dt,method,alpha) result(out)
        complex(dp), intent(in) :: b(:,:) !! numerator rows for single-input/multiple-output system
        complex(dp), intent(in) :: a(:) !! common denominator, descending powers of s
        real(dp), intent(in) :: dt !! positive sampling interval
        character(*), intent(in), optional :: method !! conversion method
        real(dp), intent(in), optional :: alpha !! gbt weight
        type(discrete_tf_result) :: out
        type(state_space_result) :: ss
        type(discrete_state_space_result) :: sd
        type(state_space_tf_result) :: tf
        out%dt = dt
        ss = tf2ss(b,a)
        if (ss%status /= state_space_success) return
        sd = cont2discrete_ss_complex(ss%a,ss%b,ss%c,ss%d,dt,method,alpha)
        if (sd%status /= cont2discrete_success) then
            out%status = sd%status
            return
        end if
        tf = ss2tf(sd%a,sd%b,sd%c,sd%d)
        if (tf%status /= state_space_success) then
            out%status = cont2discrete_numerical_failure
            return
        end if
        out%b = tf%b
        out%a = tf%a
        out%status = cont2discrete_success
    end function cont2discrete_tf_complex_multi

    function cont2discrete_zpk_real(z,p,k,dt,method,alpha) result(out)
        real(dp), intent(in) :: z(:) !! analog zeros
        real(dp), intent(in) :: p(:) !! analog poles
        real(dp), intent(in) :: k !! analog gain
        real(dp), intent(in) :: dt !! positive sampling interval
        character(*), intent(in), optional :: method !! conversion method
        real(dp), intent(in), optional :: alpha !! gbt weight
        type(discrete_zpk_result) :: out
        out = cont2discrete_zpk_complex(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp), &
            cmplx(k,0.0_dp,dp),dt,method,alpha)
    end function cont2discrete_zpk_real

    function cont2discrete_zpk_real_complex_gain(z,p,k,dt,method,alpha) result(out)
        real(dp), intent(in) :: z(:) !! analog zeros
        real(dp), intent(in) :: p(:) !! analog poles
        complex(dp), intent(in) :: k !! complex analog gain
        real(dp), intent(in) :: dt !! positive sampling interval
        character(*), intent(in), optional :: method !! conversion method
        real(dp), intent(in), optional :: alpha !! gbt weight
        type(discrete_zpk_result) :: out
        out = cont2discrete_zpk_complex(cmplx(z,0.0_dp,dp),cmplx(p,0.0_dp,dp),k,dt,method,alpha)
    end function cont2discrete_zpk_real_complex_gain

    function cont2discrete_zpk_complex_real_gain(z,p,k,dt,method,alpha) result(out)
        complex(dp), intent(in) :: z(:) !! analog zeros
        complex(dp), intent(in) :: p(:) !! analog poles
        real(dp), intent(in) :: k !! analog gain
        real(dp), intent(in) :: dt !! positive sampling interval
        character(*), intent(in), optional :: method !! conversion method
        real(dp), intent(in), optional :: alpha !! gbt weight
        type(discrete_zpk_result) :: out
        out = cont2discrete_zpk_complex(z,p,cmplx(k,0.0_dp,dp),dt,method,alpha)
    end function cont2discrete_zpk_complex_real_gain

    function cont2discrete_zpk_complex(z,p,k,dt,method,alpha) result(out)
        complex(dp), intent(in) :: z(:) !! analog zeros
        complex(dp), intent(in) :: p(:) !! analog poles
        complex(dp), intent(in) :: k !! analog gain
        real(dp), intent(in) :: dt !! positive sampling interval
        character(*), intent(in), optional :: method !! conversion method
        real(dp), intent(in), optional :: alpha !! gbt weight
        type(discrete_zpk_result) :: out
        type(state_space_result) :: ss
        type(discrete_state_space_result) :: sd
        type(filter_zpk_result) :: zp
        complex(dp), allocatable :: nf(:), df(:)
        real(dp) :: weight
        character(:), allocatable :: selected
        integer :: nz, np
        out%dt = dt
        ss = zpk2ss(z,p,k)
        if (ss%status /= state_space_success) return
        sd = cont2discrete_ss_complex(ss%a,ss%b,ss%c,ss%d,dt,method,alpha)
        if (sd%status /= cont2discrete_success) then
            out%status = sd%status
            return
        end if

        ! Rational substitution gives the zeros, poles and gain directly
        ! for generalized bilinear methods. This avoids the ill-conditioned
        ! polynomial root extraction of repeated artificial zeros (e.g. the
        ! n-fold z=-1 zero of an all-pole bilinear filter).
        selected = 'zoh'
        if (present(method)) selected = trim(adjustl(lower_case(method)))
        weight = 0.0_dp
        select case(selected)
        case('gbt')
            weight = alpha
        case('bilinear','tustin')
            weight = 0.5_dp
        case('backward_diff')
            weight = 1.0_dp
        case('euler','forward_diff')
            weight = 0.0_dp
        end select
        if (selected == 'gbt' .or. selected == 'bilinear' .or. selected == 'tustin' .or. &
            selected == 'backward_diff' .or. selected == 'euler' .or. selected == 'forward_diff') then
            np = size(p)
            nz = size(z)
            nf = 1.0_dp-weight*dt*z
            df = 1.0_dp-weight*dt*p
            if (all(abs(nf)>1.0e-12_dp) .and. all(abs(df)>1.0e-12_dp)) then
                out%poles = (1.0_dp+(1.0_dp-weight)*dt*p)/df
                if (weight > 0.0_dp) then
                    allocate(out%zeros(np))
                    out%zeros(1:nz) = (1.0_dp+(1.0_dp-weight)*dt*z)/nf
                    out%zeros(nz+1:) = cmplx(-(1.0_dp-weight)/weight,0.0_dp,dp)
                    out%gain = k*(weight*dt)**(np-nz)*product(nf)/product(df)
                else
                    allocate(out%zeros(nz))
                    out%zeros = 1.0_dp+dt*z
                    out%gain = k*dt**(np-nz)
                end if
                if (all(ieee_is_finite(real(out%zeros,dp))) .and. &
                    all(ieee_is_finite(aimag(out%zeros))) .and. &
                    all(ieee_is_finite(real(out%poles,dp))) .and. &
                    all(ieee_is_finite(aimag(out%poles))) .and. &
                    ieee_is_finite(real(out%gain,dp)) .and. ieee_is_finite(aimag(out%gain))) then
                    out%status = cont2discrete_success
                    return
                end if
            end if
            if (allocated(out%zeros)) deallocate(out%zeros)
            if (allocated(out%poles)) deallocate(out%poles)
        end if
        zp = ss2zpk(sd%a,sd%b,sd%c,sd%d)
        if (zp%status /= zpk_conversion_success) then
            out%status = cont2discrete_numerical_failure
            return
        end if
        out%zeros = zp%zeros
        out%poles = zp%poles
        out%gain = zp%gain
        out%status = cont2discrete_success
    end function cont2discrete_zpk_complex

    subroutine solve_matrix(matrix,rhs,x,info)
        ! Dense complex Gaussian elimination with partial pivoting.
        complex(dp), intent(in) :: matrix(:,:) !! square coefficient matrix
        complex(dp), intent(in) :: rhs(:,:) !! right hand sides (n,nrhs)
        complex(dp), allocatable, intent(out) :: x(:,:) !! solved vectors
        integer, intent(out) :: info !! 0 on success, 1 on failure
        complex(dp), allocatable :: lu(:,:), temp(:)
        complex(dp) :: ratio
        integer :: n, nrhs, i, j, k, pivot
        info = 1
        n = size(matrix,1)
        nrhs = size(rhs,2)
        if (size(matrix,2) /= n .or. size(rhs,1) /= n) return
        if (.not. finite_matrix(matrix) .or. .not. finite_matrix(rhs)) return
        lu = matrix
        x = rhs
        if (n == 0) then
            info = 0
            return
        end if
        allocate(temp(max(n,nrhs)))
        do k = 1, n
            pivot = k-1+maxloc(abs(lu(k:n,k)),dim=1)
            if (abs(lu(pivot,k)) <= tiny(1.0_dp)) return
            if (pivot /= k) then
                temp(1:n) = lu(pivot,:)
                lu(pivot,:) = lu(k,:)
                lu(k,:) = temp(1:n)
                temp(1:nrhs) = x(pivot,:)
                x(pivot,:) = x(k,:)
                x(k,:) = temp(1:nrhs)
            end if
            do i = k+1, n
                ratio = lu(i,k)/lu(k,k)
                lu(i,k+1:n) = lu(i,k+1:n)-ratio*lu(k,k+1:n)
                x(i,:) = x(i,:)-ratio*x(k,:)
                lu(i,k) = (0.0_dp,0.0_dp)
            end do
        end do
        do i = n, 1, -1
            do j = i+1, n
                x(i,:) = x(i,:)-lu(i,j)*x(j,:)
            end do
            x(i,:) = x(i,:)/lu(i,i)
        end do
        if (.not. finite_matrix(x)) return
        info = 0
    end subroutine solve_matrix

    subroutine matrix_exp_pade(matrix,exp_matrix,info)
        ! Order-13 Pade approximation, scaling and squaring (Higham).
        complex(dp), intent(in) :: matrix(:,:) !! matrix to exponentiate
        complex(dp), allocatable, intent(out) :: exp_matrix(:,:) !! matrix exponential
        integer, intent(out) :: info !! 0 on success, 1 on failure
        real(dp), parameter :: theta13 = 5.371920351148152_dp
        real(dp), parameter :: coeff(0:13) = [ &
            64764752532480000.0_dp, 32382376266240000.0_dp, &
            7771770303897600.0_dp, 1187353796428800.0_dp, &
            129060195264000.0_dp, 10559470521600.0_dp, &
            670442572800.0_dp, 33522128640.0_dp, &
            1323241920.0_dp, 40840800.0_dp, &
            960960.0_dp, 16380.0_dp, 182.0_dp, 1.0_dp ]
        complex(dp), allocatable :: a(:,:), a2(:,:), a4(:,:), a6(:,:), identity(:,:), u(:,:), v(:,:)
        real(dp) :: norm_one
        integer :: n, j, squarings
        info = 1
        n = size(matrix,1)
        if (size(matrix,2) /= n) return
        if (.not. finite_matrix(matrix)) return
        if (n == 0) then
            allocate(exp_matrix(0,0))
            info = 0
            return
        end if
        norm_one = maxval(sum(abs(matrix),dim=1))
        if (.not. ieee_is_finite(norm_one)) return
        squarings = 0
        if (norm_one > theta13) then
            squarings = max(0,ceiling(log(norm_one/theta13)/log(2.0_dp)))
        end if
        a = matrix*scale(1.0_dp,-squarings)
        allocate(identity(n,n))
        identity = (0.0_dp,0.0_dp)
        do j = 1, n
            identity(j,j) = (1.0_dp,0.0_dp)
        end do
        a2 = matmul(a,a)
        a4 = matmul(a2,a2)
        a6 = matmul(a2,a4)
        u = matmul(a, &
            matmul(a6,coeff(13)*a6+coeff(11)*a4+coeff(9)*a2) + &
            coeff(7)*a6+coeff(5)*a4+coeff(3)*a2+coeff(1)*identity)
        v = matmul(a6,coeff(12)*a6+coeff(10)*a4+coeff(8)*a2) + &
            coeff(6)*a6+coeff(4)*a4+coeff(2)*a2+coeff(0)*identity
        call solve_matrix(v-u,v+u,exp_matrix,info)
        if (info /= 0) return
        do j = 1, squarings
            exp_matrix = matmul(exp_matrix,exp_matrix)
            if (.not. finite_matrix(exp_matrix)) then
                info = 1
                return
            end if
        end do
    end subroutine matrix_exp_pade

end module scifort_signal_cont2discrete
