! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Smooth cubic interpolation of 3D orientations on SO(3).
! The nonlinear angular-rate continuity equations follow the published
! RotationSpline mathematical formulation (see docs/spatial_rotation_spline.md).
! Uses a self-contained block-tridiagonal 3x3 solver; no BLAS/LAPACK.
module scifort_spatial_rotation_spline
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_spatial_transform, only : rotation_type, rotation_from_rotvec, rotation_from_quat, &
        spatial_rotation_success, spatial_rotation_invalid_input
    implicit none
    private

    public :: rotation_spline_type, rotation_spline, spatial_rotation_no_convergence

    integer, parameter :: spatial_rotation_no_convergence = 2

    integer, parameter :: max_newton_iterations = 10
    real(dp), parameter :: iteration_tolerance = 1.0e-9_dp

    type :: rotation_spline_type
        private
        real(dp), allocatable :: key_times(:)
        type(rotation_type), allocatable :: key_rotations(:)
        ! Rotation-vector cubics, coefficients of s**3, s**2, s, 1.
        real(dp), allocatable :: coefficients(:,:,:)
        logical :: defined = .false.
    contains
        procedure :: valid => spline_valid
        procedure :: at => spline_at
        procedure :: angular_rate => spline_angular_rate
        procedure :: angular_acceleration => spline_angular_acceleration
        procedure :: evaluate_path => spline_evaluate_path
    end type rotation_spline_type

contains

    pure subroutine set_status(status, code)
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        integer, intent(in) :: code !! status value to return
        if (present(status)) status = code
    end subroutine set_status

    pure logical function spline_valid(self) result(ok)
        class(rotation_spline_type), intent(in) :: self !! spline to inspect
        ok = self%defined
    end function spline_valid

    pure function skew(v) result(m)
        real(dp), intent(in) :: v(3) !! vector to cross with a second vector
        real(dp) :: m(3,3)
        m(1,:) = [0.0_dp, -v(3), v(2)]
        m(2,:) = [v(3), 0.0_dp, -v(1)]
        m(3,:) = [-v(2), v(1), 0.0_dp]
    end function skew

    pure function cross(a,b) result(c)
        real(dp), intent(in) :: a(3) !! first input vector
        real(dp), intent(in) :: b(3) !! second input vector
        real(dp) :: c(3)
        c = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), &
             a(1)*b(2)-a(2)*b(1)]
    end function cross

    pure function identity3() result(m)
        real(dp) :: m(3,3)
        integer :: i
        m = 0.0_dp
        do i = 1,3
            m(i,i) = 1.0_dp
        end do
    end function identity3

    pure function angular_to_rotvec_derivative(rotvec) result(m)
        real(dp), intent(in) :: rotvec(3) !! relative rotation vector, radians
        real(dp) :: m(3,3), s(3,3), theta, theta2, k
        theta2 = sum(rotvec**2)
        theta = sqrt(theta2)
        if (theta > 1.0e-4_dp) then
            k = (1.0_dp - 0.5_dp*theta/tan(0.5_dp*theta)) / theta2
        else
            k = 1.0_dp/12.0_dp + theta2/720.0_dp + theta2**2/30240.0_dp
        end if
        s = skew(rotvec)
        m = identity3() + 0.5_dp*s + k*matmul(s,s)
    end function angular_to_rotvec_derivative

    pure function rotvec_derivative_to_angular(rotvec) result(m)
        real(dp), intent(in) :: rotvec(3) !! relative rotation vector, radians
        real(dp) :: m(3,3), s(3,3), theta, theta2, k1, k2
        theta2 = sum(rotvec**2)
        theta = sqrt(theta2)
        if (theta > 1.0e-4_dp) then
            k1 = (1.0_dp-cos(theta))/theta2
            k2 = (theta-sin(theta))/(theta2*theta)
        else
            k1 = 0.5_dp - theta2/24.0_dp + theta2**2/720.0_dp
            k2 = 1.0_dp/6.0_dp - theta2/120.0_dp + theta2**2/5040.0_dp
        end if
        s = skew(rotvec)
        m = identity3() - k1*s + k2*matmul(s,s)
    end function rotvec_derivative_to_angular

    pure function nonlinear_angular_acceleration(rotvec, rotvec_dot) result(term)
        real(dp), intent(in) :: rotvec(3) !! relative rotation vector, radians
        real(dp), intent(in) :: rotvec_dot(3) !! time derivative of relative rotation vector
        real(dp) :: term(3), theta, theta2, k1, k2, k3, cp(3), ccp(3), dccp(3)
        theta2 = sum(rotvec**2)
        theta = sqrt(theta2)
        if (theta > 1.0e-4_dp) then
            k1 = (-theta*sin(theta)-2.0_dp*(cos(theta)-1.0_dp)) / theta**4
            k2 = (-2.0_dp*theta+3.0_dp*sin(theta)-theta*cos(theta)) / theta**5
            k3 = (theta-sin(theta))/theta**3
        else
            k1 = 1.0_dp/12.0_dp - theta2/180.0_dp + theta2**2/6720.0_dp
            k2 = -1.0_dp/60.0_dp + theta2/1260.0_dp
            k3 = 1.0_dp/6.0_dp - theta2/120.0_dp + theta2**2/5040.0_dp
        end if
        cp = cross(rotvec,rotvec_dot)
        ccp = cross(rotvec,cp)
        dccp = cross(rotvec_dot,cp)
        term = dot_product(rotvec,rotvec_dot)*(k1*cp + k2*ccp) + k3*dccp
    end function nonlinear_angular_acceleration

    ! Solve a 3x3 system with four right-hand sides, pivoting rows.
    pure subroutine solve3(a, rhs, result, ok)
        real(dp), intent(in) :: a(3,3) !! nonsingular 3x3 matrix
        real(dp), intent(in) :: rhs(3,4) !! right sides, 3 matrix columns and one vector
        real(dp), intent(out) :: result(3,4) !! solutions to the four systems
        logical, intent(out) :: ok !! true when all pivots are usable
        real(dp) :: m(3,3), b(3,4), tmp(4), row(3), fac, tol, scale
        integer :: i, j, pivot
        result = 0.0_dp
        ok = .false.
        m = a
        b = rhs
        if (.not. all(ieee_is_finite(m)) .or. .not. all(ieee_is_finite(b))) return
        scale = maxval(abs(m))
        if (scale <= 0.0_dp) return
        tol = 16.0_dp*epsilon(scale)*scale
        do i = 1,3
            pivot = i-1+maxloc(abs(m(i:3,i)),dim=1)
            if (abs(m(pivot,i)) <= tol) return
            if (pivot /= i) then
                row = m(pivot,:)
                m(pivot,:) = m(i,:)
                m(i,:) = row
                tmp = b(pivot,:)
                b(pivot,:) = b(i,:)
                b(i,:) = tmp
            end if
            do j = i+1,3
                fac = m(j,i)/m(i,i)
                m(j,i:3) = m(j,i:3)-fac*m(i,i:3)
                b(j,:) = b(j,:)-fac*b(i,:)
            end do
        end do
        do i = 3,1,-1
            result(i,:) = (b(i,:)-matmul(m(i,i+1:3),result(i+1:3,:))) / m(i,i)
        end do
        ok = all(ieee_is_finite(result))
    end subroutine solve3

    pure subroutine solve_block_tridiagonal(lower, diagonal, upper, rhs, solution, ok)
        real(dp), intent(in) :: lower(:,:,:) !! (3,3,n), first unused
        real(dp), intent(in) :: diagonal(:,:,:) !! (3,3,n) block diagonal
        real(dp), intent(in) :: upper(:,:,:) !! (3,3,n), last unused
        real(dp), intent(in) :: rhs(:,:) !! (3,n) right-hand side
        real(dp), intent(out) :: solution(:,:) !! (3,n) block solution
        logical, intent(out) :: ok !! success flag
        real(dp) :: ubar(3,3,size(rhs,2)), ybar(3,size(rhs,2))
        real(dp) :: a(3,3), input(3,4), output(3,4)
        integer :: i, n
        n = size(rhs,2)
        solution = 0.0_dp
        ok = .false.
        do i = 1,n
            a = diagonal(:,:,i)
            input(:,1:3) = upper(:,:,i)
            input(:,4) = rhs(:,i)
            if (i > 1) then
                a = a-matmul(lower(:,:,i),ubar(:,:,i-1))
                input(:,4) = input(:,4)-matmul(lower(:,:,i),ybar(:,i-1))
            end if
            call solve3(a,input,output,ok)
            if (.not. ok) return
            ubar(:,:,i) = output(:,1:3)
            ybar(:,i) = output(:,4)
        end do
        solution(:,n) = ybar(:,n)
        do i = n-1,1,-1
            solution(:,i) = ybar(:,i)-matmul(ubar(:,:,i),solution(:,i+1))
        end do
        ok = all(ieee_is_finite(solution))
    end subroutine solve_block_tridiagonal

    function rotation_spline(times, rotations, status) result(spline)
        real(dp), intent(in) :: times(:) !! finite strictly increasing keyframe times, at least two
        type(rotation_type), intent(in) :: rotations(:) !! corresponding orientations
        integer, intent(out), optional :: status !! zero on success, one for invalid input
        type(rotation_spline_type) :: spline
        integer :: n, m, num_interior, i, iteration
        type(rotation_type) :: rel
        real(dp), allocatable :: dt(:), rv(:,:), rates(:,:), rdot_end(:,:), starts(:,:)
        real(dp), allocatable :: a(:,:,:), a_inv(:,:,:), lower(:,:,:), diagonal(:,:,:), upper(:,:,:)
        real(dp), allocatable :: rhs0(:,:), rhs(:,:), solution(:,:)
        real(dp) :: first(3)
        logical :: ok, converged

        call set_status(status,spatial_rotation_invalid_input)
        n = size(times)
        if (n < 2 .or. size(rotations) /= n) return
        if (.not. all(ieee_is_finite(times))) return
        if (any(times(2:n) <= times(1:n-1))) return
        do i = 1,n
            if (.not. rotations(i)%valid()) return
        end do
        m = n-1
        dt = times(2:n)-times(1:n-1)
        if (.not. all(ieee_is_finite(dt))) return
        if (any(dt <= 0.0_dp)) return
        allocate(rv(3,m),rates(3,m),rdot_end(3,m),starts(3,m))
        do i = 1,m
            rel = rotations(i)%inverse()
            rel = rel%compose(rotations(i+1))
            rv(:,i) = rel%as_rotvec()
            rates(:,i) = rv(:,i)/dt(i)
        end do
        if (.not. all(ieee_is_finite(rates))) return
        first = rates(:,1)
        if (m == 1) then
            rdot_end(:,1) = rates(:,1)
        else
            allocate(a(3,3,m),a_inv(3,3,m))
            do i = 1,m
                a(:,:,i) = angular_to_rotvec_derivative(rv(:,i))
                a_inv(:,:,i) = rotvec_derivative_to_angular(rv(:,i))
            end do
            num_interior = m-1
            allocate(lower(3,3,num_interior),diagonal(3,3,num_interior), &
                     upper(3,3,num_interior),rhs0(3,num_interior), &
                     rhs(3,num_interior),solution(3,num_interior))
            lower = 0.0_dp
            upper = 0.0_dp
            do i = 1,num_interior
                diagonal(:,:,i) = 4.0_dp*(1.0_dp/dt(i)+1.0_dp/dt(i+1))*identity3()
                rhs0(:,i) = 6.0_dp*(rv(:,i)/dt(i)**2 + rv(:,i+1)/dt(i+1)**2)
                if (i > 1) lower(:,:,i) = (2.0_dp/dt(i))*a_inv(:,:,i)
                if (i < num_interior) upper(:,:,i) = (2.0_dp/dt(i+1))*a(:,:,i+1)
            end do
            rhs0(:,1) = rhs0(:,1) - (2.0_dp/dt(1))*matmul(a_inv(:,:,1),first)
            rhs0(:,num_interior) = rhs0(:,num_interior) - &
                (2.0_dp/dt(m))*matmul(a(:,:,m),rates(:,m))
            if (.not. all(ieee_is_finite(rhs0))) return
            do iteration = 1,max_newton_iterations
                do i = 1,num_interior
                    rhs(:,i) = rhs0(:,i)-nonlinear_angular_acceleration( &
                        rv(:,i),matmul(a(:,:,i),rates(:,i)))
                end do
                call solve_block_tridiagonal(lower,diagonal,upper,rhs,solution,ok)
                if (.not. ok) return
                converged = all(abs(solution-rates(:,1:num_interior)) < &
                                iteration_tolerance*(1.0_dp+abs(solution)))
                rates(:,1:num_interior) = solution
                if (converged) exit
            end do
            ! As in the SciPy algorithm, use the last iterate when the fixed
            ! iteration budget is reached. Below, reject results whose huge
            ! coefficients destroy the required keyframe interpolation.
            do i = 1,m
                rdot_end(:,i) = matmul(a(:,:,i),rates(:,i))
            end do
        end if
        starts(:,1) = first
        starts(:,2:m) = rates(:,1:m-1)
        if (.not. all(ieee_is_finite(rdot_end))) return
        allocate(spline%coefficients(3,4,m))
        do i = 1,m
            spline%coefficients(:,1,i) = (-2.0_dp*rv(:,i)+dt(i)*(starts(:,i)+rdot_end(:,i)))/dt(i)**3
            spline%coefficients(:,2,i) = (3.0_dp*rv(:,i)-dt(i)*(2.0_dp*starts(:,i)+rdot_end(:,i)))/dt(i)**2
            spline%coefficients(:,3,i) = starts(:,i)
            spline%coefficients(:,4,i) = 0.0_dp
        end do
        if (.not. all(ieee_is_finite(spline%coefficients))) then
            deallocate(spline%coefficients)
            call set_status(status,spatial_rotation_no_convergence)
            return
        end if
        ! Extremely poorly conditioned near-pi sequences can diverge in
        ! the nonlinear iteration. Reject if coefficient cancellation no
        ! longer reconstructs the prescribed relative rotation vectors.
        do i = 1,m
            if (any(abs(((spline%coefficients(:,1,i)*dt(i) + &
                spline%coefficients(:,2,i))*dt(i) + &
                spline%coefficients(:,3,i))*dt(i)-rv(:,i)) > 5.0e-8_dp)) then
                deallocate(spline%coefficients)
                call set_status(status,spatial_rotation_no_convergence)
                return
            end if
        end do
        spline%key_times = times
        spline%key_rotations = rotations
        spline%defined = .true.
        call set_status(status,spatial_rotation_success)
    end function rotation_spline

    pure integer function segment_index(self, time) result(segment)
        class(rotation_spline_type), intent(in) :: self !! initialized rotation spline
        real(dp), intent(in) :: time !! query timestamp
        integer :: lo, hi, mid, m
        m = size(self%key_times)-1
        lo = 1
        hi = m+1
        ! Right-continuous interval choice, including endpoint extrapolation.
        do while (hi-lo > 1)
            mid = (lo+hi)/2
            if (time >= self%key_times(mid)) then
                lo = mid
            else
                hi = mid
            end if
        end do
        segment = min(m,lo)
    end function segment_index

    pure subroutine evaluate_rotvec(self,time,interval,rv,rdot,rddot)
        class(rotation_spline_type), intent(in) :: self !! initialized spline
        real(dp), intent(in) :: time !! finite evaluation time
        integer, intent(out) :: interval !! selected segment, 1-based
        real(dp), intent(out) :: rv(3) !! relative rotation vector
        real(dp), intent(out) :: rdot(3) !! first time derivative
        real(dp), intent(out) :: rddot(3) !! second time derivative
        real(dp) :: s
        interval = segment_index(self,time)
        s = time-self%key_times(interval)
        rv = ((self%coefficients(:,1,interval)*s + self%coefficients(:,2,interval))*s + &
            self%coefficients(:,3,interval))*s
        rdot = (3.0_dp*self%coefficients(:,1,interval)*s + &
                2.0_dp*self%coefficients(:,2,interval))*s + self%coefficients(:,3,interval)
        rddot = 6.0_dp*self%coefficients(:,1,interval)*s + 2.0_dp*self%coefficients(:,2,interval)
    end subroutine evaluate_rotvec

    function spline_at(self,time,status) result(r)
        class(rotation_spline_type), intent(in) :: self !! constructed spline
        real(dp), intent(in) :: time !! sample timestamp; extrapolation allowed
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        type(rotation_type) :: r, rel
        real(dp) :: rv(3),rdot(3),rddot(3)
        integer :: i
        r = rotation_from_quat([0.0_dp,0.0_dp,0.0_dp,0.0_dp])
        call set_status(status,spatial_rotation_invalid_input)
        if (.not. self%valid() .or. .not. ieee_is_finite(time)) return
        call evaluate_rotvec(self,time,i,rv,rdot,rddot)
        rel = rotation_from_rotvec(rv)
        if (.not. rel%valid()) return
        r = self%key_rotations(i)%compose(rel)
        if (r%valid()) call set_status(status,spatial_rotation_success)
    end function spline_at

    function spline_angular_rate(self,time,status) result(rate)
        class(rotation_spline_type), intent(in) :: self !! constructed spline
        real(dp), intent(in) :: time !! evaluation timestamp, extrapolation allowed
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        real(dp) :: rate(3),rv(3),rdot(3),rddot(3)
        integer :: i
        rate = ieee_value(0.0_dp,ieee_quiet_nan)
        call set_status(status,spatial_rotation_invalid_input)
        if (.not. self%valid() .or. .not. ieee_is_finite(time)) return
        call evaluate_rotvec(self,time,i,rv,rdot,rddot)
        rate = matmul(rotvec_derivative_to_angular(rv),rdot)
        if (all(ieee_is_finite(rate))) call set_status(status,spatial_rotation_success)
    end function spline_angular_rate

    function spline_angular_acceleration(self,time,status) result(acc)
        class(rotation_spline_type), intent(in) :: self !! constructed spline
        real(dp), intent(in) :: time !! evaluation timestamp, extrapolation allowed
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        real(dp) :: acc(3),rv(3),rdot(3),rddot(3)
        integer :: i
        acc = ieee_value(0.0_dp,ieee_quiet_nan)
        call set_status(status,spatial_rotation_invalid_input)
        if (.not. self%valid() .or. .not. ieee_is_finite(time)) return
        call evaluate_rotvec(self,time,i,rv,rdot,rddot)
        acc = matmul(rotvec_derivative_to_angular(rv),rddot) + nonlinear_angular_acceleration(rv,rdot)
        if (all(ieee_is_finite(acc))) call set_status(status,spatial_rotation_success)
    end function spline_angular_acceleration

    subroutine spline_evaluate_path(self, sample_times, values, angular_rates, angular_accelerations, status)
        class(rotation_spline_type), intent(in) :: self !! constructed spline
        real(dp), intent(in) :: sample_times(:) !! evaluation timestamps; unordered, extrapolation allowed
        type(rotation_type), allocatable, intent(out) :: values(:) !! interpolated orientations
        real(dp), allocatable, intent(out), optional :: angular_rates(:,:) !! (n,3), rad/time
        real(dp), allocatable, intent(out), optional :: angular_accelerations(:,:) !! (n,3), rad/time**2
        integer, intent(out), optional :: status !! zero on success, one on invalid input
        integer :: i,n,s
        call set_status(status,spatial_rotation_invalid_input)
        allocate(values(0))
        if (present(angular_rates)) allocate(angular_rates(0,3))
        if (present(angular_accelerations)) allocate(angular_accelerations(0,3))
        if (.not. self%valid() .or. .not. all(ieee_is_finite(sample_times))) return
        n = size(sample_times)
        deallocate(values)
        allocate(values(n))
        if (present(angular_rates)) then
            deallocate(angular_rates)
            allocate(angular_rates(n,3))
        end if
        if (present(angular_accelerations)) then
            deallocate(angular_accelerations)
            allocate(angular_accelerations(n,3))
        end if
        do i = 1,n
            values(i) = self%at(sample_times(i),status=s)
            if (s /= spatial_rotation_success) exit
            if (present(angular_rates)) then
                angular_rates(i,:) = self%angular_rate(sample_times(i),status=s)
                if (s /= spatial_rotation_success) exit
            end if
            if (present(angular_accelerations)) then
                angular_accelerations(i,:) = self%angular_acceleration(sample_times(i),status=s)
                if (s /= spatial_rotation_success) exit
            end if
        end do
        if (i <= n) then
            deallocate(values)
            allocate(values(0))
            if (present(angular_rates)) then
                deallocate(angular_rates)
                allocate(angular_rates(0,3))
            end if
            if (present(angular_accelerations)) then
                deallocate(angular_accelerations)
                allocate(angular_accelerations(0,3))
            end if
            return
        end if
        call set_status(status,spatial_rotation_success)
    end subroutine spline_evaluate_path
end module scifort_spatial_rotation_spline
