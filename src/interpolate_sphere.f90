! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2022, Federico Perini. All rights reserved.
! Copyright (c) 2026 SciFort contributors
!
! Spherical coefficient constraints adapted from fpsphe/fprpsp in perazz/fitpack,
! src/fitpack_core.F90, commit 1907c7bad9f5c4557970e906bdd4909eb52b2d4c.
! Modernized upstream by Federico Perini from Paul Dierckx's FITPACK.
! Full terms: LICENSES/perazz-fitpack-BSD-3-Clause.txt.
! SciFort uses its own B-spline evaluators and dense minimum-norm SVD solve;
! no adaptive knot selection or FITPACK banded solver is imported.
module scifort_interpolate_sphere
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, &
        ieee_quiet_nan, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_constants, only : scifort_pi
    use scifort_interpolate_bspline, only : bspline_design_matrix
    use scifort_interpolate_fitpack_bivariate, only : fitpack_bivariate_tck_type, bisplev
    use scifort_linalg_svd, only : svd_lstsq
    implicit none
    private
    public :: lsq_sphere_bivariate_spline, sphere_bivariate_spline_type
    integer, parameter, public :: sphere_spline_success = 0
    integer, parameter, public :: sphere_spline_invalid_input = 1
    integer, parameter, public :: sphere_spline_rank_deficient = 2
    integer, parameter, public :: sphere_spline_numerical_failure = 3
    integer, parameter, public :: sphere_spline_size_limit = 4
    real(dp), parameter :: two_pi = 2.0_dp*scifort_pi
    integer, parameter :: max_dimension = 256
    integer(int64), parameter :: max_design_entries = 2000000_int64

    type :: sphere_bivariate_spline_type
        private
        type(fitpack_bivariate_tck_type) :: tck
        integer :: numerical_rank = 0
        integer :: dimension = 0
        real(dp) :: weighted_residual = 0.0_dp
        real(dp) :: value_scale = 1.0_dp
    contains
        procedure :: valid => sphere_valid
        procedure :: value => sphere_value
        procedure :: values => sphere_values
        procedure :: grid => sphere_grid
        procedure :: get_knots => sphere_knots
        procedure :: coefficients => sphere_coefficients
        procedure :: rank => sphere_rank
        procedure :: degrees_of_freedom => sphere_dimension
        procedure :: residual => sphere_residual
    end type sphere_bivariate_spline_type

contains

    function lsq_sphere_bivariate_spline(theta, phi, r, tt, tp, w, eps, status) result(fit)
        implicit none
        real(dp), intent(in) :: theta(:) !! Colatitudes in [0,pi], radians; at least two samples.
        real(dp), intent(in) :: phi(:) !! Longitudes in [0,2*pi], radians; paired with theta.
        real(dp), intent(in) :: r(:) !! Finite scalar observations, in arbitrary sample order.
        real(dp), intent(in) :: tt(:) !! Strictly increasing interior colatitude knots in (0,pi).
        real(dp), intent(in) :: tp(:) !! At least three increasing longitude knots in (0,2*pi).
        real(dp), intent(in), optional :: w(:) !! Positive residual weights; default all one.
        real(dp), intent(in), optional :: eps !! Squared SVD cutoff in (0,1); default 1e-16.
        integer, intent(out), optional :: status !! 0 success; 2 valid rank warning; 1/3/4 failure.
        type(sphere_bivariate_spline_type) :: fit
        real(dp), allocatable :: tx(:), ty(:), nodes(:), phase(:,:), trig(:,:), cs(:,:), &
            bt(:,:), bp(:,:), periodic(:,:), weights(:), design(:,:), rhs(:,:), solution(:,:), &
            raw(:,:), predicted(:)
        real(dp) :: cutoff, rscale, wscale, co, si, q, logfp, maxcoef
        integer :: m, nx, ny, np, nf, code, rank_phase, rank_fit, i, j, k, base
        integer(int64) :: nf64

        code = sphere_spline_invalid_input
        m = size(theta)
        if (m < 2 .or. size(phi) /= m .or. size(r) /= m) goto 900
        if (size(tp) < 3) goto 900
        if (.not. all(ieee_is_finite(theta)) .or. .not. all(ieee_is_finite(phi)) &
            .or. .not. all(ieee_is_finite(r))) goto 900
        if (.not. all(ieee_is_finite(tt)) .or. .not. all(ieee_is_finite(tp))) goto 900
        if (any(theta < 0.0_dp) .or. any(theta > scifort_pi)) goto 900
        if (any(phi < 0.0_dp) .or. any(phi > two_pi)) goto 900
        if (any(tt <= 0.0_dp) .or. any(tt >= scifort_pi)) goto 900
        if (any(tp <= 0.0_dp) .or. any(tp >= two_pi)) goto 900
        if (any(tt(2:) <= tt(:size(tt)-1)) .or. any(tp(2:) <= tp(:size(tp)-1))) goto 900
        cutoff = 1.0e-16_dp
        if (present(eps)) cutoff = eps
        if (.not. ieee_is_finite(cutoff)) goto 900
        if (cutoff <= 0.0_dp .or. cutoff >= 1.0_dp) goto 900
        if (present(w)) then
            if (size(w) /= m) goto 900
            if (.not. all(ieee_is_finite(w))) goto 900
            if (any(w <= 0.0_dp)) goto 900
        end if
        code = sphere_spline_size_limit
        if (int(size(tp),int64)+1_int64 > int(max_dimension,int64)) goto 900
        nf64 = 6_int64+(int(size(tp),int64)+1_int64)*int(size(tt),int64)
        if (nf64 > int(max_dimension,int64)) goto 900
        nf = int(nf64)
        np = size(tp)+1
        if (int(m,int64)*int(nf,int64) > max_design_entries) goto 900
        ! Longitude collocation and both coordinate design matrices also need
        ! bounded storage, even when there are no interior latitude knots.
        if (np > max_dimension) goto 900
        nx = size(tt)+4
        ny = np+3
        if (int(m,int64)*int(nf+nx+ny+np,int64) > max_design_entries) goto 900
        allocate(tx(nx+4),ty(ny+4),nodes(np),weights(m))
        tx(:4) = 0.0_dp
        tx(5:nx) = tt
        tx(nx+1:) = scifort_pi
        nodes(1) = 0.0_dp
        nodes(2:) = tp
        ty(4:np+3) = nodes
        ty(np+4) = two_pi
        ty(:3) = tp(size(tp)-2:)-two_pi
        ty(np+5:) = tp(:3)+two_pi

        code = sphere_spline_numerical_failure
        phase = bspline_design_matrix(nodes,ty,3,status=i)
        if (i /= 0) goto 900
        allocate(periodic(np,np),trig(np,2),cs(np,2))
        periodic = 0.0_dp
        do j = 1, ny
            k = modulo(j-1,np)+1
            periodic(:,k) = periodic(:,k)+phase(:,j)
        end do
        trig(:,1) = cos(nodes)
        trig(:,2) = sin(nodes)
        call svd_lstsq(periodic,trig,cs,i,rank=rank_phase, &
            rcond=real(np,dp)*epsilon(1.0_dp))
        if (i /= 0 .or. rank_phase /= np) goto 900
        deallocate(periodic)
        bt = bspline_design_matrix(theta,tx,3,status=i)
        if (i /= 0) goto 900
        bp = bspline_design_matrix(phi,ty,3,status=i)
        if (i /= 0) goto 900
        allocate(periodic(m,np),design(m,nf),rhs(m,1),solution(nf,1),raw(nx,ny))
        periodic = 0.0_dp
        do j = 1, ny
            k = modulo(j-1,np)+1
            periodic(:,k) = periodic(:,k)+bp(:,j)
        end do
        design = 0.0_dp
        do i = 1, m
            co = dot_product(periodic(i,:),cs(:,1))
            si = dot_product(periodic(i,:),cs(:,2))
            design(i,1:3) = [bt(i,1)+bt(i,2),bt(i,2)*co,bt(i,2)*si]
            do j = 3, nx-2
                base = 3+(j-3)*np
                design(i,base+1:base+np) = bt(i,j)*periodic(i,:)
            end do
            design(i,nf-2:nf) = [bt(i,nx-1)*co,bt(i,nx-1)*si,bt(i,nx-1)+bt(i,nx)]
        end do
        weights = 1.0_dp
        if (present(w)) weights = w
        wscale = maxval(weights)
        weights = weights/wscale
        if (any(weights == 0.0_dp)) goto 900
        rscale = maxval(abs(r))
        if (rscale == 0.0_dp) rscale = 1.0_dp
        rhs(:,1) = (r/rscale)*weights
        do i = 1, m
            design(i,:) = design(i,:)*weights(i)
        end do
        call svd_lstsq(design,rhs,solution,i,rank=rank_fit, &
            rcond=max(sqrt(cutoff),real(max(m,nf),dp)*epsilon(1.0_dp)))
        if (i /= 0) goto 900
        ! Reconstruct the standard latitude-major, longitude-fastest FITPACK
        ! tensor coefficients from the reduced spherical coordinates (fprpsp).
        raw(1,:) = solution(1,1)
        raw(nx,:) = solution(nf,1)
        do j = 1, ny
            k = modulo(j-1,np)+1
            raw(2,j) = solution(1,1)+solution(2,1)*cs(k,1)+solution(3,1)*cs(k,2)
            raw(nx-1,j) = solution(nf,1)+solution(nf-2,1)*cs(k,1)+solution(nf-1,1)*cs(k,2)
            do i = 3, nx-2
                base = 3+(i-3)*np
                raw(i,j) = solution(base+k,1)
            end do
        end do
        maxcoef = maxval(abs(raw))
        if (maxcoef > 1.0_dp) then
            if (rscale > huge(1.0_dp)/maxcoef) goto 900
        end if
        predicted = matmul(design,solution(:,1))-rhs(:,1)
        q = dot_product(predicted,predicted)
        fit%weighted_residual = 0.0_dp
        if (q > 0.0_dp) then
            logfp = log(q)+2.0_dp*log(rscale)+2.0_dp*log(wscale)
            if (logfp > log(huge(1.0_dp))) then
                fit%weighted_residual = ieee_value(0.0_dp,ieee_positive_inf)
            else if (logfp >= log(tiny(1.0_dp))+log(epsilon(1.0_dp))) then
                fit%weighted_residual = exp(logfp)
            end if
        end if
        fit%tck%tx = tx
        fit%tck%ty = ty
        ! transpose before reshape makes longitude vary fastest.
        fit%tck%c = reshape(transpose(raw),[nx*ny])
        fit%value_scale = rscale
        fit%tck%kx = 3
        fit%tck%ky = 3
        fit%numerical_rank = rank_fit
        fit%dimension = nf
        code = sphere_spline_success
        if (rank_fit < nf) code = sphere_spline_rank_deficient
900     if (present(status)) status = code
    end function lsq_sphere_bivariate_spline

    logical function sphere_valid(self) result(ok)
        implicit none
        class(sphere_bivariate_spline_type), intent(in) :: self !! Spherical spline object.
        ok = self%dimension > 0 .and. self%tck%valid()
    end function sphere_valid

    function sphere_value(self, theta, phi, dtheta, dphi, status) result(value)
        implicit none
        class(sphere_bivariate_spline_type), intent(in) :: self !! Valid spherical spline.
        real(dp), intent(in) :: theta !! Colatitude in [0,pi], radians.
        real(dp), intent(in) :: phi !! Finite longitude, periodically wrapped into [0,2*pi).
        integer, intent(in), optional :: dtheta !! Colatitude derivative order, 0..2; default 0.
        integer, intent(in), optional :: dphi !! Longitude derivative order, 0..2; default 0.
        integer, intent(out), optional :: status !! 0 success, 1 invalid, 3 unrepresentable value.
        real(dp) :: value, longitude
        integer :: dt, dp_order, code
        value = ieee_value(0.0_dp,ieee_quiet_nan)
        code = sphere_spline_invalid_input
        if (.not. self%valid()) goto 900
        if (.not. ieee_is_finite(theta) .or. .not. ieee_is_finite(phi)) goto 900
        if (theta < 0.0_dp .or. theta > scifort_pi) goto 900
        dt = 0
        dp_order = 0
        if (present(dtheta)) dt = dtheta
        if (present(dphi)) dp_order = dphi
        if (dt < 0 .or. dt > 2 .or. dp_order < 0 .or. dp_order > 2) goto 900
        longitude = modulo(phi,two_pi)
        value = bisplev(theta,longitude,self%tck,dx=dt,dy=dp_order,status=code)
        if (code /= 0) goto 900
        code = sphere_spline_numerical_failure
        if (.not. ieee_is_finite(value)) then
            value = ieee_value(0.0_dp,ieee_quiet_nan)
            goto 900
        end if
        if (abs(value) > 1.0_dp) then
            if (self%value_scale > huge(1.0_dp)/abs(value)) then
                value = ieee_value(0.0_dp,ieee_quiet_nan)
                goto 900
            end if
        end if
        value = value*self%value_scale
        code = sphere_spline_success
900     if (present(status)) status = code
    end function sphere_value

    function sphere_values(self, theta, phi, dtheta, dphi, status) result(values)
        implicit none
        class(sphere_bivariate_spline_type), intent(in) :: self !! Valid spherical spline.
        real(dp), intent(in) :: theta(:) !! Paired colatitudes in [0,pi], radians.
        real(dp), intent(in) :: phi(:) !! Same-length finite longitudes, radians.
        integer, intent(in), optional :: dtheta !! Colatitude derivative order, 0..2.
        integer, intent(in), optional :: dphi !! Longitude derivative order, 0..2.
        integer, intent(out), optional :: status !! 0 success; 1 invalid query, 3 numerical failure.
        real(dp), allocatable :: values(:)
        integer :: i, code, overall
        allocate(values(size(theta)))
        values = ieee_value(0.0_dp,ieee_quiet_nan)
        overall = sphere_spline_invalid_input
        if (size(phi) /= size(theta) .or. .not. self%valid()) goto 900
        if (present(dtheta)) then
            if (dtheta < 0 .or. dtheta > 2) goto 900
        end if
        if (present(dphi)) then
            if (dphi < 0 .or. dphi > 2) goto 900
        end if
        overall = sphere_spline_success
        do i = 1, size(theta)
            values(i) = self%value(theta(i),phi(i),dtheta,dphi,code)
            if (code /= sphere_spline_success) overall = code
        end do
900     if (present(status)) status = overall
    end function sphere_values

    function sphere_grid(self, theta, phi, dtheta, dphi, status) result(values)
        implicit none
        class(sphere_bivariate_spline_type), intent(in) :: self !! Valid spherical spline.
        real(dp), intent(in) :: theta(:) !! Grid colatitudes in [0,pi], radians; arbitrary order.
        real(dp), intent(in) :: phi(:) !! Grid longitudes, radians; arbitrary order.
        integer, intent(in), optional :: dtheta !! Colatitude derivative order, 0..2.
        integer, intent(in), optional :: dphi !! Longitude derivative order, 0..2.
        integer, intent(out), optional :: status !! 0 success; 1 invalid query, 3 numerical failure.
        real(dp), allocatable :: values(:,:)
        integer :: i, j, code, overall
        allocate(values(size(theta),size(phi)))
        values = ieee_value(0.0_dp,ieee_quiet_nan)
        overall = sphere_spline_invalid_input
        if (.not. self%valid()) goto 900
        if (present(dtheta)) then
            if (dtheta < 0 .or. dtheta > 2) goto 900
        end if
        if (present(dphi)) then
            if (dphi < 0 .or. dphi > 2) goto 900
        end if
        overall = sphere_spline_success
        do j = 1, size(phi)
            do i = 1, size(theta)
                values(i,j) = self%value(theta(i),phi(j),dtheta,dphi,code)
                if (code /= sphere_spline_success) overall = code
            end do
        end do
900     if (present(status)) status = overall
    end function sphere_grid

    subroutine sphere_knots(self, tt, tp)
        implicit none
        class(sphere_bivariate_spline_type), intent(in) :: self !! Spherical spline object.
        real(dp), allocatable, intent(out) :: tt(:) !! Full latitude knots; empty if invalid.
        real(dp), allocatable, intent(out) :: tp(:) !! Periodic longitude knots; empty if invalid.
        allocate(tt(0),tp(0))
        if (.not. self%valid()) return
        tt = self%tck%tx
        tp = self%tck%ty
    end subroutine sphere_knots

    function sphere_coefficients(self) result(c)
        implicit none
        class(sphere_bivariate_spline_type), intent(in) :: self !! Spherical spline object.
        real(dp), allocatable :: c(:)
        allocate(c(0))
        if (self%valid()) c = self%tck%c*self%value_scale
    end function sphere_coefficients

    integer function sphere_rank(self) result(rank)
        implicit none
        class(sphere_bivariate_spline_type), intent(in) :: self !! Spherical spline object.
        rank = self%numerical_rank
    end function sphere_rank

    integer function sphere_dimension(self) result(n)
        implicit none
        class(sphere_bivariate_spline_type), intent(in) :: self !! Spherical spline object.
        n = self%dimension
    end function sphere_dimension

    real(dp) function sphere_residual(self) result(fp)
        implicit none
        class(sphere_bivariate_spline_type), intent(in) :: self !! Spherical spline object.
        fp = ieee_value(0.0_dp,ieee_quiet_nan)
        if (self%valid()) fp = self%weighted_residual
    end function sphere_residual
end module scifort_interpolate_sphere
