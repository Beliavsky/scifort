! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2022, Federico Perini. All rights reserved.
! Copyright (c) 2026 SciFort contributors
!
! Spherical coefficient constraints adapted from fpsphe/fprpsp in perazz/fitpack,
! src/fitpack_core.F90, commit 1907c7bad9f5c4557970e906bdd4909eb52b2d4c.
! Modernized upstream by Federico Perini from Paul Dierckx's FITPACK.
! Full terms: LICENSES/perazz-fitpack-BSD-3-Clause.txt.
! SciFort uses its own B-spline evaluators and dense minimum-norm SVD solve;
! the separate grid core adapts FITPACK's banded grid interpolation/smoothing.
module scifort_interpolate_sphere
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_value, &
        ieee_quiet_nan, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_constants, only : scifort_pi
    use scifort_interpolate_bspline, only : bspline_design_matrix
    use scifort_interpolate_fitpack_bivariate, only : fitpack_bivariate_tck_type, bisplev
    use scifort_linalg_svd, only : svd_lstsq
    use scifort_interpolate_sphere_grid_core, only : fitpack_sphere_grid
    use scifort_interpolate_sphere_scattered_core, only : fitpack_sphere_scattered
    implicit none
    private
    public :: lsq_sphere_bivariate_spline, sphere_bivariate_spline_type
    public :: rect_sphere_bivariate_spline,smooth_sphere_bivariate_spline
    integer, parameter, public :: sphere_spline_success = 0
    integer, parameter, public :: sphere_spline_invalid_input = 1
    integer, parameter, public :: sphere_spline_rank_deficient = 2
    integer, parameter, public :: sphere_spline_numerical_failure = 3
    integer, parameter, public :: sphere_spline_size_limit = 4
    integer, parameter, public :: sphere_spline_unsupported_smoothing = 5
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
        real(dp) :: longitude_origin = 0.0_dp
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

    function smooth_sphere_bivariate_spline(theta,phi,r,w,s,eps,ntest,npest,status) result(fit)
        implicit none
        real(dp), intent(in) :: theta(:) !! Finite scattered colatitudes in [0,pi], radians.
        real(dp), intent(in) :: phi(:) !! Paired finite longitudes in [0,2*pi], radians.
        real(dp), intent(in) :: r(:) !! Finite observations; at least two, in arbitrary order.
        real(dp), intent(in), optional :: w(:) !! Positive residual weights; default all one.
        real(dp), intent(in), optional :: s !! Nonnegative weighted residual target; default zero.
        real(dp), intent(in), optional :: eps !! Squared relative rank threshold; default 1e-16.
        integer, intent(in), optional :: ntest !! Latitude knot capacity, >=8; default 8+sqrt(m/2).
        integer, intent(in), optional :: npest !! Longitude knot capacity, >=8; same default.
        integer, intent(out), optional :: status !! 0 success, 2 rank warning; 1/3/4 failure.
        type(sphere_bivariate_spline_type) :: fit
        real(dp), allocatable :: weights(:),tt(:),tp(:),coef(:)
        real(dp) :: target,threshold,rscale,wscale,scaled_target,logfp,fp,maxcoef
        integer :: m,ntcap,npcap,code,ier,rank
        integer(int64) :: ntt,npp,ncc,ib1,ib3,work1,work2,entries
        code = sphere_spline_invalid_input
        m = size(theta)
        if (m < 2 .or. size(phi) /= m .or. size(r) /= m) goto 900
        if (.not. all(ieee_is_finite(theta)) .or. .not. all(ieee_is_finite(phi)) &
            .or. .not. all(ieee_is_finite(r))) goto 900
        if (any(theta < 0.0_dp) .or. any(theta > scifort_pi)) goto 900
        if (any(phi < 0.0_dp) .or. any(phi > two_pi)) goto 900
        target = 0.0_dp
        if (present(s)) target = s
        if (.not. ieee_is_finite(target)) goto 900
        if (target < 0.0_dp) goto 900
        threshold = 1e-16_dp
        if (present(eps)) threshold = eps
        if (.not. ieee_is_finite(threshold)) goto 900
        if (threshold <= 0.0_dp .or. threshold >= 1.0_dp) goto 900
        if (present(w)) then
            if (size(w) /= m) goto 900
            if (.not. all(ieee_is_finite(w))) goto 900
            if (any(w <= 0.0_dp)) goto 900
        end if
        ntcap = 8+int(sqrt(real(m/2,dp)))
        npcap = ntcap
        if (present(ntest)) ntcap = ntest
        if (present(npest)) npcap = npest
        if (ntcap < 8 .or. npcap < 8) goto 900
        code = sphere_spline_size_limit
        if (int(ntcap,int64)*int(npcap,int64) > max_design_entries) goto 900
        ntt = int(ntcap,int64)-7
        npp = int(npcap,int64)-7
        ncc = 6+npp*(ntt-1)
        ib1 = max(4*npp,6+3*npp)
        ib3 = max(ib1+3,6+3*npp)
        work1 = 185+52*npp+10*ntt+14*ntt*npp+8*(int(m,int64)+(ntt-1)*npp**2)
        work2 = ncc*ib3+ncc+ib3
        entries = work1+work2+int(ntcap-4,int64)*int(npcap-4,int64) &
            +5_int64*int(m,int64)+ntt*npp+ntcap+npcap
        if (entries > max_design_entries) goto 900
        code = sphere_spline_numerical_failure
        allocate(weights(m))
        weights = 1.0_dp
        wscale = 1.0_dp
        if (present(w)) then
            wscale = maxval(w)
            weights = w/wscale
        end if
        if (any(weights <= 0.0_dp)) goto 900
        rscale = maxval(abs(r))
        if (rscale == 0.0_dp) rscale = 1.0_dp
        scaled_target = 0.0_dp
        if (target > 0.0_dp) then
            logfp = log(target)-2.0_dp*(log(rscale)+log(wscale))
            if (logfp < log(tiny(1.0_dp))) goto 900
            scaled_target = exp(min(logfp,log(huge(1.0_dp))-1.0_dp))
        end if
        call fitpack_sphere_scattered(theta,phi,r/rscale,weights,scaled_target,threshold, &
            ntcap,npcap,int(work1),int(work2),tt,tp,coef,fp,rank,ier)
        if (ier == 1 .or. ier > 10) then
            code = sphere_spline_size_limit
            goto 900
        end if
        if (ier > 0) goto 900
        maxcoef = maxval(abs(coef))
        if (maxcoef > 1.0_dp) then
            if (rscale > huge(1.0_dp)/maxcoef) goto 900
        end if
        fit%tck%tx = tt
        fit%tck%ty = tp
        fit%tck%c = coef
        fit%tck%kx = 3
        fit%tck%ky = 3
        fit%value_scale = rscale
        fit%weighted_residual = 0.0_dp
        if (fp > 0.0_dp) then
            logfp = log(fp)+2.0_dp*(log(rscale)+log(wscale))
            if (logfp > log(huge(1.0_dp))) then
                fit%weighted_residual = ieee_value(0.0_dp,ieee_positive_inf)
            else
                fit%weighted_residual = exp(logfp)
            end if
        end if
        fit%dimension = 6+(size(tt)-8)*(size(tp)-7)
        if (size(tt) == 8 .and. size(tp) == 8) fit%dimension = 2
        fit%numerical_rank = rank
        code = sphere_spline_success
        if (rank >= 0 .and. rank < fit%dimension) code = sphere_spline_rank_deficient
900     if (present(status)) status = code
    end function smooth_sphere_bivariate_spline

    function rect_sphere_bivariate_spline(u,v,r,s,pole_continuity,pole_values, &
        pole_exact,pole_flat,status) result(fit)
        implicit none
        real(dp), intent(in) :: u(:) !! Strictly increasing colatitudes in (0,pi); at least four.
        real(dp), intent(in) :: v(:) !! Increasing longitudes; v(1) in [-pi,pi), span below 2*pi.
        real(dp), intent(in) :: r(:,:) !! Finite grid samples, shape (size(u),size(v)).
        real(dp), intent(in), optional :: s !! Nonnegative squared-residual target; default zero.
        logical, intent(in), optional :: pole_continuity(:) !! Two C1 pole flags; default false.
        real(dp), intent(in), optional :: pole_values(:) !! Two pole values; NaN means unspecified.
        logical, intent(in), optional :: pole_exact(:) !! Two exactness flags, default false.
        logical, intent(in), optional :: pole_flat(:) !! Two zero-slope flags; require continuity.
        integer, intent(out), optional :: status !! 0 success; 1/3/4 input/numerical/size failure.
        type(sphere_bivariate_spline_type) :: fit
        logical :: continuity(2),given(2),exact(2),flat(2)
        real(dp) :: poles(2),rscale,maxcoef,separation,target,scaled_target,fp,logfp
        real(dp), allocatable :: tx(:),ty(:),coef(:)
        integer :: mu,mv,i,code
        integer(int64) :: nu64,nv64,work64

        code = sphere_spline_invalid_input
        mu = size(u)
        mv = size(v)
        if (mu < 4 .or. mv < 4) goto 900
        if (size(r,1) /= mu .or. size(r,2) /= mv) goto 900
        if (.not. all(ieee_is_finite(u)) .or. .not. all(ieee_is_finite(v))) goto 900
        if (.not. all(ieee_is_finite(r))) goto 900
        if (u(1) <= 0.0_dp .or. u(mu) >= scifort_pi) goto 900
        if (any(u(2:) <= u(:mu-1)) .or. any(v(2:) <= v(:mv-1))) goto 900
        if (v(1) < -scifort_pi .or. v(1) >= scifort_pi) goto 900
        if (v(mv) >= v(1)+two_pi) goto 900
        target = 0.0_dp
        if (present(s)) then
            if (.not. ieee_is_finite(s)) goto 900
            if (s < 0.0_dp) goto 900
            target = s
        end if
        continuity = .false.
        exact = .false.
        flat = .false.
        poles = 0.0_dp
        given = .false.
        if (present(pole_continuity)) then
            if (size(pole_continuity) /= 2) goto 900
            continuity = pole_continuity
        end if
        if (present(pole_exact)) then
            if (size(pole_exact) /= 2) goto 900
            exact = pole_exact
        end if
        if (present(pole_flat)) then
            if (size(pole_flat) /= 2) goto 900
            flat = pole_flat
        end if
        if (any(flat .and. .not. continuity)) goto 900
        if (present(pole_values)) then
            if (size(pole_values) /= 2) goto 900
            given = .not. ieee_is_nan(pole_values)
            do i = 1, 2
                if (.not. given(i)) cycle
                if (.not. ieee_is_finite(pole_values(i))) goto 900
                poles(i) = pole_values(i)
            end do
        end if
        code = sphere_spline_size_limit
        nu64 = int(mu,int64)+8_int64
        nv64 = int(mv,int64)+7_int64
        work64 = 12+nu64*(int(mv,int64)+nv64+3)+24*nv64+4*int(mu,int64) &
            +8*int(mv,int64)+max(nu64,int(mv,int64)+nv64)
        if (work64+nu64*nv64+int(mu,int64)*int(mv,int64) > max_design_entries) goto 900
        ! Keep arithmetic away from nearly coincident grid sites/endpoints.
        code = sphere_spline_numerical_failure
        separation = 128.0_dp*epsilon(1.0_dp)*two_pi
        if (min(u(1),scifort_pi-u(mu),minval(u(2:)-u(:mu-1))) <= separation) goto 900
        if (min(two_pi-(v(mv)-v(1)),minval(v(2:)-v(:mv-1))) <= separation) goto 900
        rscale = maxval(abs(r))
        do i = 1, 2
            if (given(i)) rscale = max(rscale,abs(poles(i)))
        end do
        if (rscale == 0.0_dp) rscale = 1.0_dp
        scaled_target = 0.0_dp
        if (target > 0.0_dp) then
            logfp = log(target)-2.0_dp*log(rscale)
            if (logfp < log(tiny(1.0_dp))) goto 900
            scaled_target = exp(min(logfp,log(huge(1.0_dp))-1.0_dp))
        end if
        call fitpack_sphere_grid(u,v,r/rscale,continuity,poles/rscale,given,exact,flat, &
            scaled_target,tx,ty,coef,fp,i)
        if (i /= 0) goto 900
        maxcoef = maxval(abs(coef))
        if (maxcoef > 1.0_dp) then
            if (rscale > huge(1.0_dp)/maxcoef) goto 900
        end if
        fit%tck%tx = tx
        fit%tck%ty = ty
        fit%tck%c = coef
        fit%tck%kx = 3
        fit%tck%ky = 3
        fit%value_scale = rscale
        fit%longitude_origin = v(1)
        fit%weighted_residual = 0.0_dp
        if (fp > 0.0_dp) then
            logfp = log(fp)+2.0_dp*log(rscale)
            if (logfp > log(huge(1.0_dp))) then
                fit%weighted_residual = ieee_value(0.0_dp,ieee_positive_inf)
            else
                fit%weighted_residual = exp(logfp)
            end if
        end if
        fit%dimension = 2+2*count(continuity) &
            +(size(tx)-6-count(continuity))*(size(ty)-7)
        if (size(ty) == 8) fit%dimension = size(tx)-4-count(continuity)
        fit%numerical_rank = -1 ! The banded FITPACK path does not estimate numerical rank.
        code = sphere_spline_success
900     if (present(status)) status = code
    end function rect_sphere_bivariate_spline

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
        longitude = self%longitude_origin+modulo(phi-self%longitude_origin,two_pi)
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
