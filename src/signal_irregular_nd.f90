! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Multichannel Lomb-Scargle for rank-1..15 real observations.  Each channel
! shares the same irregular observation times and frequency grid. Shared
! weights use cached trigonometric bases; channel-specific weights and
! optionally omitted nonfinite samples use the scalar reference kernel.
! Compute all channels into temporary storage before committing the output.
module scifort_signal_irregular_nd
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_irregular, only : lombscargle, lombscargle_amplitude, &
        irregular_status_success, irregular_status_invalid_input, irregular_status_singular
    implicit none
    private

    public :: lombscargle_into, lombscargle_amplitude_into
    public :: lombscargle_plan_into, lombscargle_plan_amplitude_into

    ! Precomputed weighted design for a *fixed* observation grid, frequency grid,
    ! and set of shared sample weights. Plan evaluation does not mutate the plan.
    type, public :: lombscargle_plan
        private
        integer :: sample_count = 0
        integer :: frequency_count = 0
        logical :: floating = .false.
        real(dp), allocatable :: wt(:)
        real(dp), allocatable :: cos_basis(:,:), sin_basis(:,:)
        real(dp), allocatable :: cos_mean(:), sin_mean(:)
        real(dp), allocatable :: cos_denom(:), sin_denom(:)
        real(dp), allocatable :: cos_tau(:), sin_tau(:)
    contains
        procedure :: initialize => lomb_plan_initialize
        procedure :: clear => lomb_plan_clear
        procedure :: is_initialized => lomb_plan_is_initialized
    end type lombscargle_plan

contains


    subroutine lomb_plan_clear(self)
        class(lombscargle_plan), intent(inout) :: self !! plan to invalidate and release

        if (allocated(self%wt)) deallocate(self%wt)
        if (allocated(self%cos_basis)) deallocate(self%cos_basis)
        if (allocated(self%sin_basis)) deallocate(self%sin_basis)
        if (allocated(self%cos_mean)) deallocate(self%cos_mean)
        if (allocated(self%sin_mean)) deallocate(self%sin_mean)
        if (allocated(self%cos_denom)) deallocate(self%cos_denom)
        if (allocated(self%sin_denom)) deallocate(self%sin_denom)
        if (allocated(self%cos_tau)) deallocate(self%cos_tau)
        if (allocated(self%sin_tau)) deallocate(self%sin_tau)
        self%sample_count = 0
        self%frequency_count = 0
        self%floating = .false.
    end subroutine lomb_plan_clear

    pure logical function lomb_plan_is_initialized(self) result(ready)
        class(lombscargle_plan), intent(in) :: self !! plan whose validity is queried

        ready = allocated(self%wt) .and. allocated(self%cos_basis) .and. &
            allocated(self%sin_basis) .and. self%sample_count > 0 .and. &
            self%frequency_count > 0
    end function lomb_plan_is_initialized

    subroutine lomb_plan_initialize(self, x, freqs, weights, floating_mean, status)
        class(lombscargle_plan), intent(inout) :: self !! existing plan; retained on failure
        real(dp), intent(in) :: x(:) !! shared observation times, all finite, nonempty
        real(dp), intent(in) :: freqs(:) !! angular-frequency grid, finite and nonempty
        real(dp), intent(in), optional :: weights(:) !! shared finite nonnegative weights; default uniform
        logical, intent(in), optional :: floating_mean !! fit an intercept at each frequency; default false
        integer, intent(out), optional :: status !! zero on success; previous plan retained on failure

        type(lombscargle_plan) :: candidate
        real(dp), allocatable :: cosine(:), sine(:)
        real(dp) :: wmax, wsum, cs, sn, cc, ss, cross, angle, tau
        integer :: n, nf, i, j

        if (present(status)) status = irregular_status_invalid_input
        n = size(x)
        nf = size(freqs)
        if (n < 1 .or. nf < 1) return
        if (any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(freqs))) return
        candidate%sample_count = n
        candidate%frequency_count = nf
        if (present(floating_mean)) candidate%floating = floating_mean
        allocate(candidate%wt(n))
        if (present(weights)) then
            if (size(weights) /= n) return
            if (any(.not. ieee_is_finite(weights))) return
            if (any(weights < 0.0_dp)) return
            wmax = maxval(weights)
            if (wmax <= 0.0_dp) return
            candidate%wt = weights/wmax
        else
            candidate%wt = 1.0_dp
        end if
        wsum = sum(candidate%wt)
        if (.not. ieee_is_finite(wsum) .or. wsum <= 0.0_dp) return
        candidate%wt = candidate%wt/wsum
        allocate(candidate%cos_basis(n,nf), candidate%sin_basis(n,nf))
        allocate(candidate%cos_mean(nf), candidate%sin_mean(nf))
        allocate(candidate%cos_denom(nf), candidate%sin_denom(nf))
        allocate(candidate%cos_tau(nf), candidate%sin_tau(nf))
        allocate(cosine(n), sine(n))
        do j = 1, nf
            do i = 1, n
                angle = freqs(j)*x(i)
                if (.not. ieee_is_finite(angle)) return
                cosine(i) = cos(angle)
                sine(i) = sin(angle)
            end do
            cs = sum(candidate%wt*cosine)
            sn = sum(candidate%wt*sine)
            cc = sum(candidate%wt*cosine**2)
            ss = 1.0_dp - cc
            cross = sum(candidate%wt*cosine*sine)
            if (candidate%floating) then
                cc = cc-cs*cs
                ss = ss-sn*sn
                cross = cross-cs*sn
            end if
            tau = 0.5_dp*atan2(2.0_dp*cross,cc-ss)
            candidate%cos_tau(j) = cos(tau)
            candidate%sin_tau(j) = sin(tau)
            candidate%cos_basis(:,j) = cosine*candidate%cos_tau(j) &
                + sine*candidate%sin_tau(j)
            candidate%sin_basis(:,j) = sine*candidate%cos_tau(j) &
                - cosine*candidate%sin_tau(j)
            candidate%cos_mean(j) = sum(candidate%wt*candidate%cos_basis(:,j))
            candidate%sin_mean(j) = sum(candidate%wt*candidate%sin_basis(:,j))
            candidate%cos_denom(j) = sum(candidate%wt*candidate%cos_basis(:,j)**2)
            candidate%sin_denom(j) = 1.0_dp - candidate%cos_denom(j)
            if (candidate%floating) then
                candidate%cos_denom(j) = candidate%cos_denom(j)-candidate%cos_mean(j)**2
                candidate%sin_denom(j) = candidate%sin_denom(j)-candidate%sin_mean(j)**2
            end if
            ! Store weighted basis columns so repeated fits need only dot products.
            candidate%cos_basis(:,j) = candidate%wt*candidate%cos_basis(:,j)
            candidate%sin_basis(:,j) = candidate%wt*candidate%sin_basis(:,j)
        end do
        ! Commit only after all validation and precomputation succeed. Move
        ! allocatables to avoid copying the potentially large frequency basis.
        call self%clear()
        self%sample_count = candidate%sample_count
        self%frequency_count = candidate%frequency_count
        self%floating = candidate%floating
        call move_alloc(candidate%wt,self%wt)
        call move_alloc(candidate%cos_basis,self%cos_basis)
        call move_alloc(candidate%sin_basis,self%sin_basis)
        call move_alloc(candidate%cos_mean,self%cos_mean)
        call move_alloc(candidate%sin_mean,self%sin_mean)
        call move_alloc(candidate%cos_denom,self%cos_denom)
        call move_alloc(candidate%sin_denom,self%sin_denom)
        call move_alloc(candidate%cos_tau,self%cos_tau)
        call move_alloc(candidate%sin_tau,self%sin_tau)
        if (present(status)) status = irregular_status_success
    end subroutine lomb_plan_initialize

    subroutine transform_planned_lanes(plan, source, dimensions, ax, destination, &
        amplitude_mode, precenter, normalize, status)
        type(lombscargle_plan), intent(in) :: plan !! fixed observation/frequency/weight design
        real(dp), intent(in) :: source(:) !! flattened finite signals in Fortran order
        integer, intent(in) :: dimensions(:) !! shapes of signal dimensions
        integer, intent(in) :: ax !! resolved one-based sample axis
        complex(dp), intent(out) :: destination(:) !! staged, flattened output
        logical, intent(in) :: amplitude_mode !! true to return fitted complex amplitudes
        logical, intent(in), optional :: precenter !! subtract unweighted mean per channel
        character(*), intent(in), optional :: normalize !! 'power' or 'normalize'
        integer, intent(out) :: status !! success or invalid/singular input status

        real(dp), parameter :: min_denom = 0.5_dp*epsilon(1.0_dp)
        real(dp), allocatable :: observations(:,:), means(:), variances(:)
        real(dp) :: yc, ys, a, b, improvement
        integer(int64) :: inner, outer, nc, nt, nf, group, ch, i, j, lane
        integer(int64) :: source_idx, output_idx
        logical :: center, normalized

        status = irregular_status_invalid_input
        if (.not. plan%is_initialized()) return
        if (any(.not. ieee_is_finite(source))) return
        normalized = .false.
        if (present(normalize)) then
            select case (trim(normalize))
            case ('power')
            case ('normalize')
                normalized = .true.
            case default
                return
            end select
        end if
        nt = int(plan%sample_count,int64)
        nf = int(plan%frequency_count,int64)
        inner = product(int(dimensions(:ax-1),int64))
        outer = product(int(dimensions(ax+1:),int64))
        nc = inner*outer
        if (nc < 1 .or. size(source,kind=int64) /= nt*nc) return
        if (size(destination,kind=int64) /= nf*nc) return
        center = .false.
        if (present(precenter)) center = precenter
        allocate(observations(nt,nc),means(nc),variances(nc))
        do group = 0_int64, outer-1_int64
            do ch = 0_int64, inner-1_int64
                lane = group*inner+ch+1_int64
                do i = 0_int64, nt-1_int64
                    source_idx = group*inner*nt+i*inner+ch+1_int64
                    observations(i+1_int64,lane) = source(source_idx)
                end do
            end do
        end do
        do lane = 1_int64, nc
            if (center) observations(:,lane) = observations(:,lane) &
                - sum(observations(:,lane))/real(nt,dp)
            means(lane) = sum(plan%wt*observations(:,lane))
            if (plan%floating) then
                variances(lane) = sum(plan%wt*(observations(:,lane)-means(lane))**2)
                if (maxval(abs(observations(:,lane)-observations(1,lane))) &
                    <= tiny(1.0_dp)) variances(lane) = 0.0_dp
            else
                variances(lane) = sum(plan%wt*observations(:,lane)**2)
            end if
            if (variances(lane) <= 0.0_dp) then
                status = irregular_status_singular
                return
            end if
        end do
        do j = 1_int64, nf
            do lane = 1_int64, nc
                yc = sum(observations(:,lane)*plan%cos_basis(:,j))
                ys = sum(observations(:,lane)*plan%sin_basis(:,j))
                if (plan%floating) then
                    yc = yc-means(lane)*plan%cos_mean(j)
                    ys = ys-means(lane)*plan%sin_mean(j)
                end if
                a = yc/max(plan%cos_denom(j),min_denom)
                b = ys/max(plan%sin_denom(j),min_denom)
                improvement = a*yc+b*ys
                group = (lane-1_int64)/inner
                ch = mod(lane-1_int64,inner)
                output_idx = group*inner*nf+(j-1_int64)*inner+ch+1_int64
                if (amplitude_mode) then
                    destination(output_idx) = cmplx(a,b,dp) &
                        * cmplx(plan%cos_tau(j),plan%sin_tau(j),dp)
                else
                    if (normalized) then
                        destination(output_idx) = cmplx(improvement/variances(lane),0.0_dp,dp)
                    else
                        destination(output_idx) = cmplx(0.5_dp*real(nt,dp)*improvement,0.0_dp,dp)
                    end if
                end if
            end do
        end do
        if (any(.not. ieee_is_finite(real(destination,dp))) .or. &
            any(.not. ieee_is_finite(aimag(destination)))) return
        status = irregular_status_success
    end subroutine transform_planned_lanes

    subroutine lombscargle_plan_into(plan, y, result, axis, precenter, normalize, status)
        type(lombscargle_plan), intent(in) :: plan !! initialized shared observation design
        real(dp), intent(in) :: y(..) !! finite real observations, rank 1 through 15
        real(dp), intent(inout) :: result(..) !! caller-owned power spectrum; atomic writes
        integer, intent(in), optional :: axis !! one-based/negative sample axis, default last
        logical, intent(in), optional :: precenter !! subtract unweighted mean, default false
        character(*), intent(in), optional :: normalize !! 'power' (default) or 'normalize'
        integer, intent(out), optional :: status !! zero on success, output otherwise unchanged

        real(dp), allocatable :: flat(:)
        complex(dp), allocatable :: staged(:)
        integer :: ax, st

        if (present(status)) status = irregular_status_invalid_input
        if (.not. plan%is_initialized()) return
        call check_layout(shape(y),shape(result),plan%sample_count, &
            plan%frequency_count,axis,ax,st)
        if (st /= irregular_status_success) return
        call flatten_real(y,flat)
        allocate(staged(size(result)))
        call transform_planned_lanes(plan,flat,shape(y),ax,staged,.false., &
            precenter,normalize,st)
        if (st == irregular_status_success) call scatter_real(result,staged)
        if (present(status)) status = st
    end subroutine lombscargle_plan_into

    subroutine lombscargle_plan_amplitude_into(plan, y, result, axis, precenter, status)
        type(lombscargle_plan), intent(in) :: plan !! initialized shared observation design
        real(dp), intent(in) :: y(..) !! finite real observations, rank 1 through 15
        complex(dp), intent(inout) :: result(..) !! caller-owned fitted amplitudes; atomic writes
        integer, intent(in), optional :: axis !! one-based/negative sample axis, default last
        logical, intent(in), optional :: precenter !! subtract unweighted mean, default false
        integer, intent(out), optional :: status !! zero on success, output otherwise unchanged

        real(dp), allocatable :: flat(:)
        complex(dp), allocatable :: staged(:)
        integer :: ax, st

        if (present(status)) status = irregular_status_invalid_input
        if (.not. plan%is_initialized()) return
        call check_layout(shape(y),shape(result),plan%sample_count, &
            plan%frequency_count,axis,ax,st)
        if (st /= irregular_status_success) return
        call flatten_real(y,flat)
        allocate(staged(size(result)))
        call transform_planned_lanes(plan,flat,shape(y),ax,staged,.true.,status=st, &
            precenter=precenter)
        if (st == irregular_status_success) call scatter_complex(result,staged)
        if (present(status)) status = st
    end subroutine lombscargle_plan_amplitude_into

    subroutine check_layout(input_shape, output_shape, sample_count, frequency_count, axis, ax, status)
        integer, intent(in) :: input_shape(:) !! extent of each observation dimension
        integer, intent(in) :: output_shape(:) !! extent of each output dimension
        integer, intent(in) :: sample_count !! length of the shared observation time vector
        integer, intent(in) :: frequency_count !! number of angular frequencies
        integer, intent(in), optional :: axis !! one-based or negative sample axis; default last
        integer, intent(out) :: ax !! resolved one-based sample axis
        integer, intent(out) :: status !! success when every extent is compatible
        integer :: rank_value, j

        status = irregular_status_invalid_input
        rank_value = size(input_shape)
        ax = rank_value
        if (rank_value < 1 .or. rank_value > 15 .or. size(output_shape) /= rank_value) return
        if (present(axis)) ax = axis
        if (ax < 0) ax = ax + rank_value + 1
        if (ax < 1 .or. ax > rank_value) return
        if (sample_count < 1 .or. frequency_count < 1) return
        if (any(input_shape <= 0) .or. any(output_shape <= 0)) return
        if (input_shape(ax) /= sample_count .or. output_shape(ax) /= frequency_count) return
        do j = 1, rank_value
            if (j == ax) cycle
            if (input_shape(j) /= output_shape(j)) return
        end do
        status = irregular_status_success
    end subroutine check_layout

    subroutine flatten_real(x, values)
        real(dp), intent(in) :: x(..) !! observations, rank 1 through 15
        real(dp), allocatable, intent(out) :: values(:) !! Fortran-order flattened copy

        select rank (x)
        rank (1)
            values = reshape(x, [size(x)])
        rank (2)
            values = reshape(x, [size(x)])
        rank (3)
            values = reshape(x, [size(x)])
        rank (4)
            values = reshape(x, [size(x)])
        rank (5)
            values = reshape(x, [size(x)])
        rank (6)
            values = reshape(x, [size(x)])
        rank (7)
            values = reshape(x, [size(x)])
        rank (8)
            values = reshape(x, [size(x)])
        rank (9)
            values = reshape(x, [size(x)])
        rank (10)
            values = reshape(x, [size(x)])
        rank (11)
            values = reshape(x, [size(x)])
        rank (12)
            values = reshape(x, [size(x)])
        rank (13)
            values = reshape(x, [size(x)])
        rank (14)
            values = reshape(x, [size(x)])
        rank (15)
            values = reshape(x, [size(x)])
        rank default
            allocate(values(0))
        end select
    end subroutine flatten_real

    subroutine scatter_real(y, values)
        real(dp), intent(inout) :: y(..) !! output committed on complete success
        complex(dp), intent(in) :: values(:) !! staged spectra in Fortran order
        select rank (y)
        rank (1)
            y = reshape(real(values, dp), shape(y))
        rank (2)
            y = reshape(real(values, dp), shape(y))
        rank (3)
            y = reshape(real(values, dp), shape(y))
        rank (4)
            y = reshape(real(values, dp), shape(y))
        rank (5)
            y = reshape(real(values, dp), shape(y))
        rank (6)
            y = reshape(real(values, dp), shape(y))
        rank (7)
            y = reshape(real(values, dp), shape(y))
        rank (8)
            y = reshape(real(values, dp), shape(y))
        rank (9)
            y = reshape(real(values, dp), shape(y))
        rank (10)
            y = reshape(real(values, dp), shape(y))
        rank (11)
            y = reshape(real(values, dp), shape(y))
        rank (12)
            y = reshape(real(values, dp), shape(y))
        rank (13)
            y = reshape(real(values, dp), shape(y))
        rank (14)
            y = reshape(real(values, dp), shape(y))
        rank (15)
            y = reshape(real(values, dp), shape(y))
        rank default
            error stop 'unreachable output rank in scatter_real'
        end select
    end subroutine scatter_real

    subroutine scatter_complex(y, values)
        complex(dp), intent(inout) :: y(..) !! output committed on complete success
        complex(dp), intent(in) :: values(:) !! staged amplitudes in Fortran order
        select rank (y)
        rank (1)
            y = reshape(values, shape(y))
        rank (2)
            y = reshape(values, shape(y))
        rank (3)
            y = reshape(values, shape(y))
        rank (4)
            y = reshape(values, shape(y))
        rank (5)
            y = reshape(values, shape(y))
        rank (6)
            y = reshape(values, shape(y))
        rank (7)
            y = reshape(values, shape(y))
        rank (8)
            y = reshape(values, shape(y))
        rank (9)
            y = reshape(values, shape(y))
        rank (10)
            y = reshape(values, shape(y))
        rank (11)
            y = reshape(values, shape(y))
        rank (12)
            y = reshape(values, shape(y))
        rank (13)
            y = reshape(values, shape(y))
        rank (14)
            y = reshape(values, shape(y))
        rank (15)
            y = reshape(values, shape(y))
        rank default
            error stop 'unreachable output rank in scatter_complex'
        end select
    end subroutine scatter_complex

    subroutine transform_shared_lanes(x, source, dimensions, ax, freqs, destination, &
        amplitude_mode, precenter, normalize, weights, floating_mean, status)
        real(dp), intent(in) :: x(:) !! common observation times for all channels
        real(dp), intent(in) :: source(:) !! flattened finite channel observations
        integer, intent(in) :: dimensions(:) !! source array extents
        integer, intent(in) :: ax !! resolved one-based observation axis
        real(dp), intent(in) :: freqs(:) !! shared angular-frequency grid
        complex(dp), intent(out) :: destination(:) !! staged output, Fortran order
        logical, intent(in) :: amplitude_mode !! emit complex sinusoid coefficients
        logical, intent(in), optional :: precenter !! remove unweighted channel means
        character(*), intent(in), optional :: normalize !! power or normalized power
        real(dp), intent(in), optional :: weights(:) !! common finite nonnegative weights
        logical, intent(in), optional :: floating_mean !! regress a separate intercept
        integer, intent(out) :: status !! success or first invalid/singular channel

        real(dp), parameter :: min_denom = 0.5_dp*epsilon(1.0_dp)
        real(dp), allocatable :: observations(:,:), wt(:), &
            means(:), variances(:), cosine(:), sine(:), rotated_cos(:), rotated_sin(:)
        real(dp) :: maxweight, tau, cos_tau, sin_tau, cs, sn, cc, cross, ss, cr, sr, ccr, ssr
        real(dp) :: yc, ys, a, b, improvement
        integer(int64) :: inner, outer, nc, nt, nf, i, j, lane, group, ch
        integer(int64) :: source_index, output_index
        logical :: floating, center, normalized

        status = irregular_status_invalid_input
        nt = int(size(x),int64)
        nf = int(size(freqs),int64)
        inner = product(int(dimensions(:ax-1),int64))
        outer = product(int(dimensions(ax+1:),int64))
        nc = inner*outer
        if (nt < 1 .or. nf < 1 .or. nc < 1) return
        if (size(source,kind=int64) /= nt*nc .or. &
            size(destination,kind=int64) /= nf*nc) return
        if (any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(freqs))) return
        if (any(.not. ieee_is_finite(source))) return
        normalized = .false.
        if (present(normalize)) then
            select case (trim(normalize))
            case ('power')
            case ('normalize')
                normalized = .true.
            case default
                return
            end select
        end if
        allocate(wt(nt))
        if (present(weights)) then
            if (size(weights,kind=int64) /= nt) return
            if (any(.not. ieee_is_finite(weights))) return
            if (any(weights < 0.0_dp)) return
            maxweight = maxval(weights)
            if (maxweight <= 0.0_dp) return
            wt = weights/maxweight
        else
            wt = 1.0_dp
        end if
        if (sum(wt) <= 0.0_dp) return
        wt = wt/sum(wt)
        floating = .false.
        if (present(floating_mean)) floating = floating_mean
        center = .false.
        if (present(precenter)) center = precenter

        allocate(observations(nt,nc), means(nc), &
            variances(nc), cosine(nt), sine(nt), rotated_cos(nt), rotated_sin(nt))
        do group = 0_int64,outer-1_int64
            do ch = 0_int64,inner-1_int64
                lane = group*inner+ch+1_int64
                do i = 0_int64,nt-1_int64
                    source_index = group*inner*nt+i*inner+ch+1_int64
                    observations(i+1_int64,lane) = source(source_index)
                end do
            end do
        end do
        do lane = 1_int64,nc
            if (center) observations(:,lane) = observations(:,lane) &
                -sum(observations(:,lane))/real(nt,dp)
            means(lane) = sum(wt*observations(:,lane))
            if (floating) then
                variances(lane) = sum(wt*(observations(:,lane)-means(lane))**2)
                if (maxval(abs(observations(:,lane)-observations(1,lane))) &
                    <= tiny(1.0_dp)) variances(lane) = 0.0_dp
            else
                variances(lane) = sum(wt*observations(:,lane)**2)
            end if
            if (variances(lane) <= 0.0_dp) then
                status = irregular_status_singular
                return
            end if
        end do

        do j = 1_int64,nf
            do i = 1_int64,nt
                cosine(i) = cos(freqs(j)*x(i))
                sine(i) = sin(freqs(j)*x(i))
            end do
            cs = sum(wt*cosine)
            sn = sum(wt*sine)
            cc = sum(wt*cosine**2)
            ss = 1.0_dp-cc
            cross = sum(wt*cosine*sine)
            if (floating) then
                cc = cc-cs*cs
                ss = ss-sn*sn
                cross = cross-cs*sn
            end if
            tau = 0.5_dp*atan2(2.0_dp*cross, cc-ss)
            cos_tau = cos(tau)
            sin_tau = sin(tau)
            rotated_cos = cosine*cos_tau+sine*sin_tau
            rotated_sin = sine*cos_tau-cosine*sin_tau
            cr = sum(wt*rotated_cos)
            sr = sum(wt*rotated_sin)
            ccr = sum(wt*rotated_cos**2)
            ssr = 1.0_dp-ccr
            if (floating) then
                ccr = ccr-cr*cr
                ssr = ssr-sr*sr
            end if
            do lane = 1_int64,nc
                yc = sum(wt*observations(:,lane)*rotated_cos)
                ys = sum(wt*observations(:,lane)*rotated_sin)
                if (floating) then
                    yc = yc-means(lane)*cr
                    ys = ys-means(lane)*sr
                end if
                a = yc/max(ccr,min_denom)
                b = ys/max(ssr,min_denom)
                improvement = a*yc+b*ys
                group = (lane-1_int64)/inner
                ch = mod(lane-1_int64,inner)
                output_index = group*inner*nf+(j-1_int64)*inner+ch+1_int64
                if (amplitude_mode) then
                    destination(output_index) = cmplx(a,b,dp)*cmplx(cos_tau,sin_tau,dp)
                else
                    if (normalized) then
                        destination(output_index) = cmplx(improvement/variances(lane),0.0_dp,dp)
                    else
                        destination(output_index) = cmplx(0.5_dp*real(nt,dp)*improvement,0.0_dp,dp)
                    end if
                end if
            end do
        end do
        status = irregular_status_success
    end subroutine transform_shared_lanes

    subroutine transform_lanes(x, source, dimensions, ax, freqs, destination, &
        amplitude_mode, precenter, normalize, weights, floating_mean, status, &
        weights_flat, nan_policy)
        real(dp), intent(in) :: x(:) !! shared irregular observation times
        real(dp), intent(in) :: source(:) !! flattened real samples in Fortran order
        integer, intent(in) :: dimensions(:) !! input extents in Fortran order
        integer, intent(in) :: ax !! selected positive sample axis
        real(dp), intent(in) :: freqs(:) !! angular frequencies, radians per time
        complex(dp), intent(out) :: destination(:) !! flat, staged complete output
        logical, intent(in) :: amplitude_mode !! true for complex fitted amplitudes
        logical, intent(in), optional :: precenter !! subtract unweighted sample mean
        character(*), intent(in), optional :: normalize !! 'power' or 'normalize' (power only)
        real(dp), intent(in), optional :: weights(:) !! shared nonnegative sample weights
        logical, intent(in), optional :: floating_mean !! fit a separate intercept per channel/bin
        integer, intent(out) :: status !! status of first failing channel
        real(dp), intent(in), optional :: weights_flat(:) !! per-channel weights, flattened
        character(*), intent(in), optional :: nan_policy !! 'raise' (default) or 'omit'

        integer(int64) :: inner, outer, nx, nf, batch, channel, j, source_idx, result_idx
        integer :: st
        real(dp), allocatable :: samples(:), powers(:), channel_weights(:)
        complex(dp), allocatable :: amplitudes(:)
        logical, allocatable :: keep(:)
        logical :: omit

        status = irregular_status_invalid_input
        omit = .false.
        if (present(nan_policy)) then
            select case (trim(nan_policy))
            case ('raise')
            case ('omit')
                omit = .true.
            case default
                return
            end select
        end if
        if (present(weights) .and. present(weights_flat)) return
        if (present(weights_flat)) then
            if (size(weights_flat) /= size(source)) return
            if (any(.not. ieee_is_finite(weights_flat))) return
            if (any(weights_flat < 0.0_dp)) return
        end if
        if (present(weights)) then
            if (size(weights) /= size(x)) return
            if (any(.not. ieee_is_finite(weights))) return
            if (any(weights < 0.0_dp)) return
        end if
        if (any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(freqs))) return
        if (.not. omit .and. .not. present(weights_flat)) then
            call transform_shared_lanes(x,source,dimensions,ax,freqs,destination, &
                amplitude_mode,precenter,normalize,weights,floating_mean,status)
            return
        end if
        inner = product(int(dimensions(:ax-1), int64))
        outer = product(int(dimensions(ax+1:), int64))
        nx = int(dimensions(ax), int64)
        nf = size(freqs, kind=int64)
        if (inner < 1 .or. outer < 1 .or. nx < 1 .or. nf < 1) return
        if (size(source, kind=int64) /= inner*nx*outer) return
        if (size(destination, kind=int64) /= inner*nf*outer) return

        allocate(samples(nx),keep(nx),channel_weights(nx))
        do batch = 0_int64, outer-1_int64
            do channel = 0_int64, inner-1_int64
                do j = 0_int64, nx-1_int64
                    source_idx = batch*inner*nx + j*inner + channel + 1_int64
                    samples(j+1_int64) = source(source_idx)
                    if (present(weights_flat)) channel_weights(j+1_int64) = weights_flat(source_idx)
                end do
                if (omit) then
                    keep = ieee_is_finite(samples)
                    if (present(weights_flat)) then
                        keep = keep .and. (channel_weights > 0.0_dp)
                    else if (present(weights)) then
                        keep = keep .and. (weights > 0.0_dp)
                    end if
                    if (.not. any(keep)) return
                    if (amplitude_mode) then
                        if (present(weights_flat)) then
                            amplitudes = lombscargle_amplitude(pack(x,keep),pack(samples,keep), &
                                freqs,precenter=precenter,weights=pack(channel_weights,keep), &
                                floating_mean=floating_mean,status=st)
                        else if (present(weights)) then
                            amplitudes = lombscargle_amplitude(pack(x,keep),pack(samples,keep), &
                                freqs,precenter=precenter,weights=pack(weights,keep), &
                                floating_mean=floating_mean,status=st)
                        else
                            amplitudes = lombscargle_amplitude(pack(x,keep),pack(samples,keep), &
                                freqs,precenter=precenter,floating_mean=floating_mean,status=st)
                        end if
                    else
                        if (present(weights_flat)) then
                            powers = lombscargle(pack(x,keep),pack(samples,keep),freqs, &
                                precenter=precenter,normalize=normalize, &
                                weights=pack(channel_weights,keep), &
                                floating_mean=floating_mean,status=st)
                        else if (present(weights)) then
                            powers = lombscargle(pack(x,keep),pack(samples,keep),freqs, &
                                precenter=precenter,normalize=normalize, &
                                weights=pack(weights,keep),floating_mean=floating_mean,status=st)
                        else
                            powers = lombscargle(pack(x,keep),pack(samples,keep),freqs, &
                                precenter=precenter,normalize=normalize, &
                                floating_mean=floating_mean,status=st)
                        end if
                    end if
                else if (amplitude_mode) then
                    if (present(weights_flat)) then
                        amplitudes = lombscargle_amplitude(x,samples,freqs, &
                            precenter=precenter,weights=channel_weights, &
                            floating_mean=floating_mean,status=st)
                    else
                        amplitudes = lombscargle_amplitude(x,samples,freqs, &
                            precenter=precenter,weights=weights, &
                            floating_mean=floating_mean,status=st)
                    end if
                else
                    if (present(weights_flat)) then
                        powers = lombscargle(x,samples,freqs,precenter=precenter, &
                            normalize=normalize,weights=channel_weights, &
                            floating_mean=floating_mean,status=st)
                    else
                        powers = lombscargle(x,samples,freqs,precenter=precenter, &
                            normalize=normalize,weights=weights, &
                            floating_mean=floating_mean,status=st)
                    end if
                end if
                if (st /= irregular_status_success) then
                    status = st
                    return
                end if
                do j = 0_int64, nf-1_int64
                    result_idx = batch*inner*nf + j*inner + channel + 1_int64
                    if (amplitude_mode) then
                        destination(result_idx) = amplitudes(j+1_int64)
                    else
                        destination(result_idx) = cmplx(powers(j+1_int64), 0.0_dp, dp)
                    end if
                end do
            end do
        end do
        status = irregular_status_success
    end subroutine transform_lanes

    subroutine lombscargle_into(x, y, freqs, result, axis, precenter, normalize, &
        weights, floating_mean, status, weights_nd, nan_policy)
        real(dp), intent(in) :: x(:) !! observation times, common to every channel
        real(dp), intent(in) :: y(..) !! real signals, rank 1..15
        real(dp), intent(in) :: freqs(:) !! angular frequencies in radians per time unit
        real(dp), intent(inout) :: result(..) !! same rank as y, with sample axis replaced by freqs
        integer, intent(in), optional :: axis !! one-based or negative sample axis; default last
        logical, intent(in), optional :: precenter !! subtract unweighted channel means
        character(*), intent(in), optional :: normalize !! 'power' (default) or 'normalize'
        real(dp), intent(in), optional :: weights(:) !! nonnegative weights, length size(x)
        logical, intent(in), optional :: floating_mean !! fit a per-frequency offset
        integer, intent(out), optional :: status !! zero on success; output unchanged on failure
        real(dp), intent(in), optional :: weights_nd(..) !! nonnegative per-channel weights; same shape as y
        character(*), intent(in), optional :: nan_policy !! 'raise' (default) or 'omit' nonfinite y and zero weights

        real(dp), allocatable :: flat(:), flat_weights(:)
        complex(dp), allocatable :: staged(:)
        integer :: ax, st

        if (present(status)) status = irregular_status_invalid_input
        call check_layout(shape(y), shape(result), size(x), size(freqs), axis, ax, st)
        if (st /= irregular_status_success) return
        if (present(weights_nd)) then
            if (size(shape(weights_nd)) /= size(shape(y))) return
            if (any(shape(weights_nd) /= shape(y))) return
            call flatten_real(weights_nd,flat_weights)
        end if
        call flatten_real(y, flat)
        allocate(staged(size(result)))
        call transform_lanes(x, flat, shape(y), ax, freqs, staged, .false., &
            precenter, normalize, weights, floating_mean, st, &
            weights_flat=flat_weights, nan_policy=nan_policy)
        if (st == irregular_status_success) call scatter_real(result, staged)
        if (present(status)) status = st
    end subroutine lombscargle_into

    subroutine lombscargle_amplitude_into(x, y, freqs, result, axis, precenter, &
        weights, floating_mean, status, weights_nd, nan_policy)
        real(dp), intent(in) :: x(:) !! observation times, common to every channel
        real(dp), intent(in) :: y(..) !! real signals, rank 1..15
        real(dp), intent(in) :: freqs(:) !! angular frequencies in radians per time unit
        complex(dp), intent(inout) :: result(..) !! same rank as y, with sample axis replaced by freqs
        integer, intent(in), optional :: axis !! one-based or negative sample axis; default last
        logical, intent(in), optional :: precenter !! subtract unweighted channel means
        real(dp), intent(in), optional :: weights(:) !! nonnegative weights, length size(x)
        logical, intent(in), optional :: floating_mean !! fit a per-frequency offset
        integer, intent(out), optional :: status !! zero on success; output unchanged on failure
        real(dp), intent(in), optional :: weights_nd(..) !! nonnegative per-channel weights; same shape as y
        character(*), intent(in), optional :: nan_policy !! 'raise' (default) or 'omit' nonfinite y and zero weights

        real(dp), allocatable :: flat(:), flat_weights(:)
        complex(dp), allocatable :: staged(:)
        integer :: ax, st

        if (present(status)) status = irregular_status_invalid_input
        call check_layout(shape(y), shape(result), size(x), size(freqs), axis, ax, st)
        if (st /= irregular_status_success) return
        if (present(weights_nd)) then
            if (size(shape(weights_nd)) /= size(shape(y))) return
            if (any(shape(weights_nd) /= shape(y))) return
            call flatten_real(weights_nd,flat_weights)
        end if
        call flatten_real(y, flat)
        allocate(staged(size(result)))
        call transform_lanes(x, flat, shape(y), ax, freqs, staged, .true., &
            precenter=precenter, weights=weights, floating_mean=floating_mean, status=st, &
            weights_flat=flat_weights, nan_policy=nan_policy)
        if (st == irregular_status_success) call scatter_complex(result, staged)
        if (present(status)) status = st
    end subroutine lombscargle_amplitude_into

end module scifort_signal_irregular_nd
