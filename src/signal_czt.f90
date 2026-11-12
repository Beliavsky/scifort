! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Chirp Z transforms and zoom FFTs. Unit-circle spirals use Bluestein's
! FFT convolution; nonunit spirals use a Horner evaluation to avoid extreme
! intermediate chirp factors. No third-party FFT dependency is required.
module scifort_signal_czt
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : forward_fft
    implicit none
    private

    integer, parameter, public :: czt_status_success = 0
    integer, parameter, public :: czt_status_invalid_input = 1
    integer, parameter, public :: czt_status_nonfinite_result = 2
    real(dp), parameter :: pi = acos(-1.0_dp)

    interface czt
        module procedure czt_complex, czt_real
    end interface czt
    interface zoom_fft
        module procedure zoom_fft_complex, zoom_fft_real
        module procedure zoom_fft_complex_scalar, zoom_fft_real_scalar
    end interface zoom_fft
    public :: czt, zoom_fft, czt_points

    ! Reusable, immutable-after-initialization plans for repeated transforms.
    ! For off-unit spirals, the plan caches evaluation points and uses Horner
    ! evaluation to avoid the potentially overflowing Bluestein chirps.
    type, public :: czt_plan
        private
        integer :: n = 0
        integer :: m = 0
        integer :: fft_len = 0
        logical :: ready = .false.
        logical :: use_fft = .false.
        complex(dp) :: a = (1.0_dp, 0.0_dp)
        complex(dp) :: w = (1.0_dp, 0.0_dp)
        complex(dp), allocatable :: input_chirp(:)
        complex(dp), allocatable :: output_chirp(:)
        complex(dp), allocatable :: kernel_fft(:)
        complex(dp), allocatable :: inverse_points(:)
    contains
        procedure :: initialize => initialize_czt_plan
        procedure, private :: execute_complex => execute_czt_complex
        procedure, private :: execute_real => execute_czt_real
        procedure, private :: execute_complex_2d => execute_czt_complex_2d
        procedure, private :: execute_real_2d => execute_czt_real_2d
        generic :: transform => execute_complex, execute_real, execute_complex_2d, execute_real_2d
        procedure :: points => czt_plan_points
        procedure :: is_valid => czt_plan_is_valid
        procedure :: input_length => czt_plan_input_length
        procedure :: output_length => czt_plan_output_length
    end type czt_plan

    type, public :: zoom_fft_plan
        private
        type(czt_plan) :: czt_kernel
        real(dp) :: start_frequency = 0.0_dp
        real(dp) :: frequency_step = 0.0_dp
    contains
        procedure, private :: initialize_pair => initialize_zoom_plan_pair
        procedure, private :: initialize_scalar => initialize_zoom_plan_scalar
        generic :: initialize => initialize_pair, initialize_scalar
        procedure, private :: execute_complex => execute_zoom_complex
        procedure, private :: execute_real => execute_zoom_real
        procedure, private :: execute_complex_2d => execute_zoom_complex_2d
        procedure, private :: execute_real_2d => execute_zoom_real_2d
        generic :: transform => execute_complex, execute_real, execute_complex_2d, execute_real_2d
        procedure :: frequencies => zoom_plan_frequencies
        procedure :: is_valid => zoom_plan_is_valid
        procedure :: input_length => zoom_plan_input_length
        procedure :: output_length => zoom_plan_output_length
    end type zoom_fft_plan

contains

    function czt_points(m, w, a, status) result(points)
        ! The CZT evaluates the z-transform at z(k) = a*w**(-k).
        integer, intent(in) :: m !! number of sample points on the spiral
        complex(dp), optional, intent(in) :: w !! successive inverse point ratio
        complex(dp), optional, intent(in) :: a !! starting point, default 1
        integer, optional, intent(out) :: status !! status code
        complex(dp), allocatable :: points(:)
        complex(dp) :: first, step
        real(dp) :: angle
        integer :: k

        allocate(points(0))
        if (present(status)) status = czt_status_invalid_input
        if (m < 1 .or. m > 2000000) return
        first = (1.0_dp, 0.0_dp)
        if (present(a)) first = a
        if (.not. finite_complex(first) .or. abs(first) == 0.0_dp) return
        if (present(w)) then
            step = w
            if (.not. finite_complex(step) .or. abs(step) == 0.0_dp) return
        end if
        deallocate(points)
        allocate(points(m))
        if (present(w)) then
            do k = 0, m-1
                points(k+1) = first * step**(-k)
            end do
        else
            do k = 0, m-1
                angle = 2.0_dp*pi*real(k,dp)/real(m,dp)
                points(k+1) = first * cmplx(cos(angle),sin(angle),dp)
            end do
        end if
        if (.not. all(ieee_is_finite(real(points,dp))) .or. &
            .not. all(ieee_is_finite(aimag(points)))) then
            deallocate(points)
            allocate(points(0))
            if (present(status)) status = czt_status_nonfinite_result
            return
        end if
        if (present(status)) status = czt_status_success
    end function czt_points

    pure logical function finite_complex(x) result(ok)
        complex(dp), intent(in) :: x !! complex number being checked
        ok = ieee_is_finite(real(x,dp)) .and. ieee_is_finite(aimag(x))
    end function finite_complex

    function czt_complex(x, m, w, a, status) result(y)
        complex(dp), intent(in) :: x(:) !! input vector
        integer, intent(in), optional :: m !! number of output points (default size(x))
        complex(dp), intent(in), optional :: w !! complex spiral ratio (default exp(-2*pi*i/m))
        complex(dp), intent(in), optional :: a !! complex starting point (default 1)
        integer, intent(out), optional :: status !! success, invalid input or nonfinite result
        complex(dp), allocatable :: y(:), f(:), g(:)
        complex(dp) :: step, first, ak, q, log_first, value
        real(dp) :: phase, t
        integer :: n, count, nfft, i, j

        allocate(y(0))
        if (present(status)) status = czt_status_invalid_input
        n = size(x)
        count = n
        if (present(m)) count = m
        if (n < 1 .or. count < 1) return
        if (n > huge(n)-count+1) return
        if (.not. all(ieee_is_finite(real(x,dp)))) return
        if (.not. all(ieee_is_finite(aimag(x)))) return
        step = cmplx(cos(-2.0_dp*pi/real(count,dp)), &
            sin(-2.0_dp*pi/real(count,dp)),dp)
        if (present(w)) step = w
        first = cmplx(1.0_dp,0.0_dp,dp)
        if (present(a)) first = a
        if (.not. finite_complex(step) .or. .not. finite_complex(first)) return
        if (abs(step) == 0.0_dp .or. abs(first) == 0.0_dp) return
        if (count > 2000000 .or. n > 2000000) return

        deallocate(y)
        allocate(y(count))
        if (n == 1) then
            y = x(1)
        else if (abs(abs(step)-1.0_dp) > 8.0_dp*epsilon(1.0_dp) .or. &
            abs(abs(first)-1.0_dp) > 8.0_dp*epsilon(1.0_dp)) then
            ! For nonunit contours the FFT chirp may overflow although the
            ! polynomial is representable. Horner evaluation avoids that.
            ak = cmplx(1.0_dp,0.0_dp,dp)
            do i = 1,count
                q = ak/first
                value = x(n)
                do j = n-1,1,-1
                    value = x(j)+q*value
                end do
                y(i) = value
                ak = ak*step
            end do
        else
            nfft = 1
            do while (nfft < n+count-1)
                if (nfft > huge(nfft)/2) then
                    deallocate(y)
                    allocate(y(0))
                    return
                end if
                nfft = 2*nfft
            end do
            allocate(f(nfft),g(nfft))
            f = cmplx(0.0_dp,0.0_dp,dp)
            g = cmplx(0.0_dp,0.0_dp,dp)
            phase = atan2(aimag(step),real(step,dp))
            log_first = log(first)
            do i = 0,n-1
                t = 0.5_dp*real(i,dp)*real(i,dp)
                f(i+1) = x(i+1)*exp(cmplx(0.0_dp,t*phase,dp)-real(i,dp)*log_first)
            end do
            do i = 0,count-1
                t = -0.5_dp*real(i,dp)*real(i,dp)
                g(i+1) = exp(cmplx(0.0_dp,t*phase,dp))
            end do
            do i = 1,n-1
                t = -0.5_dp*real(i,dp)*real(i,dp)
                g(nfft-i+1) = exp(cmplx(0.0_dp,t*phase,dp))
            end do
            call forward_fft(f)
            call forward_fft(g)
            f = conjg(f*g)
            call forward_fft(f)
            f = conjg(f)/real(nfft,dp)
            do i = 0,count-1
                t = 0.5_dp*real(i,dp)*real(i,dp)
                y(i+1) = f(i+1)*exp(cmplx(0.0_dp,t*phase,dp))
            end do
        end if
        if (.not. all(ieee_is_finite(real(y,dp))) .or. &
            .not. all(ieee_is_finite(aimag(y)))) then
            deallocate(y)
            allocate(y(0))
            if (present(status)) status = czt_status_nonfinite_result
            return
        end if
        if (present(status)) status = czt_status_success
    end function czt_complex

    function czt_real(x, m, w, a, status) result(y)
        real(dp), intent(in) :: x(:) !! real input vector
        integer, intent(in), optional :: m !! output count
        complex(dp), intent(in), optional :: w !! spiral ratio
        complex(dp), intent(in), optional :: a !! starting point
        integer, intent(out), optional :: status !! status code
        complex(dp), allocatable :: y(:)
        y = czt_complex(cmplx(x,0.0_dp,dp),m,w,a,status)
    end function czt_real

    function zoom_fft_complex(x, fn, m, fs, endpoint, status) result(y)
        complex(dp), intent(in) :: x(:) !! input vector
        real(dp), intent(in) :: fn(2) !! first frequency and upper frequency limit
        integer, intent(in), optional :: m !! output count, default size(x)
        real(dp), intent(in), optional :: fs !! sampling frequency, default 2
        logical, intent(in), optional :: endpoint !! include fn(2), default false
        integer, intent(out), optional :: status !! status code
        complex(dp), allocatable :: y(:)
        complex(dp) :: w, a
        real(dp) :: rate, width, denom
        integer :: count
        logical :: include_endpoint

        allocate(y(0))
        if (present(status)) status = czt_status_invalid_input
        count = size(x)
        if (present(m)) count = m
        rate = 2.0_dp
        if (present(fs)) rate = fs
        if (count < 1 .or. size(x) < 1) return
        if (rate <= 0.0_dp .or. .not. ieee_is_finite(rate)) return
        if (.not. all(ieee_is_finite(fn))) return
        include_endpoint = .false.
        if (present(endpoint)) include_endpoint = endpoint
        denom = real(count,dp)
        if (include_endpoint .and. count > 1) denom = real(count-1,dp)
        width = (fn(2)-fn(1))/rate
        a = exp(cmplx(0.0_dp, 2.0_dp*pi*fn(1)/rate,dp))
        w = exp(cmplx(0.0_dp, -2.0_dp*pi*width/denom,dp))
        deallocate(y)
        y = czt_complex(x,m=count,w=w,a=a,status=status)
    end function zoom_fft_complex

    function zoom_fft_real(x, fn, m, fs, endpoint, status) result(y)
        real(dp), intent(in) :: x(:) !! input vector
        real(dp), intent(in) :: fn(2) !! lower and upper frequencies
        integer, intent(in), optional :: m !! output count
        real(dp), intent(in), optional :: fs !! sampling frequency
        logical, intent(in), optional :: endpoint !! include upper endpoint
        integer, intent(out), optional :: status !! status code
        complex(dp), allocatable :: y(:)
        y = zoom_fft_complex(cmplx(x,0.0_dp,dp),fn,m,fs,endpoint,status)
    end function zoom_fft_real

    function zoom_fft_complex_scalar(x, fn, m, fs, endpoint, status) result(y)
        complex(dp), intent(in) :: x(:) !! input vector
        real(dp), intent(in) :: fn !! upper frequency with starting frequency zero
        integer, intent(in), optional :: m !! output count
        real(dp), intent(in), optional :: fs !! sampling frequency
        logical, intent(in), optional :: endpoint !! include upper endpoint
        integer, intent(out), optional :: status !! status code
        complex(dp), allocatable :: y(:)
        y = zoom_fft_complex(x,[0.0_dp,fn],m,fs,endpoint,status)
    end function zoom_fft_complex_scalar

    function zoom_fft_real_scalar(x, fn, m, fs, endpoint, status) result(y)
        real(dp), intent(in) :: x(:) !! input vector
        real(dp), intent(in) :: fn !! upper frequency with starting frequency zero
        integer, intent(in), optional :: m !! output count
        real(dp), intent(in), optional :: fs !! sampling frequency
        logical, intent(in), optional :: endpoint !! include upper endpoint
        integer, intent(out), optional :: status !! status code
        complex(dp), allocatable :: y(:)
        y = zoom_fft_complex(cmplx(x,0.0_dp,dp),[0.0_dp,fn],m,fs,endpoint,status)
    end function zoom_fft_real_scalar


    subroutine initialize_czt_plan(self, n, m, w, a, status)
        class(czt_plan), intent(inout) :: self !! transform plan being initialized
        integer, intent(in) :: n !! number of input samples
        integer, intent(in), optional :: m !! number of output samples (default n)
        complex(dp), intent(in), optional :: w !! spiral ratio
        complex(dp), intent(in), optional :: a !! spiral starting point
        integer, intent(out), optional :: status !! status code
        integer :: count, fft_length, k
        real(dp) :: angle, t
        complex(dp) :: first, step, log_first
        complex(dp), allocatable :: g(:)

        self%ready = .false.
        self%n = 0
        self%m = 0
        self%fft_len = 0
        if (allocated(self%input_chirp)) deallocate(self%input_chirp)
        if (allocated(self%output_chirp)) deallocate(self%output_chirp)
        if (allocated(self%kernel_fft)) deallocate(self%kernel_fft)
        if (allocated(self%inverse_points)) deallocate(self%inverse_points)
        if (present(status)) status = czt_status_invalid_input
        count = n
        if (present(m)) count = m
        if (n < 1 .or. count < 1 .or. n > 2000000 .or. count > 2000000) return
        if (n > huge(n)-count+1) return

        step = cmplx(cos(-2.0_dp*pi/real(count,dp)), &
            sin(-2.0_dp*pi/real(count,dp)), dp)
        if (present(w)) step = w
        first = cmplx(1.0_dp,0.0_dp,dp)
        if (present(a)) first = a
        if (.not. finite_complex(step) .or. .not. finite_complex(first)) return
        if (abs(step) == 0.0_dp .or. abs(first) == 0.0_dp) return

        self%n = n
        self%m = count
        self%a = first
        self%w = step
        self%use_fft = abs(abs(step)-1.0_dp) <= 8.0_dp*epsilon(1.0_dp) .and. &
            abs(abs(first)-1.0_dp) <= 8.0_dp*epsilon(1.0_dp)
        if (n > 1 .and. self%use_fft) then
            fft_length = 1
            do while (fft_length < n+count-1)
                if (fft_length > huge(fft_length)/2) return
                fft_length = 2*fft_length
            end do
            self%fft_len = fft_length
            allocate(self%input_chirp(n),self%output_chirp(count),g(fft_length))
            g = cmplx(0.0_dp,0.0_dp,dp)
            angle = atan2(aimag(step),real(step,dp))
            log_first = log(first)
            do k = 0,n-1
                t = 0.5_dp*real(k,dp)*real(k,dp)
                self%input_chirp(k+1) = exp(cmplx(0.0_dp,t*angle,dp)-real(k,dp)*log_first)
            end do
            do k = 0,count-1
                t = 0.5_dp*real(k,dp)*real(k,dp)
                self%output_chirp(k+1) = exp(cmplx(0.0_dp,t*angle,dp))
                g(k+1) = exp(cmplx(0.0_dp,-t*angle,dp))
            end do
            do k = 1,n-1
                t = 0.5_dp*real(k,dp)*real(k,dp)
                g(fft_length-k+1) = exp(cmplx(0.0_dp,-t*angle,dp))
            end do
            call forward_fft(g)
            self%kernel_fft = g
        else if (n > 1) then
            allocate(self%inverse_points(count))
            self%inverse_points(1) = 1.0_dp/first
            do k = 2,count
                self%inverse_points(k) = self%inverse_points(k-1)*step
            end do
            ! Off-unit spirals can exceed the representable range.
            if (.not. all(ieee_is_finite(real(self%inverse_points,dp))) .or. &
                .not. all(ieee_is_finite(aimag(self%inverse_points)))) then
                deallocate(self%inverse_points)
                if (present(status)) status = czt_status_nonfinite_result
                return
            end if
        end if
        self%ready = .true.
        if (present(status)) status = czt_status_success
    end subroutine initialize_czt_plan

    logical function czt_plan_is_valid(self) result(ok)
        class(czt_plan), intent(in) :: self !! transform plan
        ok = self%ready
    end function czt_plan_is_valid

    integer function czt_plan_input_length(self) result(n)
        class(czt_plan), intent(in) :: self !! transform plan
        n = self%n
    end function czt_plan_input_length

    integer function czt_plan_output_length(self) result(m)
        class(czt_plan), intent(in) :: self !! transform plan
        m = self%m
    end function czt_plan_output_length

    function czt_plan_points(self, status) result(z)
        class(czt_plan), intent(in) :: self !! transform plan
        integer, intent(out), optional :: status !! status code
        complex(dp), allocatable :: z(:)
        allocate(z(0))
        if (present(status)) status = czt_status_invalid_input
        if (.not. self%ready) return
        z = czt_points(self%m, self%w, self%a, status)
    end function czt_plan_points

    function execute_czt_complex(self, x, status) result(y)
        class(czt_plan), intent(in) :: self !! initialized transform plan
        complex(dp), intent(in) :: x(:) !! input samples
        integer, intent(out), optional :: status !! status code
        complex(dp), allocatable :: y(:), f(:)
        complex(dp) :: val, q
        integer :: k, j

        allocate(y(0))
        if (present(status)) status = czt_status_invalid_input
        if (.not. self%ready) return
        if (size(x) /= self%n) return
        if (.not. all(ieee_is_finite(real(x,dp))) .or. &
            .not. all(ieee_is_finite(aimag(x)))) return
        deallocate(y)
        allocate(y(self%m))
        if (self%n == 1) then
            y = x(1)
        else if (self%use_fft) then
            allocate(f(self%fft_len))
            f = cmplx(0.0_dp,0.0_dp,dp)
            f(1:self%n) = x*self%input_chirp
            call forward_fft(f)
            f = conjg(f*self%kernel_fft)
            call forward_fft(f)
            y = (conjg(f(1:self%m))/real(self%fft_len,dp))*self%output_chirp
        else
            do k = 1,self%m
                q = self%inverse_points(k)
                val = x(self%n)
                do j = self%n-1,1,-1
                    val = x(j)+q*val
                end do
                y(k) = val
            end do
        end if
        if (.not. all(ieee_is_finite(real(y,dp))) .or. &
            .not. all(ieee_is_finite(aimag(y)))) then
            deallocate(y)
            allocate(y(0))
            if (present(status)) status = czt_status_nonfinite_result
            return
        end if
        if (present(status)) status = czt_status_success
    end function execute_czt_complex

    function execute_czt_real(self, x, status) result(y)
        class(czt_plan), intent(in) :: self !! initialized transform plan
        real(dp), intent(in) :: x(:) !! input samples
        integer, intent(out), optional :: status !! status code
        complex(dp), allocatable :: y(:)
        y = self%execute_complex(cmplx(x,0.0_dp,dp),status)
    end function execute_czt_real

    function execute_czt_complex_2d(self, x, axis, status) result(y)
        class(czt_plan), intent(in) :: self !! initialized transform plan
        complex(dp), intent(in) :: x(:,:) !! input channels, axis 1 or 2
        integer, intent(in), optional :: axis !! transform axis, Fortran one-based
        integer, intent(out), optional :: status !! status code
        complex(dp), allocatable :: y(:,:), segment(:)
        integer :: ax, j, stat

        allocate(y(0,0))
        if (present(status)) status = czt_status_invalid_input
        if (.not. self%ready) return
        ax = 1
        stat = czt_status_success
        if (present(axis)) ax = axis
        if (ax == 1) then
            if (size(x,1) /= self%n) return
            deallocate(y)
            allocate(y(self%m,size(x,2)))
            do j = 1,size(x,2)
                segment = self%execute_complex(x(:,j),stat)
                if (stat /= czt_status_success) exit
                y(:,j) = segment
            end do
        else if (ax == 2) then
            if (size(x,2) /= self%n) return
            deallocate(y)
            allocate(y(size(x,1),self%m))
            do j = 1,size(x,1)
                segment = self%execute_complex(x(j,:),stat)
                if (stat /= czt_status_success) exit
                y(j,:) = segment
            end do
        else
            return
        end if
        if (stat /= czt_status_success) then
            deallocate(y)
            allocate(y(0,0))
            if (present(status)) status = stat
            return
        end if
        if (present(status)) status = czt_status_success
    end function execute_czt_complex_2d

    function execute_czt_real_2d(self, x, axis, status) result(y)
        class(czt_plan), intent(in) :: self !! initialized transform plan
        real(dp), intent(in) :: x(:,:) !! input channels, axis 1 or 2
        integer, intent(in), optional :: axis !! transform axis, Fortran one-based
        integer, intent(out), optional :: status !! status code
        complex(dp), allocatable :: y(:,:)
        y = self%execute_complex_2d(cmplx(x,0.0_dp,dp),axis,status)
    end function execute_czt_real_2d

    subroutine initialize_zoom_plan_pair(self, n, fn, m, fs, endpoint, status)
        class(zoom_fft_plan), intent(inout) :: self !! zoom plan
        integer, intent(in) :: n !! input length
        real(dp), intent(in) :: fn(2) !! start and end frequencies
        integer, intent(in), optional :: m !! output length
        real(dp), intent(in), optional :: fs !! sampling frequency
        logical, intent(in), optional :: endpoint !! include final frequency
        integer, intent(out), optional :: status !! status code
        integer :: count, stat
        real(dp) :: rate, denom
        complex(dp) :: step, first
        logical :: with_endpoint

        if (present(status)) status = czt_status_invalid_input
        count = n
        if (present(m)) count = m
        rate = 2.0_dp
        if (present(fs)) rate = fs
        ! Reset the internal plan even when arguments are invalid.
        call self%czt_kernel%initialize(0)
        if (n < 1 .or. count < 1 .or. .not. ieee_is_finite(rate)) return
        if (rate <= 0.0_dp .or. .not. all(ieee_is_finite(fn))) return
        with_endpoint = .false.
        if (present(endpoint)) with_endpoint = endpoint
        denom = real(count,dp)
        if (with_endpoint .and. count > 1) denom = real(count-1,dp)
        self%start_frequency = fn(1)
        self%frequency_step = (fn(2)-fn(1))/denom
        if (.not. ieee_is_finite(self%frequency_step)) return
        if (.not. ieee_is_finite(2.0_dp*pi*fn(1)/rate)) return
        if (.not. ieee_is_finite(2.0_dp*pi*self%frequency_step/rate)) return
        first = exp(cmplx(0.0_dp,2.0_dp*pi*fn(1)/rate,dp))
        step = exp(cmplx(0.0_dp,-2.0_dp*pi*self%frequency_step/rate,dp))
        call self%czt_kernel%initialize(n,m=count,w=step,a=first,status=stat)
        if (present(status)) status = stat
    end subroutine initialize_zoom_plan_pair

    subroutine initialize_zoom_plan_scalar(self, n, fn, m, fs, endpoint, status)
        class(zoom_fft_plan), intent(inout) :: self !! zoom plan
        integer, intent(in) :: n !! input length
        real(dp), intent(in) :: fn !! upper frequency, lower is zero
        integer, intent(in), optional :: m !! output length
        real(dp), intent(in), optional :: fs !! sampling frequency
        logical, intent(in), optional :: endpoint !! include final frequency
        integer, intent(out), optional :: status !! status code
        call self%initialize_pair(n,[0.0_dp,fn],m,fs,endpoint,status)
    end subroutine initialize_zoom_plan_scalar

    logical function zoom_plan_is_valid(self) result(ok)
        class(zoom_fft_plan), intent(in) :: self !! zoom plan
        ok = self%czt_kernel%is_valid()
    end function zoom_plan_is_valid

    integer function zoom_plan_input_length(self) result(n)
        class(zoom_fft_plan), intent(in) :: self !! zoom plan
        n = self%czt_kernel%input_length()
    end function zoom_plan_input_length

    integer function zoom_plan_output_length(self) result(m)
        class(zoom_fft_plan), intent(in) :: self !! zoom plan
        m = self%czt_kernel%output_length()
    end function zoom_plan_output_length

    function zoom_plan_frequencies(self) result(frequencies)
        class(zoom_fft_plan), intent(in) :: self !! zoom plan
        real(dp), allocatable :: frequencies(:)
        integer :: k, n
        n = 0
        if (self%is_valid()) n = self%output_length()
        allocate(frequencies(n))
        do k = 1,n
            frequencies(k) = self%start_frequency+real(k-1,dp)*self%frequency_step
        end do
    end function zoom_plan_frequencies

    function execute_zoom_complex(self, x, status) result(y)
        class(zoom_fft_plan), intent(in) :: self !! zoom plan
        complex(dp), intent(in) :: x(:) !! input samples
        integer, intent(out), optional :: status !! status code
        complex(dp), allocatable :: y(:)
        y = self%czt_kernel%transform(x,status)
    end function execute_zoom_complex

    function execute_zoom_real(self, x, status) result(y)
        class(zoom_fft_plan), intent(in) :: self !! zoom plan
        real(dp), intent(in) :: x(:) !! input samples
        integer, intent(out), optional :: status !! status code
        complex(dp), allocatable :: y(:)
        y = self%czt_kernel%transform(x,status)
    end function execute_zoom_real

    function execute_zoom_complex_2d(self, x, axis, status) result(y)
        class(zoom_fft_plan), intent(in) :: self !! zoom plan
        complex(dp), intent(in) :: x(:,:) !! complex input channels
        integer, intent(in), optional :: axis !! Fortran one-based transform axis
        integer, intent(out), optional :: status !! status code
        complex(dp), allocatable :: y(:,:)
        y = self%czt_kernel%transform(x,axis,status)
    end function execute_zoom_complex_2d

    function execute_zoom_real_2d(self, x, axis, status) result(y)
        class(zoom_fft_plan), intent(in) :: self !! zoom plan
        real(dp), intent(in) :: x(:,:) !! real input channels
        integer, intent(in), optional :: axis !! Fortran one-based transform axis
        integer, intent(out), optional :: status !! status code
        complex(dp), allocatable :: y(:,:)
        y = self%czt_kernel%transform(x,axis,status)
    end function execute_zoom_real_2d

end module scifort_signal_czt
