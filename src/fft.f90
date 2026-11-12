! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Standalone binary64 transforms from the public SciPy mathematical definitions.
! FFT-backed even/odd extensions implement DCT/DST types I-IV; see provenance.
module scifort_fft
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : fft_plan
    implicit none
    private
    integer, parameter, public :: fft_status_success = 0
    integer, parameter, public :: fft_status_invalid_input = 1
    integer, parameter, public :: fft_status_setup_failure = 2
    public :: dp, fft_plan, fft, ifft, rfft, irfft, hfft, ihfft
    public :: dct, idct, dst, idst, fftfreq, rfftfreq, fftshift, ifftshift
    interface fft
        module procedure fft_real, fft_complex
    end interface fft
    interface ifft
        module procedure ifft_real, ifft_complex
    end interface ifft
    interface irfft
        module procedure irfft_real, irfft_complex
    end interface irfft
    interface hfft
        module procedure hfft_real, hfft_complex
    end interface hfft
    interface dct
        module procedure dct_real, dct_complex
    end interface dct
    interface idct
        module procedure idct_real, idct_complex
    end interface idct
    interface dst
        module procedure dst_real, dst_complex
    end interface dst
    interface idst
        module procedure idst_real, idst_complex
    end interface idst
    interface fftshift
        module procedure fftshift_real, fftshift_complex
    end interface fftshift
    interface ifftshift
        module procedure ifftshift_real, ifftshift_complex
    end interface ifftshift

contains

    function fft_real(x,n,norm,status) result(y)
        implicit none
        real(dp), intent(in) :: x(:) !! finite input samples; strided sections supported
        integer, intent(in), optional :: n !! positive transform length; crop or zero-pad
        character(*), intent(in), optional :: norm !! backward (default), forward, or ortho
        integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 setup failure
        complex(dp), allocatable :: y(:)
        integer :: code
        y = complex_transform(cmplx(x,0.0_dp,dp),n,norm,.false.,code)
        if (present(status)) status = code
    end function fft_real

    function fft_complex(x,n,norm,status) result(y)
        implicit none
        complex(dp), intent(in) :: x(:) !! finite input samples; strided sections supported
        integer, intent(in), optional :: n !! positive transform length; crop or zero-pad
        character(*), intent(in), optional :: norm !! backward (default), forward, or ortho
        integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 setup failure
        complex(dp), allocatable :: y(:)
        integer :: code
        y = complex_transform(x,n,norm,.false.,code)
        if (present(status)) status = code
    end function fft_complex

    function ifft_real(x,n,norm,status) result(y)
        implicit none
        real(dp), intent(in) :: x(:) !! finite input samples; strided sections supported
        integer, intent(in), optional :: n !! positive transform length; crop or zero-pad
        character(*), intent(in), optional :: norm !! backward (default), forward, or ortho
        integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 setup failure
        complex(dp), allocatable :: y(:)
        integer :: code
        y = complex_transform(cmplx(x,0.0_dp,dp),n,norm,.true.,code)
        if (present(status)) status = code
    end function ifft_real

    function ifft_complex(x,n,norm,status) result(y)
        implicit none
        complex(dp), intent(in) :: x(:) !! finite input samples; strided sections supported
        integer, intent(in), optional :: n !! positive transform length; crop or zero-pad
        character(*), intent(in), optional :: norm !! backward (default), forward, or ortho
        integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 setup failure
        complex(dp), allocatable :: y(:)
        integer :: code
        y = complex_transform(x,n,norm,.true.,code)
        if (present(status)) status = code
    end function ifft_complex

    function complex_transform(x,n,norm,inverse,code) result(y)
        implicit none
        complex(dp), intent(in) :: x(:) !! finite complex input
        integer, intent(in), optional :: n !! requested positive output length
        character(*), intent(in), optional :: norm !! normalization name
        logical, intent(in) :: inverse !! select the positive-exponent transform
        integer, intent(out) :: code !! completion status
        complex(dp), allocatable :: y(:), work(:)
        integer :: count, mode, ierr, copied
        allocate(y(0))
        code = fft_status_invalid_input
        count = size(x)
        if (present(n)) count = n
        mode = normalization(norm)
        if (count < 1 .or. mode < 0) return
        if (.not. finite_complex(x)) return
        allocate(work(count),stat=ierr)
        code = fft_status_setup_failure
        if (ierr /= 0) return
        work = (0.0_dp,0.0_dp)
        copied = min(count,size(x))
        work(:copied) = x(:copied)
        call apply_fft(work,inverse,mode,code)
        if (code /= fft_status_success) return
        call move_alloc(work,y)
    end function complex_transform

    subroutine apply_fft(values,inverse,mode,code)
        implicit none
        complex(dp), intent(inout) :: values(:) !! input samples, transformed in place
        logical, intent(in) :: inverse !! positive exponent when true
        integer, intent(in) :: mode !! 0 backward, 1 ortho, 2 forward
        integer, intent(out) :: code !! completion status
        type(fft_plan) :: plan
        real(dp) :: factor
        integer :: setup
        code = fft_status_setup_failure
        call plan%initialize(size(values),setup)
        if (setup /= 0) return
        factor = transform_scale(real(size(values),dp),mode,inverse)
        values = values * factor
        if (inverse) values = conjg(values)
        call plan%execute(values,setup)
        if (setup /= 0) return
        if (inverse) values = conjg(values)
        code = fft_status_success
    end subroutine apply_fft

    pure integer function normalization(norm) result(mode)
        implicit none
        character(*), intent(in), optional :: norm !! normalization name, default backward
        mode = 0
        if (.not. present(norm)) return
        select case (norm)
        case ('backward')
            mode = 0
        case ('ortho')
            mode = 1
        case ('forward')
            mode = 2
        case default
            mode = -1
        end select
    end function normalization

    pure function transform_scale(length,mode,inverse) result(factor)
        implicit none
        real(dp), intent(in) :: length !! positive logical transform length
        integer, intent(in) :: mode !! validated normalization code
        logical, intent(in) :: inverse !! select inverse normalization
        real(dp) :: factor
        factor = 1.0_dp
        if (mode == 1) then
            factor = 1.0_dp/sqrt(length)
        else if ((mode == 0 .and. inverse) .or. (mode == 2 .and. .not. inverse)) then
            factor = 1.0_dp/length
        end if
    end function transform_scale

    pure logical function finite_complex(x) result(valid)
        implicit none
        complex(dp), intent(in) :: x(:) !! complex samples to validate
        valid = all(ieee_is_finite(real(x,dp))) .and. all(ieee_is_finite(aimag(x)))
    end function finite_complex

    function rfft(x,n,norm,status) result(y)
        implicit none
        real(dp), intent(in) :: x(:) !! finite real input samples
        integer, intent(in), optional :: n !! positive full transform length; crop or zero-pad
        character(*), intent(in), optional :: norm !! backward (default), forward, or ortho
        integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 setup failure
        complex(dp), allocatable :: y(:), full(:)
        integer :: code
        allocate(y(0))
        full = complex_transform(cmplx(x,0.0_dp,dp),n,norm,.false.,code)
        if (code == fft_status_success) y = full(:size(full)/2+1)
        if (present(status)) status = code
    end function rfft

    function ihfft(x,n,norm,status) result(y)
        implicit none
        real(dp), intent(in) :: x(:) !! finite real input samples
        integer, intent(in), optional :: n !! positive full transform length; crop or zero-pad
        character(*), intent(in), optional :: norm !! backward (default), forward, or ortho
        integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 setup failure
        complex(dp), allocatable :: y(:), full(:)
        integer :: code
        allocate(y(0))
        full = complex_transform(cmplx(x,0.0_dp,dp),n,norm,.true.,code)
        if (code == fft_status_success) y = full(:size(full)/2+1)
        if (present(status)) status = code
    end function ihfft

    function irfft_real(x,n,norm,status) result(y)
        implicit none
        real(dp), intent(in) :: x(:) !! finite nonnegative-index half of a Hermitian sequence
        integer, intent(in), optional :: n !! positive full output length; default 2*(size(x)-1)
        character(*), intent(in), optional :: norm !! backward (default), forward, or ortho
        integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 setup failure
        real(dp), allocatable :: y(:)
        integer :: code
        y = hermitian_transform(cmplx(x,0.0_dp,dp),n,norm,.true.,code)
        if (present(status)) status = code
    end function irfft_real

    function irfft_complex(x,n,norm,status) result(y)
        implicit none
        complex(dp), intent(in) :: x(:) !! finite nonnegative-index half of a Hermitian sequence
        integer, intent(in), optional :: n !! positive full output length; default 2*(size(x)-1)
        character(*), intent(in), optional :: norm !! backward (default), forward, or ortho
        integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 setup failure
        real(dp), allocatable :: y(:)
        integer :: code
        y = hermitian_transform(x,n,norm,.true.,code)
        if (present(status)) status = code
    end function irfft_complex

    function hfft_real(x,n,norm,status) result(y)
        implicit none
        real(dp), intent(in) :: x(:) !! finite nonnegative-index half of a Hermitian sequence
        integer, intent(in), optional :: n !! positive full output length; default 2*(size(x)-1)
        character(*), intent(in), optional :: norm !! backward (default), forward, or ortho
        integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 setup failure
        real(dp), allocatable :: y(:)
        integer :: code
        y = hermitian_transform(cmplx(x,0.0_dp,dp),n,norm,.false.,code)
        if (present(status)) status = code
    end function hfft_real

    function hfft_complex(x,n,norm,status) result(y)
        implicit none
        complex(dp), intent(in) :: x(:) !! finite nonnegative-index half of a Hermitian sequence
        integer, intent(in), optional :: n !! positive full output length; default 2*(size(x)-1)
        character(*), intent(in), optional :: norm !! backward (default), forward, or ortho
        integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 setup failure
        real(dp), allocatable :: y(:)
        integer :: code
        y = hermitian_transform(x,n,norm,.false.,code)
        if (present(status)) status = code
    end function hfft_complex

    function hermitian_transform(x,n,norm,inverse,code) result(y)
        implicit none
        complex(dp), intent(in) :: x(:) !! supplied half sequence, cropped or padded
        integer, intent(in), optional :: n !! full output length; required for odd lengths
        character(*), intent(in), optional :: norm !! normalization name
        logical, intent(in) :: inverse !! select irfft rather than hfft
        integer, intent(out) :: code !! completion status
        real(dp), allocatable :: y(:)
        complex(dp), allocatable :: work(:)
        integer :: count, mode, copied, j, ierr
        allocate(y(0))
        code = fft_status_invalid_input
        if (present(n)) then
            count = n
        else
            if (size(x) < 2) return
            if (size(x)-1 > huge(count)/2) return
            count = 2*(size(x)-1)
        end if
        mode = normalization(norm)
        if (count < 1 .or. mode < 0) return
        if (.not. finite_complex(x)) return
        allocate(work(count),stat=ierr)
        code = fft_status_setup_failure
        if (ierr /= 0) return
        work = (0.0_dp,0.0_dp)
        copied = min(size(x),count/2+1)
        work(:copied) = x(:copied)
        ! DC, and even-length Nyquist, are self-conjugate and therefore real.
        work(1) = cmplx(real(work(1),dp),0.0_dp,dp)
        if (mod(count,2) == 0) work(count/2+1) = cmplx(real(work(count/2+1),dp),0.0_dp,dp)
        do j = 2, (count+1)/2
            work(count-j+2) = conjg(work(j))
        end do
        call apply_fft(work,inverse,mode,code)
        if (code == fft_status_success) y = real(work,dp)
    end function hermitian_transform

    function dct_real(x,type,n,norm,orthogonalize,status) result(y)
        implicit none
        real(dp), intent(in) :: x(:) !! finite real or complex input samples
        integer, intent(in), optional :: type !! transform type 1, 2 (default), 3, or 4
        integer, intent(in), optional :: n !! positive output length; crop or zero-pad
        character(*), intent(in), optional :: norm !! backward (default), forward, or ortho
        logical, intent(in), optional :: orthogonalize !! defaults true only for norm=ortho
        integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 setup failure
        real(dp), allocatable :: y(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code
        transformed = r2r_transform(cmplx(x,0.0_dp,dp),type,n,norm,orthogonalize,.false.,.false.,code)
        y = real(transformed,dp)
        if (present(status)) status = code
    end function dct_real

    function dct_complex(x,type,n,norm,orthogonalize,status) result(y)
        implicit none
        complex(dp), intent(in) :: x(:) !! finite real or complex input samples
        integer, intent(in), optional :: type !! transform type 1, 2 (default), 3, or 4
        integer, intent(in), optional :: n !! positive output length; crop or zero-pad
        character(*), intent(in), optional :: norm !! backward (default), forward, or ortho
        logical, intent(in), optional :: orthogonalize !! defaults true only for norm=ortho
        integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 setup failure
        complex(dp), allocatable :: y(:)
        integer :: code
        y = r2r_transform(x,type,n,norm,orthogonalize,.false.,.false.,code)
        if (present(status)) status = code
    end function dct_complex

    function idct_real(x,type,n,norm,orthogonalize,status) result(y)
        implicit none
        real(dp), intent(in) :: x(:) !! finite real or complex input samples
        integer, intent(in), optional :: type !! transform type 1, 2 (default), 3, or 4
        integer, intent(in), optional :: n !! positive output length; crop or zero-pad
        character(*), intent(in), optional :: norm !! backward (default), forward, or ortho
        logical, intent(in), optional :: orthogonalize !! defaults true only for norm=ortho
        integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 setup failure
        real(dp), allocatable :: y(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code
        transformed = r2r_transform(cmplx(x,0.0_dp,dp),type,n,norm,orthogonalize,.false.,.true.,code)
        y = real(transformed,dp)
        if (present(status)) status = code
    end function idct_real

    function idct_complex(x,type,n,norm,orthogonalize,status) result(y)
        implicit none
        complex(dp), intent(in) :: x(:) !! finite real or complex input samples
        integer, intent(in), optional :: type !! transform type 1, 2 (default), 3, or 4
        integer, intent(in), optional :: n !! positive output length; crop or zero-pad
        character(*), intent(in), optional :: norm !! backward (default), forward, or ortho
        logical, intent(in), optional :: orthogonalize !! defaults true only for norm=ortho
        integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 setup failure
        complex(dp), allocatable :: y(:)
        integer :: code
        y = r2r_transform(x,type,n,norm,orthogonalize,.false.,.true.,code)
        if (present(status)) status = code
    end function idct_complex

    function dst_real(x,type,n,norm,orthogonalize,status) result(y)
        implicit none
        real(dp), intent(in) :: x(:) !! finite real or complex input samples
        integer, intent(in), optional :: type !! transform type 1, 2 (default), 3, or 4
        integer, intent(in), optional :: n !! positive output length; crop or zero-pad
        character(*), intent(in), optional :: norm !! backward (default), forward, or ortho
        logical, intent(in), optional :: orthogonalize !! defaults true only for norm=ortho
        integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 setup failure
        real(dp), allocatable :: y(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code
        transformed = r2r_transform(cmplx(x,0.0_dp,dp),type,n,norm,orthogonalize,.true.,.false.,code)
        y = real(transformed,dp)
        if (present(status)) status = code
    end function dst_real

    function dst_complex(x,type,n,norm,orthogonalize,status) result(y)
        implicit none
        complex(dp), intent(in) :: x(:) !! finite real or complex input samples
        integer, intent(in), optional :: type !! transform type 1, 2 (default), 3, or 4
        integer, intent(in), optional :: n !! positive output length; crop or zero-pad
        character(*), intent(in), optional :: norm !! backward (default), forward, or ortho
        logical, intent(in), optional :: orthogonalize !! defaults true only for norm=ortho
        integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 setup failure
        complex(dp), allocatable :: y(:)
        integer :: code
        y = r2r_transform(x,type,n,norm,orthogonalize,.true.,.false.,code)
        if (present(status)) status = code
    end function dst_complex

    function idst_real(x,type,n,norm,orthogonalize,status) result(y)
        implicit none
        real(dp), intent(in) :: x(:) !! finite real or complex input samples
        integer, intent(in), optional :: type !! transform type 1, 2 (default), 3, or 4
        integer, intent(in), optional :: n !! positive output length; crop or zero-pad
        character(*), intent(in), optional :: norm !! backward (default), forward, or ortho
        logical, intent(in), optional :: orthogonalize !! defaults true only for norm=ortho
        integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 setup failure
        real(dp), allocatable :: y(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code
        transformed = r2r_transform(cmplx(x,0.0_dp,dp),type,n,norm,orthogonalize,.true.,.true.,code)
        y = real(transformed,dp)
        if (present(status)) status = code
    end function idst_real

    function idst_complex(x,type,n,norm,orthogonalize,status) result(y)
        implicit none
        complex(dp), intent(in) :: x(:) !! finite real or complex input samples
        integer, intent(in), optional :: type !! transform type 1, 2 (default), 3, or 4
        integer, intent(in), optional :: n !! positive output length; crop or zero-pad
        character(*), intent(in), optional :: norm !! backward (default), forward, or ortho
        logical, intent(in), optional :: orthogonalize !! defaults true only for norm=ortho
        integer, intent(out), optional :: status !! 0 success, 1 invalid input, 2 setup failure
        complex(dp), allocatable :: y(:)
        integer :: code
        y = r2r_transform(x,type,n,norm,orthogonalize,.true.,.true.,code)
        if (present(status)) status = code
    end function idst_complex

    function r2r_transform(x,type,n,norm,orthogonalize,sine,inverse,code) result(y)
        implicit none
        complex(dp), intent(in) :: x(:) !! finite input, padded or cropped before weighting
        integer, intent(in), optional :: type !! requested transform type, default 2
        integer, intent(in), optional :: n !! requested positive output length
        character(*), intent(in), optional :: norm !! normalization name
        logical, intent(in), optional :: orthogonalize !! endpoint weighting selection
        logical, intent(in) :: sine !! select DST rather than DCT
        logical, intent(in) :: inverse !! select inverse for the requested type
        integer, intent(out) :: code !! completion status
        complex(dp), allocatable :: y(:), data(:), work(:), output(:)
        complex(dp) :: value
        integer :: count, actual, mode, logical_size, length, j, index, frequency, ierr
        logical :: orthogonal
        real(dp) :: factor, parity
        allocate(y(0))
        code = fft_status_invalid_input
        count = size(x)
        if (present(n)) count = n
        actual = 2
        if (present(type)) actual = type
        mode = normalization(norm)
        if (count < 1 .or. count > huge(count)/8) return
        if (actual < 1 .or. actual > 4 .or. mode < 0) return
        if (.not. sine .and. actual == 1 .and. count == 1) return
        if (.not. finite_complex(x)) return
        orthogonal = mode == 1
        if (present(orthogonalize)) orthogonal = orthogonalize
        if (actual == 1) then
            if (sine) then
                logical_size = 2*(count+1)
            else
                logical_size = 2*(count-1)
            end if
        else
            logical_size = 2*count
        end if
        factor = transform_scale(real(logical_size,dp),mode,inverse)
        if (inverse .and. (actual == 2 .or. actual == 3)) actual = 5-actual
        select case (actual)
        case (1)
            length = logical_size
        case (2,3)
            length = 4*count
        case (4)
            length = 8*count
        end select
        allocate(data(count),work(length),output(count),stat=ierr)
        code = fft_status_setup_failure
        if (ierr /= 0) return
        data = (0.0_dp,0.0_dp)
        data(:min(count,size(x))) = x(:min(count,size(x))) * factor
        if (orthogonal) then
            if (.not. sine .and. actual == 1) then
                data(1) = data(1)*sqrt(2.0_dp)
                data(count) = data(count)*sqrt(2.0_dp)
            else if (.not. sine .and. actual == 3) then
                data(1) = data(1)*sqrt(2.0_dp)
            else if (sine .and. actual == 3) then
                data(count) = data(count)*sqrt(2.0_dp)
            end if
        end if
        work = (0.0_dp,0.0_dp)
        parity = 1.0_dp
        if (sine) parity = -1.0_dp
        do j = 0, count-1
            select case (actual)
            case (1,3)
                index = j
                if (sine) index = j+1
            case (2,4)
                index = 2*j+1
            end select
            value = data(j+1)
            if (sine .and. actual == 3 .and. j == count-1) value = value*0.5_dp
            work(index+1) = value
            if (index /= 0 .and. 2*index /= length) work(length-index+1) = parity*value
        end do
        call apply_fft(work,.false.,0,code)
        if (code /= fft_status_success) return
        do j = 0, count-1
            select case (actual)
            case (1,2)
                frequency = j
                if (sine) frequency = j+1
            case (3,4)
                frequency = 2*j+1
            end select
            output(j+1) = work(frequency+1)
            if (sine) output(j+1) = cmplx(0.0_dp,1.0_dp,dp)*output(j+1)
        end do
        if (orthogonal) then
            if (.not. sine .and. actual == 1) then
                output(1) = output(1)/sqrt(2.0_dp)
                output(count) = output(count)/sqrt(2.0_dp)
            else if (.not. sine .and. actual == 2) then
                output(1) = output(1)/sqrt(2.0_dp)
            else if (sine .and. actual == 2) then
                output(count) = output(count)/sqrt(2.0_dp)
            end if
        end if
        call move_alloc(output,y)
    end function r2r_transform

    function fftfreq(n,d,status) result(y)
        implicit none
        integer, intent(in) :: n !! positive full transform length
        real(dp), intent(in), optional :: d !! finite nonzero sample spacing, default 1
        integer, intent(out), optional :: status !! 0 success, 1 invalid/unrepresentable, 2 setup failure
        real(dp), allocatable :: y(:)
        y = frequency_bins(n,d,.false.,status)
    end function fftfreq

    function rfftfreq(n,d,status) result(y)
        implicit none
        integer, intent(in) :: n !! positive full transform length
        real(dp), intent(in), optional :: d !! finite nonzero sample spacing, default 1
        integer, intent(out), optional :: status !! 0 success, 1 invalid/unrepresentable, 2 setup failure
        real(dp), allocatable :: y(:)
        y = frequency_bins(n,d,.true.,status)
    end function rfftfreq

    function frequency_bins(n,d,half,status) result(y)
        implicit none
        integer, intent(in) :: n !! positive full transform length
        real(dp), intent(in), optional :: d !! finite nonzero sample spacing
        logical, intent(in) :: half !! keep the nonnegative-frequency half
        integer, intent(out), optional :: status !! completion status
        real(dp), allocatable :: y(:), work(:)
        real(dp) :: spacing, unit, ratio
        integer :: count, j, bin, power, ierr
        allocate(y(0))
        if (present(status)) status = fft_status_invalid_input
        if (n < 1) return
        spacing = 1.0_dp
        if (present(d)) spacing = d
        if (.not. ieee_is_finite(spacing)) return
        if (spacing == 0.0_dp) return
        count = n
        if (half) count = n/2+1
        allocate(work(count),stat=ierr)
        if (present(status)) status = fft_status_setup_failure
        if (ierr /= 0) return
        do j = 0, count-1
            bin = j
            if (.not. half) then
                if (j > (n-1)/2) bin = j-n
            end if
            unit = real(bin,dp)/real(n,dp)
            work(j+1) = 0.0_dp
            if (unit == 0.0_dp) cycle
            ratio = fraction(unit)/fraction(spacing)
            power = exponent(unit)-exponent(spacing)+exponent(ratio)
            if (power > maxexponent(unit)) then
                if (present(status)) status = fft_status_invalid_input
                return
            end if
            work(j+1) = scale(fraction(ratio),power)
        end do
        call move_alloc(work,y)
        if (present(status)) status = fft_status_success
    end function frequency_bins

    function fftshift_real(x) result(y)
        implicit none
        real(dp), intent(in) :: x(:) !! sequence to permute; nonfinite entries are preserved
        real(dp), allocatable :: y(:)
        y = cshift(x,-(size(x)/2))
    end function fftshift_real

    function fftshift_complex(x) result(y)
        implicit none
        complex(dp), intent(in) :: x(:) !! sequence to permute; nonfinite entries are preserved
        complex(dp), allocatable :: y(:)
        y = cshift(x,-(size(x)/2))
    end function fftshift_complex

    function ifftshift_real(x) result(y)
        implicit none
        real(dp), intent(in) :: x(:) !! sequence to permute; nonfinite entries are preserved
        real(dp), allocatable :: y(:)
        y = cshift(x,(size(x)/2))
    end function ifftshift_real

    function ifftshift_complex(x) result(y)
        implicit none
        complex(dp), intent(in) :: x(:) !! sequence to permute; nonfinite entries are preserved
        complex(dp), allocatable :: y(:)
        y = cshift(x,(size(x)/2))
    end function ifftshift_complex

end module scifort_fft
