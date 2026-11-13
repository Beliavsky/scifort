! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Standalone binary64 transforms from the public SciPy mathematical definitions.
! FFT-backed even/odd extensions implement DCT/DST types I-IV; see provenance.
module scifort_fft
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
    use, intrinsic :: iso_c_binding, only : c_loc, c_f_pointer
    use, intrinsic :: iso_fortran_env, only : int64
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

    public :: fftn, ifftn, rfftn, irfftn, hfftn, ihfftn, dctn, idctn, dstn, idstn
    interface fftn
        module procedure fftn_real, fftn_complex
    end interface fftn
    interface ifftn
        module procedure ifftn_real, ifftn_complex
    end interface ifftn
    interface rfftn
        module procedure rfftn_real
    end interface rfftn
    interface irfftn
        module procedure irfftn_real, irfftn_complex
    end interface irfftn
    interface hfftn
        module procedure hfftn_real, hfftn_complex
    end interface hfftn
    interface ihfftn
        module procedure ihfftn_real
    end interface ihfftn
    interface dctn
        module procedure dctn_real, dctn_complex
    end interface dctn
    interface idctn
        module procedure idctn_real, idctn_complex
    end interface idctn
    interface dstn
        module procedure dstn_real, dstn_complex
    end interface dstn
    interface idstn
        module procedure idstn_real, idstn_complex
    end interface idstn

contains

    subroutine fftn_real(x,y,s,axes,norm,status)
        implicit none
        real(dp), contiguous, target, intent(in) :: x(..) !! finite rank 1-15 input
        complex(dp), contiguous, target, intent(out) :: y(..) !! preallocated expected result shape
        integer, optional, intent(in) :: s(:) !! positive full lengths in axes order
        integer, optional, intent(in) :: axes(:) !! distinct one-based axes; default last size(s) or all
        character(*), optional, intent(in) :: norm !! backward (default), forward, or ortho
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 setup failure
        real(dp), pointer :: input(:)
        complex(dp), pointer :: output(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code

        if (present(status)) status = fft_status_invalid_input
        if (rank(y) < 1 .or. rank(y) > 15) return
        if (size(y,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(y),output,[size(y)])
        output = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
        if (rank(x) /= rank(y)) return
        if (size(x,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(x),input,[size(x)])
        call nd_transform(cmplx(input,0.0_dp,dp),shape(x),shape(y),transformed, &
            1,s,axes,norm,code)
        if (code == fft_status_success) output = transformed
        if (present(status)) status = code
    end subroutine fftn_real

    subroutine fftn_complex(x,y,s,axes,norm,status)
        implicit none
        complex(dp), contiguous, target, intent(in) :: x(..) !! finite rank 1-15 input
        complex(dp), contiguous, target, intent(out) :: y(..) !! preallocated expected result shape
        integer, optional, intent(in) :: s(:) !! positive full lengths in axes order
        integer, optional, intent(in) :: axes(:) !! distinct one-based axes; default last size(s) or all
        character(*), optional, intent(in) :: norm !! backward (default), forward, or ortho
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 setup failure
        complex(dp), pointer :: input(:)
        complex(dp), pointer :: output(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code

        if (present(status)) status = fft_status_invalid_input
        if (rank(y) < 1 .or. rank(y) > 15) return
        if (size(y,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(y),output,[size(y)])
        output = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
        if (rank(x) /= rank(y)) return
        if (size(x,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(x),input,[size(x)])
        call nd_transform(input,shape(x),shape(y),transformed, &
            1,s,axes,norm,code)
        if (code == fft_status_success) output = transformed
        if (present(status)) status = code
    end subroutine fftn_complex

    subroutine ifftn_real(x,y,s,axes,norm,status)
        implicit none
        real(dp), contiguous, target, intent(in) :: x(..) !! finite rank 1-15 input
        complex(dp), contiguous, target, intent(out) :: y(..) !! preallocated expected result shape
        integer, optional, intent(in) :: s(:) !! positive full lengths in axes order
        integer, optional, intent(in) :: axes(:) !! distinct one-based axes; default last size(s) or all
        character(*), optional, intent(in) :: norm !! backward (default), forward, or ortho
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 setup failure
        real(dp), pointer :: input(:)
        complex(dp), pointer :: output(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code

        if (present(status)) status = fft_status_invalid_input
        if (rank(y) < 1 .or. rank(y) > 15) return
        if (size(y,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(y),output,[size(y)])
        output = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
        if (rank(x) /= rank(y)) return
        if (size(x,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(x),input,[size(x)])
        call nd_transform(cmplx(input,0.0_dp,dp),shape(x),shape(y),transformed, &
            2,s,axes,norm,code)
        if (code == fft_status_success) output = transformed
        if (present(status)) status = code
    end subroutine ifftn_real

    subroutine ifftn_complex(x,y,s,axes,norm,status)
        implicit none
        complex(dp), contiguous, target, intent(in) :: x(..) !! finite rank 1-15 input
        complex(dp), contiguous, target, intent(out) :: y(..) !! preallocated expected result shape
        integer, optional, intent(in) :: s(:) !! positive full lengths in axes order
        integer, optional, intent(in) :: axes(:) !! distinct one-based axes; default last size(s) or all
        character(*), optional, intent(in) :: norm !! backward (default), forward, or ortho
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 setup failure
        complex(dp), pointer :: input(:)
        complex(dp), pointer :: output(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code

        if (present(status)) status = fft_status_invalid_input
        if (rank(y) < 1 .or. rank(y) > 15) return
        if (size(y,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(y),output,[size(y)])
        output = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
        if (rank(x) /= rank(y)) return
        if (size(x,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(x),input,[size(x)])
        call nd_transform(input,shape(x),shape(y),transformed, &
            2,s,axes,norm,code)
        if (code == fft_status_success) output = transformed
        if (present(status)) status = code
    end subroutine ifftn_complex

    subroutine rfftn_real(x,y,s,axes,norm,status)
        implicit none
        real(dp), contiguous, target, intent(in) :: x(..) !! finite rank 1-15 input
        complex(dp), contiguous, target, intent(out) :: y(..) !! preallocated expected result shape
        integer, optional, intent(in) :: s(:) !! positive full lengths in axes order
        integer, optional, intent(in) :: axes(:) !! distinct one-based axes; default last size(s) or all
        character(*), optional, intent(in) :: norm !! backward (default), forward, or ortho
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 setup failure
        real(dp), pointer :: input(:)
        complex(dp), pointer :: output(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code

        if (present(status)) status = fft_status_invalid_input
        if (rank(y) < 1 .or. rank(y) > 15) return
        if (size(y,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(y),output,[size(y)])
        output = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
        if (rank(x) /= rank(y)) return
        if (size(x,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(x),input,[size(x)])
        call nd_transform(cmplx(input,0.0_dp,dp),shape(x),shape(y),transformed, &
            3,s,axes,norm,code)
        if (code == fft_status_success) output = transformed
        if (present(status)) status = code
    end subroutine rfftn_real

    subroutine irfftn_real(x,y,s,axes,norm,status)
        implicit none
        real(dp), contiguous, target, intent(in) :: x(..) !! finite rank 1-15 input
        real(dp), contiguous, target, intent(out) :: y(..) !! preallocated expected result shape
        integer, optional, intent(in) :: s(:) !! positive full lengths in axes order
        integer, optional, intent(in) :: axes(:) !! distinct one-based axes; default last size(s) or all
        character(*), optional, intent(in) :: norm !! backward (default), forward, or ortho
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 setup failure
        real(dp), pointer :: input(:)
        real(dp), pointer :: output(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code

        if (present(status)) status = fft_status_invalid_input
        if (rank(y) < 1 .or. rank(y) > 15) return
        if (size(y,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(y),output,[size(y)])
        output = ieee_value(0.0_dp,ieee_quiet_nan)
        if (rank(x) /= rank(y)) return
        if (size(x,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(x),input,[size(x)])
        call nd_transform(cmplx(input,0.0_dp,dp),shape(x),shape(y),transformed, &
            4,s,axes,norm,code)
        if (code == fft_status_success) output = real(transformed,dp)
        if (present(status)) status = code
    end subroutine irfftn_real

    subroutine irfftn_complex(x,y,s,axes,norm,status)
        implicit none
        complex(dp), contiguous, target, intent(in) :: x(..) !! finite rank 1-15 input
        real(dp), contiguous, target, intent(out) :: y(..) !! preallocated expected result shape
        integer, optional, intent(in) :: s(:) !! positive full lengths in axes order
        integer, optional, intent(in) :: axes(:) !! distinct one-based axes; default last size(s) or all
        character(*), optional, intent(in) :: norm !! backward (default), forward, or ortho
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 setup failure
        complex(dp), pointer :: input(:)
        real(dp), pointer :: output(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code

        if (present(status)) status = fft_status_invalid_input
        if (rank(y) < 1 .or. rank(y) > 15) return
        if (size(y,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(y),output,[size(y)])
        output = ieee_value(0.0_dp,ieee_quiet_nan)
        if (rank(x) /= rank(y)) return
        if (size(x,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(x),input,[size(x)])
        call nd_transform(input,shape(x),shape(y),transformed, &
            4,s,axes,norm,code)
        if (code == fft_status_success) output = real(transformed,dp)
        if (present(status)) status = code
    end subroutine irfftn_complex

    subroutine hfftn_real(x,y,s,axes,norm,status)
        implicit none
        real(dp), contiguous, target, intent(in) :: x(..) !! finite rank 1-15 input
        real(dp), contiguous, target, intent(out) :: y(..) !! preallocated expected result shape
        integer, optional, intent(in) :: s(:) !! positive full lengths in axes order
        integer, optional, intent(in) :: axes(:) !! distinct one-based axes; default last size(s) or all
        character(*), optional, intent(in) :: norm !! backward (default), forward, or ortho
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 setup failure
        real(dp), pointer :: input(:)
        real(dp), pointer :: output(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code

        if (present(status)) status = fft_status_invalid_input
        if (rank(y) < 1 .or. rank(y) > 15) return
        if (size(y,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(y),output,[size(y)])
        output = ieee_value(0.0_dp,ieee_quiet_nan)
        if (rank(x) /= rank(y)) return
        if (size(x,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(x),input,[size(x)])
        call nd_transform(cmplx(input,0.0_dp,dp),shape(x),shape(y),transformed, &
            5,s,axes,norm,code)
        if (code == fft_status_success) output = real(transformed,dp)
        if (present(status)) status = code
    end subroutine hfftn_real

    subroutine hfftn_complex(x,y,s,axes,norm,status)
        implicit none
        complex(dp), contiguous, target, intent(in) :: x(..) !! finite rank 1-15 input
        real(dp), contiguous, target, intent(out) :: y(..) !! preallocated expected result shape
        integer, optional, intent(in) :: s(:) !! positive full lengths in axes order
        integer, optional, intent(in) :: axes(:) !! distinct one-based axes; default last size(s) or all
        character(*), optional, intent(in) :: norm !! backward (default), forward, or ortho
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 setup failure
        complex(dp), pointer :: input(:)
        real(dp), pointer :: output(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code

        if (present(status)) status = fft_status_invalid_input
        if (rank(y) < 1 .or. rank(y) > 15) return
        if (size(y,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(y),output,[size(y)])
        output = ieee_value(0.0_dp,ieee_quiet_nan)
        if (rank(x) /= rank(y)) return
        if (size(x,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(x),input,[size(x)])
        call nd_transform(input,shape(x),shape(y),transformed, &
            5,s,axes,norm,code)
        if (code == fft_status_success) output = real(transformed,dp)
        if (present(status)) status = code
    end subroutine hfftn_complex

    subroutine ihfftn_real(x,y,s,axes,norm,status)
        implicit none
        real(dp), contiguous, target, intent(in) :: x(..) !! finite rank 1-15 input
        complex(dp), contiguous, target, intent(out) :: y(..) !! preallocated expected result shape
        integer, optional, intent(in) :: s(:) !! positive full lengths in axes order
        integer, optional, intent(in) :: axes(:) !! distinct one-based axes; default last size(s) or all
        character(*), optional, intent(in) :: norm !! backward (default), forward, or ortho
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 setup failure
        real(dp), pointer :: input(:)
        complex(dp), pointer :: output(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code

        if (present(status)) status = fft_status_invalid_input
        if (rank(y) < 1 .or. rank(y) > 15) return
        if (size(y,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(y),output,[size(y)])
        output = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
        if (rank(x) /= rank(y)) return
        if (size(x,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(x),input,[size(x)])
        call nd_transform(cmplx(input,0.0_dp,dp),shape(x),shape(y),transformed, &
            6,s,axes,norm,code)
        if (code == fft_status_success) output = transformed
        if (present(status)) status = code
    end subroutine ihfftn_real

    subroutine dctn_real(x,y,type,s,axes,norm,orthogonalize,status)
        implicit none
        real(dp), contiguous, target, intent(in) :: x(..) !! finite rank 1-15 input
        real(dp), contiguous, target, intent(out) :: y(..) !! preallocated expected result shape
        integer, optional, intent(in) :: type !! type 1, 2 (default), 3, or 4
        integer, optional, intent(in) :: s(:) !! positive full lengths in axes order
        integer, optional, intent(in) :: axes(:) !! distinct one-based axes; default last size(s) or all
        character(*), optional, intent(in) :: norm !! backward (default), forward, or ortho
        logical, optional, intent(in) :: orthogonalize !! default true only for norm=ortho
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 setup failure
        real(dp), pointer :: input(:)
        real(dp), pointer :: output(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code

        if (present(status)) status = fft_status_invalid_input
        if (rank(y) < 1 .or. rank(y) > 15) return
        if (size(y,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(y),output,[size(y)])
        output = ieee_value(0.0_dp,ieee_quiet_nan)
        if (rank(x) /= rank(y)) return
        if (size(x,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(x),input,[size(x)])
        call nd_transform(cmplx(input,0.0_dp,dp),shape(x),shape(y),transformed, &
            7,s,axes,norm,code,type,orthogonalize)
        if (code == fft_status_success) output = real(transformed,dp)
        if (present(status)) status = code
    end subroutine dctn_real

    subroutine dctn_complex(x,y,type,s,axes,norm,orthogonalize,status)
        implicit none
        complex(dp), contiguous, target, intent(in) :: x(..) !! finite rank 1-15 input
        complex(dp), contiguous, target, intent(out) :: y(..) !! preallocated expected result shape
        integer, optional, intent(in) :: type !! type 1, 2 (default), 3, or 4
        integer, optional, intent(in) :: s(:) !! positive full lengths in axes order
        integer, optional, intent(in) :: axes(:) !! distinct one-based axes; default last size(s) or all
        character(*), optional, intent(in) :: norm !! backward (default), forward, or ortho
        logical, optional, intent(in) :: orthogonalize !! default true only for norm=ortho
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 setup failure
        complex(dp), pointer :: input(:)
        complex(dp), pointer :: output(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code

        if (present(status)) status = fft_status_invalid_input
        if (rank(y) < 1 .or. rank(y) > 15) return
        if (size(y,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(y),output,[size(y)])
        output = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
        if (rank(x) /= rank(y)) return
        if (size(x,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(x),input,[size(x)])
        call nd_transform(input,shape(x),shape(y),transformed, &
            7,s,axes,norm,code,type,orthogonalize)
        if (code == fft_status_success) output = transformed
        if (present(status)) status = code
    end subroutine dctn_complex

    subroutine idctn_real(x,y,type,s,axes,norm,orthogonalize,status)
        implicit none
        real(dp), contiguous, target, intent(in) :: x(..) !! finite rank 1-15 input
        real(dp), contiguous, target, intent(out) :: y(..) !! preallocated expected result shape
        integer, optional, intent(in) :: type !! type 1, 2 (default), 3, or 4
        integer, optional, intent(in) :: s(:) !! positive full lengths in axes order
        integer, optional, intent(in) :: axes(:) !! distinct one-based axes; default last size(s) or all
        character(*), optional, intent(in) :: norm !! backward (default), forward, or ortho
        logical, optional, intent(in) :: orthogonalize !! default true only for norm=ortho
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 setup failure
        real(dp), pointer :: input(:)
        real(dp), pointer :: output(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code

        if (present(status)) status = fft_status_invalid_input
        if (rank(y) < 1 .or. rank(y) > 15) return
        if (size(y,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(y),output,[size(y)])
        output = ieee_value(0.0_dp,ieee_quiet_nan)
        if (rank(x) /= rank(y)) return
        if (size(x,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(x),input,[size(x)])
        call nd_transform(cmplx(input,0.0_dp,dp),shape(x),shape(y),transformed, &
            8,s,axes,norm,code,type,orthogonalize)
        if (code == fft_status_success) output = real(transformed,dp)
        if (present(status)) status = code
    end subroutine idctn_real

    subroutine idctn_complex(x,y,type,s,axes,norm,orthogonalize,status)
        implicit none
        complex(dp), contiguous, target, intent(in) :: x(..) !! finite rank 1-15 input
        complex(dp), contiguous, target, intent(out) :: y(..) !! preallocated expected result shape
        integer, optional, intent(in) :: type !! type 1, 2 (default), 3, or 4
        integer, optional, intent(in) :: s(:) !! positive full lengths in axes order
        integer, optional, intent(in) :: axes(:) !! distinct one-based axes; default last size(s) or all
        character(*), optional, intent(in) :: norm !! backward (default), forward, or ortho
        logical, optional, intent(in) :: orthogonalize !! default true only for norm=ortho
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 setup failure
        complex(dp), pointer :: input(:)
        complex(dp), pointer :: output(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code

        if (present(status)) status = fft_status_invalid_input
        if (rank(y) < 1 .or. rank(y) > 15) return
        if (size(y,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(y),output,[size(y)])
        output = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
        if (rank(x) /= rank(y)) return
        if (size(x,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(x),input,[size(x)])
        call nd_transform(input,shape(x),shape(y),transformed, &
            8,s,axes,norm,code,type,orthogonalize)
        if (code == fft_status_success) output = transformed
        if (present(status)) status = code
    end subroutine idctn_complex

    subroutine dstn_real(x,y,type,s,axes,norm,orthogonalize,status)
        implicit none
        real(dp), contiguous, target, intent(in) :: x(..) !! finite rank 1-15 input
        real(dp), contiguous, target, intent(out) :: y(..) !! preallocated expected result shape
        integer, optional, intent(in) :: type !! type 1, 2 (default), 3, or 4
        integer, optional, intent(in) :: s(:) !! positive full lengths in axes order
        integer, optional, intent(in) :: axes(:) !! distinct one-based axes; default last size(s) or all
        character(*), optional, intent(in) :: norm !! backward (default), forward, or ortho
        logical, optional, intent(in) :: orthogonalize !! default true only for norm=ortho
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 setup failure
        real(dp), pointer :: input(:)
        real(dp), pointer :: output(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code

        if (present(status)) status = fft_status_invalid_input
        if (rank(y) < 1 .or. rank(y) > 15) return
        if (size(y,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(y),output,[size(y)])
        output = ieee_value(0.0_dp,ieee_quiet_nan)
        if (rank(x) /= rank(y)) return
        if (size(x,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(x),input,[size(x)])
        call nd_transform(cmplx(input,0.0_dp,dp),shape(x),shape(y),transformed, &
            9,s,axes,norm,code,type,orthogonalize)
        if (code == fft_status_success) output = real(transformed,dp)
        if (present(status)) status = code
    end subroutine dstn_real

    subroutine dstn_complex(x,y,type,s,axes,norm,orthogonalize,status)
        implicit none
        complex(dp), contiguous, target, intent(in) :: x(..) !! finite rank 1-15 input
        complex(dp), contiguous, target, intent(out) :: y(..) !! preallocated expected result shape
        integer, optional, intent(in) :: type !! type 1, 2 (default), 3, or 4
        integer, optional, intent(in) :: s(:) !! positive full lengths in axes order
        integer, optional, intent(in) :: axes(:) !! distinct one-based axes; default last size(s) or all
        character(*), optional, intent(in) :: norm !! backward (default), forward, or ortho
        logical, optional, intent(in) :: orthogonalize !! default true only for norm=ortho
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 setup failure
        complex(dp), pointer :: input(:)
        complex(dp), pointer :: output(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code

        if (present(status)) status = fft_status_invalid_input
        if (rank(y) < 1 .or. rank(y) > 15) return
        if (size(y,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(y),output,[size(y)])
        output = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
        if (rank(x) /= rank(y)) return
        if (size(x,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(x),input,[size(x)])
        call nd_transform(input,shape(x),shape(y),transformed, &
            9,s,axes,norm,code,type,orthogonalize)
        if (code == fft_status_success) output = transformed
        if (present(status)) status = code
    end subroutine dstn_complex

    subroutine idstn_real(x,y,type,s,axes,norm,orthogonalize,status)
        implicit none
        real(dp), contiguous, target, intent(in) :: x(..) !! finite rank 1-15 input
        real(dp), contiguous, target, intent(out) :: y(..) !! preallocated expected result shape
        integer, optional, intent(in) :: type !! type 1, 2 (default), 3, or 4
        integer, optional, intent(in) :: s(:) !! positive full lengths in axes order
        integer, optional, intent(in) :: axes(:) !! distinct one-based axes; default last size(s) or all
        character(*), optional, intent(in) :: norm !! backward (default), forward, or ortho
        logical, optional, intent(in) :: orthogonalize !! default true only for norm=ortho
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 setup failure
        real(dp), pointer :: input(:)
        real(dp), pointer :: output(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code

        if (present(status)) status = fft_status_invalid_input
        if (rank(y) < 1 .or. rank(y) > 15) return
        if (size(y,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(y),output,[size(y)])
        output = ieee_value(0.0_dp,ieee_quiet_nan)
        if (rank(x) /= rank(y)) return
        if (size(x,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(x),input,[size(x)])
        call nd_transform(cmplx(input,0.0_dp,dp),shape(x),shape(y),transformed, &
            10,s,axes,norm,code,type,orthogonalize)
        if (code == fft_status_success) output = real(transformed,dp)
        if (present(status)) status = code
    end subroutine idstn_real

    subroutine idstn_complex(x,y,type,s,axes,norm,orthogonalize,status)
        implicit none
        complex(dp), contiguous, target, intent(in) :: x(..) !! finite rank 1-15 input
        complex(dp), contiguous, target, intent(out) :: y(..) !! preallocated expected result shape
        integer, optional, intent(in) :: type !! type 1, 2 (default), 3, or 4
        integer, optional, intent(in) :: s(:) !! positive full lengths in axes order
        integer, optional, intent(in) :: axes(:) !! distinct one-based axes; default last size(s) or all
        character(*), optional, intent(in) :: norm !! backward (default), forward, or ortho
        logical, optional, intent(in) :: orthogonalize !! default true only for norm=ortho
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 setup failure
        complex(dp), pointer :: input(:)
        complex(dp), pointer :: output(:)
        complex(dp), allocatable :: transformed(:)
        integer :: code

        if (present(status)) status = fft_status_invalid_input
        if (rank(y) < 1 .or. rank(y) > 15) return
        if (size(y,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(y),output,[size(y)])
        output = cmplx(ieee_value(0.0_dp,ieee_quiet_nan),0.0_dp,dp)
        if (rank(x) /= rank(y)) return
        if (size(x,kind=int64) > int(huge(0),int64)) return
        call c_f_pointer(c_loc(x),input,[size(x)])
        call nd_transform(input,shape(x),shape(y),transformed, &
            10,s,axes,norm,code,type,orthogonalize)
        if (code == fft_status_success) output = transformed
        if (present(status)) status = code
    end subroutine idstn_complex


    subroutine nd_transform(x,xdims,ydims,y,op,s,axes,norm,code,type,orthogonalize)
        implicit none
        complex(dp), intent(in) :: x(:) !! flattened finite input in Fortran element order
        integer, intent(in) :: xdims(:) !! input extents
        integer, intent(in) :: ydims(:) !! caller-provided output extents
        complex(dp), allocatable, intent(out) :: y(:) !! flattened result on success
        integer, intent(in) :: op !! 1-6 Fourier pairs, 7-10 cosine/sine pairs
        integer, optional, intent(in) :: s(:) !! positive full transform lengths
        integer, optional, intent(in) :: axes(:) !! distinct one-based axes
        character(*), optional, intent(in) :: norm !! normalization name
        integer, intent(out) :: code !! completion status
        integer, optional, intent(in) :: type !! cosine/sine type, default 2
        logical, optional, intent(in) :: orthogonalize !! endpoint orthogonalization
        integer :: selected(15), lengths(15), dims(15), expected(15)
        integer :: r, count, i, j, axis, mode, kind, ierr, axis_op
        logical :: short_half, full_half

        code = fft_status_invalid_input
        r = size(xdims)
        if (r < 1 .or. r > 15 .or. size(ydims) /= r) return
        mode = normalization(norm)
        if (mode < 0 .or. .not. finite_complex(x)) return
        count = r
        if (present(s)) count = size(s)
        if (present(axes)) count = size(axes)
        if (count > r) return
        if (present(s)) then
            if (size(s) /= count) return
            if (any(s < 1)) return
        end if
        do i = 1, count
            selected(i) = r-count+i
            if (present(axes)) selected(i) = axes(i)
            axis = selected(i)
            if (axis < 1 .or. axis > r) return
            if (any(selected(:i-1) == axis)) return
            lengths(i) = xdims(axis)
            if (present(s)) lengths(i) = s(i)
        end do
        short_half = op == 3 .or. op == 6
        full_half = op == 4 .or. op == 5
        if (short_half .or. full_half) then
            if (count == 0) return
        end if
        if (full_half .and. .not. present(s)) then
            axis = selected(count)
            if (xdims(axis) < 2 .or. xdims(axis)-1 > huge(0)/2) return
            lengths(count) = 2*(xdims(axis)-1)
        end if
        if (any(lengths(:count) < 1)) return
        kind = 2
        if (present(type)) kind = type
        if (op >= 7) then
            if (kind < 1 .or. kind > 4) return
            if (any(lengths(:count) > huge(0)/8)) return
            if (kind == 1 .and. (op == 7 .or. op == 8)) then
                if (any(lengths(:count) < 2)) return
            end if
        end if
        expected(:r) = xdims
        do i = 1, count
            expected(selected(i)) = lengths(i)
        end do
        if (short_half) expected(selected(count)) = lengths(count)/2+1
        if (any(expected(:r) /= ydims)) return
        ! Limit all intermediate arrays as well as input/output to default-integer indexing.
        dims(:r) = xdims
        do i = 1, count
            axis = selected(i)
            dims(axis) = max(dims(axis),expected(axis))
        end do
        if (.not. valid_nd_size(dims(:r))) return
        code = fft_status_setup_failure
        allocate(y(size(x)),stat=ierr)
        if (ierr /= 0) return
        y = x
        dims(:r) = xdims
        ! Real-to-half transforms must precede the complex axes. Half-to-real
        ! transforms must follow them: only then are their endpoints real.
        if (short_half) then
            call nd_axis(y,dims(:r),selected(count),lengths(count),op, &
                mode,norm,kind,orthogonalize,code)
            if (code /= fft_status_success) return
        end if
        do j = 1, count
            i = j
            if (short_half .and. i == count) cycle
            if (full_half .and. i == count) cycle
            axis_op = op
            if (op == 3 .or. op == 5) axis_op = 1
            if (op == 4 .or. op == 6) axis_op = 2
            call nd_axis(y,dims(:r),selected(i),lengths(i),axis_op, &
                mode,norm,kind,orthogonalize,code)
            if (code /= fft_status_success) return
        end do
        if (full_half) then
            call nd_axis(y,dims(:r),selected(count),lengths(count),op, &
                mode,norm,kind,orthogonalize,code)
            if (code /= fft_status_success) return
        end if
        code = fft_status_success
    end subroutine nd_transform

    pure logical function valid_nd_size(dims) result(valid)
        implicit none
        integer, intent(in) :: dims(:) !! nonnegative extents
        integer :: total, i
        valid = .false.
        total = 1
        do i = 1, size(dims)
            if (dims(i) < 0) return
            if (dims(i) == 0) then
                total = 0
            else
                if (total > huge(0)/dims(i)) return
                total = total*dims(i)
            end if
        end do
        valid = .true.
    end function valid_nd_size

    subroutine nd_axis(values,dims,axis,n,op,mode,norm,kind,orthogonalize,code)
        implicit none
        complex(dp), allocatable, intent(inout) :: values(:) !! flattened input and result
        integer, intent(inout) :: dims(:) !! extents updated after this transform
        integer, intent(in) :: axis !! one-based transform dimension
        integer, intent(in) :: n !! positive full transform length
        integer, intent(in) :: op !! transform selector
        integer, intent(in) :: mode !! validated normalization code
        character(*), optional, intent(in) :: norm !! normalization name
        integer, intent(in) :: kind !! validated cosine/sine type
        logical, optional, intent(in) :: orthogonalize !! endpoint orthogonalization
        integer, intent(out) :: code !! completion status
        complex(dp), allocatable :: output(:), line(:), transformed(:)
        type(fft_plan) :: plan
        integer :: stride, blocks, oldn, newn, copied, block, offset, j, base, dest
        integer :: ierr, setup
        real(dp) :: factor
        logical :: inverse, short_half, full_half

        code = fft_status_setup_failure
        oldn = dims(axis)
        newn = n
        short_half = op == 3 .or. op == 6
        full_half = op == 4 .or. op == 5
        if (short_half) newn = n/2+1
        ! An unselected zero extent is an empty batch. Avoid stride products
        ! and plan setup altogether, even when other extents are very large.
        if (any(dims(:axis-1) == 0) .or. any(dims(axis+1:) == 0)) then
            dims(axis) = newn
            allocate(output(0),stat=ierr)
            if (ierr /= 0) return
            call move_alloc(output,values)
            code = fft_status_success
            return
        end if
        stride = product(dims(:axis-1))
        blocks = product(dims(axis+1:))
        dims(axis) = newn
        allocate(output(product(dims)),line(n),stat=ierr)
        if (ierr /= 0) return
        inverse = op == 2 .or. op == 4 .or. op == 6 .or. op == 8 .or. op == 10
        if (op <= 6) then
            call plan%initialize(n,setup)
            if (setup /= 0) return
        end if
        factor = transform_scale(real(n,dp),mode,inverse)
        copied = min(oldn,n)
        if (full_half) copied = min(oldn,n/2+1)
        do block = 0, blocks-1
            do offset = 1, stride
                line = (0.0_dp,0.0_dp)
                base = block*stride*oldn+offset
                do j = 1, copied
                    line(j) = values(base+(j-1)*stride)
                end do
                if (full_half) then
                    line(1) = cmplx(real(line(1),dp),0.0_dp,dp)
                    if (mod(n,2) == 0) line(n/2+1) = cmplx(real(line(n/2+1),dp),0.0_dp,dp)
                    do j = 2, n/2+mod(n,2)
                        line(n-j+2) = conjg(line(j))
                    end do
                end if
                if (op <= 6) then
                    line = line*factor
                    if (inverse) line = conjg(line)
                    call plan%execute(line,setup)
                    if (setup /= 0) return
                    if (inverse) line = conjg(line)
                    if (full_half) line = cmplx(real(line,dp),0.0_dp,dp)
                else
                    transformed = r2r_transform(line,kind,n,norm,orthogonalize,op >= 9,inverse,code)
                    if (code /= fft_status_success) return
                    line = transformed
                end if
                dest = block*stride*newn+offset
                do j = 1, newn
                    output(dest+(j-1)*stride) = line(j)
                end do
            end do
        end do
        call move_alloc(output,values)
        code = fft_status_success
    end subroutine nd_axis

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
