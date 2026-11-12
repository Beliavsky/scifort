! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2021 Jacob Williams
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
!
! Scalar minimization for the scipy.optimize translation.
!
! The bounded Brent kernel is adapted from fmin 1.1.2 by Jacob Williams
! (BSD-3-Clause) and reconciled with SciPy 1.17.0 bounded minimization
! semantics. The downhill bracketing, unbounded Brent, golden-section,
! dispatch, result, and evaluation-count semantics are adapted from SciPy
! 1.17.0 scipy.optimize._optimize and scipy.optimize._minimize (BSD-3-Clause).
! The scalar bracket_minimum expansion, finite-limit interpolation, and
! status accounting are adapted from SciPy 1.17.0 scipy.optimize._bracket.
! Scalar Chandrupatla minimization is adapted from SciPy 1.17.0
! scipy.optimize._chandrupatla (BSD-3-Clause). The rank-one
! find_minimum_vectorized active-lane controller is original SciFort code.

module scifort_optimize_minimize_scalar
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_negative_inf, ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_optimize_elementwise_batch, only : elementwise_broadcast, &
        elementwise_broadcast_optional, elementwise_flatten, &
        elementwise_broadcast_arguments, &
        elementwise_argument, elementwise_broadcast_argument_list, &
        elementwise_typed_argument, elementwise_broadcast_typed_argument_list, &
        elementwise_pack_typed_arguments, elementwise_gather_typed_indices
    implicit none
    private

    integer, parameter, public :: minimize_status_success = 0
    integer, parameter, public :: minimize_status_bracket_error = -1
    integer, parameter, public :: minimize_status_convergence_error = -2
    integer, parameter, public :: minimize_status_value_error = -3
    integer, parameter, public :: minimize_status_invalid_input = -5
    integer, parameter, public :: bracket_minimum_status_limit = -1
    integer, parameter, public :: find_minimum_status_callback = -4

    real(dp), parameter :: default_brent_xtol = 1.48e-8_dp
    real(dp), parameter :: default_bounded_xatol = 1.0e-5_dp
    real(dp), parameter :: default_golden_xtol = sqrt(epsilon(1.0_dp))
    integer, parameter :: default_brent_maxiter = 500
    integer, parameter :: default_bounded_maxiter = 500
    integer, parameter :: default_golden_maxiter = 5000
    integer, parameter :: default_bracket_maxiter = 1000

    type, public :: minimize_result
        real(dp) :: x = 0.0_dp !! approximate minimizer
        real(dp) :: fun = 0.0_dp !! objective value at x
        integer :: iterations = 0 !! solver iterations performed
        integer :: function_calls = 0 !! objective evaluations performed
        integer :: status = minimize_status_invalid_input !! termination code
        logical :: success = .false. !! true when the stopping rule was met
        character(len=48) :: flag = 'invalid input' !! human-readable termination flag
        character(len=16) :: method = '' !! solver method actually used
    end type minimize_result

    type :: bracket_state
        real(dp) :: xa = 0.0_dp
        real(dp) :: xb = 0.0_dp
        real(dp) :: xc = 0.0_dp
        real(dp) :: fa = 0.0_dp
        real(dp) :: fb = 0.0_dp
        real(dp) :: fc = 0.0_dp
        integer :: function_calls = 0
        integer :: iterations = 0
        logical :: valid = .false.
        integer :: status = minimize_status_bracket_error
    end type bracket_state

    type, public :: bracket_result
        real(dp) :: xa = 0.0_dp !! first point in the downhill bracket
        real(dp) :: xb = 0.0_dp !! middle point with the lowest sampled objective
        real(dp) :: xc = 0.0_dp !! third point in the downhill bracket
        real(dp) :: fa = 0.0_dp !! objective at xa
        real(dp) :: fb = 0.0_dp !! objective at xb
        real(dp) :: fc = 0.0_dp !! objective at xc
        integer :: function_calls = 0 !! number of objective evaluations
        integer :: iterations = 0 !! iterations of bracket expansion
        integer :: status = minimize_status_invalid_input !! zero for a valid bracket
        logical :: success = .false. !! true when the three points form a valid bracket
    end type bracket_result

    type, public :: bracket_minimum_result
        real(dp) :: xl = 0.0_dp !! leftmost bracket endpoint
        real(dp) :: xm = 0.0_dp !! middle bracket point
        real(dp) :: xr = 0.0_dp !! rightmost bracket endpoint
        real(dp) :: fl = 0.0_dp !! objective value at xl
        real(dp) :: fm = 0.0_dp !! objective value at xm
        real(dp) :: fr = 0.0_dp !! objective value at xr
        integer :: nit = 0 !! number of bracket expansion iterations
        integer :: nfev = 0 !! number of scalar objective evaluations
        integer :: status = minimize_status_invalid_input !! SciPy-style termination code
        logical :: success = .false. !! true if the triple brackets a minimum
    end type bracket_minimum_result

    type, public :: find_minimum_result
        real(dp) :: x = 0.0_dp !! current best abscissa
        real(dp) :: f_x = 0.0_dp !! value at current best abscissa
        real(dp) :: xl = 0.0_dp !! left bracket endpoint
        real(dp) :: xm = 0.0_dp !! middle bracket point
        real(dp) :: xr = 0.0_dp !! right bracket endpoint
        real(dp) :: fl = 0.0_dp !! function at left endpoint
        real(dp) :: fm = 0.0_dp !! function at middle point
        real(dp) :: fr = 0.0_dp !! function at right endpoint
        integer :: nit = 0 !! number of iterations
        integer :: nfev = 0 !! total scalar function evaluations
        integer :: status = minimize_status_invalid_input !! 0 success, -1/-2/-3/-4 failure
        logical :: success = .false. !! convergence to requested tolerances
    end type find_minimum_result

    ! State for one downhill search in vectorized minimum bracketing.
    type :: vector_bracket_minimum_lane
        real(dp) :: xl = 0.0_dp
        real(dp) :: xm = 0.0_dp
        real(dp) :: xr = 0.0_dp
        real(dp) :: fl = 0.0_dp
        real(dp) :: fm = 0.0_dp
        real(dp) :: fr = 0.0_dp
        real(dp) :: limit = 0.0_dp
        real(dp) :: origin = 0.0_dp
        real(dp) :: step = 0.0_dp
        logical :: limited = .false.
    end type vector_bracket_minimum_lane

    abstract interface
        function scalar_objective_function(x) result(fx)
            import dp
            implicit none
            real(dp), intent(in) :: x !! scalar abscissa
            real(dp) :: fx !! scalar objective value
        end function scalar_objective_function

        subroutine elementwise_minimum_vector_function(x, fx)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! packed objective evaluation points
            real(dp), intent(out) :: fx(:) !! objective values, same length as x
        end subroutine elementwise_minimum_vector_function
        subroutine elementwise_minimum_args_function(x, args, fx)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! packed active objective points
            real(dp), intent(in) :: args(:,:) !! packed parameters (nargs, size(x))
            real(dp), intent(out) :: fx(:) !! values for each packed point
        end subroutine elementwise_minimum_args_function
        subroutine elementwise_minimum_typed_args_function(x, args, fx)
            import dp, elementwise_typed_argument
            implicit none
            real(dp), intent(in) :: x(:) !! packed active abscissae
            type(elementwise_typed_argument), intent(in) :: args(:) !! typed values, per active lane
            real(dp), intent(out) :: fx(:) !! returned objective values
        end subroutine elementwise_minimum_typed_args_function
        subroutine elementwise_minimum_indexed_args_function(x, args, indices, fx)
            import dp, elementwise_typed_argument
            implicit none
            real(dp), intent(in) :: x(:) !! active abscissae in packed callback order
            type(elementwise_typed_argument), intent(in) :: args(:) !! original typed input snapshots
            integer, intent(in) :: indices(:,:) !! (nargs,nactive) offsets into args(i)%values
            real(dp), intent(out) :: fx(:) !! objective at each active abscissa
        end subroutine elementwise_minimum_indexed_args_function

        subroutine find_minimum_callback(state, stop)
            import find_minimum_result
            implicit none
            type(find_minimum_result), intent(in) :: state !! current minimizer and bracket
            logical, intent(out) :: stop !! true to stop iteration
        end subroutine find_minimum_callback
    end interface

    public :: scalar_objective_function
    public :: minimize_scalar
    public :: bracket
    public :: bracket_minimum
    public :: find_minimum
    public :: find_minimum_vectorized
    public :: find_minimum_vectorized_args_into
    public :: find_minimum_vectorized_indexed_args_into
    public :: bracket_minimum_vectorized_args_into
    public :: bracket_minimum_vectorized_indexed_args_into
    public :: elementwise_minimum_args_function
    public :: elementwise_minimum_typed_args_function
    public :: elementwise_minimum_indexed_args_function
    public :: find_minimum_vectorized_into
    public :: bracket_minimum_vectorized
    public :: bracket_minimum_vectorized_into
    public :: elementwise_minimum_vector_function
    public :: find_minimum_callback
    public :: brent
    public :: bounded
    public :: golden

    ! Packed vector callbacks accept independent grids of rank one through four.
    ! The objective receives active points, not the entire problem grid.
    interface find_minimum_vectorized_args_into
        module procedure find_minimum_vectorized_args_into_stacked
        module procedure find_minimum_vectorized_args_into_list
        module procedure find_minimum_vectorized_args_into_typed
    end interface find_minimum_vectorized_args_into

    interface bracket_minimum_vectorized_args_into
        module procedure bracket_minimum_vectorized_args_into_stacked
        module procedure bracket_minimum_vectorized_args_into_list
        module procedure bracket_minimum_vectorized_args_into_typed
    end interface bracket_minimum_vectorized_args_into

    interface find_minimum_vectorized
        module procedure find_minimum_vectorized_rank1
        module procedure find_minimum_vectorized_rank2
        module procedure find_minimum_vectorized_rank3
        module procedure find_minimum_vectorized_rank4
    end interface find_minimum_vectorized

    interface bracket_minimum_vectorized
        module procedure bracket_minimum_vectorized_rank1
        module procedure bracket_minimum_vectorized_rank2
        module procedure bracket_minimum_vectorized_rank3
        module procedure bracket_minimum_vectorized_rank4
    end interface bracket_minimum_vectorized

    interface find_minimum
        module procedure find_minimum_scalar
        module procedure find_minimum_vector
        module procedure find_minimum_matrix
        module procedure find_minimum_tensor
        module procedure find_minimum_rank14
        module procedure find_minimum_rank13
        module procedure find_minimum_rank12
        module procedure find_minimum_rank11
        module procedure find_minimum_rank10
        module procedure find_minimum_rank9
        module procedure find_minimum_rank8
        module procedure find_minimum_rank7
        module procedure find_minimum_rank6
        module procedure find_minimum_rank5
        module procedure find_minimum_rank4
    end interface find_minimum

    interface bracket_minimum
        module procedure bracket_minimum_scalar
        module procedure bracket_minimum_vector
        module procedure bracket_minimum_matrix
        module procedure bracket_minimum_tensor
        module procedure bracket_minimum_rank15
        module procedure bracket_minimum_rank14
        module procedure bracket_minimum_rank13
        module procedure bracket_minimum_rank12
        module procedure bracket_minimum_rank11
        module procedure bracket_minimum_rank10
        module procedure bracket_minimum_rank9
        module procedure bracket_minimum_rank8
        module procedure bracket_minimum_rank7
        module procedure bracket_minimum_rank6
        module procedure bracket_minimum_rank5
        module procedure bracket_minimum_rank4
    end interface bracket_minimum

contains

    subroutine find_minimum_vectorized_args_into_stacked(f, bracket, args, res, &
            xatol, xrtol, fatol, frtol, maxiter, callback)
        ! Explicit args entry: retain the legacy (x, fx) callback ABI.
        procedure(elementwise_minimum_args_function) :: f !! receives x, packed args, fx
        real(dp), intent(in) :: bracket(..) !! starting problem grid
        real(dp), intent(in) :: args(..) !! stacked parameters, first axis = nargs
        type(find_minimum_result), intent(inout) :: res(..) !! initialized grid results
        real(dp), intent(in), optional :: xatol(..) !! broadcast tolerance
        real(dp), intent(in), optional :: xrtol(..) !! broadcast tolerance
        real(dp), intent(in), optional :: fatol(..) !! broadcast tolerance
        real(dp), intent(in), optional :: frtol(..) !! broadcast tolerance
        integer, intent(in), optional :: maxiter !! iteration limit
        procedure(find_minimum_callback), optional :: callback !! lane stop callback
        call find_minimum_vectorized_into(bracket=bracket, res=res, &
            xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, &
            callback=callback, args=args, f_args=f)
    end subroutine find_minimum_vectorized_args_into_stacked

    subroutine bracket_minimum_vectorized_args_into_stacked(f, xm0, args, res, &
            xl0, xr0, xmin, xmax, factor, maxiter)
        ! Explicit args entry: retain the legacy (x, fx) callback ABI.
        procedure(elementwise_minimum_args_function) :: f !! receives x, packed args, fx
        real(dp), intent(in) :: xm0(..) !! starting problem grid
        real(dp), intent(in) :: args(..) !! stacked parameters, first axis = nargs
        type(bracket_minimum_result), intent(inout) :: res(..) !! initialized grid results
        real(dp), intent(in), optional :: xl0(..) !! broadcast left endpoint
        real(dp), intent(in), optional :: xr0(..) !! broadcast right endpoint
        real(dp), intent(in), optional :: xmin(..) !! broadcast lower limit
        real(dp), intent(in), optional :: xmax(..) !! broadcast upper limit
        real(dp), intent(in), optional :: factor(..) !! broadcast expansion factor
        integer, intent(in), optional :: maxiter !! iteration limit
        call bracket_minimum_vectorized_into(xm0=xm0, res=res, &
            xl0=xl0, xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, &
            maxiter=maxiter, args=args, f_args=f)
    end subroutine bracket_minimum_vectorized_args_into_stacked

    subroutine find_minimum_vectorized_args_into_list(f, bracket, args, res, &
            xatol, xrtol, fatol, frtol, maxiter, callback)
        ! Alternative to a stacked args array: every parameter may have its
        ! own scalar/array rank and shape.  The callback ABI is unchanged.
        procedure(elementwise_minimum_args_function) :: f !! receives x, packed args, fx
        real(dp), intent(in) :: bracket(..) !! starting problem grid
        type(elementwise_argument), intent(in) :: args(:) !! independent parameters
        type(find_minimum_result), intent(inout) :: res(..) !! result grid
        real(dp), intent(in), optional :: xatol(..) !! broadcast tolerance
        real(dp), intent(in), optional :: xrtol(..) !! broadcast tolerance
        real(dp), intent(in), optional :: fatol(..) !! broadcast tolerance
        real(dp), intent(in), optional :: frtol(..) !! broadcast tolerance
        integer, intent(in), optional :: maxiter !! iteration limit
        procedure(find_minimum_callback), optional :: callback !! lane stop callback
        call find_minimum_vectorized_into(bracket=bracket, res=res, &
            xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, &
            callback=callback, args_list=args, f_args=f)
    end subroutine find_minimum_vectorized_args_into_list

    subroutine bracket_minimum_vectorized_args_into_list(f, xm0, args, res, &
            xl0, xr0, xmin, xmax, factor, maxiter)
        procedure(elementwise_minimum_args_function) :: f !! receives x, packed args, fx
        real(dp), intent(in) :: xm0(..) !! middle starting grid
        type(bracket_minimum_result), intent(inout) :: res(..) !! result grid
        type(elementwise_argument), intent(in) :: args(:) !! independent parameters
        real(dp), intent(in), optional :: xl0(..) !! left starting points
        real(dp), intent(in), optional :: xr0(..) !! right starting points
        real(dp), intent(in), optional :: xmin(..) !! lower bounds
        real(dp), intent(in), optional :: xmax(..) !! upper bounds
        real(dp), intent(in), optional :: factor(..) !! expansion factor
        integer, intent(in), optional :: maxiter !! iteration limit
        call bracket_minimum_vectorized_into(xm0=xm0, res=res, &
            xl0=xl0, xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, &
            maxiter=maxiter, args_list=args, f_args=f)
    end subroutine bracket_minimum_vectorized_args_into_list



    subroutine find_minimum_vectorized_args_into_typed(f, bracket, args, res, &
            xatol, xrtol, fatol, frtol, maxiter, callback)
        ! No conversion to real(dp): forward independently broadcast typed parameters.
        procedure(elementwise_minimum_typed_args_function) :: f !! receives typed, packed active parameters
        real(dp), intent(in) :: bracket(..) !! problem-grid inputs
        type(elementwise_typed_argument), intent(in) :: args(:) !! typed argument list
        type(find_minimum_result), intent(inout) :: res(..) !! grid-shaped output results
        real(dp), intent(in), optional :: xatol(..) !! broadcast absolute x tolerance
        real(dp), intent(in), optional :: xrtol(..) !! broadcast relative x tolerance
        real(dp), intent(in), optional :: fatol(..) !! broadcast absolute function tolerance
        real(dp), intent(in), optional :: frtol(..) !! broadcast relative function tolerance
        integer, intent(in), optional :: maxiter !! iteration limit
        procedure(find_minimum_callback), optional :: callback !! per-lane callback
        call find_minimum_vectorized_into(bracket=bracket, res=res, &
            xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, &
            maxiter=maxiter, callback=callback, typed_args=args, f_typed=f)
    end subroutine find_minimum_vectorized_args_into_typed

    subroutine bracket_minimum_vectorized_args_into_typed(f, xm0, args, res, &
            xl0, xr0, xmin, xmax, factor, maxiter)
        ! No conversion to real(dp): forward independently broadcast typed parameters.
        procedure(elementwise_minimum_typed_args_function) :: f !! receives typed, packed active parameters
        real(dp), intent(in) :: xm0(..) !! problem-grid inputs
        type(elementwise_typed_argument), intent(in) :: args(:) !! typed argument list
        type(bracket_minimum_result), intent(inout) :: res(..) !! grid-shaped output results
        real(dp), intent(in), optional :: xl0(..) !! optional left endpoint
        real(dp), intent(in), optional :: xr0(..) !! optional right endpoint
        real(dp), intent(in), optional :: xmin(..) !! broadcast lower bound
        real(dp), intent(in), optional :: xmax(..) !! broadcast upper bound
        real(dp), intent(in), optional :: factor(..) !! broadcast expansion factor
        integer, intent(in), optional :: maxiter !! iteration limit
        call bracket_minimum_vectorized_into(xm0=xm0, res=res, &
            xl0=xl0, xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, &
            maxiter=maxiter, typed_args=args, f_typed=f)
    end subroutine bracket_minimum_vectorized_args_into_typed

    subroutine find_minimum_vectorized_indexed_args_into(f, bracket, args, res, &
            xatol, xrtol, fatol, frtol, maxiter, callback)
        ! Zero-copy typed callback: source boxes are never replicated per evaluation.
        procedure(elementwise_minimum_indexed_args_function) :: f !! receives source values and active-lane offsets
        real(dp), intent(in) :: bracket(..) !! initial three-point brackets
        type(elementwise_typed_argument), intent(in) :: args(:) !! mixed-type argument list
        type(find_minimum_result), intent(inout) :: res(..) !! grid of minimizers
        real(dp), intent(in), optional :: xatol(..) !! absolute abscissa tolerance
        real(dp), intent(in), optional :: xrtol(..) !! relative abscissa tolerance
        real(dp), intent(in), optional :: fatol(..) !! absolute objective tolerance
        real(dp), intent(in), optional :: frtol(..) !! relative objective tolerance
        integer, intent(in), optional :: maxiter !! iteration limit
        procedure(find_minimum_callback), optional :: callback !! per-lane stopping callback
        call find_minimum_vectorized_into(bracket=bracket, res=res, &
            xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, &
            maxiter=maxiter, callback=callback, typed_args=args, f_indexed=f)
    end subroutine find_minimum_vectorized_indexed_args_into

    subroutine bracket_minimum_vectorized_indexed_args_into(f, xm0, args, res, &
            xl0, xr0, xmin, xmax, factor, maxiter)
        procedure(elementwise_minimum_indexed_args_function) :: f !! receives source values and active-lane offsets
        real(dp), intent(in) :: xm0(..) !! initial middle points
        type(elementwise_typed_argument), intent(in) :: args(:) !! mixed-type argument list
        type(bracket_minimum_result), intent(inout) :: res(..) !! grid of minimum brackets
        real(dp), intent(in), optional :: xl0(..) !! initial left points
        real(dp), intent(in), optional :: xr0(..) !! initial right points
        real(dp), intent(in), optional :: xmin(..) !! lower search limits
        real(dp), intent(in), optional :: xmax(..) !! upper search limits
        real(dp), intent(in), optional :: factor(..) !! expansion factors
        integer, intent(in), optional :: maxiter !! iteration limit
        call bracket_minimum_vectorized_into(xm0=xm0, res=res, &
            xl0=xl0, xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, &
            maxiter=maxiter, typed_args=args, f_indexed=f)
    end subroutine bracket_minimum_vectorized_indexed_args_into

    subroutine find_minimum_vectorized_into(f, bracket, res, xatol, xrtol, fatol, frtol, &
            maxiter, &
            callback, args, f_args, args_list, typed_args, f_typed, f_indexed)
        ! Assumed-rank grid dispatcher. Result rank and shape must match input.
        ! Only this subroutine broadcasts and reshapes; the rank-one solver
        ! remains responsible for active-lane callback ordering and statuses.
        procedure(elementwise_minimum_vector_function), optional :: f !! packed vector objective
        real(dp), intent(in), optional :: args(..) !! leading-axis objective parameters
        type(elementwise_argument), intent(in), optional :: args_list(:) !! independently shaped args
        procedure(elementwise_minimum_args_function), optional :: f_args !! packed parameterized objective
        type(elementwise_typed_argument), intent(in), optional :: typed_args(:) !! typed arguments
        procedure(elementwise_minimum_typed_args_function), optional :: f_typed !! typed packed callback
        procedure(elementwise_minimum_indexed_args_function), optional :: f_indexed !! zero-copy typed callback
        real(dp), intent(in) :: bracket(..) !! grid inputs with a leading bracket axis
        type(find_minimum_result), intent(inout) :: res(..) !! results, same problem-grid shape
        real(dp), intent(in), optional :: xatol(..) !! broadcast tolerance argument
        real(dp), intent(in), optional :: xrtol(..) !! broadcast tolerance argument
        real(dp), intent(in), optional :: fatol(..) !! broadcast tolerance argument
        real(dp), intent(in), optional :: frtol(..) !! broadcast tolerance argument
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_minimum_callback), optional :: callback !! scalar per-lane stop callback
        type(find_minimum_result), allocatable :: flat(:)
        real(dp), allocatable :: linear(:)
        real(dp), allocatable :: pairs(:,:)
        real(dp), allocatable :: ax(:), rx(:), af(:), rf(:)
        integer, allocatable :: dims(:), grid(:)
        real(dp), allocatable :: packed_args(:,:)
        integer, allocatable :: typed_indices(:,:)
        logical :: valid

        ! Initialize even when shapes or broadcast parameters are invalid.
        select rank(res)
        rank(0)
            res = find_minimum_result()
        rank(1)
            res = find_minimum_result()
        rank(2)
            res = find_minimum_result()
        rank(3)
            res = find_minimum_result()
        rank(4)
            res = find_minimum_result()
        rank(5)
            res = find_minimum_result()
        rank(6)
            res = find_minimum_result()
        rank(7)
            res = find_minimum_result()
        rank(8)
            res = find_minimum_result()
        rank(9)
            res = find_minimum_result()
        rank(10)
            res = find_minimum_result()
        rank(11)
            res = find_minimum_result()
        rank(12)
            res = find_minimum_result()
        rank(13)
            res = find_minimum_result()
        rank(14)
            res = find_minimum_result()
        rank(15)
            res = find_minimum_result()
        rank default
            return
        end select

        if (present(args) .and. present(args_list)) return
        if (present(f_args) .neqv. (present(args) .or. present(args_list))) return
        if (present(f_typed) .and. present(f_indexed)) return
        if ((present(f_typed) .or. present(f_indexed)) .neqv. present(typed_args)) return
        if (present(typed_args) .and. &
                (present(f) .or. present(f_args) .or. present(args) .or. present(args_list))) return
        if (.not. present(f) .and. .not. present(f_args) .and. &
                .not. present(f_typed) .and. .not. present(f_indexed)) return
        dims = shape(bracket)
        grid = shape(res)
        if (size(dims) /= size(grid) + 1) return
        if (dims(1) /= 3) return
        if (any(dims(2:) /= grid)) return
        if (size(res) == 0) return
        valid = .true.
        call elementwise_broadcast_optional(xatol, grid, ax, valid)
        call elementwise_broadcast_optional(xrtol, grid, rx, valid)
        call elementwise_broadcast_optional(fatol, grid, af, valid)
        call elementwise_broadcast_optional(frtol, grid, rf, valid)
        if (present(args)) call elementwise_broadcast_arguments(args, grid, packed_args, valid)
        if (present(args_list)) call elementwise_broadcast_argument_list(args_list, grid, packed_args, valid)
        if (present(typed_args)) then
            call elementwise_broadcast_typed_argument_list(typed_args, grid, typed_indices, valid)
        end if
        if (.not. valid) return
        call elementwise_flatten(bracket, linear)
        pairs = reshape(linear, [3,size(res)])
        flat = find_minimum_vectorized_rank1(f, pairs, &
            xatol=ax, xrtol=rx, fatol=af, frtol=rf, maxiter=maxiter, callback=callback, &
            args_data=packed_args, f_args=f_args, &
            typed_args=typed_args, typed_indices=typed_indices, f_typed=f_typed, f_indexed=f_indexed)
        ! SELECT RANK is required because Fortran 2018 has no assumed-rank
        ! function result or rank-independent RESHAPE assignment.
        select rank(res)
        rank(0)
            res = flat(1)
        rank(1)
            res = reshape(flat, shape(res))
        rank(2)
            res = reshape(flat, shape(res))
        rank(3)
            res = reshape(flat, shape(res))
        rank(4)
            res = reshape(flat, shape(res))
        rank(5)
            res = reshape(flat, shape(res))
        rank(6)
            res = reshape(flat, shape(res))
        rank(7)
            res = reshape(flat, shape(res))
        rank(8)
            res = reshape(flat, shape(res))
        rank(9)
            res = reshape(flat, shape(res))
        rank(10)
            res = reshape(flat, shape(res))
        rank(11)
            res = reshape(flat, shape(res))
        rank(12)
            res = reshape(flat, shape(res))
        rank(13)
            res = reshape(flat, shape(res))
        rank(14)
            res = reshape(flat, shape(res))
        rank(15)
            res = reshape(flat, shape(res))
        rank default
            return
        end select
    end subroutine find_minimum_vectorized_into

    subroutine bracket_minimum_vectorized_into(f, xm0, res, xl0, xr0, xmin, xmax, factor, &
            maxiter, args, f_args, args_list, typed_args, f_typed, f_indexed)
        ! Assumed-rank grid dispatcher. Result rank and shape must match input.
        ! Only this subroutine broadcasts and reshapes; the rank-one solver
        ! remains responsible for active-lane callback ordering and statuses.
        procedure(elementwise_minimum_vector_function), optional :: f !! packed vector objective
        real(dp), intent(in), optional :: args(..) !! leading-axis objective parameters
        type(elementwise_argument), intent(in), optional :: args_list(:) !! independently shaped args
        procedure(elementwise_minimum_args_function), optional :: f_args !! packed parameterized objective
        type(elementwise_typed_argument), intent(in), optional :: typed_args(:) !! typed arguments
        procedure(elementwise_minimum_typed_args_function), optional :: f_typed !! typed packed callback
        procedure(elementwise_minimum_indexed_args_function), optional :: f_indexed !! zero-copy typed callback
        real(dp), intent(in) :: xm0(..) !! grid inputs with arbitrary grid rank
        type(bracket_minimum_result), intent(inout) :: res(..) !! results, same problem-grid shape
        real(dp), intent(in), optional :: xl0(..) !! broadcast bracketing argument
        real(dp), intent(in), optional :: xr0(..) !! broadcast bracketing argument
        real(dp), intent(in), optional :: xmin(..) !! broadcast bracketing argument
        real(dp), intent(in), optional :: xmax(..) !! broadcast bracketing argument
        real(dp), intent(in), optional :: factor(..) !! broadcast bracketing argument
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_minimum_result), allocatable :: flat(:)
        real(dp), allocatable :: linear(:)
        real(dp), allocatable :: left(:), right(:), lo(:), hi(:), grow(:)
        integer, allocatable :: dims(:), grid(:)
        real(dp), allocatable :: packed_args(:,:)
        integer, allocatable :: typed_indices(:,:)
        logical :: valid

        ! Initialize even when shapes or broadcast parameters are invalid.
        select rank(res)
        rank(0)
            res = bracket_minimum_result()
        rank(1)
            res = bracket_minimum_result()
        rank(2)
            res = bracket_minimum_result()
        rank(3)
            res = bracket_minimum_result()
        rank(4)
            res = bracket_minimum_result()
        rank(5)
            res = bracket_minimum_result()
        rank(6)
            res = bracket_minimum_result()
        rank(7)
            res = bracket_minimum_result()
        rank(8)
            res = bracket_minimum_result()
        rank(9)
            res = bracket_minimum_result()
        rank(10)
            res = bracket_minimum_result()
        rank(11)
            res = bracket_minimum_result()
        rank(12)
            res = bracket_minimum_result()
        rank(13)
            res = bracket_minimum_result()
        rank(14)
            res = bracket_minimum_result()
        rank(15)
            res = bracket_minimum_result()
        rank default
            return
        end select

        if (present(args) .and. present(args_list)) return
        if (present(f_args) .neqv. (present(args) .or. present(args_list))) return
        if (present(f_typed) .and. present(f_indexed)) return
        if ((present(f_typed) .or. present(f_indexed)) .neqv. present(typed_args)) return
        if (present(typed_args) .and. &
                (present(f) .or. present(f_args) .or. present(args) .or. present(args_list))) return
        if (.not. present(f) .and. .not. present(f_args) .and. &
                .not. present(f_typed) .and. .not. present(f_indexed)) return
        dims = shape(xm0)
        grid = shape(res)
        if (size(dims) /= size(grid)) return
        if (any(dims /= grid)) return
        if (size(res) == 0) return
        valid = .true.
        call elementwise_broadcast_optional(xl0, grid, left, valid)
        call elementwise_broadcast_optional(xr0, grid, right, valid)
        call elementwise_broadcast_optional(xmin, grid, lo, valid)
        call elementwise_broadcast_optional(xmax, grid, hi, valid)
        call elementwise_broadcast_optional(factor, grid, grow, valid)
        if (present(args)) call elementwise_broadcast_arguments(args, grid, packed_args, valid)
        if (present(args_list)) call elementwise_broadcast_argument_list(args_list, grid, packed_args, valid)
        if (present(typed_args)) then
            call elementwise_broadcast_typed_argument_list(typed_args, grid, typed_indices, valid)
        end if
        if (.not. valid) return
        call elementwise_flatten(xm0, linear)
        flat = bracket_minimum_vectorized_rank1(f, linear, &
            xl0=left, xr0=right, xmin=lo, xmax=hi, factor=grow, maxiter=maxiter, args_data=packed_args, f_args=f_args, &
            typed_args=typed_args, typed_indices=typed_indices, f_typed=f_typed, f_indexed=f_indexed)
        ! SELECT RANK is required because Fortran 2018 has no assumed-rank
        ! function result or rank-independent RESHAPE assignment.
        select rank(res)
        rank(0)
            res = flat(1)
        rank(1)
            res = reshape(flat, shape(res))
        rank(2)
            res = reshape(flat, shape(res))
        rank(3)
            res = reshape(flat, shape(res))
        rank(4)
            res = reshape(flat, shape(res))
        rank(5)
            res = reshape(flat, shape(res))
        rank(6)
            res = reshape(flat, shape(res))
        rank(7)
            res = reshape(flat, shape(res))
        rank(8)
            res = reshape(flat, shape(res))
        rank(9)
            res = reshape(flat, shape(res))
        rank(10)
            res = reshape(flat, shape(res))
        rank(11)
            res = reshape(flat, shape(res))
        rank(12)
            res = reshape(flat, shape(res))
        rank(13)
            res = reshape(flat, shape(res))
        rank(14)
            res = reshape(flat, shape(res))
        rank(15)
            res = reshape(flat, shape(res))
        rank default
            return
        end select
    end subroutine bracket_minimum_vectorized_into


    function find_minimum_vectorized_rank3(f, bracket, xatol, xrtol, fatol, frtol, maxiter, &
            callback) result(res)
        ! Compatibility function; numerical work lives in the assumed-rank dispatcher.
        procedure(elementwise_minimum_vector_function) :: f !! packed objective callback
        real(dp), intent(in) :: bracket(:,:,:,:) !! input grid
        real(dp), intent(in), optional :: xatol(..) !! broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_minimum_callback), optional :: callback !! scalar per-lane callback
        type(find_minimum_result) :: res(size(bracket,2),size(bracket,3),size(bracket,4))

        call find_minimum_vectorized_into(f, bracket, res, &
            xatol, xrtol, fatol, frtol, maxiter, callback)
    end function find_minimum_vectorized_rank3


    function find_minimum_vectorized_rank4(f, bracket, xatol, xrtol, fatol, frtol, maxiter, &
            callback) result(res)
        ! Compatibility function; numerical work lives in the assumed-rank dispatcher.
        procedure(elementwise_minimum_vector_function) :: f !! packed objective callback
        real(dp), intent(in) :: bracket(:,:,:,:,:) !! input grid
        real(dp), intent(in), optional :: xatol(..) !! broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_minimum_callback), optional :: callback !! scalar per-lane callback
        type(find_minimum_result) :: res(size(bracket,2),size(bracket,3),size(bracket,4),size(bracket,5))

        call find_minimum_vectorized_into(f, bracket, res, &
            xatol, xrtol, fatol, frtol, maxiter, callback)
    end function find_minimum_vectorized_rank4


    function bracket_minimum_vectorized_rank3(f, xm0, xl0, xr0, xmin, xmax, factor, &
            maxiter) result(res)
        ! Compatibility function; numerical work lives in the assumed-rank dispatcher.
        procedure(elementwise_minimum_vector_function) :: f !! packed objective callback
        real(dp), intent(in) :: xm0(:,:,:) !! input grid
        real(dp), intent(in), optional :: xl0(..) !! broadcast parameter
        real(dp), intent(in), optional :: xr0(..) !! broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! broadcast parameter
        real(dp), intent(in), optional :: factor(..) !! broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_minimum_result) :: res(size(xm0,1),size(xm0,2),size(xm0,3))

        call bracket_minimum_vectorized_into(f, xm0, res, &
            xl0, xr0, xmin, xmax, factor, maxiter)
    end function bracket_minimum_vectorized_rank3


    function bracket_minimum_vectorized_rank4(f, xm0, xl0, xr0, xmin, xmax, factor, &
            maxiter) result(res)
        ! Compatibility function; numerical work lives in the assumed-rank dispatcher.
        procedure(elementwise_minimum_vector_function) :: f !! packed objective callback
        real(dp), intent(in) :: xm0(:,:,:,:) !! input grid
        real(dp), intent(in), optional :: xl0(..) !! broadcast parameter
        real(dp), intent(in), optional :: xr0(..) !! broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! broadcast parameter
        real(dp), intent(in), optional :: factor(..) !! broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_minimum_result) :: res(size(xm0,1),size(xm0,2),size(xm0,3),size(xm0,4))

        call bracket_minimum_vectorized_into(f, xm0, res, &
            xl0, xr0, xmin, xmax, factor, maxiter)
    end function bracket_minimum_vectorized_rank4

    function find_minimum_vectorized_rank2(f, bracket, xatol, xrtol, fatol, frtol, maxiter, &
            callback) result(res)
        ! Compatibility function; numerical work lives in the assumed-rank dispatcher.
        procedure(elementwise_minimum_vector_function) :: f !! packed objective callback
        real(dp), intent(in) :: bracket(:,:,:) !! input grid
        real(dp), intent(in), optional :: xatol(..) !! broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_minimum_callback), optional :: callback !! scalar per-lane callback
        type(find_minimum_result) :: res(size(bracket,2),size(bracket,3))

        call find_minimum_vectorized_into(f, bracket, res, &
            xatol, xrtol, fatol, frtol, maxiter, callback)
    end function find_minimum_vectorized_rank2


    function bracket_minimum_vectorized_rank2(f, xm0, xl0, xr0, xmin, xmax, factor, &
            maxiter) result(res)
        ! Compatibility function; numerical work lives in the assumed-rank dispatcher.
        procedure(elementwise_minimum_vector_function) :: f !! packed objective callback
        real(dp), intent(in) :: xm0(:,:) !! input grid
        real(dp), intent(in), optional :: xl0(..) !! broadcast parameter
        real(dp), intent(in), optional :: xr0(..) !! broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! broadcast parameter
        real(dp), intent(in), optional :: factor(..) !! broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_minimum_result) :: res(size(xm0,1),size(xm0,2))

        call bracket_minimum_vectorized_into(f, xm0, res, &
            xl0, xr0, xmin, xmax, factor, maxiter)
    end function bracket_minimum_vectorized_rank2


    function find_minimum_vectorized_rank1(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback, args_data, f_args, &
            typed_args, typed_indices, f_typed, f_indexed) result(res)
        ! Chandrupatla minimum refinement with packed vector evaluations.
        ! Each lane runs the same scalar state machine; active lanes are packed
        ! before each call, so the callback must be elementwise and lane-free.
        procedure(elementwise_minimum_vector_function), optional :: f !! packed vector objective
        real(dp), intent(in), optional :: args_data(:,:) !! prebroadcast (nargs,nproblems)
        procedure(elementwise_minimum_args_function), optional :: f_args !! packed parameterized objective
        type(elementwise_typed_argument), intent(in), optional :: typed_args(:) !! typed input snapshots
        integer, intent(in), optional :: typed_indices(:,:) !! source offsets for every lane
        procedure(elementwise_minimum_typed_args_function), optional :: f_typed !! typed packed callback
        procedure(elementwise_minimum_indexed_args_function), optional :: f_indexed !! zero-copy typed callback
        real(dp), intent(in) :: bracket(:,:) !! starting triples (3,n)
        real(dp), intent(in), optional :: xatol(..) !! broadcast absolute x tolerance
        real(dp), intent(in), optional :: xrtol(..) !! broadcast relative x tolerance
        real(dp), intent(in), optional :: fatol(..) !! broadcast absolute f tolerance
        real(dp), intent(in), optional :: frtol(..) !! broadcast relative f tolerance
        integer, intent(in), optional :: maxiter !! common maximum iterations
        procedure(find_minimum_callback), optional :: callback !! scalar state callback
        type(find_minimum_result) :: res(size(bracket,2))

        real(dp), parameter :: phi = (1.0_dp + sqrt(5.0_dp)) / 2.0_dp
        real(dp), allocatable :: ax(:), rx(:), af(:), rf(:)
        real(dp), allocatable :: x1(:), x2(:), x3(:), f1(:), f2(:), f3(:)
        real(dp), allocatable :: q0(:), xtol(:), points(:), values(:)
        type(elementwise_typed_argument), allocatable :: active_typed(:)
        integer, allocatable :: packed_typed_indices(:,:)
        logical, allocatable :: active(:)
        integer, allocatable :: ids(:)
        integer, allocatable :: packed_ids(:)
        real(dp) :: x(3), fv(3), tmp, xa, xr, a, b, c, q1, trial
        integer :: n, niter, i, j, k, h, m, it, state
        logical :: valid, stopped

        if (present(f_args) .neqv. present(args_data)) return
        if (present(f_typed) .and. present(f_indexed)) return
        if ((present(f_typed) .or. present(f_indexed)) .neqv. present(typed_args)) return
        if ((present(f_typed) .or. present(f_indexed)) .neqv. present(typed_indices)) return
        if ((present(f_typed) .or. present(f_indexed)) .and. &
                (present(f) .or. present(f_args))) return
        if (.not. present(f) .and. .not. present(f_args) .and. &
                .not. present(f_typed) .and. .not. present(f_indexed)) return
        if (present(args_data)) then
            if (size(args_data,1) < 1) return
            if (size(args_data,2) /= size(bracket,2)) return
        end if
        if (present(typed_args)) then
            if (size(typed_args) == 0) return
            if (size(typed_indices,1) /= size(typed_args)) return
            if (size(typed_indices,2) /= size(bracket,2)) return
        end if
        n = size(bracket,2)
        if (size(bracket,1) /= 3 .or. n == 0) return
        niter = 100
        if (present(maxiter)) niter = maxiter
        if (niter < 0) return

        allocate(ax(n), rx(n), af(n), rf(n))
        allocate(x1(n), x2(n), x3(n), f1(n), f2(n), f3(n))
        allocate(q0(n), xtol(n), active(n), ids(n))
        allocate(points(3*n), values(3*n))
        if (present(f_indexed)) allocate(packed_typed_indices(size(typed_indices,1), 3*n))
        if (present(f_args) .or. present(f_typed) .or. present(f_indexed)) allocate(packed_ids(3*n))
        ax = tiny(1.0_dp)
        rx = sqrt(epsilon(1.0_dp))
        af = tiny(1.0_dp)
        rf = tiny(1.0_dp)
        valid = .true.
        call elementwise_broadcast(xatol, [n], ax, valid)
        call elementwise_broadcast(xrtol, [n], rx, valid)
        call elementwise_broadcast(fatol, [n], af, valid)
        call elementwise_broadcast(frtol, [n], rf, valid)
        if (.not. valid) return

        ! Exclude invalid lanes before the first callback. Check finiteness
        ! before comparisons to avoid NaN comparisons under IEEE traps.
        active = .false.
        m = 0
        do i = 1, n
            if (.not. all(ieee_is_finite(bracket(:,i)))) cycle
            if (.not. ieee_is_finite(ax(i))) cycle
            if (.not. ieee_is_finite(rx(i))) cycle
            if (.not. ieee_is_finite(af(i))) cycle
            if (.not. ieee_is_finite(rf(i))) cycle
            if (ax(i) < 0.0_dp .or. rx(i) < 0.0_dp) cycle
            if (af(i) < 0.0_dp .or. rf(i) < 0.0_dp) cycle
            active(i) = .true.
            m = m + 1
            ids(m) = i
            points(3*m-2:3*m) = bracket(:,i)
            if (present(f_args) .or. present(f_typed) .or. present(f_indexed)) packed_ids(3*m-2:3*m) = i
        end do
        if (m == 0) return

        ! Evaluate all three endpoints for each valid lane in supplied order.
        if (present(f_indexed)) then
            call elementwise_gather_typed_indices(typed_indices, packed_ids(:3*m), &
                packed_typed_indices(:, :3*m))
            call f_indexed(points(:3*m), typed_args, &
                packed_typed_indices(:, :3*m), values(:3*m))
        else if (present(f_typed)) then
            call elementwise_pack_typed_arguments(typed_args, &
                typed_indices(:,packed_ids(:3*m)), active_typed)
            call f_typed(points(:3*m), active_typed, values(:3*m))
        else if (present(f_args)) then
            call f_args(points(:3*m), args_data(:,packed_ids(:3*m)), values(:3*m))
        else
            call f(points(:3*m), values(:3*m))
        end if
        do j = 1, m
            i = ids(j)
            x = bracket(:,i)
            fv = values(3*j-2:3*j)
            res(i)%nfev = 3
            ! Match the scalar solver's stable sort of x and f(x) together.
            do k = 2, 3
                do h = k, 2, -1
                    if (x(h) >= x(h-1)) exit
                    tmp = x(h)
                    x(h) = x(h-1)
                    x(h-1) = tmp
                    tmp = fv(h)
                    fv(h) = fv(h-1)
                    fv(h-1) = tmp
                end do
            end do
            x1(i) = x(1)
            x2(i) = x(2)
            x3(i) = x(3)
            f1(i) = fv(1)
            f2(i) = fv(2)
            f3(i) = fv(3)
            q0(i) = x3(i)
            if (.not. all(ieee_is_finite(fv))) then
                res(i)%status = minimize_status_value_error
                call store_minimum_state(res(i), x1(i), x2(i), x3(i), f1(i), f2(i), f3(i))
                call invalidate_minimum(res(i))
                active(i) = .false.
                cycle
            end if
            if (x1(i) >= x2(i) .or. x2(i) >= x3(i) .or. &
                    f2(i) > f1(i) .or. f2(i) > f3(i)) then
                res(i)%status = minimize_status_bracket_error
                call store_minimum_state(res(i), x1(i), x2(i), x3(i), f1(i), f2(i), f3(i))
                call invalidate_minimum(res(i))
                active(i) = .false.
                cycle
            end if
            call check_find_minimum(x1(i), x2(i), x3(i), f1(i), f2(i), f3(i), &
                ax(i), rx(i), af(i), rf(i), xtol(i), state)
            res(i)%status = state
            call store_minimum_state(res(i), x1(i), x2(i), x3(i), f1(i), f2(i), f3(i))
            if (present(callback)) then
                stopped = .false.
                call callback(res(i), stopped)
                if (stopped .and. res(i)%status == 1) res(i)%status = find_minimum_status_callback
            end if
            active(i) = res(i)%status == 1
        end do

        do it = 1, niter
            m = 0
            do i = 1, n
                if (.not. active(i)) cycle
                xa = x2(i) - x1(i)
                xr = x3(i) - x2(i)
                a = xa * (f3(i) - f2(i))
                b = xr * (f1(i) - f2(i))
                if (a+b == 0.0_dp) then
                    q1 = x2(i)
                else
                    c = a / (a+b)
                    q1 = 0.5_dp * (c * (x1(i)-x3(i)) + x2(i) + x3(i))
                end if
                trial = x2(i) + (2.0_dp-phi)*xr
                if (ieee_is_finite(q1)) then
                    if (abs(q1-q0(i)) < 0.5_dp*abs(xa)) then
                        trial = q1
                        if (abs(trial-x2(i)) <= xtol(i)) trial = x2(i) + sign(xtol(i),xr)
                    end if
                end if
                q0(i) = q1
                if (trial == x2(i) .or. trial == x1(i) .or. trial == x3(i)) then
                    trial = x2(i) + (2.0_dp-phi)*xr
                    if (trial == x2(i) .or. trial == x3(i)) then
                        res(i)%status = minimize_status_convergence_error
                        active(i) = .false.
                        cycle
                    end if
                end if
                m = m + 1
                ids(m) = i
                points(m) = trial
            end do
            if (m == 0) exit
            if (present(f_indexed)) then
                call elementwise_gather_typed_indices(typed_indices, ids(:m), &
                    packed_typed_indices(:, :m))
                call f_indexed(points(:m), typed_args, &
                    packed_typed_indices(:, :m), values(:m))
            else if (present(f_typed)) then
                    call elementwise_pack_typed_arguments(typed_args, &
                        typed_indices(:,ids(:m)), active_typed)
                    call f_typed(points(:m), active_typed, values(:m))
            else if (present(f_args)) then
                    call f_args(points(:m), args_data(:,ids(:m)), values(:m))
            else
                    call f(points(:m), values(:m))
                end if
            do j = 1, m
                i = ids(j)
                trial = points(j)
                res(i)%nit = res(i)%nit + 1
                res(i)%nfev = res(i)%nfev + 1
                if (.not. ieee_is_finite(values(j))) then
                    res(i)%status = minimize_status_value_error
                    call invalidate_minimum(res(i))
                    active(i) = .false.
                    cycle
                end if
                if ((trial > x2(i)) .eqv. (x3(i) > x2(i))) then
                    if (values(j) > f2(i)) then
                        x3(i) = trial
                        f3(i) = values(j)
                    else
                        x1(i) = x2(i)
                        f1(i) = f2(i)
                        x2(i) = trial
                        f2(i) = values(j)
                    end if
                else
                    if (values(j) > f2(i)) then
                        x1(i) = trial
                        f1(i) = values(j)
                    else
                        x3(i) = x2(i)
                        f3(i) = f2(i)
                        x2(i) = trial
                        f2(i) = values(j)
                    end if
                end if
                call check_find_minimum(x1(i), x2(i), x3(i), f1(i), f2(i), f3(i), &
                    ax(i), rx(i), af(i), rf(i), xtol(i), state)
                res(i)%status = state
                call store_minimum_state(res(i), x1(i), x2(i), x3(i), f1(i), f2(i), f3(i))
                if (present(callback)) then
                    stopped = .false.
                    call callback(res(i), stopped)
                    if (stopped .and. res(i)%status == 1) res(i)%status = find_minimum_status_callback
                end if
                active(i) = res(i)%status == 1
            end do
        end do
        do i = 1, n
            if (res(i)%status == 1) res(i)%status = minimize_status_convergence_error
            res(i)%success = res(i)%status == minimize_status_success
        end do
    end function find_minimum_vectorized_rank1

    subroutine find_minimum_batch(f, bracket, batch_shape, xatol, xrtol, fatol, frtol, maxiter, callback, &
            res)
        implicit none
        ! Subroutine dispatch avoids temporary descriptors for absent assumed-rank
        ! arguments in array-valued function calls on affected compilers.
        ! One independent scalar solve per element; reshape by thin rank wrappers.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:, :) !! flattened problem inputs
        integer, intent(in) :: batch_shape(:) !! desired output dimensions
        real(dp), intent(in), optional :: xatol(..) !! broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! broadcast parameter
        integer, intent(in), optional :: maxiter !! shared iteration limit
        procedure(find_minimum_callback), optional :: callback !! scalar callback
        type(find_minimum_result), allocatable, intent(out) :: res(:) !! flattened independent solver results
        logical :: valid
        integer :: i, n
        real(dp), allocatable :: ax(:), rx(:), af(:), rf(:)

        n = size(bracket, 2)
        allocate(res(n), ax(n), rx(n), af(n), rf(n))
        if (size(bracket, 1) /= 3) return
        if (n /= product(batch_shape)) return
        ax = tiny(1.0_dp)
        rx = sqrt(epsilon(1.0_dp))
        af = tiny(1.0_dp)
        rf = tiny(1.0_dp)
        valid = .true.
        call elementwise_broadcast(xatol, batch_shape, ax, valid)
        call elementwise_broadcast(xrtol, batch_shape, rx, valid)
        call elementwise_broadcast(fatol, batch_shape, af, valid)
        call elementwise_broadcast(frtol, batch_shape, rf, valid)
        if (.not. valid) return
        do i = 1, n
            res(i) = find_minimum_scalar(f, bracket(:, i), xatol=ax(i), xrtol=rx(i), &
                fatol=af(i), frtol=rf(i), maxiter=maxiter, callback=callback)
        end do
    end subroutine find_minimum_batch

    function find_minimum_vector(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 1-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_minimum_callback), optional :: callback !! per-element callback
        type(find_minimum_result), allocatable :: flat_result(:)
        type(find_minimum_result) :: res(size(bracket, 2))

        call find_minimum_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_minimum_vector

    function find_minimum_matrix(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 2-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_minimum_callback), optional :: callback !! per-element callback
        type(find_minimum_result), allocatable :: flat_result(:)
        type(find_minimum_result) :: res(size(bracket, 2), size(bracket, 3))

        call find_minimum_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_minimum_matrix

    function find_minimum_tensor(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 3-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_minimum_callback), optional :: callback !! per-element callback
        type(find_minimum_result), allocatable :: flat_result(:)
        type(find_minimum_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4))

        call find_minimum_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_minimum_tensor

    function find_minimum_rank4(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 4-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_minimum_callback), optional :: callback !! per-element callback
        type(find_minimum_result), allocatable :: flat_result(:)
        type(find_minimum_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), size(bracket, 5))

        call find_minimum_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_minimum_rank4

    function find_minimum_rank5(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 5-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_minimum_callback), optional :: callback !! per-element callback
        type(find_minimum_result), allocatable :: flat_result(:)
        type(find_minimum_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), size(bracket, 5), size(bracket, 6))

        call find_minimum_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_minimum_rank5

    function find_minimum_rank6(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 6-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_minimum_callback), optional :: callback !! per-element callback
        type(find_minimum_result), allocatable :: flat_result(:)
        type(find_minimum_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), &
            size(bracket, 5), size(bracket, 6), size(bracket, 7))

        call find_minimum_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_minimum_rank6

    function find_minimum_rank7(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 7-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_minimum_callback), optional :: callback !! per-element callback
        type(find_minimum_result), allocatable :: flat_result(:)
        type(find_minimum_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), &
            size(bracket, 5), size(bracket, 6), size(bracket, 7), size(bracket, 8))

        call find_minimum_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_minimum_rank7

    function find_minimum_rank8(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 8-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_minimum_callback), optional :: callback !! per-element callback
        type(find_minimum_result), allocatable :: flat_result(:)
        type(find_minimum_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), &
            size(bracket, 5), size(bracket, 6), size(bracket, 7), size(bracket, 8), size(bracket, 9))

        call find_minimum_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_minimum_rank8

    function find_minimum_rank9(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 9-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:,:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_minimum_callback), optional :: callback !! per-element callback
        type(find_minimum_result), allocatable :: flat_result(:)
        type(find_minimum_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), &
            size(bracket, 5), size(bracket, 6), size(bracket, 7), size(bracket, 8), size(bracket, 9), &
            size(bracket, 10))

        call find_minimum_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_minimum_rank9

    function find_minimum_rank10(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 10-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:,:,:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_minimum_callback), optional :: callback !! per-element callback
        type(find_minimum_result), allocatable :: flat_result(:)
        type(find_minimum_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), &
            size(bracket, 5), size(bracket, 6), size(bracket, 7), size(bracket, 8), size(bracket, 9), &
            size(bracket, 10), size(bracket, 11))

        call find_minimum_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_minimum_rank10

    function find_minimum_rank11(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 11-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:,:,:,:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_minimum_callback), optional :: callback !! per-element callback
        type(find_minimum_result), allocatable :: flat_result(:)
        type(find_minimum_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), &
            size(bracket, 5), size(bracket, 6), size(bracket, 7), size(bracket, 8), size(bracket, 9), &
            size(bracket, 10), size(bracket, 11), size(bracket, 12))

        call find_minimum_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_minimum_rank11

    function find_minimum_rank12(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 12-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:,:,:,:,:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_minimum_callback), optional :: callback !! per-element callback
        type(find_minimum_result), allocatable :: flat_result(:)
        type(find_minimum_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), &
            size(bracket, 5), size(bracket, 6), size(bracket, 7), size(bracket, 8), size(bracket, 9), &
            size(bracket, 10), size(bracket, 11), size(bracket, 12), size(bracket, 13))

        call find_minimum_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_minimum_rank12

    function find_minimum_rank13(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 13-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_minimum_callback), optional :: callback !! per-element callback
        type(find_minimum_result), allocatable :: flat_result(:)
        type(find_minimum_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), &
            size(bracket, 5), size(bracket, 6), size(bracket, 7), size(bracket, 8), size(bracket, 9), &
            size(bracket, 10), size(bracket, 11), size(bracket, 12), size(bracket, 13), size(bracket, 14))

        call find_minimum_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_minimum_rank13

    function find_minimum_rank14(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 14-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_minimum_callback), optional :: callback !! per-element callback
        type(find_minimum_result), allocatable :: flat_result(:)
        type(find_minimum_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), &
            size(bracket, 5), size(bracket, 6), size(bracket, 7), size(bracket, 8), size(bracket, 9), &
            size(bracket, 10), size(bracket, 11), size(bracket, 12), size(bracket, 13), size(bracket, 14), &
            size(bracket, 15))

        call find_minimum_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_minimum_rank14

    subroutine bracket_minimum_batch(f, xm0, batch_shape, xl0, xr0, xmin, xmax, factor, maxiter, &
            res)
        implicit none
        ! Subroutine dispatch avoids temporary descriptors for absent assumed-rank
        ! arguments in array-valued function calls on affected compilers.
        ! One independent scalar solve per element; reshape by thin rank wrappers.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: xm0(:) !! flattened problem inputs
        integer, intent(in) :: batch_shape(:) !! desired output dimensions
        real(dp), intent(in), optional :: xl0(..) !! broadcast parameter
        real(dp), intent(in), optional :: xr0(..) !! broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! broadcast parameter
        real(dp), intent(in), optional :: factor !! shared growth factor
        integer, intent(in), optional :: maxiter !! shared iteration limit
        type(bracket_minimum_result), allocatable, intent(out) :: res(:) !! flattened independent solver results
        logical :: valid
        integer :: i, n
        real(dp), allocatable :: left(:), right(:), lo(:), hi(:)

        n = size(xm0)
        allocate(res(n), left(n), right(n), lo(n), hi(n))
        if (n /= product(batch_shape)) return
        left = 0.0_dp
        right = 0.0_dp
        lo = ieee_value(0.0_dp, ieee_negative_inf)
        hi = ieee_value(0.0_dp, ieee_positive_inf)
        valid = .true.
        call elementwise_broadcast(xl0, batch_shape, left, valid)
        call elementwise_broadcast(xr0, batch_shape, right, valid)
        call elementwise_broadcast(xmin, batch_shape, lo, valid)
        call elementwise_broadcast(xmax, batch_shape, hi, valid)
        if (.not. valid) return
        do i = 1, n
            if (present(xl0)) then
                if (present(xr0)) then
                    res(i) = bracket_minimum_scalar(f, xm0(i), xl0=left(i), xr0=right(i), &
                        xmin=lo(i), xmax=hi(i), factor=factor, maxiter=maxiter)
                else
                    res(i) = bracket_minimum_scalar(f, xm0(i), xl0=left(i), &
                        xmin=lo(i), xmax=hi(i), factor=factor, maxiter=maxiter)
                end if
            else if (present(xr0)) then
                res(i) = bracket_minimum_scalar(f, xm0(i), xr0=right(i), &
                    xmin=lo(i), xmax=hi(i), factor=factor, maxiter=maxiter)
            else
                res(i) = bracket_minimum_scalar(f, xm0(i), &
                    xmin=lo(i), xmax=hi(i), factor=factor, maxiter=maxiter)
            end if
        end do
    end subroutine bracket_minimum_batch

    function bracket_minimum_vector(f, xm0, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 1-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: xm0(:) !! problem starting points
        real(dp), intent(in), optional :: xl0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_minimum_result), allocatable :: flat_result(:)
        type(bracket_minimum_result) :: res(size(xm0, 1))

        call bracket_minimum_batch(f, reshape(xm0, [size(xm0)]), &
            shape(res), xl0=xl0, xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_minimum_vector

    function bracket_minimum_matrix(f, xm0, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 2-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: xm0(:,:) !! problem starting points
        real(dp), intent(in), optional :: xl0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_minimum_result), allocatable :: flat_result(:)
        type(bracket_minimum_result) :: res(size(xm0, 1), size(xm0, 2))

        call bracket_minimum_batch(f, reshape(xm0, [size(xm0)]), &
            shape(res), xl0=xl0, xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_minimum_matrix

    function bracket_minimum_tensor(f, xm0, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 3-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: xm0(:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xl0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_minimum_result), allocatable :: flat_result(:)
        type(bracket_minimum_result) :: res(size(xm0, 1), size(xm0, 2), size(xm0, 3))

        call bracket_minimum_batch(f, reshape(xm0, [size(xm0)]), &
            shape(res), xl0=xl0, xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_minimum_tensor

    function bracket_minimum_rank4(f, xm0, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 4-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: xm0(:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xl0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_minimum_result), allocatable :: flat_result(:)
        type(bracket_minimum_result) :: res(size(xm0, 1), size(xm0, 2), size(xm0, 3), size(xm0, 4))

        call bracket_minimum_batch(f, reshape(xm0, [size(xm0)]), &
            shape(res), xl0=xl0, xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_minimum_rank4

    function bracket_minimum_rank5(f, xm0, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 5-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: xm0(:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xl0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_minimum_result), allocatable :: flat_result(:)
        type(bracket_minimum_result) :: res(size(xm0, 1), size(xm0, 2), size(xm0, 3), size(xm0, 4), size(xm0, 5))

        call bracket_minimum_batch(f, reshape(xm0, [size(xm0)]), &
            shape(res), xl0=xl0, xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_minimum_rank5

    function bracket_minimum_rank6(f, xm0, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 6-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: xm0(:,:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xl0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_minimum_result), allocatable :: flat_result(:)
        type(bracket_minimum_result) :: res(size(xm0, 1), size(xm0, 2), size(xm0, 3), size(xm0, 4), size(xm0, 5), size(xm0, 6))

        call bracket_minimum_batch(f, reshape(xm0, [size(xm0)]), &
            shape(res), xl0=xl0, xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_minimum_rank6

    function bracket_minimum_rank7(f, xm0, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 7-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: xm0(:,:,:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xl0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_minimum_result), allocatable :: flat_result(:)
        type(bracket_minimum_result) :: res(size(xm0, 1), size(xm0, 2), size(xm0, 3), size(xm0, 4), &
            size(xm0, 5), size(xm0, 6), size(xm0, 7))

        call bracket_minimum_batch(f, reshape(xm0, [size(xm0)]), &
            shape(res), xl0=xl0, xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_minimum_rank7

    function bracket_minimum_rank8(f, xm0, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 8-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: xm0(:,:,:,:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xl0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_minimum_result), allocatable :: flat_result(:)
        type(bracket_minimum_result) :: res(size(xm0, 1), size(xm0, 2), size(xm0, 3), size(xm0, 4), &
            size(xm0, 5), size(xm0, 6), size(xm0, 7), size(xm0, 8))

        call bracket_minimum_batch(f, reshape(xm0, [size(xm0)]), &
            shape(res), xl0=xl0, xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_minimum_rank8

    function bracket_minimum_rank9(f, xm0, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 9-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: xm0(:,:,:,:,:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xl0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_minimum_result), allocatable :: flat_result(:)
        type(bracket_minimum_result) :: res(size(xm0, 1), size(xm0, 2), size(xm0, 3), size(xm0, 4), &
            size(xm0, 5), size(xm0, 6), size(xm0, 7), size(xm0, 8), size(xm0, 9))

        call bracket_minimum_batch(f, reshape(xm0, [size(xm0)]), &
            shape(res), xl0=xl0, xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_minimum_rank9

    function bracket_minimum_rank10(f, xm0, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 10-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: xm0(:,:,:,:,:,:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xl0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_minimum_result), allocatable :: flat_result(:)
        type(bracket_minimum_result) :: res(size(xm0, 1), size(xm0, 2), size(xm0, 3), size(xm0, 4), &
            size(xm0, 5), size(xm0, 6), size(xm0, 7), size(xm0, 8), size(xm0, 9), size(xm0, 10))

        call bracket_minimum_batch(f, reshape(xm0, [size(xm0)]), &
            shape(res), xl0=xl0, xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_minimum_rank10

    function bracket_minimum_rank11(f, xm0, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 11-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: xm0(:,:,:,:,:,:,:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xl0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_minimum_result), allocatable :: flat_result(:)
        type(bracket_minimum_result) :: res(size(xm0, 1), size(xm0, 2), size(xm0, 3), size(xm0, 4), &
            size(xm0, 5), size(xm0, 6), size(xm0, 7), size(xm0, 8), size(xm0, 9), size(xm0, 10), &
            size(xm0, 11))

        call bracket_minimum_batch(f, reshape(xm0, [size(xm0)]), &
            shape(res), xl0=xl0, xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_minimum_rank11

    function bracket_minimum_rank12(f, xm0, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 12-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: xm0(:,:,:,:,:,:,:,:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xl0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_minimum_result), allocatable :: flat_result(:)
        type(bracket_minimum_result) :: res(size(xm0, 1), size(xm0, 2), size(xm0, 3), size(xm0, 4), &
            size(xm0, 5), size(xm0, 6), size(xm0, 7), size(xm0, 8), size(xm0, 9), size(xm0, 10), &
            size(xm0, 11), size(xm0, 12))

        call bracket_minimum_batch(f, reshape(xm0, [size(xm0)]), &
            shape(res), xl0=xl0, xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_minimum_rank12

    function bracket_minimum_rank13(f, xm0, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 13-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: xm0(:,:,:,:,:,:,:,:,:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xl0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_minimum_result), allocatable :: flat_result(:)
        type(bracket_minimum_result) :: res(size(xm0, 1), size(xm0, 2), size(xm0, 3), size(xm0, 4), &
            size(xm0, 5), size(xm0, 6), size(xm0, 7), size(xm0, 8), size(xm0, 9), size(xm0, 10), &
            size(xm0, 11), size(xm0, 12), size(xm0, 13))

        call bracket_minimum_batch(f, reshape(xm0, [size(xm0)]), &
            shape(res), xl0=xl0, xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_minimum_rank13

    function bracket_minimum_rank14(f, xm0, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 14-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: xm0(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xl0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_minimum_result), allocatable :: flat_result(:)
        type(bracket_minimum_result) :: res(size(xm0, 1), size(xm0, 2), size(xm0, 3), size(xm0, 4), &
            size(xm0, 5), size(xm0, 6), size(xm0, 7), size(xm0, 8), size(xm0, 9), size(xm0, 10), &
            size(xm0, 11), size(xm0, 12), size(xm0, 13), size(xm0, 14))

        call bracket_minimum_batch(f, reshape(xm0, [size(xm0)]), &
            shape(res), xl0=xl0, xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_minimum_rank14

    function bracket_minimum_rank15(f, xm0, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 15-dimensional batch, via the shared flat engine.
        procedure(scalar_objective_function) :: f !! common scalar function
        real(dp), intent(in) :: xm0(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xl0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_minimum_result), allocatable :: flat_result(:)
        type(bracket_minimum_result) :: res(size(xm0, 1), size(xm0, 2), size(xm0, 3), size(xm0, 4), &
            size(xm0, 5), size(xm0, 6), size(xm0, 7), size(xm0, 8), size(xm0, 9), size(xm0, 10), &
            size(xm0, 11), size(xm0, 12), size(xm0, 13), size(xm0, 14), size(xm0, 15))

        call bracket_minimum_batch(f, reshape(xm0, [size(xm0)]), &
            shape(res), xl0=xl0, xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_minimum_rank15

    function find_minimum_scalar(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Scalar version of SciPy elementwise.find_minimum (Chandrupatla).
        ! The initial bracket must have three abscissae around a local minimum.
        procedure(scalar_objective_function) :: f !! scalar objective
        real(dp), intent(in) :: bracket(:) !! three initial points, any order
        real(dp), intent(in), optional :: xatol !! absolute abscissa tolerance
        real(dp), intent(in), optional :: xrtol !! relative abscissa tolerance
        real(dp), intent(in), optional :: fatol !! absolute function tolerance
        real(dp), intent(in), optional :: frtol !! relative function tolerance
        integer, intent(in), optional :: maxiter !! maximum iterations (default 100)
        procedure(find_minimum_callback), optional :: callback !! optional stopping callback
        type(find_minimum_result) :: res

        real(dp) :: x(3), fval(3), tmp, q0, q1, c, a, b
        real(dp) :: x1, x2, x3, f1, f2, f3, xtol
        real(dp) :: xa, xr, trial, f_trial
        real(dp) :: atol_x, rtol_x, atol_f, rtol_f
        real(dp), parameter :: phi = (1.0_dp + sqrt(5.0_dp)) / 2.0_dp
        integer :: i, j, niter, status
        logical :: stopped

        if (size(bracket) /= 3) return
        if (.not. all(ieee_is_finite(bracket))) return
        atol_x = tiny(1.0_dp)
        rtol_x = sqrt(epsilon(1.0_dp))
        atol_f = tiny(1.0_dp)
        rtol_f = tiny(1.0_dp)
        if (present(xatol)) atol_x = xatol
        if (present(xrtol)) rtol_x = xrtol
        if (present(fatol)) atol_f = fatol
        if (present(frtol)) rtol_f = frtol
        niter = 100
        if (present(maxiter)) niter = maxiter
        if (niter < 0) return
        if (.not. all(ieee_is_finite([atol_x, rtol_x, atol_f, rtol_f]))) return
        if (min(atol_x, rtol_x, atol_f, rtol_f) < 0.0_dp) return

        ! SciPy evaluates in the supplied order, then sorts the abscissae.
        x = bracket
        do i = 1, 3
            fval(i) = f(x(i))
        end do
        res%nfev = 3
        do i = 2, 3
            do j = i, 2, -1
                if (x(j) >= x(j-1)) exit
                tmp = x(j)
                x(j) = x(j-1)
                x(j-1) = tmp
                tmp = fval(j)
                fval(j) = fval(j-1)
                fval(j-1) = tmp
            end do
        end do

        x1 = x(1)
        x2 = x(2)
        x3 = x(3)
        f1 = fval(1)
        f2 = fval(2)
        f3 = fval(3)
        q0 = x3
        if (.not. all(ieee_is_finite(fval))) then
            res%status = minimize_status_value_error
            call store_minimum_state(res, x1, x2, x3, f1, f2, f3)
            call invalidate_minimum(res)
            return
        end if
        if (x1 >= x2 .or. x2 >= x3 .or. f2 > f1 .or. f2 > f3) then
            res%status = minimize_status_bracket_error
            call store_minimum_state(res, x1, x2, x3, f1, f2, f3)
            call invalidate_minimum(res)
            return
        end if

        call check_find_minimum(x1, x2, x3, f1, f2, f3, atol_x, rtol_x, &
            atol_f, rtol_f, xtol, status)
        res%status = status
        call store_minimum_state(res, x1, x2, x3, f1, f2, f3)
        if (present(callback)) then
            stopped = .false.
            call callback(res, stopped)
            if (stopped .and. res%status == 1) res%status = find_minimum_status_callback
        end if

        do i = 1, niter
            if (res%status /= 1) exit
            ! The check above makes the interval (x2, x3) the larger one.
            xa = x2 - x1
            xr = x3 - x2
            a = xa * (f3 - f2)
            b = xr * (f1 - f2)
            if (a+b == 0.0_dp) then
                q1 = x2
            else
                c = a / (a+b)
                q1 = 0.5_dp * (c * (x1-x3) + x2 + x3)
            end if
            trial = x2 + (2.0_dp-phi)*xr
            if (ieee_is_finite(q1)) then
                if (abs(q1-q0) < 0.5_dp*abs(xa)) then
                    trial = q1
                    if (abs(trial-x2) <= xtol) trial = x2 + sign(xtol, xr)
                end if
            end if
            q0 = q1

            ! Avoid evaluating endpoints after rounding a proposed step.
            ! A fully collapsed bracket has already met the x tolerance.
            if (trial == x2 .or. trial == x1 .or. trial == x3) then
                trial = x2 + (2.0_dp-phi)*xr
                if (trial == x2 .or. trial == x3) then
                    res%status = minimize_status_convergence_error
                    exit
                end if
            end if
            f_trial = f(trial)
            res%nit = res%nit + 1
            res%nfev = res%nfev + 1
            if (.not. ieee_is_finite(f_trial)) then
                res%status = minimize_status_value_error
                exit
            end if

            ! Maintain f(x2) <= f(x1), f(x3); x3 may lie to the left.
            if ((trial > x2) .eqv. (x3 > x2)) then
                if (f_trial > f2) then
                    x3 = trial
                    f3 = f_trial
                else
                    x1 = x2
                    f1 = f2
                    x2 = trial
                    f2 = f_trial
                end if
            else
                if (f_trial > f2) then
                    x1 = trial
                    f1 = f_trial
                else
                    x3 = x2
                    f3 = f2
                    x2 = trial
                    f2 = f_trial
                end if
            end if

            call check_find_minimum(x1, x2, x3, f1, f2, f3, atol_x, rtol_x, &
                atol_f, rtol_f, xtol, status)
            res%status = status
            call store_minimum_state(res, x1, x2, x3, f1, f2, f3)
            if (present(callback)) then
                stopped = .false.
                call callback(res, stopped)
                if (stopped .and. res%status == 1) res%status = find_minimum_status_callback
            end if
        end do
        if (res%status == 1) res%status = minimize_status_convergence_error
        res%success = res%status == minimize_status_success
        if (res%status == minimize_status_value_error) call invalidate_minimum(res)
    end function find_minimum_scalar

    subroutine check_find_minimum(x1, x2, x3, f1, f2, f3, &
            atol_x, rtol_x, atol_f, rtol_f, xtol, status)
        ! Chandrupatla's test: either position or function bracket converged.
        real(dp), intent(inout) :: x1 !! first outer bracket point
        real(dp), intent(inout) :: x2 !! middle point
        real(dp), intent(inout) :: x3 !! second outer bracket point
        real(dp), intent(inout) :: f1 !! value at first outer point
        real(dp), intent(inout) :: f2 !! value at middle point
        real(dp), intent(inout) :: f3 !! value at second outer point
        real(dp), intent(in) :: atol_x !! absolute x tolerance
        real(dp), intent(in) :: rtol_x !! relative x tolerance
        real(dp), intent(in) :: atol_f !! absolute function tolerance
        real(dp), intent(in) :: rtol_f !! relative function tolerance
        real(dp), intent(out) :: xtol !! current absolute x tolerance
        integer, intent(out) :: status !! 0 converged or 1 in progress
        real(dp) :: tmp, ftol

        status = 1
        if (abs(x3-x2) < abs(x2-x1)) then
            tmp = x1
            x1 = x3
            x3 = tmp
            tmp = f1
            f1 = f3
            f3 = tmp
        end if
        xtol = abs(x2)*rtol_x + atol_x
        ftol = abs(f2)*rtol_f + atol_f
        if (abs(x3-x2) <= 2.0_dp*xtol .or. &
                f1 - 2.0_dp*f2 + f3 <= 2.0_dp*ftol) status = minimize_status_success
    end subroutine check_find_minimum

    subroutine store_minimum_state(res, x1, x2, x3, f1, f2, f3)
        type(find_minimum_result), intent(inout) :: res !! updated result
        real(dp), intent(in) :: x1 !! first outer bracket point
        real(dp), intent(in) :: x2 !! middle point
        real(dp), intent(in) :: x3 !! second outer bracket point
        real(dp), intent(in) :: f1 !! value at first point
        real(dp), intent(in) :: f2 !! value at middle point
        real(dp), intent(in) :: f3 !! value at second point

        res%x = x2
        res%f_x = f2
        res%xm = x2
        res%fm = f2
        if (x1 <= x3) then
            res%xl = x1
            res%fl = f1
            res%xr = x3
            res%fr = f3
        else
            res%xl = x3
            res%fl = f3
            res%xr = x1
            res%fr = f1
        end if
        res%success = res%status == minimize_status_success
    end subroutine store_minimum_state

    subroutine invalidate_minimum(res)
        type(find_minimum_result), intent(inout) :: res !! result with invalid minimum
        res%x = ieee_value(0.0_dp, ieee_quiet_nan)
        res%f_x = ieee_value(0.0_dp, ieee_quiet_nan)
        res%xm = res%x
        res%fm = res%f_x
        res%success = .false.
    end subroutine invalidate_minimum

    function bracket_minimum_vectorized_rank1(f, xm0, xl0, xr0, xmin, xmax, factor, maxiter, args_data, f_args, &
            typed_args, typed_indices, f_typed, f_indexed) result(res)
        ! Evaluate initial triples together, then pack only active trial points.
        ! Independent lanes reproduce bracket_minimum_scalar's search rules.
        procedure(elementwise_minimum_vector_function), optional :: f !! packed vector objective
        real(dp), intent(in), optional :: args_data(:,:) !! prebroadcast (nargs,nproblems)
        procedure(elementwise_minimum_args_function), optional :: f_args !! packed parameterized objective
        type(elementwise_typed_argument), intent(in), optional :: typed_args(:) !! typed input snapshots
        integer, intent(in), optional :: typed_indices(:,:) !! source offsets for every lane
        procedure(elementwise_minimum_typed_args_function), optional :: f_typed !! typed packed callback
        procedure(elementwise_minimum_indexed_args_function), optional :: f_indexed !! zero-copy typed callback
        real(dp), intent(in) :: xm0(:) !! initial middle abscissae
        real(dp), intent(in), optional :: xl0(..) !! broadcast initial left points
        real(dp), intent(in), optional :: xr0(..) !! broadcast initial right points
        real(dp), intent(in), optional :: xmin(..) !! broadcast lower search limits
        real(dp), intent(in), optional :: xmax(..) !! broadcast upper search limits
        real(dp), intent(in), optional :: factor(..) !! broadcast growth factor
        integer, intent(in), optional :: maxiter !! common number of expansions
        type(bracket_minimum_result) :: res(size(xm0))

        type(vector_bracket_minimum_lane), allocatable :: lane(:)
        real(dp), allocatable :: left(:), right(:), lo(:), hi(:), grow(:)
        real(dp), allocatable :: points(:), values(:)
        type(elementwise_typed_argument), allocatable :: active_typed(:)
        integer, allocatable :: packed_typed_indices(:,:)
        integer, allocatable :: ids(:), sides(:)
        logical, allocatable :: work(:)
        real(dp) :: space, scratch, next_x
        integer :: n, niter, i, j, k, m, it
        logical :: valid

        if (present(f_args) .neqv. present(args_data)) return
        if (present(f_typed) .and. present(f_indexed)) return
        if ((present(f_typed) .or. present(f_indexed)) .neqv. present(typed_args)) return
        if ((present(f_typed) .or. present(f_indexed)) .neqv. present(typed_indices)) return
        if ((present(f_typed) .or. present(f_indexed)) .and. &
                (present(f) .or. present(f_args))) return
        if (.not. present(f) .and. .not. present(f_args) .and. &
                .not. present(f_typed) .and. .not. present(f_indexed)) return
        if (present(args_data)) then
            if (size(args_data,1) < 1) return
            if (size(args_data,2) /= size(xm0)) return
        end if
        if (present(typed_args)) then
            if (size(typed_args) == 0) return
            if (size(typed_indices,1) /= size(typed_args)) return
            if (size(typed_indices,2) /= size(xm0)) return
        end if
        n = size(xm0)
        if (n == 0) return
        niter = 1000
        if (present(maxiter)) niter = maxiter
        if (niter < 0) return
        allocate(lane(n), work(n), left(n), right(n), lo(n), hi(n), grow(n))
        allocate(points(3*n), values(3*n), ids(3*n), sides(3*n))
        if (present(f_indexed)) allocate(packed_typed_indices(size(typed_indices,1), 3*n))
        left = 0.0_dp
        right = 0.0_dp
        lo = ieee_value(0.0_dp, ieee_negative_inf)
        hi = ieee_value(0.0_dp, ieee_positive_inf)
        grow = 2.0_dp
        valid = .true.
        call elementwise_broadcast(xl0, [n], left, valid)
        call elementwise_broadcast(xr0, [n], right, valid)
        call elementwise_broadcast(xmin, [n], lo, valid)
        call elementwise_broadcast(xmax, [n], hi, valid)
        call elementwise_broadcast(factor, [n], grow, valid)
        if (.not. valid) return

        work = .false.
        do i = 1, n
            if (.not. ieee_is_finite(xm0(i))) cycle
            if (ieee_is_nan(lo(i)) .or. ieee_is_nan(hi(i))) cycle
            if (.not. ieee_is_finite(grow(i))) cycle
            if (grow(i) <= 1.0_dp) cycle
            lane(i)%xm = xm0(i)
            lane(i)%xl = xm0(i) - 0.5_dp
            if (ieee_is_finite(lo(i))) then
                space = max(0.0_dp, xm0(i)/16.0_dp - lo(i)/16.0_dp)
                lane(i)%xl = xm0(i) - min(0.5_dp, space)
            end if
            lane(i)%xr = xm0(i) + 0.5_dp
            if (ieee_is_finite(hi(i))) then
                space = max(0.0_dp, hi(i)/16.0_dp - xm0(i)/16.0_dp)
                lane(i)%xr = xm0(i) + min(0.5_dp, space)
            end if
            if (present(xl0)) lane(i)%xl = left(i)
            if (present(xr0)) lane(i)%xr = right(i)
            if (.not. ieee_is_finite(lane(i)%xl) .or. .not. ieee_is_finite(lane(i)%xr)) cycle
            if (lane(i)%xl < lo(i) .or. lane(i)%xl >= xm0(i)) cycle
            if (lane(i)%xr > hi(i) .or. lane(i)%xr <= xm0(i)) cycle
            work(i) = .true.
        end do

        m = 0
        do i = 1, n
            if (.not. work(i)) cycle
            do j = 1, 3
                m = m + 1
                ids(m) = i
                sides(m) = j
                select case (j)
                case (1)
                    points(m) = lane(i)%xl
                case (2)
                    points(m) = lane(i)%xm
                case (3)
                    points(m) = lane(i)%xr
                end select
            end do
        end do
        if (m > 0) then
            if (present(f_indexed)) then
                call elementwise_gather_typed_indices(typed_indices, ids(:m), &
                    packed_typed_indices(:, :m))
                call f_indexed(points(:m), typed_args, &
                    packed_typed_indices(:, :m), values(:m))
            else if (present(f_typed)) then
                call elementwise_pack_typed_arguments(typed_args, &
                    typed_indices(:,ids(:m)), active_typed)
                call f_typed(points(:m), active_typed, values(:m))
            else if (present(f_args)) then
                call f_args(points(:m), args_data(:,ids(:m)), values(:m))
            else
                call f(points(:m), values(:m))
            end if
            do k = 1, m
                i = ids(k)
                select case (sides(k))
                case (1)
                    lane(i)%fl = values(k)
                case (2)
                    lane(i)%fm = values(k)
                case (3)
                    lane(i)%fr = values(k)
                end select
                res(i)%nfev = res(i)%nfev + 1
            end do
        end if
        do i = 1, n
            if (.not. work(i)) cycle
            if (ieee_is_finite(lane(i)%fl) .and. ieee_is_finite(lane(i)%fr)) then
                if (lane(i)%fl < lane(i)%fr) then
                    scratch = lane(i)%xl
                    lane(i)%xl = lane(i)%xr
                    lane(i)%xr = scratch
                    scratch = lane(i)%fl
                    lane(i)%fl = lane(i)%fr
                    lane(i)%fr = scratch
                    lane(i)%limit = lo(i)
                else
                    lane(i)%limit = hi(i)
                end if
            else
                lane(i)%limit = hi(i)
            end if
            lane(i)%limited = ieee_is_finite(lane(i)%limit)
            lane(i)%origin = lane(i)%xr
            lane(i)%step = lane(i)%xr - lane(i)%xm
            if (lane(i)%limited) then
                if (abs(lane(i)%limit / 2.0_dp - lane(i)%origin / 2.0_dp) > &
                        huge(1.0_dp) / 2.0_dp) then
                    res(i)%status = minimize_status_value_error
                    work(i) = .false.
                else
                    lane(i)%step = lane(i)%limit - lane(i)%origin
                end if
            end if
        end do

        do it = 0, niter
            if (.not. any(work)) exit
            m = 0
            do i = 1, n
                if (.not. work(i)) cycle
                if (.not. ieee_is_finite(lane(i)%fl) .or. &
                        .not. ieee_is_finite(lane(i)%fm) .or. &
                        .not. ieee_is_finite(lane(i)%fr)) then
                    res(i)%status = minimize_status_value_error
                    work(i) = .false.
                    cycle
                end if
                if ((lane(i)%fl >= lane(i)%fm .and. lane(i)%fr > lane(i)%fm) .or. &
                        (lane(i)%fl > lane(i)%fm .and. lane(i)%fr >= lane(i)%fm)) then
                    res(i)%status = minimize_status_success
                    res(i)%success = .true.
                    work(i) = .false.
                    cycle
                end if
                if (lane(i)%limited) then
                    if (lane(i)%xr == lane(i)%limit) then
                        res(i)%status = bracket_minimum_status_limit
                        work(i) = .false.
                        cycle
                    end if
                end if
                if (it == niter) then
                    res(i)%status = minimize_status_convergence_error
                    work(i) = .false.
                    cycle
                end if
                if (lane(i)%limited) then
                    lane(i)%step = lane(i)%step / grow(i)
                    next_x = lane(i)%limit - lane(i)%step
                    if (next_x == lane(i)%xr) next_x = lane(i)%limit
                else
                    if (abs(lane(i)%step) > huge(1.0_dp) / grow(i)) then
                        res(i)%status = minimize_status_value_error
                        work(i) = .false.
                        cycle
                    end if
                    lane(i)%step = lane(i)%step * grow(i)
                    if (abs(lane(i)%origin / 2.0_dp + lane(i)%step / 2.0_dp) > &
                            huge(1.0_dp) / 2.0_dp) then
                        res(i)%status = minimize_status_value_error
                        work(i) = .false.
                        cycle
                    end if
                    next_x = lane(i)%origin + lane(i)%step
                end if
                if (.not. ieee_is_finite(next_x)) then
                    res(i)%status = minimize_status_value_error
                    work(i) = .false.
                    cycle
                end if
                m = m + 1
                ids(m) = i
                points(m) = next_x
            end do
            if (m > 0) then
                if (present(f_indexed)) then
                    call elementwise_gather_typed_indices(typed_indices, ids(:m), &
                        packed_typed_indices(:, :m))
                    call f_indexed(points(:m), typed_args, &
                        packed_typed_indices(:, :m), values(:m))
                else if (present(f_typed)) then
                call elementwise_pack_typed_arguments(typed_args, &
                    typed_indices(:,ids(:m)), active_typed)
                call f_typed(points(:m), active_typed, values(:m))
                else if (present(f_args)) then
                call f_args(points(:m), args_data(:,ids(:m)), values(:m))
                else
                call f(points(:m), values(:m))
            end if
                do k = 1, m
                    i = ids(k)
                    lane(i)%xl = lane(i)%xm
                    lane(i)%fl = lane(i)%fm
                    lane(i)%xm = lane(i)%xr
                    lane(i)%fm = lane(i)%fr
                    lane(i)%xr = points(k)
                    lane(i)%fr = values(k)
                    res(i)%nfev = res(i)%nfev + 1
                    res(i)%nit = res(i)%nit + 1
                end do
            end if
        end do

        do i = 1, n
            if (res(i)%nfev == 0) cycle
            if (lane(i)%xl > lane(i)%xr) then
                scratch = lane(i)%xl
                lane(i)%xl = lane(i)%xr
                lane(i)%xr = scratch
                scratch = lane(i)%fl
                lane(i)%fl = lane(i)%fr
                lane(i)%fr = scratch
            end if
            res(i)%xl = lane(i)%xl
            res(i)%xm = lane(i)%xm
            res(i)%xr = lane(i)%xr
            res(i)%fl = lane(i)%fl
            res(i)%fm = lane(i)%fm
            res(i)%fr = lane(i)%fr
        end do
    end function bracket_minimum_vectorized_rank1

    function bracket_minimum_scalar(f, xm0, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Scalar adaptation of scipy.optimize.elementwise.bracket_minimum.
        ! The moving endpoint expands geometrically or approaches a finite limit.
        procedure(scalar_objective_function) :: f !! real scalar objective to bracket
        real(dp), intent(in) :: xm0 !! initial middle point
        real(dp), intent(in), optional :: xl0 !! initial left point (default xm0-0.5)
        real(dp), intent(in), optional :: xr0 !! initial right point (default xm0+0.5)
        real(dp), intent(in), optional :: xmin !! left search limit, inclusive
        real(dp), intent(in), optional :: xmax !! right search limit, inclusive
        real(dp), intent(in), optional :: factor !! expansion factor greater than 1 (default 2)
        integer, intent(in), optional :: maxiter !! maximum expansion steps (default 1000)
        type(bracket_minimum_result) :: res

        real(dp) :: xl, xm, xr, fl, fm, fr
        real(dp) :: left_limit, right_limit, limit, grow, step, origin
        real(dp) :: next_x, next_f, scratch, space
        integer :: maxiter_use, iter
        logical :: limited, can_expand

        left_limit = ieee_value(0.0_dp, ieee_negative_inf)
        right_limit = ieee_value(0.0_dp, ieee_positive_inf)
        if (present(xmin)) left_limit = xmin
        if (present(xmax)) right_limit = xmax
        grow = 2.0_dp
        if (present(factor)) grow = factor
        maxiter_use = 1000
        if (present(maxiter)) maxiter_use = maxiter

        if (.not. ieee_is_finite(xm0)) return
        if (ieee_is_nan(left_limit) .or. ieee_is_nan(right_limit)) return
        if (.not. ieee_is_finite(grow)) return
        if (grow <= 1.0_dp .or. maxiter_use < 0) return

        ! Scale bound distances before subtraction to avoid overflow when
        ! bounds have large opposite signs. SciPy uses 1/16 of each distance.
        xl = xm0 - 0.5_dp
        if (ieee_is_finite(left_limit)) then
            space = max(0.0_dp, xm0 / 16.0_dp - left_limit / 16.0_dp)
            xl = xm0 - min(0.5_dp, space)
        end if
        xr = xm0 + 0.5_dp
        if (ieee_is_finite(right_limit)) then
            space = max(0.0_dp, right_limit / 16.0_dp - xm0 / 16.0_dp)
            xr = xm0 + min(0.5_dp, space)
        end if
        if (present(xl0)) xl = xl0
        if (present(xr0)) xr = xr0
        if (.not. ieee_is_finite(xl) .or. .not. ieee_is_finite(xr)) return
        if (xl < left_limit .or. xl >= xm0 .or. xm0 >= xr .or. xr > right_limit) return

        xm = xm0
        fl = f(xl)
        fm = f(xm)
        fr = f(xr)
        res%nfev = 3

        ! Keep the working triple in downhill order. Leftward searches are
        ! reversed again before returning the public, increasing triple.
        if (ieee_is_finite(fl) .and. ieee_is_finite(fr)) then
            if (fl < fr) then
                scratch = xl
                xl = xr
                xr = scratch
                scratch = fl
                fl = fr
                fr = scratch
                limit = left_limit
            else
                limit = right_limit
            end if
        else
            limit = right_limit
        end if

        limited = ieee_is_finite(limit)
        origin = xr
        step = xr - xm
        can_expand = .true.
        if (limited) then
            ! The interval remaining before the limit shrinks by 1/factor.
            if (abs(limit / 2.0_dp - origin / 2.0_dp) > huge(1.0_dp) / 2.0_dp) then
                res%status = minimize_status_value_error
                can_expand = .false.
            else
                step = limit - origin
            end if
        end if

        if (can_expand) then
            do iter = 0, maxiter_use
            if (.not. ieee_is_finite(fl) .or. .not. ieee_is_finite(fm) .or. &
                    .not. ieee_is_finite(fr)) then
                res%status = minimize_status_value_error
                exit
            end if
            if ((fl >= fm .and. fr > fm) .or. (fl > fm .and. fr >= fm)) then
                res%status = minimize_status_success
                res%success = .true.
                exit
            end if
            if (limited) then
                if (xr == limit) then
                    res%status = bracket_minimum_status_limit
                    exit
                end if
            end if
            if (iter == maxiter_use) then
                res%status = minimize_status_convergence_error
                exit
            end if

            if (limited) then
                step = step / grow
                next_x = limit - step
                if (next_x == xr) next_x = limit
            else
                if (abs(step) > huge(1.0_dp) / grow) then
                    res%status = minimize_status_value_error
                    exit
                end if
                step = step * grow
                if (abs(origin / 2.0_dp + step / 2.0_dp) > huge(1.0_dp) / 2.0_dp) then
                    res%status = minimize_status_value_error
                    exit
                end if
                next_x = origin + step
            end if
            if (.not. ieee_is_finite(next_x)) then
                res%status = minimize_status_value_error
                exit
            end if

            next_f = f(next_x)
            res%nfev = res%nfev + 1
            res%nit = res%nit + 1
            xl = xm
            fl = fm
            xm = xr
            fm = fr
            xr = next_x
            fr = next_f
            end do
        end if

        if (xl > xr) then
            scratch = xl
            xl = xr
            xr = scratch
            scratch = fl
            fl = fr
            fr = scratch
        end if
        res%xl = xl
        res%xm = xm
        res%xr = xr
        res%fl = fl
        res%fm = fm
        res%fr = fr
    end function bracket_minimum_scalar

    function bracket(f, xa, xb, grow_limit, maxiter) result(res)
        ! Locate a downhill triple, returning SciPy's bracket data and work count.
        ! Bracketing can legitimately return in decreasing abscissa order.
        procedure(scalar_objective_function) :: f !! scalar objective function
        real(dp), intent(in), optional :: xa !! first trial point, default 0
        real(dp), intent(in), optional :: xb !! second trial point, default 1
        real(dp), intent(in), optional :: grow_limit !! maximum parabolic expansion factor
        integer, intent(in), optional :: maxiter !! maximum expansion iterations
        type(bracket_result) :: res

        type(bracket_state) :: state
        real(dp) :: xa_use
        real(dp) :: xb_use
        real(dp) :: limit_use
        integer :: maxiter_use

        xa_use = 0.0_dp
        xb_use = 1.0_dp
        limit_use = 110.0_dp
        maxiter_use = default_bracket_maxiter
        if (present(xa)) xa_use = xa
        if (present(xb)) xb_use = xb
        if (present(grow_limit)) limit_use = grow_limit
        if (present(maxiter)) maxiter_use = maxiter
        if (.not. ieee_is_finite(xa_use) .or. .not. ieee_is_finite(xb_use) .or. &
                .not. ieee_is_finite(limit_use) .or. &
                maxiter_use < 0) return

        call bracket_search(f, xa_use, xb_use, maxiter_use, state, limit_use)
        res%xa = state%xa
        res%xb = state%xb
        res%xc = state%xc
        res%fa = state%fa
        res%fb = state%fb
        res%fc = state%fc
        res%function_calls = state%function_calls
        res%iterations = state%iterations
        res%status = state%status
        res%success = state%valid
    end function bracket

    function minimize_scalar(f, method, bracket, bounds, tol, xtol, xatol, &
            maxiter) result(res)
        procedure(scalar_objective_function) :: f !! scalar objective function
        character(len=*), intent(in), optional :: method !! brent, bounded, or golden
        real(dp), intent(in), optional :: bracket(:) !! two- or three-point bracket data
        real(dp), intent(in), optional :: bounds(:) !! two finite bounds for bounded method
        real(dp), intent(in), optional :: tol !! SciPy-style generic x tolerance
        real(dp), intent(in), optional :: xtol !! relative x tolerance for brent or golden
        real(dp), intent(in), optional :: xatol !! absolute x tolerance for bounded method
        integer, intent(in), optional :: maxiter !! maximum solver iterations or calls
        type(minimize_result) :: res

        character(len=16) :: meth
        real(dp) :: tol_use

        call init_result(res, '')
        meth = ''
        if (present(method)) meth = lowercase(trim(method))
        if (len_trim(meth) == 0) then
            if (present(bounds)) then
                meth = 'bounded'
            else
                meth = 'brent'
            end if
        end if

        if (present(bounds)) then
            if (trim(meth) == 'brent' .or. trim(meth) == 'golden') return
        end if

        select case (trim(meth))
        case ('bounded')
            if (.not. present(bounds)) return
            if (size(bounds) /= 2) return
            tol_use = default_bounded_xatol
            if (present(tol)) tol_use = tol
            if (present(xatol)) tol_use = xatol
            res = bounded(f, bounds(1), bounds(2), xatol=tol_use, maxiter=maxiter)
        case ('brent')
            tol_use = default_brent_xtol
            if (present(tol)) tol_use = tol
            if (present(xtol)) tol_use = xtol
            if (present(bracket)) then
                res = brent(f, bracket=bracket, xtol=tol_use, maxiter=maxiter)
            else
                res = brent(f, xtol=tol_use, maxiter=maxiter)
            end if
        case ('golden')
            tol_use = default_golden_xtol
            if (present(tol)) tol_use = tol
            if (present(xtol)) tol_use = xtol
            if (present(bracket)) then
                res = golden(f, bracket=bracket, xtol=tol_use, maxiter=maxiter)
            else
                res = golden(f, xtol=tol_use, maxiter=maxiter)
            end if
        case default
            return
        end select
    end function minimize_scalar

    function bounded(f, lower, upper, xatol, maxiter) result(res)
        procedure(scalar_objective_function) :: f !! scalar objective function
        real(dp), intent(in) :: lower !! finite lower optimization bound
        real(dp), intent(in) :: upper !! finite upper optimization bound
        real(dp), intent(in), optional :: xatol !! absolute x tolerance, default 1e-5
        integer, intent(in), optional :: maxiter !! function-call cap, default 500
        type(minimize_result) :: res

        real(dp), parameter :: golden_mean = 0.5_dp * (3.0_dp - sqrt(5.0_dp))
        real(dp), parameter :: sqrt_eps = sqrt(2.2e-16_dp)
        real(dp) :: a
        real(dp) :: b
        real(dp) :: e
        real(dp) :: ffulc
        real(dp) :: fnfc
        real(dp) :: fu
        real(dp) :: fulc
        real(dp) :: fx
        real(dp) :: nfc
        real(dp) :: p
        real(dp) :: q
        real(dp) :: r
        real(dp) :: rat
        real(dp) :: si
        real(dp) :: tol1
        real(dp) :: tol2
        real(dp) :: tol_use
        real(dp) :: x
        real(dp) :: xf
        real(dp) :: xm
        integer :: maxfun
        integer :: num
        logical :: use_golden

        call init_result(res, 'bounded')
        tol_use = default_bounded_xatol
        if (present(xatol)) tol_use = xatol
        maxfun = default_bounded_maxiter
        if (present(maxiter)) maxfun = maxiter
        if (.not. ieee_is_finite(lower) .or. .not. ieee_is_finite(upper)) return
        if (lower > upper) return
        if (.not. ieee_is_finite(tol_use) .or. tol_use < 0.0_dp) return

        a = lower
        b = upper
        fulc = a + golden_mean * (b - a)
        nfc = fulc
        xf = fulc
        rat = 0.0_dp
        e = 0.0_dp
        x = xf
        fx = f(x)
        num = 1
        ffulc = fx
        fnfc = fx
        fu = fx
        if (ieee_is_nan(fx)) then
            res%x = xf
            res%fun = fx
            res%function_calls = num
            res%iterations = num
            call fail_result(res, minimize_status_value_error)
            return
        end if

        xm = 0.5_dp * (a + b)
        tol1 = sqrt_eps * abs(xf) + tol_use / 3.0_dp
        tol2 = 2.0_dp * tol1

        do while (abs(xf - xm) > (tol2 - 0.5_dp * (b - a)))
            use_golden = .true.
            if (abs(e) > tol1) then
                use_golden = .false.
                r = (xf - nfc) * (fx - ffulc)
                q = (xf - fulc) * (fx - fnfc)
                p = (xf - fulc) * q - (xf - nfc) * r
                q = 2.0_dp * (q - r)
                if (q > 0.0_dp) p = -p
                q = abs(q)
                r = e
                e = rat
                if (abs(p) < abs(0.5_dp * q * r) .and. &
                        p > q * (a - xf) .and. p < q * (b - xf)) then
                    rat = p / q
                    x = xf + rat
                    if ((x - a) < tol2 .or. (b - x) < tol2) then
                        si = signed_one(xm - xf)
                        rat = tol1 * si
                    end if
                else
                    use_golden = .true.
                end if
            end if

            if (use_golden) then
                if (xf >= xm) then
                    e = a - xf
                else
                    e = b - xf
                end if
                rat = golden_mean * e
            end if

            si = signed_one(rat)
            x = xf + si * max(abs(rat), tol1)
            fu = f(x)
            num = num + 1
            if (ieee_is_nan(fu)) then
                res%x = x
                res%fun = fu
                res%function_calls = num
                res%iterations = num
                call fail_result(res, minimize_status_value_error)
                return
            end if

            if (fu <= fx) then
                if (x >= xf) then
                    a = xf
                else
                    b = xf
                end if
                fulc = nfc
                ffulc = fnfc
                nfc = xf
                fnfc = fx
                xf = x
                fx = fu
            else
                if (x < xf) then
                    a = x
                else
                    b = x
                end if
                if (fu <= fnfc .or. same_real(nfc, xf)) then
                    fulc = nfc
                    ffulc = fnfc
                    nfc = x
                    fnfc = fu
                else if (fu <= ffulc .or. same_real(fulc, xf) .or. same_real(fulc, nfc)) then
                    fulc = x
                    ffulc = fu
                end if
            end if

            xm = 0.5_dp * (a + b)
            tol1 = sqrt_eps * abs(xf) + tol_use / 3.0_dp
            tol2 = 2.0_dp * tol1
            if (num >= maxfun) exit
        end do

        res%x = xf
        res%fun = fx
        res%function_calls = num
        res%iterations = num
        if (num >= maxfun .and. &
                abs(xf - xm) > (tol2 - 0.5_dp * (b - a))) then
            call fail_result(res, minimize_status_convergence_error)
        else
            call succeed_result(res)
        end if
    end function bounded

    function brent(f, bracket, xtol, maxiter) result(res)
        procedure(scalar_objective_function) :: f !! scalar objective function
        real(dp), intent(in), optional :: bracket(:) !! optional two- or three-point bracket
        real(dp), intent(in), optional :: xtol !! relative x tolerance, default 1.48e-8
        integer, intent(in), optional :: maxiter !! iteration cap, default 500
        type(minimize_result) :: res

        real(dp), parameter :: cg = 0.3819660_dp
        real(dp), parameter :: mintol = 1.0e-11_dp
        type(bracket_state) :: br
        real(dp) :: a
        real(dp) :: b
        real(dp) :: deltax
        real(dp) :: dx_temp
        real(dp) :: fu
        real(dp) :: fv
        real(dp) :: fw
        real(dp) :: fx
        real(dp) :: p
        real(dp) :: rat
        real(dp) :: tmp1
        real(dp) :: tmp2
        real(dp) :: tol1
        real(dp) :: tol2
        real(dp) :: tol_use
        real(dp) :: u
        real(dp) :: v
        real(dp) :: w
        real(dp) :: x
        real(dp) :: xmid
        integer :: iter
        integer :: mit

        call init_result(res, 'brent')
        tol_use = default_brent_xtol
        if (present(xtol)) tol_use = xtol
        mit = default_brent_maxiter
        if (present(maxiter)) mit = maxiter
        if (.not. ieee_is_finite(tol_use) .or. tol_use < 0.0_dp) return

        if (present(bracket)) then
            call prepare_bracket(f, bracket, br)
        else
            call bracket_search(f, 0.0_dp, 1.0_dp, default_bracket_maxiter, br)
        end if
        if (.not. br%valid) then
            res%function_calls = br%function_calls
            if (br%status == minimize_status_value_error) then
                res%x = ieee_value(0.0_dp, ieee_quiet_nan)
                res%fun = ieee_value(0.0_dp, ieee_quiet_nan)
            else
                res%x = br%xb
                res%fun = br%fb
            end if
            call fail_result(res, br%status)
            return
        end if

        x = br%xb
        w = x
        v = x
        fw = br%fb
        fv = fw
        fx = fw
        a = min(br%xa, br%xc)
        b = max(br%xa, br%xc)
        deltax = 0.0_dp
        rat = 0.0_dp
        iter = 0

        do while (iter < mit)
            tol1 = tol_use * abs(x) + mintol
            tol2 = 2.0_dp * tol1
            xmid = 0.5_dp * (a + b)
            if (abs(x - xmid) < (tol2 - 0.5_dp * (b - a))) exit

            if (abs(deltax) <= tol1) then
                if (x >= xmid) then
                    deltax = a - x
                else
                    deltax = b - x
                end if
                rat = cg * deltax
            else
                tmp1 = (x - w) * (fx - fv)
                tmp2 = (x - v) * (fx - fw)
                p = (x - v) * tmp2 - (x - w) * tmp1
                tmp2 = 2.0_dp * (tmp2 - tmp1)
                if (tmp2 > 0.0_dp) p = -p
                tmp2 = abs(tmp2)
                dx_temp = deltax
                deltax = rat
                if (p > tmp2 * (a - x) .and. p < tmp2 * (b - x) .and. &
                        abs(p) < abs(0.5_dp * tmp2 * dx_temp)) then
                    rat = p / tmp2
                    u = x + rat
                    if ((u - a) < tol2 .or. (b - u) < tol2) then
                        rat = sign(tol1, xmid - x)
                        if (same_real(xmid, x)) rat = tol1
                    end if
                else
                    if (x >= xmid) then
                        deltax = a - x
                    else
                        deltax = b - x
                    end if
                    rat = cg * deltax
                end if
            end if

            if (abs(rat) < tol1) then
                if (rat >= 0.0_dp) then
                    u = x + tol1
                else
                    u = x - tol1
                end if
            else
                u = x + rat
            end if
            fu = f(u)
            br%function_calls = br%function_calls + 1
            if (ieee_is_nan(fu)) then
                res%x = u
                res%fun = fu
                res%iterations = iter
                res%function_calls = br%function_calls
                call fail_result(res, minimize_status_value_error)
                return
            end if

            if (fu > fx) then
                if (u < x) then
                    a = u
                else
                    b = u
                end if
                if (fu <= fw .or. same_real(w, x)) then
                    v = w
                    w = u
                    fv = fw
                    fw = fu
                else if (fu <= fv .or. same_real(v, x) .or. same_real(v, w)) then
                    v = u
                    fv = fu
                end if
            else
                if (u >= x) then
                    a = x
                else
                    b = x
                end if
                v = w
                w = x
                x = u
                fv = fw
                fw = fx
                fx = fu
            end if
            iter = iter + 1
        end do

        res%x = x
        res%fun = fx
        res%iterations = iter
        res%function_calls = br%function_calls
        if (ieee_is_nan(x) .or. ieee_is_nan(fx)) then
            call fail_result(res, minimize_status_value_error)
        else if (iter >= mit) then
            call fail_result(res, minimize_status_convergence_error)
        else
            call succeed_result(res)
        end if
    end function brent

    function golden(f, bracket, xtol, maxiter) result(res)
        procedure(scalar_objective_function) :: f !! scalar objective function
        real(dp), intent(in), optional :: bracket(:) !! optional two- or three-point bracket
        real(dp), intent(in), optional :: xtol !! relative x tolerance, default epsilon
        integer, intent(in), optional :: maxiter !! iteration cap, default 5000
        type(minimize_result) :: res

        real(dp), parameter :: gr = 0.61803399_dp
        real(dp), parameter :: gc = 1.0_dp - gr
        type(bracket_state) :: br
        real(dp) :: f1
        real(dp) :: f2
        real(dp) :: tol_use
        real(dp) :: x0
        real(dp) :: x1
        real(dp) :: x2
        real(dp) :: x3
        integer :: i
        integer :: mit
        integer :: nit

        call init_result(res, 'golden')
        tol_use = default_golden_xtol
        if (present(xtol)) tol_use = xtol
        mit = default_golden_maxiter
        if (present(maxiter)) mit = maxiter
        if (.not. ieee_is_finite(tol_use) .or. tol_use < 0.0_dp) return

        if (present(bracket)) then
            call prepare_bracket(f, bracket, br)
        else
            call bracket_search(f, 0.0_dp, 1.0_dp, default_bracket_maxiter, br)
        end if
        if (.not. br%valid) then
            res%function_calls = br%function_calls
            if (br%status == minimize_status_value_error) then
                res%x = ieee_value(0.0_dp, ieee_quiet_nan)
                res%fun = ieee_value(0.0_dp, ieee_quiet_nan)
            else
                res%x = br%xb
                res%fun = br%fb
            end if
            call fail_result(res, br%status)
            return
        end if

        x3 = br%xc
        x0 = br%xa
        if (abs(br%xc - br%xb) > abs(br%xb - br%xa)) then
            x1 = br%xb
            x2 = br%xb + gc * (br%xc - br%xb)
        else
            x2 = br%xb
            x1 = br%xb - gc * (br%xb - br%xa)
        end if
        f1 = f(x1)
        f2 = f(x2)
        br%function_calls = br%function_calls + 2
        if (ieee_is_nan(f1) .or. ieee_is_nan(f2)) then
            res%x = merge(x1, x2, ieee_is_nan(f1))
            res%fun = merge(f1, f2, ieee_is_nan(f1))
            res%function_calls = br%function_calls
            call fail_result(res, minimize_status_value_error)
            return
        end if

        nit = 0
        do i = 1, mit
            if (abs(x3 - x0) <= tol_use * (abs(x1) + abs(x2))) exit
            if (f2 < f1) then
                x0 = x1
                x1 = x2
                x2 = gr * x1 + gc * x3
                f1 = f2
                f2 = f(x2)
                if (ieee_is_nan(f2)) then
                    res%x = x2
                    res%fun = f2
                    res%iterations = nit
                    res%function_calls = br%function_calls + 1
                    call fail_result(res, minimize_status_value_error)
                    return
                end if
            else
                x3 = x2
                x2 = x1
                x1 = gr * x2 + gc * x0
                f2 = f1
                f1 = f(x1)
                if (ieee_is_nan(f1)) then
                    res%x = x1
                    res%fun = f1
                    res%iterations = nit
                    res%function_calls = br%function_calls + 1
                    call fail_result(res, minimize_status_value_error)
                    return
                end if
            end if
            br%function_calls = br%function_calls + 1
            nit = nit + 1
        end do

        if (f1 < f2) then
            res%x = x1
            res%fun = f1
        else
            res%x = x2
            res%fun = f2
        end if
        res%iterations = nit
        res%function_calls = br%function_calls
        if (nit >= mit) then
            call fail_result(res, minimize_status_convergence_error)
        else
            call succeed_result(res)
        end if
    end function golden

    subroutine prepare_bracket(f, points, br)
        procedure(scalar_objective_function) :: f !! scalar objective function
        real(dp), intent(in) :: points(:) !! two initial points or a valid three-point bracket
        type(bracket_state), intent(out) :: br !! prepared bracket and objective values

        real(dp) :: xa
        real(dp) :: xb
        real(dp) :: xc
        real(dp) :: tmp

        br = bracket_state()
        if (size(points) == 2) then
            call bracket_search(f, points(1), points(2), default_bracket_maxiter, br)
            return
        end if
        if (size(points) /= 3) then
            br%status = minimize_status_invalid_input
            return
        end if

        xa = points(1)
        xb = points(2)
        xc = points(3)
        if (xa > xc) then
            tmp = xa
            xa = xc
            xc = tmp
        end if
        br%xa = xa
        br%xb = xb
        br%xc = xc
        if (.not. (xa < xb .and. xb < xc)) then
            br%status = minimize_status_invalid_input
            return
        end if
        if (.not. ieee_is_finite(xa) .or. .not. ieee_is_finite(xb) .or. &
                .not. ieee_is_finite(xc)) then
            br%status = minimize_status_invalid_input
            return
        end if

        br%fa = f(xa)
        br%fb = f(xb)
        br%fc = f(xc)
        br%function_calls = 3
        if (ieee_is_nan(br%fa) .or. ieee_is_nan(br%fb) .or. ieee_is_nan(br%fc)) then
            br%status = minimize_status_value_error
            return
        end if
        if (.not. (br%fb < br%fa .and. br%fb < br%fc)) then
            br%status = minimize_status_bracket_error
            return
        end if
        br%valid = .true.
        br%status = minimize_status_success
    end subroutine prepare_bracket

    subroutine bracket_search(f, xa_in, xb_in, maxiter, br, limit)
        procedure(scalar_objective_function) :: f !! scalar objective function
        real(dp), intent(in) :: xa_in !! first initial point
        real(dp), intent(in) :: xb_in !! second initial point
        integer, intent(in) :: maxiter !! maximum expansion iterations
        type(bracket_state), intent(out) :: br !! resulting three-point bracket
        real(dp), intent(in), optional :: limit !! maximum growth factor, default 110

        real(dp), parameter :: gold = 1.618034_dp
        real(dp), parameter :: verysmall = 1.0e-21_dp
        real(dp) :: grow_limit
        real(dp) :: denom
        real(dp) :: fa
        real(dp) :: fb
        real(dp) :: fc
        real(dp) :: fw
        real(dp) :: tmp
        real(dp) :: tmp1
        real(dp) :: tmp2
        real(dp) :: val
        real(dp) :: w
        real(dp) :: wlim
        real(dp) :: xa
        real(dp) :: xb
        real(dp) :: xc
        integer :: iter
        logical :: cond1
        logical :: cond2
        logical :: cond3

        br = bracket_state()
        grow_limit = 110.0_dp
        if (present(limit)) grow_limit = limit
        if (.not. ieee_is_finite(xa_in) .or. .not. ieee_is_finite(xb_in)) then
            br%status = minimize_status_invalid_input
            return
        end if
        xa = xa_in
        xb = xb_in
        fa = f(xa)
        fb = f(xb)
        br%function_calls = 2
        ! Quiet NaNs must not enter an ordered comparison when IEEE invalid
        ! exceptions are trapped. Keep the normal three-call failure accounting.
        if (.not. ieee_is_nan(fa) .and. .not. ieee_is_nan(fb)) then
            if (fa < fb) then
                tmp = xa
                xa = xb
                xb = tmp
                tmp = fa
                fa = fb
                fb = tmp
            end if
        end if
        xc = xb + gold * (xb - xa)
        fc = f(xc)
        br%function_calls = 3
        if (ieee_is_nan(fa) .or. ieee_is_nan(fb) .or. ieee_is_nan(fc)) then
            br%xa = xa
            br%xb = xb
            br%xc = xc
            br%fa = fa
            br%fb = fb
            br%fc = fc
            br%status = minimize_status_value_error
            return
        end if

        iter = 0
        do while (fc < fb)
            tmp1 = (xb - xa) * (fb - fc)
            tmp2 = (xb - xc) * (fb - fa)
            val = tmp2 - tmp1
            if (abs(val) < verysmall) then
                denom = 2.0_dp * verysmall
            else
                denom = 2.0_dp * val
            end if
            w = xb - ((xb - xc) * tmp2 - (xb - xa) * tmp1) / denom
            wlim = xb + grow_limit * (xc - xb)
            if (iter > maxiter) then
                br%status = minimize_status_convergence_error
                exit
            end if
            iter = iter + 1

            if ((w - xc) * (xb - w) > 0.0_dp) then
                fw = f(w)
                br%function_calls = br%function_calls + 1
                if (ieee_is_nan(fw)) then
                    br%status = minimize_status_value_error
                    exit
                end if
                if (fw < fc) then
                    xa = xb
                    xb = w
                    fa = fb
                    fb = fw
                    exit
                else if (fw > fb) then
                    xc = w
                    fc = fw
                    exit
                end if
                w = xc + gold * (xc - xb)
                fw = f(w)
                br%function_calls = br%function_calls + 1
            else if ((w - wlim) * (wlim - xc) >= 0.0_dp) then
                w = wlim
                fw = f(w)
                br%function_calls = br%function_calls + 1
            else if ((w - wlim) * (xc - w) > 0.0_dp) then
                fw = f(w)
                br%function_calls = br%function_calls + 1
                if (ieee_is_nan(fw)) then
                    br%status = minimize_status_value_error
                    exit
                end if
                if (fw < fc) then
                    xb = xc
                    xc = w
                    w = xc + gold * (xc - xb)
                    fb = fc
                    fc = fw
                    fw = f(w)
                    br%function_calls = br%function_calls + 1
                end if
            else
                w = xc + gold * (xc - xb)
                fw = f(w)
                br%function_calls = br%function_calls + 1
            end if
            if (ieee_is_nan(fw)) then
                br%status = minimize_status_value_error
                exit
            end if
            xa = xb
            xb = xc
            xc = w
            fa = fb
            fb = fc
            fc = fw
        end do

        br%iterations = iter

        br%xa = xa
        br%xb = xb
        br%xc = xc
        br%fa = fa
        br%fb = fb
        br%fc = fc
        if (br%status == minimize_status_value_error .or. &
                br%status == minimize_status_convergence_error) return

        cond1 = (fb < fc .and. fb <= fa) .or. (fb < fa .and. fb <= fc)
        cond2 = (xa < xb .and. xb < xc) .or. (xc < xb .and. xb < xa)
        cond3 = ieee_is_finite(xa) .and. ieee_is_finite(xb) .and. ieee_is_finite(xc)
        if (cond1 .and. cond2 .and. cond3) then
            br%valid = .true.
            br%status = minimize_status_success
        else
            br%valid = .false.
            br%status = minimize_status_bracket_error
        end if
    end subroutine bracket_search

    subroutine init_result(res, method)
        type(minimize_result), intent(out) :: res !! result object to initialize
        character(len=*), intent(in) :: method !! solver method label

        res%x = ieee_value(0.0_dp, ieee_quiet_nan)
        res%fun = ieee_value(0.0_dp, ieee_quiet_nan)
        res%iterations = 0
        res%function_calls = 0
        res%status = minimize_status_invalid_input
        res%success = .false.
        res%flag = 'invalid input'
        res%method = method
    end subroutine init_result

    subroutine succeed_result(res)
        type(minimize_result), intent(inout) :: res !! result object to mark successful

        res%status = minimize_status_success
        res%success = .true.
        res%flag = 'converged'
    end subroutine succeed_result

    subroutine fail_result(res, status)
        type(minimize_result), intent(inout) :: res !! result object to mark unsuccessful
        integer, intent(in) :: status !! failure status code

        res%status = status
        res%success = .false.
        select case (status)
        case (minimize_status_bracket_error)
            res%flag = 'bracket error'
        case (minimize_status_convergence_error)
            res%flag = 'convergence error'
        case (minimize_status_value_error)
            res%flag = 'value error'
        case default
            res%flag = 'invalid input'
        end select
    end subroutine fail_result

    pure function same_real(a, b) result(equal)
        real(dp), intent(in) :: a !! first finite value to compare exactly
        real(dp), intent(in) :: b !! second finite value to compare exactly
        logical :: equal

        equal = a <= b .and. a >= b
    end function same_real

    pure function signed_one(x) result(s)
        real(dp), intent(in) :: x !! value whose sign selects plus or minus one
        real(dp) :: s

        if (x < 0.0_dp) then
            s = -1.0_dp
        else
            s = 1.0_dp
        end if
    end function signed_one

    pure function lowercase(text) result(out)
        character(len=*), intent(in) :: text !! text to convert to lowercase ASCII
        character(len=len(text)) :: out

        integer :: code
        integer :: i

        out = text
        do i = 1, len(text)
            code = iachar(text(i:i))
            if (code >= iachar('A') .and. code <= iachar('Z')) then
                out(i:i) = achar(code + iachar('a') - iachar('A'))
            end if
        end do
    end function lowercase

end module scifort_optimize_minimize_scalar
