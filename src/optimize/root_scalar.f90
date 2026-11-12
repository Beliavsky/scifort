! SPDX-License-Identifier: BSD-3-Clause
! Copyright (c) 2026 SciFort contributors
! Copyright (c) 2021 Jacob Williams
! Copyright (c) 2001-2002 Enthought, Inc. 2003, SciPy Developers.
!
! Scalar root finding for the scipy.optimize translation.
!
! The bisection, Brent brentq/brenth, and Ridders kernels are adapted from
! roots-fortran 1.0.0 by Jacob Williams (BSD-3-Clause). roots-fortran's
! brentq/brenth routines in turn document their relationship to SciPy's
! implementations. The root_scalar dispatch/result semantics, Newton/secant/
! Halley iteration, and TOMS748 organization are adapted from SciPy 1.17.0
! (BSD-3-Clause). SciFort-specific validation and IEEE-trap-safe sign tests
! replace multiplication-based sign checks.
! Scalar root bracketing follows SciPy 1.17.0 scipy.optimize._bracket
! (BSD-3-Clause). Generic rank-one through rank-15 batches use scalar kernels.
! The vectorized refinement and bracketing routines use packed array callbacks.

module scifort_optimize_root_scalar
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_negative_inf, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
    use scifort_differentiate, only : derivative, derivative_result
    use scifort_kinds, only : dp
    use scifort_optimize_elementwise_batch, only : elementwise_broadcast, &
        elementwise_broadcast_optional, elementwise_flatten, &
        elementwise_broadcast_arguments, &
        elementwise_argument, elementwise_broadcast_argument_list, &
        elementwise_typed_argument, elementwise_broadcast_typed_argument_list, &
        elementwise_pack_typed_arguments, elementwise_gather_typed_indices
    implicit none
    private

    integer, parameter, public :: root_status_success = 0
    integer, parameter, public :: root_status_sign_error = -1
    integer, parameter, public :: root_status_convergence_error = -2
    integer, parameter, public :: root_status_value_error = -3
    integer, parameter, public :: root_status_invalid_input = -5
    integer, parameter, public :: bracket_root_status_limit = -1
    integer, parameter, public :: find_root_status_callback = -4

    real(dp), parameter :: default_xtol = 2.0e-12_dp
    real(dp), parameter :: default_rtol = 4.0_dp * epsilon(1.0_dp)
    integer, parameter :: default_bracket_maxiter = 100
    real(dp), parameter :: default_newton_tol = 1.48e-8_dp
    integer, parameter :: default_newton_maxiter = 50

    type, public :: root_result
        real(dp) :: root = 0.0_dp !! approximate root
        integer :: iterations = 0 !! iterations performed
        integer :: function_calls = 0 !! objective/derivative evaluations
        integer :: status = root_status_invalid_input !! numerical termination code
        logical :: converged = .false. !! true when the requested stopping rule was met
        character(len=32) :: flag = 'invalid input' !! human-readable termination flag
        character(len=16) :: method = '' !! solver method actually used
    end type root_result

    type, public :: bracket_root_result
        real(dp) :: xl = 0.0_dp !! left endpoint of the discovered or searched interval
        real(dp) :: xr = 0.0_dp !! right endpoint of the discovered or searched interval
        real(dp) :: fl = 0.0_dp !! objective value at xl
        real(dp) :: fr = 0.0_dp !! objective value at xr
        integer :: nfev = 0 !! number of function evaluations
        integer :: nit = 0 !! number of expansion iterations
        integer :: status = root_status_invalid_input !! 0, -1, -2, -3, or -5
        logical :: success = .false. !! true if an endpoint vanishes or signs differ
    end type bracket_root_result

    type, public :: find_root_result
        real(dp) :: x = 0.0_dp !! best current root estimate
        real(dp) :: f_x = 0.0_dp !! function value at x
        real(dp) :: xl = 0.0_dp !! left endpoint of current bracket
        real(dp) :: xr = 0.0_dp !! right endpoint of current bracket
        real(dp) :: fl = 0.0_dp !! function value at xl
        real(dp) :: fr = 0.0_dp !! function value at xr
        integer :: nit = 0 !! refinement iterations
        integer :: nfev = 0 !! scalar function evaluations
        integer :: status = root_status_invalid_input !! 0/-1/-2/-3/-4 or -5 invalid input
        logical :: success = .false. !! true when convergence criterion is met
    end type find_root_result

    ! The two independently expanding sides of one vectorized root bracket.
    type :: vector_bracket_root_lane
        real(dp) :: x(2) = 0.0_dp
        real(dp) :: fx(2) = 0.0_dp
        real(dp) :: xlast(2) = 0.0_dp
        real(dp) :: flast(2) = 0.0_dp
        real(dp) :: anchor(2) = 0.0_dp
        real(dp) :: bound(2) = 0.0_dp
        real(dp) :: distance(2) = 0.0_dp
        logical :: bounded(2) = .false.
        logical :: active(2) = .false.
        logical :: found(2) = .false.
        logical :: nonfinite = .false.
    end type vector_bracket_root_lane

    abstract interface
        function scalar_root_function(x) result(fx)
            import dp
            implicit none
            real(dp), intent(in) :: x !! scalar abscissa
            real(dp) :: fx !! scalar function value
        end function scalar_root_function
        subroutine find_root_callback(state, stop)
            import find_root_result
            implicit none
            type(find_root_result), intent(in) :: state !! current root and bracket
            logical, intent(out) :: stop !! true to terminate refinement
        end subroutine find_root_callback
        subroutine elementwise_root_vector_function(x, fx)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! active abscissae, tightly packed
            real(dp), intent(out) :: fx(:) !! function values, same size as x
        end subroutine elementwise_root_vector_function
        subroutine elementwise_root_args_function(x, args, fx)
            import dp
            implicit none
            real(dp), intent(in) :: x(:) !! packed active objective points
            real(dp), intent(in) :: args(:,:) !! packed parameters (nargs, size(x))
            real(dp), intent(out) :: fx(:) !! values for each packed point
        end subroutine elementwise_root_args_function
        subroutine elementwise_root_typed_args_function(x, args, fx)
            import dp, elementwise_typed_argument
            implicit none
            real(dp), intent(in) :: x(:) !! packed active abscissae
            type(elementwise_typed_argument), intent(in) :: args(:) !! typed values, per active lane
            real(dp), intent(out) :: fx(:) !! returned objective values
        end subroutine elementwise_root_typed_args_function
        subroutine elementwise_root_indexed_args_function(x, args, indices, fx)
            import dp, elementwise_typed_argument
            implicit none
            real(dp), intent(in) :: x(:) !! active abscissae in packed callback order
            type(elementwise_typed_argument), intent(in) :: args(:) !! original typed input snapshots
            integer, intent(in) :: indices(:,:) !! (nargs,nactive) offsets into args(i)%values
            real(dp), intent(out) :: fx(:) !! objective at each active abscissa
        end subroutine elementwise_root_indexed_args_function
    end interface

    public :: find_root_vectorized
    public :: find_root_vectorized_args_into
    public :: find_root_vectorized_indexed_args_into
    public :: bracket_root_vectorized_args_into
    public :: bracket_root_vectorized_indexed_args_into
    public :: elementwise_root_args_function
    public :: elementwise_root_typed_args_function
    public :: elementwise_root_indexed_args_function
    public :: find_root_vectorized_into
    public :: bracket_root_vectorized
    public :: bracket_root_vectorized_into
    public :: elementwise_root_vector_function
    public :: scalar_root_function
    public :: root_scalar
    public :: bracket_root
    public :: find_root
    public :: find_root_callback
    public :: bisect
    public :: brentq
    public :: brenth
    public :: ridder
    public :: toms748
    public :: newton

    ! Batch overloads through rank 15 preserve the scalar API. Each result
    ! reports its own work counts and status; callbacks remain scalar.
    ! The packed callback always receives a rank-one array of active points.
    ! These generics accept rank-one through rank-four problem grids.
    interface find_root_vectorized_args_into
        module procedure find_root_vectorized_args_into_stacked
        module procedure find_root_vectorized_args_into_list
        module procedure find_root_vectorized_args_into_typed
    end interface find_root_vectorized_args_into

    interface bracket_root_vectorized_args_into
        module procedure bracket_root_vectorized_args_into_stacked
        module procedure bracket_root_vectorized_args_into_list
        module procedure bracket_root_vectorized_args_into_typed
    end interface bracket_root_vectorized_args_into

    interface find_root_vectorized
        module procedure find_root_vectorized_rank1
        module procedure find_root_vectorized_rank2
        module procedure find_root_vectorized_rank3
        module procedure find_root_vectorized_rank4
    end interface find_root_vectorized

    interface bracket_root_vectorized
        module procedure bracket_root_vectorized_rank1
        module procedure bracket_root_vectorized_rank2
        module procedure bracket_root_vectorized_rank3
        module procedure bracket_root_vectorized_rank4
    end interface bracket_root_vectorized

    interface find_root
        module procedure find_root_scalar
        module procedure find_root_vector
        module procedure find_root_matrix
        module procedure find_root_tensor
        module procedure find_root_rank14
        module procedure find_root_rank13
        module procedure find_root_rank12
        module procedure find_root_rank11
        module procedure find_root_rank10
        module procedure find_root_rank9
        module procedure find_root_rank8
        module procedure find_root_rank7
        module procedure find_root_rank6
        module procedure find_root_rank5
        module procedure find_root_rank4
    end interface find_root

    interface bracket_root
        module procedure bracket_root_scalar
        module procedure bracket_root_vector
        module procedure bracket_root_matrix
        module procedure bracket_root_tensor
        module procedure bracket_root_rank15
        module procedure bracket_root_rank14
        module procedure bracket_root_rank13
        module procedure bracket_root_rank12
        module procedure bracket_root_rank11
        module procedure bracket_root_rank10
        module procedure bracket_root_rank9
        module procedure bracket_root_rank8
        module procedure bracket_root_rank7
        module procedure bracket_root_rank6
        module procedure bracket_root_rank5
        module procedure bracket_root_rank4
    end interface bracket_root

contains

    subroutine find_root_vectorized_args_into_stacked(f, bracket, args, res, &
            xatol, xrtol, fatol, frtol, maxiter, callback)
        ! Explicit args entry: retain the legacy (x, fx) callback ABI.
        procedure(elementwise_root_args_function) :: f !! receives x, packed args, fx
        real(dp), intent(in) :: bracket(..) !! starting problem grid
        real(dp), intent(in) :: args(..) !! stacked parameters, first axis = nargs
        type(find_root_result), intent(inout) :: res(..) !! initialized grid results
        real(dp), intent(in), optional :: xatol(..) !! broadcast tolerance
        real(dp), intent(in), optional :: xrtol(..) !! broadcast tolerance
        real(dp), intent(in), optional :: fatol(..) !! broadcast tolerance
        real(dp), intent(in), optional :: frtol(..) !! broadcast tolerance
        integer, intent(in), optional :: maxiter !! iteration limit
        procedure(find_root_callback), optional :: callback !! lane stop callback
        call find_root_vectorized_into(bracket=bracket, res=res, &
            xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, &
            callback=callback, args=args, f_args=f)
    end subroutine find_root_vectorized_args_into_stacked

    subroutine bracket_root_vectorized_args_into_stacked(f, xl0, args, res, &
            xr0, xmin, xmax, factor, maxiter)
        ! Explicit args entry: retain the legacy (x, fx) callback ABI.
        procedure(elementwise_root_args_function) :: f !! receives x, packed args, fx
        real(dp), intent(in) :: xl0(..) !! starting problem grid
        real(dp), intent(in) :: args(..) !! stacked parameters, first axis = nargs
        type(bracket_root_result), intent(inout) :: res(..) !! initialized grid results
        real(dp), intent(in), optional :: xr0(..) !! broadcast right endpoint
        real(dp), intent(in), optional :: xmin(..) !! broadcast lower limit
        real(dp), intent(in), optional :: xmax(..) !! broadcast upper limit
        real(dp), intent(in), optional :: factor(..) !! broadcast expansion factor
        integer, intent(in), optional :: maxiter !! iteration limit
        call bracket_root_vectorized_into(xl0=xl0, res=res, &
            xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, &
            args=args, f_args=f)
    end subroutine bracket_root_vectorized_args_into_stacked

    subroutine find_root_vectorized_args_into_list(f, bracket, args, res, &
            xatol, xrtol, fatol, frtol, maxiter, callback)
        ! Alternative to a stacked args array: every parameter may have its
        ! own scalar/array rank and shape.  The callback ABI is unchanged.
        procedure(elementwise_root_args_function) :: f !! receives x, packed args, fx
        real(dp), intent(in) :: bracket(..) !! starting problem grid
        type(elementwise_argument), intent(in) :: args(:) !! independent parameters
        type(find_root_result), intent(inout) :: res(..) !! result grid
        real(dp), intent(in), optional :: xatol(..) !! broadcast tolerance
        real(dp), intent(in), optional :: xrtol(..) !! broadcast tolerance
        real(dp), intent(in), optional :: fatol(..) !! broadcast tolerance
        real(dp), intent(in), optional :: frtol(..) !! broadcast tolerance
        integer, intent(in), optional :: maxiter !! iteration limit
        procedure(find_root_callback), optional :: callback !! lane stop callback
        call find_root_vectorized_into(bracket=bracket, res=res, &
            xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, &
            callback=callback, args_list=args, f_args=f)
    end subroutine find_root_vectorized_args_into_list

    subroutine bracket_root_vectorized_args_into_list(f, xl0, args, res, &
            xr0, xmin, xmax, factor, maxiter)
        procedure(elementwise_root_args_function) :: f !! receives x, packed args, fx
        real(dp), intent(in) :: xl0(..) !! left starting grid
        type(elementwise_argument), intent(in) :: args(:) !! independent parameters
        type(bracket_root_result), intent(inout) :: res(..) !! result grid
        real(dp), intent(in), optional :: xr0(..) !! right starting points
        real(dp), intent(in), optional :: xmin(..) !! lower bounds
        real(dp), intent(in), optional :: xmax(..) !! upper bounds
        real(dp), intent(in), optional :: factor(..) !! expansion factor
        integer, intent(in), optional :: maxiter !! iteration limit
        call bracket_root_vectorized_into(xl0=xl0, res=res, &
            xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, &
            args_list=args, f_args=f)
    end subroutine bracket_root_vectorized_args_into_list



    subroutine find_root_vectorized_args_into_typed(f, bracket, args, res, &
            xatol, xrtol, fatol, frtol, maxiter, callback)
        ! No conversion to real(dp): forward independently broadcast typed parameters.
        procedure(elementwise_root_typed_args_function) :: f !! receives typed, packed active parameters
        real(dp), intent(in) :: bracket(..) !! problem-grid inputs
        type(elementwise_typed_argument), intent(in) :: args(:) !! typed argument list
        type(find_root_result), intent(inout) :: res(..) !! grid-shaped output results
        real(dp), intent(in), optional :: xatol(..) !! broadcast absolute x tolerance
        real(dp), intent(in), optional :: xrtol(..) !! broadcast relative x tolerance
        real(dp), intent(in), optional :: fatol(..) !! broadcast absolute function tolerance
        real(dp), intent(in), optional :: frtol(..) !! broadcast relative function tolerance
        integer, intent(in), optional :: maxiter !! iteration limit
        procedure(find_root_callback), optional :: callback !! per-lane callback
        call find_root_vectorized_into(bracket=bracket, res=res, &
            xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, &
            maxiter=maxiter, callback=callback, typed_args=args, f_typed=f)
    end subroutine find_root_vectorized_args_into_typed

    subroutine bracket_root_vectorized_args_into_typed(f, xl0, args, res, &
            xr0, xmin, xmax, factor, maxiter)
        ! No conversion to real(dp): forward independently broadcast typed parameters.
        procedure(elementwise_root_typed_args_function) :: f !! receives typed, packed active parameters
        real(dp), intent(in) :: xl0(..) !! problem-grid inputs
        type(elementwise_typed_argument), intent(in) :: args(:) !! typed argument list
        type(bracket_root_result), intent(inout) :: res(..) !! grid-shaped output results
        real(dp), intent(in), optional :: xr0(..) !! optional right endpoint
        real(dp), intent(in), optional :: xmin(..) !! broadcast lower bound
        real(dp), intent(in), optional :: xmax(..) !! broadcast upper bound
        real(dp), intent(in), optional :: factor(..) !! broadcast expansion factor
        integer, intent(in), optional :: maxiter !! iteration limit
        call bracket_root_vectorized_into(xl0=xl0, res=res, &
            xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, &
            maxiter=maxiter, typed_args=args, f_typed=f)
    end subroutine bracket_root_vectorized_args_into_typed

    subroutine find_root_vectorized_indexed_args_into(f, bracket, args, res, &
            xatol, xrtol, fatol, frtol, maxiter, callback)
        ! Zero-copy typed callback: source boxes are never replicated per evaluation.
        procedure(elementwise_root_indexed_args_function) :: f !! receives source values and active-lane offsets
        real(dp), intent(in) :: bracket(..) !! initial two-point brackets
        type(elementwise_typed_argument), intent(in) :: args(:) !! mixed-type argument list
        type(find_root_result), intent(inout) :: res(..) !! grid of root estimates
        real(dp), intent(in), optional :: xatol(..) !! absolute root tolerance
        real(dp), intent(in), optional :: xrtol(..) !! relative root tolerance
        real(dp), intent(in), optional :: fatol(..) !! absolute residual tolerance
        real(dp), intent(in), optional :: frtol(..) !! relative residual tolerance
        integer, intent(in), optional :: maxiter !! iteration limit
        procedure(find_root_callback), optional :: callback !! per-lane stopping callback
        call find_root_vectorized_into(bracket=bracket, res=res, &
            xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, &
            maxiter=maxiter, callback=callback, typed_args=args, f_indexed=f)
    end subroutine find_root_vectorized_indexed_args_into

    subroutine bracket_root_vectorized_indexed_args_into(f, xl0, args, res, &
            xr0, xmin, xmax, factor, maxiter)
        procedure(elementwise_root_indexed_args_function) :: f !! receives source values and active-lane offsets
        real(dp), intent(in) :: xl0(..) !! left initial points
        type(elementwise_typed_argument), intent(in) :: args(:) !! mixed-type argument list
        type(bracket_root_result), intent(inout) :: res(..) !! grid of root brackets
        real(dp), intent(in), optional :: xr0(..) !! right initial points
        real(dp), intent(in), optional :: xmin(..) !! lower search limits
        real(dp), intent(in), optional :: xmax(..) !! upper search limits
        real(dp), intent(in), optional :: factor(..) !! expansion factors
        integer, intent(in), optional :: maxiter !! iteration limit
        call bracket_root_vectorized_into(xl0=xl0, res=res, &
            xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, &
            maxiter=maxiter, typed_args=args, f_indexed=f)
    end subroutine bracket_root_vectorized_indexed_args_into

    subroutine find_root_vectorized_into(f, bracket, res, xatol, xrtol, fatol, frtol, maxiter, &
            callback, args, f_args, args_list, typed_args, f_typed, f_indexed)
        ! Assumed-rank grid dispatcher. Result rank and shape must match input.
        ! Only this subroutine broadcasts and reshapes; the rank-one solver
        ! remains responsible for active-lane callback ordering and statuses.
        procedure(elementwise_root_vector_function), optional :: f !! packed vector objective
        real(dp), intent(in), optional :: args(..) !! leading-axis objective parameters
        type(elementwise_argument), intent(in), optional :: args_list(:) !! independently shaped args
        procedure(elementwise_root_args_function), optional :: f_args !! packed parameterized objective
        type(elementwise_typed_argument), intent(in), optional :: typed_args(:) !! typed arguments
        procedure(elementwise_root_typed_args_function), optional :: f_typed !! typed packed callback
        procedure(elementwise_root_indexed_args_function), optional :: f_indexed !! zero-copy typed callback
        real(dp), intent(in) :: bracket(..) !! grid inputs with a leading bracket axis
        type(find_root_result), intent(inout) :: res(..) !! results, same problem-grid shape
        real(dp), intent(in), optional :: xatol(..) !! broadcast tolerance argument
        real(dp), intent(in), optional :: xrtol(..) !! broadcast tolerance argument
        real(dp), intent(in), optional :: fatol(..) !! broadcast tolerance argument
        real(dp), intent(in), optional :: frtol(..) !! broadcast tolerance argument
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_root_callback), optional :: callback !! scalar per-lane stop callback
        type(find_root_result), allocatable :: flat(:)
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
            res = find_root_result()
        rank(1)
            res = find_root_result()
        rank(2)
            res = find_root_result()
        rank(3)
            res = find_root_result()
        rank(4)
            res = find_root_result()
        rank(5)
            res = find_root_result()
        rank(6)
            res = find_root_result()
        rank(7)
            res = find_root_result()
        rank(8)
            res = find_root_result()
        rank(9)
            res = find_root_result()
        rank(10)
            res = find_root_result()
        rank(11)
            res = find_root_result()
        rank(12)
            res = find_root_result()
        rank(13)
            res = find_root_result()
        rank(14)
            res = find_root_result()
        rank(15)
            res = find_root_result()
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
        if (dims(1) /= 2) return
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
        pairs = reshape(linear, [2,size(res)])
        flat = find_root_vectorized_rank1(f, pairs, &
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
    end subroutine find_root_vectorized_into

    subroutine bracket_root_vectorized_into(f, xl0, res, xr0, xmin, xmax, factor, maxiter, &
            args, f_args, args_list, typed_args, f_typed, f_indexed)
        ! Assumed-rank grid dispatcher. Result rank and shape must match input.
        ! Only this subroutine broadcasts and reshapes; the rank-one solver
        ! remains responsible for active-lane callback ordering and statuses.
        procedure(elementwise_root_vector_function), optional :: f !! packed vector objective
        real(dp), intent(in), optional :: args(..) !! leading-axis objective parameters
        type(elementwise_argument), intent(in), optional :: args_list(:) !! independently shaped args
        procedure(elementwise_root_args_function), optional :: f_args !! packed parameterized objective
        type(elementwise_typed_argument), intent(in), optional :: typed_args(:) !! typed arguments
        procedure(elementwise_root_typed_args_function), optional :: f_typed !! typed packed callback
        procedure(elementwise_root_indexed_args_function), optional :: f_indexed !! zero-copy typed callback
        real(dp), intent(in) :: xl0(..) !! grid inputs with arbitrary grid rank
        type(bracket_root_result), intent(inout) :: res(..) !! results, same problem-grid shape
        real(dp), intent(in), optional :: xr0(..) !! broadcast bracketing argument
        real(dp), intent(in), optional :: xmin(..) !! broadcast bracketing argument
        real(dp), intent(in), optional :: xmax(..) !! broadcast bracketing argument
        real(dp), intent(in), optional :: factor(..) !! broadcast bracketing argument
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_root_result), allocatable :: flat(:)
        real(dp), allocatable :: linear(:)
        real(dp), allocatable :: right(:), lo(:), hi(:), grow(:)
        integer, allocatable :: dims(:), grid(:)
        real(dp), allocatable :: packed_args(:,:)
        integer, allocatable :: typed_indices(:,:)
        logical :: valid

        ! Initialize even when shapes or broadcast parameters are invalid.
        select rank(res)
        rank(0)
            res = bracket_root_result()
        rank(1)
            res = bracket_root_result()
        rank(2)
            res = bracket_root_result()
        rank(3)
            res = bracket_root_result()
        rank(4)
            res = bracket_root_result()
        rank(5)
            res = bracket_root_result()
        rank(6)
            res = bracket_root_result()
        rank(7)
            res = bracket_root_result()
        rank(8)
            res = bracket_root_result()
        rank(9)
            res = bracket_root_result()
        rank(10)
            res = bracket_root_result()
        rank(11)
            res = bracket_root_result()
        rank(12)
            res = bracket_root_result()
        rank(13)
            res = bracket_root_result()
        rank(14)
            res = bracket_root_result()
        rank(15)
            res = bracket_root_result()
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
        dims = shape(xl0)
        grid = shape(res)
        if (size(dims) /= size(grid)) return
        if (any(dims /= grid)) return
        if (size(res) == 0) return
        valid = .true.
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
        call elementwise_flatten(xl0, linear)
        flat = bracket_root_vectorized_rank1(f, linear, &
            xr0=right, xmin=lo, xmax=hi, factor=grow, maxiter=maxiter, args_data=packed_args, f_args=f_args, &
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
    end subroutine bracket_root_vectorized_into


    function find_root_vectorized_rank3(f, bracket, xatol, xrtol, fatol, frtol, maxiter, &
            callback) result(res)
        ! Compatibility function; numerical work lives in the assumed-rank dispatcher.
        procedure(elementwise_root_vector_function) :: f !! packed objective callback
        real(dp), intent(in) :: bracket(:,:,:,:) !! input grid
        real(dp), intent(in), optional :: xatol(..) !! broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_root_callback), optional :: callback !! scalar per-lane callback
        type(find_root_result) :: res(size(bracket,2),size(bracket,3),size(bracket,4))

        call find_root_vectorized_into(f, bracket, res, &
            xatol, xrtol, fatol, frtol, maxiter, callback)
    end function find_root_vectorized_rank3


    function find_root_vectorized_rank4(f, bracket, xatol, xrtol, fatol, frtol, maxiter, &
            callback) result(res)
        ! Compatibility function; numerical work lives in the assumed-rank dispatcher.
        procedure(elementwise_root_vector_function) :: f !! packed objective callback
        real(dp), intent(in) :: bracket(:,:,:,:,:) !! input grid
        real(dp), intent(in), optional :: xatol(..) !! broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_root_callback), optional :: callback !! scalar per-lane callback
        type(find_root_result) :: res(size(bracket,2),size(bracket,3),size(bracket,4),size(bracket,5))

        call find_root_vectorized_into(f, bracket, res, &
            xatol, xrtol, fatol, frtol, maxiter, callback)
    end function find_root_vectorized_rank4


    function bracket_root_vectorized_rank3(f, xl0, xr0, xmin, xmax, factor, &
            maxiter) result(res)
        ! Compatibility function; numerical work lives in the assumed-rank dispatcher.
        procedure(elementwise_root_vector_function) :: f !! packed objective callback
        real(dp), intent(in) :: xl0(:,:,:) !! input grid
        real(dp), intent(in), optional :: xr0(..) !! broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! broadcast parameter
        real(dp), intent(in), optional :: factor(..) !! broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_root_result) :: res(size(xl0,1),size(xl0,2),size(xl0,3))

        call bracket_root_vectorized_into(f, xl0, res, &
            xr0, xmin, xmax, factor, maxiter)
    end function bracket_root_vectorized_rank3


    function bracket_root_vectorized_rank4(f, xl0, xr0, xmin, xmax, factor, &
            maxiter) result(res)
        ! Compatibility function; numerical work lives in the assumed-rank dispatcher.
        procedure(elementwise_root_vector_function) :: f !! packed objective callback
        real(dp), intent(in) :: xl0(:,:,:,:) !! input grid
        real(dp), intent(in), optional :: xr0(..) !! broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! broadcast parameter
        real(dp), intent(in), optional :: factor(..) !! broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_root_result) :: res(size(xl0,1),size(xl0,2),size(xl0,3),size(xl0,4))

        call bracket_root_vectorized_into(f, xl0, res, &
            xr0, xmin, xmax, factor, maxiter)
    end function bracket_root_vectorized_rank4

    function find_root_vectorized_rank2(f, bracket, xatol, xrtol, fatol, frtol, maxiter, &
            callback) result(res)
        ! Compatibility function; numerical work lives in the assumed-rank dispatcher.
        procedure(elementwise_root_vector_function) :: f !! packed objective callback
        real(dp), intent(in) :: bracket(:,:,:) !! input grid
        real(dp), intent(in), optional :: xatol(..) !! broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_root_callback), optional :: callback !! scalar per-lane callback
        type(find_root_result) :: res(size(bracket,2),size(bracket,3))

        call find_root_vectorized_into(f, bracket, res, &
            xatol, xrtol, fatol, frtol, maxiter, callback)
    end function find_root_vectorized_rank2


    function bracket_root_vectorized_rank2(f, xl0, xr0, xmin, xmax, factor, &
            maxiter) result(res)
        ! Compatibility function; numerical work lives in the assumed-rank dispatcher.
        procedure(elementwise_root_vector_function) :: f !! packed objective callback
        real(dp), intent(in) :: xl0(:,:) !! input grid
        real(dp), intent(in), optional :: xr0(..) !! broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! broadcast parameter
        real(dp), intent(in), optional :: factor(..) !! broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_root_result) :: res(size(xl0,1),size(xl0,2))

        call bracket_root_vectorized_into(f, xl0, res, &
            xr0, xmin, xmax, factor, maxiter)
    end function bracket_root_vectorized_rank2


    subroutine find_root_batch(f, bracket, batch_shape, xatol, xrtol, fatol, frtol, maxiter, callback, &
            res)
        implicit none
        ! Subroutine dispatch avoids temporary descriptors for absent assumed-rank
        ! arguments in array-valued function calls on affected compilers.
        ! One independent scalar solve per element; reshape by thin rank wrappers.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:, :) !! flattened problem inputs
        integer, intent(in) :: batch_shape(:) !! desired output dimensions
        real(dp), intent(in), optional :: xatol(..) !! broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! broadcast parameter
        integer, intent(in), optional :: maxiter !! shared iteration limit
        procedure(find_root_callback), optional :: callback !! scalar callback
        type(find_root_result), allocatable, intent(out) :: res(:) !! flattened independent solver results
        logical :: valid
        integer :: i, n
        real(dp), allocatable :: ax(:), rx(:), af(:), rf(:)

        n = size(bracket, 2)
        allocate(res(n), ax(n), rx(n), af(n), rf(n))
        if (size(bracket, 1) /= 2) return
        if (n /= product(batch_shape)) return
        ax = 4.0_dp * tiny(1.0_dp)
        rx = 4.0_dp * epsilon(1.0_dp)
        af = tiny(1.0_dp)
        rf = 0.0_dp
        valid = .true.
        call elementwise_broadcast(xatol, batch_shape, ax, valid)
        call elementwise_broadcast(xrtol, batch_shape, rx, valid)
        call elementwise_broadcast(fatol, batch_shape, af, valid)
        call elementwise_broadcast(frtol, batch_shape, rf, valid)
        if (.not. valid) return
        do i = 1, n
            res(i) = find_root_scalar(f, bracket(:, i), xatol=ax(i), xrtol=rx(i), &
                fatol=af(i), frtol=rf(i), maxiter=maxiter, callback=callback)
        end do
    end subroutine find_root_batch

    function find_root_vectorized_rank1(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback, args_data, f_args, &
            typed_args, typed_indices, f_typed, f_indexed) result(res)
        ! Solve independent roots with one vector call per active iteration.
        ! Unlike find_root(f, bracket), f here takes arrays of active points.
        ! The scalar Chandrupatla controller is retained for each element.
        procedure(elementwise_root_vector_function), optional :: f !! active-point vector callback
        real(dp), intent(in), optional :: args_data(:,:) !! prebroadcast (nargs,nproblems)
        procedure(elementwise_root_args_function), optional :: f_args !! packed parameterized objective
        type(elementwise_typed_argument), intent(in), optional :: typed_args(:) !! typed input snapshots
        integer, intent(in), optional :: typed_indices(:,:) !! source offsets for every lane
        procedure(elementwise_root_typed_args_function), optional :: f_typed !! typed packed callback
        procedure(elementwise_root_indexed_args_function), optional :: f_indexed !! zero-copy typed callback
        real(dp), intent(in) :: bracket(:,:) !! endpoints (2,n)
        real(dp), intent(in), optional :: xatol(..) !! scalar or length-n tolerance
        real(dp), intent(in), optional :: xrtol(..) !! scalar or length-n tolerance
        real(dp), intent(in), optional :: fatol(..) !! scalar or length-n tolerance
        real(dp), intent(in), optional :: frtol(..) !! scalar or length-n tolerance
        integer, intent(in), optional :: maxiter !! common maximum iteration count
        procedure(find_root_callback), optional :: callback !! scalar state callback
        type(find_root_result) :: res(size(bracket, 2))

        integer :: n, niter, i, j, m, it
        integer, allocatable :: ids(:)
        integer, allocatable :: packed_ids(:)
        real(dp), allocatable :: ax(:), rx(:), af(:), rf(:), ftarget(:)
        real(dp), allocatable :: x1(:), x2(:), f1(:), f2(:), t(:), tol(:), dx(:)
        real(dp), allocatable :: points(:), values(:)
        type(elementwise_typed_argument), allocatable :: active_typed(:)
        integer, allocatable :: packed_typed_indices(:,:)
        real(dp) :: x3, f3, fnew, xnew, fscale, next_t, xi1, phi1
        real(dp) :: alpha, denom, p1, p2, tl
        logical, allocatable :: active(:)
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
        n = size(bracket, 2)
        if (size(bracket, 1) /= 2) return
        if (n == 0) return
        niter = 2046
        if (present(maxiter)) niter = maxiter
        if (niter < 0) return

        allocate(ax(n), rx(n), af(n), rf(n), ftarget(n))
        allocate(x1(n), x2(n), f1(n), f2(n), t(n), tol(n), dx(n))
        allocate(ids(n), points(2*n), values(2*n), active(n))
        if (present(f_indexed)) allocate(packed_typed_indices(size(typed_indices,1), 2*n))
        if (present(f_args) .or. present(f_typed) .or. present(f_indexed)) allocate(packed_ids(2*n))
        ax = 4.0_dp * tiny(1.0_dp)
        rx = 4.0_dp * epsilon(1.0_dp)
        af = tiny(1.0_dp)
        rf = 0.0_dp
        valid = .true.
        call elementwise_broadcast(xatol, [n], ax, valid)
        call elementwise_broadcast(xrtol, [n], rx, valid)
        call elementwise_broadcast(fatol, [n], af, valid)
        call elementwise_broadcast(frtol, [n], rf, valid)
        if (.not. valid) return

        ! Validate each lane before forming a function-evaluation batch.
        ! In particular, check finiteness before numerical comparisons so
        ! that a NaN input is safe with IEEE invalid-operation traps enabled.
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
            points(2*m-1:2*m) = bracket(:,i)
            if (present(f_args) .or. present(f_typed) .or. present(f_indexed)) packed_ids(2*m-1:2*m) = i
        end do
        if (m == 0) return

        ! A single callback evaluates both endpoints of every valid lane.
        if (present(f_indexed)) then
            call elementwise_gather_typed_indices(typed_indices, packed_ids(:2*m), &
                packed_typed_indices(:, :2*m))
            call f_indexed(points(:2*m), typed_args, &
                packed_typed_indices(:, :2*m), values(:2*m))
        else if (present(f_typed)) then
            call elementwise_pack_typed_arguments(typed_args, &
                typed_indices(:,packed_ids(:2*m)), active_typed)
            call f_typed(points(:2*m), active_typed, values(:2*m))
        else if (present(f_args)) then
            call f_args(points(:2*m), args_data(:,packed_ids(:2*m)), values(:2*m))
        else
            call f(points(:2*m), values(:2*m))
        end if
        do j = 1, m
            i = ids(j)
            x1(i) = bracket(1,i)
            x2(i) = bracket(2,i)
            f1(i) = values(2*j-1)
            f2(i) = values(2*j)
            res(i)%nfev = 2
            call store_find_root_state(res(i), x1(i), x2(i), f1(i), f2(i))
            if (.not. ieee_is_finite(f1(i)) .or. .not. ieee_is_finite(f2(i))) then
                res(i)%status = root_status_value_error
                call invalidate_find_root(res(i))
                active(i) = .false.
                cycle
            end if
            ftarget(i) = 0.0_dp
            if (rf(i) > 0.0_dp) then
                fscale = min(abs(f1(i)), abs(f2(i)))
                if (rf(i) > 1.0_dp) then
                    if (fscale > huge(1.0_dp) / rf(i)) then
                        ftarget(i) = huge(1.0_dp)
                    else
                        ftarget(i) = rf(i) * fscale
                    end if
                else
                    ftarget(i) = rf(i) * fscale
                end if
            end if
            t(i) = 0.5_dp
            call check_find_root(res(i), x1(i), x2(i), f1(i), f2(i), &
                ax(i), rx(i), af(i), ftarget(i), tol(i), dx(i))
            if (present(callback)) then
                stopped = .false.
                call callback(res(i), stopped)
                if (stopped .and. res(i)%status == 1) res(i)%status = find_root_status_callback
            end if
            active(i) = (res(i)%status == 1)
        end do

        do it = 1, niter
            m = 0
            do i = 1, n
                if (.not. active(i)) cycle
                xnew = x1(i) + t(i) * (x2(i) - x1(i))
                if (.not. ieee_is_finite(xnew)) then
                    res(i)%status = root_status_value_error
                    call invalidate_find_root(res(i))
                    active(i) = .false.
                    cycle
                end if
                m = m + 1
                ids(m) = i
                points(m) = xnew
            end do
            if (m == 0) exit
            ! Shrink the callback input as individual problems converge.
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
                xnew = points(j)
                fnew = values(j)
                res(i)%nfev = res(i)%nfev + 1
                res(i)%nit = it
                x3 = x2(i)
                f3 = f2(i)
                if (ieee_is_finite(fnew)) then
                    if (same_root_sign(fnew, f1(i))) then
                        x3 = x1(i)
                        f3 = f1(i)
                    else
                        x2(i) = x1(i)
                        f2(i) = f1(i)
                    end if
                end if
                x1(i) = xnew
                f1(i) = fnew
                if (.not. ieee_is_finite(fnew)) then
                    res(i)%status = root_status_value_error
                    call store_find_root_state(res(i), x1(i), x2(i), f1(i), f2(i))
                    call invalidate_find_root(res(i))
                    active(i) = .false.
                    cycle
                end if
                call check_find_root(res(i), x1(i), x2(i), f1(i), f2(i), &
                    ax(i), rx(i), af(i), ftarget(i), tol(i), dx(i))
                if (present(callback)) then
                    stopped = .false.
                    call callback(res(i), stopped)
                    if (stopped .and. res(i)%status == 1) res(i)%status = find_root_status_callback
                end if
                active(i) = (res(i)%status == 1)
                if (.not. active(i)) cycle

                ! Same safeguarded inverse-quadratic step as the scalar kernel.
                next_t = 0.5_dp
                if (x3 /= x2(i) .and. f3 /= f2(i) .and. x2(i) /= x1(i)) then
                    xi1 = (x1(i) - x2(i)) / (x3 - x2(i))
                    if (xi1 > 0.0_dp .and. xi1 < 1.0_dp) then
                        phi1 = (f1(i) - f2(i)) / (f3 - f2(i))
                        if (phi1 > 1.0_dp - sqrt(1.0_dp - xi1) .and. phi1 < sqrt(xi1)) then
                            alpha = (x3 - x1(i)) / (x2(i) - x1(i))
                            denom = f3 - f1(i)
                            if (f1(i) /= f2(i) .and. f3 /= f2(i) .and. &
                                denom /= 0.0_dp .and. f2(i) /= f3) then
                                p1 = (f1(i) / (f1(i)-f2(i))) * (f3 / (f3-f2(i)))
                                p2 = alpha * (f1(i) / denom) * (f2(i) / (f2(i)-f3))
                                if (ieee_is_finite(p1) .and. ieee_is_finite(p2)) next_t = p1-p2
                            end if
                        end if
                    end if
                end if
                if (dx(i) > 0.0_dp) then
                    tl = 0.5_dp * tol(i) / dx(i)
                    t(i) = max(tl, min(1.0_dp-tl, next_t))
                else
                    t(i) = 0.5_dp
                end if
            end do
        end do
        do i = 1, n
            if (res(i)%status == 1) res(i)%status = root_status_convergence_error
            res(i)%success = (res(i)%status == root_status_success)
        end do
    end function find_root_vectorized_rank1

    function find_root_vector(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 1-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_root_callback), optional :: callback !! per-element callback
        type(find_root_result), allocatable :: flat_result(:)
        type(find_root_result) :: res(size(bracket, 2))

        call find_root_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_root_vector

    function find_root_matrix(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 2-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_root_callback), optional :: callback !! per-element callback
        type(find_root_result), allocatable :: flat_result(:)
        type(find_root_result) :: res(size(bracket, 2), size(bracket, 3))

        call find_root_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_root_matrix

    function find_root_tensor(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 3-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_root_callback), optional :: callback !! per-element callback
        type(find_root_result), allocatable :: flat_result(:)
        type(find_root_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4))

        call find_root_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_root_tensor

    function find_root_rank4(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 4-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_root_callback), optional :: callback !! per-element callback
        type(find_root_result), allocatable :: flat_result(:)
        type(find_root_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), size(bracket, 5))

        call find_root_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_root_rank4

    function find_root_rank5(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 5-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_root_callback), optional :: callback !! per-element callback
        type(find_root_result), allocatable :: flat_result(:)
        type(find_root_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), size(bracket, 5), size(bracket, 6))

        call find_root_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_root_rank5

    function find_root_rank6(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 6-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_root_callback), optional :: callback !! per-element callback
        type(find_root_result), allocatable :: flat_result(:)
        type(find_root_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), &
            size(bracket, 5), size(bracket, 6), size(bracket, 7))

        call find_root_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_root_rank6

    function find_root_rank7(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 7-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_root_callback), optional :: callback !! per-element callback
        type(find_root_result), allocatable :: flat_result(:)
        type(find_root_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), &
            size(bracket, 5), size(bracket, 6), size(bracket, 7), size(bracket, 8))

        call find_root_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_root_rank7

    function find_root_rank8(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 8-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_root_callback), optional :: callback !! per-element callback
        type(find_root_result), allocatable :: flat_result(:)
        type(find_root_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), &
            size(bracket, 5), size(bracket, 6), size(bracket, 7), size(bracket, 8), size(bracket, 9))

        call find_root_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_root_rank8

    function find_root_rank9(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 9-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:,:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_root_callback), optional :: callback !! per-element callback
        type(find_root_result), allocatable :: flat_result(:)
        type(find_root_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), &
            size(bracket, 5), size(bracket, 6), size(bracket, 7), size(bracket, 8), size(bracket, 9), &
            size(bracket, 10))

        call find_root_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_root_rank9

    function find_root_rank10(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 10-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:,:,:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_root_callback), optional :: callback !! per-element callback
        type(find_root_result), allocatable :: flat_result(:)
        type(find_root_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), &
            size(bracket, 5), size(bracket, 6), size(bracket, 7), size(bracket, 8), size(bracket, 9), &
            size(bracket, 10), size(bracket, 11))

        call find_root_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_root_rank10

    function find_root_rank11(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 11-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:,:,:,:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_root_callback), optional :: callback !! per-element callback
        type(find_root_result), allocatable :: flat_result(:)
        type(find_root_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), &
            size(bracket, 5), size(bracket, 6), size(bracket, 7), size(bracket, 8), size(bracket, 9), &
            size(bracket, 10), size(bracket, 11), size(bracket, 12))

        call find_root_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_root_rank11

    function find_root_rank12(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 12-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:,:,:,:,:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_root_callback), optional :: callback !! per-element callback
        type(find_root_result), allocatable :: flat_result(:)
        type(find_root_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), &
            size(bracket, 5), size(bracket, 6), size(bracket, 7), size(bracket, 8), size(bracket, 9), &
            size(bracket, 10), size(bracket, 11), size(bracket, 12), size(bracket, 13))

        call find_root_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_root_rank12

    function find_root_rank13(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 13-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_root_callback), optional :: callback !! per-element callback
        type(find_root_result), allocatable :: flat_result(:)
        type(find_root_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), &
            size(bracket, 5), size(bracket, 6), size(bracket, 7), size(bracket, 8), size(bracket, 9), &
            size(bracket, 10), size(bracket, 11), size(bracket, 12), size(bracket, 13), size(bracket, 14))

        call find_root_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_root_rank13

    function find_root_rank14(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Shape-preserving independent 14-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: bracket(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! problem brackets
        real(dp), intent(in), optional :: xatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xrtol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: fatol(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: frtol(..) !! trailing-axis broadcast parameter
        integer, intent(in), optional :: maxiter !! common iteration limit
        procedure(find_root_callback), optional :: callback !! per-element callback
        type(find_root_result), allocatable :: flat_result(:)
        type(find_root_result) :: res(size(bracket, 2), size(bracket, 3), size(bracket, 4), &
            size(bracket, 5), size(bracket, 6), size(bracket, 7), size(bracket, 8), size(bracket, 9), &
            size(bracket, 10), size(bracket, 11), size(bracket, 12), size(bracket, 13), size(bracket, 14), &
            size(bracket, 15))

        call find_root_batch(f, reshape(bracket, [size(bracket, 1), size(res)]), &
            shape(res), xatol=xatol, xrtol=xrtol, fatol=fatol, frtol=frtol, maxiter=maxiter, callback=callback, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function find_root_rank14

    subroutine bracket_root_batch(f, xl0, batch_shape, xr0, xmin, xmax, factor, maxiter, &
            res)
        implicit none
        ! Subroutine dispatch avoids temporary descriptors for absent assumed-rank
        ! arguments in array-valued function calls on affected compilers.
        ! One independent scalar solve per element; reshape by thin rank wrappers.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: xl0(:) !! flattened problem inputs
        integer, intent(in) :: batch_shape(:) !! desired output dimensions
        real(dp), intent(in), optional :: xr0(..) !! broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! broadcast parameter
        real(dp), intent(in), optional :: factor !! shared growth factor
        integer, intent(in), optional :: maxiter !! shared iteration limit
        type(bracket_root_result), allocatable, intent(out) :: res(:) !! flattened independent solver results
        logical :: valid
        integer :: i, n
        real(dp), allocatable :: right(:), lo(:), hi(:)

        n = size(xl0)
        allocate(res(n), right(n), lo(n), hi(n))
        if (n /= product(batch_shape)) return
        right = 0.0_dp
        lo = ieee_value(0.0_dp, ieee_negative_inf)
        hi = ieee_value(0.0_dp, ieee_positive_inf)
        valid = .true.
        call elementwise_broadcast(xr0, batch_shape, right, valid)
        call elementwise_broadcast(xmin, batch_shape, lo, valid)
        call elementwise_broadcast(xmax, batch_shape, hi, valid)
        if (.not. valid) return
        do i = 1, n
            if (present(xr0)) then
                res(i) = bracket_root_scalar(f, xl0(i), xr0=right(i), &
                    xmin=lo(i), xmax=hi(i), factor=factor, maxiter=maxiter)
            else
                res(i) = bracket_root_scalar(f, xl0(i), &
                    xmin=lo(i), xmax=hi(i), factor=factor, maxiter=maxiter)
            end if
        end do
    end subroutine bracket_root_batch

    function bracket_root_vectorized_rank1(f, xl0, xr0, xmin, xmax, factor, maxiter, args_data, f_args, &
            typed_args, typed_indices, f_typed, f_indexed) result(res)
        ! Pack initial endpoint evaluations and subsequent active expansion
        ! points into vector calls. Each lane keeps SciPy's scalar rules.
        procedure(elementwise_root_vector_function), optional :: f !! packed vector callback
        real(dp), intent(in), optional :: args_data(:,:) !! prebroadcast (nargs,nproblems)
        procedure(elementwise_root_args_function), optional :: f_args !! packed parameterized objective
        type(elementwise_typed_argument), intent(in), optional :: typed_args(:) !! typed input snapshots
        integer, intent(in), optional :: typed_indices(:,:) !! source offsets for every lane
        procedure(elementwise_root_typed_args_function), optional :: f_typed !! typed packed callback
        procedure(elementwise_root_indexed_args_function), optional :: f_indexed !! zero-copy typed callback
        real(dp), intent(in) :: xl0(:) !! initial left endpoint for each lane
        real(dp), intent(in), optional :: xr0(..) !! broadcast right endpoint
        real(dp), intent(in), optional :: xmin(..) !! broadcast lower search limit
        real(dp), intent(in), optional :: xmax(..) !! broadcast upper search limit
        real(dp), intent(in), optional :: factor(..) !! broadcast growth factors
        integer, intent(in), optional :: maxiter !! common expansion limit
        type(bracket_root_result) :: res(size(xl0))

        type(vector_bracket_root_lane), allocatable :: lane(:)
        real(dp), allocatable :: right(:), lo(:), hi(:), grow(:)
        real(dp), allocatable :: points(:), values(:)
        type(elementwise_typed_argument), allocatable :: active_typed(:)
        integer, allocatable :: packed_typed_indices(:,:)
        integer, allocatable :: ids(:), sides(:)
        logical, allocatable :: work(:), selected(:)
        real(dp) :: xr, space, next_x
        integer :: n, i, j, k, m, it, niter
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
            if (size(args_data,2) /= size(xl0)) return
        end if
        if (present(typed_args)) then
            if (size(typed_args) == 0) return
            if (size(typed_indices,1) /= size(typed_args)) return
            if (size(typed_indices,2) /= size(xl0)) return
        end if
        n = size(xl0)
        if (n == 0) return
        niter = 1000
        if (present(maxiter)) niter = maxiter
        if (niter < 0) return
        allocate(lane(n), work(n), selected(n), right(n), lo(n), hi(n), grow(n))
        allocate(points(2*n), values(2*n), ids(2*n), sides(2*n))
        if (present(f_indexed)) allocate(packed_typed_indices(size(typed_indices,1), 2*n))
        right = 0.0_dp
        lo = ieee_value(0.0_dp, ieee_negative_inf)
        hi = ieee_value(0.0_dp, ieee_positive_inf)
        grow = 2.0_dp
        valid = .true.
        call elementwise_broadcast(xr0, [n], right, valid)
        call elementwise_broadcast(xmin, [n], lo, valid)
        call elementwise_broadcast(xmax, [n], hi, valid)
        call elementwise_broadcast(factor, [n], grow, valid)
        if (.not. valid) return

        selected = .false.
        do i = 1, n
            res(i)%xl = xl0(i)
            res(i)%xr = xl0(i)
            if (.not. ieee_is_finite(xl0(i))) cycle
            if (ieee_is_nan(lo(i)) .or. ieee_is_nan(hi(i))) cycle
            if (.not. ieee_is_finite(grow(i))) cycle
            if (grow(i) <= 1.0_dp .or. lo(i) > xl0(i)) cycle
            xr = right(i)
            if (.not. present(xr0)) then
                if (hi(i) <= xl0(i)) cycle
                if (ieee_is_finite(hi(i))) then
                    if (xl0(i) < 0.0_dp .and. hi(i) > 0.0_dp) then
                        space = hi(i) / 8.0_dp - xl0(i) / 8.0_dp
                    else
                        space = (hi(i) - xl0(i)) / 8.0_dp
                    end if
                    xr = xl0(i) + min(space, 1.0_dp)
                else
                    if (xl0(i) > huge(1.0_dp) - 1.0_dp) cycle
                    xr = xl0(i) + 1.0_dp
                end if
            end if
            res(i)%xr = xr
            if (.not. ieee_is_finite(xr)) cycle
            if (xl0(i) >= xr .or. hi(i) < xr) cycle
            lane(i)%x = [xl0(i), xr]
            lane(i)%bounded = [ieee_is_finite(lo(i)), ieee_is_finite(hi(i))]
            lane(i)%bound = [lo(i), hi(i)]
            selected(i) = .true.
        end do

        ! One packed call for both endpoints of all valid lanes.
        m = 0
        do i = 1, n
            if (.not. selected(i)) cycle
            do j = 1, 2
                m = m + 1
                points(m) = lane(i)%x(j)
                ids(m) = i
                sides(m) = j
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
                j = sides(k)
                lane(i)%fx(j) = values(k)
                res(i)%nfev = res(i)%nfev + 1
            end do
        end if
        work = .false.
        do i = 1, n
            if (.not. selected(i)) cycle
            res(i)%fl = lane(i)%fx(1)
            res(i)%fr = lane(i)%fx(2)
            if (.not. all(ieee_is_finite(lane(i)%fx))) then
                res(i)%status = root_status_value_error
                cycle
            end if
            if (bracket_has_root(lane(i)%fx(1), lane(i)%fx(2))) then
                res(i)%status = root_status_success
                res(i)%success = .true.
                cycle
            end if
            if (xl0(i) < 0.0_dp .and. lane(i)%x(2) > 0.0_dp) then
                if (lane(i)%x(2) > huge(1.0_dp) + xl0(i)) then
                    res(i)%status = root_status_value_error
                    cycle
                end if
            end if
            xr = lane(i)%x(2)
            lane(i)%anchor = [xr, xl0(i)]
            lane(i)%xlast = [xr, xl0(i)]
            lane(i)%flast = [lane(i)%fx(2), lane(i)%fx(1)]
            lane(i)%distance = [xl0(i) - xr, xr - xl0(i)]
            if (lane(i)%bounded(1)) lane(i)%distance(1) = lo(i) - xl0(i)
            if (lane(i)%bounded(2)) lane(i)%distance(2) = hi(i) - xr
            lane(i)%active = .true.
            if (lane(i)%bounded(1)) lane(i)%active(1) = xl0(i) > lo(i)
            if (lane(i)%bounded(2)) lane(i)%active(2) = xr < hi(i)
            if (.not. any(lane(i)%active)) then
                res(i)%status = bracket_root_status_limit
            else
                work(i) = .true.
            end if
        end do

        do it = 1, niter
            if (.not. any(work)) exit
            m = 0
            do i = 1, n
                if (.not. work(i)) cycle
                res(i)%nit = it
                do j = 1, 2
                    if (.not. lane(i)%active(j)) cycle
                    if (lane(i)%bounded(j)) then
                        lane(i)%distance(j) = lane(i)%distance(j) / grow(i)
                        next_x = lane(i)%bound(j) - lane(i)%distance(j)
                    else
                        if (abs(lane(i)%distance(j)) > huge(1.0_dp) / grow(i)) then
                            lane(i)%nonfinite = .true.
                            lane(i)%active(j) = .false.
                            cycle
                        end if
                        lane(i)%distance(j) = lane(i)%distance(j) * grow(i)
                        if (abs(lane(i)%anchor(j)) > huge(1.0_dp) - abs(lane(i)%distance(j))) then
                            lane(i)%nonfinite = .true.
                            lane(i)%active(j) = .false.
                            cycle
                        end if
                        next_x = lane(i)%anchor(j) + lane(i)%distance(j)
                    end if
                    if (.not. ieee_is_finite(next_x)) then
                        lane(i)%nonfinite = .true.
                        lane(i)%active(j) = .false.
                        cycle
                    end if
                    if (next_x == lane(i)%x(j)) then
                        lane(i)%active(j) = .false.
                        cycle
                    end if
                    lane(i)%xlast(j) = lane(i)%x(j)
                    lane(i)%flast(j) = lane(i)%fx(j)
                    lane(i)%x(j) = next_x
                    m = m + 1
                    ids(m) = i
                    sides(m) = j
                    points(m) = next_x
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
                    j = sides(k)
                    lane(i)%fx(j) = values(k)
                    res(i)%nfev = res(i)%nfev + 1
                    if (.not. ieee_is_finite(values(k))) then
                        lane(i)%nonfinite = .true.
                        lane(i)%active(j) = .false.
                        cycle
                    end if
                    if (bracket_has_root(lane(i)%fx(j), lane(i)%flast(j))) then
                        lane(i)%found(j) = .true.
                        lane(i)%active(j) = .false.
                    else if (lane(i)%bounded(j)) then
                        if (lane(i)%x(j) == lane(i)%bound(j)) lane(i)%active(j) = .false.
                    end if
                end do
            end if
            do i = 1, n
                if (.not. work(i)) cycle
                if (any(lane(i)%found) .or. .not. any(lane(i)%active)) work(i) = .false.
            end do
        end do

        do i = 1, n
            if (.not. selected(i)) cycle
            if (res(i)%status /= root_status_invalid_input) cycle
            ! An initially valid but subsequently invalid objective is already
            ! marked above. Remaining lanes have visited the expansion loop.
            if (res(i)%nfev == 0) cycle
            res(i)%xl = lane(i)%x(1)
            res(i)%xr = lane(i)%x(2)
            res(i)%fl = lane(i)%fx(1)
            res(i)%fr = lane(i)%fx(2)
            if (any(lane(i)%found)) then
                if (all(lane(i)%found)) then
                    if (lane(i)%xlast(1) - lane(i)%x(1) <= lane(i)%x(2) - lane(i)%xlast(2)) then
                        res(i)%xr = lane(i)%xlast(1)
                        res(i)%fr = lane(i)%flast(1)
                    else
                        res(i)%xl = lane(i)%xlast(2)
                        res(i)%fl = lane(i)%flast(2)
                    end if
                else if (lane(i)%found(1)) then
                    res(i)%xr = lane(i)%xlast(1)
                    res(i)%fr = lane(i)%flast(1)
                else
                    res(i)%xl = lane(i)%xlast(2)
                    res(i)%fl = lane(i)%flast(2)
                end if
                res(i)%status = root_status_success
                res(i)%success = .true.
            else if (lane(i)%nonfinite) then
                res(i)%status = root_status_value_error
            else if (.not. any(lane(i)%active)) then
                res(i)%status = bracket_root_status_limit
            else
                res(i)%status = root_status_convergence_error
            end if
        end do
    end function bracket_root_vectorized_rank1

    function bracket_root_vector(f, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 1-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: xl0(:) !! problem starting points
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_root_result), allocatable :: flat_result(:)
        type(bracket_root_result) :: res(size(xl0, 1))

        call bracket_root_batch(f, reshape(xl0, [size(xl0)]), &
            shape(res), xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_root_vector

    function bracket_root_matrix(f, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 2-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: xl0(:,:) !! problem starting points
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_root_result), allocatable :: flat_result(:)
        type(bracket_root_result) :: res(size(xl0, 1), size(xl0, 2))

        call bracket_root_batch(f, reshape(xl0, [size(xl0)]), &
            shape(res), xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_root_matrix

    function bracket_root_tensor(f, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 3-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: xl0(:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_root_result), allocatable :: flat_result(:)
        type(bracket_root_result) :: res(size(xl0, 1), size(xl0, 2), size(xl0, 3))

        call bracket_root_batch(f, reshape(xl0, [size(xl0)]), &
            shape(res), xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_root_tensor

    function bracket_root_rank4(f, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 4-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: xl0(:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_root_result), allocatable :: flat_result(:)
        type(bracket_root_result) :: res(size(xl0, 1), size(xl0, 2), size(xl0, 3), size(xl0, 4))

        call bracket_root_batch(f, reshape(xl0, [size(xl0)]), &
            shape(res), xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_root_rank4

    function bracket_root_rank5(f, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 5-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: xl0(:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_root_result), allocatable :: flat_result(:)
        type(bracket_root_result) :: res(size(xl0, 1), size(xl0, 2), size(xl0, 3), size(xl0, 4), size(xl0, 5))

        call bracket_root_batch(f, reshape(xl0, [size(xl0)]), &
            shape(res), xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_root_rank5

    function bracket_root_rank6(f, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 6-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: xl0(:,:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_root_result), allocatable :: flat_result(:)
        type(bracket_root_result) :: res(size(xl0, 1), size(xl0, 2), size(xl0, 3), size(xl0, 4), size(xl0, 5), size(xl0, 6))

        call bracket_root_batch(f, reshape(xl0, [size(xl0)]), &
            shape(res), xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_root_rank6

    function bracket_root_rank7(f, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 7-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: xl0(:,:,:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_root_result), allocatable :: flat_result(:)
        type(bracket_root_result) :: res(size(xl0, 1), size(xl0, 2), size(xl0, 3), size(xl0, 4), &
            size(xl0, 5), size(xl0, 6), size(xl0, 7))

        call bracket_root_batch(f, reshape(xl0, [size(xl0)]), &
            shape(res), xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_root_rank7

    function bracket_root_rank8(f, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 8-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: xl0(:,:,:,:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_root_result), allocatable :: flat_result(:)
        type(bracket_root_result) :: res(size(xl0, 1), size(xl0, 2), size(xl0, 3), size(xl0, 4), &
            size(xl0, 5), size(xl0, 6), size(xl0, 7), size(xl0, 8))

        call bracket_root_batch(f, reshape(xl0, [size(xl0)]), &
            shape(res), xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_root_rank8

    function bracket_root_rank9(f, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 9-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: xl0(:,:,:,:,:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_root_result), allocatable :: flat_result(:)
        type(bracket_root_result) :: res(size(xl0, 1), size(xl0, 2), size(xl0, 3), size(xl0, 4), &
            size(xl0, 5), size(xl0, 6), size(xl0, 7), size(xl0, 8), size(xl0, 9))

        call bracket_root_batch(f, reshape(xl0, [size(xl0)]), &
            shape(res), xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_root_rank9

    function bracket_root_rank10(f, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 10-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: xl0(:,:,:,:,:,:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_root_result), allocatable :: flat_result(:)
        type(bracket_root_result) :: res(size(xl0, 1), size(xl0, 2), size(xl0, 3), size(xl0, 4), &
            size(xl0, 5), size(xl0, 6), size(xl0, 7), size(xl0, 8), size(xl0, 9), size(xl0, 10))

        call bracket_root_batch(f, reshape(xl0, [size(xl0)]), &
            shape(res), xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_root_rank10

    function bracket_root_rank11(f, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 11-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: xl0(:,:,:,:,:,:,:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_root_result), allocatable :: flat_result(:)
        type(bracket_root_result) :: res(size(xl0, 1), size(xl0, 2), size(xl0, 3), size(xl0, 4), &
            size(xl0, 5), size(xl0, 6), size(xl0, 7), size(xl0, 8), size(xl0, 9), size(xl0, 10), &
            size(xl0, 11))

        call bracket_root_batch(f, reshape(xl0, [size(xl0)]), &
            shape(res), xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_root_rank11

    function bracket_root_rank12(f, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 12-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: xl0(:,:,:,:,:,:,:,:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_root_result), allocatable :: flat_result(:)
        type(bracket_root_result) :: res(size(xl0, 1), size(xl0, 2), size(xl0, 3), size(xl0, 4), &
            size(xl0, 5), size(xl0, 6), size(xl0, 7), size(xl0, 8), size(xl0, 9), size(xl0, 10), &
            size(xl0, 11), size(xl0, 12))

        call bracket_root_batch(f, reshape(xl0, [size(xl0)]), &
            shape(res), xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_root_rank12

    function bracket_root_rank13(f, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 13-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: xl0(:,:,:,:,:,:,:,:,:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_root_result), allocatable :: flat_result(:)
        type(bracket_root_result) :: res(size(xl0, 1), size(xl0, 2), size(xl0, 3), size(xl0, 4), &
            size(xl0, 5), size(xl0, 6), size(xl0, 7), size(xl0, 8), size(xl0, 9), size(xl0, 10), &
            size(xl0, 11), size(xl0, 12), size(xl0, 13))

        call bracket_root_batch(f, reshape(xl0, [size(xl0)]), &
            shape(res), xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_root_rank13

    function bracket_root_rank14(f, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 14-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: xl0(:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_root_result), allocatable :: flat_result(:)
        type(bracket_root_result) :: res(size(xl0, 1), size(xl0, 2), size(xl0, 3), size(xl0, 4), &
            size(xl0, 5), size(xl0, 6), size(xl0, 7), size(xl0, 8), size(xl0, 9), size(xl0, 10), &
            size(xl0, 11), size(xl0, 12), size(xl0, 13), size(xl0, 14))

        call bracket_root_batch(f, reshape(xl0, [size(xl0)]), &
            shape(res), xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_root_rank14

    function bracket_root_rank15(f, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Shape-preserving independent 15-dimensional batch, via the shared flat engine.
        procedure(scalar_root_function) :: f !! common scalar function
        real(dp), intent(in) :: xl0(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:) !! problem starting points
        real(dp), intent(in), optional :: xr0(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmin(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: xmax(..) !! trailing-axis broadcast parameter
        real(dp), intent(in), optional :: factor !! common expansion factor
        integer, intent(in), optional :: maxiter !! common iteration limit
        type(bracket_root_result), allocatable :: flat_result(:)
        type(bracket_root_result) :: res(size(xl0, 1), size(xl0, 2), size(xl0, 3), size(xl0, 4), &
            size(xl0, 5), size(xl0, 6), size(xl0, 7), size(xl0, 8), size(xl0, 9), size(xl0, 10), &
            size(xl0, 11), size(xl0, 12), size(xl0, 13), size(xl0, 14), size(xl0, 15))

        call bracket_root_batch(f, reshape(xl0, [size(xl0)]), &
            shape(res), xr0=xr0, xmin=xmin, xmax=xmax, factor=factor, maxiter=maxiter, res=flat_result)
        res = reshape(flat_result, shape(res))
    end function bracket_root_rank15

    function find_root_scalar(f, bracket, xatol, xrtol, fatol, frtol, maxiter, callback) result(res)
        ! Scalar form of SciPy elementwise.find_root, using Chandrupatla's method.
        ! A two-point bracket with opposite function signs is required.
        procedure(scalar_root_function) :: f !! function of one real variable
        real(dp), intent(in) :: bracket(:) !! two initial endpoints (either order)
        real(dp), intent(in), optional :: xatol !! absolute abscissa tolerance
        real(dp), intent(in), optional :: xrtol !! relative abscissa tolerance
        real(dp), intent(in), optional :: fatol !! absolute function tolerance
        real(dp), intent(in), optional :: frtol !! relative function tolerance
        integer, intent(in), optional :: maxiter !! maximum refinement iterations
        procedure(find_root_callback), optional :: callback !! optional stop hook
        type(find_root_result) :: res

        real(dp) :: x1, x2, x3, f1, f2, f3, xnew, fnew
        real(dp) :: atol_x, rtol_x, atol_f, rtol_f, ftarget
        real(dp) :: t, tl, tol, dx, xi1, phi1, alpha, next_t
        real(dp) :: denom, fscale, p1, p2
        integer :: i, niter
        logical :: stopped

        if (size(bracket) /= 2) return
        if (.not. all(ieee_is_finite(bracket))) return
        atol_x = 4.0_dp * tiny(1.0_dp)
        rtol_x = 4.0_dp * epsilon(1.0_dp)
        atol_f = tiny(1.0_dp)
        rtol_f = 0.0_dp
        if (present(xatol)) atol_x = xatol
        if (present(xrtol)) rtol_x = xrtol
        if (present(fatol)) atol_f = fatol
        if (present(frtol)) rtol_f = frtol
        ! Check tolerances separately before arithmetic: an optimized MIN over
        ! a NaN may raise IEEE_INVALID when exception trapping is enabled.
        if (.not. ieee_is_finite(atol_x)) return
        if (.not. ieee_is_finite(rtol_x)) return
        if (.not. ieee_is_finite(atol_f)) return
        if (.not. ieee_is_finite(rtol_f)) return
        if (atol_x < 0.0_dp) return
        if (rtol_x < 0.0_dp) return
        if (atol_f < 0.0_dp) return
        if (rtol_f < 0.0_dp) return
        niter = 2046
        if (present(maxiter)) niter = maxiter
        if (niter < 0) return

        x1 = bracket(1)
        x2 = bracket(2)
        f1 = f(x1)
        f2 = f(x2)
        res%nfev = 2
        call store_find_root_state(res, x1, x2, f1, f2)
        if (.not. ieee_is_finite(f1) .or. .not. ieee_is_finite(f2)) then
            res%status = root_status_value_error
            call invalidate_find_root(res)
            return
        end if

        ! SciPy scales frtol with the smaller initial absolute function value.
        ftarget = 0.0_dp
        if (rtol_f > 0.0_dp) then
            fscale = min(abs(f1), abs(f2))
            if (rtol_f > 1.0_dp) then
                if (fscale > huge(1.0_dp)/rtol_f) then
                    ftarget = huge(1.0_dp)
                else
                    ftarget = rtol_f * fscale
                end if
            else
                ftarget = rtol_f * fscale
            end if
        end if
        t = 0.5_dp
        call check_find_root(res, x1, x2, f1, f2, atol_x, rtol_x, atol_f, ftarget, tol, dx)
        if (present(callback)) then
            stopped = .false.
            call callback(res, stopped)
            if (stopped .and. res%status == 1) res%status = find_root_status_callback
        end if

        do i = 1, niter
            if (res%status /= 1) exit
            ! Use the current bracket, preserving the sign-changing endpoints.
            xnew = x1 + t * (x2 - x1)
            if (.not. ieee_is_finite(xnew)) then
                res%status = root_status_value_error
                call invalidate_find_root(res)
                exit
            end if
            fnew = f(xnew)
            res%nfev = res%nfev + 1
            res%nit = i
            x3 = x2
            f3 = f2
            if (ieee_is_finite(fnew)) then
                if (same_root_sign(fnew, f1)) then
                    x3 = x1
                    f3 = f1
                else
                    x2 = x1
                    f2 = f1
                end if
            end if
            x1 = xnew
            f1 = fnew
            if (.not. ieee_is_finite(fnew)) then
                res%status = root_status_value_error
                call store_find_root_state(res, x1, x2, f1, f2)
                call invalidate_find_root(res)
                exit
            end if
            call check_find_root(res, x1, x2, f1, f2, atol_x, rtol_x, atol_f, ftarget, tol, dx)
            if (present(callback)) then
                stopped = .false.
                call callback(res, stopped)
                if (stopped .and. res%status == 1) res%status = find_root_status_callback
            end if
            if (res%status /= 1) exit

            ! Chandrupatla's inverse-quadratic interpolation criterion.
            ! Bisection is the safe fallback for degenerate triples.
            next_t = 0.5_dp
            if (x3 /= x2 .and. f3 /= f2 .and. x2 /= x1) then
                xi1 = (x1 - x2) / (x3 - x2)
                if (xi1 > 0.0_dp .and. xi1 < 1.0_dp) then
                    phi1 = (f1 - f2) / (f3 - f2)
                    if (phi1 > 1.0_dp - sqrt(1.0_dp - xi1) .and. phi1 < sqrt(xi1)) then
                        alpha = (x3 - x1) / (x2 - x1)
                        denom = f3 - f1
                        if (f1 /= f2 .and. f3 /= f2 .and. denom /= 0.0_dp .and. f2 /= f3) then
                            p1 = (f1 / (f1-f2)) * (f3 / (f3-f2))
                            p2 = alpha * (f1 / denom) * (f2 / (f2-f3))
                            if (ieee_is_finite(p1) .and. ieee_is_finite(p2)) next_t = p1-p2
                        end if
                    end if
                end if
            end if
            if (dx > 0.0_dp) then
                tl = 0.5_dp * tol / dx
                t = max(tl, min(1.0_dp-tl, next_t))
            else
                t = 0.5_dp
            end if
        end do

        if (res%status == 1) res%status = root_status_convergence_error
        res%success = (res%status == root_status_success)
    end function find_root_scalar

    subroutine check_find_root(res, x1, x2, f1, f2, atol_x, rtol_x, atol_f, ftarget, tol, dx)
        ! Preserve SciPy's function-first, bracket-second stopping semantics.
        type(find_root_result), intent(inout) :: res !! iteration result
        real(dp), intent(in) :: x1 !! first current endpoint
        real(dp), intent(in) :: x2 !! second current endpoint
        real(dp), intent(in) :: f1 !! function value at x1
        real(dp), intent(in) :: f2 !! function value at x2
        real(dp), intent(in) :: atol_x !! absolute root tolerance
        real(dp), intent(in) :: rtol_x !! relative root tolerance
        real(dp), intent(in) :: atol_f !! absolute function tolerance
        real(dp), intent(in) :: ftarget !! scaled relative function tolerance
        real(dp), intent(out) :: tol !! current root-width tolerance
        real(dp), intent(out) :: dx !! current bracket width
        real(dp) :: best

        call store_find_root_state(res, x1, x2, f1, f2)
        res%status = 1
        best = abs(res%x)
        tol = atol_x + best * rtol_x
        dx = abs(x2-x1)
        if (abs(res%f_x) <= atol_f + ftarget) then
            res%status = root_status_success
        else if (same_root_sign(f1, f2)) then
            res%status = root_status_sign_error
            call invalidate_find_root(res)
        else if (dx < tol) then
            res%status = root_status_success
        end if
        res%success = (res%status == root_status_success)
    end subroutine check_find_root

    subroutine store_find_root_state(res, x1, x2, f1, f2)
        ! Sort the reported bracket and retain the smaller-residual endpoint.
        type(find_root_result), intent(inout) :: res !! iteration result
        real(dp), intent(in) :: x1 !! first current endpoint
        real(dp), intent(in) :: x2 !! second current endpoint
        real(dp), intent(in) :: f1 !! function value at x1
        real(dp), intent(in) :: f2 !! function value at x2

        if (x1 < x2) then
            res%xl = x1
            res%xr = x2
            res%fl = f1
            res%fr = f2
        else
            res%xl = x2
            res%xr = x1
            res%fl = f2
            res%fr = f1
        end if
        if (.not. ieee_is_finite(f1) .or. .not. ieee_is_finite(f2)) then
            call invalidate_find_root(res)
            return
        end if
        if (abs(f1) < abs(f2)) then
            res%x = x1
            res%f_x = f1
        else
            res%x = x2
            res%f_x = f2
        end if
    end subroutine store_find_root_state

    subroutine invalidate_find_root(res)
        ! Report invalid/nonfinite root estimates without losing bracket data.
        type(find_root_result), intent(inout) :: res !! iteration result
        res%x = ieee_value(0.0_dp, ieee_quiet_nan)
        res%f_x = ieee_value(0.0_dp, ieee_quiet_nan)
        res%success = .false.
    end subroutine invalidate_find_root

    pure logical function same_root_sign(a, b)
        ! Compare signs without multiplying (safe under IEEE overflow traps).
        real(dp), intent(in) :: a !! first finite function value
        real(dp), intent(in) :: b !! second finite function value
        same_root_sign = (a > 0.0_dp .and. b > 0.0_dp) .or. &
            (a < 0.0_dp .and. b < 0.0_dp) .or. &
            (a == 0.0_dp .and. b == 0.0_dp)
    end function same_root_sign

    function bracket_root_scalar(f, xl0, xr0, xmin, xmax, factor, maxiter) result(res)
        ! Find a sign-changing interval by expanding to the left and right.
        ! This is the scalar counterpart of scipy.optimize.elementwise.bracket_root.
        procedure(scalar_root_function) :: f !! scalar real-valued function
        real(dp), intent(in) :: xl0 !! initial left search endpoint
        real(dp), intent(in), optional :: xr0 !! initial right endpoint, default xl0 + 1
        real(dp), intent(in), optional :: xmin !! allowable lower search bound
        real(dp), intent(in), optional :: xmax !! allowable upper search bound
        real(dp), intent(in), optional :: factor !! expansion factor, default 2, must exceed 1
        integer, intent(in), optional :: maxiter !! maximum expansion iterations, default 1000
        type(bracket_root_result) :: res

        real(dp) :: x(2)
        real(dp) :: fx(2)
        real(dp) :: xlast(2)
        real(dp) :: flast(2)
        real(dp) :: anchor(2)
        real(dp) :: bound(2)
        real(dp) :: distance(2)
        real(dp) :: xnew
        real(dp) :: grow
        real(dp) :: space
        real(dp) :: xr
        integer :: i
        integer :: j
        integer :: limit_iter
        logical :: bounded_side(2)
        logical :: active(2)
        logical :: found(2)
        logical :: nonfinite

        res%xl = xl0
        res%xr = xl0
        if (.not. ieee_is_finite(xl0)) return
        if (present(xr0)) then
            xr = xr0
        else
            ! SciPy shrinks its default initial gap if xmax is supplied.
            if (present(xmax)) then
                if (ieee_is_nan(xmax)) return
                if (xmax <= xl0) return
                if (ieee_is_finite(xmax)) then
                    ! Divide first to avoid overflow of xmax - xl0.
                    if (xl0 < 0.0_dp .and. xmax > 0.0_dp) then
                        space = xmax / 8.0_dp - xl0 / 8.0_dp
                    else
                        space = (xmax - xl0) / 8.0_dp
                    end if
                    xr = xl0 + min(space, 1.0_dp)
                else
                    if (xl0 > huge(1.0_dp) - 1.0_dp) return
                    xr = xl0 + 1.0_dp
                end if
            else
                if (xl0 > huge(1.0_dp) - 1.0_dp) return
                xr = xl0 + 1.0_dp
            end if
        end if
        res%xr = xr
        if (.not. ieee_is_finite(xr)) return
        if (xl0 >= xr) return
        grow = 2.0_dp
        if (present(factor)) grow = factor
        if (.not. ieee_is_finite(grow)) return
        if (grow <= 1.0_dp) return
        limit_iter = 1000
        if (present(maxiter)) limit_iter = maxiter
        if (limit_iter < 0) return
        bounded_side = [present(xmin), present(xmax)]
        bound = 0.0_dp
        if (present(xmin)) then
            if (ieee_is_nan(xmin)) return
            if (xmin > xl0) return
            if (ieee_is_finite(xmin)) then
                bound(1) = xmin
            else
                bounded_side(1) = .false.  ! Negative infinity is unbounded.
            end if
        end if
        if (present(xmax)) then
            if (ieee_is_nan(xmax)) return
            if (xmax < xr) return
            if (ieee_is_finite(xmax)) then
                bound(2) = xmax
            else
                bounded_side(2) = .false.  ! Positive infinity is unbounded.
            end if
        end if

        x = [xl0, xr]
        anchor = [xr, xl0]
        xlast = [xr, xl0]
        fx(1) = f(x(1))
        fx(2) = f(x(2))
        res%nfev = 2
        res%fl = fx(1)
        res%fr = fx(2)
        if (.not. ieee_is_finite(fx(1)) .or. .not. ieee_is_finite(fx(2))) then
            res%status = root_status_value_error
            return
        end if
        if (bracket_has_root(fx(1), fx(2))) then
            res%success = .true.
            res%status = root_status_success
            return
        end if
        ! The initial endpoints may straddle the root even when their
        ! separation would overflow in floating-point arithmetic.
        if (xl0 < 0.0_dp .and. xr > 0.0_dp) then
            if (xr > huge(1.0_dp) + xl0) then
                res%status = root_status_value_error
                return
            end if
        end if
        distance = [xl0 - xr, xr - xl0]
        if (bounded_side(1)) distance(1) = bound(1) - xl0
        if (bounded_side(2)) distance(2) = bound(2) - xr
        flast = [fx(2), fx(1)]
        active = .true.
        if (bounded_side(1)) active(1) = xl0 > bound(1)
        if (bounded_side(2)) active(2) = xr < bound(2)
        found = .false.
        nonfinite = .false.
        if (.not. any(active)) then
            res%status = bracket_root_status_limit
            return
        end if

        do i = 1, limit_iter
            do j = 1, 2
                if (.not. active(j)) cycle
                if (bounded_side(j)) then
                    distance(j) = distance(j) / grow
                    xnew = bound(j) - distance(j)
                else
                    ! Avoid overflows with floating-point exception traps enabled.
                    if (abs(distance(j)) > huge(1.0_dp) / grow) then
                        nonfinite = .true.
                        active(j) = .false.
                        cycle
                    end if
                    distance(j) = distance(j) * grow
                    if (abs(anchor(j)) > huge(1.0_dp) - abs(distance(j))) then
                        nonfinite = .true.
                        active(j) = .false.
                        cycle
                    end if
                    xnew = anchor(j) + distance(j)
                end if
                if (.not. ieee_is_finite(xnew)) then
                    nonfinite = .true.
                    active(j) = .false.
                    cycle
                end if
                if (xnew == x(j)) then
                    active(j) = .false.
                    cycle
                end if
                xlast(j) = x(j)
                flast(j) = fx(j)
                x(j) = xnew
                fx(j) = f(xnew)
                res%nfev = res%nfev + 1
                if (.not. ieee_is_finite(fx(j))) then
                    nonfinite = .true.
                    active(j) = .false.
                    cycle
                end if
                if (bracket_has_root(fx(j), flast(j))) then
                    found(j) = .true.
                    active(j) = .false.
                else if (bounded_side(j)) then
                    if (x(j) == bound(j)) active(j) = .false.
                end if
            end do
            res%nit = i
            if (any(found)) exit
            if (.not. any(active)) exit
        end do

        res%xl = x(1)
        res%fl = fx(1)
        res%xr = x(2)
        res%fr = fx(2)
        if (any(found)) then
            ! Select the shorter bracket when both searches succeed together.
            if (found(1) .and. found(2)) then
                if (xlast(1) - x(1) <= x(2) - xlast(2)) then
                    res%xr = xlast(1)
                    res%fr = flast(1)
                else
                    res%xl = xlast(2)
                    res%fl = flast(2)
                end if
            else if (found(1)) then
                res%xr = xlast(1)
                res%fr = flast(1)
            else
                res%xl = xlast(2)
                res%fl = flast(2)
            end if
            res%status = root_status_success
            res%success = .true.
        else if (nonfinite) then
            res%status = root_status_value_error
        else if (.not. any(active)) then
            res%status = bracket_root_status_limit
        else
            res%status = root_status_convergence_error
        end if
    end function bracket_root_scalar

    pure logical function bracket_has_root(fa, fb) result(ok)
        real(dp), intent(in) :: fa !! first finite function value
        real(dp), intent(in) :: fb !! second finite function value
        ok = fa == 0.0_dp .or. fb == 0.0_dp .or. opposite_signs(fa, fb)
    end function bracket_has_root

    function root_scalar(f, method, bracket, x0, x1, fprime, fprime2, xtol, rtol, &
            maxiter, k) result(res)
        procedure(scalar_root_function) :: f !! scalar objective function
        character(len=*), intent(in), optional :: method !! solver method
        real(dp), intent(in), optional :: bracket(:) !! two-element sign-changing interval
        real(dp), intent(in), optional :: x0 !! initial point for open methods
        real(dp), intent(in), optional :: x1 !! second point for the secant method
        procedure(scalar_root_function), optional :: fprime !! first derivative
        procedure(scalar_root_function), optional :: fprime2 !! second derivative
        real(dp), intent(in), optional :: xtol !! absolute tolerance
        real(dp), intent(in), optional :: rtol !! relative tolerance
        integer, intent(in), optional :: maxiter !! maximum iterations
        integer, intent(in), optional :: k !! TOMS748 Newton-quadratic steps
        type(root_result) :: res

        character(len=16) :: meth

        call init_result(res, '')
        meth = ''
        if (present(method)) meth = lowercase(trim(method))

        if (len_trim(meth) == 0) then
            if (present(bracket)) then
                meth = 'brentq'
            else if (present(x0)) then
                if (present(fprime)) then
                    if (present(fprime2)) then
                        meth = 'halley'
                    else
                        meth = 'newton'
                    end if
                else if (present(x1)) then
                    meth = 'secant'
                else
                    meth = 'newton'
                end if
            else
                return
            end if
        end if

        select case (trim(meth))
        case ('bisect')
            if (.not. present(bracket)) return
            if (size(bracket) < 2) return
            res = bisect(f, bracket(1), bracket(2), xtol, rtol, maxiter)
        case ('brentq')
            if (.not. present(bracket)) return
            if (size(bracket) < 2) return
            res = brentq(f, bracket(1), bracket(2), xtol, rtol, maxiter)
        case ('brenth')
            if (.not. present(bracket)) return
            if (size(bracket) < 2) return
            res = brenth(f, bracket(1), bracket(2), xtol, rtol, maxiter)
        case ('ridder')
            if (.not. present(bracket)) return
            if (size(bracket) < 2) return
            res = ridder(f, bracket(1), bracket(2), xtol, rtol, maxiter)
        case ('toms748')
            if (.not. present(bracket)) return
            if (size(bracket) < 2) return
            res = toms748(f, bracket(1), bracket(2), xtol, rtol, maxiter, k)
        case ('secant')
            if (.not. present(x0)) return
            if (present(x1)) then
                res = secant_solver(f, x0, x1=x1, tol=xtol, rtol=rtol, maxiter=maxiter)
            else
                res = secant_solver(f, x0, tol=xtol, rtol=rtol, maxiter=maxiter)
            end if
        case ('newton')
            if (.not. present(x0)) return
            if (present(fprime)) then
                res = newton_solver(f, x0, fprime=fprime, tol=xtol, rtol=rtol, maxiter=maxiter)
            else
                res = numerical_newton_solver(f, x0, tol=xtol, rtol=rtol, maxiter=maxiter)
            end if
        case ('halley')
            if (.not. present(x0)) return
            if (.not. present(fprime)) return
            if (.not. present(fprime2)) return
            res = newton_solver(f, x0, fprime=fprime, fprime2=fprime2, tol=xtol, &
                rtol=rtol, maxiter=maxiter)
        case default
            return
        end select
    end function root_scalar

    function bisect(f, a, b, xtol, rtol, maxiter) result(res)
        procedure(scalar_root_function) :: f !! objective function
        real(dp), intent(in) :: a !! lower endpoint
        real(dp), intent(in) :: b !! upper endpoint
        real(dp), intent(in), optional :: xtol !! absolute tolerance
        real(dp), intent(in), optional :: rtol !! relative tolerance
        integer, intent(in), optional :: maxiter !! maximum iterations
        type(root_result) :: res

        integer :: i
        integer :: mit
        real(dp) :: aa
        real(dp) :: bb
        real(dp) :: fa
        real(dp) :: fb
        real(dp) :: fc
        real(dp) :: tol_abs
        real(dp) :: tol_rel
        real(dp) :: c

        call init_result(res, 'bisect')
        if (.not. bracket_options(xtol, rtol, maxiter, tol_abs, tol_rel, mit)) return
        if (.not. valid_ab(a, b)) return

        aa = a
        bb = b
        fa = f(aa)
        fb = f(bb)
        res%function_calls = 2
        if (.not. finite_value(fa) .or. .not. finite_value(fb)) then
            call fail_result(res, root_status_value_error)
            res%root = merge(aa, bb, .not. finite_value(fa))
            return
        end if
        if (near_zero(fa)) then
            call succeed_result(res, aa, 0)
            return
        end if
        if (near_zero(fb)) then
            call succeed_result(res, bb, 0)
            return
        end if
        if (.not. opposite_signs(fa, fb)) then
            call fail_result(res, root_status_sign_error)
            res%root = aa
            return
        end if

        do i = 1, mit
            c = aa + 0.5_dp * (bb - aa)
            fc = f(c)
            res%function_calls = res%function_calls + 1
            res%iterations = i
            if (.not. finite_value(fc)) then
                res%root = c
                call fail_result(res, root_status_value_error)
                return
            end if
            if (near_zero(fc)) then
                call succeed_result(res, c, i)
                return
            end if
            if (isclose(aa, bb, tol_abs, tol_rel)) then
                call succeed_result(res, c, i)
                return
            end if
            if (opposite_signs(fc, fb)) then
                aa = c
                fa = fc
            else
                bb = c
                fb = fc
            end if
            if (isclose(aa, bb, tol_abs, tol_rel)) then
                if (abs(fa) <= abs(fb)) then
                    call succeed_result(res, aa, i)
                else
                    call succeed_result(res, bb, i)
                end if
                return
            end if
        end do

        res%root = aa + 0.5_dp * (bb - aa)
        call fail_result(res, root_status_convergence_error)
    end function bisect

    function brentq(f, a, b, xtol, rtol, maxiter) result(res)
        procedure(scalar_root_function) :: f !! objective function
        real(dp), intent(in) :: a !! lower endpoint
        real(dp), intent(in) :: b !! upper endpoint
        real(dp), intent(in), optional :: xtol !! absolute tolerance
        real(dp), intent(in), optional :: rtol !! relative tolerance
        integer, intent(in), optional :: maxiter !! maximum iterations
        type(root_result) :: res

        res = brent_kernel(f, a, b, xtol, rtol, maxiter, .false.)
        res%method = 'brentq'
    end function brentq

    function brenth(f, a, b, xtol, rtol, maxiter) result(res)
        procedure(scalar_root_function) :: f !! objective function
        real(dp), intent(in) :: a !! lower endpoint
        real(dp), intent(in) :: b !! upper endpoint
        real(dp), intent(in), optional :: xtol !! absolute tolerance
        real(dp), intent(in), optional :: rtol !! relative tolerance
        integer, intent(in), optional :: maxiter !! maximum iterations
        type(root_result) :: res

        res = brent_kernel(f, a, b, xtol, rtol, maxiter, .true.)
        res%method = 'brenth'
    end function brenth

    function brent_kernel(f, a, b, xtol, rtol, maxiter, hyperbolic) result(res)
        procedure(scalar_root_function) :: f
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp), intent(in), optional :: xtol
        real(dp), intent(in), optional :: rtol
        integer, intent(in), optional :: maxiter
        logical, intent(in) :: hyperbolic
        type(root_result) :: res

        integer :: i
        integer :: mit
        real(dp) :: dblk
        real(dp) :: delta
        real(dp) :: dpre
        real(dp) :: fblk
        real(dp) :: fcur
        real(dp) :: fpre
        real(dp) :: sbis
        real(dp) :: scur
        real(dp) :: spre
        real(dp) :: stry
        real(dp) :: tmpf
        real(dp) :: tmpx
        real(dp) :: tol_abs
        real(dp) :: tol_rel
        real(dp) :: xblk
        real(dp) :: xcur
        real(dp) :: xdelta
        real(dp) :: xpre

        call init_result(res, 'brentq')
        if (hyperbolic) res%method = 'brenth'
        if (.not. bracket_options(xtol, rtol, maxiter, tol_abs, tol_rel, mit)) return
        if (.not. valid_ab(a, b)) return

        xpre = a
        xcur = b
        fpre = f(xpre)
        fcur = f(xcur)
        res%function_calls = 2
        if (.not. finite_value(fpre) .or. .not. finite_value(fcur)) then
            res%root = merge(xpre, xcur, .not. finite_value(fpre))
            call fail_result(res, root_status_value_error)
            return
        end if
        if (near_zero(fpre)) then
            call succeed_result(res, xpre, 0)
            return
        end if
        if (near_zero(fcur)) then
            call succeed_result(res, xcur, 0)
            return
        end if
        if (.not. opposite_signs(fpre, fcur)) then
            res%root = xpre
            call fail_result(res, root_status_sign_error)
            return
        end if

        xblk = xpre
        fblk = fpre
        scur = xcur - xpre
        spre = scur

        do i = 1, mit
            res%iterations = i
            if (opposite_signs(fpre, fcur)) then
                xblk = xpre
                fblk = fpre
                scur = xcur - xpre
                spre = scur
            end if
            if (abs(fblk) < abs(fcur)) then
                tmpx = xcur
                xcur = xblk
                xpre = tmpx
                xblk = tmpx
                tmpf = fcur
                fcur = fblk
                fpre = tmpf
                fblk = tmpf
            end if

            delta = 0.5_dp * (tol_abs + tol_rel * abs(xcur))
            sbis = 0.5_dp * (xblk - xcur)
            if (abs(sbis) < delta) then
                call succeed_result(res, xcur, i)
                return
            end if

            if (abs(spre) > delta .and. abs(fcur) < abs(fpre)) then
                if (near_equal(xpre, xblk)) then
                    stry = -fcur * (xcur - xpre) / (fcur - fpre)
                else
                    dpre = (fpre - fcur) / (xpre - xcur)
                    dblk = (fblk - fcur) / (xblk - xcur)
                    if (hyperbolic) then
                        stry = -fcur * (fblk - fpre) / (fblk * dpre - fpre * dblk)
                    else
                        stry = -fcur * (fblk * dblk - fpre * dpre) / &
                            (dblk * dpre * (fblk - fpre))
                    end if
                end if
                if (2.0_dp * abs(stry) < min(abs(spre), 3.0_dp * abs(sbis) - delta)) then
                    spre = scur
                    scur = stry
                else
                    spre = sbis
                    scur = sbis
                end if
            else
                spre = sbis
                scur = sbis
            end if

            xpre = xcur
            fpre = fcur
            if (abs(scur) > delta) then
                xcur = xcur + scur
            else
                xdelta = sign(delta, sbis)
                xcur = xcur + xdelta
            end if
            fcur = f(xcur)
            res%function_calls = res%function_calls + 1
            if (.not. finite_value(fcur)) then
                res%root = xcur
                call fail_result(res, root_status_value_error)
                return
            end if
            if (near_zero(fcur)) then
                call succeed_result(res, xcur, i + 1)
                return
            end if
        end do

        res%root = xcur
        call fail_result(res, root_status_convergence_error)
    end function brent_kernel

    function ridder(f, a, b, xtol, rtol, maxiter) result(res)
        procedure(scalar_root_function) :: f !! objective function
        real(dp), intent(in) :: a !! lower endpoint
        real(dp), intent(in) :: b !! upper endpoint
        real(dp), intent(in), optional :: xtol !! absolute tolerance
        real(dp), intent(in), optional :: rtol !! relative tolerance
        integer, intent(in), optional :: maxiter !! maximum iterations
        type(root_result) :: res

        integer :: i
        integer :: mit
        real(dp) :: denom
        real(dp) :: disc
        real(dp) :: fh
        real(dp) :: fl
        real(dp) :: fm
        real(dp) :: fnew
        real(dp) :: tol_abs
        real(dp) :: tol_rel
        real(dp) :: xh
        real(dp) :: xl
        real(dp) :: xm
        real(dp) :: xnew
        real(dp) :: xold

        call init_result(res, 'ridder')
        if (.not. bracket_options(xtol, rtol, maxiter, tol_abs, tol_rel, mit)) return
        if (.not. valid_ab(a, b)) return

        xl = a
        xh = b
        fl = f(xl)
        fh = f(xh)
        res%function_calls = 2
        if (.not. finite_value(fl) .or. .not. finite_value(fh)) then
            res%root = merge(xl, xh, .not. finite_value(fl))
            call fail_result(res, root_status_value_error)
            return
        end if
        if (near_zero(fl)) then
            call succeed_result(res, xl, 0)
            return
        end if
        if (near_zero(fh)) then
            call succeed_result(res, xh, 0)
            return
        end if
        if (.not. opposite_signs(fl, fh)) then
            res%root = xl
            call fail_result(res, root_status_sign_error)
            return
        end if

        xold = huge(1.0_dp)
        do i = 1, mit
            res%iterations = i
            xm = xl + 0.5_dp * (xh - xl)
            fm = f(xm)
            res%function_calls = res%function_calls + 1
            if (.not. finite_value(fm)) then
                res%root = xm
                call fail_result(res, root_status_value_error)
                return
            end if

            disc = fm * fm - fl * fh
            if (disc < 0.0_dp) disc = 0.0_dp
            denom = sqrt(disc)
            if (near_zero(denom)) then
                res%root = xm
                call fail_result(res, root_status_convergence_error)
                return
            end if

            xnew = xm + (xm - xl) * sign(1.0_dp, fl - fh) * fm / denom
            fnew = f(xnew)
            res%function_calls = res%function_calls + 1
            if (.not. finite_value(fnew)) then
                res%root = xnew
                call fail_result(res, root_status_value_error)
                return
            end if
            if (near_zero(fnew)) then
                call succeed_result(res, xnew, i)
                return
            end if
            if (isclose(xold, xnew, tol_abs, tol_rel)) then
                call succeed_result(res, xnew, i)
                return
            end if
            xold = xnew

            if (opposite_signs(fm, fnew)) then
                xl = xm
                fl = fm
                xh = xnew
                fh = fnew
            else if (opposite_signs(fl, fnew)) then
                xh = xnew
                fh = fnew
            else if (opposite_signs(fh, fnew)) then
                xl = xnew
                fl = fnew
            else
                res%root = xnew
                call fail_result(res, root_status_convergence_error)
                return
            end if

            if (isclose(xl, xh, tol_abs, tol_rel)) then
                call succeed_result(res, xnew, i)
                return
            end if
        end do

        res%root = xold
        call fail_result(res, root_status_convergence_error)
    end function ridder

    function toms748(f, a, b, xtol, rtol, maxiter, k) result(res)
        procedure(scalar_root_function) :: f !! objective function
        real(dp), intent(in) :: a !! lower endpoint
        real(dp), intent(in) :: b !! upper endpoint
        real(dp), intent(in), optional :: xtol !! absolute tolerance
        real(dp), intent(in), optional :: rtol !! relative tolerance
        integer, intent(in), optional :: maxiter !! maximum iterations
        integer, intent(in), optional :: k !! Newton-quadratic steps per substep
        type(root_result) :: res

        integer :: kk
        integer :: mit
        integer :: nsteps
        real(dp) :: aa
        real(dp) :: bb
        real(dp) :: c
        real(dp) :: c0
        real(dp) :: d
        real(dp) :: e
        real(dp) :: fa
        real(dp) :: fb
        real(dp) :: fc
        real(dp) :: fd
        real(dp) :: fe
        real(dp) :: fu
        real(dp) :: old_width
        real(dp) :: sec
        real(dp) :: tol_abs
        real(dp) :: tol_rel
        real(dp) :: u
        real(dp) :: width
        real(dp) :: z
        real(dp) :: fz
        logical :: have_e
        logical :: interp_ok

        call init_result(res, 'toms748')
        if (.not. bracket_options(xtol, rtol, maxiter, tol_abs, tol_rel, mit)) return
        if (.not. valid_ab(a, b)) return
        kk = 1
        if (present(k)) kk = max(1, min(100, k))

        aa = a
        bb = b
        fa = f(aa)
        res%function_calls = 1
        if (.not. finite_value(fa)) then
            res%root = aa
            call fail_result(res, root_status_value_error)
            return
        end if
        if (near_zero(fa)) then
            call succeed_result(res, aa, 0)
            return
        end if
        fb = f(bb)
        res%function_calls = 2
        if (.not. finite_value(fb)) then
            res%root = bb
            call fail_result(res, root_status_value_error)
            return
        end if
        if (near_zero(fb)) then
            call succeed_result(res, bb, 0)
            return
        end if
        if (.not. opposite_signs(fa, fb)) then
            res%root = aa
            call fail_result(res, root_status_sign_error)
            return
        end if
        if (isclose(aa, bb, tol_abs, tol_rel)) then
            call succeed_result(res, aa + 0.5_dp * (bb - aa), 0)
            return
        end if

        ! First secant step.
        c = safe_secant(aa, bb, fa, fb)
        if (.not. inside_open(c, aa, bb)) c = aa + 0.5_dp * (bb - aa)
        fc = f(c)
        res%function_calls = res%function_calls + 1
        if (.not. finite_value(fc)) then
            res%root = c
            call fail_result(res, root_status_value_error)
            return
        end if
        if (near_zero(fc)) then
            call succeed_result(res, c, 0)
            return
        end if
        call update_bracket(aa, bb, fa, fb, c, fc, d, fd)
        e = 0.0_dp
        fe = 0.0_dp
        have_e = .false.
        res%iterations = 1

        do while (res%iterations < mit)
            if (isclose(aa, bb, tol_abs, tol_rel)) then
                call succeed_result(res, aa + 0.5_dp * (bb - aa), res%iterations)
                return
            end if

            res%iterations = res%iterations + 1
            old_width = bb - aa
            c = quiet_nan()

            do nsteps = 2, kk + 1
                interp_ok = .false.
                if (have_e) then
                    if (four_values_separated(fa, fb, fd, fe)) then
                        c0 = inverse_cubic_zero(aa, bb, d, e, fa, fb, fd, fe)
                        if (inside_open(c0, aa, bb)) then
                            c = c0
                            interp_ok = .true.
                        end if
                    end if
                end if
                if (.not. interp_ok) c = newton_quadratic(aa, bb, fa, fb, d, fd, nsteps)
                if (.not. inside_open(c, aa, bb)) c = aa + 0.5_dp * (bb - aa)

                fc = f(c)
                res%function_calls = res%function_calls + 1
                if (.not. finite_value(fc)) then
                    res%root = c
                    call fail_result(res, root_status_value_error)
                    return
                end if
                if (near_zero(fc)) then
                    call succeed_result(res, c, res%iterations)
                    return
                end if
                e = d
                fe = fd
                have_e = .true.
                call update_bracket(aa, bb, fa, fb, c, fc, d, fd)
            end do

            if (abs(fa) < abs(fb)) then
                u = aa
                fu = fa
            else
                u = bb
                fu = fb
            end if
            sec = divided_difference_two(aa, bb, fa, fb)
            if (near_zero(sec)) then
                c = aa + 0.5_dp * (bb - aa)
            else
                c = u - 2.0_dp * fu / sec
            end if
            if (abs(c - u) > 0.5_dp * (bb - aa)) c = aa + 0.5_dp * (bb - aa)
            if (.not. inside_open(c, aa, bb)) c = aa + 0.5_dp * (bb - aa)

            fc = f(c)
            res%function_calls = res%function_calls + 1
            if (.not. finite_value(fc)) then
                res%root = c
                call fail_result(res, root_status_value_error)
                return
            end if
            if (near_zero(fc)) then
                call succeed_result(res, c, res%iterations)
                return
            end if
            e = d
            fe = fd
            have_e = .true.
            call update_bracket(aa, bb, fa, fb, c, fc, d, fd)

            width = bb - aa
            if (width > 0.5_dp * old_width) then
                z = aa + 0.5_dp * (bb - aa)
                fz = f(z)
                res%function_calls = res%function_calls + 1
                if (.not. finite_value(fz)) then
                    res%root = z
                    call fail_result(res, root_status_value_error)
                    return
                end if
                if (near_zero(fz)) then
                    call succeed_result(res, z, res%iterations)
                    return
                end if
                e = d
                fe = fd
                have_e = .true.
                call update_bracket(aa, bb, fa, fb, z, fz, d, fd)
            end if
        end do

        res%root = aa + 0.5_dp * (bb - aa)
        call fail_result(res, root_status_convergence_error)
    end function toms748

    function newton(f, x0, fprime, fprime2, x1, tol, rtol, maxiter) result(res)
        procedure(scalar_root_function) :: f !! objective function
        real(dp), intent(in) :: x0 !! starting point
        procedure(scalar_root_function), optional :: fprime !! first derivative
        procedure(scalar_root_function), optional :: fprime2 !! second derivative
        real(dp), intent(in), optional :: x1 !! second point for secant mode
        real(dp), intent(in), optional :: tol !! absolute step tolerance
        real(dp), intent(in), optional :: rtol !! relative step tolerance
        integer, intent(in), optional :: maxiter !! maximum iterations
        type(root_result) :: res

        if (present(fprime)) then
            if (present(fprime2)) then
                res = newton_solver(f, x0, fprime=fprime, fprime2=fprime2, &
                    tol=tol, rtol=rtol, maxiter=maxiter)
            else
                res = newton_solver(f, x0, fprime=fprime, tol=tol, rtol=rtol, &
                    maxiter=maxiter)
            end if
        else
            if (present(x1)) then
                res = secant_solver(f, x0, x1=x1, tol=tol, rtol=rtol, maxiter=maxiter)
            else
                res = secant_solver(f, x0, tol=tol, rtol=rtol, maxiter=maxiter)
            end if
        end if
    end function newton

    function newton_solver(f, x0, fprime, fprime2, tol, rtol, maxiter) result(res)
        procedure(scalar_root_function) :: f
        real(dp), intent(in) :: x0
        procedure(scalar_root_function) :: fprime
        procedure(scalar_root_function), optional :: fprime2
        real(dp), intent(in), optional :: tol
        real(dp), intent(in), optional :: rtol
        integer, intent(in), optional :: maxiter
        type(root_result) :: res

        integer :: i
        integer :: mit
        real(dp) :: abs_tol
        real(dp) :: adj
        real(dp) :: fd
        real(dp) :: fd2
        real(dp) :: fv
        real(dp) :: p
        real(dp) :: p0
        real(dp) :: rel_tol
        real(dp) :: step

        call init_result(res, 'newton')
        if (.not. newton_options(tol, rtol, maxiter, abs_tol, rel_tol, mit)) return
        if (.not. finite_value(x0)) return

        p0 = x0
        do i = 1, mit
            fv = f(p0)
            res%function_calls = res%function_calls + 1
            if (.not. finite_value(fv)) then
                res%root = p0
                call fail_result(res, root_status_value_error)
                return
            end if
            if (near_zero(fv)) then
                call succeed_result(res, p0, i - 1)
                return
            end if
            fd = fprime(p0)
            res%function_calls = res%function_calls + 1
            if (.not. finite_value(fd)) then
                res%root = p0
                call fail_result(res, root_status_value_error)
                return
            end if
            if (near_zero(fd)) then
                res%root = p0
                res%iterations = i
                call fail_result(res, root_status_convergence_error)
                return
            end if

            step = fv / fd
            if (present(fprime2)) then
                fd2 = fprime2(p0)
                res%function_calls = res%function_calls + 1
                res%method = 'halley'
                if (.not. finite_value(fd2)) then
                    res%root = p0
                    call fail_result(res, root_status_value_error)
                    return
                end if
                adj = 0.5_dp * step * fd2 / fd
                if (abs(adj) < 1.0_dp) step = step / (1.0_dp - adj)
            end if
            p = p0 - step
            res%iterations = i
            if (isclose(p, p0, abs_tol, rel_tol)) then
                call succeed_result(res, p, i)
                return
            end if
            if (.not. finite_value(p)) then
                res%root = p0
                call fail_result(res, root_status_value_error)
                return
            end if
            p0 = p
        end do

        res%root = p0
        call fail_result(res, root_status_convergence_error)
    end function newton_solver

    function numerical_newton_solver(f, x0, tol, rtol, maxiter) result(res)
        procedure(scalar_root_function) :: f
        real(dp), intent(in) :: x0
        real(dp), intent(in), optional :: tol
        real(dp), intent(in), optional :: rtol
        integer, intent(in), optional :: maxiter
        type(root_result) :: res

        integer :: i
        integer :: mit
        integer :: before
        real(dp) :: abs_tol
        real(dp) :: fv
        real(dp) :: p
        real(dp) :: p0
        real(dp) :: rel_tol
        real(dp) :: step
        type(root_result) :: tmp

        call init_result(res, 'newton')
        if (.not. newton_options(tol, rtol, maxiter, abs_tol, rel_tol, mit)) return
        if (.not. finite_value(x0)) return
        p0 = x0

        do i = 1, mit
            fv = f(p0)
            res%function_calls = res%function_calls + 1
            if (.not. finite_value(fv)) then
                res%root = p0
                call fail_result(res, root_status_value_error)
                return
            end if
            if (near_zero(fv)) then
                call succeed_result(res, p0, i - 1)
                return
            end if
            before = res%function_calls
            tmp = finite_difference_newton_step(p0, fv)
            res%function_calls = res%function_calls + 1
            if (.not. tmp%converged) then
                res%root = p0
                res%iterations = i
                call fail_result(res, tmp%status)
                return
            end if
            step = tmp%root
            if (res%function_calls == before) then
                res%root = p0
                call fail_result(res, root_status_value_error)
                return
            end if
            p = p0 - step
            res%iterations = i
            if (isclose(p, p0, abs_tol, rel_tol)) then
                call succeed_result(res, p, i)
                return
            end if
            if (.not. finite_value(p)) then
                res%root = p0
                call fail_result(res, root_status_value_error)
                return
            end if
            p0 = p
        end do

        res%root = p0
        call fail_result(res, root_status_convergence_error)

    contains

        function finite_difference_newton_step(x, fx) result(step_res)
            real(dp), intent(in) :: x
            real(dp), intent(in) :: fx
            type(root_result) :: step_res
            type(derivative_result) :: dres
            real(dp) :: h

            call init_result(step_res, 'newton')
            h = sqrt(epsilon(1.0_dp)) * max(1.0_dp, abs(x))
            dres = derivative(f, x, order=2, maxiter=1, initial_step=h)
            if (.not. finite_value(dres%df)) then
                call fail_result(step_res, root_status_value_error)
                return
            end if
            if (near_zero(dres%df)) then
                call fail_result(step_res, root_status_convergence_error)
                return
            end if
            step_res%root = fx / dres%df
            step_res%converged = .true.
            step_res%status = root_status_success
            step_res%flag = 'converged'
        end function finite_difference_newton_step

    end function numerical_newton_solver

    function secant_solver(f, x0, x1, tol, rtol, maxiter) result(res)
        procedure(scalar_root_function) :: f
        real(dp), intent(in) :: x0
        real(dp), intent(in), optional :: x1
        real(dp), intent(in), optional :: tol
        real(dp), intent(in), optional :: rtol
        integer, intent(in), optional :: maxiter
        type(root_result) :: res

        integer :: i
        integer :: mit
        real(dp) :: abs_tol
        real(dp) :: eps0
        real(dp) :: p
        real(dp) :: p0
        real(dp) :: p1
        real(dp) :: q0
        real(dp) :: q1
        real(dp) :: rel_tol
        real(dp) :: tmpp
        real(dp) :: tmpq

        call init_result(res, 'secant')
        if (.not. newton_options(tol, rtol, maxiter, abs_tol, rel_tol, mit)) return
        if (.not. finite_value(x0)) return
        p0 = x0
        if (present(x1)) then
            if (.not. finite_value(x1)) return
            if (near_equal(x1, x0)) return
            p1 = x1
        else
            eps0 = 1.0e-4_dp
            p1 = x0 * (1.0_dp + eps0)
            if (p1 >= 0.0_dp) then
                p1 = p1 + eps0
            else
                p1 = p1 - eps0
            end if
        end if

        q0 = f(p0)
        q1 = f(p1)
        res%function_calls = 2
        if (.not. finite_value(q0) .or. .not. finite_value(q1)) then
            res%root = merge(p0, p1, .not. finite_value(q0))
            call fail_result(res, root_status_value_error)
            return
        end if
        if (abs(q1) < abs(q0)) then
            tmpp = p0
            p0 = p1
            p1 = tmpp
            tmpq = q0
            q0 = q1
            q1 = tmpq
        end if

        do i = 1, mit
            res%iterations = i
            if (near_equal(q1, q0)) then
                p = 0.5_dp * (p1 + p0)
                res%root = p
                call fail_result(res, root_status_convergence_error)
                return
            end if
            if (abs(q1) > abs(q0)) then
                p = (-q0 / q1 * p1 + p0) / (1.0_dp - q0 / q1)
            else
                p = (-q1 / q0 * p0 + p1) / (1.0_dp - q1 / q0)
            end if
            if (isclose(p, p1, abs_tol, rel_tol)) then
                call succeed_result(res, p, i)
                return
            end if
            p0 = p1
            q0 = q1
            p1 = p
            q1 = f(p1)
            res%function_calls = res%function_calls + 1
            if (.not. finite_value(q1)) then
                res%root = p1
                call fail_result(res, root_status_value_error)
                return
            end if
        end do

        res%root = p1
        call fail_result(res, root_status_convergence_error)
    end function secant_solver

    logical function bracket_options(xtol, rtol, maxiter, tol_abs, tol_rel, mit) result(ok)
        real(dp), intent(in), optional :: xtol
        real(dp), intent(in), optional :: rtol
        integer, intent(in), optional :: maxiter
        real(dp), intent(out) :: tol_abs
        real(dp), intent(out) :: tol_rel
        integer, intent(out) :: mit

        tol_abs = default_xtol
        if (present(xtol)) tol_abs = xtol
        tol_rel = default_rtol
        if (present(rtol)) tol_rel = rtol
        mit = default_bracket_maxiter
        if (present(maxiter)) mit = maxiter
        ok = finite_value(tol_abs)
        if (.not. ok) return
        ok = finite_value(tol_rel)
        if (.not. ok) return
        ok = tol_abs > 0.0_dp
        if (.not. ok) return
        ok = tol_rel >= default_rtol
        if (.not. ok) return
        ok = mit >= 0
    end function bracket_options

    logical function newton_options(tol, rtol, maxiter, tol_abs, tol_rel, mit) result(ok)
        real(dp), intent(in), optional :: tol
        real(dp), intent(in), optional :: rtol
        integer, intent(in), optional :: maxiter
        real(dp), intent(out) :: tol_abs
        real(dp), intent(out) :: tol_rel
        integer, intent(out) :: mit

        tol_abs = default_newton_tol
        if (present(tol)) tol_abs = tol
        tol_rel = 0.0_dp
        if (present(rtol)) tol_rel = rtol
        mit = default_newton_maxiter
        if (present(maxiter)) mit = maxiter
        ok = finite_value(tol_abs)
        if (.not. ok) return
        ok = finite_value(tol_rel)
        if (.not. ok) return
        ok = tol_abs > 0.0_dp
        if (.not. ok) return
        ok = tol_rel >= 0.0_dp
        if (.not. ok) return
        ok = mit >= 1
    end function newton_options

    pure logical function valid_ab(a, b) result(ok)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        ok = finite_value(a)
        if (.not. ok) return
        ok = finite_value(b)
        if (.not. ok) return
        ok = .not. near_equal(a, b)
    end function valid_ab

    pure logical function near_zero(x) result(ok)
        real(dp), intent(in) :: x
        ok = abs(x) <= tiny(1.0_dp)
    end function near_zero

    pure logical function near_equal(a, b) result(ok)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp) :: scale
        scale = max(1.0_dp, abs(a), abs(b))
        ok = abs(a - b) <= tiny(1.0_dp) * scale
    end function near_equal

    pure logical function finite_value(x) result(ok)
        real(dp), intent(in) :: x
        ok = ieee_is_finite(x)
    end function finite_value

    pure logical function opposite_signs(a, b) result(ok)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        ok = .false.
        if (a < 0.0_dp) then
            if (b > 0.0_dp) ok = .true.
        else if (a > 0.0_dp) then
            if (b < 0.0_dp) ok = .true.
        end if
    end function opposite_signs

    pure logical function isclose(a, b, atol, rtol) result(ok)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp), intent(in) :: atol
        real(dp), intent(in) :: rtol
        real(dp) :: scale
        scale = max(abs(a), abs(b))
        ok = abs(a - b) <= atol + rtol * scale
    end function isclose

    pure logical function inside_open(x, a, b) result(ok)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        if (a < b) then
            ok = x > a .and. x < b
        else
            ok = x > b .and. x < a
        end if
    end function inside_open

    pure real(dp) function safe_secant(x0, x1, f0, f1) result(x2)
        real(dp), intent(in) :: x0
        real(dp), intent(in) :: x1
        real(dp), intent(in) :: f0
        real(dp), intent(in) :: f1
        if (near_equal(f0, f1)) then
            x2 = quiet_nan()
        else if (abs(f1) > abs(f0)) then
            x2 = (-f0 / f1 * x1 + x0) / (1.0_dp - f0 / f1)
        else
            x2 = (-f1 / f0 * x0 + x1) / (1.0_dp - f1 / f0)
        end if
    end function safe_secant

    subroutine update_bracket(a, b, fa, fb, c, fc, discarded_x, discarded_f)
        real(dp), intent(inout) :: a
        real(dp), intent(inout) :: b
        real(dp), intent(inout) :: fa
        real(dp), intent(inout) :: fb
        real(dp), intent(in) :: c
        real(dp), intent(in) :: fc
        real(dp), intent(out) :: discarded_x
        real(dp), intent(out) :: discarded_f

        if (same_sign(fa, fc)) then
            discarded_x = a
            discarded_f = fa
            a = c
            fa = fc
        else
            discarded_x = b
            discarded_f = fb
            b = c
            fb = fc
        end if
        if (a > b) call swap_pairs(a, b, fa, fb)
    end subroutine update_bracket

    pure logical function same_sign(a, b) result(ok)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        ok = (a > 0.0_dp .and. b > 0.0_dp) .or. (a < 0.0_dp .and. b < 0.0_dp)
    end function same_sign

    subroutine swap_pairs(a, b, fa, fb)
        real(dp), intent(inout) :: a
        real(dp), intent(inout) :: b
        real(dp), intent(inout) :: fa
        real(dp), intent(inout) :: fb
        real(dp) :: t
        t = a
        a = b
        b = t
        t = fa
        fa = fb
        fb = t
    end subroutine swap_pairs

    pure logical function four_values_separated(f1, f2, f3, f4) result(ok)
        real(dp), intent(in) :: f1
        real(dp), intent(in) :: f2
        real(dp), intent(in) :: f3
        real(dp), intent(in) :: f4
        real(dp) :: fs(4)
        integer :: i
        integer :: j
        real(dp) :: scale
        real(dp) :: threshold

        fs = [f1, f2, f3, f4]
        ok = .false.
        do i = 1, 4
            if (.not. finite_value(fs(i))) return
            if (near_zero(fs(i))) return
        end do
        threshold = 32.0_dp * epsilon(1.0_dp)
        do i = 1, 3
            do j = i + 1, 4
                scale = max(1.0_dp, abs(fs(i)), abs(fs(j)))
                if (abs(fs(i) / scale - fs(j) / scale) <= threshold / scale) return
            end do
        end do
        ok = .true.
    end function four_values_separated

    pure real(dp) function inverse_cubic_zero(a, b, c, d, fa, fb, fc, fd) result(x)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp), intent(in) :: c
        real(dp), intent(in) :: d
        real(dp), intent(in) :: fa
        real(dp), intent(in) :: fb
        real(dp), intent(in) :: fc
        real(dp), intent(in) :: fd
        real(dp) :: xv(4)
        real(dp) :: fv(4)
        real(dp) :: q(4, 4)
        real(dp) :: dd(4, 4)
        real(dp) :: alpha
        real(dp) :: diffik
        integer :: i
        integer :: j
        integer :: k
        integer :: n

        xv = [a, b, c, d]
        fv = [fa, fb, fc, fd]
        q = 0.0_dp
        dd = 0.0_dp
        q(:, 1) = xv
        dd(:, 1) = xv
        n = 4
        do k = 2, n
            do i = k, n
                j = i - k + 1
                alpha = dd(i, k - 1) - q(i - 1, k - 1)
                diffik = fv(j) - fv(i)
                q(i, k) = fv(i) / diffik * alpha
                dd(i, k) = fv(j) / diffik * alpha
            end do
        end do
        x = sum(q(n, 1:n))
    end function inverse_cubic_zero

    pure real(dp) function divided_difference_two(a, b, fa, fb) result(dd)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp), intent(in) :: fa
        real(dp), intent(in) :: fb
        dd = (fb - fa) / (b - a)
    end function divided_difference_two

    pure real(dp) function newton_quadratic(a, b, fa, fb, d, fd, k) result(r)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: b
        real(dp), intent(in) :: fa
        real(dp), intent(in) :: fb
        real(dp), intent(in) :: d
        real(dp), intent(in) :: fd
        integer, intent(in) :: k
        real(dp) :: aa
        real(dp) :: bb
        real(dp) :: den
        real(dp) :: p
        real(dp) :: r1
        integer :: i

        bb = (fb - fa) / (b - a)
        den = d - b
        if (near_zero(den)) then
            r = a + 0.5_dp * (b - a)
            return
        end if
        aa = ((fd - fb) / den - bb) / (d - a)
        if (near_zero(aa)) then
            if (near_zero(bb)) then
                r = a + 0.5_dp * (b - a)
            else
                r = a - fa / bb
            end if
            return
        end if
        if (same_sign(aa, fa)) then
            r = a
        else
            r = b
        end if
        do i = 1, k
            p = (aa * (r - b) + bb) * (r - a) + fa
            den = bb + aa * (2.0_dp * r - a - b)
            if (near_zero(den)) exit
            r1 = r - p / den
            if (.not. inside_open(r1, a, b)) then
                if (inside_open(r, a, b)) return
                r = a + 0.5_dp * (b - a)
                return
            end if
            r = r1
        end do
    end function newton_quadratic

    subroutine init_result(res, method)
        type(root_result), intent(out) :: res
        character(len=*), intent(in) :: method
        res%root = quiet_nan()
        res%iterations = 0
        res%function_calls = 0
        res%status = root_status_invalid_input
        res%converged = .false.
        res%flag = 'invalid input'
        res%method = method
    end subroutine init_result

    subroutine succeed_result(res, root, iterations)
        type(root_result), intent(inout) :: res
        real(dp), intent(in) :: root
        integer, intent(in) :: iterations
        res%root = root
        res%iterations = iterations
        res%status = root_status_success
        res%converged = .true.
        res%flag = 'converged'
    end subroutine succeed_result

    subroutine fail_result(res, status)
        type(root_result), intent(inout) :: res
        integer, intent(in) :: status
        res%status = status
        res%converged = .false.
        select case (status)
        case (root_status_sign_error)
            res%flag = 'sign error'
        case (root_status_convergence_error)
            res%flag = 'convergence error'
        case (root_status_value_error)
            res%flag = 'value error'
        case default
            res%flag = 'invalid input'
        end select
    end subroutine fail_result

    pure real(dp) function quiet_nan() result(x)
        x = ieee_value(0.0_dp, ieee_quiet_nan)
    end function quiet_nan

    pure function lowercase(s) result(out)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: out
        integer :: c
        integer :: i
        do i = 1, len(s)
            c = iachar(s(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) then
                out(i:i) = achar(c + iachar('a') - iachar('A'))
            else
                out(i:i) = s(i:i)
            end if
        end do
    end function lowercase

end module scifort_optimize_root_scalar
