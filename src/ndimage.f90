! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Separable Gaussian and neighborhood filters for arrays of ranks 1 through 15.
! Independent finite-kernel implementation; SciPy is a reference, not source.
module scifort_ndimage
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_ndimage_watershed, only : watershed_ift_nd
    use scifort_ndimage_distance, only : distance_transform_edt_nd, distance_transform_bf_nd
    use scifort_ndimage_chamfer, only : distance_transform_cdt_nd, distance_transform_weighted_cdt_nd
    use scifort_ndimage_measurements, only : sum_labels_nd, mean_labels_nd, &
        variance_labels_nd, center_of_mass_nd, standard_deviation_labels_nd, &
        count_labels_nd, histogram_labels_nd, extrema_labels_nd, &
        minimum_labels_nd, maximum_labels_nd, minimum_position_nd, maximum_position_nd, &
        bounding_boxes_nd, find_objects_nd, value_indices_nd, &
        labeled_comprehension_nd, ndimage_region_reduce
    use scifort_ndimage_morphology, only : binary_dilation_nd, binary_erosion_nd, &
        label_nd, binary_opening_nd, binary_closing_nd, binary_fill_holes_nd, &
        binary_propagation_nd, binary_hit_or_miss_nd, generate_binary_structure_nd, &
        iterate_structure_nd, morphology_status_success, morphology_status_invalid_input, &
        morphology_status_too_expensive
    use scifort_signal_spline_nd, only : spline_filter_nd, spline_filter1d, &
        map_coordinates_nd, map_coordinates_grid_nd, &
        spline_interpolation_plan, prepare_spline_plan_nd, &
        evaluate_spline_plan_nd, evaluate_spline_plan_grid_nd, &
        evaluate_spline_plan_with_gradient_nd, &
        evaluate_spline_plan_grid_with_gradient_nd, &
        evaluate_spline_plan_affine_with_gradient_nd, &
        evaluate_spline_plan_with_hessian_nd, evaluate_spline_plan_grid_with_hessian_nd, &
        evaluate_spline_plan_displacement_grid_with_gradient_nd, &
        evaluate_spline_plan_affine_with_hessian_nd, &
        evaluate_spline_plan_displacement_with_hessian_nd, &
        evaluate_spline_plan_displacement_grid_with_hessian_nd, &
        spline_axis_grid, evaluate_spline_plan_tensor_nd, &
        evaluate_spline_plan_affine_nd, evaluate_spline_plan_geometric_nd, &
        evaluate_spline_plan_geometric_batch_nd, &
        evaluate_spline_plan_displacement_nd, evaluate_spline_plan_displacement_grid_nd, &
        evaluate_spline_plan_displacement_with_gradient_nd, &
        ndimage_coordinate_batch_map, &
        geometric_transform_nd, ndimage_coordinate_map, &
        affine_transform_nd, shift_nd, zoom_nd, rotate_nd, &
        zoom_shape_nd, rotate_shape_nd
    implicit none
    private

    integer, parameter, public :: ndimage_status_success = 0
    integer, parameter, public :: ndimage_status_invalid_input = 1
    integer, parameter, public :: ndimage_status_too_expensive = 2
    integer, parameter :: max_kernel_radius=100000

    public :: spline_filter_nd, spline_filter1d, map_coordinates_nd, map_coordinates_grid_nd
    public :: spline_interpolation_plan, prepare_spline_plan_nd
    public :: evaluate_spline_plan_nd, evaluate_spline_plan_grid_nd
    public :: evaluate_spline_plan_with_gradient_nd
    public :: evaluate_spline_plan_grid_with_gradient_nd
    public :: evaluate_spline_plan_affine_with_gradient_nd
    public :: evaluate_spline_plan_with_hessian_nd, evaluate_spline_plan_grid_with_hessian_nd
    public :: evaluate_spline_plan_displacement_grid_with_gradient_nd
    public :: evaluate_spline_plan_affine_with_hessian_nd
    public :: evaluate_spline_plan_displacement_with_hessian_nd
    public :: evaluate_spline_plan_displacement_grid_with_hessian_nd
    public :: spline_axis_grid, evaluate_spline_plan_tensor_nd
    public :: evaluate_spline_plan_affine_nd, evaluate_spline_plan_geometric_nd
    public :: evaluate_spline_plan_geometric_batch_nd, ndimage_coordinate_batch_map
    public :: evaluate_spline_plan_displacement_nd, evaluate_spline_plan_displacement_grid_nd
    public :: evaluate_spline_plan_displacement_with_gradient_nd
    public :: geometric_transform_nd, ndimage_coordinate_map
    public :: affine_transform_nd, shift_nd, zoom_nd, rotate_nd
    public :: zoom_shape_nd, rotate_shape_nd
    public :: gaussian_filter_nd, gaussian_filter1d_nd
    public :: fourier_gaussian_nd, fourier_uniform_nd, fourier_shift_nd, fourier_ellipsoid_nd
    public :: gaussian_gradient_magnitude_nd, gaussian_laplace_nd
    public :: generic_gradient_magnitude_nd, generic_laplace_nd, ndimage_derivative_func
    public :: ndimage_complex_derivative_func
    public :: uniform_filter_nd, minimum_filter_nd, maximum_filter_nd
    public :: uniform_filter1d_nd, minimum_filter1d_nd, maximum_filter1d_nd
    public :: sobel_nd, prewitt_nd, laplace_nd
    public :: correlate_nd, convolve_nd, correlate1d_nd, convolve1d_nd
    public :: generic_filter_nd, generic_filter1d_nd, ndimage_line_func, ndimage_complex_line_func
    public :: vectorized_filter_nd, vectorized_filter_multi_nd
    public :: ndimage_batch_real_func, ndimage_batch_complex_func
    public :: ndimage_batch_real_multi_func, ndimage_batch_complex_multi_func
    public :: ndimage_batch_real_to_complex_func, ndimage_batch_complex_to_real_func
    public :: ndimage_batch_real_to_complex_multi_func, ndimage_batch_complex_to_real_multi_func
    public :: ndimage_neighborhood_func, ndimage_complex_neighborhood_func, median_filter_nd
    public :: rank_filter_nd, percentile_filter_nd
    public :: watershed_ift_nd
    public :: distance_transform_edt_nd, distance_transform_bf_nd, distance_transform_cdt_nd, distance_transform_weighted_cdt_nd
    public :: sum_labels_nd, mean_labels_nd, variance_labels_nd, center_of_mass_nd
    public :: standard_deviation_labels_nd, count_labels_nd, histogram_labels_nd
    public :: extrema_labels_nd, minimum_labels_nd, maximum_labels_nd
    public :: minimum_position_nd, maximum_position_nd
    public :: bounding_boxes_nd, find_objects_nd, value_indices_nd, &
        labeled_comprehension_nd, ndimage_region_reduce
    public :: grey_dilation_nd, grey_erosion_nd, grey_opening_nd, grey_closing_nd
    public :: morphological_gradient_nd, morphological_laplace_nd
    public :: white_tophat_nd, black_tophat_nd
    public :: binary_dilation_nd, binary_erosion_nd, label_nd
    public :: binary_opening_nd, binary_closing_nd, binary_fill_holes_nd
    public :: binary_propagation_nd, binary_hit_or_miss_nd
    public :: generate_binary_structure_nd, iterate_structure_nd
    public :: morphology_status_success, morphology_status_invalid_input
    public :: morphology_status_too_expensive

    abstract interface
        subroutine ndimage_derivative_func(x,axis,y,mode,cval,status)
            import dp
            real(dp), intent(in) :: x(..) !! caller's image, ranks 1 through 15
            integer, intent(in) :: axis !! one-based differentiation axis
            real(dp), intent(out) :: y(..) !! derivative output, same shape as x
            character(len=*), intent(in) :: mode !! extension for selected axis
            real(dp), intent(in) :: cval !! constant extension value
            integer, intent(out) :: status !! 0 for success
        end subroutine ndimage_derivative_func
        subroutine ndimage_complex_derivative_func(x,axis,y,mode,cval,status)
            import dp
            complex(dp), intent(in) :: x(..) !! complex input, ranks 1 through 15
            integer, intent(in) :: axis !! one-based differentiation axis
            complex(dp), intent(out) :: y(..) !! complex axis derivative, same shape
            character(len=*), intent(in) :: mode !! extension for selected axis
            complex(dp), intent(in) :: cval !! complex constant extension value
            integer, intent(out) :: status !! 0 for success
        end subroutine ndimage_complex_derivative_func
        function ndimage_neighborhood_func(values) result(value)
            import dp
            real(dp), intent(in) :: values(:) !! active neighbors in C-order (SciPy convention)
            real(dp) :: value !! caller-provided scalar reduction
        end function ndimage_neighborhood_func
        subroutine ndimage_line_func(input_line,output_line)
            import dp
            real(dp), intent(in) :: input_line(:) !! padded source line
            real(dp), intent(out) :: output_line(:) !! output along requested axis
        end subroutine ndimage_line_func
        subroutine ndimage_complex_line_func(input_line,output_line)
            import dp
            complex(dp), intent(in) :: input_line(:) !! padded complex source line
            complex(dp), intent(out) :: output_line(:) !! output along requested axis
        end subroutine ndimage_complex_line_func
        function ndimage_complex_neighborhood_func(values) result(value)
            import dp
            complex(dp), intent(in) :: values(:) !! active complex neighbors in C-order
            complex(dp) :: value !! complex-valued neighborhood reduction
        end function ndimage_complex_neighborhood_func
        subroutine ndimage_batch_real_func(windows,results)
            import dp
            real(dp), intent(in) :: windows(:,:) !! C-order neighbors by Fortran-order output batch
            real(dp), intent(out) :: results(:) !! one scalar result per batch column
        end subroutine ndimage_batch_real_func
        subroutine ndimage_batch_complex_func(windows,results)
            import dp
            complex(dp), intent(in) :: windows(:,:) !! complex neighbors by output batch
            complex(dp), intent(out) :: results(:) !! one scalar result per batch column
        end subroutine ndimage_batch_complex_func
        subroutine ndimage_batch_real_multi_func(windows,results)
            import dp
            real(dp), intent(in) :: windows(:,:) !! active neighbors x batch
            real(dp), intent(out) :: results(:,:) !! channels x batch
        end subroutine ndimage_batch_real_multi_func
        subroutine ndimage_batch_complex_multi_func(windows,results)
            import dp
            complex(dp), intent(in) :: windows(:,:) !! complex neighbors x batch
            complex(dp), intent(out) :: results(:,:) !! channels x batch
        end subroutine ndimage_batch_complex_multi_func
        subroutine ndimage_batch_real_to_complex_func(windows,results)
            import dp
            real(dp), intent(in) :: windows(:,:) !! real neighbors by batch
            complex(dp), intent(out) :: results(:) !! complex result per neighborhood
        end subroutine ndimage_batch_real_to_complex_func
        subroutine ndimage_batch_complex_to_real_func(windows,results)
            import dp
            complex(dp), intent(in) :: windows(:,:) !! complex neighbors by batch
            real(dp), intent(out) :: results(:) !! real result per neighborhood
        end subroutine ndimage_batch_complex_to_real_func
        subroutine ndimage_batch_real_to_complex_multi_func(windows,results)
            import dp
            real(dp), intent(in) :: windows(:,:) !! real neighbors by batch
            complex(dp), intent(out) :: results(:,:) !! complex channels by batch
        end subroutine ndimage_batch_real_to_complex_multi_func
        subroutine ndimage_batch_complex_to_real_multi_func(windows,results)
            import dp
            complex(dp), intent(in) :: windows(:,:) !! complex neighbors by batch
            real(dp), intent(out) :: results(:,:) !! real channels by batch
        end subroutine ndimage_batch_complex_to_real_multi_func
    end interface

    interface vectorized_filter_multi_nd
        module procedure vectorized_filter_multi_real_nd, vectorized_filter_multi_complex_nd
        module procedure vectorized_filter_multi_r2c_nd, vectorized_filter_multi_c2r_nd
    end interface

    interface vectorized_filter_nd
        module procedure vectorized_filter_real_nd, vectorized_filter_complex_nd
        module procedure vectorized_filter_r2c_nd, vectorized_filter_c2r_nd
    end interface

    interface generic_gradient_magnitude_nd
        module procedure generic_gradient_magnitude_real_nd, generic_gradient_magnitude_complex_nd
    end interface
    interface generic_laplace_nd
        module procedure generic_laplace_real_nd, generic_laplace_complex_nd
    end interface
    interface generic_filter_nd
        module procedure generic_filter_real_nd, generic_filter_complex_nd
    end interface
    interface generic_filter1d_nd
        module procedure generic_filter1d_real_nd, generic_filter1d_complex_nd
    end interface
    interface correlate1d_nd
        module procedure correlate1d_real_nd, correlate1d_complex_nd
    end interface
    interface convolve1d_nd
        module procedure convolve1d_real_nd, convolve1d_complex_nd
    end interface
    interface correlate_nd
        module procedure correlate_real_nd, correlate_complex_nd
    end interface
    interface convolve_nd
        module procedure convolve_real_nd, convolve_complex_nd
    end interface
    interface median_filter_nd
        module procedure median_filter_scalar, median_filter_vector
    end interface
    interface rank_filter_nd
        module procedure rank_filter_scalar, rank_filter_vector
    end interface
    interface percentile_filter_nd
        module procedure percentile_filter_scalar, percentile_filter_vector
    end interface

    interface fourier_ellipsoid_nd
        module procedure fourier_ellipsoid_vector, fourier_ellipsoid_scalar
    end interface
    interface fourier_gaussian_nd
        module procedure fourier_gaussian_vector, fourier_gaussian_scalar
    end interface
    interface fourier_uniform_nd
        module procedure fourier_uniform_vector, fourier_uniform_scalar
    end interface
    interface fourier_shift_nd
        module procedure fourier_shift_vector, fourier_shift_scalar
    end interface
    interface gaussian_filter_nd
        module procedure gaussian_real_vector, gaussian_real_scalar
        module procedure gaussian_complex_vector, gaussian_complex_scalar
    end interface
    interface gaussian_filter1d_nd
        module procedure gaussian1d_real, gaussian1d_complex
    end interface
    interface gaussian_gradient_magnitude_nd
        module procedure gradient_magnitude_vector, gradient_magnitude_scalar
    end interface
    interface gaussian_laplace_nd
        module procedure gaussian_laplace_vector, gaussian_laplace_scalar
    end interface
    interface uniform_filter1d_nd
        module procedure uniform1d_real, uniform1d_complex
    end interface
    interface uniform_filter_nd
        module procedure uniform_real_vector, uniform_real_scalar
        module procedure uniform_complex_vector, uniform_complex_scalar
    end interface
    interface minimum_filter_nd
        module procedure minimum_vector, minimum_scalar
    end interface
    interface maximum_filter_nd
        module procedure maximum_vector, maximum_scalar
    end interface

contains
    pure logical function mode_valid(mode) result(ok)
        character(len=*), intent(in) :: mode !! boundary extension string
        select case(mode)
        case('reflect','grid-mirror','mirror','nearest','wrap','grid-wrap', &
             'constant','grid-constant')
            ok=.true.
        case default
            ok=.false.
        end select
    end function mode_valid

    pure integer function extension_index(i,n,mode) result(idx)
        integer, intent(in) :: i !! zero-based unbounded position
        integer, intent(in) :: n !! axis extent, positive
        character(len=*), intent(in) :: mode !! extension type
        integer :: z
        idx=0
        select case(mode)
        case('constant','grid-constant')
            if (i>=0 .and. i<n) idx=i+1
        case('nearest')
            idx=max(0,min(n-1,i))+1
        case('wrap','grid-wrap')
            idx=modulo(i,n)+1
        case('reflect','grid-mirror')
            z=modulo(i,2*n)
            if (z>=n) z=2*n-1-z
            idx=z+1
        case('mirror')
            if (n==1) then
                idx=1
            else
                z=modulo(i,2*(n-1))
                if (z>=n) z=2*(n-1)-z
                idx=z+1
            end if
        end select
    end function extension_index

    subroutine resolve_axis_modes(nd,mode,modes,boundaries,ok)
        integer, intent(in) :: nd !! number of image dimensions
        character(len=*), optional, intent(in) :: mode !! common extension mode
        character(len=*), optional, intent(in) :: modes(:) !! optional per-axis modes
        character(len=32), allocatable, intent(out) :: boundaries(:) !! resolved modes
        logical, intent(out) :: ok !! all modes recognized
        integer :: d

        ok=.false.
        if (nd<1 .or. nd>15) return
        if (present(mode) .and. present(modes)) return
        if (present(modes)) then
            if (size(modes)/=nd) return
        end if
        allocate(boundaries(nd))
        boundaries='reflect'
        if (present(mode)) boundaries=trim(mode)
        if (present(modes)) then
            do d=1,nd
                if (len_trim(modes(d))>len(boundaries)) return
                boundaries(d)=trim(modes(d))
            end do
        end if
        do d=1,nd
            if (.not.mode_valid(trim(boundaries(d)))) return
        end do
        ok=.true.
    end subroutine resolve_axis_modes

    ! Axis-selective filters accept compact parameters in the same order as
    ! axes(:).  Outside the selected axes the size is one or sigma zero.
    subroutine selected_axes_mask(nd,axes,selected,ok)
        integer, intent(in) :: nd !! rank of source array
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        logical, intent(out) :: selected(nd) !! selected dimensions
        logical, intent(out) :: ok !! uniqueness and range validation
        integer :: k,d
        selected=.true.
        ok=.true.
        if (.not.present(axes)) return
        selected=.false.
        if (size(axes)>nd) then
            ok=.false.
            return
        end if
        do k=1,size(axes)
            d=axes(k)
            if (d<1 .or. d>nd) then
                ok=.false.
                return
            end if
            if (selected(d)) then
                ok=.false.
                return
            end if
            selected(d)=.true.
        end do
    end subroutine selected_axes_mask

    subroutine expand_selected_int(values,nd,axes,identity,expanded,ok)
        integer, intent(in) :: values(:) !! compact parameters in axes order
        integer, intent(in) :: nd !! array rank
        integer, optional, intent(in) :: axes(:) !! selected axes
        integer, intent(in) :: identity !! value on unselected axes
        integer, allocatable, intent(out) :: expanded(:) !! full-rank parameters
        logical, intent(out) :: ok !! valid axis selection and parameter count
        logical :: selected(nd)
        integer :: k
        call selected_axes_mask(nd,axes,selected,ok)
        if (.not.ok) return
        if (present(axes)) then
            if (size(values)/=size(axes)) then
                ok=.false.
                return
            end if
            allocate(expanded(nd),source=identity)
            do k=1,size(axes)
                expanded(axes(k))=values(k)
            end do
        else
            if (size(values)/=nd) then
                ok=.false.
                return
            end if
            expanded=values
        end if
    end subroutine expand_selected_int

    subroutine expand_selected_real(values,nd,axes,identity,expanded,ok)
        real(dp), intent(in) :: values(:) !! compact parameters in axes order
        integer, intent(in) :: nd !! array rank
        integer, optional, intent(in) :: axes(:) !! selected axes
        real(dp), intent(in) :: identity !! value on unselected axes
        real(dp), allocatable, intent(out) :: expanded(:) !! full-rank parameters
        logical, intent(out) :: ok !! valid selection
        logical :: selected(nd)
        integer :: k
        call selected_axes_mask(nd,axes,selected,ok)
        if (.not.ok) return
        if (present(axes)) then
            if (size(values)/=size(axes)) then
                ok=.false.
                return
            end if
            allocate(expanded(nd),source=identity)
            do k=1,size(axes)
                expanded(axes(k))=values(k)
            end do
        else
            if (size(values)/=nd) then
                ok=.false.
                return
            end if
            expanded=values
        end if
    end subroutine expand_selected_real

    subroutine selected_filter_modes(nd,axes,mode,modes,boundaries,ok)
        integer, intent(in) :: nd !! rank of image
        integer, optional, intent(in) :: axes(:) !! one-based axes
        character(len=*), optional, intent(in) :: mode !! common boundary
        character(len=*), optional, intent(in) :: modes(:) !! modes for selected axes
        character(len=32), allocatable, intent(out) :: boundaries(:) !! full-rank modes
        logical, intent(out) :: ok !! validation result
        logical :: selected(nd)
        integer :: d,k
        call selected_axes_mask(nd,axes,selected,ok)
        if (.not.ok) return
        if (.not.present(axes)) then
            call resolve_axis_modes(nd,mode,modes,boundaries,ok)
            return
        end if
        ok=.false.
        if (present(mode) .and. present(modes)) return
        allocate(boundaries(nd))
        boundaries='reflect'
        if (present(mode)) then
            if (len_trim(mode)>len(boundaries)) return
            boundaries=trim(mode)
        end if
        if (present(modes)) then
            if (size(modes)/=size(axes)) return
            do k=1,size(axes)
                if (len_trim(modes(k))>len(boundaries)) return
                boundaries(axes(k))=trim(modes(k))
            end do
        end if
        do d=1,nd
            if (.not.mode_valid(trim(boundaries(d)))) return
        end do
        ok=.true.
    end subroutine selected_filter_modes

    ! Monotone deque for a one-dimensional minimum/maximum over an extended
    ! line.  Each sample is inserted and removed at most once, so the work is
    ! O(n + width) rather than O(n * width).  Boundary samples are generated
    ! on demand; the window is not materialized as a padded image.
    subroutine extrema_axis_fast(values,dims,axis,width,mode,fill,maximum,origin)
        real(dp), intent(inout) :: values(:) !! packed values, modified along one axis
        integer, intent(in) :: dims(:) !! array extents
        integer, intent(in) :: axis !! one-based axis
        integer, intent(in) :: width !! positive window width
        character(len=*), intent(in) :: mode !! valid extension mode
        real(dp), intent(in) :: fill !! fill for constant modes
        logical, intent(in) :: maximum !! select maximum instead of minimum
        integer, intent(in) :: origin !! window displacement
        real(dp), allocatable :: result(:),queue_values(:)
        integer, allocatable :: queue_positions(:)
        integer(int64) :: stride,outer,block,o,j,pos
        integer :: n,p,t,front,back,left,center
        real(dp) :: v

        n=dims(axis)
        stride=product(int(dims(:axis-1),int64))
        outer=product(int(dims(axis+1:),int64))
        block=stride*int(n,int64)
        center=width/2+origin
        allocate(result(n),queue_values(n+width-1),queue_positions(n+width-1))
        do o=0_int64,outer-1_int64
            do j=0_int64,stride-1_int64
                front=1
                back=0
                do p=0,n+width-2
                    t=extension_index(p-center,n,mode)
                    v=fill
                    if (t/=0) then
                        pos=o*block+j+int(t-1,int64)*stride+1_int64
                        v=values(pos)
                    end if
                    if (maximum) then
                        do while(back>=front)
                            if (queue_values(back)>v) exit
                            back=back-1
                        end do
                    else
                        do while(back>=front)
                            if (queue_values(back)<v) exit
                            back=back-1
                        end do
                    end if
                    back=back+1
                    queue_values(back)=v
                    queue_positions(back)=p
                    left=p-width+1
                    if (left>=0) then
                        do while(queue_positions(front)<left)
                            front=front+1
                        end do
                        result(left+1)=queue_values(front)
                    end if
                end do
                do p=0,n-1
                    pos=o*block+j+int(p,int64)*stride+1_int64
                    values(pos)=result(p+1)
                end do
            end do
        end do
    end subroutine extrema_axis_fast

    ! Linear-time box filter for one axis.  Each new output window removes
    ! its leftmost value and adds its next rightmost value.  Like the deque
    ! extrema implementation, extension samples are generated on demand.
    ! The optional origin is used by uniform_filter1d_nd only.
    subroutine uniform_axis_running(values,dims,axis,width,mode,fill,origin)
        complex(dp), intent(inout) :: values(:) !! packed array, modified along one axis
        integer, intent(in) :: dims(:) !! full image dimensions
        integer, intent(in) :: axis !! one-based axis to filter
        integer, intent(in) :: width !! positive averaging window length
        character(len=*), intent(in) :: mode !! extension mode
        complex(dp), intent(in) :: fill !! constant boundary value
        integer, intent(in) :: origin !! SciPy center displacement
        complex(dp), allocatable :: result(:)
        integer(int64) :: stride,outer,block,o,j,pos
        integer :: n,i,k,t,center
        complex(dp) :: running,left_value,right_value

        if (width==1) return
        n=dims(axis)
        stride=product(int(dims(:axis-1),int64))
        outer=product(int(dims(axis+1:),int64))
        block=stride*int(n,int64)
        center=width/2+origin
        allocate(result(n))
        do o=0_int64,outer-1_int64
            do j=0_int64,stride-1_int64
                running=(0.0_dp,0.0_dp)
                do k=0,width-1
                    t=extension_index(k-center,n,mode)
                    right_value=fill
                    if (t/=0) then
                        pos=o*block+j+int(t-1,int64)*stride+1_int64
                        right_value=values(pos)
                    end if
                    running=running+right_value
                end do
                result(1)=running/real(width,dp)
                do i=1,n-1
                    t=extension_index(i-1-center,n,mode)
                    left_value=fill
                    if (t/=0) then
                        pos=o*block+j+int(t-1,int64)*stride+1_int64
                        left_value=values(pos)
                    end if
                    t=extension_index(i+width-1-center,n,mode)
                    right_value=fill
                    if (t/=0) then
                        pos=o*block+j+int(t-1,int64)*stride+1_int64
                        right_value=values(pos)
                    end if
                    running=running+(right_value-left_value)
                    result(i+1)=running/real(width,dp)
                end do
                do i=0,n-1
                    pos=o*block+j+int(i,int64)*stride+1_int64
                    values(pos)=result(i+1)
                end do
            end do
        end do
    end subroutine uniform_axis_running

    subroutine axis_convolution(values,dims,axis,offsets,weights,mode,fill)
        complex(dp), intent(inout) :: values(:) !! packed input and output
        integer, intent(in) :: dims(:) !! packed array extents
        integer, intent(in) :: axis !! one-based filtering axis
        integer, intent(in) :: offsets(:) !! offsets relative to output point
        real(dp), intent(in) :: weights(:) !! correlation weights
        character(len=*), intent(in) :: mode !! boundary mode
        complex(dp), intent(in) :: fill !! constant extension value
        complex(dp), allocatable :: result(:)
        integer(int64) :: stride,outer,block,o,j,pos
        integer :: n,i,k,t
        complex(dp) :: acc

        n=dims(axis)
        stride=product(int(dims(:axis-1),int64))
        outer=product(int(dims(axis+1:),int64))
        block=stride*int(n,int64)
        allocate(result(n))
        do o=0_int64,outer-1_int64
            do j=0_int64,stride-1_int64
                do i=0,n-1
                    acc=(0.0_dp,0.0_dp)
                    do k=1,size(offsets)
                        t=extension_index(i+offsets(k),n,mode)
                        if (t==0) then
                            acc=acc+weights(k)*fill
                        else
                            pos=o*block+j+int(t-1,int64)*stride+1_int64
                            acc=acc+weights(k)*values(pos)
                        end if
                    end do
                    result(i+1)=acc
                end do
                do i=0,n-1
                    pos=o*block+j+int(i,int64)*stride+1_int64
                    values(pos)=result(i+1)
                end do
            end do
        end do
    end subroutine axis_convolution

    subroutine gaussian_kernel(sigma,order,radius,offsets,weights)
        real(dp), intent(in) :: sigma !! Gaussian standard deviation, positive
        integer, intent(in) :: order !! derivative order 0 through 3
        integer, intent(in) :: radius !! nonnegative kernel radius
        integer, allocatable, intent(out) :: offsets(:) !! symmetric lags
        real(dp), allocatable, intent(out) :: weights(:) !! correlation kernel
        real(dp) :: x,t,scale,scale2
        integer :: k

        allocate(offsets(2*radius+1),weights(2*radius+1))
        scale=1.0_dp/(sigma*sigma)
        scale2=scale*scale
        do k=-radius,radius
            x=real(k,dp)
            offsets(k+radius+1)=k
            weights(k+radius+1)=exp(-0.5_dp*(x/sigma)**2)
        end do
        weights=weights/sum(weights)
        do k=1,size(weights)
            x=real(offsets(k),dp)
            select case(order)
            case(0)
                t=1.0_dp
            case(1)
                t=x*scale
            case(2)
                t=x*x*scale2-scale
            case(3)
                t=x*x*x*scale2*scale-3.0_dp*x*scale2
            end select
            weights(k)=weights(k)*t
        end do
    end subroutine gaussian_kernel

    subroutine gaussian_process(values,dims,sigma,orders,mode,fill,truncate,radii,status,modes,axes)
        complex(dp), intent(inout) :: values(:) !! packed image, updated in place
        integer, intent(in) :: dims(:) !! image extents
        real(dp), intent(in) :: sigma(:) !! per-axis standard deviations
        integer, intent(in) :: orders(:) !! derivative orders 0..3 per axis
        character(len=*), intent(in) :: mode !! common boundary mode
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary modes
        integer, optional, intent(in) :: axes(:) !! original filtering order, one-based
        complex(dp), intent(in) :: fill !! constant extension value
        real(dp), intent(in) :: truncate !! kernel cutoff in sigmas
        integer, intent(in) :: radii(:) !! negative: calculate from truncate
        integer, intent(out) :: status !! zero on success
        integer, allocatable :: offsets(:)
        real(dp), allocatable :: weights(:)
        real(dp) :: bound
        integer :: d,r,k
        character(len=32), allocatable :: boundaries(:)
        logical :: valid_modes

        status=ndimage_status_invalid_input
        if (size(dims)<1 .or. size(dims)>15) return
        if (size(sigma)/=size(dims) .or. size(orders)/=size(dims) .or. &
            size(radii)/=size(dims)) return
        if (any(dims<1) .or. any(sigma<0.0_dp)) return
        if (.not.all(ieee_is_finite(sigma))) return
        if (any(orders<0) .or. any(orders>3) .or. any(radii< -1)) return
        call resolve_axis_modes(size(dims),mode=mode,boundaries=boundaries,ok=valid_modes)
        if (present(modes)) call resolve_axis_modes(size(dims),modes=modes,boundaries=boundaries,ok=valid_modes)
        if (.not.valid_modes) return
        if (.not.ieee_is_finite(truncate) .or. truncate<0.0_dp) return
        if (.not.ieee_is_finite(real(fill,dp)) .or. .not.ieee_is_finite(aimag(fill))) return
        do k=1,size(dims)
            d=k
            if (present(axes)) then
                if (k>size(axes)) exit
                d=axes(k)
            end if
            if (sigma(d)==0.0_dp) then
                if (orders(d)>0) return
                cycle
            end if
            bound=truncate*sigma(d)
            if (.not.ieee_is_finite(bound)) return
            if (radii(d)>=0) then
                r=radii(d)
            else
                if (bound>real(max_kernel_radius,dp)) then
                    status=ndimage_status_too_expensive
                    return
                end if
                r=int(bound+0.5_dp)
            end if
            if (r>max_kernel_radius) then
                status=ndimage_status_too_expensive
                return
            end if
            call gaussian_kernel(sigma(d),orders(d),r,offsets,weights)
            call axis_convolution(values,dims,d,offsets,weights,trim(boundaries(d)),fill)
            deallocate(offsets,weights)
        end do
        status=ndimage_status_success
    end subroutine gaussian_process

    subroutine uniform_process(values,dims,widths,mode,fill,status,modes,axes)
        complex(dp), intent(inout) :: values(:) !! packed image, updated in place
        integer, intent(in) :: dims(:) !! extents
        integer, intent(in) :: widths(:) !! rectangular window widths
        character(len=*), intent(in) :: mode !! boundary extension
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary modes
        integer, optional, intent(in) :: axes(:) !! original filtering order, one-based
        complex(dp), intent(in) :: fill !! constant boundary value
        integer, intent(out) :: status !! zero on success
        integer :: d,n,step
        character(len=32), allocatable :: boundaries(:)
        logical :: valid_modes

        status=ndimage_status_invalid_input
        if (size(widths)/=size(dims)) return
        if (any(widths<1)) return
        call resolve_axis_modes(size(dims),mode=mode,boundaries=boundaries,ok=valid_modes)
        if (present(modes)) call resolve_axis_modes(size(dims),modes=modes,boundaries=boundaries,ok=valid_modes)
        if (.not.valid_modes) return
        if (.not.ieee_is_finite(real(fill,dp)) .or. .not.ieee_is_finite(aimag(fill))) return
        do step=1,size(dims)
            d=step
            if (present(axes)) then
                if (step>size(axes)) exit
                d=axes(step)
            end if
            n=widths(d)
            if (n>2*max_kernel_radius+1 .or. &
                int(n,int64)+int(dims(d),int64)-1_int64>int(huge(1),int64)) then
                status=ndimage_status_too_expensive
                return
            end if
            if (n>1) call uniform_axis_running(values,dims,d,n,trim(boundaries(d)),fill,0)
        end do
        status=ndimage_status_success
    end subroutine uniform_process

    subroutine extrema_process(values,dims,widths,mode,fill,maximum,status,modes,axes)
        real(dp), intent(inout) :: values(:) !! packed image, modified in place
        integer, intent(in) :: dims(:) !! original image extents
        integer, intent(in) :: widths(:) !! width along each axis
        character(len=*), intent(in) :: mode !! common boundary mode
        real(dp), intent(in) :: fill !! constant boundary value
        logical, intent(in) :: maximum !! calculate maximum instead of minimum
        integer, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: modes(:) !! per-axis modes
        integer, optional, intent(in) :: axes(:) !! original filtering order, one-based
        integer :: d,step
        character(len=:), allocatable :: boundary

        status=ndimage_status_invalid_input
        if (size(widths)/=size(dims)) return
        if (any(dims<1) .or. any(widths<1)) return
        if (.not.mode_valid(mode) .or. .not.ieee_is_finite(fill)) return
        if (present(modes)) then
            if (size(modes)/=size(dims)) return
            do d=1,size(dims)
                if (.not.mode_valid(trim(modes(d)))) return
            end do
        end if
        if (any(widths>2*max_kernel_radius+1)) then
            status=ndimage_status_too_expensive
            return
        end if
        do d=1,size(dims)
            ! Guard the index range used by the sliding deque.
            if (int(dims(d),int64)+int(widths(d),int64)-1_int64>int(huge(1),int64)) then
                status=ndimage_status_too_expensive
                return
            end if
        end do
        do step=1,size(dims)
            d=step
            if (present(axes)) then
                if (step>size(axes)) exit
                d=axes(step)
            end if
            if (widths(d)==1) cycle
            boundary=mode
            if (present(modes)) boundary=trim(modes(d))
            call extrema_axis_fast(values,dims,d,widths(d),boundary,fill,maximum,0)
        end do
        status=ndimage_status_success
    end subroutine extrema_process

    subroutine pack_real(x,buf)
        real(dp), intent(in) :: x(..) !! input array rank one through fifteen
        real(dp), allocatable, intent(out) :: buf(:) !! column-major linearized data
        select rank(x)
        rank(1)
            buf=reshape(x,[size(x)])
        rank(2)
            buf=reshape(x,[size(x)])
        rank(3)
            buf=reshape(x,[size(x)])
        rank(4)
            buf=reshape(x,[size(x)])
        rank(5)
            buf=reshape(x,[size(x)])
        rank(6)
            buf=reshape(x,[size(x)])
        rank(7)
            buf=reshape(x,[size(x)])
        rank(8)
            buf=reshape(x,[size(x)])
        rank(9)
            buf=reshape(x,[size(x)])
        rank(10)
            buf=reshape(x,[size(x)])
        rank(11)
            buf=reshape(x,[size(x)])
        rank(12)
            buf=reshape(x,[size(x)])
        rank(13)
            buf=reshape(x,[size(x)])
        rank(14)
            buf=reshape(x,[size(x)])
        rank(15)
            buf=reshape(x,[size(x)])
        rank default
            allocate(buf(0))
        end select
    end subroutine pack_real

    subroutine scatter_real(buf,y)
        real(dp), intent(in) :: buf(:) !! column-major linearized data
        real(dp), intent(out) :: y(..) !! destination, same shape as original
        select rank(y)
        rank(1)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(2)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(3)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(4)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(5)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(6)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(7)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(8)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(9)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(10)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(11)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(12)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(13)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(14)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(15)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank default
        end select
    end subroutine scatter_real

    subroutine pack_complex(x,buf)
        complex(dp), intent(in) :: x(..) !! input array rank one through fifteen
        complex(dp), allocatable, intent(out) :: buf(:) !! column-major linearized data
        select rank(x)
        rank(1)
            buf=reshape(x,[size(x)])
        rank(2)
            buf=reshape(x,[size(x)])
        rank(3)
            buf=reshape(x,[size(x)])
        rank(4)
            buf=reshape(x,[size(x)])
        rank(5)
            buf=reshape(x,[size(x)])
        rank(6)
            buf=reshape(x,[size(x)])
        rank(7)
            buf=reshape(x,[size(x)])
        rank(8)
            buf=reshape(x,[size(x)])
        rank(9)
            buf=reshape(x,[size(x)])
        rank(10)
            buf=reshape(x,[size(x)])
        rank(11)
            buf=reshape(x,[size(x)])
        rank(12)
            buf=reshape(x,[size(x)])
        rank(13)
            buf=reshape(x,[size(x)])
        rank(14)
            buf=reshape(x,[size(x)])
        rank(15)
            buf=reshape(x,[size(x)])
        rank default
            allocate(buf(0))
        end select
    end subroutine pack_complex

    subroutine scatter_complex(buf,y)
        complex(dp), intent(in) :: buf(:) !! column-major linearized data
        complex(dp), intent(out) :: y(..) !! destination, same shape as original
        select rank(y)
        rank(1)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(2)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(3)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(4)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(5)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(6)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(7)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(8)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(9)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(10)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(11)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(12)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(13)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(14)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank(15)
            if (size(buf)==size(y)) y=reshape(buf,shape(y))
        rank default
        end select
    end subroutine scatter_complex

    subroutine gaussian_real_work(x,sigma,y,orders,status,mode,cval,truncate,radius,modes,axes)
        real(dp), intent(in) :: x(..) !! source image, rank 1..15
        real(dp), intent(in) :: sigma(:) !! per-axis Gaussian widths
        real(dp), intent(out) :: y(..) !! caller-allocated filtered image
        integer, intent(in) :: orders(:) !! per-axis derivative orders
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! default reflect
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary modes
        real(dp), optional, intent(in) :: cval !! default zero
        real(dp), optional, intent(in) :: truncate !! default four sigmas
        integer, optional, intent(in) :: radius(:) !! explicit radii per axis
        integer, optional, intent(in) :: axes(:) !! selected filtering axes
        real(dp), allocatable :: sigmas(:),data(:)
        complex(dp), allocatable :: work(:)
        integer, allocatable :: rad(:),derivatives(:)
        logical :: ok
        character(len=32), allocatable :: boundaries(:)
        integer :: code
        character(len=:), allocatable :: boundary
        real(dp) :: cutoff
        complex(dp) :: fill

        code=ndimage_status_invalid_input
        if (present(status)) status=code
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        call pack_real(x,data)
        if (.not.all(ieee_is_finite(data))) return
        work=cmplx(data,0.0_dp,dp)
        call expand_selected_real(sigma,rank(x),axes,0.0_dp,sigmas,ok)
        if (.not.ok) return
        call expand_selected_int(orders,rank(x),axes,0,derivatives,ok)
        if (.not.ok) return
        rad=spread(-1,1,rank(x))
        if (present(radius)) then
            call expand_selected_int(radius,rank(x),axes,-1,rad,ok)
            if (.not.ok) return
        end if
        call selected_filter_modes(rank(x),axes,mode,modes,boundaries,ok)
        if (.not.ok) return
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cmplx(cval,0.0_dp,dp)
        cutoff=4.0_dp
        if (present(truncate)) cutoff=truncate
        call gaussian_process(work,shape(x),sigmas,derivatives,'reflect',fill,cutoff,rad,code,boundaries,axes)
        if (code==ndimage_status_success) then
            data=real(work,dp)
            call scatter_real(data,y)
        end if
        if (present(status)) status=code
    end subroutine gaussian_real_work

    subroutine gaussian_complex_work(x,sigma,y,orders,status,mode,cval,truncate,radius,modes,axes)
        complex(dp), intent(in) :: x(..) !! complex input, rank 1..15
        real(dp), intent(in) :: sigma(:) !! widths in each dimension
        complex(dp), intent(out) :: y(..) !! filtered complex image
        integer, intent(in) :: orders(:) !! derivative orders 0..3
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! default reflect
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary modes
        complex(dp), optional, intent(in) :: cval !! default zero
        real(dp), optional, intent(in) :: truncate !! default four sigmas
        integer, optional, intent(in) :: radius(:) !! per-axis radius override
        integer, optional, intent(in) :: axes(:) !! selected filtering axes
        real(dp), allocatable :: sigmas(:)
        complex(dp), allocatable :: work(:)
        integer, allocatable :: rad(:),derivatives(:)
        logical :: ok
        character(len=32), allocatable :: boundaries(:)
        integer :: code
        character(len=:), allocatable :: boundary
        real(dp) :: cutoff
        complex(dp) :: fill

        code=ndimage_status_invalid_input
        if (present(status)) status=code
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        call pack_complex(x,work)
        if (.not.all(ieee_is_finite(real(work,dp))) .or. &
            .not.all(ieee_is_finite(aimag(work)))) return
        call expand_selected_real(sigma,rank(x),axes,0.0_dp,sigmas,ok)
        if (.not.ok) return
        call expand_selected_int(orders,rank(x),axes,0,derivatives,ok)
        if (.not.ok) return
        rad=spread(-1,1,rank(x))
        if (present(radius)) then
            call expand_selected_int(radius,rank(x),axes,-1,rad,ok)
            if (.not.ok) return
        end if
        call selected_filter_modes(rank(x),axes,mode,modes,boundaries,ok)
        if (.not.ok) return
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cval
        cutoff=4.0_dp
        if (present(truncate)) cutoff=truncate
        call gaussian_process(work,shape(x),sigmas,derivatives,'reflect',fill,cutoff,rad,code,boundaries,axes)
        if (code==ndimage_status_success) call scatter_complex(work,y)
        if (present(status)) status=code
    end subroutine gaussian_complex_work

    subroutine gaussian_real_vector(x,sigma,y,order,status,mode,cval,truncate,radius,modes,axes)
        real(dp), intent(in) :: x(..) !! source array
        real(dp), intent(in) :: sigma(:) !! per-axis standard deviations
        real(dp), intent(out) :: y(..) !! filtered array, same shape
        integer, optional, intent(in) :: order(:) !! derivative orders; default zero
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary extension
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary modes
        real(dp), optional, intent(in) :: cval !! constant outside value
        real(dp), optional, intent(in) :: truncate !! cutoff in sigmas
        integer, optional, intent(in) :: radius(:) !! per-axis kernel radii
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        integer, allocatable :: orders(:)
        if (present(axes)) then
            orders=spread(0,1,size(axes))
        else
            orders=spread(0,1,rank(x))
        end if
        if (present(order)) orders=order
        call gaussian_real_work(x,sigma,y,orders,status,mode,cval,truncate,radius,modes,axes)
    end subroutine gaussian_real_vector

    subroutine gaussian_complex_vector(x,sigma,y,order,status,mode,cval,truncate,radius,modes,axes)
        complex(dp), intent(in) :: x(..) !! source image
        real(dp), intent(in) :: sigma(:) !! widths per axis
        complex(dp), intent(out) :: y(..) !! filtered image
        integer, optional, intent(in) :: order(:) !! derivative orders
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary extension
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary modes
        complex(dp), optional, intent(in) :: cval !! outside value
        real(dp), optional, intent(in) :: truncate !! Gaussian cutoff
        integer, optional, intent(in) :: radius(:) !! kernel radii per axis
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        integer, allocatable :: orders(:)
        if (present(axes)) then
            orders=spread(0,1,size(axes))
        else
            orders=spread(0,1,rank(x))
        end if
        if (present(order)) orders=order
        call gaussian_complex_work(x,sigma,y,orders,status,mode,cval,truncate,radius,modes,axes)
    end subroutine gaussian_complex_vector

    subroutine gaussian_real_scalar(x,sigma,y,order,status,mode,cval,truncate,radius,modes,axes)
        real(dp), intent(in) :: x(..) !! source array
        real(dp), intent(in) :: sigma !! shared Gaussian width
        real(dp), intent(out) :: y(..) !! filtered image
        integer, optional, intent(in) :: order !! common derivative order
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary extension
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary modes
        real(dp), optional, intent(in) :: cval !! constant boundary value
        real(dp), optional, intent(in) :: truncate !! cutoff in sigmas
        integer, optional, intent(in) :: radius !! common kernel radius
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        integer :: derivative,r
        integer, allocatable :: rad(:)
        derivative=0
        if (present(order)) derivative=order
        r=-1
        if (present(radius)) r=radius
        if (present(axes)) then
            rad=spread(r,1,size(axes))
            call gaussian_real_work(x,spread(sigma,1,size(axes)),y, &
                spread(derivative,1,size(axes)),status,mode,cval,truncate,rad,modes,axes)
        else
            rad=spread(r,1,rank(x))
            call gaussian_real_work(x,spread(sigma,1,rank(x)),y, &
                spread(derivative,1,rank(x)),status,mode,cval,truncate,rad,modes)
        end if
    end subroutine gaussian_real_scalar

    subroutine gaussian_complex_scalar(x,sigma,y,order,status,mode,cval,truncate,radius,modes,axes)
        complex(dp), intent(in) :: x(..) !! source complex array
        real(dp), intent(in) :: sigma !! shared Gaussian width
        complex(dp), intent(out) :: y(..) !! filtered image
        integer, optional, intent(in) :: order !! common derivative order
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary extension
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary modes
        complex(dp), optional, intent(in) :: cval !! constant boundary value
        real(dp), optional, intent(in) :: truncate !! cutoff in sigmas
        integer, optional, intent(in) :: radius !! common kernel radius
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        integer :: derivative,r
        integer, allocatable :: rad(:)
        derivative=0
        if (present(order)) derivative=order
        r=-1
        if (present(radius)) r=radius
        if (present(axes)) then
            rad=spread(r,1,size(axes))
            call gaussian_complex_work(x,spread(sigma,1,size(axes)),y, &
                spread(derivative,1,size(axes)),status,mode,cval,truncate,rad,modes,axes)
        else
            rad=spread(r,1,rank(x))
            call gaussian_complex_work(x,spread(sigma,1,rank(x)),y, &
                spread(derivative,1,rank(x)),status,mode,cval,truncate,rad,modes)
        end if
    end subroutine gaussian_complex_scalar

    subroutine gaussian1d_real(x,sigma,axis,y,order,status,mode,cval,truncate,radius)
        real(dp), intent(in) :: x(..) !! original data
        real(dp), intent(in) :: sigma !! Gaussian width on selected axis
        integer, intent(in) :: axis !! one-based axis
        real(dp), intent(out) :: y(..) !! filtered image
        integer, optional, intent(in) :: order !! derivative order 0..3
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary extension
        real(dp), optional, intent(in) :: cval !! constant outside value
        real(dp), optional, intent(in) :: truncate !! Gaussian cutoff
        integer, optional, intent(in) :: radius !! selected-axis kernel radius
        integer, allocatable :: rad(:),orders(:)
        real(dp), allocatable :: sigmas(:)
        integer :: derivative
        if (present(status)) status=ndimage_status_invalid_input
        if (axis<1 .or. axis>rank(x)) return
        sigmas=spread(0.0_dp,1,rank(x))
        orders=spread(0,1,rank(x))
        rad=spread(-1,1,rank(x))
        sigmas(axis)=sigma
        derivative=0
        if (present(order)) derivative=order
        orders(axis)=derivative
        if (present(radius)) rad(axis)=radius
        call gaussian_real_work(x,sigmas,y,orders,status,mode,cval,truncate,rad)
    end subroutine gaussian1d_real

    subroutine gaussian1d_complex(x,sigma,axis,y,order,status,mode,cval,truncate,radius)
        complex(dp), intent(in) :: x(..) !! complex original data
        real(dp), intent(in) :: sigma !! Gaussian width on selected axis
        integer, intent(in) :: axis !! one-based axis
        complex(dp), intent(out) :: y(..) !! filtered image
        integer, optional, intent(in) :: order !! derivative order 0..3
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary extension
        complex(dp), optional, intent(in) :: cval !! constant outside value
        real(dp), optional, intent(in) :: truncate !! Gaussian cutoff
        integer, optional, intent(in) :: radius !! selected-axis kernel radius
        integer, allocatable :: rad(:),orders(:)
        real(dp), allocatable :: sigmas(:)
        integer :: derivative
        if (present(status)) status=ndimage_status_invalid_input
        if (axis<1 .or. axis>rank(x)) return
        sigmas=spread(0.0_dp,1,rank(x))
        orders=spread(0,1,rank(x))
        rad=spread(-1,1,rank(x))
        sigmas(axis)=sigma
        derivative=0
        if (present(order)) derivative=order
        orders(axis)=derivative
        if (present(radius)) rad(axis)=radius
        call gaussian_complex_work(x,sigmas,y,orders,status,mode,cval,truncate,rad)
    end subroutine gaussian1d_complex

    subroutine uniform_real_work(x,widths,y,status,mode,cval,modes,axes)
        real(dp), intent(in) :: x(..) !! real source array
        integer, intent(in) :: widths(:) !! rectangular neighborhood dimensions
        real(dp), intent(out) :: y(..) !! filtered image
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary extension
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary modes
        real(dp), optional, intent(in) :: cval !! constant fill value
        integer, optional, intent(in) :: axes(:) !! axes to filter
        real(dp), allocatable :: data(:)
        complex(dp), allocatable :: work(:)
        complex(dp) :: fill
        character(len=:), allocatable :: boundary
        integer, allocatable :: full_widths(:)
        character(len=32), allocatable :: boundaries(:)
        logical :: ok
        integer :: code
        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(x)/=rank(y)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        call pack_real(x,data)
        if (.not.all(ieee_is_finite(data))) return
        work=cmplx(data,0.0_dp,dp)
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cmplx(cval,0.0_dp,dp)
        call expand_selected_int(widths,rank(x),axes,1,full_widths,ok)
        if (.not.ok) return
        call selected_filter_modes(rank(x),axes,mode,modes,boundaries,ok)
        if (.not.ok) return
        call uniform_process(work,shape(x),full_widths,'reflect',fill,code,boundaries,axes)
        if (code==ndimage_status_success) then
            data=real(work,dp)
            call scatter_real(data,y)
        end if
        if (present(status)) status=code
    end subroutine uniform_real_work

    subroutine uniform_complex_work(x,widths,y,status,mode,cval,modes,axes)
        complex(dp), intent(in) :: x(..) !! complex source image
        integer, intent(in) :: widths(:) !! dimensions of rectangular filter
        complex(dp), intent(out) :: y(..) !! filtered image
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary extension
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary modes
        complex(dp), optional, intent(in) :: cval !! constant fill value
        integer, optional, intent(in) :: axes(:) !! axes to filter
        complex(dp), allocatable :: data(:)
        complex(dp) :: fill
        character(len=:), allocatable :: boundary
        integer, allocatable :: full_widths(:)
        character(len=32), allocatable :: boundaries(:)
        logical :: ok
        integer :: code
        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(x)/=rank(y)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        call pack_complex(x,data)
        if (.not.all(ieee_is_finite(real(data,dp))) .or. &
            .not.all(ieee_is_finite(aimag(data)))) return
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cval
        call expand_selected_int(widths,rank(x),axes,1,full_widths,ok)
        if (.not.ok) return
        call selected_filter_modes(rank(x),axes,mode,modes,boundaries,ok)
        if (.not.ok) return
        call uniform_process(data,shape(x),full_widths,'reflect',fill,code,boundaries,axes)
        if (code==ndimage_status_success) call scatter_complex(data,y)
        if (present(status)) status=code
    end subroutine uniform_complex_work

    subroutine uniform_real_vector(x,sizes,y,status,mode,cval,modes,axes)
        real(dp), intent(in) :: x(..) !! source array
        integer, intent(in) :: sizes(:) !! per-axis filter sizes
        real(dp), intent(out) :: y(..) !! filtered array
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! default reflect
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary modes
        real(dp), optional, intent(in) :: cval !! constant outside value
        integer, optional, intent(in) :: axes(:) !! one-based axes to filter
        call uniform_real_work(x,sizes,y,status,mode,cval,modes,axes)
    end subroutine uniform_real_vector

    subroutine uniform_real_scalar(x,sizes,y,status,mode,cval,modes,axes)
        real(dp), intent(in) :: x(..) !! real input
        integer, intent(in) :: sizes !! width shared by all axes
        real(dp), intent(out) :: y(..) !! real output
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary extension
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary modes
        real(dp), optional, intent(in) :: cval !! constant fill value
        integer, optional, intent(in) :: axes(:) !! one-based axes to filter
        if (present(axes)) then
            call uniform_real_work(x,spread(sizes,1,size(axes)),y,status,mode,cval,modes,axes)
        else
            call uniform_real_work(x,spread(sizes,1,rank(x)),y,status,mode,cval,modes)
        end if
    end subroutine uniform_real_scalar

    subroutine uniform_complex_vector(x,sizes,y,status,mode,cval,modes,axes)
        complex(dp), intent(in) :: x(..) !! complex input
        integer, intent(in) :: sizes(:) !! widths along each axis
        complex(dp), intent(out) :: y(..) !! complex output
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary extension
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary modes
        complex(dp), optional, intent(in) :: cval !! constant outside value
        integer, optional, intent(in) :: axes(:) !! one-based axes to filter
        call uniform_complex_work(x,sizes,y,status,mode,cval,modes,axes)
    end subroutine uniform_complex_vector

    subroutine uniform_complex_scalar(x,sizes,y,status,mode,cval,modes,axes)
        complex(dp), intent(in) :: x(..) !! complex input
        integer, intent(in) :: sizes !! common rectangular window width
        complex(dp), intent(out) :: y(..) !! complex output
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary extension
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary modes
        complex(dp), optional, intent(in) :: cval !! constant fill value
        integer, optional, intent(in) :: axes(:) !! one-based axes to filter
        if (present(axes)) then
            call uniform_complex_work(x,spread(sizes,1,size(axes)),y,status,mode,cval,modes,axes)
        else
            call uniform_complex_work(x,spread(sizes,1,rank(x)),y,status,mode,cval,modes)
        end if
    end subroutine uniform_complex_scalar

    subroutine extrema_footprint_process(values,dims,offsets,mode,fill,maximum,result,modes)
        real(dp), intent(in) :: values(:) !! packed source image
        integer, intent(in) :: dims(:) !! image dimensions
        integer, intent(in) :: offsets(:,:) !! selected offsets (rank,active samples)
        character(len=*), intent(in) :: mode !! boundary mode
        real(dp), intent(in) :: fill !! constant-extension value
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary modes
        logical, intent(in) :: maximum !! true to compute maximum
        real(dp), intent(out) :: result(:) !! packed neighborhood extrema
        integer(int64) :: strides(size(dims)),pos,rem,q
        integer :: loc(size(dims)),d,k,t
        real(dp) :: best,value
        logical :: inside
        character(len=:), allocatable :: boundary

        strides(1)=1_int64
        do d=2,size(dims)
            strides(d)=strides(d-1)*int(dims(d-1),int64)
        end do
        do pos=1_int64,size(values,kind=int64)
            rem=pos-1_int64
            do d=1,size(dims)
                loc(d)=int(modulo(rem,int(dims(d),int64)))
                rem=rem/int(dims(d),int64)
            end do
            if (maximum) then
                best=-huge(1.0_dp)
            else
                best=huge(1.0_dp)
            end if
            do k=1,size(offsets,2)
                inside=.true.
                q=1_int64
                do d=1,size(dims)
                    boundary=mode
                    if (present(modes)) boundary=trim(modes(d))
                    t=extension_index(loc(d)+offsets(d,k),dims(d),boundary)
                    if (t==0) then
                        inside=.false.
                        exit
                    end if
                    q=q+int(t-1,int64)*strides(d)
                end do
                value=fill
                if (inside) value=values(q)
                if (maximum) then
                    best=max(best,value)
                else
                    best=min(best,value)
                end if
            end do
            result(pos)=best
        end do
    end subroutine extrema_footprint_process

    subroutine extremum_real_work(x,sizes,y,maximum,status,mode,cval,footprint,origin,modes,axes)
        real(dp), intent(in) :: x(..) !! real source data
        integer, intent(in) :: sizes(:) !! neighborhood widths along axes
        real(dp), intent(out) :: y(..) !! local maxima or minima
        logical, intent(in) :: maximum !! select maximum or minimum
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary extension
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary modes
        real(dp), optional, intent(in) :: cval !! constant outside value
        logical, optional, intent(in) :: footprint(..) !! active locations, overrides sizes
        integer, optional, intent(in) :: origin(:) !! per-axis offsets, SciPy convention
        integer, optional, intent(in) :: axes(:) !! selected rectangular axes
        real(dp), allocatable :: data(:),result(:),heights(:)
        integer, allocatable :: widths(:),offsets(:,:)
        real(dp) :: fill
        integer :: code
        character(len=:), allocatable :: boundary
        character(len=32), allocatable :: boundaries(:)
        logical :: valid_modes

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        if (present(axes)) then
            if (present(footprint) .or. present(origin)) return
        end if
        call expand_selected_int(sizes,rank(x),axes,1,widths,valid_modes)
        if (.not.valid_modes) return
        call pack_real(x,data)
        if (.not.all(ieee_is_finite(data))) return
        fill=0.0_dp
        if (present(cval)) fill=cval
        call selected_filter_modes(rank(x),axes,mode,modes,boundaries,valid_modes)
        if (.not.valid_modes .or. .not.ieee_is_finite(fill)) return
        boundary='reflect' 
        if (present(footprint) .or. present(origin)) then
            call gray_stencil(rank(x),sizes=widths,footprint=footprint,origin=origin, &
                offsets=offsets,heights=heights,status=code)
            if (code==ndimage_status_success) then
                allocate(result(size(data)))
                call extrema_footprint_process(data,shape(x),offsets,boundary,fill,maximum,result,boundaries)
                call scatter_real(result,y)
            end if
        else
            call extrema_process(data,shape(x),widths,boundary,fill,maximum,code,boundaries,axes)
            if (code==ndimage_status_success) call scatter_real(data,y)
        end if
        if (present(status)) status=code
    end subroutine extremum_real_work

    subroutine minimum_vector(x,sizes,y,status,mode,cval,footprint,origin,modes,axes)
        real(dp), intent(in) :: x(..) !! source data
        integer, optional, intent(in) :: sizes(:) !! per-axis widths, default 3
        real(dp), intent(out) :: y(..) !! local minimum image
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary mode
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary modes
        real(dp), optional, intent(in) :: cval !! constant fill value
        logical, optional, intent(in) :: footprint(..) !! active filter window
        integer, optional, intent(in) :: origin(:) !! per-axis stencil origin
        integer, optional, intent(in) :: axes(:) !! selected rectangular axes
        integer, allocatable :: widths(:)
        if (present(axes)) then
            widths=spread(3,1,size(axes))
        else
            widths=spread(3,1,rank(x))
        end if
        if (present(sizes)) widths=sizes
        call extremum_real_work(x,widths,y,.false.,status,mode,cval,footprint,origin,modes,axes)
    end subroutine minimum_vector

    subroutine minimum_scalar(x,sizes,y,status,mode,cval,footprint,origin,modes,axes)
        real(dp), intent(in) :: x(..) !! source data
        integer, intent(in) :: sizes !! common neighborhood width
        real(dp), intent(out) :: y(..) !! local minimum image
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary mode
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary modes
        real(dp), optional, intent(in) :: cval !! constant outside value
        logical, optional, intent(in) :: footprint(..) !! active filter window
        integer, optional, intent(in) :: origin(:) !! per-axis stencil origin
        integer, optional, intent(in) :: axes(:) !! selected rectangular axes
        if (present(axes)) then
            call extremum_real_work(x,spread(sizes,1,size(axes)),y,.false., &
                status,mode,cval,footprint,origin,modes,axes)
        else
            call extremum_real_work(x,spread(sizes,1,rank(x)),y,.false., &
                status,mode,cval,footprint,origin,modes)
        end if
    end subroutine minimum_scalar

    subroutine maximum_vector(x,sizes,y,status,mode,cval,footprint,origin,modes,axes)
        real(dp), intent(in) :: x(..) !! source data
        integer, optional, intent(in) :: sizes(:) !! per-axis widths, default 3
        real(dp), intent(out) :: y(..) !! local maximum image
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary mode
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary modes
        real(dp), optional, intent(in) :: cval !! constant fill value
        logical, optional, intent(in) :: footprint(..) !! active filter window
        integer, optional, intent(in) :: origin(:) !! per-axis stencil origin
        integer, optional, intent(in) :: axes(:) !! selected rectangular axes
        integer, allocatable :: widths(:)
        if (present(axes)) then
            widths=spread(3,1,size(axes))
        else
            widths=spread(3,1,rank(x))
        end if
        if (present(sizes)) widths=sizes
        call extremum_real_work(x,widths,y,.true.,status,mode,cval,footprint,origin,modes,axes)
    end subroutine maximum_vector

    subroutine maximum_scalar(x,sizes,y,status,mode,cval,footprint,origin,modes,axes)
        real(dp), intent(in) :: x(..) !! source data
        integer, intent(in) :: sizes !! common neighborhood width
        real(dp), intent(out) :: y(..) !! local maximum image
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary mode
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary modes
        real(dp), optional, intent(in) :: cval !! constant outside value
        logical, optional, intent(in) :: footprint(..) !! active filter window
        integer, optional, intent(in) :: origin(:) !! per-axis stencil origin
        integer, optional, intent(in) :: axes(:) !! selected rectangular axes
        if (present(axes)) then
            call extremum_real_work(x,spread(sizes,1,size(axes)),y,.true., &
                status,mode,cval,footprint,origin,modes,axes)
        else
            call extremum_real_work(x,spread(sizes,1,rank(x)),y,.true., &
                status,mode,cval,footprint,origin,modes)
        end if
    end subroutine maximum_scalar

    subroutine derivative_image(x,sigma,y,laplacian,status,mode,cval,truncate)
        real(dp), intent(in) :: x(..) !! source real image
        real(dp), intent(in) :: sigma(:) !! Gaussian widths for all axes
        real(dp), intent(out) :: y(..) !! gradient magnitude or Laplacian
        logical, intent(in) :: laplacian !! true: sum 2nd derivatives
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary extension
        real(dp), optional, intent(in) :: cval !! constant fill value
        real(dp), optional, intent(in) :: truncate !! radius in sigmas
        complex(dp), allocatable :: base(:), work(:)
        real(dp), allocatable :: data(:),accum(:)
        integer, allocatable :: orders(:),radii(:)
        complex(dp) :: fill
        real(dp) :: cutoff
        character(len=:), allocatable :: boundary
        integer :: d,code

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(x)/=rank(y)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        if (size(sigma)/=rank(x) .or. any(sigma<=0.0_dp)) return
        if (.not.all(ieee_is_finite(sigma))) return
        call pack_real(x,data)
        if (.not.all(ieee_is_finite(data))) return
        base=cmplx(data,0.0_dp,dp)
        allocate(accum(size(data)))
        accum=0.0_dp
        orders=spread(0,1,rank(x))
        radii=spread(-1,1,rank(x))
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cmplx(cval,0.0_dp,dp)
        boundary='reflect'
        if (present(mode)) boundary=trim(mode)
        cutoff=4.0_dp
        if (present(truncate)) cutoff=truncate
        do d=1,rank(x)
            orders=0
            if (laplacian) then
                orders(d)=2
            else
                orders(d)=1
            end if
            work=base
            call gaussian_process(work,shape(x),sigma,orders,boundary,fill,cutoff,radii,code)
            if (code/=ndimage_status_success) then
                if (present(status)) status=code
                return
            end if
            if (laplacian) then
                accum=accum+real(work,dp)
            else
                accum=accum+real(work,dp)**2
            end if
        end do
        if (.not.laplacian) accum=sqrt(accum)
        call scatter_real(accum,y)
        if (present(status)) status=ndimage_status_success
    end subroutine derivative_image

    subroutine gradient_magnitude_vector(x,sigma,y,status,mode,cval,truncate)
        real(dp), intent(in) :: x(..) !! source image
        real(dp), intent(in) :: sigma(:) !! per-axis standard deviations
        real(dp), intent(out) :: y(..) !! Gaussian gradient magnitude
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary mode
        real(dp), optional, intent(in) :: cval !! constant outside value
        real(dp), optional, intent(in) :: truncate !! cutoff in sigmas
        call derivative_image(x,sigma,y,.false.,status,mode,cval,truncate)
    end subroutine gradient_magnitude_vector

    subroutine gradient_magnitude_scalar(x,sigma,y,status,mode,cval,truncate)
        real(dp), intent(in) :: x(..) !! source image
        real(dp), intent(in) :: sigma !! common Gaussian width
        real(dp), intent(out) :: y(..) !! gradient magnitude
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary mode
        real(dp), optional, intent(in) :: cval !! constant outside value
        real(dp), optional, intent(in) :: truncate !! cutoff in sigmas
        call derivative_image(x,spread(sigma,1,rank(x)),y,.false., &
                              status,mode,cval,truncate)
    end subroutine gradient_magnitude_scalar

    subroutine gaussian_laplace_vector(x,sigma,y,status,mode,cval,truncate)
        real(dp), intent(in) :: x(..) !! source image
        real(dp), intent(in) :: sigma(:) !! per-axis Gaussian widths
        real(dp), intent(out) :: y(..) !! Laplacian of Gaussian
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary mode
        real(dp), optional, intent(in) :: cval !! constant outside value
        real(dp), optional, intent(in) :: truncate !! cutoff in sigmas
        call derivative_image(x,sigma,y,.true.,status,mode,cval,truncate)
    end subroutine gaussian_laplace_vector

    subroutine gaussian_laplace_scalar(x,sigma,y,status,mode,cval,truncate)
        real(dp), intent(in) :: x(..) !! source image
        real(dp), intent(in) :: sigma !! common Gaussian width
        real(dp), intent(out) :: y(..) !! Gaussian Laplacian
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary mode
        real(dp), optional, intent(in) :: cval !! constant outside value
        real(dp), optional, intent(in) :: truncate !! cutoff in sigmas
        call derivative_image(x,spread(sigma,1,rank(x)),y,.true., &
                              status,mode,cval,truncate)
    end subroutine gaussian_laplace_scalar

    ! The caller supplies an arbitrary axis derivative.  The user callback
    ! writes into y as scratch for each axis; the combined result is written
    ! to y only after all callbacks finish. This avoids rank-specific buffers.
    subroutine generic_derivatives(x,derivative,y,status,mode,cval,modes,laplacian,axes)
        real(dp), intent(in) :: x(..) !! sampled real image
        procedure(ndimage_derivative_func) :: derivative !! axis derivative callback
        real(dp), intent(out) :: y(..) !! gradient magnitude or Laplacian
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! shared boundary extension
        real(dp), optional, intent(in) :: cval !! fill value
        character(len=*), optional, intent(in) :: modes(:) !! optional per-axis modes
        logical, intent(in) :: laplacian !! sum terms versus square and square-root
        integer, optional, intent(in) :: axes(:) !! one-based axes to differentiate, without duplicates
        character(len=32), allocatable :: boundaries(:)
        real(dp), allocatable :: accum(:),term(:),input_flat(:)
        integer, allocatable :: selected(:)
        real(dp) :: fill
        logical :: valid
        integer :: d,code,i

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(x)/=rank(y)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        fill=0.0_dp
        if (present(cval)) fill=cval
        if (.not.ieee_is_finite(fill)) return
        selected=[(d,d=1,rank(x))]
        if (present(axes)) then
            if (size(axes)==0) return
            if (any(axes<1) .or. any(axes>rank(x))) return
            do i=1,size(axes)
                if (count(axes==axes(i))/=1) return
            end do
            selected=axes
        end if
        if (present(modes) .and. present(axes)) then
            ! Accept a compact mode sequence in the order of the chosen axes.
            ! Full-rank modes(:) retains the original SciFort axis indexing.
            if (size(modes)==size(selected) .and. size(selected)/=rank(x)) then
                if (present(mode)) return
                allocate(boundaries(rank(x)))
                boundaries='reflect'
                do i=1,size(selected)
                    if (.not.mode_valid(trim(modes(i)))) return
                    boundaries(selected(i))=trim(modes(i))
                end do
            else
                call resolve_axis_modes(rank(x),mode,modes,boundaries,valid)
                if (.not.valid) return
            end if
        else
            call resolve_axis_modes(rank(x),mode,modes,boundaries,valid)
            if (.not.valid) return
        end if
        call pack_real(x,input_flat)
        if (.not.all(ieee_is_finite(input_flat))) return
        allocate(accum(size(input_flat)))
        accum=0.0_dp
        do i=1,size(selected)
            d=selected(i)
            call derivative(x,d,y,trim(boundaries(d)),fill,code)
            if (code/=ndimage_status_success) then
                if (present(status)) status=code
                return
            end if
            call pack_real(y,term)
            if (.not.all(ieee_is_finite(term))) return
            if (laplacian) then
                accum=accum+term
            else
                accum=accum+term**2
            end if
        end do
        if (.not.laplacian) accum=sqrt(accum)
        call scatter_real(accum,y)
        if (present(status)) status=ndimage_status_success
    end subroutine generic_derivatives

    subroutine generic_gradient_magnitude_real_nd(x,derivative,y,status,mode,cval,modes,axes)
        real(dp), intent(in) :: x(..) !! source image
        procedure(ndimage_derivative_func) :: derivative !! first derivative along one axis
        real(dp), intent(out) :: y(..) !! Euclidean norm of first derivatives
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! common boundary extension
        real(dp), optional, intent(in) :: cval !! constant outside value
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary extensions
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        call generic_derivatives(x,derivative,y,status,mode,cval,modes,.false.,axes)
    end subroutine generic_gradient_magnitude_real_nd

    subroutine generic_laplace_real_nd(x,derivative2,y,status,mode,cval,modes,axes)
        real(dp), intent(in) :: x(..) !! source image
        procedure(ndimage_derivative_func) :: derivative2 !! second axis derivative callback
        real(dp), intent(out) :: y(..) !! sum of second derivatives
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! common boundary extension
        real(dp), optional, intent(in) :: cval !! constant outside value
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundary extensions
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        call generic_derivatives(x,derivative2,y,status,mode,cval,modes,.true.,axes)
    end subroutine generic_laplace_real_nd

    ! Complex counterpart to generic_derivatives.  Following scipy.ndimage,
    ! gradient magnitude is sqrt(sum(derivative**2)), without conjugation.
    ! Consequently the output is complex rather than a Hermitian norm.
    subroutine generic_complex_derivatives(x,derivative,y,status,mode,cval,modes,laplacian,axes)
        complex(dp), intent(in) :: x(..) !! complex image, ranks 1 through 15
        procedure(ndimage_complex_derivative_func) :: derivative !! complex callback
        complex(dp), intent(out) :: y(..) !! complex gradient magnitude or Laplacian
        integer, optional, intent(out) :: status !! 0 success, 1 invalid, callback error passed through
        character(len=*), optional, intent(in) :: mode !! common boundary extension
        complex(dp), optional, intent(in) :: cval !! complex outside value
        character(len=*), optional, intent(in) :: modes(:) !! per-axis modes
        logical, intent(in) :: laplacian !! sum derivatives instead of their squares
        integer, optional, intent(in) :: axes(:) !! selected one-based axes, no duplicates
        character(len=32), allocatable :: boundaries(:)
        complex(dp), allocatable :: accum(:),term(:),input_flat(:)
        integer, allocatable :: selected(:)
        complex(dp) :: fill
        logical :: valid
        integer :: d,code,i

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(x)/=rank(y)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        fill=cmplx(0.0_dp,0.0_dp,dp)
        if (present(cval)) fill=cval
        if (.not.ieee_is_finite(real(fill,dp)) .or. .not.ieee_is_finite(aimag(fill))) return
        selected=[(d,d=1,rank(x))]
        if (present(axes)) then
            if (size(axes)==0) return
            if (any(axes<1) .or. any(axes>rank(x))) return
            do i=1,size(axes)
                if (count(axes==axes(i))/=1) return
            end do
            selected=axes
        end if
        if (present(modes) .and. present(axes)) then
            if (size(modes)==size(selected) .and. size(selected)/=rank(x)) then
                if (present(mode)) return
                allocate(boundaries(rank(x)))
                boundaries='reflect'
                do i=1,size(selected)
                    if (.not.mode_valid(trim(modes(i)))) return
                    boundaries(selected(i))=trim(modes(i))
                end do
            else
                call resolve_axis_modes(rank(x),mode,modes,boundaries,valid)
                if (.not.valid) return
            end if
        else
            call resolve_axis_modes(rank(x),mode,modes,boundaries,valid)
            if (.not.valid) return
        end if
        call pack_complex(x,input_flat)
        if (.not.all(ieee_is_finite(real(input_flat,dp)))) return
        if (.not.all(ieee_is_finite(aimag(input_flat)))) return
        allocate(accum(size(input_flat)))
        accum=cmplx(0.0_dp,0.0_dp,dp)
        do i=1,size(selected)
            d=selected(i)
            call derivative(x,d,y,trim(boundaries(d)),fill,code)
            if (code/=ndimage_status_success) then
                if (present(status)) status=code
                return
            end if
            call pack_complex(y,term)
            if (.not.all(ieee_is_finite(real(term,dp)))) return
            if (.not.all(ieee_is_finite(aimag(term)))) return
            if (laplacian) then
                accum=accum+term
            else
                accum=accum+term**2
            end if
        end do
        if (.not.laplacian) accum=sqrt(accum)
        call scatter_complex(accum,y)
        if (present(status)) status=ndimage_status_success
    end subroutine generic_complex_derivatives

    subroutine generic_gradient_magnitude_complex_nd(x,derivative,y,status,mode,cval,modes,axes)
        complex(dp), intent(in) :: x(..) !! complex image
        procedure(ndimage_complex_derivative_func) :: derivative !! first derivative callback
        complex(dp), intent(out) :: y(..) !! SciPy-style complex square-root of squared derivatives
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! common mode
        complex(dp), optional, intent(in) :: cval !! constant outside value
        character(len=*), optional, intent(in) :: modes(:) !! per-axis modes
        integer, optional, intent(in) :: axes(:) !! one-based axes
        call generic_complex_derivatives(x,derivative,y,status,mode,cval,modes,.false.,axes)
    end subroutine generic_gradient_magnitude_complex_nd

    subroutine generic_laplace_complex_nd(x,derivative2,y,status,mode,cval,modes,axes)
        complex(dp), intent(in) :: x(..) !! complex image
        procedure(ndimage_complex_derivative_func) :: derivative2 !! second derivative callback
        complex(dp), intent(out) :: y(..) !! sum of second derivatives
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! common mode
        complex(dp), optional, intent(in) :: cval !! constant outside value
        character(len=*), optional, intent(in) :: modes(:) !! per-axis modes
        integer, optional, intent(in) :: axes(:) !! one-based axes
        call generic_complex_derivatives(x,derivative2,y,status,mode,cval,modes,.true.,axes)
    end subroutine generic_laplace_complex_nd

    subroutine edge_filter_work(x,axis,y,sobel,status,mode,cval)
        real(dp), intent(in) :: x(..) !! sampled real image
        integer, intent(in) :: axis !! one-based derivative axis
        real(dp), intent(out) :: y(..) !! weighted first difference
        logical, intent(in) :: sobel !! select Sobel versus Prewitt smoothing
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! default reflect
        real(dp), optional, intent(in) :: cval !! constant outside value
        real(dp), allocatable :: data(:)
        complex(dp), allocatable :: work(:)
        complex(dp) :: fill
        character(len=:), allocatable :: boundary
        integer :: d

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(x)/=rank(y)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (axis<1 .or. axis>rank(x)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        boundary='reflect'
        if (present(mode)) boundary=trim(mode)
        if (.not.mode_valid(boundary)) return
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cmplx(cval,0.0_dp,dp)
        if (.not.ieee_is_finite(real(fill,dp))) return
        call pack_real(x,data)
        if (.not.all(ieee_is_finite(data))) return
        work=cmplx(data,0.0_dp,dp)
        call axis_convolution(work,shape(x),axis,[-1,0,1], &
                              [-1.0_dp,0.0_dp,1.0_dp],boundary,fill)
        do d=1,rank(x)
            if (d==axis) cycle
            if (sobel) then
                call axis_convolution(work,shape(x),d,[-1,0,1], &
                                      [1.0_dp,2.0_dp,1.0_dp],boundary,fill)
            else
                call axis_convolution(work,shape(x),d,[-1,0,1], &
                                      [1.0_dp,1.0_dp,1.0_dp],boundary,fill)
            end if
        end do
        data=real(work,dp)
        call scatter_real(data,y)
        if (present(status)) status=ndimage_status_success
    end subroutine edge_filter_work

    subroutine sobel_nd(x,axis,y,status,mode,cval)
        real(dp), intent(in) :: x(..) !! sampled real image
        integer, intent(in) :: axis !! one-based derivative axis
        real(dp), intent(out) :: y(..) !! Sobel edge response
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary mode
        real(dp), optional, intent(in) :: cval !! outside-image value
        call edge_filter_work(x,axis,y,.true.,status,mode,cval)
    end subroutine sobel_nd

    subroutine prewitt_nd(x,axis,y,status,mode,cval)
        real(dp), intent(in) :: x(..) !! sampled real image
        integer, intent(in) :: axis !! one-based derivative axis
        real(dp), intent(out) :: y(..) !! Prewitt edge response
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary mode
        real(dp), optional, intent(in) :: cval !! outside-image value
        call edge_filter_work(x,axis,y,.false.,status,mode,cval)
    end subroutine prewitt_nd

    subroutine laplace_nd(x,y,status,mode,cval)
        real(dp), intent(in) :: x(..) !! sampled real image
        real(dp), intent(out) :: y(..) !! sum of axis second differences
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! default reflect
        real(dp), optional, intent(in) :: cval !! constant outside value
        real(dp), allocatable :: data(:),accum(:)
        complex(dp), allocatable :: work(:),base(:)
        complex(dp) :: fill
        character(len=:), allocatable :: boundary
        integer :: d

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(x)/=rank(y)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        boundary='reflect'
        if (present(mode)) boundary=trim(mode)
        if (.not.mode_valid(boundary)) return
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cmplx(cval,0.0_dp,dp)
        if (.not.ieee_is_finite(real(fill,dp))) return
        call pack_real(x,data)
        if (.not.all(ieee_is_finite(data))) return
        base=cmplx(data,0.0_dp,dp)
        allocate(accum(size(data)))
        accum=0.0_dp
        do d=1,rank(x)
            work=base
            call axis_convolution(work,shape(x),d,[-1,0,1], &
                                  [1.0_dp,-2.0_dp,1.0_dp],boundary,fill)
            accum=accum+real(work,dp)
        end do
        call scatter_real(accum,y)
        if (present(status)) status=ndimage_status_success
    end subroutine laplace_nd

    ! A full-dimensional correlation/convolution; the kernel is not required
    ! to be separable, symmetric, or of odd length along any axis.
    subroutine correlate_process(values,dims,weights,wdims,mode,fill,convolution,result,modes)
        complex(dp), intent(in) :: values(:) !! flattened source image
        integer, intent(in) :: dims(:) !! source extents
        complex(dp), intent(in) :: weights(:) !! flattened kernel
        integer, intent(in) :: wdims(:) !! kernel extents
        character(len=*), intent(in) :: mode !! boundary mode
        complex(dp), intent(in) :: fill !! value outside image
        character(len=*), optional, intent(in) :: modes(:) !! per-axis modes
        logical, intent(in) :: convolution !! use convolution instead of correlation
        complex(dp), intent(out) :: result(:) !! flattened filtered output
        integer(int64) :: stride(size(dims)),pos,q,rem,kw
        integer :: coords(size(dims)),d,t,k,offset
        complex(dp) :: acc,w
        logical :: valid
        character(len=:), allocatable :: boundary

        stride(1)=1_int64
        do d=2,size(dims)
            stride(d)=stride(d-1)*int(dims(d-1),int64)
        end do
        do pos=1_int64,size(values,kind=int64)
            rem=pos-1_int64
            do d=1,size(dims)
                coords(d)=int(modulo(rem,int(dims(d),int64)))
                rem=rem/int(dims(d),int64)
            end do
            acc=(0.0_dp,0.0_dp)
            do k=1,size(weights)
                kw=int(k-1,int64)
                valid=.true.
                q=1_int64
                do d=1,size(dims)
                    t=int(modulo(kw,int(wdims(d),int64)))
                    kw=kw/int(wdims(d),int64)
                    if (convolution) then
                        offset=wdims(d)/2-t
                    else
                        offset=t-wdims(d)/2
                    end if
                    boundary=mode
                    if (present(modes)) boundary=trim(modes(d))
                    t=extension_index(coords(d)+offset,dims(d),boundary)
                    if (t==0) then
                        valid=.false.
                        exit
                    end if
                    q=q+int(t-1,int64)*stride(d)
                end do
                w=weights(k)
                if (.not.convolution) w=conjg(w)
                if (valid) then
                    acc=acc+w*values(q)
                else
                    acc=acc+w*fill
                end if
            end do
            result(pos)=acc
        end do
    end subroutine correlate_process

    subroutine correlate_real_work(x,weights,y,convolution,status,mode,cval,modes)
        real(dp), intent(in) :: x(..) !! real source of ranks 1 through 15
        real(dp), intent(in) :: weights(..) !! real full-rank kernel
        real(dp), intent(out) :: y(..) !! same-shape real result
        logical, intent(in) :: convolution !! reverse kernel sense
        integer, optional, intent(out) :: status !! validation status
        character(len=*), optional, intent(in) :: mode !! common boundary mode
        character(len=*), optional, intent(in) :: modes(:) !! per-axis modes
        real(dp), optional, intent(in) :: cval !! outside-image value
        real(dp), allocatable :: flat(:),coeff(:)
        complex(dp), allocatable :: res(:)
        real(dp) :: fill
        character(len=:), allocatable :: boundary
        character(len=32), allocatable :: boundaries(:)
        logical :: valid_modes

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x) .or. rank(weights)/=rank(x)) return
        if (any(shape(y)/=shape(x)) .or. any(shape(weights)<1) .or. any(shape(x)<1)) return
        if (size(weights,kind=int64)>50000_int64) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        call resolve_axis_modes(rank(x),mode,modes,boundaries,valid_modes)
        if (.not.valid_modes) return
        boundary='reflect' 
        fill=0.0_dp
        if (present(cval)) fill=cval
        if (.not.ieee_is_finite(fill)) return
        call pack_real(x,flat)
        call pack_real(weights,coeff)
        if (.not.all(ieee_is_finite(flat)) .or. .not.all(ieee_is_finite(coeff))) return
        allocate(res(size(flat)))
        call correlate_process(cmplx(flat,0.0_dp,dp),shape(x),cmplx(coeff,0.0_dp,dp), &
            shape(weights),boundary,cmplx(fill,0.0_dp,dp),convolution,res,boundaries)
        call scatter_real(real(res,dp),y)
        if (present(status)) status=ndimage_status_success
    end subroutine correlate_real_work

    subroutine correlate_complex_work(x,weights,y,convolution,status,mode,cval,modes)
        complex(dp), intent(in) :: x(..) !! complex source of ranks 1 through 15
        complex(dp), intent(in) :: weights(..) !! complex full-rank kernel
        complex(dp), intent(out) :: y(..) !! same-shape complex result
        logical, intent(in) :: convolution !! reverse kernel sense
        integer, optional, intent(out) :: status !! validation status
        character(len=*), optional, intent(in) :: mode !! boundary mode
        character(len=*), optional, intent(in) :: modes(:) !! per-axis modes
        complex(dp), optional, intent(in) :: cval !! outside-image value
        complex(dp), allocatable :: flat(:),coeff(:),res(:)
        complex(dp) :: fill
        character(len=:), allocatable :: boundary
        character(len=32), allocatable :: boundaries(:)
        logical :: valid_modes

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x) .or. rank(weights)/=rank(x)) return
        if (any(shape(y)/=shape(x)) .or. any(shape(weights)<1) .or. any(shape(x)<1)) return
        if (size(weights,kind=int64)>50000_int64) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        call resolve_axis_modes(rank(x),mode,modes,boundaries,valid_modes)
        if (.not.valid_modes) return
        boundary='reflect' 
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cval
        if (.not.mode_valid(boundary) .or. .not.ieee_is_finite(real(fill,dp)) .or. &
            .not.ieee_is_finite(aimag(fill))) return
        call pack_complex(x,flat)
        call pack_complex(weights,coeff)
        if (.not.all(ieee_is_finite(real(flat,dp))) .or. &
            .not.all(ieee_is_finite(aimag(flat)))) return
        if (.not.all(ieee_is_finite(real(coeff,dp))) .or. &
            .not.all(ieee_is_finite(aimag(coeff)))) return
        allocate(res(size(flat)))
        call correlate_process(flat,shape(x),coeff,shape(weights),boundary,fill,convolution,res,boundaries)
        call scatter_complex(res,y)
        if (present(status)) status=ndimage_status_success
    end subroutine correlate_complex_work

    subroutine correlate_real_nd(x,weights,y,status,mode,cval,modes)
        real(dp), intent(in) :: x(..) !! source image
        real(dp), intent(in) :: weights(..) !! correlation kernel
        real(dp), intent(out) :: y(..) !! result
        integer, optional, intent(out) :: status !! status code
        character(len=*), optional, intent(in) :: mode !! boundary mode
        character(len=*), optional, intent(in) :: modes(:) !! per-axis modes
        real(dp), optional, intent(in) :: cval !! constant outside value
        call correlate_real_work(x,weights,y,.false.,status,mode,cval,modes)
    end subroutine correlate_real_nd

    subroutine convolve_real_nd(x,weights,y,status,mode,cval,modes)
        real(dp), intent(in) :: x(..) !! source image
        real(dp), intent(in) :: weights(..) !! convolution kernel
        real(dp), intent(out) :: y(..) !! result
        integer, optional, intent(out) :: status !! status code
        character(len=*), optional, intent(in) :: mode !! boundary mode
        character(len=*), optional, intent(in) :: modes(:) !! per-axis modes
        real(dp), optional, intent(in) :: cval !! constant outside value
        call correlate_real_work(x,weights,y,.true.,status,mode,cval,modes)
    end subroutine convolve_real_nd

    subroutine correlate_complex_nd(x,weights,y,status,mode,cval,modes)
        complex(dp), intent(in) :: x(..) !! source image
        complex(dp), intent(in) :: weights(..) !! complex correlation kernel
        complex(dp), intent(out) :: y(..) !! result
        integer, optional, intent(out) :: status !! status code
        character(len=*), optional, intent(in) :: mode !! boundary mode
        character(len=*), optional, intent(in) :: modes(:) !! per-axis modes
        complex(dp), optional, intent(in) :: cval !! complex outside value
        call correlate_complex_work(x,weights,y,.false.,status,mode,cval,modes)
    end subroutine correlate_complex_nd

    subroutine convolve_complex_nd(x,weights,y,status,mode,cval,modes)
        complex(dp), intent(in) :: x(..) !! source image
        complex(dp), intent(in) :: weights(..) !! complex convolution kernel
        complex(dp), intent(out) :: y(..) !! result
        integer, optional, intent(out) :: status !! status code
        character(len=*), optional, intent(in) :: mode !! boundary mode
        character(len=*), optional, intent(in) :: modes(:) !! per-axis modes
        complex(dp), optional, intent(in) :: cval !! complex outside value
        call correlate_complex_work(x,weights,y,.true.,status,mode,cval,modes)
    end subroutine convolve_complex_nd

    ! One-axis arbitrary real or complex correlation/convolution. SciPy has
    ! different even-kernel centering conventions for the two operations.
    subroutine filter1d_process(x,dims,weights,axis,origin,mode,fill,convolution,y)
        complex(dp), intent(in) :: x(:) !! flattened complex image
        integer, intent(in) :: dims(:) !! image dimensions
        complex(dp), intent(in) :: weights(:) !! complex coefficients
        integer, intent(in) :: axis !! one-based processing axis
        integer, intent(in) :: origin !! kernel-center offset
        character(len=*), intent(in) :: mode !! border extension
        complex(dp), intent(in) :: fill !! outside value
        logical, intent(in) :: convolution !! apply convolution rather than correlation
        complex(dp), intent(out) :: y(:) !! flattened result
        integer(int64) :: stride,block,outer,o,j,pos,q
        integer :: i,k,offset,t,n,width
        complex(dp) :: acc,coeff

        n=dims(axis)
        width=size(weights)
        stride=product(int(dims(:axis-1),int64))
        block=stride*int(n,int64)
        outer=product(int(dims(axis+1:),int64))
        do o=0_int64,outer-1_int64
            do j=0_int64,stride-1_int64
                do i=0,n-1
                    acc=(0.0_dp,0.0_dp)
                    do k=1,width
                        if (convolution) then
                            offset=width/2-(k-1)+origin
                        else
                            offset=(k-1)-width/2-origin
                        end if
                        t=extension_index(i+offset,n,mode)
                        coeff=weights(k)
                        if (.not.convolution) coeff=conjg(coeff)
                        if (t==0) then
                            acc=acc+coeff*fill
                        else
                            q=o*block+j+int(t-1,int64)*stride+1_int64
                            acc=acc+coeff*x(q)
                        end if
                    end do
                    pos=o*block+j+int(i,int64)*stride+1_int64
                    y(pos)=acc
                end do
            end do
        end do
    end subroutine filter1d_process

    subroutine filter1d_complex_work(x,weights,axis,y,convolution,status,mode,cval,origin)
        complex(dp), intent(in) :: x(..) !! complex input image
        complex(dp), intent(in) :: weights(:) !! filter coefficients
        integer, intent(in) :: axis !! one-based axis (1..rank)
        complex(dp), intent(out) :: y(..) !! same-shaped output
        logical, intent(in) :: convolution !! reverse input/kernel offset
        integer, optional, intent(out) :: status !! 0 successful, 1 invalid, 2 resource limit
        character(len=*), optional, intent(in) :: mode !! extension type
        complex(dp), optional, intent(in) :: cval !! constant fill value
        integer, optional, intent(in) :: origin !! shifted center, SciPy convention
        complex(dp), allocatable :: flat(:),res(:)
        complex(dp) :: fill
        integer :: shift
        character(len=:), allocatable :: boundary

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (axis<1 .or. axis>rank(x) .or. size(weights)<1) return
        if (size(weights)>50000) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        shift=0
        if (present(origin)) shift=origin
        if (shift<-(size(weights)/2) .or. shift>(size(weights)-1)/2) return
        boundary='reflect'
        if (present(mode)) boundary=trim(mode)
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cval
        if (.not.mode_valid(boundary)) return
        if (.not.all(ieee_is_finite(real(weights,dp))) .or. &
            .not.all(ieee_is_finite(aimag(weights)))) return
        if (.not.ieee_is_finite(real(fill,dp)) .or. .not.ieee_is_finite(aimag(fill))) return
        call pack_complex(x,flat)
        if (.not.all(ieee_is_finite(real(flat,dp))) .or. &
            .not.all(ieee_is_finite(aimag(flat)))) return
        allocate(res(size(flat)))
        call filter1d_process(flat,shape(x),weights,axis,shift,boundary,fill,convolution,res)
        call scatter_complex(res,y)
        if (present(status)) status=ndimage_status_success
    end subroutine filter1d_complex_work

    subroutine filter1d_real_work(x,weights,axis,y,convolution,status,mode,cval,origin)
        real(dp), intent(in) :: x(..) !! real input image
        real(dp), intent(in) :: weights(:) !! real filter coefficients
        integer, intent(in) :: axis !! one-based axis
        real(dp), intent(out) :: y(..) !! same-shaped result
        logical, intent(in) :: convolution !! use convolution conventions
        integer, optional, intent(out) :: status !! 0 success, 1 invalid, 2 oversized kernel
        character(len=*), optional, intent(in) :: mode !! extension
        real(dp), optional, intent(in) :: cval !! constant exterior value
        integer, optional, intent(in) :: origin !! center displacement
        real(dp), allocatable :: flat(:)
        complex(dp), allocatable :: res(:)
        real(dp) :: fill
        integer :: shift
        character(len=:), allocatable :: boundary

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (axis<1 .or. axis>rank(x) .or. size(weights)<1) return
        if (size(weights)>50000) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        shift=0
        if (present(origin)) shift=origin
        if (shift<-(size(weights)/2) .or. shift>(size(weights)-1)/2) return
        boundary='reflect'
        if (present(mode)) boundary=trim(mode)
        fill=0.0_dp
        if (present(cval)) fill=cval
        if (.not.mode_valid(boundary)) return
        if (.not.all(ieee_is_finite(weights)) .or. .not.ieee_is_finite(fill)) return
        call pack_real(x,flat)
        if (.not.all(ieee_is_finite(flat))) return
        allocate(res(size(flat)))
        call filter1d_process(cmplx(flat,0.0_dp,dp),shape(x),cmplx(weights,0.0_dp,dp), &
            axis,shift,boundary,cmplx(fill,0.0_dp,dp),convolution,res)
        call scatter_real(real(res,dp),y)
        if (present(status)) status=ndimage_status_success
    end subroutine filter1d_real_work

    subroutine correlate1d_real_nd(x,weights,axis,y,status,mode,cval,origin)
        real(dp), intent(in) :: x(..) !! real input
        real(dp), intent(in) :: weights(:) !! correlation coefficients
        integer, intent(in) :: axis !! one-based axis
        real(dp), intent(out) :: y(..) !! output
        integer, optional, intent(out) :: status !! validation status
        character(len=*), optional, intent(in) :: mode !! border mode
        real(dp), optional, intent(in) :: cval !! constant value
        integer, optional, intent(in) :: origin !! shifted center
        call filter1d_real_work(x,weights,axis,y,.false.,status,mode,cval,origin)
    end subroutine correlate1d_real_nd

    subroutine convolve1d_real_nd(x,weights,axis,y,status,mode,cval,origin)
        real(dp), intent(in) :: x(..) !! real input
        real(dp), intent(in) :: weights(:) !! convolution coefficients
        integer, intent(in) :: axis !! one-based axis
        real(dp), intent(out) :: y(..) !! output
        integer, optional, intent(out) :: status !! validation status
        character(len=*), optional, intent(in) :: mode !! border mode
        real(dp), optional, intent(in) :: cval !! constant value
        integer, optional, intent(in) :: origin !! shifted center
        call filter1d_real_work(x,weights,axis,y,.true.,status,mode,cval,origin)
    end subroutine convolve1d_real_nd

    subroutine correlate1d_complex_nd(x,weights,axis,y,status,mode,cval,origin)
        complex(dp), intent(in) :: x(..) !! complex input
        complex(dp), intent(in) :: weights(:) !! correlation coefficients
        integer, intent(in) :: axis !! one-based axis
        complex(dp), intent(out) :: y(..) !! output
        integer, optional, intent(out) :: status !! validation status
        character(len=*), optional, intent(in) :: mode !! border mode
        complex(dp), optional, intent(in) :: cval !! complex constant value
        integer, optional, intent(in) :: origin !! shifted center
        call filter1d_complex_work(x,weights,axis,y,.false.,status,mode,cval,origin)
    end subroutine correlate1d_complex_nd

    subroutine convolve1d_complex_nd(x,weights,axis,y,status,mode,cval,origin)
        complex(dp), intent(in) :: x(..) !! complex input
        complex(dp), intent(in) :: weights(:) !! convolution coefficients
        integer, intent(in) :: axis !! one-based axis
        complex(dp), intent(out) :: y(..) !! output
        integer, optional, intent(out) :: status !! validation status
        character(len=*), optional, intent(in) :: mode !! border mode
        complex(dp), optional, intent(in) :: cval !! complex constant value
        integer, optional, intent(in) :: origin !! shifted center
        call filter1d_complex_work(x,weights,axis,y,.true.,status,mode,cval,origin)
    end subroutine convolve1d_complex_nd

    ! One-axis uniform averaging uses the same correlation kernel and origin
    ! convention as scipy.ndimage.uniform_filter1d.
    subroutine uniform1d_real(x,filter_size,axis,y,status,mode,cval,origin)
        real(dp), intent(in) :: x(..) !! source real array
        integer, intent(in) :: filter_size !! window width, positive
        integer, intent(in) :: axis !! one-based selected axis
        real(dp), intent(out) :: y(..) !! filtered result
        integer, optional, intent(out) :: status !! status code
        character(len=*), optional, intent(in) :: mode !! extension mode
        real(dp), optional, intent(in) :: cval !! constant outside value
        integer, optional, intent(in) :: origin !! window-center displacement
        real(dp), allocatable :: flat(:)
        complex(dp), allocatable :: data(:)
        character(len=:), allocatable :: boundary
        complex(dp) :: fill
        integer :: shift

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (axis<1 .or. axis>rank(x) .or. filter_size<1) return
        shift=0
        if (present(origin)) shift=origin
        if (shift<-(filter_size/2) .or. shift>(filter_size-1)/2) return
        if (filter_size>50000 .or. size(x,kind=int64)>int(huge(1),int64) .or. &
            int(size(x,axis),int64)+int(filter_size,int64)-1_int64>int(huge(1),int64)) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        boundary='reflect'
        if (present(mode)) boundary=trim(mode)
        if (.not.mode_valid(boundary)) return
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cmplx(cval,0.0_dp,dp)
        if (.not.ieee_is_finite(real(fill,dp))) return
        call pack_real(x,flat)
        if (.not.all(ieee_is_finite(flat))) return
        data=cmplx(flat,0.0_dp,dp)
        call uniform_axis_running(data,shape(x),axis,filter_size,boundary,fill,shift)
        call scatter_real(real(data,dp),y)
        if (present(status)) status=ndimage_status_success
    end subroutine uniform1d_real

    subroutine uniform1d_complex(x,filter_size,axis,y,status,mode,cval,origin)
        complex(dp), intent(in) :: x(..) !! source complex array
        integer, intent(in) :: filter_size !! positive window size
        integer, intent(in) :: axis !! one-based selected axis
        complex(dp), intent(out) :: y(..) !! filtered result
        integer, optional, intent(out) :: status !! status code
        character(len=*), optional, intent(in) :: mode !! extension mode
        complex(dp), optional, intent(in) :: cval !! constant outside value
        integer, optional, intent(in) :: origin !! window-center displacement
        complex(dp), allocatable :: data(:)
        character(len=:), allocatable :: boundary
        complex(dp) :: fill
        integer :: shift

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (axis<1 .or. axis>rank(x) .or. filter_size<1) return
        shift=0
        if (present(origin)) shift=origin
        if (shift<-(filter_size/2) .or. shift>(filter_size-1)/2) return
        if (filter_size>50000 .or. size(x,kind=int64)>int(huge(1),int64) .or. &
            int(size(x,axis),int64)+int(filter_size,int64)-1_int64>int(huge(1),int64)) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        boundary='reflect'
        if (present(mode)) boundary=trim(mode)
        if (.not.mode_valid(boundary)) return
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cval
        if (.not.ieee_is_finite(real(fill,dp)) .or. .not.ieee_is_finite(aimag(fill))) return
        call pack_complex(x,data)
        if (.not.all(ieee_is_finite(real(data,dp))) .or. .not.all(ieee_is_finite(aimag(data)))) return
        call uniform_axis_running(data,shape(x),axis,filter_size,boundary,fill,shift)
        call scatter_complex(data,y)
        if (present(status)) status=ndimage_status_success
    end subroutine uniform1d_complex

    subroutine extrema1d_work(x,filter_size,axis,y,maximum,status,mode,cval,origin)
        real(dp), intent(in) :: x(..) !! source array
        integer, intent(in) :: filter_size !! neighborhood width
        integer, intent(in) :: axis !! one-based axis
        real(dp), intent(out) :: y(..) !! filtered output
        logical, intent(in) :: maximum !! max instead of min
        integer, optional, intent(out) :: status !! status code
        character(len=*), optional, intent(in) :: mode !! boundary mode
        real(dp), optional, intent(in) :: cval !! outside fill
        integer, optional, intent(in) :: origin !! displaced window origin
        real(dp), allocatable :: data(:)
        real(dp) :: fill
        integer :: shift,center
        character(len=:), allocatable :: boundary

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (axis<1 .or. axis>rank(x) .or. filter_size<1) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        if (filter_size>50000 .or. &
            int(size(x,axis),int64)+int(filter_size,int64)-1_int64>int(huge(1),int64)) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        shift=0
        if (present(origin)) shift=origin
        center=filter_size/2
        if (shift<(-center) .or. shift>(filter_size-1)/2) return
        boundary='reflect'
        if (present(mode)) boundary=trim(mode)
        if (.not.mode_valid(boundary)) return
        fill=0.0_dp
        if (present(cval)) fill=cval
        if (.not.ieee_is_finite(fill)) return
        call pack_real(x,data)
        if (.not.all(ieee_is_finite(data))) return
        if (filter_size>1) call extrema_axis_fast(data,shape(x),axis,filter_size,boundary,fill,maximum,shift)
        call scatter_real(data,y)
        if (present(status)) status=ndimage_status_success
    end subroutine extrema1d_work

    subroutine minimum_filter1d_nd(x,filter_size,axis,y,status,mode,cval,origin)
        real(dp), intent(in) :: x(..) !! source real array
        integer, intent(in) :: filter_size !! width of neighborhood
        integer, intent(in) :: axis !! one-based axis
        real(dp), intent(out) :: y(..) !! minima of each window
        integer, optional, intent(out) :: status !! success/invalid/resource guard
        character(len=*), optional, intent(in) :: mode !! boundary extension
        real(dp), optional, intent(in) :: cval !! constant outside value
        integer, optional, intent(in) :: origin !! shifted window center
        call extrema1d_work(x,filter_size,axis,y,.false.,status,mode,cval,origin)
    end subroutine minimum_filter1d_nd

    subroutine maximum_filter1d_nd(x,filter_size,axis,y,status,mode,cval,origin)
        real(dp), intent(in) :: x(..) !! source real array
        integer, intent(in) :: filter_size !! width of neighborhood
        integer, intent(in) :: axis !! one-based axis
        real(dp), intent(out) :: y(..) !! maxima of each window
        integer, optional, intent(out) :: status !! success/invalid/resource guard
        character(len=*), optional, intent(in) :: mode !! boundary extension
        real(dp), optional, intent(in) :: cval !! constant outside value
        integer, optional, intent(in) :: origin !! shifted window center
        call extrema1d_work(x,filter_size,axis,y,.true.,status,mode,cval,origin)
    end subroutine maximum_filter1d_nd

    ! Iterative quickselect: a rank filter needs only one sorted order statistic,
    ! not an entire sorted neighborhood.
    real(dp) function select_order(values,chosen) result(value)
        real(dp), intent(inout) :: values(:) !! scratch neighborhood, rearranged
        integer, intent(in) :: chosen !! one-based order statistic
        integer :: low,high,i,j
        real(dp) :: pivot,tmp

        low=1
        high=size(values)
        do while(low<high)
            pivot=values(low+(high-low)/2)
            i=low
            j=high
            do
                do while(values(i)<pivot)
                    i=i+1
                end do
                do while(values(j)>pivot)
                    j=j-1
                end do
                if (i>j) exit
                tmp=values(i)
                values(i)=values(j)
                values(j)=tmp
                i=i+1
                j=j-1
                if (i>j) exit
            end do
            if (chosen<=j) then
                high=j
            else if (chosen>=i) then
                low=i
            else
                exit
            end if
        end do
        value=values(chosen)
    end function select_order

    subroutine order_stat_process(values,dims,offsets,chosen,boundaries,fill,result)
        real(dp), intent(in) :: values(:) !! flattened real image
        integer, intent(in) :: dims(:) !! image extents
        integer, intent(in) :: offsets(:,:) !! active window offsets, (rank, count)
        integer, intent(in) :: chosen !! zero-based sorted rank
        character(len=*), intent(in) :: boundaries(:) !! per-axis boundary extensions
        real(dp), intent(in) :: fill !! constant exterior value
        real(dp), intent(out) :: result(:) !! packed output
        integer(int64) :: stride(size(dims)),pos,q,rem
        real(dp), allocatable :: samples(:)
        integer :: coords(size(dims)),d,t,k
        logical :: valid

        allocate(samples(size(offsets,2)))
        stride(1)=1_int64
        do d=2,size(dims)
            stride(d)=stride(d-1)*int(dims(d-1),int64)
        end do
        do pos=1_int64,size(values,kind=int64)
            rem=pos-1_int64
            do d=1,size(dims)
                coords(d)=int(modulo(rem,int(dims(d),int64)))
                rem=rem/int(dims(d),int64)
            end do
            do k=1,size(samples)
                valid=.true.
                q=1_int64
                do d=1,size(dims)
                    t=extension_index(coords(d)+offsets(d,k),dims(d),trim(boundaries(d)))
                    if (t==0) then
                        valid=.false.
                        exit
                    end if
                    q=q+int(t-1,int64)*stride(d)
                end do
                if (valid) then
                    samples(k)=values(q)
                else
                    samples(k)=fill
                end if
            end do
            result(pos)=select_order(samples,chosen+1)
        end do
    end subroutine order_stat_process

    subroutine order_stat_work(x,widths,y,rank_value,percentile,status,mode,cval,footprint,origin,modes,axes)
        real(dp), intent(in) :: x(..) !! source image, 1 to 15 dimensions
        integer, intent(in) :: widths(:) !! per-axis window widths
        real(dp), intent(in) :: rank_value !! zero-based order rank, or percentile
        logical, intent(in) :: percentile !! interpret rank_value as percentage
        real(dp), intent(out) :: y(..) !! same-shape image
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary mode
        real(dp), optional, intent(in) :: cval !! fill value
        logical, optional, intent(in) :: footprint(..) !! Boolean custom window, overrides sizes
        integer, optional, intent(in) :: origin(:) !! per-axis window displacement
        integer, optional, intent(in) :: axes(:) !! selected axes, in one-based parameter order
        real(dp), allocatable :: flat(:),result(:)
        character(len=32), allocatable :: boundaries(:)
        character(len=*), optional, intent(in) :: modes(:) !! optional per-axis boundary modes
        logical :: valid_modes
        integer(int64) :: n
        integer, allocatable :: offsets(:,:),compact_offsets(:,:),canonical_axes(:),compact_origin(:)
        real(dp), allocatable :: heights(:)
        integer :: requested,code,k,j,tmp
        logical :: selected(rank(x))
        real(dp) :: fill

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        call selected_axes_mask(rank(x),axes,selected,valid_modes)
        if (.not.valid_modes) return
        if (present(axes)) then
            ! The SciPy rank/median/percentile filters associate compact
            ! parameters with axes in increasing numerical order, even if
            ! the caller lists the selected axes in another order.
            canonical_axes=axes
            do k=1,size(canonical_axes)
                do j=k+1,size(canonical_axes)
                    if (canonical_axes(j)>=canonical_axes(k)) cycle
                    tmp=canonical_axes(j)
                    canonical_axes(j)=canonical_axes(k)
                    canonical_axes(k)=tmp
                end do
            end do
            if (size(widths)/=size(axes)) return
        else
            if (size(widths)/=rank(x)) return
        end if
        if (.not.present(footprint)) then
            if (any(widths<1)) return
        end if
        if (present(axes)) then
            if (size(axes)==0) then
                if (present(footprint)) return
                if (present(origin)) then
                    if (size(origin)/=0) return
                end if
                allocate(offsets(rank(x),1),source=0)
                code=ndimage_status_success
            else
                ! Build a compact selected-axis stencil, then map its offsets
                ! into the full image rank.  This avoids materializing a
                ! rank-15 footprint with singleton unselected dimensions.
                if (present(origin)) then
                    if (size(origin)/=size(axes)) return
                    allocate(compact_origin(size(axes)))
                    ! Footprint/size axes are canonicalized as in SciPy,
                    ! but origin entries refer to the axes as supplied.
                    do k=1,size(axes)
                        do j=1,size(axes)
                            if (canonical_axes(k)==axes(j)) compact_origin(k)=origin(j)
                        end do
                    end do
                end if
                call gray_stencil(size(axes),sizes=widths,footprint=footprint,origin=compact_origin, &
                    offsets=compact_offsets,heights=heights,status=code)
                if (code==ndimage_status_success) then
                    allocate(offsets(rank(x),size(compact_offsets,2)),source=0)
                    do k=1,size(axes)
                        offsets(canonical_axes(k),:)=compact_offsets(k,:)
                    end do
                end if
            end if
        else
            call gray_stencil(rank(x),sizes=widths,footprint=footprint,origin=origin, &
                offsets=offsets,heights=heights,status=code)
        end if
        if (code/=ndimage_status_success) then
            if (present(status)) status=code
            return
        end if
        n=int(size(offsets,2),int64)
        if (.not.ieee_is_finite(rank_value)) return
        if (present(axes)) then
            call selected_filter_modes(rank(x),axes,mode,modes,boundaries,valid_modes)
        else
            call resolve_axis_modes(rank(x),mode,modes,boundaries,valid_modes)
        end if
        if (.not.valid_modes) return
        fill=0.0_dp
        if (present(cval)) fill=cval
        if (.not.ieee_is_finite(fill)) return
        if (percentile) then
            if (rank_value<0.0_dp .or. rank_value>100.0_dp) return
            requested=min(int(n)-1,int(real(n,dp)*rank_value/100.0_dp))
        else
            if (rank_value>=real(n,dp) .or. rank_value< -real(n,dp)) return
            if (rank_value/=real(int(rank_value),dp)) return
            requested=int(rank_value)
            if (requested<0) requested=requested+int(n)
        end if
        call pack_real(x,flat)
        if (.not.all(ieee_is_finite(flat))) return
        allocate(result(size(flat)))
        call order_stat_process(flat,shape(x),offsets,requested,boundaries,fill,result)
        call scatter_real(result,y)
        if (present(status)) status=ndimage_status_success
    end subroutine order_stat_work

    ! Batched counterpart of generic_filter_nd.  A callback receives
    ! windows(neighbor,batch), with active neighbors in C-order and batch
    ! columns in Fortran storage order.  The largest intermediate window
    ! matrix contains at most two million values, independent of image size.
    ! SciPy's 'valid' mode keeps only complete windows, reducing selected
    ! output extents; it is incompatible with nonzero origins/mixed modes.
    subroutine vectorized_layout(dims,output_dims,sizes,footprint,origin,mode,modes,axes, &
        offsets,boundaries,centers,valid,code)
        integer, intent(in) :: dims(:) !! input image dimensions
        integer, intent(in) :: output_dims(:) !! caller's output array dimensions
        integer, optional, intent(in) :: sizes(:) !! compact neighborhood dimensions
        logical, optional, intent(in) :: footprint(..) !! compact mask
        integer, optional, intent(in) :: origin(:) !! compact origin offsets
        character(len=*), optional, intent(in) :: mode !! common boundary or valid
        character(len=*), optional, intent(in) :: modes(:) !! mixed non-valid boundaries
        integer, optional, intent(in) :: axes(:) !! chosen dimensions in stencil order
        integer, allocatable, intent(out) :: offsets(:,:) !! active full-rank offsets
        character(len=32), allocatable, intent(out) :: boundaries(:) !! per-axis boundary
        integer, allocatable, intent(out) :: centers(:) !! shifts for valid output
        logical, intent(out) :: valid !! use reduced output dimensions
        integer, intent(out) :: code !! success, invalid, too expensive
        integer, allocatable :: widths(:),axis_order(:),expected(:)
        real(dp), allocatable :: heights(:)
        logical :: selected(size(dims)),ok
        integer :: nd,i,d

        code=ndimage_status_invalid_input
        nd=size(dims)
        if (nd<1 .or. nd>15 .or. size(output_dims)/=nd) return
        if (any(dims<1) .or. any(output_dims<1)) return
        if (.not.present(sizes) .and. .not.present(footprint)) return
        call selected_axes_mask(nd,axes,selected,ok)
        if (.not.ok) return
        if (present(axes)) then
            axis_order=axes
        else
            axis_order=[(i,i=1,nd)]
        end if
        if (present(footprint)) then
            if (rank(footprint)/=size(axis_order)) return
            widths=shape(footprint)
        else
            if (size(sizes)/=size(axis_order)) return
            widths=sizes
        end if
        if (any(widths<1)) return
        if (present(origin)) then
            if (size(origin)/=size(axis_order)) return
        end if
        valid=.false.
        if (present(mode)) valid=(trim(mode)=='valid')
        allocate(centers(nd),source=0)
        expected=dims
        if (valid) then
            if (present(modes)) return
            if (present(origin)) then
                if (any(origin/=0)) return
            end if
            do i=1,size(axis_order)
                d=axis_order(i)
                expected(d)=dims(d)-widths(i)+1
                centers(d)=widths(i)/2
            end do
            if (any(expected<1) .or. any(expected/=output_dims)) return
            allocate(boundaries(nd))
            boundaries='constant' ! valid windows cannot actually extend past the input
        else
            if (any(output_dims/=dims)) return
            call selected_filter_modes(nd,axes,mode,modes,boundaries,ok)
            if (.not.ok) return
        end if
        ! When footprint is provided SciPy ignores size entirely.
        if (present(footprint)) then
            call selected_gray_stencil(nd,footprint=footprint,origin=origin,axes=axes, &
                offsets=offsets,heights=heights,status=code,c_order=.true.)
        else
            call selected_gray_stencil(nd,sizes=sizes,origin=origin,axes=axes, &
                offsets=offsets,heights=heights,status=code,c_order=.true.)
        end if
    end subroutine vectorized_layout

    ! Determine the ordinary output extents for a multi-result batch. For
    ! rank-15 input, the flat (pixels,channels) representation avoids a
    ! prohibited rank-16 output. Both forms preserve Fortran pixel order.
    subroutine vectorized_multi_shape(dims,ydims,sizes,footprint,mode,axes,odims,channels,code)
        integer, intent(in) :: dims(:) !! shape of source array
        integer, intent(in) :: ydims(:) !! shape of channel-last or flat output
        integer, optional, intent(in) :: sizes(:) !! selected-axis widths
        logical, optional, intent(in) :: footprint(..) !! selected-axis Boolean stencil
        character(len=*), optional, intent(in) :: mode !! extension including valid
        integer, optional, intent(in) :: axes(:) !! one-based axes
        integer, allocatable, intent(out) :: odims(:) !! output pixel extents
        integer, intent(out) :: channels !! results per neighborhood
        integer, intent(out) :: code !! validation status
        integer, allocatable :: axis_order(:),widths(:)
        integer :: i,nd
        integer(int64) :: npixels
        logical :: selected(size(dims)),ok,valid,shaped

        code=ndimage_status_invalid_input
        channels=0
        nd=size(dims)
        if (nd<1 .or. nd>15 .or. any(dims<1) .or. any(ydims<1)) return
        shaped=size(ydims)==nd+1
        if (.not.shaped .and. size(ydims)/=2) return
        channels=ydims(size(ydims))
        if (channels<1) return
        call selected_axes_mask(nd,axes,selected,ok)
        if (.not.ok) return
        if (present(axes)) then
            axis_order=axes
        else
            axis_order=[(i,i=1,nd)]
        end if
        if (present(footprint)) then
            if (rank(footprint)/=size(axis_order)) return
            widths=shape(footprint)
        else
            if (.not.present(sizes)) return
            if (size(sizes)/=size(axis_order)) return
            widths=sizes
        end if
        if (any(widths<1)) return
        odims=dims
        valid=.false.
        if (present(mode)) valid=trim(mode)=='valid'
        if (valid) then
            do i=1,size(axis_order)
                odims(axis_order(i))=dims(axis_order(i))-widths(i)+1
            end do
            if (any(odims<1)) return
        end if
        npixels=1_int64
        do i=1,nd
            npixels=npixels*int(odims(i),int64)
            if (npixels>int(huge(1),int64)) then
                code=ndimage_status_too_expensive
                return
            end if
        end do
        if (shaped) then
            if (any(ydims(:nd)/=odims)) return
        else
            if (int(ydims(1),int64)/=npixels) return
        end if
        code=ndimage_status_success
    end subroutine vectorized_multi_shape

    ! Fast path: if the entire footprint lies inside the source image, all
    ! offsets can be converted to linear displacements once per batch.
    ! Otherwise use the ordinary per-axis extension. This helper is shared
    ! by scalar and multi-result batch filters.
    subroutine gather_vectorized_real(flat,dims,odims,offsets,boundaries,centers,valid,fill,first,windows)
        real(dp), intent(in) :: flat(:) !! packed source in Fortran order
        integer, intent(in) :: dims(:) !! input shape
        integer, intent(in) :: odims(:) !! filtered pixel shape
        integer, intent(in) :: offsets(:,:) !! full-rank active neighbors
        character(len=*), intent(in) :: boundaries(:) !! resolved border modes
        integer, intent(in) :: centers(:) !! valid-mode centering offsets
        logical, intent(in) :: valid !! no extended neighbors in valid mode
        real(dp), intent(in) :: fill !! constant extension value
        integer(int64), intent(in) :: first !! first pixel, one-based flattened index
        real(dp), intent(out) :: windows(:,:) !! extracted windows, one per column
        integer(int64) :: strides(size(dims)),delta(size(offsets,2)),rem,base,q
        integer :: coords(size(dims)),lo(size(dims)),hi(size(dims))
        integer :: d,j,k,t
        logical :: interior,inside

        strides(1)=1_int64
        do d=2,size(dims)
            strides(d)=strides(d-1)*int(dims(d-1),int64)
        end do
        do d=1,size(dims)
            lo(d)=minval(offsets(d,:))
            hi(d)=maxval(offsets(d,:))
        end do
        do k=1,size(offsets,2)
            delta(k)=sum(int(offsets(:,k),int64)*strides)
        end do
        do j=1,size(windows,2)
            rem=first+int(j,int64)-2_int64
            do d=1,size(dims)
                coords(d)=int(modulo(rem,int(odims(d),int64)))+centers(d)
                rem=rem/int(odims(d),int64)
            end do
            interior=all(coords+lo>=0 .and. coords+hi<dims)
            if (interior) then
                base=1_int64+sum(int(coords,int64)*strides)
                do k=1,size(offsets,2)
                    windows(k,j)=flat(base+delta(k))
                end do
            else
                do k=1,size(offsets,2)
                    inside=.true.
                    q=1_int64
                    do d=1,size(dims)
                        if (valid) then
                            t=coords(d)+offsets(d,k)+1
                        else
                            t=extension_index(coords(d)+offsets(d,k),dims(d),trim(boundaries(d)))
                        end if
                        if (t<1 .or. t>dims(d)) then
                            inside=.false.
                            exit
                        end if
                        q=q+int(t-1,int64)*strides(d)
                    end do
                    windows(k,j)=fill
                    if (inside) windows(k,j)=flat(q)
                end do
            end if
        end do
    end subroutine gather_vectorized_real

    subroutine gather_vectorized_complex(flat,dims,odims,offsets,boundaries,centers,valid,fill,first,windows)
        complex(dp), intent(in) :: flat(:) !! packed complex source
        integer, intent(in) :: dims(:) !! source shape
        integer, intent(in) :: odims(:) !! output pixel shape
        integer, intent(in) :: offsets(:,:) !! active full-rank neighbors
        character(len=*), intent(in) :: boundaries(:) !! axis boundary modes
        integer, intent(in) :: centers(:) !! valid-mode center offsets
        logical, intent(in) :: valid !! unextended valid region
        complex(dp), intent(in) :: fill !! exterior value
        integer(int64), intent(in) :: first !! first output pixel
        complex(dp), intent(out) :: windows(:,:) !! gathered windows
        integer(int64) :: strides(size(dims)),delta(size(offsets,2)),rem,base,q
        integer :: coords(size(dims)),lo(size(dims)),hi(size(dims))
        integer :: d,j,k,t
        logical :: interior,inside

        strides(1)=1_int64
        do d=2,size(dims)
            strides(d)=strides(d-1)*int(dims(d-1),int64)
        end do
        do d=1,size(dims)
            lo(d)=minval(offsets(d,:))
            hi(d)=maxval(offsets(d,:))
        end do
        do k=1,size(offsets,2)
            delta(k)=sum(int(offsets(:,k),int64)*strides)
        end do
        do j=1,size(windows,2)
            rem=first+int(j,int64)-2_int64
            do d=1,size(dims)
                coords(d)=int(modulo(rem,int(odims(d),int64)))+centers(d)
                rem=rem/int(odims(d),int64)
            end do
            interior=all(coords+lo>=0 .and. coords+hi<dims)
            if (interior) then
                base=1_int64+sum(int(coords,int64)*strides)
                do k=1,size(offsets,2)
                    windows(k,j)=flat(base+delta(k))
                end do
            else
                do k=1,size(offsets,2)
                    inside=.true.
                    q=1_int64
                    do d=1,size(dims)
                        if (valid) then
                            t=coords(d)+offsets(d,k)+1
                        else
                            t=extension_index(coords(d)+offsets(d,k),dims(d),trim(boundaries(d)))
                        end if
                        if (t<1 .or. t>dims(d)) then
                            inside=.false.
                            exit
                        end if
                        q=q+int(t-1,int64)*strides(d)
                    end do
                    windows(k,j)=fill
                    if (inside) windows(k,j)=flat(q)
                end do
            end if
        end do
    end subroutine gather_vectorized_complex

    subroutine vectorized_filter_real_nd(x,func,y,status,sizes,footprint,origin,mode,cval,modes,axes,batch_size)
        real(dp), intent(in) :: x(..) !! source array of ranks 1 through 15
        procedure(ndimage_batch_real_func) :: func !! batched neighborhood callback
        real(dp), intent(out) :: y(..) !! filtered data; smaller shape in valid mode
        integer, optional, intent(out) :: status !! 0 success, 1 invalid, 2 too large
        integer, optional, intent(in) :: sizes(:) !! footprint widths in axes order
        logical, optional, intent(in) :: footprint(..) !! active stencil cells
        integer, optional, intent(in) :: origin(:) !! origins in axes order
        character(len=*), optional, intent(in) :: mode !! boundary extension or valid
        real(dp), optional, intent(in) :: cval !! exterior value, may be NaN
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundaries
        integer, optional, intent(in) :: axes(:) !! selected axes, one-based
        integer, optional, intent(in) :: batch_size !! upper bound on output pixels per call
        integer, allocatable :: offsets(:,:),centers(:)
        character(len=32), allocatable :: boundaries(:)
        real(dp), allocatable :: flat(:),result(:),windows(:,:),batch_results(:)
        integer(int64) :: pos,nout
        integer :: dims(rank(x)),odims(rank(x)),batch,code,n,nw,cap
        logical :: valid
        real(dp) :: fill

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x)) return
        dims=shape(x)
        odims=shape(y)
        if (size(x,kind=int64)>int(huge(1),int64) .or. &
            size(y,kind=int64)>int(huge(1),int64)) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        call vectorized_layout(dims,odims,sizes,footprint,origin,mode,modes,axes, &
            offsets,boundaries,centers,valid,code)
        if (code/=ndimage_status_success) then
            if (present(status)) status=code
            return
        end if
        nw=size(offsets,2)
        cap=max(1,2000000/nw)
        batch=min(4096,cap)
        if (present(batch_size)) then
            if (batch_size<1) return
            batch=min(batch_size,cap)
        end if
        fill=0.0_dp
        if (present(cval)) fill=cval
        call pack_real(x,flat)
        nout=size(y,kind=int64)
        allocate(result(int(nout)),windows(nw,batch),batch_results(batch))
        do pos=1_int64,nout,int(batch,int64)
            n=int(min(int(batch,int64),nout-pos+1_int64))
            call gather_vectorized_real(flat,dims,odims,offsets,boundaries,centers, &
                valid,fill,pos,windows(:,:n))
            call func(windows(:,:n),batch_results(:n))
            result(pos:pos+int(n,int64)-1_int64)=batch_results(:n)
        end do
        call scatter_real(result,y)
        if (present(status)) status=ndimage_status_success
    end subroutine vectorized_filter_real_nd

    subroutine vectorized_filter_complex_nd(x,func,y,status,sizes,footprint,origin,mode,cval,modes,axes,batch_size)
        complex(dp), intent(in) :: x(..) !! complex source array, ranks 1 through 15
        procedure(ndimage_batch_complex_func) :: func !! batched complex reduction
        complex(dp), intent(out) :: y(..) !! complex output, reduced for valid mode
        integer, optional, intent(out) :: status !! 0 success, 1 invalid, 2 too large
        integer, optional, intent(in) :: sizes(:) !! rectangular neighborhood widths
        logical, optional, intent(in) :: footprint(..) !! active stencil cells
        integer, optional, intent(in) :: origin(:) !! per-axis origins
        character(len=*), optional, intent(in) :: mode !! boundary extension or valid
        complex(dp), optional, intent(in) :: cval !! complex constant exterior
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundaries
        integer, optional, intent(in) :: axes(:) !! chosen image dimensions
        integer, optional, intent(in) :: batch_size !! maximum output pixels per batch
        integer, allocatable :: offsets(:,:),centers(:)
        character(len=32), allocatable :: boundaries(:)
        complex(dp), allocatable :: flat(:),result(:),windows(:,:),batch_results(:)
        integer(int64) :: pos,nout
        integer :: dims(rank(x)),odims(rank(x)),batch,code,n,nw,cap
        logical :: valid
        complex(dp) :: fill

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x)) return
        dims=shape(x)
        odims=shape(y)
        if (size(x,kind=int64)>int(huge(1),int64) .or. &
            size(y,kind=int64)>int(huge(1),int64)) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        call vectorized_layout(dims,odims,sizes,footprint,origin,mode,modes,axes, &
            offsets,boundaries,centers,valid,code)
        if (code/=ndimage_status_success) then
            if (present(status)) status=code
            return
        end if
        nw=size(offsets,2)
        cap=max(1,2000000/nw)
        batch=min(4096,cap)
        if (present(batch_size)) then
            if (batch_size<1) return
            batch=min(batch_size,cap)
        end if
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cval
        call pack_complex(x,flat)
        nout=size(y,kind=int64)
        allocate(result(int(nout)),windows(nw,batch),batch_results(batch))
        do pos=1_int64,nout,int(batch,int64)
            n=int(min(int(batch,int64),nout-pos+1_int64))
            call gather_vectorized_complex(flat,dims,odims,offsets,boundaries,centers, &
                valid,fill,pos,windows(:,:n))
            call func(windows(:,:n),batch_results(:n))
            result(pos:pos+int(n,int64)-1_int64)=batch_results(:n)
        end do
        call scatter_complex(result,y)
        if (present(status)) status=ndimage_status_success
    end subroutine vectorized_filter_complex_nd

    ! Multi-valued batched filtering. The callback emits results(channel,pixel)
    ! and the result array stores channels as its final axis. At rank 15,
    ! y(pixels,channels) is an alternative that avoids rank-16 arrays.
    subroutine vectorized_filter_multi_real_nd(x,func,y,status,sizes,footprint,origin,mode,cval,modes,axes,batch_size)
        real(dp), intent(in) :: x(..) !! real image, rank 1 through 15
        procedure(ndimage_batch_real_multi_func) :: func !! multi-result neighborhood callback
        real(dp), intent(out) :: y(..) !! shape [output shape,channels] or [pixels,channels]
        integer, optional, intent(out) :: status !! 0 success, 1 invalid, 2 resource limit
        integer, optional, intent(in) :: sizes(:) !! compact window widths
        logical, optional, intent(in) :: footprint(..) !! compact Boolean stencil
        integer, optional, intent(in) :: origin(:) !! selected-axis origins
        character(len=*), optional, intent(in) :: mode !! border mode, including valid
        real(dp), optional, intent(in) :: cval !! constant exterior value, may be NaN
        character(len=*), optional, intent(in) :: modes(:) !! mixed axis boundary modes
        integer, optional, intent(in) :: axes(:) !! selected one-based axes
        integer, optional, intent(in) :: batch_size !! maximum windows per callback
        integer, allocatable :: offsets(:,:),centers(:),odims(:)
        character(len=32), allocatable :: boundaries(:)
        real(dp), allocatable :: flat(:),result(:),windows(:,:),batch_results(:,:)
        integer(int64) :: pos,nout
        integer :: dims(rank(x)),channels,code,nw,batch,cap,n,c
        logical :: valid
        real(dp) :: fill

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)<2 .or. rank(y)>15) return
        if (size(x,kind=int64)>int(huge(1),int64) .or. &
            size(y,kind=int64)>int(huge(1),int64)) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        dims=shape(x)
        call vectorized_multi_shape(dims,shape(y),sizes,footprint,mode,axes,odims,channels,code)
        if (code/=ndimage_status_success) then
            if (present(status)) status=code
            return
        end if
        call vectorized_layout(dims,odims,sizes,footprint,origin,mode,modes,axes, &
            offsets,boundaries,centers,valid,code)
        if (code/=ndimage_status_success) then
            if (present(status)) status=code
            return
        end if
        nw=size(offsets,2)
        if (int(nw,int64)+int(channels,int64)>2000000_int64) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        cap=2000000/(nw+channels)
        batch=min(4096,cap)
        if (present(batch_size)) then
            if (batch_size<1) return
            batch=min(batch_size,cap)
        end if
        fill=0.0_dp
        if (present(cval)) fill=cval
        call pack_real(x,flat)
        nout=product(int(odims,int64))
        allocate(result(int(nout*int(channels,int64))),windows(nw,batch),batch_results(channels,batch))
        do pos=1_int64,nout,int(batch,int64)
            n=int(min(int(batch,int64),nout-pos+1_int64))
            call gather_vectorized_real(flat,dims,odims,offsets,boundaries,centers, &
                valid,fill,pos,windows(:,:n))
            call func(windows(:,:n),batch_results(:,:n))
            do c=1,channels
                result((int(c,int64)-1_int64)*nout+pos:(int(c,int64)-1_int64)*nout+pos+n-1)=batch_results(c,:n)
            end do
        end do
        call scatter_real(result,y)
        if (present(status)) status=ndimage_status_success
    end subroutine vectorized_filter_multi_real_nd

    subroutine vectorized_filter_multi_complex_nd(x,func,y,status,sizes,footprint,origin,mode,cval,modes,axes,batch_size)
        complex(dp), intent(in) :: x(..) !! complex source, rank 1 through 15
        procedure(ndimage_batch_complex_multi_func) :: func !! multi-result complex callback
        complex(dp), intent(out) :: y(..) !! output with trailing channel axis, or flat
        integer, optional, intent(out) :: status !! 0 success, 1 invalid, 2 resource limit
        integer, optional, intent(in) :: sizes(:) !! selected-axis window widths
        logical, optional, intent(in) :: footprint(..) !! compact active footprint
        integer, optional, intent(in) :: origin(:) !! shifted window origins
        character(len=*), optional, intent(in) :: mode !! border extension mode or valid
        complex(dp), optional, intent(in) :: cval !! exterior value
        character(len=*), optional, intent(in) :: modes(:) !! per-axis boundaries
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        integer, optional, intent(in) :: batch_size !! maximum windows per callback
        integer, allocatable :: offsets(:,:),centers(:),odims(:)
        character(len=32), allocatable :: boundaries(:)
        complex(dp), allocatable :: flat(:),result(:),windows(:,:),batch_results(:,:)
        integer(int64) :: pos,nout
        integer :: dims(rank(x)),channels,code,nw,batch,cap,n,c
        logical :: valid
        complex(dp) :: fill

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)<2 .or. rank(y)>15) return
        if (size(x,kind=int64)>int(huge(1),int64) .or. &
            size(y,kind=int64)>int(huge(1),int64)) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        dims=shape(x)
        call vectorized_multi_shape(dims,shape(y),sizes,footprint,mode,axes,odims,channels,code)
        if (code/=ndimage_status_success) then
            if (present(status)) status=code
            return
        end if
        call vectorized_layout(dims,odims,sizes,footprint,origin,mode,modes,axes, &
            offsets,boundaries,centers,valid,code)
        if (code/=ndimage_status_success) then
            if (present(status)) status=code
            return
        end if
        nw=size(offsets,2)
        if (int(nw,int64)+int(channels,int64)>2000000_int64) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        cap=2000000/(nw+channels)
        batch=min(4096,cap)
        if (present(batch_size)) then
            if (batch_size<1) return
            batch=min(batch_size,cap)
        end if
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cval
        call pack_complex(x,flat)
        nout=product(int(odims,int64))
        allocate(result(int(nout*int(channels,int64))),windows(nw,batch),batch_results(channels,batch))
        do pos=1_int64,nout,int(batch,int64)
            n=int(min(int(batch,int64),nout-pos+1_int64))
            call gather_vectorized_complex(flat,dims,odims,offsets,boundaries,centers, &
                valid,fill,pos,windows(:,:n))
            call func(windows(:,:n),batch_results(:,:n))
            do c=1,channels
                result((int(c,int64)-1_int64)*nout+pos:(int(c,int64)-1_int64)*nout+pos+n-1)=batch_results(c,:n)
            end do
        end do
        call scatter_complex(result,y)
        if (present(status)) status=ndimage_status_success
    end subroutine vectorized_filter_multi_complex_nd

    ! Mixed-type vectorized filtering uses the same footprint extraction
    ! and batch limits as the homogeneous real/complex overloads.
    subroutine vectorized_filter_r2c_nd(x,func,y,status,sizes,footprint,origin,mode,cval,modes,axes,batch_size)
        real(dp), intent(in) :: x(..) !! input image of ranks 1 through 15
        procedure(ndimage_batch_real_to_complex_func) :: func !! user-supplied mixed-type callback
        complex(dp), intent(out) :: y(..) !! image-shaped output
        integer, optional, intent(out) :: status !! zero for success
        integer, optional, intent(in) :: sizes(:) !! selected-axis widths
        logical, optional, intent(in) :: footprint(..) !! selected active neighbors
        integer, optional, intent(in) :: origin(:) !! selected-axis offsets
        character(len=*), optional, intent(in) :: mode !! boundary mode or valid
        real(dp), optional, intent(in) :: cval !! outside value, same type as input
        character(len=*), optional, intent(in) :: modes(:) !! axis-specific modes
        integer, optional, intent(in) :: axes(:) !! selected axes
        integer, optional, intent(in) :: batch_size !! maximum batch width
        integer, allocatable :: offsets(:,:),centers(:)
        character(len=32), allocatable :: boundaries(:)
        real(dp), allocatable :: flat(:),windows(:,:)
        complex(dp), allocatable :: result(:),batch_results(:)
        integer(int64) :: pos,nout
        integer :: dims(rank(x)),odims(rank(x)),batch,code,n,nw,cap
        logical :: valid
        real(dp) :: fill

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x)) return
        dims=shape(x)
        odims=shape(y)
        if (size(x,kind=int64)>int(huge(1),int64) .or. &
            size(y,kind=int64)>int(huge(1),int64)) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        call vectorized_layout(dims,odims,sizes,footprint,origin,mode,modes,axes, &
            offsets,boundaries,centers,valid,code)
        if (code/=ndimage_status_success) then
            if (present(status)) status=code
            return
        end if
        nw=size(offsets,2)
        cap=max(1,2000000/nw)
        batch=min(4096,cap)
        if (present(batch_size)) then
            if (batch_size<1) return
            batch=min(batch_size,cap)
        end if
        fill=0.0_dp
        if (present(cval)) fill=cval
        call pack_real(x,flat)
        nout=size(y,kind=int64)
        allocate(result(int(nout)),windows(nw,batch),batch_results(batch))
        do pos=1_int64,nout,int(batch,int64)
            n=int(min(int(batch,int64),nout-pos+1_int64))
            call gather_vectorized_real(flat,dims,odims,offsets,boundaries,centers, &
                valid,fill,pos,windows(:,:n))
            call func(windows(:,:n),batch_results(:n))
            result(pos:pos+int(n,int64)-1_int64)=batch_results(:n)
        end do
        call scatter_complex(result,y)
        if (present(status)) status=ndimage_status_success
    end subroutine vectorized_filter_r2c_nd

    subroutine vectorized_filter_multi_r2c_nd(x,func,y,status,sizes,footprint,origin,mode,cval,modes,axes,batch_size)
        real(dp), intent(in) :: x(..) !! input image of ranks 1 through 15
        procedure(ndimage_batch_real_to_complex_multi_func) :: func !! user-supplied mixed-type multi-result callback
        complex(dp), intent(out) :: y(..) !! last-axis channels or rank-15 flat form
        integer, optional, intent(out) :: status !! zero for success
        integer, optional, intent(in) :: sizes(:) !! selected-axis widths
        logical, optional, intent(in) :: footprint(..) !! active compact stencil
        integer, optional, intent(in) :: origin(:) !! selected-axis origins
        character(len=*), optional, intent(in) :: mode !! boundary mode or valid
        real(dp), optional, intent(in) :: cval !! outside input value
        character(len=*), optional, intent(in) :: modes(:) !! selected axis boundaries
        integer, optional, intent(in) :: axes(:) !! selected axes
        integer, optional, intent(in) :: batch_size !! maximum batch width
        integer, allocatable :: offsets(:,:),centers(:),odims(:)
        character(len=32), allocatable :: boundaries(:)
        real(dp), allocatable :: flat(:),windows(:,:)
        complex(dp), allocatable :: result(:),batch_results(:,:)
        integer(int64) :: pos,nout
        integer :: dims(rank(x)),channels,code,nw,batch,cap,n,c
        logical :: valid
        real(dp) :: fill

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)<2 .or. rank(y)>15) return
        if (size(x,kind=int64)>int(huge(1),int64) .or. &
            size(y,kind=int64)>int(huge(1),int64)) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        dims=shape(x)
        call vectorized_multi_shape(dims,shape(y),sizes,footprint,mode,axes,odims,channels,code)
        if (code/=ndimage_status_success) then
            if (present(status)) status=code
            return
        end if
        call vectorized_layout(dims,odims,sizes,footprint,origin,mode,modes,axes, &
            offsets,boundaries,centers,valid,code)
        if (code/=ndimage_status_success) then
            if (present(status)) status=code
            return
        end if
        nw=size(offsets,2)
        if (int(nw,int64)+int(channels,int64)>2000000_int64) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        cap=2000000/(nw+channels)
        batch=min(4096,cap)
        if (present(batch_size)) then
            if (batch_size<1) return
            batch=min(batch_size,cap)
        end if
        fill=0.0_dp
        if (present(cval)) fill=cval
        call pack_real(x,flat)
        nout=product(int(odims,int64))
        allocate(result(int(nout*int(channels,int64))),windows(nw,batch),batch_results(channels,batch))
        do pos=1_int64,nout,int(batch,int64)
            n=int(min(int(batch,int64),nout-pos+1_int64))
            call gather_vectorized_real(flat,dims,odims,offsets,boundaries,centers, &
                valid,fill,pos,windows(:,:n))
            call func(windows(:,:n),batch_results(:,:n))
            do c=1,channels
                result((int(c,int64)-1_int64)*nout+pos: &
                    (int(c,int64)-1_int64)*nout+pos+n-1)=batch_results(c,:n)
            end do
        end do
        call scatter_complex(result,y)
        if (present(status)) status=ndimage_status_success
    end subroutine vectorized_filter_multi_r2c_nd

    ! Mixed-type vectorized filtering uses the same footprint extraction
    ! and batch limits as the homogeneous real/complex overloads.
    subroutine vectorized_filter_c2r_nd(x,func,y,status,sizes,footprint,origin,mode,cval,modes,axes,batch_size)
        complex(dp), intent(in) :: x(..) !! input image of ranks 1 through 15
        procedure(ndimage_batch_complex_to_real_func) :: func !! user-supplied mixed-type callback
        real(dp), intent(out) :: y(..) !! image-shaped output
        integer, optional, intent(out) :: status !! zero for success
        integer, optional, intent(in) :: sizes(:) !! selected-axis widths
        logical, optional, intent(in) :: footprint(..) !! selected active neighbors
        integer, optional, intent(in) :: origin(:) !! selected-axis offsets
        character(len=*), optional, intent(in) :: mode !! boundary mode or valid
        complex(dp), optional, intent(in) :: cval !! outside value, same type as input
        character(len=*), optional, intent(in) :: modes(:) !! axis-specific modes
        integer, optional, intent(in) :: axes(:) !! selected axes
        integer, optional, intent(in) :: batch_size !! maximum batch width
        integer, allocatable :: offsets(:,:),centers(:)
        character(len=32), allocatable :: boundaries(:)
        complex(dp), allocatable :: flat(:),windows(:,:)
        real(dp), allocatable :: result(:),batch_results(:)
        integer(int64) :: pos,nout
        integer :: dims(rank(x)),odims(rank(x)),batch,code,n,nw,cap
        logical :: valid
        complex(dp) :: fill

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x)) return
        dims=shape(x)
        odims=shape(y)
        if (size(x,kind=int64)>int(huge(1),int64) .or. &
            size(y,kind=int64)>int(huge(1),int64)) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        call vectorized_layout(dims,odims,sizes,footprint,origin,mode,modes,axes, &
            offsets,boundaries,centers,valid,code)
        if (code/=ndimage_status_success) then
            if (present(status)) status=code
            return
        end if
        nw=size(offsets,2)
        cap=max(1,2000000/nw)
        batch=min(4096,cap)
        if (present(batch_size)) then
            if (batch_size<1) return
            batch=min(batch_size,cap)
        end if
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cval
        call pack_complex(x,flat)
        nout=size(y,kind=int64)
        allocate(result(int(nout)),windows(nw,batch),batch_results(batch))
        do pos=1_int64,nout,int(batch,int64)
            n=int(min(int(batch,int64),nout-pos+1_int64))
            call gather_vectorized_complex(flat,dims,odims,offsets,boundaries,centers, &
                valid,fill,pos,windows(:,:n))
            call func(windows(:,:n),batch_results(:n))
            result(pos:pos+int(n,int64)-1_int64)=batch_results(:n)
        end do
        call scatter_real(result,y)
        if (present(status)) status=ndimage_status_success
    end subroutine vectorized_filter_c2r_nd

    subroutine vectorized_filter_multi_c2r_nd(x,func,y,status,sizes,footprint,origin,mode,cval,modes,axes,batch_size)
        complex(dp), intent(in) :: x(..) !! input image of ranks 1 through 15
        procedure(ndimage_batch_complex_to_real_multi_func) :: func !! user-supplied mixed-type multi-result callback
        real(dp), intent(out) :: y(..) !! last-axis channels or rank-15 flat form
        integer, optional, intent(out) :: status !! zero for success
        integer, optional, intent(in) :: sizes(:) !! selected-axis widths
        logical, optional, intent(in) :: footprint(..) !! active compact stencil
        integer, optional, intent(in) :: origin(:) !! selected-axis origins
        character(len=*), optional, intent(in) :: mode !! boundary mode or valid
        complex(dp), optional, intent(in) :: cval !! outside input value
        character(len=*), optional, intent(in) :: modes(:) !! selected axis boundaries
        integer, optional, intent(in) :: axes(:) !! selected axes
        integer, optional, intent(in) :: batch_size !! maximum batch width
        integer, allocatable :: offsets(:,:),centers(:),odims(:)
        character(len=32), allocatable :: boundaries(:)
        complex(dp), allocatable :: flat(:),windows(:,:)
        real(dp), allocatable :: result(:),batch_results(:,:)
        integer(int64) :: pos,nout
        integer :: dims(rank(x)),channels,code,nw,batch,cap,n,c
        logical :: valid
        complex(dp) :: fill

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)<2 .or. rank(y)>15) return
        if (size(x,kind=int64)>int(huge(1),int64) .or. &
            size(y,kind=int64)>int(huge(1),int64)) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        dims=shape(x)
        call vectorized_multi_shape(dims,shape(y),sizes,footprint,mode,axes,odims,channels,code)
        if (code/=ndimage_status_success) then
            if (present(status)) status=code
            return
        end if
        call vectorized_layout(dims,odims,sizes,footprint,origin,mode,modes,axes, &
            offsets,boundaries,centers,valid,code)
        if (code/=ndimage_status_success) then
            if (present(status)) status=code
            return
        end if
        nw=size(offsets,2)
        if (int(nw,int64)+int(channels,int64)>2000000_int64) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        cap=2000000/(nw+channels)
        batch=min(4096,cap)
        if (present(batch_size)) then
            if (batch_size<1) return
            batch=min(batch_size,cap)
        end if
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cval
        call pack_complex(x,flat)
        nout=product(int(odims,int64))
        allocate(result(int(nout*int(channels,int64))),windows(nw,batch),batch_results(channels,batch))
        do pos=1_int64,nout,int(batch,int64)
            n=int(min(int(batch,int64),nout-pos+1_int64))
            call gather_vectorized_complex(flat,dims,odims,offsets,boundaries,centers, &
                valid,fill,pos,windows(:,:n))
            call func(windows(:,:n),batch_results(:,:n))
            do c=1,channels
                result((int(c,int64)-1_int64)*nout+pos: &
                    (int(c,int64)-1_int64)*nout+pos+n-1)=batch_results(c,:n)
            end do
        end do
        call scatter_real(result,y)
        if (present(status)) status=ndimage_status_success
    end subroutine vectorized_filter_multi_c2r_nd

    ! Apply a user-supplied real reduction to each active neighborhood.
    ! The callback receives active footprint elements in SciPy's C-order.
    subroutine generic_filter_real_nd(x,func,y,status,sizes,footprint,origin,mode,cval,modes,axes)
        real(dp), intent(in) :: x(..) !! real input image, ranks 1 through 15
        procedure(ndimage_neighborhood_func) :: func !! neighborhood reduction callback
        real(dp), intent(out) :: y(..) !! filtered output with shape(x)
        integer, optional, intent(out) :: status !! 0 success, 1 invalid, 2 oversized stencil
        integer, optional, intent(in) :: sizes(:) !! rectangular footprint sizes
        logical, optional, intent(in) :: footprint(..) !! custom active neighborhood
        integer, optional, intent(in) :: origin(:) !! offset of each axis window center
        character(len=*), optional, intent(in) :: mode !! border extension mode
        real(dp), optional, intent(in) :: cval !! constant extension value
        character(len=*), optional, intent(in) :: modes(:) !! optional per-axis modes
        integer, optional, intent(in) :: axes(:) !! one-based selected dimensions
        integer, allocatable :: offsets(:,:)
        real(dp), allocatable :: heights(:),flat(:),result(:),samples(:)
        character(len=32), allocatable :: boundaries(:)
        integer(int64) :: strides(rank(x)),pos,rem,q
        integer :: code,coords(rank(x)),d,k,t
        logical :: valid,valid_modes
        real(dp) :: fill

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(x)/=rank(y)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        if (.not.present(sizes) .and. .not.present(footprint)) return
        call selected_filter_modes(rank(x),axes,mode,modes,boundaries,valid_modes)
        if (.not.valid_modes) return
        fill=0.0_dp
        if (present(cval)) fill=cval
        if (.not.ieee_is_finite(fill)) return
        call selected_gray_stencil(rank(x),sizes=sizes,footprint=footprint,origin=origin, &
            axes=axes,offsets=offsets,heights=heights,status=code,c_order=.true.,canonical=.true.)
        if (code/=ndimage_status_success) then
            if (present(status)) status=code
            return
        end if
        call pack_real(x,flat)
        if (.not.all(ieee_is_finite(flat))) return
        allocate(result(size(flat)),samples(size(offsets,2)))
        strides(1)=1_int64
        do d=2,rank(x)
            strides(d)=strides(d-1)*int(size(x,d-1),int64)
        end do
        do pos=1_int64,size(flat,kind=int64)
            rem=pos-1_int64
            do d=1,rank(x)
                coords(d)=int(modulo(rem,int(size(x,d),int64)))
                rem=rem/int(size(x,d),int64)
            end do
            do k=1,size(samples)
                valid=.true.
                q=1_int64
                do d=1,rank(x)
                    t=extension_index(coords(d)+offsets(d,k),size(x,d),trim(boundaries(d)))
                    if (t==0) then
                        valid=.false.
                        exit
                    end if
                    q=q+int(t-1,int64)*strides(d)
                end do
                samples(k)=fill
                if (valid) samples(k)=flat(q)
            end do
            result(pos)=func(samples)
        end do
        call scatter_real(result,y)
        if (present(status)) status=ndimage_status_success
    end subroutine generic_filter_real_nd

    subroutine generic_filter_complex_nd(x,func,y,status,sizes,footprint,origin,mode,cval,modes,axes)
        complex(dp), intent(in) :: x(..) !! complex input image, ranks 1 through 15
        procedure(ndimage_complex_neighborhood_func) :: func !! complex neighborhood callback
        complex(dp), intent(out) :: y(..) !! complex output with shape(x)
        integer, optional, intent(out) :: status !! 0 success, 1 invalid, 2 oversized stencil
        integer, optional, intent(in) :: sizes(:) !! rectangular footprint sizes
        logical, optional, intent(in) :: footprint(..) !! custom active neighborhood
        integer, optional, intent(in) :: origin(:) !! offset of each axis window center
        character(len=*), optional, intent(in) :: mode !! border extension mode
        complex(dp), optional, intent(in) :: cval !! constant extension value
        character(len=*), optional, intent(in) :: modes(:) !! optional per-axis modes
        integer, optional, intent(in) :: axes(:) !! one-based selected dimensions
        integer, allocatable :: offsets(:,:)
        real(dp), allocatable :: heights(:)
        complex(dp), allocatable :: flat(:),result(:),samples(:)
        character(len=32), allocatable :: boundaries(:)
        integer(int64) :: strides(rank(x)),pos,rem,q
        integer :: code,coords(rank(x)),d,k,t
        logical :: valid,valid_modes
        complex(dp) :: fill

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(x)/=rank(y)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        if (.not.present(sizes) .and. .not.present(footprint)) return
        call selected_filter_modes(rank(x),axes,mode,modes,boundaries,valid_modes)
        if (.not.valid_modes) return
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cval
        if (.not.ieee_is_finite(real(fill,dp)) .or. .not.ieee_is_finite(aimag(fill))) return
        call selected_gray_stencil(rank(x),sizes=sizes,footprint=footprint,origin=origin, &
            axes=axes,offsets=offsets,heights=heights,status=code,c_order=.true.,canonical=.true.)
        if (code/=ndimage_status_success) then
            if (present(status)) status=code
            return
        end if
        call pack_complex(x,flat)
        if (.not.all(ieee_is_finite(real(flat,dp))) .or. &
            .not.all(ieee_is_finite(aimag(flat)))) return
        allocate(result(size(flat)),samples(size(offsets,2)))
        strides(1)=1_int64
        do d=2,rank(x)
            strides(d)=strides(d-1)*int(size(x,d-1),int64)
        end do
        do pos=1_int64,size(flat,kind=int64)
            rem=pos-1_int64
            do d=1,rank(x)
                coords(d)=int(modulo(rem,int(size(x,d),int64)))
                rem=rem/int(size(x,d),int64)
            end do
            do k=1,size(samples)
                valid=.true.
                q=1_int64
                do d=1,rank(x)
                    t=extension_index(coords(d)+offsets(d,k),size(x,d),trim(boundaries(d)))
                    if (t==0) then
                        valid=.false.
                        exit
                    end if
                    q=q+int(t-1,int64)*strides(d)
                end do
                samples(k)=fill
                if (valid) samples(k)=flat(q)
            end do
            result(pos)=func(samples)
        end do
        call scatter_complex(result,y)
        if (present(status)) status=ndimage_status_success
    end subroutine generic_filter_complex_nd

    ! SciPy-style line callback: the caller receives a padded one-dimensional
    ! input of extent axis_length+filter_size-1 and fills axis_length outputs.
    ! Origin shifts the window center: input_line(j:j+filter_size-1) spans
    ! input(j-filter_size/2-origin : j+filter_size/2-origin) in zero-based
    ! image coordinates, allowing the callback to reuse several offsets.
    subroutine generic_filter1d_real_nd(x,func,filter_size,axis,y,status,mode,cval,origin)
        real(dp), intent(in) :: x(..) !! real image with ranks 1 through 15
        procedure(ndimage_line_func) :: func !! one call for each padded axis line
        integer, intent(in) :: filter_size !! positive filter width
        integer, intent(in) :: axis !! one-based line axis
        real(dp), intent(out) :: y(..) !! image of same shape as x
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 resource guard
        character(len=*), optional, intent(in) :: mode !! axis boundary extension
        real(dp), optional, intent(in) :: cval !! constant exterior value
        integer, optional, intent(in) :: origin !! center shift
        real(dp), allocatable :: flat(:),result(:),input_line(:),output_line(:)
        character(len=:), allocatable :: boundary
        integer :: dims(rank(x)),nd,n,shift,i,j,t
        integer(int64) :: stride,block,outer,o,inner,base,pos,total
        real(dp) :: fill

        if (present(status)) status=ndimage_status_invalid_input
        nd=rank(x)
        if (nd<1 .or. nd>15 .or. rank(y)/=nd) return
        dims=shape(x)
        if (any(dims<1) .or. any(dims/=shape(y))) return
        if (axis<1 .or. axis>nd .or. filter_size<1) return
        shift=0
        if (present(origin)) shift=origin
        if (shift<-(filter_size/2) .or. shift>(filter_size-1)/2) return
        if (filter_size>max_kernel_radius) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        if (dims(axis)>huge(1)-filter_size+1) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        total=size(x,kind=int64)
        if (total>int(huge(1),int64)) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        boundary='reflect'
        if (present(mode)) boundary=mode
        if (.not.mode_valid(boundary)) return
        fill=0.0_dp
        if (present(cval)) fill=cval
        if (.not.ieee_is_finite(fill)) return
        call pack_real(x,flat)
        if (.not.all(ieee_is_finite(flat))) return
        allocate(result(size(flat)))
        n=dims(axis)
        allocate(input_line(n+filter_size-1),output_line(n))
        stride=1_int64
        do i=1,axis-1
            stride=stride*int(dims(i),int64)
        end do
        block=stride*int(n,int64)
        outer=total/block
        do o=0_int64,outer-1_int64
            do inner=0_int64,stride-1_int64
                base=o*block+inner
                do j=1,size(input_line)
                    t=extension_index(j-1-filter_size/2-shift,n,boundary)
                    input_line(j)=fill
                    if (t>0) input_line(j)=flat(base+int(t-1,int64)*stride+1_int64)
                end do
                call func(input_line,output_line)
                do j=1,n
                    pos=base+int(j-1,int64)*stride+1_int64
                    result(pos)=output_line(j)
                end do
            end do
        end do
        call scatter_real(result,y)
        if (present(status)) status=ndimage_status_success
    end subroutine generic_filter1d_real_nd
    ! SciPy's generic_filter1d currently accepts only real inputs; this
    ! additional overload extends the same padded-line contract to complex
    ! arrays without changing the existing real callback interface.
    subroutine generic_filter1d_complex_nd(x,func,filter_size,axis,y,status,mode,cval,origin)
        complex(dp), intent(in) :: x(..) !! real image with ranks 1 through 15
        procedure(ndimage_complex_line_func) :: func !! one call for each padded axis line
        integer, intent(in) :: filter_size !! positive filter width
        integer, intent(in) :: axis !! one-based line axis
        complex(dp), intent(out) :: y(..) !! image of same shape as x
        integer, optional, intent(out) :: status !! 0 success, 1 invalid input, 2 resource guard
        character(len=*), optional, intent(in) :: mode !! axis boundary extension
        complex(dp), optional, intent(in) :: cval !! constant exterior value
        integer, optional, intent(in) :: origin !! center shift
        complex(dp), allocatable :: flat(:),result(:),input_line(:),output_line(:)
        character(len=:), allocatable :: boundary
        integer :: dims(rank(x)),nd,n,shift,i,j,t
        integer(int64) :: stride,block,outer,o,inner,base,pos,total
        complex(dp) :: fill

        if (present(status)) status=ndimage_status_invalid_input
        nd=rank(x)
        if (nd<1 .or. nd>15 .or. rank(y)/=nd) return
        dims=shape(x)
        if (any(dims<1) .or. any(dims/=shape(y))) return
        if (axis<1 .or. axis>nd .or. filter_size<1) return
        shift=0
        if (present(origin)) shift=origin
        if (shift<-(filter_size/2) .or. shift>(filter_size-1)/2) return
        if (filter_size>max_kernel_radius) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        if (dims(axis)>huge(1)-filter_size+1) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        total=size(x,kind=int64)
        if (total>int(huge(1),int64)) then
            if (present(status)) status=ndimage_status_too_expensive
            return
        end if
        boundary='reflect'
        if (present(mode)) boundary=mode
        if (.not.mode_valid(boundary)) return
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cval
        if (.not.ieee_is_finite(real(fill,dp)) .or. .not.ieee_is_finite(aimag(fill))) return
        call pack_complex(x,flat)
        if (.not.all(ieee_is_finite(real(flat,dp))) .or. &
            .not.all(ieee_is_finite(aimag(flat)))) return
        allocate(result(size(flat)))
        n=dims(axis)
        allocate(input_line(n+filter_size-1),output_line(n))
        stride=1_int64
        do i=1,axis-1
            stride=stride*int(dims(i),int64)
        end do
        block=stride*int(n,int64)
        outer=total/block
        do o=0_int64,outer-1_int64
            do inner=0_int64,stride-1_int64
                base=o*block+inner
                do j=1,size(input_line)
                    t=extension_index(j-1-filter_size/2-shift,n,boundary)
                    input_line(j)=fill
                    if (t>0) input_line(j)=flat(base+int(t-1,int64)*stride+1_int64)
                end do
                call func(input_line,output_line)
                do j=1,n
                    pos=base+int(j-1,int64)*stride+1_int64
                    result(pos)=output_line(j)
                end do
            end do
        end do
        call scatter_complex(result,y)
        if (present(status)) status=ndimage_status_success
    end subroutine generic_filter1d_complex_nd

    subroutine median_filter_scalar(x,sizes,y,status,mode,cval,footprint,origin,modes,axes)
        real(dp), intent(in) :: x(..) !! image
        integer, intent(in) :: sizes !! common window width
        real(dp), intent(out) :: y(..) !! filtered image
        integer, optional, intent(out) :: status !! code
        character(len=*), optional, intent(in) :: mode !! extension mode
        real(dp), optional, intent(in) :: cval !! fill value
        logical, optional, intent(in) :: footprint(..) !! selected neighborhood entries
        integer, optional, intent(in) :: origin(:) !! per-axis shifted center
        character(len=*), optional, intent(in) :: modes(:) !! per-axis extension modes
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        integer, allocatable :: widths(:)
        if (present(axes)) then
            widths=spread(sizes,1,size(axes))
        else
            widths=spread(sizes,1,rank(x))
        end if
        call order_stat_work(x,widths,y,50.0_dp,.true.,status,mode,cval,footprint,origin,modes,axes)
    end subroutine median_filter_scalar

    subroutine median_filter_vector(x,sizes,y,status,mode,cval,footprint,origin,modes,axes)
        real(dp), intent(in) :: x(..) !! image
        integer, optional, intent(in) :: sizes(:) !! per-axis widths; may omit when footprint is provided
        real(dp), intent(out) :: y(..) !! filtered image
        integer, optional, intent(out) :: status !! code
        character(len=*), optional, intent(in) :: mode !! extension mode
        real(dp), optional, intent(in) :: cval !! fill value
        logical, optional, intent(in) :: footprint(..) !! selected neighborhood entries
        integer, optional, intent(in) :: origin(:) !! per-axis shifted center
        character(len=*), optional, intent(in) :: modes(:) !! per-axis extension modes
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        integer, allocatable :: widths(:)
        if (present(axes)) then
            widths=spread(1,1,size(axes))
        else
            widths=spread(1,1,rank(x))
        end if
        if (present(sizes)) then
            widths=sizes
        else if (.not.present(footprint)) then
            if (present(status)) status=ndimage_status_invalid_input
            return
        end if
        call order_stat_work(x,widths,y,50.0_dp,.true.,status,mode,cval,footprint,origin,modes,axes)
    end subroutine median_filter_vector

    subroutine rank_filter_scalar(x,rank_value,sizes,y,status,mode,cval,footprint,origin,modes,axes)
        real(dp), intent(in) :: x(..) !! image
        integer, intent(in) :: rank_value !! zero-based rank, negative from largest
        integer, intent(in) :: sizes !! common window width
        real(dp), intent(out) :: y(..) !! filtered image
        integer, optional, intent(out) :: status !! code
        character(len=*), optional, intent(in) :: mode !! extension mode
        real(dp), optional, intent(in) :: cval !! fill value
        logical, optional, intent(in) :: footprint(..) !! selected neighborhood entries
        integer, optional, intent(in) :: origin(:) !! per-axis shifted center
        character(len=*), optional, intent(in) :: modes(:) !! per-axis extension modes
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        integer, allocatable :: widths(:)
        if (present(axes)) then
            widths=spread(sizes,1,size(axes))
        else
            widths=spread(sizes,1,rank(x))
        end if
        call order_stat_work(x,widths,y,real(rank_value,dp),.false.,status,mode,cval,footprint,origin,modes,axes)
    end subroutine rank_filter_scalar

    subroutine rank_filter_vector(x,rank_value,sizes,y,status,mode,cval,footprint,origin,modes,axes)
        real(dp), intent(in) :: x(..) !! image
        integer, intent(in) :: rank_value !! zero-based rank, negative from largest
        integer, optional, intent(in) :: sizes(:) !! per-axis widths; may omit when footprint is provided
        real(dp), intent(out) :: y(..) !! filtered image
        integer, optional, intent(out) :: status !! code
        character(len=*), optional, intent(in) :: mode !! extension mode
        real(dp), optional, intent(in) :: cval !! fill value
        logical, optional, intent(in) :: footprint(..) !! selected neighborhood entries
        integer, optional, intent(in) :: origin(:) !! per-axis shifted center
        character(len=*), optional, intent(in) :: modes(:) !! per-axis extension modes
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        integer, allocatable :: widths(:)
        if (present(axes)) then
            widths=spread(1,1,size(axes))
        else
            widths=spread(1,1,rank(x))
        end if
        if (present(sizes)) then
            widths=sizes
        else if (.not.present(footprint)) then
            if (present(status)) status=ndimage_status_invalid_input
            return
        end if
        call order_stat_work(x,widths,y,real(rank_value,dp),.false.,status,mode,cval,footprint,origin,modes,axes)
    end subroutine rank_filter_vector

    subroutine percentile_filter_scalar(x,percentile,sizes,y,status,mode,cval,footprint,origin,modes,axes)
        real(dp), intent(in) :: x(..) !! image
        real(dp), intent(in) :: percentile !! percentile in [0,100]
        integer, intent(in) :: sizes !! common window width
        real(dp), intent(out) :: y(..) !! filtered image
        integer, optional, intent(out) :: status !! code
        character(len=*), optional, intent(in) :: mode !! extension mode
        real(dp), optional, intent(in) :: cval !! fill value
        logical, optional, intent(in) :: footprint(..) !! selected neighborhood entries
        integer, optional, intent(in) :: origin(:) !! per-axis shifted center
        character(len=*), optional, intent(in) :: modes(:) !! per-axis extension modes
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        integer, allocatable :: widths(:)
        if (present(axes)) then
            widths=spread(sizes,1,size(axes))
        else
            widths=spread(sizes,1,rank(x))
        end if
        call order_stat_work(x,widths,y,percentile,.true.,status,mode,cval,footprint,origin,modes,axes)
    end subroutine percentile_filter_scalar

    subroutine percentile_filter_vector(x,percentile,sizes,y,status,mode,cval,footprint,origin,modes,axes)
        real(dp), intent(in) :: x(..) !! image
        real(dp), intent(in) :: percentile !! percentile in [0,100]
        integer, optional, intent(in) :: sizes(:) !! per-axis widths; may omit when footprint is provided
        real(dp), intent(out) :: y(..) !! filtered image
        integer, optional, intent(out) :: status !! code
        character(len=*), optional, intent(in) :: mode !! extension mode
        real(dp), optional, intent(in) :: cval !! fill value
        logical, optional, intent(in) :: footprint(..) !! selected neighborhood entries
        integer, optional, intent(in) :: origin(:) !! per-axis shifted center
        character(len=*), optional, intent(in) :: modes(:) !! per-axis extension modes
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        integer, allocatable :: widths(:)
        if (present(axes)) then
            widths=spread(1,1,size(axes))
        else
            widths=spread(1,1,rank(x))
        end if
        if (present(sizes)) then
            widths=sizes
        else if (.not.present(footprint)) then
            if (present(status)) status=ndimage_status_invalid_input
            return
        end if
        call order_stat_work(x,widths,y,percentile,.true.,status,mode,cval,footprint,origin,modes,axes)
    end subroutine percentile_filter_vector

    subroutine pack_footprint(x,buf)
        logical, intent(in) :: x(..) !! arbitrary-rank Boolean footprint
        logical, allocatable, intent(out) :: buf(:) !! Fortran-order flattened mask
        select rank(x)
        rank(1)
            buf=reshape(x,[size(x)])
        rank(2)
            buf=reshape(x,[size(x)])
        rank(3)
            buf=reshape(x,[size(x)])
        rank(4)
            buf=reshape(x,[size(x)])
        rank(5)
            buf=reshape(x,[size(x)])
        rank(6)
            buf=reshape(x,[size(x)])
        rank(7)
            buf=reshape(x,[size(x)])
        rank(8)
            buf=reshape(x,[size(x)])
        rank(9)
            buf=reshape(x,[size(x)])
        rank(10)
            buf=reshape(x,[size(x)])
        rank(11)
            buf=reshape(x,[size(x)])
        rank(12)
            buf=reshape(x,[size(x)])
        rank(13)
            buf=reshape(x,[size(x)])
        rank(14)
            buf=reshape(x,[size(x)])
        rank(15)
            buf=reshape(x,[size(x)])
        rank default
            allocate(buf(0))
        end select
    end subroutine pack_footprint

    subroutine gray_stencil(nd,sizes,footprint,structure,origin,offsets,heights,status,c_order)
        integer, intent(in) :: nd !! input rank
        integer, optional, intent(in) :: sizes(:) !! rectangular window widths
        logical, optional, intent(in) :: footprint(..) !! active entries, same rank as input
        real(dp), optional, intent(in) :: structure(..) !! additive non-flat offsets
        integer, optional, intent(in) :: origin(:) !! per-axis SciPy origin offsets
        integer, allocatable, intent(out) :: offsets(:,:) !! active window offsets
        real(dp), allocatable, intent(out) :: heights(:) !! additive offsets for each active entry
        integer, intent(out) :: status !! zero on success
        logical, optional, intent(in) :: c_order !! visit footprint pixels in C order
        integer :: dims(nd),loc(nd),origins(nd),n,j,k,d,flat_index
        integer(int64) :: total,rem,findex,stride
        logical, allocatable :: active(:)
        real(dp), allocatable :: weight(:)

        status=ndimage_status_invalid_input
        dims=3
        if (present(sizes)) then
            if (size(sizes)/=nd) return
            dims=sizes
        end if
        if (present(structure)) then
            if (rank(structure)/=nd) return
            dims=shape(structure)
        end if
        if (present(footprint)) then
            if (rank(footprint)/=nd) return
            dims=shape(footprint)
            if (present(structure)) then
                if (any(dims/=shape(structure))) return
            end if
        end if
        if (any(dims<1)) return
        origins=0
        if (present(origin)) then
            if (size(origin)/=nd) return
            origins=origin
        end if
        if (any(origins<-(dims/2)) .or. any(origins>(dims-1)/2)) return
        total=1_int64
        do d=1,nd
            if (int(dims(d),int64)>50000_int64/total) then
                status=ndimage_status_too_expensive
                return
            end if
            total=total*int(dims(d),int64)
        end do
        if (present(footprint)) then
            call pack_footprint(footprint,active)
        else
            allocate(active(int(total)))
            active=.true.
        end if
        if (present(structure)) then
            call pack_real(structure,weight)
            if (.not.all(ieee_is_finite(weight))) return
        else
            allocate(weight(int(total)))
            weight=0.0_dp
        end if
        n=count(active)
        if (n==0) return
        allocate(offsets(nd,n),heights(n))
        k=0
        do j=1,size(active)
            flat_index=j
            if (present(c_order)) then
                if (c_order) then
                    rem=int(j-1,int64)
                    findex=1_int64
                    do d=nd,1,-1
                        stride=product(int(dims(:d-1),int64))
                        findex=findex+modulo(rem,int(dims(d),int64))*stride
                        rem=rem/int(dims(d),int64)
                    end do
                    flat_index=int(findex)
                end if
            end if
            if (.not.active(flat_index)) cycle
            rem=int(flat_index-1,int64)
            do d=1,nd
                loc(d)=int(modulo(rem,int(dims(d),int64)))-dims(d)/2-origins(d)
                rem=rem/int(dims(d),int64)
            end do
            k=k+1
            offsets(:,k)=loc
            heights(k)=weight(flat_index)
        end do
        status=ndimage_status_success
    end subroutine gray_stencil


    ! Construct a compact footprint on the selected axes without materializing
    ! a rank-15 singleton-extended array.  Some SciPy order-statistic and
    ! generic callbacks canonicalize footprint axes, whereas grayscale
    ! morphology associates its compact widths with the axes as supplied.
    subroutine selected_gray_stencil(nd,sizes,footprint,structure,origin,axes, &
        offsets,heights,status,c_order,canonical)
        integer, intent(in) :: nd !! source image rank
        integer, optional, intent(in) :: sizes(:) !! compact window extents
        logical, optional, intent(in) :: footprint(..) !! compact Boolean footprint
        real(dp), optional, intent(in) :: structure(..) !! compact non-flat stencil
        integer, optional, intent(in) :: origin(:) !! origins in supplied axes order
        integer, optional, intent(in) :: axes(:) !! selected one-based dimensions
        integer, allocatable, intent(out) :: offsets(:,:) !! full-rank active offsets
        real(dp), allocatable, intent(out) :: heights(:) !! corresponding stencil heights
        integer, intent(out) :: status !! success or invalid/resource error
        logical, optional, intent(in) :: c_order !! sample order for callbacks
        logical, optional, intent(in) :: canonical !! sort selected axes for compact footprint
        integer, allocatable :: reduced(:,:),axis_order(:),reduced_origin(:)
        integer :: i,j,tmp
        logical :: selected(nd),ok

        status=ndimage_status_invalid_input
        call selected_axes_mask(nd,axes,selected,ok)
        if (.not.ok) return
        if (.not.present(axes)) then
            call gray_stencil(nd,sizes,footprint,structure,origin,offsets,heights,status,c_order)
            return
        end if
        if (present(sizes)) then
            if (size(sizes)/=size(axes)) return
        end if
        if (present(origin)) then
            if (size(origin)/=size(axes)) return
        end if
        if (size(axes)==0) then
            if (present(footprint) .or. present(structure)) return
            if (present(sizes)) then
                if (size(sizes)/=0) return
            end if
            allocate(offsets(nd,1),source=0)
            allocate(heights(1),source=0.0_dp)
            status=ndimage_status_success
            return
        end if
        axis_order=axes
        if (present(canonical)) then
            if (canonical) then
                do i=1,size(axis_order)
                    do j=i+1,size(axis_order)
                        if (axis_order(j)>=axis_order(i)) cycle
                        tmp=axis_order(j)
                        axis_order(j)=axis_order(i)
                        axis_order(i)=tmp
                    end do
                end do
            end if
        end if
        if (present(origin)) then
            allocate(reduced_origin(size(axes)))
            do i=1,size(axes)
                do j=1,size(axes)
                    if (axis_order(i)==axes(j)) reduced_origin(i)=origin(j)
                end do
            end do
        end if
        call gray_stencil(size(axes),sizes,footprint,structure,reduced_origin, &
            reduced,heights,status,c_order)
        if (status/=ndimage_status_success) return
        allocate(offsets(nd,size(reduced,2)),source=0)
        do i=1,size(axis_order)
            offsets(axis_order(i),:)=reduced(i,:)
        end do
    end subroutine selected_gray_stencil

    subroutine gray_kernel(input,dims,offsets,heights,boundaries,cval,dilate,output)
        real(dp), intent(in) :: input(:) !! Fortran packed scalar input
        integer, intent(in) :: dims(:) !! extents
        integer, intent(in) :: offsets(:,:) !! active window offsets
        real(dp), intent(in) :: heights(:) !! non-flat values matching offsets
        character(len=*), intent(in) :: boundaries(:) !! boundary mode per axis
        real(dp), intent(in) :: cval !! constant exterior value
        logical, intent(in) :: dilate !! true for max, false for min
        real(dp), intent(out) :: output(:) !! filtered image
        integer(int64) :: strides(size(dims)),pos,rem,q
        integer :: loc(size(dims)),d,k,ix,side
        real(dp) :: value,best
        logical :: outside

        strides(1)=1_int64
        do d=2,size(dims)
            strides(d)=strides(d-1)*int(dims(d-1),int64)
        end do
        do pos=1_int64,size(input,kind=int64)
            rem=pos-1_int64
            do d=1,size(dims)
                loc(d)=int(modulo(rem,int(dims(d),int64)))
                rem=rem/int(dims(d),int64)
            end do
            if (dilate) then
                best=-huge(1.0_dp)
                side=-1
            else
                best=huge(1.0_dp)
                side=1
            end if
            do k=1,size(offsets,2)
                outside=.false.
                q=1_int64
                do d=1,size(dims)
                    ix=extension_index(loc(d)+side*offsets(d,k),dims(d),trim(boundaries(d)))
                    if (ix==0) then
                        outside=.true.
                        exit
                    end if
                    q=q+int(ix-1,int64)*strides(d)
                end do
                value=cval
                if (.not.outside) value=input(q)
                if (dilate) then
                    value=value+heights(k)
                else
                    value=value-heights(k)
                end if
                if (dilate) then
                    best=max(best,value)
                else
                    best=min(best,value)
                end if
            end do
            ! An empty footprint leaves the input unchanged.
            if (size(offsets,2)==0) best=input(pos)
            output(pos)=best
        end do
    end subroutine gray_kernel

    subroutine gray_operation(x,y,operation,status,sizes,footprint,mode,cval,structure,origin,axes,modes)
        real(dp), intent(in) :: x(..) !! image, ranks 1 through 15
        real(dp), intent(out) :: y(..) !! processed image of the same shape
        integer, intent(in) :: operation !! 1 dilation, 2 erosion, 3 opening, 4 closing,
                                       !! 5 gradient, 6 laplace, 7 white top-hat, 8 black top-hat
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: sizes(:) !! footprint size per axis
        logical, optional, intent(in) :: footprint(..) !! active neighbors
        character(len=*), optional, intent(in) :: mode !! boundary extension
        real(dp), optional, intent(in) :: cval !! value outside for constant modes
        real(dp), optional, intent(in) :: structure(..) !! non-flat additive stencil
        integer, optional, intent(in) :: origin(:) !! per-axis origin displacement
        integer, optional, intent(in) :: axes(:) !! selected axes, supplied order
        character(len=*), optional, intent(in) :: modes(:) !! modes on selected axes
        integer, allocatable :: offsets(:,:),dilation_offsets(:,:),compact_dims(:), &
            axis_order(:),adjusted_origin(:),requested_origin(:)
        real(dp), allocatable :: heights(:),dilation_heights(:),flat(:),r(:),t(:)
        logical, allocatable :: mask(:)
        logical :: nonseparable
        integer :: d,j,k,tmp,neven,meven
        character(len=32), allocatable :: boundaries(:)
        logical :: valid_modes
        integer :: code
        real(dp) :: fill

        if (present(status)) status=ndimage_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x)) return
        if (any(shape(x)/=shape(y)) .or. any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        if (operation<1 .or. operation>8) return
        call selected_filter_modes(rank(x),axes,mode,modes,boundaries,valid_modes)
        if (.not.valid_modes) return
        fill=0.0_dp
        if (present(cval)) fill=cval
        if (.not.ieee_is_finite(fill)) return
        ! SciPy's non-separable footprint expansion inserts singleton axes
        ! in increasing dimension order, whereas separable window sizes are
        ! associated with the user-supplied axes order.  Dilation additionally
        ! reverses even-width origin shifts before the expansion.  Build the
        ! erosion and dilation stencils separately when these differ.
        nonseparable=present(structure)
        if (present(footprint) .and. .not.present(structure)) then
            call pack_footprint(footprint,mask)
            nonseparable=.not.all(mask)
        end if
        if (present(axes) .and. nonseparable) then
            if (present(footprint)) then
                compact_dims=shape(footprint)
            else
                compact_dims=shape(structure)
            end if
            if (size(compact_dims)/=size(axes)) return
            allocate(requested_origin(size(axes)),source=0)
            if (present(origin)) then
                if (size(origin)/=size(axes)) return
                requested_origin=origin
            end if
            if (any(requested_origin < -(compact_dims/2)) .or. &
                any(requested_origin > (compact_dims-1)/2)) return
            axis_order=axes
            do d=1,size(axis_order)
                do j=d+1,size(axis_order)
                    if (axis_order(j)>=axis_order(d)) cycle
                    tmp=axis_order(j)
                    axis_order(j)=axis_order(d)
                    axis_order(d)=tmp
                end do
            end do
            adjusted_origin=requested_origin
            do j=1,size(axes)
                do k=1,size(axis_order)
                    if (axes(j)/=axis_order(k)) cycle
                    neven=modulo(compact_dims(j)+1,2)
                    meven=modulo(compact_dims(k)+1,2)
                    adjusted_origin(j)=requested_origin(j)+neven-meven
                end do
            end do
            call selected_gray_stencil(rank(x),sizes,footprint,structure,requested_origin,axes, &
                offsets,heights,code,canonical=.true.)
            if (code==ndimage_status_success) then
                call selected_gray_stencil(rank(x),sizes,footprint,structure,adjusted_origin,axes, &
                    dilation_offsets,dilation_heights,code,canonical=.true.)
            end if
        else
            call selected_gray_stencil(rank(x),sizes,footprint,structure,origin,axes, &
                offsets,heights,code)
            if (code==ndimage_status_success) then
                dilation_offsets=offsets
                dilation_heights=heights
            end if
        end if
        if (code/=ndimage_status_success) then
            if (present(status)) status=code
            return
        end if
        call pack_real(x,flat)
        if (.not.all(ieee_is_finite(flat))) return
        allocate(r(size(flat)),t(size(flat)))
        select case(operation)
        case(1)
            call gray_kernel(flat,shape(x),dilation_offsets,dilation_heights,boundaries,fill,.true.,r)
        case(2)
            call gray_kernel(flat,shape(x),offsets,heights,boundaries,fill,.false.,r)
        case(3,7)
            call gray_kernel(flat,shape(x),offsets,heights,boundaries,fill,.false.,t)
            call gray_kernel(t,shape(x),dilation_offsets,dilation_heights,boundaries,fill,.true.,r)
            if (operation==7) r=flat-r
        case(4,8)
            call gray_kernel(flat,shape(x),dilation_offsets,dilation_heights,boundaries,fill,.true.,t)
            call gray_kernel(t,shape(x),offsets,heights,boundaries,fill,.false.,r)
            if (operation==8) r=r-flat
        case(5,6)
            call gray_kernel(flat,shape(x),dilation_offsets,dilation_heights,boundaries,fill,.true.,r)
            call gray_kernel(flat,shape(x),offsets,heights,boundaries,fill,.false.,t)
            if (operation==5) then
                r=r-t
            else
                r=r+t-2.0_dp*flat
            end if
        end select
        call scatter_real(r,y)
        if (present(status)) status=ndimage_status_success
    end subroutine gray_operation

    subroutine grey_dilation_nd(x,y,status,sizes,footprint,mode,cval,structure,origin,axes,modes)
        real(dp), intent(in) :: x(..) !! real-valued multidimensional image
        real(dp), intent(out) :: y(..) !! caller-allocated output
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: sizes(:) !! per-axis rectangular sizes
        logical, optional, intent(in) :: footprint(..) !! custom active neighborhood
        character(len=*), optional, intent(in) :: mode !! boundary mode
        real(dp), optional, intent(in) :: cval !! constant boundary value
        real(dp), optional, intent(in) :: structure(..) !! non-flat structuring element
        integer, optional, intent(in) :: origin(:) !! per-axis origin offsets
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        character(len=*), optional, intent(in) :: modes(:) !! boundary mode per selected axis
        call gray_operation(x,y,1,status,sizes,footprint,mode,cval,structure,origin,axes,modes)
    end subroutine grey_dilation_nd

    subroutine grey_erosion_nd(x,y,status,sizes,footprint,mode,cval,structure,origin,axes,modes)
        real(dp), intent(in) :: x(..) !! real-valued multidimensional image
        real(dp), intent(out) :: y(..) !! caller-allocated output
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: sizes(:) !! per-axis rectangular sizes
        logical, optional, intent(in) :: footprint(..) !! custom active neighborhood
        character(len=*), optional, intent(in) :: mode !! boundary mode
        real(dp), optional, intent(in) :: cval !! constant boundary value
        real(dp), optional, intent(in) :: structure(..) !! non-flat structuring element
        integer, optional, intent(in) :: origin(:) !! per-axis origin offsets
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        character(len=*), optional, intent(in) :: modes(:) !! boundary mode per selected axis
        call gray_operation(x,y,2,status,sizes,footprint,mode,cval,structure,origin,axes,modes)
    end subroutine grey_erosion_nd

    subroutine grey_opening_nd(x,y,status,sizes,footprint,mode,cval,structure,origin,axes,modes)
        real(dp), intent(in) :: x(..) !! real-valued multidimensional image
        real(dp), intent(out) :: y(..) !! caller-allocated output
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: sizes(:) !! per-axis rectangular sizes
        logical, optional, intent(in) :: footprint(..) !! custom active neighborhood
        character(len=*), optional, intent(in) :: mode !! boundary mode
        real(dp), optional, intent(in) :: cval !! constant boundary value
        real(dp), optional, intent(in) :: structure(..) !! non-flat structuring element
        integer, optional, intent(in) :: origin(:) !! per-axis origin offsets
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        character(len=*), optional, intent(in) :: modes(:) !! boundary mode per selected axis
        call gray_operation(x,y,3,status,sizes,footprint,mode,cval,structure,origin,axes,modes)
    end subroutine grey_opening_nd

    subroutine grey_closing_nd(x,y,status,sizes,footprint,mode,cval,structure,origin,axes,modes)
        real(dp), intent(in) :: x(..) !! real-valued multidimensional image
        real(dp), intent(out) :: y(..) !! caller-allocated output
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: sizes(:) !! per-axis rectangular sizes
        logical, optional, intent(in) :: footprint(..) !! custom active neighborhood
        character(len=*), optional, intent(in) :: mode !! boundary mode
        real(dp), optional, intent(in) :: cval !! constant boundary value
        real(dp), optional, intent(in) :: structure(..) !! non-flat structuring element
        integer, optional, intent(in) :: origin(:) !! per-axis origin offsets
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        character(len=*), optional, intent(in) :: modes(:) !! boundary mode per selected axis
        call gray_operation(x,y,4,status,sizes,footprint,mode,cval,structure,origin,axes,modes)
    end subroutine grey_closing_nd

    subroutine morphological_gradient_nd(x,y,status,sizes,footprint,mode,cval,structure,origin,axes,modes)
        real(dp), intent(in) :: x(..) !! real-valued multidimensional image
        real(dp), intent(out) :: y(..) !! caller-allocated output
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: sizes(:) !! per-axis rectangular sizes
        logical, optional, intent(in) :: footprint(..) !! custom active neighborhood
        character(len=*), optional, intent(in) :: mode !! boundary mode
        real(dp), optional, intent(in) :: cval !! constant boundary value
        real(dp), optional, intent(in) :: structure(..) !! non-flat structuring element
        integer, optional, intent(in) :: origin(:) !! per-axis origin offsets
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        character(len=*), optional, intent(in) :: modes(:) !! boundary mode per selected axis
        call gray_operation(x,y,5,status,sizes,footprint,mode,cval,structure,origin,axes,modes)
    end subroutine morphological_gradient_nd

    subroutine morphological_laplace_nd(x,y,status,sizes,footprint,mode,cval,structure,origin,axes,modes)
        real(dp), intent(in) :: x(..) !! real-valued multidimensional image
        real(dp), intent(out) :: y(..) !! caller-allocated output
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: sizes(:) !! per-axis rectangular sizes
        logical, optional, intent(in) :: footprint(..) !! custom active neighborhood
        character(len=*), optional, intent(in) :: mode !! boundary mode
        real(dp), optional, intent(in) :: cval !! constant boundary value
        real(dp), optional, intent(in) :: structure(..) !! non-flat structuring element
        integer, optional, intent(in) :: origin(:) !! per-axis origin offsets
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        character(len=*), optional, intent(in) :: modes(:) !! boundary mode per selected axis
        call gray_operation(x,y,6,status,sizes,footprint,mode,cval,structure,origin,axes,modes)
    end subroutine morphological_laplace_nd

    subroutine white_tophat_nd(x,y,status,sizes,footprint,mode,cval,structure,origin,axes,modes)
        real(dp), intent(in) :: x(..) !! real-valued multidimensional image
        real(dp), intent(out) :: y(..) !! caller-allocated output
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: sizes(:) !! per-axis rectangular sizes
        logical, optional, intent(in) :: footprint(..) !! custom active neighborhood
        character(len=*), optional, intent(in) :: mode !! boundary mode
        real(dp), optional, intent(in) :: cval !! constant boundary value
        real(dp), optional, intent(in) :: structure(..) !! non-flat structuring element
        integer, optional, intent(in) :: origin(:) !! per-axis origin offsets
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        character(len=*), optional, intent(in) :: modes(:) !! boundary mode per selected axis
        call gray_operation(x,y,7,status,sizes,footprint,mode,cval,structure,origin,axes,modes)
    end subroutine white_tophat_nd

    subroutine black_tophat_nd(x,y,status,sizes,footprint,mode,cval,structure,origin,axes,modes)
        real(dp), intent(in) :: x(..) !! real-valued multidimensional image
        real(dp), intent(out) :: y(..) !! caller-allocated output
        integer, optional, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: sizes(:) !! per-axis rectangular sizes
        logical, optional, intent(in) :: footprint(..) !! custom active neighborhood
        character(len=*), optional, intent(in) :: mode !! boundary mode
        real(dp), optional, intent(in) :: cval !! constant boundary value
        real(dp), optional, intent(in) :: structure(..) !! non-flat structuring element
        integer, optional, intent(in) :: origin(:) !! per-axis origin offsets
        integer, optional, intent(in) :: axes(:) !! one-based selected axes
        character(len=*), optional, intent(in) :: modes(:) !! boundary mode per selected axis
        call gray_operation(x,y,8,status,sizes,footprint,mode,cval,structure,origin,axes,modes)
    end subroutine black_tophat_nd

    ! Frequency-domain operations multiply the supplied FFT spectrum, without
    ! computing a forward or inverse transform. For n>0, axis identifies the
    ! one-based real FFT axis, whose stored size must be n/2+1.
    subroutine fourier_process(x,params,y,operation,status,n,axis)
        complex(dp), intent(in) :: x(..) !! complete FFT or RFFT spectrum
        real(dp), intent(in) :: params(:) !! per-axis sigma, width or translation
        complex(dp), intent(out) :: y(..) !! filtered frequency-domain array
        integer, intent(in) :: operation !! 1 Gaussian, 2 uniform, 3 shift
        integer, optional, intent(out) :: status !! zero success, one invalid
        integer, optional, intent(in) :: n !! -1 complex FFT, >0 real FFT length
        integer, optional, intent(in) :: axis !! one-based RFFT axis, default last
        complex(dp), allocatable :: flat(:),res(:)
        integer :: nd,d,fft_length,real_axis,rfft_size
        integer :: dims(rank(x))
        integer(int64) :: pos,rem,k,ksigned
        integer :: nfft
        real(dp) :: freq,omega,acc,factor,pi
        complex(dp) :: multiplier

        if (present(status)) status=ndimage_status_invalid_input
        nd=rank(x)
        if (nd<1 .or. nd>15 .or. rank(y)/=nd) return
        dims=shape(x)
        if (any(dims<1) .or. any(dims/=shape(y))) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        if (size(params)/=nd) return
        if (.not.all(ieee_is_finite(params))) return
        if (operation==1 .and. any(params<0.0_dp)) return
        if (operation<1 .or. operation>3) return
        nfft=-1
        if (present(n)) nfft=n
        if (nfft==0 .or. nfft< -1) return
        real_axis=nd
        if (present(axis)) real_axis=axis
        if (real_axis<0) real_axis=nd+real_axis+1
        if (real_axis<1 .or. real_axis>nd) return
        if (nfft>0) then
            rfft_size=nfft/2+1
            if (dims(real_axis)/=rfft_size) return
        end if
        call pack_complex(x,flat)
        if (.not.all(ieee_is_finite(real(flat,dp))) .or. &
            .not.all(ieee_is_finite(aimag(flat)))) return
        allocate(res(size(flat)))
        pi=acos(-1.0_dp)
        do pos=1_int64,size(flat,kind=int64)
            rem=pos-1_int64
            acc=0.0_dp
            factor=1.0_dp
            do d=1,nd
                k=modulo(rem,int(dims(d),int64))
                rem=rem/int(dims(d),int64)
                fft_length=dims(d)
                ksigned=k
                if (nfft>0 .and. d==real_axis) then
                    fft_length=nfft
                else
                    if (2_int64*k>=int(fft_length,int64)) ksigned=k-int(fft_length,int64)
                end if
                freq=real(ksigned,dp)/real(fft_length,dp)
                select case(operation)
                case(1)
                    acc=acc+(freq*params(d))**2
                case(2)
                    omega=pi*freq*params(d)
                    if (omega/=0.0_dp) factor=factor*sin(omega)/omega
                case(3)
                    acc=acc+freq*params(d)
                end select
            end do
            select case(operation)
            case(1)
                multiplier=cmplx(exp(-2.0_dp*pi*pi*acc),0.0_dp,dp)
            case(2)
                multiplier=cmplx(factor,0.0_dp,dp)
            case(3)
                omega=2.0_dp*pi*acc
                multiplier=cmplx(cos(omega),-sin(omega),dp)
            end select
            res(pos)=flat(pos)*multiplier
        end do
        call scatter_complex(res,y)
        if (present(status)) status=ndimage_status_success
    end subroutine fourier_process

    subroutine fourier_gaussian_vector(x,sigma,y,status,n,axis)
        complex(dp), intent(in) :: x(..) !! frequency-domain image
        real(dp), intent(in) :: sigma(:) !! Gaussian sigmas per axis
        complex(dp), intent(out) :: y(..) !! filtered spectrum
        integer, optional, intent(out) :: status !! status code
        integer, optional, intent(in) :: n !! real FFT length, or -1
        integer, optional, intent(in) :: axis !! real FFT axis, one-based
        call fourier_process(x,sigma,y,1,status,n,axis)
    end subroutine fourier_gaussian_vector

    subroutine fourier_gaussian_scalar(x,sigma,y,status,n,axis)
        complex(dp), intent(in) :: x(..) !! frequency-domain image
        real(dp), intent(in) :: sigma !! shared Gaussian sigma
        complex(dp), intent(out) :: y(..) !! filtered spectrum
        integer, optional, intent(out) :: status !! status code
        integer, optional, intent(in) :: n !! real FFT length, or -1
        integer, optional, intent(in) :: axis !! real FFT axis, one-based
        call fourier_process(x,spread(sigma,1,rank(x)),y,1,status,n,axis)
    end subroutine fourier_gaussian_scalar

    subroutine fourier_uniform_vector(x,sizes,y,status,n,axis)
        complex(dp), intent(in) :: x(..) !! frequency-domain image
        real(dp), intent(in) :: sizes(:) !! continuous uniform widths per axis
        complex(dp), intent(out) :: y(..) !! filtered spectrum
        integer, optional, intent(out) :: status !! status code
        integer, optional, intent(in) :: n !! real FFT length, or -1
        integer, optional, intent(in) :: axis !! real FFT axis, one-based
        call fourier_process(x,sizes,y,2,status,n,axis)
    end subroutine fourier_uniform_vector

    subroutine fourier_uniform_scalar(x,sizes,y,status,n,axis)
        complex(dp), intent(in) :: x(..) !! frequency-domain image
        real(dp), intent(in) :: sizes !! uniform width for each axis
        complex(dp), intent(out) :: y(..) !! filtered spectrum
        integer, optional, intent(out) :: status !! status code
        integer, optional, intent(in) :: n !! real FFT length, or -1
        integer, optional, intent(in) :: axis !! real FFT axis, one-based
        call fourier_process(x,spread(sizes,1,rank(x)),y,2,status,n,axis)
    end subroutine fourier_uniform_scalar

    subroutine fourier_shift_vector(x,shift,y,status,n,axis)
        complex(dp), intent(in) :: x(..) !! frequency-domain image
        real(dp), intent(in) :: shift(:) !! translations per axis, in pixels
        complex(dp), intent(out) :: y(..) !! shifted spectrum
        integer, optional, intent(out) :: status !! status code
        integer, optional, intent(in) :: n !! real FFT length, or -1
        integer, optional, intent(in) :: axis !! real FFT axis, one-based
        call fourier_process(x,shift,y,3,status,n,axis)
    end subroutine fourier_shift_vector

    subroutine fourier_shift_scalar(x,shift,y,status,n,axis)
        complex(dp), intent(in) :: x(..) !! frequency-domain image
        real(dp), intent(in) :: shift !! translation on each axis
        complex(dp), intent(out) :: y(..) !! shifted spectrum
        integer, optional, intent(out) :: status !! status code
        integer, optional, intent(in) :: n !! real FFT length, or -1
        integer, optional, intent(in) :: axis !! real FFT axis, one-based
        call fourier_process(x,spread(shift,1,rank(x)),y,3,status,n,axis)
    end subroutine fourier_shift_scalar

    ! Normalized Fourier transform of a uniform ellipsoid (an N-D ball after
    ! per-axis coordinate scaling). SciPy implements ranks 1--3; the same
    ! mathematical expression extends to higher ranks in SciFort.
    ! For dimension D and radius parameter q=pi*||size*frequency||,
    !     F_D(q) = Gamma(D/2+1) * (2/q)**(D/2) * J_(D/2)(q).
    ! A power series avoids cancellation around q=0. For larger q,
    ! integer- and half-integer-order Bessel functions are evaluated directly.
    pure real(dp) function ellipsoid_response(q,nd) result(value)
        real(dp), intent(in) :: q !! dimensionless Fourier radius
        integer, intent(in) :: nd !! array rank, 1 through 15
        real(dp) :: term,total,jnm1,jn,jnnext,scale
        integer :: k,m

        if (q==0.0_dp) then
            value=1.0_dp
        else if (q<8.0_dp) then
            term=1.0_dp
            total=term
            do k=1,256
                term=term*(-0.25_dp*q*q)/(real(k,dp)*(0.5_dp*real(nd,dp)+real(k,dp)))
                total=total+term
                if (abs(term)<=4.0_dp*epsilon(total)*max(1.0_dp,abs(total))) exit
            end do
            value=total
        else if (mod(nd,2)==0) then
            m=nd/2
            scale=1.0_dp
            do k=1,m
                scale=scale*2.0_dp*real(k,dp)/q
            end do
            value=scale*bessel_jn(m,q)
        else
            m=(nd-1)/2
            jnm1=sin(q)/q
            if (m==0) then
                jn=jnm1
            else
                jn=sin(q)/(q*q)-cos(q)/q
                do k=1,m-1
                    jnnext=real(2*k+1,dp)*jn/q-jnm1
                    jnm1=jn
                    jn=jnnext
                end do
            end if
            scale=1.0_dp
            do k=1,m
                scale=scale*real(2*k+1,dp)/q
            end do
            value=scale*jn
        end if
    end function ellipsoid_response

    subroutine fourier_ellipsoid_vector(x,sizes,y,status,n,axis)
        complex(dp), intent(in) :: x(..) !! Fourier coefficients, ranks 1 to 15
        real(dp), intent(in) :: sizes(:) !! diameter of ellipsoid along each axis
        complex(dp), intent(out) :: y(..) !! ellipsoid-filtered spectrum
        integer, optional, intent(out) :: status !! success or invalid input
        integer, optional, intent(in) :: n !! packed real-FFT original extent or -1
        integer, optional, intent(in) :: axis !! packed real-FFT axis (one-based)
        complex(dp), allocatable :: flat(:),res(:)
        integer :: nd,d,fft_length,real_axis,nfft,dims(rank(x))
        integer(int64) :: pos,rem,k,ksigned
        real(dp) :: freq,radius2,pi,q

        if (present(status)) status=ndimage_status_invalid_input
        nd=rank(x)
        if (nd<1 .or. nd>15 .or. rank(y)/=nd) return
        dims=shape(x)
        if (any(dims<1) .or. any(dims/=shape(y))) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        if (size(sizes)/=nd) return
        if (.not.all(ieee_is_finite(sizes))) return
        if (any(abs(sizes)>1.0e6_dp)) return
        nfft=-1
        if (present(n)) nfft=n
        if (nfft==0 .or. nfft< -1) return
        real_axis=nd
        if (present(axis)) real_axis=axis
        if (real_axis<0) real_axis=nd+real_axis+1
        if (real_axis<1 .or. real_axis>nd) return
        if (nfft>0) then
            if (dims(real_axis)/=nfft/2+1) return
        end if
        call pack_complex(x,flat)
        if (.not.all(ieee_is_finite(real(flat,dp))) .or. &
            .not.all(ieee_is_finite(aimag(flat)))) return
        allocate(res(size(flat)))
        pi=acos(-1.0_dp)
        do pos=1_int64,size(flat,kind=int64)
            rem=pos-1_int64
            radius2=0.0_dp
            do d=1,nd
                k=modulo(rem,int(dims(d),int64))
                rem=rem/int(dims(d),int64)
                fft_length=dims(d)
                ksigned=k
                if (nfft>0 .and. d==real_axis) then
                    fft_length=nfft
                else
                    if (2_int64*k>=int(fft_length,int64)) ksigned=k-int(fft_length,int64)
                end if
                freq=real(ksigned,dp)*sizes(d)/real(fft_length,dp)
                radius2=radius2+freq*freq
            end do
            q=pi*sqrt(radius2)
            res(pos)=flat(pos)*ellipsoid_response(q,nd)
        end do
        call scatter_complex(res,y)
        if (present(status)) status=ndimage_status_success
    end subroutine fourier_ellipsoid_vector

    subroutine fourier_ellipsoid_scalar(x,sizes,y,status,n,axis)
        complex(dp), intent(in) :: x(..) !! Fourier coefficients
        real(dp), intent(in) :: sizes !! common diameter for each axis
        complex(dp), intent(out) :: y(..) !! filtered spectrum
        integer, optional, intent(out) :: status !! success or invalid input
        integer, optional, intent(in) :: n !! packed real-FFT original extent
        integer, optional, intent(in) :: axis !! packed axis (one-based)
        call fourier_ellipsoid_vector(x,spread(sizes,1,rank(x)),y,status,n,axis)
    end subroutine fourier_ellipsoid_scalar

end module scifort_ndimage
