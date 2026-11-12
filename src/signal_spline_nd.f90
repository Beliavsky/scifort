! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Rank-independent spline prefiltering and coordinate interpolation.
! Independent implementation based on cardinal B-spline collocation.
! SciPy ndimage boundary conventions; all coordinates are one-based.
! Public coordinates are ONE-BASED. See docs/signal_spline_nd.md.
module scifort_signal_spline_nd
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_spline, only : spline_status_success, spline_status_invalid_input
    implicit none
    private

    ! The high-rank output remains defined on numerical-complexity failure.
    integer, parameter, public :: spline_nd_status_too_expensive = 2
    integer(int64), parameter :: max_evaluation_terms = 2000000_int64
    integer(int64), parameter :: max_padded_elements = 8000000_int64
    integer, parameter :: spline_halo = 12
    integer(int64), parameter :: max_geometric_coordinates = 4000000_int64

    public :: spline_filter_nd, spline_filter1d, cspline_nd, qspline_nd
    public :: map_coordinates_nd, map_coordinates_grid_nd
    public :: geometric_transform_nd, ndimage_coordinate_map
    public :: affine_transform_nd, shift_nd
    public :: zoom_nd, rotate_nd, zoom_shape_nd, rotate_shape_nd
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
    public :: ndimage_coordinate_batch_map, evaluate_spline_plan_geometric_batch_nd
    public :: evaluate_spline_plan_displacement_nd, evaluate_spline_plan_displacement_grid_nd
    public :: evaluate_spline_plan_displacement_with_gradient_nd

    ! Reuse filtered coefficients across many unrelated coordinate queries.
    ! Components are private so callers cannot invalidate cached metadata.
    type :: spline_interpolation_plan
        private
        complex(dp), allocatable :: coefficients(:)
        integer, allocatable :: original_shape(:), working_shape(:), halo(:)
        integer :: degree = -1
        character(len=16) :: boundary = 'mirror'
        complex(dp) :: fill = (0.0_dp,0.0_dp)
        logical :: filtered = .false.
        logical :: complex_input = .false.
    end type spline_interpolation_plan

    ! Coordinates along each output axis for a tensor-product interpolation.
    ! The result has one dimension per input axis, avoiding a full coordinate mesh.
    type :: spline_axis_grid
        real(dp), allocatable :: coordinates(:)
    end type spline_axis_grid

    interface evaluate_spline_plan_with_hessian_nd
        module procedure evaluate_plan_hessian_real, evaluate_plan_hessian_complex
    end interface
    interface evaluate_spline_plan_grid_with_hessian_nd
        module procedure evaluate_plan_grid_hessian_real, evaluate_plan_grid_hessian_complex
    end interface
    interface evaluate_spline_plan_affine_with_hessian_nd
        module procedure evaluate_plan_affine_hessian_real
        module procedure evaluate_plan_affine_hessian_complex
    end interface
    interface evaluate_spline_plan_displacement_with_hessian_nd
        module procedure evaluate_plan_displacement_hessian_real
        module procedure evaluate_plan_displacement_hessian_complex
    end interface
    interface evaluate_spline_plan_displacement_grid_with_hessian_nd
        module procedure evaluate_plan_displacement_grid_hessian_real
        module procedure evaluate_plan_displacement_grid_hessian_complex
    end interface
    interface evaluate_spline_plan_displacement_grid_with_gradient_nd
        module procedure evaluate_plan_displacement_grid_gradient_real
        module procedure evaluate_plan_displacement_grid_gradient_complex
    end interface
    interface evaluate_spline_plan_grid_with_gradient_nd
        module procedure evaluate_plan_grid_gradient_real, evaluate_plan_grid_gradient_complex
    end interface
    interface evaluate_spline_plan_affine_with_gradient_nd
        module procedure evaluate_plan_affine_gradient_real, evaluate_plan_affine_gradient_complex
    end interface
    interface evaluate_spline_plan_tensor_nd
        module procedure evaluate_plan_tensor_real, evaluate_plan_tensor_complex
    end interface

    interface evaluate_spline_plan_displacement_with_gradient_nd
        module procedure evaluate_plan_displacement_gradient_real
        module procedure evaluate_plan_displacement_gradient_complex
    end interface
    interface evaluate_spline_plan_displacement_nd
        module procedure evaluate_plan_displacement_real, evaluate_plan_displacement_complex
    end interface
    interface evaluate_spline_plan_displacement_grid_nd
        module procedure evaluate_plan_displacement_grid_real, evaluate_plan_displacement_grid_complex
    end interface

    interface evaluate_spline_plan_geometric_batch_nd
        module procedure evaluate_plan_geometric_batch_real, evaluate_plan_geometric_batch_complex
    end interface

    interface evaluate_spline_plan_affine_nd
        module procedure evaluate_plan_affine_real, evaluate_plan_affine_complex
    end interface
    interface evaluate_spline_plan_geometric_nd
        module procedure evaluate_plan_geometric_real, evaluate_plan_geometric_complex
    end interface

    interface prepare_spline_plan_nd
        module procedure prepare_plan_real, prepare_plan_complex
    end interface
    interface evaluate_spline_plan_nd
        module procedure evaluate_plan_real, evaluate_plan_complex
    end interface
    interface evaluate_spline_plan_with_gradient_nd
        module procedure evaluate_plan_gradient_real, evaluate_plan_gradient_complex
    end interface
    interface evaluate_spline_plan_grid_nd
        module procedure evaluate_plan_grid_real, evaluate_plan_grid_complex
    end interface

    interface spline_filter_nd
        module procedure filter_all_real, filter_all_complex
    end interface
    interface spline_filter1d
        module procedure filter_axis_real, filter_axis_complex
    end interface
    interface cspline_nd
        module procedure cubic_real, cubic_complex
    end interface
    interface qspline_nd
        module procedure quadratic_real, quadratic_complex
    end interface
    interface map_coordinates_nd
        module procedure coordinates_real, coordinates_complex
    end interface
    interface map_coordinates_grid_nd
        module procedure coordinate_grid_real, coordinate_grid_complex
    end interface
    abstract interface
        subroutine ndimage_coordinate_map(output_coords,input_coords)
            import dp
            integer, intent(in) :: output_coords(:) !! one-based destination coordinates
            real(dp), intent(out) :: input_coords(:) !! one-based source coordinates
        end subroutine ndimage_coordinate_map

        ! A whole block of one-based coordinates per callback invocation.
        ! The number of columns may be smaller in the final block.
        subroutine ndimage_coordinate_batch_map(output_coords,input_coords)
            import dp
            integer, intent(in) :: output_coords(:,:) !! (output rank, number of points)
            real(dp), intent(out) :: input_coords(:,:) !! (input rank, number of points)
        end subroutine ndimage_coordinate_batch_map
    end interface
    interface geometric_transform_nd
        module procedure geometric_real, geometric_complex
    end interface
    interface affine_transform_nd
        module procedure affine_real, affine_complex
    end interface
    interface shift_nd
        module procedure shift_real, shift_complex
    end interface
    interface zoom_nd
        module procedure zoom_real, zoom_complex, zoom_real_scalar, zoom_complex_scalar
    end interface
    interface zoom_shape_nd
        module procedure zoom_shape_vector, zoom_shape_scalar
    end interface
    interface rotate_nd
        module procedure rotate_real, rotate_complex
    end interface

contains
    pure integer function mirrored(raw, n) result(idx)
        integer, intent(in) :: raw !! zero-based integer coordinate
        integer, intent(in) :: n !! extent of axis
        integer :: t
        idx = 1
        if (n <= 1) return
        t = modulo(raw, 2*(n-1))
        if (t >= n) t = 2*(n-1)-t
        idx = t+1
    end function mirrored

    pure integer function boundary_index(raw,n,mode) result(idx)
        integer, intent(in) :: raw !! zero-based index
        integer, intent(in) :: n !! axis extent
        character(len=*), intent(in) :: mode !! boundary extension
        integer :: t

        idx=0
        select case (mode)
        case ('reflect','grid-mirror')
            t=modulo(raw,2*n)
            if (t>=n) t=2*n-1-t
            idx=t+1
        case ('nearest')
            idx=max(1,min(n,raw+1))
        case ('grid-wrap')
            idx=modulo(raw,n)+1
        case ('grid-constant')
            if (raw>=0 .and. raw<n) idx=raw+1
        case default
            idx=mirrored(raw,n)
        end select
    end function boundary_index

    pure logical function valid_boundary(mode) result(ok)
        character(len=*), intent(in) :: mode
        select case (mode)
        case ('mirror','constant','reflect','grid-mirror','nearest', &
              'wrap','grid-wrap','grid-constant')
            ok=.true.
        case default
            ok=.false.
        end select
    end function valid_boundary

    pure real(dp) function cardinal(degree, t) result(v)
        integer, intent(in) :: degree !! cardinal spline degree, zero through five
        real(dp), intent(in) :: t !! coordinate measured from knot
        real(dp) :: a
        integer :: k
        real(dp) :: vals(0:5, 0:5)

        ! Recurrence of centered B-splines. The half-value at a zeroth-
        ! order discontinuity gives correct integer-knot evaluations.
        vals = 0.0_dp
        do k = 0, degree
            a = t + 0.5_dp*real(degree,dp) - real(k,dp)
            if (abs(a) < 0.5_dp) vals(0,k) = 1.0_dp
            if (abs(abs(a)-0.5_dp) <= 4*epsilon(a)) vals(0,k) = 0.5_dp
        end do
        do k = 1, degree
            block
                integer :: i
                real(dp) :: u
                do i = 0, degree-k
                    u = t + real(degree-k,dp)*0.5_dp-real(i,dp)
                    vals(k,i) = ((real(k+1,dp)*0.5_dp+u)*vals(k-1,i) + &
                        (real(k+1,dp)*0.5_dp-u)*vals(k-1,i+1)) / real(k,dp)
                end do
            end block
        end do
        v = vals(degree,0)
    end function cardinal

    ! Periodic B-spline prefilter using the causal/anticausal pole recurrences.
    ! Unlike a direct DFT, this requires O(n) arithmetic per line for orders
    ! 2 through 5. The poles follow the cardinal spline characteristic roots.
    subroutine periodic_line_coeff(samples,degree,output,status)
        complex(dp), intent(in) :: samples(:) !! periodic samples
        integer, intent(in) :: degree !! spline order 2 through 5
        complex(dp), intent(out) :: output(:) !! periodic spline coefficients
        integer, intent(out) :: status !! numerical status
        real(dp) :: poles(2), pole, gain, power
        integer :: n, i, j, count

        status=spline_status_invalid_input
        n=size(samples)
        if (size(output)/=n .or. n<1) return
        poles=0.0_dp
        select case(degree)
        case(2)
            poles(1)=-0.171572875253809902396622551580603843_dp
            count=1
        case(3)
            poles(1)=-0.267949192431122706472553658494127633_dp
            count=1
        case(4)
            poles(1)=-0.361341225900220177092212841325675255_dp
            poles(2)=-0.013725429297339121360331226939128204_dp
            count=2
        case(5)
            poles(1)=-0.430575347099973791851434783493520110_dp
            poles(2)=-0.043096288203264653822712376822550182_dp
            count=2
        case default
            return
        end select
        output=samples
        if (n==1) then
            status=spline_status_success
            return
        end if
        gain=1.0_dp
        do j=1,count
            pole=poles(j)
            gain=gain*(1.0_dp-pole)*(1.0_dp-1.0_dp/pole)
        end do
        output=output*gain
        do j=1,count
            pole=poles(j)
            power=pole
            do i=1,n-1
                output(1)=output(1)+power*output(n-i+1)
                power=power*pole
            end do
            output(1)=output(1)/(1.0_dp-power)
            do i=2,n
                output(i)=output(i)+pole*output(i-1)
            end do
            power=pole
            do i=1,n-1
                output(n)=output(n)+power*output(i)
                power=power*pole
            end do
            output(n)=output(n)*pole/(power-1.0_dp)
            do i=n-1,1,-1
                output(i)=pole*(output(i+1)-output(i))
            end do
        end do
        status=spline_status_success
    end subroutine periodic_line_coeff

    subroutine line_coeff(samples, degree, output, status, mode)
        complex(dp), intent(in) :: samples(:) !! equally spaced samples
        integer, intent(in) :: degree !! 2, 3, 4, or 5
        complex(dp), intent(out) :: output(:) !! spline coefficients
        integer, intent(out) :: status !! numerical status
        character(len=*), optional, intent(in) :: mode !! spline boundary mode
        complex(dp), allocatable :: rhs(:)
        real(dp), allocatable :: a(:,:)
        real(dp) :: pivot, factor, wt
        integer :: n, i, j, k, dest
        character(len=:), allocatable :: boundary

        status = spline_status_invalid_input
        n = size(samples)
        if (size(output) /= n .or. n < 1) return
        if (degree < 2 .or. degree > 5) return
        boundary='mirror'
        if (present(mode)) boundary=mode
        if (boundary=='grid-wrap') then
            call periodic_line_coeff(samples,degree,output,status)
            return
        end if
        ! SciPy's spline_filter1d modes differ from interpolation padding:
        ! nearest/reflect use half-sample symmetry; constant/wrap/mirror
        ! and grid-constant use whole-sample symmetry.
        select case(boundary)
        case('reflect','grid-mirror','nearest')
            boundary='reflect'
        case default
            boundary='mirror'
        end select
        if (n == 1) then
            output = samples
            status = spline_status_success
            return
        end if
        allocate(a(5,n), rhs(n))
        a = 0.0_dp
        rhs = samples
        do i = 1, n
            do k = -2, 2
                wt = cardinal(degree,real(k,dp))
                if (wt == 0.0_dp) cycle
                dest = boundary_index(i-1+k,n,boundary)
                if (dest==0) cycle
                a(dest-i+3,i) = a(dest-i+3,i) + wt
            end do
        end do
        ! Banded Gaussian elimination without pivoting. These cardinal
        ! spline collocation matrices are strictly row-diagonally dominant.
        do i = 1, n
            pivot = a(3,i)
            if (abs(pivot) <= tiny(pivot)) return
            do j = i+1, min(n,i+2)
                factor = a(i-j+3,j)/pivot
                a(i-j+3,j) = 0.0_dp
                if (factor == 0.0_dp) cycle
                do k = i+1, min(n,i+2,j+2)
                    a(k-j+3,j) = a(k-j+3,j) - factor*a(k-i+3,i)
                end do
                rhs(j) = rhs(j) - factor*rhs(i)
            end do
        end do
        do i = n, 1, -1
            do j = i+1, min(n,i+2)
                rhs(i) = rhs(i) - a(j-i+3,i)*output(j)
            end do
            output(i) = rhs(i)/a(3,i)
        end do
        status = spline_status_success
    end subroutine line_coeff

    subroutine filter_flat(values, dims, degree, axis, status, mode)
        complex(dp), intent(inout) :: values(:) !! column-major packed data
        integer, intent(in) :: dims(:) !! array extents
        integer, intent(in) :: degree !! order 0 through 5
        integer, intent(in) :: axis !! zero means all axes
        integer, intent(out) :: status !! numerical status
        character(len=*), optional, intent(in) :: mode !! spline boundary
        complex(dp), allocatable :: line(:), result(:)
        integer(int64) :: stride, outer, block, i, j, pos
        integer :: d, first, last, n, k, code

        status = spline_status_invalid_input
        if (degree < 0 .or. degree > 5) return
        if (axis < 0 .or. axis > size(dims)) return
        if (size(dims) == 0 .or. any(dims < 1)) return
        if (degree <= 1) then
            status = spline_status_success
            return
        end if
        first = 1
        last = size(dims)
        if (axis /= 0) then
            first = axis
            last = axis
        end if
        do d = first, last
            n = dims(d)
            if (n == 1) cycle
            allocate(line(n),result(n))
            stride = product(int(dims(:d-1),int64))
            block = stride*int(n,int64)
            outer = product(int(dims(d+1:),int64))
            do j = 0_int64, outer-1_int64
                do i = 0_int64, stride-1_int64
                    do k = 1, n
                        pos = j*block+i+1_int64+int(k-1,int64)*stride
                        line(k) = values(pos)
                    end do
                    call line_coeff(line,degree,result,code,mode)
                    if (code /= spline_status_success) return
                    do k = 1, n
                        pos = j*block+i+1_int64+int(k-1,int64)*stride
                        values(pos) = result(k)
                    end do
                end do
            end do
            deallocate(line,result)
        end do
        status = spline_status_success
    end subroutine filter_flat

    ! SciPy uses a 12-sample data halo for high-order nearest and
    ! grid-constant interpolation. Padding every singleton axis would make
    ! rank-15 data impractical. Singleton nearest axes are constant, while
    ! a singleton grid-constant axis contributes a separable impulse factor.
    subroutine prepare_interpolation_coeff(data,dims,degree,mode,fill,coeff, &
                                             expanded,halo,status)
        complex(dp), intent(in) :: data(:) !! packed original samples
        integer, intent(in) :: dims(:) !! input extents
        integer, intent(in) :: degree !! B-spline order
        character(len=*), intent(in) :: mode !! interpolation boundary mode
        complex(dp), intent(in) :: fill !! constant extension value
        complex(dp), allocatable, intent(out) :: coeff(:) !! packed spline coefficients
        integer, allocatable, intent(out) :: expanded(:) !! working extents
        integer, allocatable, intent(out) :: halo(:) !! original-to-working offsets
        integer, intent(out) :: status !! spline status
        integer(int64) :: total, p, q, source, multiplier
        integer :: d, coordinate, original, n
        logical :: outside

        status=spline_status_invalid_input
        n=size(dims)
        expanded=dims
        allocate(halo(n))
        halo=0
        if (n<1 .or. any(dims<1)) return
        if (mode=='nearest' .or. mode=='grid-constant') then
            where (dims>1) halo=spline_halo
            if (any(dims>huge(1)-2*spline_halo)) return
            expanded=dims+2*halo
            total=1_int64
            do d=1,n
                if (total>max_padded_elements/int(expanded(d),int64)) then
                    status=spline_nd_status_too_expensive
                    return
                end if
                total=total*int(expanded(d),int64)
            end do
            allocate(coeff(total))
            coeff=(0.0_dp,0.0_dp)
            do p=0_int64,total-1_int64
                q=p
                source=1_int64
                multiplier=1_int64
                outside=.false.
                do d=1,n
                    coordinate=int(modulo(q,int(expanded(d),int64)))-halo(d)
                    q=q/int(expanded(d),int64)
                    if (mode=='grid-constant' .and. &
                        (coordinate<0 .or. coordinate>=dims(d))) then
                        outside=.true.
                        exit
                    end if
                    original=max(0,min(dims(d)-1,coordinate))
                    source=source+int(original,int64)*multiplier
                    multiplier=multiplier*int(dims(d),int64)
                end do
                if (outside) cycle
                if (mode=='grid-constant') then
                    coeff(p+1_int64)=data(source)-fill
                else
                    coeff(p+1_int64)=data(source)
                end if
            end do
            ! Filtering the padded input at the selected boundary mode
            ! also produces the coefficients used by SciPy outside the
            ! original data rectangle.
            if (mode=='nearest') then
                call filter_flat(coeff,expanded,degree,0,status,'reflect')
            else
                call filter_flat(coeff,expanded,degree,0,status,'mirror')
            end if
        else
            coeff=data
            if (mode=='reflect' .or. mode=='grid-mirror' .or. &
                mode=='grid-wrap') then
                call filter_flat(coeff,expanded,degree,0,status,mode)
            else
                call filter_flat(coeff,expanded,degree,0,status)
            end if
        end if
    end subroutine prepare_interpolation_coeff

    subroutine evaluate_interpolation(coeff,expanded,halo,original,coords, &
                                       degree,mode,fill,values,status)
        complex(dp), intent(in) :: coeff(:) !! prefiltered packed coefficients
        integer, intent(in) :: expanded(:) !! padded array shape
        integer, intent(in) :: halo(:) !! coordinate shifts
        integer, intent(in) :: original(:) !! original input extents
        real(dp), intent(in) :: coords(:,:) !! one-based sample positions
        integer, intent(in) :: degree !! spline degree
        character(len=*), intent(in) :: mode !! chosen boundary mode
        complex(dp), intent(in) :: fill !! constant fill value
        complex(dp), intent(out) :: values(:) !! interpolated samples
        integer, intent(out) :: status !! zero on success
        real(dp), allocatable :: shifted(:,:)
        complex(dp) :: impulse(25), impulse_coeff(25)
        real(dp) :: weights(6), weight_factor
        integer :: indices(6), nterms, d, k, j, code

        shifted=coords
        do d=1,size(expanded)
            shifted(d,:)=shifted(d,:)+real(halo(d),dp)
        end do
        if (mode/='grid-constant' .or. .not.any(original==1)) then
            if (mode=='grid-constant') then
                call interpolate(coeff,expanded,shifted,values,degree,mode, &
                                 (0.0_dp,0.0_dp),status)
                if (status==spline_status_success) values=values+fill
            else
                call interpolate(coeff,expanded,shifted,values,degree,mode,fill,status)
            end if
            return
        end if
        ! An axis of length one contributes a tensor-product interpolation
        ! factor from the padded unit impulse, with zero exterior samples.
        ! This avoids a 25-fold allocation penalty per singleton dimension.
        impulse=(0.0_dp,0.0_dp)
        impulse(13)=(1.0_dp,0.0_dp)
        call line_coeff(impulse,degree,impulse_coeff,code,'mirror')
        status=code
        if (code/=spline_status_success) return
        call interpolate(coeff,expanded,shifted,values,degree,mode, &
                         (0.0_dp,0.0_dp),status,collapse_singletons=.true.)
        if (status/=spline_status_success) return
        do k=1,size(values)
            weight_factor=1.0_dp
            do d=1,size(original)
                if (original(d)/=1) cycle
                call support(coords(d,k)+12.0_dp,25,degree,mode,indices,weights,nterms)
                block
                    real(dp) :: w
                    w=0.0_dp
                    do j=1,nterms
                        if (indices(j)==0) cycle
                        w=w+weights(j)*real(impulse_coeff(indices(j)),dp)
                    end do
                    weight_factor=weight_factor*w
                end block
            end do
            values(k)=fill+weight_factor*values(k)
        end do
    end subroutine evaluate_interpolation

    subroutine support(x,n,degree,mode,indices,weights,num,dweights,ddweights)
        real(dp), intent(in) :: x !! one-based coordinate
        integer, intent(in) :: n !! input axis extent
        integer, intent(in) :: degree !! 0 through 5
        character(len=*), intent(in) :: mode !! boundary extension
        integer, intent(out) :: indices(6) !! one-based coefficient indices, zero for cval
        real(dp), intent(out) :: weights(6) !! combined B-spline weights
        integer, intent(out) :: num !! number of distinct neighbors
        real(dp), optional, intent(out) :: dweights(6) !! derivatives with respect to x
        real(dp), optional, intent(out) :: ddweights(6) !! second derivatives with respect to x
        real(dp) :: u,period,wt,z,slope,dwt,ddwt,t
        integer :: start,raw,k,idx,j,nsupport
        logical :: linear_knot

        indices=0
        weights=0.0_dp
        if (present(dweights)) dweights=0.0_dp
        if (present(ddweights)) ddweights=0.0_dp
        num=0
        z=x-1.0_dp
        u=z
        slope=1.0_dp
        select case(mode)
        case('mirror','constant')
            if (n>1) then
                period=2.0_dp*real(n-1,dp)
                u=modulo(z,period)
                if (u>real(n-1,dp)) then
                    u=period-u
                    slope=-1.0_dp
                end if
            else
                u=0.0_dp
                slope=0.0_dp
            end if
        case('wrap')
            if (n>1) then
                period=real(n-1,dp)
                if (z<0.0_dp .or. z>period) then
                    u=modulo(z,period)
                    if (z<0.0_dp .and. u==0.0_dp) u=period
                end if
            else
                u=0.0_dp
                slope=0.0_dp
            end if
        case('reflect','grid-mirror')
            period=2.0_dp*real(n,dp)
            u=modulo(z+0.5_dp,period)-0.5_dp
            if (u>=real(n,dp)-0.5_dp) then
                u=period-1.0_dp-u
                slope=-1.0_dp
            end if
        case('grid-wrap')
            u=modulo(z,real(n,dp))
        case('nearest','grid-constant')
            ! Avoid converting far-outside coordinates to an overflowing INTEGER.
            if (z<-real(degree+3,dp) .or. &
                z>real(n+degree+3,dp)) then
                u=max(-real(degree+3,dp), &
                    min(real(n+degree+3,dp),z))
                slope=0.0_dp
            end if
        end select
        if (mod(degree,2)==0) then
            start=floor(u+0.5_dp)-degree/2
            if (degree==0 .and. mode=='mirror' .and. n>1) then
                period=2.0_dp*real(n-1,dp)
                if (z>=0.0_dp .and. modulo(z,period)>real(n-1,dp)) then
                    if (abs(u-real(floor(u),dp)-0.5_dp)<=8*epsilon(u)) &
                        start=start-1
                end if
            end if
        else
            start=floor(u)-degree/2
        end if
        ! The first-order spline has a corner at integral coordinates.
        ! Average its one-sided derivatives there by including both
        ! neighboring endpoint weights; value-only support is unchanged.
        linear_knot=.false.
        if ((present(dweights) .or. present(ddweights)) .and. degree==1) then
            linear_knot=u==real(floor(u),dp)
        end if
        nsupport=degree+1
        if (linear_knot) then
            start=start-1
            nsupport=3
        end if
        do k=0,nsupport-1
            raw=start+k
            idx=boundary_index(raw,n,mode)
            t=u-real(raw,dp)
            wt=cardinal(degree,t)
            if (degree==0) wt=1.0_dp
            dwt=0.0_dp
            if (present(dweights) .and. degree>0) then
                dwt=slope*(cardinal(degree-1,t+0.5_dp)- &
                    cardinal(degree-1,t-0.5_dp))
            end if
            ddwt=0.0_dp
            if (present(ddweights) .and. degree>=2) then
                ddwt=slope*slope*(cardinal(degree-2,t+1.0_dp)- &
                    2.0_dp*cardinal(degree-2,t)+cardinal(degree-2,t-1.0_dp))
            end if
            if (wt==0.0_dp .and. dwt==0.0_dp .and. ddwt==0.0_dp) cycle
            do j=1,num
                if (indices(j)==idx) exit
            end do
            if (j<=num) then
                weights(j)=weights(j)+wt
                if (present(dweights)) dweights(j)=dweights(j)+dwt
                if (present(ddweights)) ddweights(j)=ddweights(j)+ddwt
            else
                num=num+1
                indices(num)=idx
                weights(num)=wt
                if (present(dweights)) dweights(num)=dwt
                if (present(ddweights)) ddweights(num)=ddwt
            end if
        end do
    end subroutine support

    subroutine interpolate(values,dims,coords,result,degree,mode,fill,status, &
                           collapse_singletons)
        complex(dp), intent(in) :: values(:) !! packed coefficients
        integer, intent(in) :: dims(:) !! original shape
        real(dp), intent(in) :: coords(:,:) !! (rank, npoints) one-based
        complex(dp), intent(out) :: result(:) !! interpolated points
        integer, intent(in) :: degree !! spline order 0..5
        character(len=*), intent(in) :: mode !! mirror or constant
        complex(dp), intent(in) :: fill !! constant out-of-bounds value
        integer, intent(out) :: status !! error code
        logical, optional, intent(in) :: collapse_singletons !! ignore singleton support
        integer :: ix(6,size(dims)), counts(size(dims)), picks(size(dims))
        integer(int64) :: strides(size(dims)), nterms, off
        real(dp) :: wt(6,size(dims)), total
        integer :: d,k,npts,kfill

        status = spline_status_invalid_input
        npts = size(coords,2)
        result = cmplx(0.0_dp,0.0_dp,dp)
        if (degree < 0 .or. degree > 5) return
        if (size(coords,1) /= size(dims) .or. size(result) /= npts) return
        if (.not. valid_boundary(mode)) return
        if (any(dims < 1)) return
        if (.not. all(ieee_is_finite(coords))) return
        if (.not. ieee_is_finite(real(fill,dp)) .or. &
            .not. ieee_is_finite(aimag(fill))) return
        strides(1)=1_int64
        do d = 2, size(dims)
            strides(d) = strides(d-1)*int(dims(d-1),int64)
        end do
        do k = 1, npts
            if (mode == 'constant') then
                if (any(coords(:,k)<1.0_dp) .or. &
                    any(coords(:,k)>real(dims,dp))) then
                    result(k) = fill
                    cycle
                end if
            end if
            nterms = 1_int64
            do d = 1, size(dims)
                if (present(collapse_singletons)) then
                    if (collapse_singletons .and. dims(d)==1) then
                        ix(:,d)=0
                        ix(1,d)=1
                        wt(:,d)=0.0_dp
                        wt(1,d)=1.0_dp
                        counts(d)=1
                    else
                        call support(coords(d,k),dims(d),degree,mode, &
                                     ix(:,d),wt(:,d),counts(d))
                    end if
                else
                    call support(coords(d,k),dims(d),degree,mode, &
                                 ix(:,d),wt(:,d),counts(d))
                end if
                if (counts(d) < 1) exit
                if (nterms > max_evaluation_terms / int(counts(d),int64)) then
                    status=spline_nd_status_too_expensive
                    return
                end if
                nterms=nterms*int(counts(d),int64)
            end do
            if (d <= size(dims)) then
                status=spline_status_invalid_input
                return
            end if
            picks=1
            do
                total=1.0_dp
                off=1_int64
                do d = 1, size(dims)
                    total=total*wt(picks(d),d)
                    if (ix(picks(d),d)==0) exit
                    off=off+int(ix(picks(d),d)-1,int64)*strides(d)
                end do
                if (d<=size(dims)) then
                    do kfill=d+1,size(dims)
                        total=total*wt(picks(kfill),kfill)
                    end do
                    result(k)=result(k)+total*fill
                else
                    result(k)=result(k)+total*values(off)
                end if
                do d = 1, size(dims)
                    picks(d)=picks(d)+1
                    if (picks(d) <= counts(d)) exit
                    picks(d)=1
                end do
                if (d > size(dims)) exit
            end do
        end do
        status=spline_status_success
    end subroutine interpolate

    ! Interpolate both the value and its coordinate derivatives in one support
    ! traversal. B-spline derivatives use B'_(p)(t) = B_(p-1)(t+1/2)
    ! - B_(p-1)(t-1/2). Prefix/suffix products avoid division by zero.
    subroutine interpolate_with_gradient(values,dims,coords,result,grad,degree,mode,fill,status, &
                                         collapse_singletons)
        complex(dp), intent(in) :: values(:) !! packed coefficient array
        integer, intent(in) :: dims(:) !! coefficient extents
        real(dp), intent(in) :: coords(:,:) !! one-based positions (rank, points)
        complex(dp), intent(out) :: result(:) !! interpolated values
        complex(dp), intent(out) :: grad(:,:) !! coordinate derivatives (rank, points)
        integer, intent(in) :: degree !! spline order 0..5
        character(len=*), intent(in) :: mode !! extension mode
        complex(dp), intent(in) :: fill !! constant boundary value
        integer, intent(out) :: status !! success or invalid/complexity error
        logical, optional, intent(in) :: collapse_singletons !! ignore singleton factors
        integer :: ix(6,size(dims)),counts(size(dims)),picks(size(dims))
        integer(int64) :: strides(size(dims)),nterms,off
        real(dp) :: wt(6,size(dims)),dwt(6,size(dims))
        real(dp) :: prefix(0:size(dims)),suffix(size(dims)+1)
        complex(dp) :: sample
        integer :: d,k,rank_input,npts

        status=spline_status_invalid_input
        result=(0.0_dp,0.0_dp)
        grad=(0.0_dp,0.0_dp)
        rank_input=size(dims)
        npts=size(coords,2)
        if (rank_input<1 .or. rank_input>15) return
        if (degree<0 .or. degree>5) return
        if (size(coords,1)/=rank_input .or. size(result)/=npts) return
        if (size(grad,1)/=rank_input .or. size(grad,2)/=npts) return
        if (.not.valid_boundary(mode) .or. any(dims<1)) return
        if (.not.all(ieee_is_finite(coords))) return
        if (.not.ieee_is_finite(real(fill,dp)) .or. .not.ieee_is_finite(aimag(fill))) return
        strides(1)=1_int64
        do d=2,rank_input
            strides(d)=strides(d-1)*int(dims(d-1),int64)
        end do
        do k=1,npts
            if (mode=='constant') then
                if (any(coords(:,k)<1.0_dp) .or. &
                    any(coords(:,k)>real(dims,dp))) then
                    result(k)=fill
                    cycle
                end if
            end if
            nterms=1_int64
            do d=1,rank_input
                if (present(collapse_singletons)) then
                    if (collapse_singletons .and. dims(d)==1) then
                        ix(:,d)=0
                        ix(1,d)=1
                        wt(:,d)=0.0_dp
                        wt(1,d)=1.0_dp
                        dwt(:,d)=0.0_dp
                        counts(d)=1
                    else
                        call support(coords(d,k),dims(d),degree,mode, &
                            ix(:,d),wt(:,d),counts(d),dwt(:,d))
                    end if
                else
                    call support(coords(d,k),dims(d),degree,mode, &
                        ix(:,d),wt(:,d),counts(d),dwt(:,d))
                end if
                if (counts(d)<1) return
                if (nterms>max_evaluation_terms/int(counts(d),int64)) then
                    status=spline_nd_status_too_expensive
                    return
                end if
                nterms=nterms*int(counts(d),int64)
            end do
            picks=1
            do
                off=1_int64
                sample=fill
                do d=1,rank_input
                    if (ix(picks(d),d)==0) exit
                    off=off+int(ix(picks(d),d)-1,int64)*strides(d)
                end do
                if (d>rank_input) sample=values(off)
                prefix(0)=1.0_dp
                do d=1,rank_input
                    prefix(d)=prefix(d-1)*wt(picks(d),d)
                end do
                suffix(rank_input+1)=1.0_dp
                do d=rank_input,1,-1
                    suffix(d)=suffix(d+1)*wt(picks(d),d)
                end do
                result(k)=result(k)+prefix(rank_input)*sample
                do d=1,rank_input
                    grad(d,k)=grad(d,k)+prefix(d-1)*suffix(d+1)*dwt(picks(d),d)*sample
                end do
                do d=1,rank_input
                    picks(d)=picks(d)+1
                    if (picks(d)<=counts(d)) exit
                    picks(d)=1
                end do
                if (d>rank_input) exit
            end do
        end do
        status=spline_status_success
    end subroutine interpolate_with_gradient

    subroutine evaluate_interpolation_with_gradient(plan,coords,values,grad,status)
        type(spline_interpolation_plan), intent(in) :: plan !! prepared spline plan
        real(dp), intent(in) :: coords(:,:) !! one-based original coordinates
        complex(dp), intent(out) :: values(:) !! interpolated values
        complex(dp), intent(out) :: grad(:,:) !! partial derivatives
        integer, intent(out) :: status !! zero for success
        real(dp), allocatable :: shifted(:,:)
        complex(dp) :: impulse(25),impulse_coeff(25)
        real(dp) :: w(6),dw(6),factors(size(plan%original_shape))
        real(dp) :: dfactors(size(plan%original_shape)),total
        integer :: idx(6),nterms,d,k,j,code,e

        shifted=coords
        do d=1,size(plan%original_shape)
            shifted(d,:)=shifted(d,:)+real(plan%halo(d),dp)
        end do
        if (plan%boundary/='grid-constant' .or. .not.any(plan%original_shape==1)) then
            if (plan%boundary=='grid-constant') then
                call interpolate_with_gradient(plan%coefficients,plan%working_shape, &
                    shifted,values,grad,plan%degree,trim(plan%boundary), &
                    (0.0_dp,0.0_dp),status)
                if (status==spline_status_success) values=values+plan%fill
            else
                call interpolate_with_gradient(plan%coefficients,plan%working_shape, &
                    shifted,values,grad,plan%degree,trim(plan%boundary),plan%fill,status)
            end if
            return
        end if
        impulse=(0.0_dp,0.0_dp)
        impulse(13)=(1.0_dp,0.0_dp)
        call line_coeff(impulse,plan%degree,impulse_coeff,code,'mirror')
        status=code
        if (code/=spline_status_success) return
        call interpolate_with_gradient(plan%coefficients,plan%working_shape, &
            shifted,values,grad,plan%degree,trim(plan%boundary), &
            (0.0_dp,0.0_dp),status,collapse_singletons=.true.)
        if (status/=spline_status_success) return
        do k=1,size(values)
            factors=1.0_dp
            dfactors=0.0_dp
            do d=1,size(plan%original_shape)
                if (plan%original_shape(d)/=1) cycle
                call support(coords(d,k)+12.0_dp,25,plan%degree, &
                    trim(plan%boundary),idx,w,nterms,dw)
                factors(d)=0.0_dp
                do j=1,nterms
                    if (idx(j)==0) cycle
                    factors(d)=factors(d)+w(j)*real(impulse_coeff(idx(j)),dp)
                    dfactors(d)=dfactors(d)+dw(j)*real(impulse_coeff(idx(j)),dp)
                end do
            end do
            total=product(factors)
            do d=1,size(plan%original_shape)
                if (plan%original_shape(d)/=1) then
                    grad(d,k)=grad(d,k)*total
                else
                    block
                        real(dp) :: g
                        g=dfactors(d)
                        do e=1,size(plan%original_shape)
                            if (e/=d) g=g*factors(e)
                        end do
                        grad(d,k)=values(k)*g
                    end block
                end if
            end do
            values(k)=plan%fill+total*values(k)
        end do
    end subroutine evaluate_interpolation_with_gradient

    subroutine evaluate_plan_gradient_core(plan,coords,values,grad,status)
        type(spline_interpolation_plan), intent(in) :: plan !! prepared plan
        real(dp), intent(in) :: coords(:,:) !! one-based positions (rank, npoints)
        complex(dp), intent(out) :: values(:) !! interpolated values
        complex(dp), intent(out) :: grad(:,:) !! coordinate derivatives
        integer, intent(out) :: status !! zero on success

        status=spline_status_invalid_input
        values=(0.0_dp,0.0_dp)
        grad=(0.0_dp,0.0_dp)
        if (.not.allocated(plan%coefficients)) return
        if (size(coords,1)/=size(plan%original_shape)) return
        if (size(coords,2)/=size(values)) return
        if (size(grad,1)/=size(plan%original_shape) .or. &
            size(grad,2)/=size(values)) return
        if (plan%filtered) then
            call evaluate_interpolation_with_gradient(plan,coords,values,grad,status)
        else
            call interpolate_with_gradient(plan%coefficients,plan%original_shape, &
                coords,values,grad,plan%degree,trim(plan%boundary),plan%fill,status)
        end if
    end subroutine evaluate_plan_gradient_core

    ! This extension returns analytic derivatives with respect to the
    ! one-based source coordinates (not derivatives of the input samples).
    ! Order zero has a zero derivative away from nearest-neighbor jumps.
    subroutine evaluate_plan_gradient_real(plan,coords,values,grad,status)
        type(spline_interpolation_plan), intent(in) :: plan !! real-valued prepared plan
        real(dp), intent(in) :: coords(:,:) !! one-based query coordinates
        real(dp), intent(out) :: values(:) !! scalar interpolated values
        real(dp), intent(out) :: grad(:,:) !! (source rank, npoints) Jacobian
        integer, intent(out) :: status !! zero on success
        complex(dp), allocatable :: work(:),derivatives(:,:)

        status=spline_status_invalid_input
        values=0.0_dp
        grad=0.0_dp
        if (.not.allocated(plan%coefficients)) return
        if (plan%complex_input) return
        allocate(work(size(values)),derivatives(size(grad,1),size(grad,2)))
        call evaluate_plan_gradient_core(plan,coords,work,derivatives,status)
        if (status==spline_status_success) then
            values=real(work,dp)
            grad=real(derivatives,dp)
        end if
    end subroutine evaluate_plan_gradient_real

    subroutine evaluate_plan_gradient_complex(plan,coords,values,grad,status)
        type(spline_interpolation_plan), intent(in) :: plan !! complex-valued prepared plan
        real(dp), intent(in) :: coords(:,:) !! one-based query coordinates
        complex(dp), intent(out) :: values(:) !! scalar interpolated values
        complex(dp), intent(out) :: grad(:,:) !! (source rank, npoints) derivatives
        integer, intent(out) :: status !! zero on success

        status=spline_status_invalid_input
        values=(0.0_dp,0.0_dp)
        grad=(0.0_dp,0.0_dp)
        if (.not.allocated(plan%coefficients)) return
        if (.not.plan%complex_input) return
        call evaluate_plan_gradient_core(plan,coords,values,grad,status)
    end subroutine evaluate_plan_gradient_complex

    ! The Hessian is analytic in the interiors of spline polynomial pieces.
    ! At knots where a derivative is discontinuous, the lower-degree cardinal
    ! convention provides a symmetric average (orders zero and one give zero
    ! classical second derivatives away from jumps).
    subroutine interpolate_with_hessian(values,dims,coords,result,grad,hess,degree,mode,fill,status, &
                                        collapse_singletons)
        complex(dp), intent(in) :: values(:) !! packed interpolation coefficients
        integer, intent(in) :: dims(:) !! coefficient extents
        real(dp), intent(in) :: coords(:,:) !! one-based (source rank, points)
        complex(dp), intent(out) :: result(:) !! interpolated values
        complex(dp), intent(out) :: grad(:,:) !! first coordinate derivatives
        complex(dp), intent(out) :: hess(:,:,:) !! symmetric second derivatives
        integer, intent(in) :: degree !! spline order from zero to five
        character(len=*), intent(in) :: mode !! extension mode
        complex(dp), intent(in) :: fill !! constant exterior value
        integer, intent(out) :: status !! success or an error code
        logical, optional, intent(in) :: collapse_singletons !! avoid singleton expansion
        integer :: ix(6,size(dims)), counts(size(dims)), picks(size(dims))
        integer(int64) :: strides(size(dims)), nterms, off
        real(dp) :: wt(6,size(dims)), dwt(6,size(dims)), ddwt(6,size(dims))
        real(dp) :: prefix(0:size(dims)), suffix(size(dims)+1), middle, term
        complex(dp) :: sample
        integer :: d,e,k,nr,np

        status=spline_status_invalid_input
        result=(0.0_dp,0.0_dp)
        grad=(0.0_dp,0.0_dp)
        hess=(0.0_dp,0.0_dp)
        nr=size(dims)
        np=size(coords,2)
        if (nr<1 .or. nr>15 .or. degree<0 .or. degree>5) return
        if (any(dims<1) .or. .not.valid_boundary(mode)) return
        if (size(coords,1)/=nr .or. size(result)/=np) return
        if (any(shape(grad)/=[nr,np]) .or. any(shape(hess)/=[nr,nr,np])) return
        if (.not.all(ieee_is_finite(coords))) return
        if (.not.ieee_is_finite(real(fill,dp)) .or. .not.ieee_is_finite(aimag(fill))) return
        strides(1)=1_int64
        do d=2,nr
            strides(d)=strides(d-1)*int(dims(d-1),int64)
        end do
        do k=1,np
            if (mode=='constant') then
                if (any(coords(:,k)<1.0_dp) .or. any(coords(:,k)>real(dims,dp))) then
                    result(k)=fill
                    cycle
                end if
            end if
            nterms=1_int64
            do d=1,nr
                if (present(collapse_singletons)) then
                    if (collapse_singletons .and. dims(d)==1) then
                        ix(:,d)=0
                        ix(1,d)=1
                        wt(:,d)=0.0_dp
                        wt(1,d)=1.0_dp
                        dwt(:,d)=0.0_dp
                        ddwt(:,d)=0.0_dp
                        counts(d)=1
                    else
                        call support(coords(d,k),dims(d),degree,mode, &
                            ix(:,d),wt(:,d),counts(d),dwt(:,d),ddwt(:,d))
                    end if
                else
                    call support(coords(d,k),dims(d),degree,mode, &
                        ix(:,d),wt(:,d),counts(d),dwt(:,d),ddwt(:,d))
                end if
                if (counts(d)<1) return
                if (nterms>max_evaluation_terms/int(counts(d),int64)) then
                    status=spline_nd_status_too_expensive
                    return
                end if
                nterms=nterms*int(counts(d),int64)
            end do
            picks=1
            do
                off=1_int64
                sample=fill
                do d=1,nr
                    if (ix(picks(d),d)==0) exit
                    off=off+int(ix(picks(d),d)-1,int64)*strides(d)
                end do
                if (d>nr) sample=values(off)
                prefix(0)=1.0_dp
                do d=1,nr
                    prefix(d)=prefix(d-1)*wt(picks(d),d)
                end do
                suffix(nr+1)=1.0_dp
                do d=nr,1,-1
                    suffix(d)=suffix(d+1)*wt(picks(d),d)
                end do
                result(k)=result(k)+prefix(nr)*sample
                do d=1,nr
                    term=prefix(d-1)*suffix(d+1)
                    grad(d,k)=grad(d,k)+term*dwt(picks(d),d)*sample
                    hess(d,d,k)=hess(d,d,k)+term*ddwt(picks(d),d)*sample
                    middle=1.0_dp
                    do e=d+1,nr
                        term=prefix(d-1)*dwt(picks(d),d)*middle* &
                            dwt(picks(e),e)*suffix(e+1)
                        hess(d,e,k)=hess(d,e,k)+term*sample
                        middle=middle*wt(picks(e),e)
                    end do
                end do
                do d=1,nr
                    picks(d)=picks(d)+1
                    if (picks(d)<=counts(d)) exit
                    picks(d)=1
                end do
                if (d>nr) exit
            end do
            do d=1,nr
                do e=d+1,nr
                    hess(e,d,k)=hess(d,e,k)
                end do
            end do
        end do
        status=spline_status_success
    end subroutine interpolate_with_hessian

    subroutine evaluate_interpolation_with_hessian(plan,coords,values,grad,hess,status)
        type(spline_interpolation_plan), intent(in) :: plan !! cached coefficients
        real(dp), intent(in) :: coords(:,:) !! original, one-based coordinates
        complex(dp), intent(out) :: values(:) !! interpolated samples
        complex(dp), intent(out) :: grad(:,:) !! first derivatives
        complex(dp), intent(out) :: hess(:,:,:) !! second derivatives
        integer, intent(out) :: status !! spline error code
        real(dp), allocatable :: shifted(:,:)
        complex(dp) :: impulse(25), impulse_coeff(25), v0
        complex(dp) :: g0(size(plan%original_shape)), h0(size(plan%original_shape),size(plan%original_shape))
        real(dp) :: factors(size(plan%original_shape)), first(size(plan%original_shape))
        real(dp) :: second(size(plan%original_shape)), factor, w(6), dw(6), ddw(6)
        real(dp) :: other
        integer :: nr,d,e,t,k,code,idx(6),nterms

        nr=size(plan%original_shape)
        shifted=coords
        do d=1,nr
            shifted(d,:)=shifted(d,:)+real(plan%halo(d),dp)
        end do
        if (plan%boundary/='grid-constant' .or. .not.any(plan%original_shape==1)) then
            if (plan%boundary=='grid-constant') then
                call interpolate_with_hessian(plan%coefficients,plan%working_shape, &
                    shifted,values,grad,hess,plan%degree,trim(plan%boundary), &
                    (0.0_dp,0.0_dp),status)
                if (status==spline_status_success) values=values+plan%fill
            else
                call interpolate_with_hessian(plan%coefficients,plan%working_shape, &
                    shifted,values,grad,hess,plan%degree,trim(plan%boundary),plan%fill,status)
            end if
            return
        end if
        impulse=(0.0_dp,0.0_dp)
        impulse(13)=(1.0_dp,0.0_dp)
        call line_coeff(impulse,plan%degree,impulse_coeff,code,'mirror')
        status=code
        if (code/=spline_status_success) return
        call interpolate_with_hessian(plan%coefficients,plan%working_shape, &
            shifted,values,grad,hess,plan%degree,trim(plan%boundary), &
            (0.0_dp,0.0_dp),status,collapse_singletons=.true.)
        if (status/=spline_status_success) return
        do k=1,size(values)
            factors=1.0_dp
            first=0.0_dp
            second=0.0_dp
            do d=1,nr
                if (plan%original_shape(d)/=1) cycle
                call support(coords(d,k)+12.0_dp,25,plan%degree, &
                    trim(plan%boundary),idx,w,nterms,dw,ddw)
                factors(d)=0.0_dp
                do t=1,nterms
                    if (idx(t)==0) cycle
                    factors(d)=factors(d)+w(t)*real(impulse_coeff(idx(t)),dp)
                    first(d)=first(d)+dw(t)*real(impulse_coeff(idx(t)),dp)
                    second(d)=second(d)+ddw(t)*real(impulse_coeff(idx(t)),dp)
                end do
            end do
            v0=values(k)
            g0=grad(:,k)
            h0=hess(:,:,k)
            factor=product(factors)
            values(k)=plan%fill+factor*v0
            do d=1,nr
                other=1.0_dp
                do t=1,nr
                    if (t/=d) other=other*factors(t)
                end do
                grad(d,k)=factor*g0(d)+other*first(d)*v0
                hess(d,d,k)=factor*h0(d,d)+2.0_dp*other*first(d)*g0(d)+ &
                    other*second(d)*v0
                do e=d+1,nr
                    block
                        real(dp) :: except_two
                        except_two=1.0_dp
                        do t=1,nr
                            if (t/=d .and. t/=e) except_two=except_two*factors(t)
                        end do
                        hess(d,e,k)=factor*h0(d,e)+other*first(d)*g0(e)+ &
                            (except_two*factors(d))*first(e)*g0(d)+ &
                            except_two*first(d)*first(e)*v0
                        hess(e,d,k)=hess(d,e,k)
                    end block
                end do
            end do
        end do
    end subroutine evaluate_interpolation_with_hessian

    subroutine evaluate_plan_hessian_core(plan,coords,values,grad,hess,status)
        type(spline_interpolation_plan), intent(in) :: plan !! prepared spline plan
        real(dp), intent(in) :: coords(:,:) !! (input rank, points), one-based
        complex(dp), intent(out) :: values(:) !! interpolated values
        complex(dp), intent(out) :: grad(:,:) !! (input rank, points)
        complex(dp), intent(out) :: hess(:,:,:) !! (input rank, input rank, points)
        integer, intent(out) :: status !! zero on success
        integer :: nr,np

        status=spline_status_invalid_input
        values=(0.0_dp,0.0_dp)
        grad=(0.0_dp,0.0_dp)
        hess=(0.0_dp,0.0_dp)
        if (.not.allocated(plan%coefficients)) return
        nr=size(plan%original_shape)
        np=size(values)
        if (size(coords,1)/=nr .or. size(coords,2)/=np) return
        if (any(shape(grad)/=[nr,np])) return
        if (any(shape(hess)/=[nr,nr,np])) return
        if (int(np,int64)>max_padded_elements/int(nr*nr,int64)) then
            status=spline_nd_status_too_expensive
            return
        end if
        if (plan%filtered) then
            call evaluate_interpolation_with_hessian(plan,coords,values,grad,hess,status)
        else
            call interpolate_with_hessian(plan%coefficients,plan%original_shape, &
                coords,values,grad,hess,plan%degree,trim(plan%boundary),plan%fill,status)
        end if
    end subroutine evaluate_plan_hessian_core

    subroutine evaluate_plan_hessian_real(plan,coords,values,gradient,hessian,status)
        type(spline_interpolation_plan), intent(in) :: plan !! cached real-valued coefficients
        real(dp), intent(in) :: coords(:,:) !! (input rank, samples), one-based
        real(dp), intent(out) :: values(:) !! samples
        real(dp), intent(out) :: gradient(:,:) !! first coordinate derivatives
        real(dp), intent(out) :: hessian(:,:,:) !! second coordinate derivatives
        integer, intent(out) :: status !! zero on success
        complex(dp), allocatable :: v(:),g(:,:),h(:,:,:)
        integer :: nr,np
        status=spline_status_invalid_input
        values=0.0_dp
        gradient=0.0_dp
        hessian=0.0_dp
        if (.not.allocated(plan%coefficients)) return
        if (plan%complex_input) return
        nr=size(plan%original_shape)
        np=size(values)
        if (int(np,int64)>max_padded_elements/int(nr*nr,int64)) then
            status=spline_nd_status_too_expensive
            return
        end if
        allocate(v(np),g(nr,np),h(nr,nr,np))
        call evaluate_plan_hessian_core(plan,coords,v,g,h,status)
        if (status==spline_status_success) then
            values=real(v,dp)
            gradient=real(g,dp)
            hessian=real(h,dp)
        end if
    end subroutine evaluate_plan_hessian_real

    subroutine evaluate_plan_hessian_complex(plan,coords,values,gradient,hessian,status)
        type(spline_interpolation_plan), intent(in) :: plan !! cached complex coefficients
        real(dp), intent(in) :: coords(:,:) !! one-based query locations
        complex(dp), intent(out) :: values(:) !! samples
        complex(dp), intent(out) :: gradient(:,:) !! first coordinate derivatives
        complex(dp), intent(out) :: hessian(:,:,:) !! second coordinate derivatives
        integer, intent(out) :: status !! zero on success
        status=spline_status_invalid_input
        values=(0.0_dp,0.0_dp)
        gradient=(0.0_dp,0.0_dp)
        hessian=(0.0_dp,0.0_dp)
        if (.not.allocated(plan%coefficients)) return
        if (.not.plan%complex_input) return
        call evaluate_plan_hessian_core(plan,coords,values,gradient,hessian,status)
    end subroutine evaluate_plan_hessian_complex

    ! Grid-shaped Hessians have two leading source-coordinate dimensions.
    ! Thus a rank-13 output fits within Fortran's maximum rank 15.
    subroutine evaluate_plan_grid_hessian_real(plan,coords,y,gradient,hessian,status)
        type(spline_interpolation_plan), intent(in) :: plan !! cached real coefficients
        real(dp), intent(in) :: coords(..) !! shape [input rank, shape(y)]
        real(dp), intent(out) :: y(..) !! interpolated values, ranks 1 through 13
        real(dp), intent(out) :: gradient(..) !! shape [input rank, shape(y)]
        real(dp), intent(out) :: hessian(..) !! shape [input rank, input rank, shape(y)]
        integer, intent(out) :: status !! zero on success
        real(dp), allocatable :: packed(:), locations(:,:),v(:),g(:,:),h(:,:,:)
        integer :: nr,np
        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (plan%complex_input) return
        if (rank(y)<1 .or. rank(y)>13) return
        nr=size(plan%original_shape)
        if (rank(coords)/=rank(y)+1 .or. rank(gradient)/=rank(y)+1) return
        if (rank(hessian)/=rank(y)+2) return
        if (any(shape(coords)/=[nr,shape(y)])) return
        if (any(shape(gradient)/=[nr,shape(y)])) return
        if (any(shape(hessian)/=[nr,nr,shape(y)])) return
        if (any(shape(y)<1)) return
        if (size(y,kind=int64)>max_padded_elements/int(nr*nr,int64)) then
            status=spline_nd_status_too_expensive
            return
        end if
        np=size(y)
        call pack_real(coords,packed)
        locations=reshape(packed,[nr,np])
        allocate(v(np),g(nr,np),h(nr,nr,np))
        call evaluate_plan_hessian_real(plan,locations,v,g,h,status)
        if (status/=spline_status_success) return
        call scatter_real(v,y)
        call scatter_real(reshape(g,[size(g)]),gradient)
        call scatter_real(reshape(h,[size(h)]),hessian)
    end subroutine evaluate_plan_grid_hessian_real

    subroutine evaluate_plan_grid_hessian_complex(plan,coords,y,gradient,hessian,status)
        type(spline_interpolation_plan), intent(in) :: plan !! cached complex coefficients
        real(dp), intent(in) :: coords(..) !! shape [input rank, shape(y)]
        complex(dp), intent(out) :: y(..) !! interpolated values
        complex(dp), intent(out) :: gradient(..) !! first derivatives
        complex(dp), intent(out) :: hessian(..) !! second derivatives
        integer, intent(out) :: status !! zero on success
        real(dp), allocatable :: packed(:), locations(:,:)
        complex(dp), allocatable :: v(:),g(:,:),h(:,:,:)
        integer :: nr,np
        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (.not.plan%complex_input) return
        if (rank(y)<1 .or. rank(y)>13) return
        nr=size(plan%original_shape)
        if (rank(coords)/=rank(y)+1 .or. rank(gradient)/=rank(y)+1) return
        if (rank(hessian)/=rank(y)+2) return
        if (any(shape(coords)/=[nr,shape(y)])) return
        if (any(shape(gradient)/=[nr,shape(y)])) return
        if (any(shape(hessian)/=[nr,nr,shape(y)])) return
        if (any(shape(y)<1)) return
        if (size(y,kind=int64)>max_padded_elements/int(nr*nr,int64)) then
            status=spline_nd_status_too_expensive
            return
        end if
        np=size(y)
        call pack_real(coords,packed)
        locations=reshape(packed,[nr,np])
        allocate(v(np),g(nr,np),h(nr,nr,np))
        call evaluate_plan_hessian_complex(plan,locations,v,g,h,status)
        if (status/=spline_status_success) return
        call scatter_complex(v,y)
        call scatter_complex(reshape(g,[size(g)]),gradient)
        call scatter_complex(reshape(h,[size(h)]),hessian)
    end subroutine evaluate_plan_grid_hessian_complex

    ! Dense displacement gradients with a grid-shaped leading source axis.
    ! Callers requiring rank-15 output should use the flat displacement API.
    subroutine evaluate_plan_displacement_grid_gradient_real(plan,field,y,gradient,status, &
                                                             matrix,offset,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! prepared real coefficients
        real(dp), intent(in) :: field(..) !! (input rank, output dimensions)
        real(dp), intent(out) :: y(..) !! warped image, rank 1 through 14
        real(dp), intent(out) :: gradient(..) !! (input rank, output dimensions)
        integer, intent(out) :: status !! zero on success
        real(dp), optional, intent(in) :: matrix(:,:) !! affine base mapping
        real(dp), optional, intent(in) :: offset(:) !! affine base offset
        integer, optional, intent(in) :: batch_size !! block size
        real(dp), allocatable :: displacements(:,:),g(:,:)
        integer :: nr
        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (plan%complex_input) return
        if (rank(y)<1 .or. rank(y)>14) return
        if (rank(gradient)/=rank(y)+1) return
        nr=size(plan%original_shape)
        if (any(shape(gradient)/=[nr,shape(y)])) return
        if (size(y,kind=int64)>max_padded_elements/int(nr,int64)) then
            status=spline_nd_status_too_expensive
            return
        end if
        call flatten_plan_displacements(field,nr,shape(y),displacements,status)
        if (status/=spline_status_success) return
        allocate(g(nr,size(y)))
        call evaluate_plan_displacement_gradient_real(plan,displacements,y,g,status, &
            matrix,offset,batch_size)
        if (status==spline_status_success) call scatter_real(reshape(g,[size(g)]),gradient)
    end subroutine evaluate_plan_displacement_grid_gradient_real

    subroutine evaluate_plan_displacement_grid_gradient_complex(plan,field,y,gradient,status, &
                                                                matrix,offset,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! prepared complex coefficients
        real(dp), intent(in) :: field(..) !! (input rank, output dimensions)
        complex(dp), intent(out) :: y(..) !! warped image, rank 1 through 14
        complex(dp), intent(out) :: gradient(..) !! (input rank, output dimensions)
        integer, intent(out) :: status !! zero on success
        real(dp), optional, intent(in) :: matrix(:,:) !! affine base mapping
        real(dp), optional, intent(in) :: offset(:) !! affine base offset
        integer, optional, intent(in) :: batch_size !! block size
        real(dp), allocatable :: displacements(:,:)
        complex(dp), allocatable :: g(:,:)
        integer :: nr
        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (.not.plan%complex_input) return
        if (rank(y)<1 .or. rank(y)>14) return
        if (rank(gradient)/=rank(y)+1) return
        nr=size(plan%original_shape)
        if (any(shape(gradient)/=[nr,shape(y)])) return
        if (size(y,kind=int64)>max_padded_elements/int(nr,int64)) then
            status=spline_nd_status_too_expensive
            return
        end if
        call flatten_plan_displacements(field,nr,shape(y),displacements,status)
        if (status/=spline_status_success) return
        allocate(g(nr,size(y)))
        call evaluate_plan_displacement_gradient_complex(plan,displacements,y,g,status, &
            matrix,offset,batch_size)
        if (status==spline_status_success) call scatter_complex(reshape(g,[size(g)]),gradient)
    end subroutine evaluate_plan_displacement_grid_gradient_complex

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

    subroutine filter_real_core(x,y,order,axis,status,mode)
        real(dp), intent(in) :: x(..) !! arbitrary-rank sampled array
        real(dp), intent(out) :: y(..) !! coefficient array, same shape
        integer, intent(in) :: order !! B-spline order 0..5
        integer, intent(in) :: axis !! 0 for all axes, otherwise 1..rank
        integer, intent(out) :: status !! zero for success
        character(len=*), optional, intent(in) :: mode !! spline boundary mode
        real(dp), allocatable :: input(:)
        complex(dp), allocatable :: work(:)
        integer, allocatable :: dims(:)

        status=spline_status_invalid_input
        if (rank(x) < 1 .or. rank(x) > 15 .or. rank(y) /= rank(x)) return
        if (any(shape(x) /= shape(y))) return
        if (size(x,kind=int64) > int(huge(1),int64)) return
        if (any(shape(x)<1)) return
        if (axis < 0 .or. axis > rank(x)) return
        if (order < 0 .or. order > 5) return
        if (present(mode)) then
            if (.not.valid_boundary(trim(mode))) return
        end if
        call pack_real(x,input)
        if (size(input) /= size(x)) return
        work=cmplx(input,0.0_dp,dp)
        if (any(.not. ieee_is_finite(real(work,dp))) .or. &
            any(.not. ieee_is_finite(aimag(work)))) return
        dims=shape(x)
        call filter_flat(work,dims,order,axis,status,mode)
        if (status /= spline_status_success) return
        input=real(work,dp)
        call scatter_real(input,y)
    end subroutine filter_real_core

    subroutine filter_all_real(x,y,order,status,mode)
        real(dp), intent(in) :: x(..) !! sampled array, rank 1..15
        real(dp), intent(out) :: y(..) !! same-shape spline coefficients
        integer, optional, intent(in) :: order !! default 3, supported 0..5
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! default mirror
        integer :: degree, code
        degree=3
        if (present(order)) degree=order
        call filter_real_core(x,y,degree,0,code,mode)
        if (present(status)) status=code
    end subroutine filter_all_real

    subroutine filter_axis_real(x,y,axis,order,status,mode)
        real(dp), intent(in) :: x(..) !! sampled array, rank 1..15
        real(dp), intent(out) :: y(..) !! same-shape filtered array
        integer, intent(in) :: axis !! one-based axis to filter
        integer, optional, intent(in) :: order !! default 3
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! default mirror
        integer :: degree, code
        degree=3
        if (present(order)) degree=order
        call filter_real_core(x,y,degree,axis,code,mode)
        if (present(status)) status=code
    end subroutine filter_axis_real

    subroutine cubic_real(x,y,status)
        real(dp), intent(in) :: x(..) !! sampled array, rank 1..15
        real(dp), intent(out) :: y(..) !! same-shape cubic spline coefficients
        integer, optional, intent(out) :: status !! zero on success
        call filter_all_real(x,y,3,status)
    end subroutine cubic_real

    subroutine quadratic_real(x,y,status)
        real(dp), intent(in) :: x(..) !! sampled array, rank 1..15
        real(dp), intent(out) :: y(..) !! same-shape quadratic spline coefficients
        integer, optional, intent(out) :: status !! zero on success
        call filter_all_real(x,y,2,status)
    end subroutine quadratic_real

    subroutine filter_complex_core(x,y,order,axis,status,mode)
        complex(dp), intent(in) :: x(..) !! arbitrary-rank sampled array
        complex(dp), intent(out) :: y(..) !! coefficient array, same shape
        integer, intent(in) :: order !! B-spline order 0..5
        integer, intent(in) :: axis !! 0 for all axes, otherwise 1..rank
        integer, intent(out) :: status !! zero for success
        character(len=*), optional, intent(in) :: mode !! spline boundary mode
        complex(dp), allocatable :: input(:)
        complex(dp), allocatable :: work(:)
        integer, allocatable :: dims(:)

        status=spline_status_invalid_input
        if (rank(x) < 1 .or. rank(x) > 15 .or. rank(y) /= rank(x)) return
        if (any(shape(x) /= shape(y))) return
        if (size(x,kind=int64) > int(huge(1),int64)) return
        if (any(shape(x)<1)) return
        if (axis < 0 .or. axis > rank(x)) return
        if (order < 0 .or. order > 5) return
        if (present(mode)) then
            if (.not.valid_boundary(trim(mode))) return
        end if
        call pack_complex(x,input)
        if (size(input) /= size(x)) return
        work=input
        if (any(.not. ieee_is_finite(real(work,dp))) .or. &
            any(.not. ieee_is_finite(aimag(work)))) return
        dims=shape(x)
        call filter_flat(work,dims,order,axis,status,mode)
        if (status /= spline_status_success) return
        input=work
        call scatter_complex(input,y)
    end subroutine filter_complex_core

    subroutine filter_all_complex(x,y,order,status,mode)
        complex(dp), intent(in) :: x(..) !! sampled array, rank 1..15
        complex(dp), intent(out) :: y(..) !! same-shape spline coefficients
        integer, optional, intent(in) :: order !! default 3, supported 0..5
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! default mirror
        integer :: degree, code
        degree=3
        if (present(order)) degree=order
        call filter_complex_core(x,y,degree,0,code,mode)
        if (present(status)) status=code
    end subroutine filter_all_complex

    subroutine filter_axis_complex(x,y,axis,order,status,mode)
        complex(dp), intent(in) :: x(..) !! sampled array, rank 1..15
        complex(dp), intent(out) :: y(..) !! same-shape filtered array
        integer, intent(in) :: axis !! one-based axis to filter
        integer, optional, intent(in) :: order !! default 3
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! default mirror
        integer :: degree, code
        degree=3
        if (present(order)) degree=order
        call filter_complex_core(x,y,degree,axis,code,mode)
        if (present(status)) status=code
    end subroutine filter_axis_complex

    subroutine cubic_complex(x,y,status)
        complex(dp), intent(in) :: x(..) !! sampled array, rank 1..15
        complex(dp), intent(out) :: y(..) !! same-shape cubic spline coefficients
        integer, optional, intent(out) :: status !! zero on success
        call filter_all_complex(x,y,3,status)
    end subroutine cubic_complex

    subroutine quadratic_complex(x,y,status)
        complex(dp), intent(in) :: x(..) !! sampled array, rank 1..15
        complex(dp), intent(out) :: y(..) !! same-shape quadratic spline coefficients
        integer, optional, intent(out) :: status !! zero on success
        call filter_all_complex(x,y,2,status)
    end subroutine quadratic_complex

    ! The plan owns a copy of its input (or its prefiltered coefficients),
    ! so it remains usable after the caller modifies or frees the input.
    subroutine prepare_plan_core(samples,dims,plan,degree,boundary,filtered,fill,status,is_complex)
        complex(dp), intent(in) :: samples(:) !! packed samples or coefficients
        integer, intent(in) :: dims(:) !! shape of source image
        type(spline_interpolation_plan), intent(inout) :: plan !! prepared cache
        integer, intent(in) :: degree !! interpolation order, 0 through 5
        character(len=*), intent(in) :: boundary !! ndimage boundary extension
        logical, intent(in) :: filtered !! prefilter the source once
        complex(dp), intent(in) :: fill !! out-of-bounds constant
        integer, intent(out) :: status !! zero for success
        logical, intent(in) :: is_complex !! identify complex source type
        complex(dp), allocatable :: result(:)
        integer, allocatable :: expanded(:), padding(:)
        integer :: rank_input

        ! Reset before validation: an invalid rebuild never retains old data.
        plan=spline_interpolation_plan()
        status=spline_status_invalid_input
        rank_input=size(dims)
        if (rank_input<1 .or. rank_input>15) return
        if (degree<0 .or. degree>5) return
        if (.not.valid_boundary(boundary)) return
        if (any(dims<1)) return
        if (size(samples,kind=int64)>int(huge(1),int64)) return
        if (any(.not.ieee_is_finite(real(samples,dp))) .or. &
            any(.not.ieee_is_finite(aimag(samples)))) return
        if (.not.ieee_is_finite(real(fill,dp)) .or. &
            .not.ieee_is_finite(aimag(fill))) return

        expanded=dims
        allocate(padding(rank_input))
        padding=0
        if (filtered .and. degree>1) then
            call prepare_interpolation_coeff(samples,dims,degree,boundary,fill, &
                result,expanded,padding,status)
            if (status/=spline_status_success) return
        else
            result=samples
        end if
        call move_alloc(result,plan%coefficients)
        plan%original_shape=dims
        plan%working_shape=expanded
        plan%halo=padding
        plan%degree=degree
        plan%boundary=boundary
        plan%filtered=filtered .and. degree>1
        plan%fill=fill
        plan%complex_input=is_complex
        status=spline_status_success
    end subroutine prepare_plan_core

    subroutine prepare_plan_real(x,plan,order,status,mode,prefilter,cval)
        real(dp), intent(in) :: x(..) !! real samples or prefiltered coefficients
        type(spline_interpolation_plan), intent(inout) :: plan !! reusable spline cache
        integer, optional, intent(in) :: order !! spline order, default 3
        integer, intent(out) :: status !! zero for success
        character(len=*), optional, intent(in) :: mode !! extension, default mirror
        logical, optional, intent(in) :: prefilter !! filter samples once, default true
        real(dp), optional, intent(in) :: cval !! constant exterior, default zero
        real(dp), allocatable :: packed(:)
        integer :: degree
        character(len=:), allocatable :: boundary
        logical :: filtered
        complex(dp) :: fill

        plan=spline_interpolation_plan()
        status=spline_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15) return
        if (any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        degree=3
        if (present(order)) degree=order
        boundary='mirror'
        if (present(mode)) boundary=trim(mode)
        filtered=.true.
        if (present(prefilter)) filtered=prefilter
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cmplx(cval,0.0_dp,dp)
        call pack_real(x,packed)
        call prepare_plan_core(cmplx(packed,0.0_dp,dp),shape(x),plan, &
            degree,boundary,filtered,fill,status,.false.)
    end subroutine prepare_plan_real

    subroutine prepare_plan_complex(x,plan,order,status,mode,prefilter,cval)
        complex(dp), intent(in) :: x(..) !! complex samples or filtered coefficients
        type(spline_interpolation_plan), intent(inout) :: plan !! reusable spline cache
        integer, optional, intent(in) :: order !! spline order, default 3
        integer, intent(out) :: status !! zero for success
        character(len=*), optional, intent(in) :: mode !! extension, default mirror
        logical, optional, intent(in) :: prefilter !! filter samples once, default true
        complex(dp), optional, intent(in) :: cval !! constant exterior, default zero
        complex(dp), allocatable :: packed(:)
        integer :: degree
        character(len=:), allocatable :: boundary
        logical :: filtered
        complex(dp) :: fill

        plan=spline_interpolation_plan()
        status=spline_status_invalid_input
        if (rank(x)<1 .or. rank(x)>15) return
        if (any(shape(x)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        degree=3
        if (present(order)) degree=order
        boundary='mirror'
        if (present(mode)) boundary=trim(mode)
        filtered=.true.
        if (present(prefilter)) filtered=prefilter
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cval
        call pack_complex(x,packed)
        call prepare_plan_core(packed,shape(x),plan,degree,boundary,filtered, &
            fill,status,.true.)
    end subroutine prepare_plan_complex

    subroutine evaluate_plan_core(plan,coords,values,status)
        type(spline_interpolation_plan), intent(in) :: plan !! prepared cache
        real(dp), intent(in) :: coords(:,:) !! one-based (rank, npoints)
        complex(dp), intent(out) :: values(:) !! complex-valued query results
        integer, intent(out) :: status !! zero for success

        status=spline_status_invalid_input
        values=(0.0_dp,0.0_dp)
        if (.not.allocated(plan%coefficients)) return
        if (size(coords,1)/=size(plan%original_shape)) return
        if (size(coords,2)/=size(values)) return
        if (plan%filtered) then
            call evaluate_interpolation(plan%coefficients,plan%working_shape, &
                plan%halo,plan%original_shape,coords,plan%degree, &
                trim(plan%boundary),plan%fill,values,status)
        else
            call interpolate(plan%coefficients,plan%original_shape,coords, &
                values,plan%degree,trim(plan%boundary),plan%fill,status)
        end if
    end subroutine evaluate_plan_core

    subroutine evaluate_plan_real(plan,coords,values,status)
        type(spline_interpolation_plan), intent(in) :: plan !! prepared real cache
        real(dp), intent(in) :: coords(:,:) !! one-based positions
        real(dp), intent(out) :: values(:) !! interpolated real samples
        integer, intent(out) :: status !! zero for success
        complex(dp), allocatable :: result(:)

        values=0.0_dp
        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (plan%complex_input) return
        allocate(result(size(values)))
        call evaluate_plan_core(plan,coords,result,status)
        if (status==spline_status_success) values=real(result,dp)
    end subroutine evaluate_plan_real

    subroutine evaluate_plan_complex(plan,coords,values,status)
        type(spline_interpolation_plan), intent(in) :: plan !! prepared complex cache
        real(dp), intent(in) :: coords(:,:) !! one-based positions
        complex(dp), intent(out) :: values(:) !! interpolated complex samples
        integer, intent(out) :: status !! zero for success

        values=(0.0_dp,0.0_dp)
        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (.not.plan%complex_input) return
        call evaluate_plan_core(plan,coords,values,status)
    end subroutine evaluate_plan_complex

    ! Shape-preserving plan evaluation. Coordinate leading axis specifies
    ! source dimensions and trailing axes describe the destination shape.
    subroutine evaluate_plan_grid_real(plan,coords,y,status)
        type(spline_interpolation_plan), intent(in) :: plan !! prepared real cache
        real(dp), intent(in) :: coords(..) !! shape (rank(source), shape(y))
        real(dp), intent(out) :: y(..) !! output rank 1 through 14
        integer, intent(out) :: status !! zero for success
        real(dp), allocatable :: packed(:), locations(:,:), values(:)
        integer, allocatable :: dims(:)

        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (plan%complex_input) return
        if (rank(coords)<2 .or. rank(coords)>15) return
        if (rank(y)/=rank(coords)-1) return
        dims=shape(coords)
        if (dims(1)/=size(plan%original_shape)) return
        if (any(dims(2:)/=shape(y))) return
        if (any(shape(y)<1)) return
        if (size(y,kind=int64)>int(huge(1),int64)) return
        call pack_real(coords,packed)
        locations=reshape(packed,[dims(1),size(y)])
        allocate(values(size(y)))
        call evaluate_plan_real(plan,locations,values,status)
        if (status==spline_status_success) call scatter_real(values,y)
    end subroutine evaluate_plan_grid_real

    subroutine evaluate_plan_grid_complex(plan,coords,y,status)
        type(spline_interpolation_plan), intent(in) :: plan !! prepared complex cache
        real(dp), intent(in) :: coords(..) !! shape (rank(source), shape(y))
        complex(dp), intent(out) :: y(..) !! output rank 1 through 14
        integer, intent(out) :: status !! zero for success
        real(dp), allocatable :: packed(:), locations(:,:)
        complex(dp), allocatable :: values(:)
        integer, allocatable :: dims(:)

        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (.not.plan%complex_input) return
        if (rank(coords)<2 .or. rank(coords)>15) return
        if (rank(y)/=rank(coords)-1) return
        dims=shape(coords)
        if (dims(1)/=size(plan%original_shape)) return
        if (any(dims(2:)/=shape(y))) return
        if (any(shape(y)<1)) return
        if (size(y,kind=int64)>int(huge(1),int64)) return
        call pack_real(coords,packed)
        locations=reshape(packed,[dims(1),size(y)])
        allocate(values(size(y)))
        call evaluate_plan_complex(plan,locations,values,status)
        if (status==spline_status_success) call scatter_complex(values,y)
    end subroutine evaluate_plan_grid_complex

    ! Grid-shaped analytic coordinate derivatives. The leading gradient axis
    ! is the input rank; the remaining axes match the image destination.
    ! Coordinates and destination may be noncontiguous; output rank <=14.
    subroutine evaluate_plan_grid_gradient_real(plan,coords,y,gradient,status)
        type(spline_interpolation_plan), intent(in) :: plan !! prepared real spline plan
        real(dp), intent(in) :: coords(..) !! shape [input rank, shape(y)]
        real(dp), intent(out) :: y(..) !! shape-preserving interpolated values
        real(dp), intent(out) :: gradient(..) !! [input rank, shape(y)] derivatives
        integer, intent(out) :: status !! success or input/resource error
        real(dp), allocatable :: packed(:),locations(:,:),v(:),g(:,:)
        integer, allocatable :: dims(:)
        integer :: nr

        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (plan%complex_input) return
        if (rank(y)<1 .or. rank(y)>14) return
        if (rank(coords)/=rank(y)+1 .or. rank(gradient)/=rank(y)+1) return
        nr=size(plan%original_shape)
        if (any(shape(coords)/=[nr,shape(y)])) return
        if (any(shape(gradient)/=[nr,shape(y)])) return
        if (any(shape(y)<1)) return
        if (size(y,kind=int64)>max_padded_elements/int(nr,int64)) then
            status=spline_nd_status_too_expensive
            return
        end if
        dims=shape(coords)
        call pack_real(coords,packed)
        locations=reshape(packed,[nr,size(y)])
        allocate(v(size(y)),g(nr,size(y)))
        call evaluate_plan_gradient_real(plan,locations,v,g,status)
        if (status/=spline_status_success) return
        call scatter_real(v,y)
        call scatter_real(reshape(g,[size(g)]),gradient)
    end subroutine evaluate_plan_grid_gradient_real

    subroutine evaluate_plan_grid_gradient_complex(plan,coords,y,gradient,status)
        type(spline_interpolation_plan), intent(in) :: plan !! prepared complex spline plan
        real(dp), intent(in) :: coords(..) !! shape [input rank, shape(y)]
        complex(dp), intent(out) :: y(..) !! shape-preserving interpolated values
        complex(dp), intent(out) :: gradient(..) !! [input rank, shape(y)] derivatives
        integer, intent(out) :: status !! success or input/resource error
        real(dp), allocatable :: packed(:),locations(:,:)
        complex(dp), allocatable :: v(:),g(:,:)
        integer :: nr

        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (.not.plan%complex_input) return
        if (rank(y)<1 .or. rank(y)>14) return
        if (rank(coords)/=rank(y)+1 .or. rank(gradient)/=rank(y)+1) return
        nr=size(plan%original_shape)
        if (any(shape(coords)/=[nr,shape(y)])) return
        if (any(shape(gradient)/=[nr,shape(y)])) return
        if (any(shape(y)<1)) return
        if (size(y,kind=int64)>max_padded_elements/int(nr,int64)) then
            status=spline_nd_status_too_expensive
            return
        end if
        call pack_real(coords,packed)
        locations=reshape(packed,[nr,size(y)])
        allocate(v(size(y)),g(nr,size(y)))
        call evaluate_plan_gradient_complex(plan,locations,v,g,status)
        if (status/=spline_status_success) return
        call scatter_complex(v,y)
        call scatter_complex(reshape(g,[size(g)]),gradient)
    end subroutine evaluate_plan_grid_gradient_complex

    ! Tensor-product grid evaluation performs one-dimensional contractions
    ! instead of repeating a rank-dimensional kernel at every output pixel.
    ! Working buffers are bounded independently of the kernel rank.
    subroutine evaluate_plan_tensor_core(plan,grid,values,status)
        type(spline_interpolation_plan), intent(in) :: plan !! prepared source and coefficients
        type(spline_axis_grid), intent(in) :: grid(:) !! one coordinate vector per source axis
        complex(dp), intent(out) :: values(:) !! packed output in Fortran array order
        integer, intent(out) :: status !! spline status
        integer(int64), parameter :: max_tensor_elements=8000000_int64
        complex(dp), allocatable :: work(:), next(:)
        complex(dp) :: impulse(25), impulse_coeff(25), accum
        integer, allocatable :: dims(:), outdims(:), ix(:,:), counts(:)
        real(dp), allocatable :: weights(:,:)
        integer :: nd, d, j, k, code, n, m, i0(6), nterms, t
        integer(int64) :: block, outside, oldbase, newbase, pos, total
        real(dp) :: w0(6), shift, original_x, factor
        logical :: use_impulse, is_constant, is_grid_constant, out_of_bounds

        status=spline_status_invalid_input
        values=(0.0_dp,0.0_dp)
        if (.not.allocated(plan%coefficients)) return
        nd=size(plan%original_shape)
        if (nd<1 .or. nd>15 .or. size(grid)/=nd) return
        if (size(values,kind=int64)>max_tensor_elements) then
            status=spline_nd_status_too_expensive
            return
        end if
        allocate(outdims(nd))
        do d=1,nd
            if (.not.allocated(grid(d)%coordinates)) return
            m=size(grid(d)%coordinates)
            if (m<1) return
            if (.not.all(ieee_is_finite(grid(d)%coordinates))) return
            outdims(d)=m
        end do
        total=1_int64
        do d=1,nd
            if (total>max_tensor_elements/int(outdims(d),int64)) then
                status=spline_nd_status_too_expensive
                return
            end if
            total=total*int(outdims(d),int64)
        end do
        if (total/=size(values,kind=int64)) return
        if (size(plan%coefficients,kind=int64)>max_tensor_elements) then
            status=spline_nd_status_too_expensive
            return
        end if

        is_constant=trim(plan%boundary)=='constant'
        is_grid_constant=trim(plan%boundary)=='grid-constant'
        use_impulse=is_grid_constant .and. plan%filtered .and. &
                    any(plan%original_shape==1)
        if (use_impulse) then
            impulse=(0.0_dp,0.0_dp)
            impulse(13)=(1.0_dp,0.0_dp)
            call line_coeff(impulse,plan%degree,impulse_coeff,code,'mirror')
            if (code/=spline_status_success) then
                status=code
                return
            end if
        end if
        work=plan%coefficients
        dims=plan%working_shape
        ! A constant contribution can be evaluated separately: zero is the
        ! correct extension of (samples - cval) for grid-constant mode.
        if (is_grid_constant .and. .not.plan%filtered) work=work-plan%fill
        do d=1,nd
            n=dims(d)
            m=outdims(d)
            allocate(ix(6,m),weights(6,m),counts(m))
            ix=0
            weights=0.0_dp
            counts=0
            do j=1,m
                original_x=grid(d)%coordinates(j)
                if (is_constant) then
                    if (original_x<1.0_dp .or. &
                        original_x>real(plan%original_shape(d),dp)) cycle
                end if
                if (use_impulse .and. plan%original_shape(d)==1) then
                    ! Nearest padded singleton axes are constant, but a
                    ! grid-constant singleton is a filtered unit impulse.
                    call support(original_x+12.0_dp,25,plan%degree, &
                        'grid-constant',i0,w0,nterms)
                    factor=0.0_dp
                    do k=1,nterms
                        if (i0(k)>0) factor=factor+w0(k)* &
                            real(impulse_coeff(i0(k)),dp)
                    end do
                    counts(j)=1
                    ix(1,j)=1
                    weights(1,j)=factor
                else
                    shift=real(plan%halo(d),dp)
                    call support(original_x+shift,n,plan%degree, &
                        trim(plan%boundary),ix(:,j),weights(:,j),counts(j))
                end if
            end do
            block=1_int64
            do t=1,d-1
                block=block*int(dims(t),int64)
            end do
            outside=int(size(work,kind=int64)/(block*int(n,int64)))
            total=outside*block*int(m,int64)
            if (total>max_tensor_elements) then
                status=spline_nd_status_too_expensive
                return
            end if
            allocate(next(total))
            do pos=0_int64,outside-1_int64
                oldbase=pos*block*int(n,int64)
                newbase=pos*block*int(m,int64)
                do j=1,m
                    do t=0,int(block)-1
                        accum=(0.0_dp,0.0_dp)
                        do k=1,counts(j)
                            if (ix(k,j)<=0) cycle
                            accum=accum+weights(k,j)* &
                                work(oldbase+int(ix(k,j)-1,int64)*block+int(t,int64)+1_int64)
                        end do
                        next(newbase+int(j-1,int64)*block+int(t,int64)+1_int64)=accum
                    end do
                end do
            end do
            call move_alloc(next,work)
            dims(d)=m
            deallocate(ix,weights,counts)
        end do
        if (is_grid_constant) work=work+plan%fill
        if (is_constant) then
            ! With constant boundaries, any outside coordinate sets the
            ! whole output pixel to cval, not one fill per outside axis.
            do pos=0_int64,size(work,kind=int64)-1_int64
                total=pos
                out_of_bounds=.false.
                do d=1,nd
                    j=int(modulo(total,int(outdims(d),int64)))+1
                    total=total/int(outdims(d),int64)
                    if (grid(d)%coordinates(j)<1.0_dp .or. &
                        grid(d)%coordinates(j)>real(plan%original_shape(d),dp)) then
                        out_of_bounds=.true.
                        exit
                    end if
                end do
                if (out_of_bounds) work(pos+1_int64)=plan%fill
            end do
        end if
        values=work
        status=spline_status_success
    end subroutine evaluate_plan_tensor_core

    subroutine evaluate_plan_tensor_real(plan,grid,y,status)
        type(spline_interpolation_plan), intent(in) :: plan !! prepared real spline source
        type(spline_axis_grid), intent(in) :: grid(:) !! one-based axis coordinate vectors
        real(dp), intent(out) :: y(..) !! tensor product output, rank matches source
        integer, intent(out) :: status !! spline status
        complex(dp), allocatable :: values(:)
        integer :: d

        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (plan%complex_input .or. rank(y)/=size(plan%original_shape)) return
        if (size(grid)/=rank(y)) return
        if (size(y,kind=int64)>8000000_int64) then
            status=spline_nd_status_too_expensive
            return
        end if
        do d=1,size(grid)
            if (.not.allocated(grid(d)%coordinates)) return
            if (size(grid(d)%coordinates)/=size(y,dim=d)) return
        end do
        allocate(values(size(y)))
        call evaluate_plan_tensor_core(plan,grid,values,status)
        if (status==spline_status_success) call scatter_real(real(values,dp),y)
    end subroutine evaluate_plan_tensor_real

    subroutine evaluate_plan_tensor_complex(plan,grid,y,status)
        type(spline_interpolation_plan), intent(in) :: plan !! prepared complex spline source
        type(spline_axis_grid), intent(in) :: grid(:) !! one-based axis coordinate vectors
        complex(dp), intent(out) :: y(..) !! tensor product output, rank matches source
        integer, intent(out) :: status !! spline status
        complex(dp), allocatable :: values(:)
        integer :: d

        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (.not.plan%complex_input .or. rank(y)/=size(plan%original_shape)) return
        if (size(grid)/=rank(y)) return
        if (size(y,kind=int64)>8000000_int64) then
            status=spline_nd_status_too_expensive
            return
        end if
        do d=1,size(grid)
            if (.not.allocated(grid(d)%coordinates)) return
            if (size(grid(d)%coordinates)/=size(y,dim=d)) return
        end do
        allocate(values(size(y)))
        call evaluate_plan_tensor_core(plan,grid,values,status)
        if (status==spline_status_success) call scatter_complex(values,y)
    end subroutine evaluate_plan_tensor_complex

    ! Evaluate arbitrary nonseparable output mappings using a fixed-size
    ! coordinate batch. This avoids allocating a (source rank, output size)
    ! coordinate table; the source plan is not modified by any evaluation.
    ! Affine offsets and matrices follow SciPy's zero-based convention.
    ! Callback coordinates, on both sides, use Fortran one-based indexing.
    subroutine evaluate_plan_transform_core(plan,out_dims,values,status,matrix,offset,mapping,batch_size,batch_mapping, &
        displacements,gradient,dmatrix,doffset,hessian,hmm,hmo,hoo)
        type(spline_interpolation_plan), intent(in) :: plan !! prepared source coefficients
        integer, intent(in) :: out_dims(:) !! destination extents
        complex(dp), intent(out) :: values(:) !! packed destination samples
        integer, intent(out) :: status !! zero on success
        real(dp), optional, intent(in) :: matrix(:,:) !! source-by-output affine matrix
        real(dp), optional, intent(in) :: offset(:) !! zero-based source offset
        procedure(ndimage_coordinate_map), optional :: mapping !! one-based scalar callback
        integer, optional, intent(in) :: batch_size !! maximum positions per interpolation block
        procedure(ndimage_coordinate_batch_map), optional :: batch_mapping !! vectorized coordinate callback
        real(dp), optional, intent(in) :: displacements(:,:) !! (input rank, flat output pixels), additive offsets
        complex(dp), optional, intent(out) :: gradient(:,:) !! (input rank, flat output), dI/d displacement
        complex(dp), optional, intent(out) :: dmatrix(:,:,:) !! (input rank, output rank, pixels)
        complex(dp), optional, intent(out) :: doffset(:,:) !! (input rank, pixels)
        complex(dp), optional, intent(out) :: hessian(:,:,:) !! displacement Hessian (input rank, input rank, pixels)
        complex(dp), optional, intent(out) :: hmm(:,:,:,:,:) !! affine matrix-matrix second derivatives
        complex(dp), optional, intent(out) :: hmo(:,:,:,:) !! affine matrix-offset second derivatives
        complex(dp), optional, intent(out) :: hoo(:,:,:) !! affine offset-offset second derivatives
        real(dp), allocatable :: coords(:,:)
        integer, allocatable :: output_block(:,:)
        complex(dp), allocatable :: block_values(:),block_gradient(:,:),block_hessian(:,:,:)
        integer :: nr,nd,i,d,n_here,capacity,output_coords(size(out_dims)),code
        integer(int64) :: base,linear,product_dims

        status=spline_status_invalid_input
        values=(0.0_dp,0.0_dp)
        if (present(gradient)) gradient=(0.0_dp,0.0_dp)
        if (present(dmatrix)) dmatrix=(0.0_dp,0.0_dp)
        if (present(doffset)) doffset=(0.0_dp,0.0_dp)
        if (present(hessian)) hessian=(0.0_dp,0.0_dp)
        if (present(hmm)) hmm=(0.0_dp,0.0_dp)
        if (present(hmo)) hmo=(0.0_dp,0.0_dp)
        if (present(hoo)) hoo=(0.0_dp,0.0_dp)
        if (present(hmm) .neqv. present(hmo)) return
        if (present(hmm) .neqv. present(hoo)) return
        if (present(hmm) .and. .not.present(dmatrix)) return
        if (present(dmatrix) .neqv. present(doffset)) return
        if (.not.allocated(plan%coefficients)) return
        nr=size(plan%original_shape)
        nd=size(out_dims)
        if (nr<1 .or. nr>15 .or. nd<1 .or. nd>15) return
        if (any(out_dims<1)) return
        capacity=256
        if (present(batch_size)) capacity=batch_size
        if (capacity<1 .or. capacity>8192) return
        if (present(mapping) .or. present(batch_mapping)) then
            if (present(mapping) .and. present(batch_mapping)) return
            if (present(matrix) .or. present(offset) .or. present(displacements)) return
        else
            if (present(matrix) .neqv. present(offset)) return
            if (.not.present(displacements) .and. .not.present(matrix)) return
            if (present(matrix)) then
                if (size(matrix,1)/=nr .or. size(matrix,2)/=nd) return
                if (size(offset)/=nr) return
                if (.not.all(ieee_is_finite(matrix))) return
                if (.not.all(ieee_is_finite(offset))) return
            else
                ! Without an affine base, the displacement is relative to
                ! the same one-based output coordinate in every dimension.
                if (nr/=nd) return
            end if
        end if
        product_dims=1_int64
        do d=1,nd
            if (product_dims>max_padded_elements/int(out_dims(d),int64)) then
                status=spline_nd_status_too_expensive
                return
            end if
            product_dims=product_dims*int(out_dims(d),int64)
        end do
        if (product_dims/=size(values,kind=int64)) return
        if (present(dmatrix)) then
            if (.not.present(matrix) .or. .not.present(offset)) return
            if (size(dmatrix,1)/=nr .or. size(dmatrix,2)/=nd) return
            if (size(dmatrix,3,kind=int64)/=product_dims) return
            if (size(doffset,1)/=nr .or. size(doffset,2,kind=int64)/=product_dims) return
            if (product_dims>max_padded_elements/int(nr*nd,int64)) then
                status=spline_nd_status_too_expensive
                return
            end if
        end if
        if (present(gradient)) then
            if (size(gradient,1)/=nr .or. size(gradient,2,kind=int64)/=product_dims) return
            if (product_dims>max_padded_elements/int(nr,int64)) then
                status=spline_nd_status_too_expensive
                return
            end if
        end if
        if (present(hessian)) then
            if (any(shape(hessian)/=[nr,nr,int(product_dims)])) return
            if (product_dims>max_padded_elements/int(nr*nr,int64)) then
                status=spline_nd_status_too_expensive
                return
            end if
        end if
        if (present(hmm)) then
            if (any(shape(hmm)/=[nr,nd,nr,nd,int(product_dims)])) return
            if (any(shape(hmo)/=[nr,nd,nr,int(product_dims)])) return
            if (any(shape(hoo)/=[nr,nr,int(product_dims)])) return
            if (product_dims>max_padded_elements/int(nr*nr*(nd+1)*(nd+1),int64)) then
                status=spline_nd_status_too_expensive
                return
            end if
        end if
        if (present(displacements)) then
            if (size(displacements,1)/=nr) return
            if (size(displacements,2,kind=int64)/=product_dims) return
            if (.not.all(ieee_is_finite(displacements))) return
        end if
        allocate(coords(nr,capacity),block_values(capacity))
        if (present(gradient) .or. present(dmatrix) .or. present(hessian)) &
            allocate(block_gradient(nr,capacity))
        if (present(hessian) .or. present(hmm)) allocate(block_hessian(nr,nr,capacity))
        if (present(batch_mapping) .or. present(dmatrix)) allocate(output_block(nd,capacity))
        base=0_int64
        do while (base<product_dims)
            n_here=int(min(int(capacity,int64),product_dims-base))
            do i=1,n_here
                linear=base+int(i-1,int64)
                do d=1,nd
                    output_coords(d)=int(modulo(linear,int(out_dims(d),int64)))+1
                    linear=linear/int(out_dims(d),int64)
                end do
                if (present(batch_mapping) .or. present(dmatrix)) output_block(:,i)=output_coords
                if (present(batch_mapping)) then
                    ! The callback fills source coordinates after this loop.
                else if (present(mapping)) then
                    call mapping(output_coords,coords(:,i))
                else if (present(matrix)) then
                    coords(:,i)=1.0_dp+offset+matmul(matrix,real(output_coords-1,dp))
                else
                    coords(:,i)=real(output_coords,dp)
                end if
                if (present(displacements)) then
                    coords(:,i)=coords(:,i)+displacements(:,base+int(i,int64))
                end if
            end do
            if (present(batch_mapping)) then
                call batch_mapping(output_block(:,:n_here),coords(:,:n_here))
            end if
            if (.not.all(ieee_is_finite(coords(:,:n_here)))) return
            if (present(hessian) .or. present(hmm)) then
                call evaluate_plan_hessian_core(plan,coords(:,:n_here), &
                    block_values(:n_here),block_gradient(:,:n_here),block_hessian(:,:,:n_here),code)
            else if (present(gradient) .or. present(dmatrix)) then
                call evaluate_plan_gradient_core(plan,coords(:,:n_here), &
                    block_values(:n_here),block_gradient(:,:n_here),code)
            else
                call evaluate_plan_core(plan,coords(:,:n_here),block_values(:n_here),code)
            end if
            if (code/=spline_status_success) then
                status=code
                return
            end if
            values(base+1_int64:base+int(n_here,int64))=block_values(:n_here)
            if (present(gradient)) then
                gradient(:,base+1_int64:base+int(n_here,int64))=block_gradient(:,:n_here)
            end if
            if (present(hessian)) then
                hessian(:,:,base+1_int64:base+int(n_here,int64))=block_hessian(:,:,:n_here)
            end if
            if (present(hmm)) then
                block
                    integer :: a,b,c,e
                    real(dp) :: xj,xl
                    do i=1,n_here
                        do a=1,nr
                            do b=1,nr
                                hoo(a,b,base+int(i,int64))=block_hessian(a,b,i)
                            end do
                        end do
                        do c=1,nd
                            xj=real(output_block(c,i)-1,dp)
                            do a=1,nr
                                do b=1,nr
                                    hmo(a,c,b,base+int(i,int64))=block_hessian(a,b,i)*xj
                                end do
                            end do
                            do e=1,nd
                                xl=real(output_block(e,i)-1,dp)
                                do a=1,nr
                                    do b=1,nr
                                        hmm(a,c,b,e,base+int(i,int64))=block_hessian(a,b,i)*xj*xl
                                    end do
                                end do
                            end do
                        end do
                    end do
                end block
            end if
            if (present(dmatrix)) then
                do i=1,n_here
                    do d=1,nd
                        dmatrix(:,d,base+int(i,int64))=block_gradient(:,i)* &
                            real(output_block(d,i)-1,dp)
                    end do
                    doffset(:,base+int(i,int64))=block_gradient(:,i)
                end do
            end if
            base=base+int(n_here,int64)
        end do
        status=spline_status_success
    end subroutine evaluate_plan_transform_core

    ! Dense displacement fields use shape [source_rank, shape(output)].
    ! A rank-15 output instead uses the flat field interface, because
    ! Fortran permits no rank-16 arrays. Source coordinates are one-based;
    ! displacements themselves are additive pixel offsets.
    subroutine flatten_plan_displacements(field,nr,out_dims,flattened,status)
        real(dp), intent(in) :: field(..) !! real displacement grid, rank 2..15
        integer, intent(in) :: nr !! source rank
        integer, intent(in) :: out_dims(:) !! destination dimensions
        real(dp), allocatable, intent(out) :: flattened(:,:) !! (nr, product(out_dims))
        integer, intent(out) :: status !! zero on success
        integer :: npixels,dd,d
        integer(int64) :: count

        status=spline_status_invalid_input
        dd=size(out_dims)
        if (dd<1 .or. dd>14) return
        if (rank(field)/=dd+1) return
        if (any(shape(field)/=[nr,out_dims])) return
        if (any(out_dims<1)) return
        count=1_int64
        do d=1,dd
            if (count>max_padded_elements/int(out_dims(d),int64)) then
                status=spline_nd_status_too_expensive
                return
            end if
            count=count*int(out_dims(d),int64)
        end do
        npixels=int(count)
        select rank (field)
        rank (2)
            flattened=reshape(field,[nr,npixels])
        rank (3)
            flattened=reshape(field,[nr,npixels])
        rank (4)
            flattened=reshape(field,[nr,npixels])
        rank (5)
            flattened=reshape(field,[nr,npixels])
        rank (6)
            flattened=reshape(field,[nr,npixels])
        rank (7)
            flattened=reshape(field,[nr,npixels])
        rank (8)
            flattened=reshape(field,[nr,npixels])
        rank (9)
            flattened=reshape(field,[nr,npixels])
        rank (10)
            flattened=reshape(field,[nr,npixels])
        rank (11)
            flattened=reshape(field,[nr,npixels])
        rank (12)
            flattened=reshape(field,[nr,npixels])
        rank (13)
            flattened=reshape(field,[nr,npixels])
        rank (14)
            flattened=reshape(field,[nr,npixels])
        rank (15)
            flattened=reshape(field,[nr,npixels])
        rank default
            return
        end select
        status=spline_status_success
    end subroutine flatten_plan_displacements

    ! Displacements are packed in Fortran column-major order: one column per
    ! output pixel, with a source-rank component along the first dimension.
    ! With no affine base, source rank must equal destination rank, and each
    ! displacement is relative to the one-based identity mapping.
    ! With matrix and offset, mapped coordinates are
    ! 1 + offset + matmul(matrix, output_index - 1) + displacements(:,pixel).
    ! Return the warped image and its partial derivatives with respect to
    ! the displacement components at each output pixel. These are derivatives
    ! with respect to source coordinates, NOT spatial derivatives of the
    ! displacement field or the Jacobian of the affine base transform.
    ! The gradient is flattened (source rank, number of output pixels).
    subroutine evaluate_plan_displacement_gradient_real(plan,displacements,y,gradient,status, &
                                                        matrix,offset,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! cached real source
        real(dp), intent(in) :: displacements(:,:) !! (source rank, output pixels)
        real(dp), intent(out) :: y(..) !! warped real image, rank 1..15
        real(dp), intent(out) :: gradient(:,:) !! derivatives of output w.r.t. displacements
        integer, intent(out) :: status !! zero for success
        real(dp), optional, intent(in) :: matrix(:,:) !! affine source-by-output map
        real(dp), optional, intent(in) :: offset(:) !! affine source offset
        integer, optional, intent(in) :: batch_size !! coordinate batch size
        complex(dp), allocatable :: flat(:),grad(:,:)
        integer :: nr

        status=spline_status_invalid_input
        gradient=0.0_dp
        if (.not.allocated(plan%coefficients)) return
        if (plan%complex_input) return
        if (rank(y)<1 .or. rank(y)>15) return
        nr=size(plan%original_shape)
        if (size(y,kind=int64)>max_padded_elements/int(nr,int64)) then
            status=spline_nd_status_too_expensive
            return
        end if
        if (size(gradient,1)/=nr .or. size(gradient,2,kind=int64)/=size(y,kind=int64)) return
        allocate(flat(size(y)),grad(nr,size(y)))
        call evaluate_plan_transform_core(plan,shape(y),flat,status, &
            matrix=matrix,offset=offset,batch_size=batch_size, &
            displacements=displacements,gradient=grad)
        if (status==spline_status_success) then
            call scatter_real(real(flat,dp),y)
            gradient=real(grad,dp)
        end if
    end subroutine evaluate_plan_displacement_gradient_real

    subroutine evaluate_plan_displacement_gradient_complex(plan,displacements,y,gradient,status, &
                                                           matrix,offset,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! cached complex source
        real(dp), intent(in) :: displacements(:,:) !! (source rank, output pixels)
        complex(dp), intent(out) :: y(..) !! warped complex image, rank 1..15
        complex(dp), intent(out) :: gradient(:,:) !! derivatives of output w.r.t. displacements
        integer, intent(out) :: status !! zero for success
        real(dp), optional, intent(in) :: matrix(:,:) !! affine source-by-output map
        real(dp), optional, intent(in) :: offset(:) !! affine source offset
        integer, optional, intent(in) :: batch_size !! coordinate batch size
        complex(dp), allocatable :: flat(:)
        integer :: nr

        status=spline_status_invalid_input
        gradient=(0.0_dp,0.0_dp)
        if (.not.allocated(plan%coefficients)) return
        if (.not.plan%complex_input) return
        if (rank(y)<1 .or. rank(y)>15) return
        nr=size(plan%original_shape)
        if (size(y,kind=int64)>max_padded_elements/int(nr,int64)) then
            status=spline_nd_status_too_expensive
            return
        end if
        if (size(gradient,1)/=nr .or. size(gradient,2,kind=int64)/=size(y,kind=int64)) return
        allocate(flat(size(y)))
        call evaluate_plan_transform_core(plan,shape(y),flat,status, &
            matrix=matrix,offset=offset,batch_size=batch_size, &
            displacements=displacements,gradient=gradient)
        if (status==spline_status_success) call scatter_complex(flat,y)
    end subroutine evaluate_plan_displacement_gradient_complex

    subroutine evaluate_plan_displacement_real(plan,displacements,y,status,matrix,offset,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! cached real image coefficients
        real(dp), intent(in) :: displacements(:,:) !! (source rank, number of output pixels)
        real(dp), intent(out) :: y(..) !! real destination, rank 1..15
        integer, intent(out) :: status !! zero on success
        real(dp), optional, intent(in) :: matrix(:,:) !! zero-based affine source-by-output map
        real(dp), optional, intent(in) :: offset(:) !! zero-based source offset
        integer, optional, intent(in) :: batch_size !! output pixels per interpolation block
        complex(dp), allocatable :: flat(:)

        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (plan%complex_input) return
        if (rank(y)<1 .or. rank(y)>15) return
        if (size(y,kind=int64)>max_padded_elements) then
            status=spline_nd_status_too_expensive
            return
        end if
        allocate(flat(size(y)))
        call evaluate_plan_transform_core(plan,shape(y),flat,status, &
            matrix=matrix,offset=offset,batch_size=batch_size,displacements=displacements)
        if (status==spline_status_success) call scatter_real(real(flat,dp),y)
    end subroutine evaluate_plan_displacement_real

    subroutine evaluate_plan_displacement_complex(plan,displacements,y,status,matrix,offset,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! cached complex image coefficients
        real(dp), intent(in) :: displacements(:,:) !! (source rank, number of output pixels)
        complex(dp), intent(out) :: y(..) !! complex destination, rank 1..15
        integer, intent(out) :: status !! zero on success
        real(dp), optional, intent(in) :: matrix(:,:) !! zero-based affine source-by-output map
        real(dp), optional, intent(in) :: offset(:) !! zero-based source offset
        integer, optional, intent(in) :: batch_size !! output pixels per interpolation block
        complex(dp), allocatable :: flat(:)

        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (.not.plan%complex_input) return
        if (rank(y)<1 .or. rank(y)>15) return
        if (size(y,kind=int64)>max_padded_elements) then
            status=spline_nd_status_too_expensive
            return
        end if
        allocate(flat(size(y)))
        call evaluate_plan_transform_core(plan,shape(y),flat,status, &
            matrix=matrix,offset=offset,batch_size=batch_size,displacements=displacements)
        if (status==spline_status_success) call scatter_complex(flat,y)
    end subroutine evaluate_plan_displacement_complex

    ! Convenience interface for a dense field with a leading source-axis
    ! component. The output may have at most 14 axes in this representation.
    subroutine evaluate_plan_displacement_grid_real(plan,field,y,status,matrix,offset,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! cached real image coefficients
        real(dp), intent(in) :: field(..) !! shape (source rank, output extents...)
        real(dp), intent(out) :: y(..) !! real destination, rank 1..14
        integer, intent(out) :: status !! zero on success
        real(dp), optional, intent(in) :: matrix(:,:) !! optional zero-based affine matrix
        real(dp), optional, intent(in) :: offset(:) !! optional zero-based affine offset
        integer, optional, intent(in) :: batch_size !! interpolation block size
        real(dp), allocatable :: flat_field(:,:)

        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (plan%complex_input) return
        call flatten_plan_displacements(field,size(plan%original_shape),shape(y),flat_field,status)
        if (status/=spline_status_success) return
        call evaluate_plan_displacement_real(plan,flat_field,y,status,matrix,offset,batch_size)
    end subroutine evaluate_plan_displacement_grid_real

    subroutine evaluate_plan_displacement_grid_complex(plan,field,y,status,matrix,offset,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! cached complex image coefficients
        real(dp), intent(in) :: field(..) !! shape (source rank, output extents...)
        complex(dp), intent(out) :: y(..) !! complex destination, rank 1..14
        integer, intent(out) :: status !! zero on success
        real(dp), optional, intent(in) :: matrix(:,:) !! optional zero-based affine matrix
        real(dp), optional, intent(in) :: offset(:) !! optional zero-based affine offset
        integer, optional, intent(in) :: batch_size !! interpolation block size
        real(dp), allocatable :: flat_field(:,:)

        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (.not.plan%complex_input) return
        call flatten_plan_displacements(field,size(plan%original_shape),shape(y),flat_field,status)
        if (status/=spline_status_success) return
        call evaluate_plan_displacement_complex(plan,flat_field,y,status,matrix,offset,batch_size)
    end subroutine evaluate_plan_displacement_grid_complex

    ! Affine parameter Hessians. The last index is a flattened output pixel.
    ! The parameter maps are x_i = 1 + offset_i + sum_j matrix_ij * (out_j-1).
    ! hmatrix(a,j,b,k,p) = d2 value_p / dmatrix(a,j) dmatrix(b,k).
    ! hcross(a,j,b,p) = d2 value_p / dmatrix(a,j) doffset(b).
    ! hoffset(a,b,p) = d2 value_p / doffset(a) doffset(b).
    subroutine evaluate_plan_affine_hessian_real(plan,matrix,offset,y,dmatrix,doffset, &
                                                hmatrix,hcross,hoffset,status,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! cached real spline coefficients
        real(dp), intent(in) :: matrix(:,:) !! source-by-destination affine map
        real(dp), intent(in) :: offset(:) !! zero-based source offset
        real(dp), intent(out) :: y(..) !! transformed image
        real(dp), intent(out) :: dmatrix(:,:,:) !! (input rank, output rank, pixels)
        real(dp), intent(out) :: doffset(:,:) !! (input rank, pixels)
        real(dp), intent(out) :: hmatrix(:,:,:,:,:) !! (input rank, output rank, input rank, output rank, pixels)
        real(dp), intent(out) :: hcross(:,:,:,:) !! (input rank, output rank, input rank, pixels)
        real(dp), intent(out) :: hoffset(:,:,:) !! (input rank, input rank, pixels)
        integer, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: batch_size !! maximum block size
        complex(dp), allocatable :: v(:),dm(:,:,:),db(:,:),hm(:,:,:,:,:),hc(:,:,:,:),hb(:,:,:)
        integer :: nr,nd,np

        status=spline_status_invalid_input
        dmatrix=0.0_dp
        doffset=0.0_dp
        hmatrix=0.0_dp
        hcross=0.0_dp
        hoffset=0.0_dp
        if (.not.allocated(plan%coefficients)) return
        if (plan%complex_input) return
        nr=size(plan%original_shape)
        nd=rank(y)
        if (nd<1 .or. nd>15) return
        if (any(shape(y)<1)) return
        np=size(y)
        if (int(np,int64)>max_padded_elements/int(nr*nr*(nd+1)*(nd+1),int64)) then
            status=spline_nd_status_too_expensive
            return
        end if
        if (any(shape(dmatrix)/=[nr,nd,np])) return
        if (any(shape(doffset)/=[nr,np])) return
        if (any(shape(hmatrix)/=[nr,nd,nr,nd,np])) return
        if (any(shape(hcross)/=[nr,nd,nr,np])) return
        if (any(shape(hoffset)/=[nr,nr,np])) return
        allocate(v(np),dm(nr,nd,np),db(nr,np),hm(nr,nd,nr,nd,np),hc(nr,nd,nr,np),hb(nr,nr,np))
        call evaluate_plan_transform_core(plan,shape(y),v,status, &
            matrix=matrix,offset=offset,batch_size=batch_size,dmatrix=dm,doffset=db, &
            hmm=hm,hmo=hc,hoo=hb)
        if (status/=spline_status_success) return
        call scatter_real(real(v,dp),y)
        dmatrix=real(dm,dp)
        doffset=real(db,dp)
        hmatrix=real(hm,dp)
        hcross=real(hc,dp)
        hoffset=real(hb,dp)
    end subroutine evaluate_plan_affine_hessian_real

    subroutine evaluate_plan_affine_hessian_complex(plan,matrix,offset,y,dmatrix,doffset, &
                                                   hmatrix,hcross,hoffset,status,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! cached complex spline coefficients
        real(dp), intent(in) :: matrix(:,:) !! source-by-destination affine map
        real(dp), intent(in) :: offset(:) !! zero-based source offset
        complex(dp), intent(out) :: y(..) !! transformed image
        complex(dp), intent(out) :: dmatrix(:,:,:) !! (input rank, output rank, pixels)
        complex(dp), intent(out) :: doffset(:,:) !! (input rank, pixels)
        complex(dp), intent(out) :: hmatrix(:,:,:,:,:) !! second matrix derivatives
        complex(dp), intent(out) :: hcross(:,:,:,:) !! matrix-offset derivatives
        complex(dp), intent(out) :: hoffset(:,:,:) !! offset-offset derivatives
        integer, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: batch_size !! maximum block size
        complex(dp), allocatable :: v(:)
        integer :: nr,nd,np

        status=spline_status_invalid_input
        dmatrix=(0.0_dp,0.0_dp)
        doffset=(0.0_dp,0.0_dp)
        hmatrix=(0.0_dp,0.0_dp)
        hcross=(0.0_dp,0.0_dp)
        hoffset=(0.0_dp,0.0_dp)
        if (.not.allocated(plan%coefficients)) return
        if (.not.plan%complex_input) return
        nr=size(plan%original_shape)
        nd=rank(y)
        if (nd<1 .or. nd>15) return
        if (any(shape(y)<1)) return
        np=size(y)
        if (int(np,int64)>max_padded_elements/int(nr*nr*(nd+1)*(nd+1),int64)) then
            status=spline_nd_status_too_expensive
            return
        end if
        if (any(shape(dmatrix)/=[nr,nd,np])) return
        if (any(shape(doffset)/=[nr,np])) return
        if (any(shape(hmatrix)/=[nr,nd,nr,nd,np])) return
        if (any(shape(hcross)/=[nr,nd,nr,np])) return
        if (any(shape(hoffset)/=[nr,nr,np])) return
        allocate(v(np))
        call evaluate_plan_transform_core(plan,shape(y),v,status, &
            matrix=matrix,offset=offset,batch_size=batch_size,dmatrix=dmatrix,doffset=doffset, &
            hmm=hmatrix,hmo=hcross,hoo=hoffset)
        if (status==spline_status_success) call scatter_complex(v,y)
    end subroutine evaluate_plan_affine_hessian_complex

    ! A displacement is additive in source-coordinate space. Its Hessian
    ! is therefore exactly the spline coordinate Hessian at each mapped site.
    subroutine evaluate_plan_displacement_hessian_real(plan,displacements,y,gradient,hessian,status, &
                                                       matrix,offset,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! cached real spline coefficients
        real(dp), intent(in) :: displacements(:,:) !! (source rank, output pixels)
        real(dp), intent(out) :: y(..) !! warped image, ranks 1 through 15
        real(dp), intent(out) :: gradient(:,:) !! (source rank, output pixels)
        real(dp), intent(out) :: hessian(:,:,:) !! (source rank, source rank, output pixels)
        integer, intent(out) :: status !! zero on success
        real(dp), optional, intent(in) :: matrix(:,:) !! zero-based affine base
        real(dp), optional, intent(in) :: offset(:) !! zero-based affine offset
        integer, optional, intent(in) :: batch_size !! maximum batch size
        complex(dp), allocatable :: v(:),g(:,:),h(:,:,:)
        integer :: nr,np

        status=spline_status_invalid_input
        gradient=0.0_dp
        hessian=0.0_dp
        if (.not.allocated(plan%coefficients)) return
        if (plan%complex_input) return
        nr=size(plan%original_shape)
        np=size(y)
        if (rank(y)<1 .or. rank(y)>15) return
        if (int(np,int64)>max_padded_elements/int(nr*nr,int64)) then
            status=spline_nd_status_too_expensive
            return
        end if
        if (any(shape(gradient)/=[nr,np]) .or. any(shape(hessian)/=[nr,nr,np])) return
        allocate(v(np),g(nr,np),h(nr,nr,np))
        call evaluate_plan_transform_core(plan,shape(y),v,status, &
            matrix=matrix,offset=offset,batch_size=batch_size,displacements=displacements, &
            gradient=g,hessian=h)
        if (status/=spline_status_success) return
        call scatter_real(real(v,dp),y)
        gradient=real(g,dp)
        hessian=real(h,dp)
    end subroutine evaluate_plan_displacement_hessian_real

    subroutine evaluate_plan_displacement_hessian_complex(plan,displacements,y,gradient,hessian,status, &
                                                          matrix,offset,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! cached complex spline coefficients
        real(dp), intent(in) :: displacements(:,:) !! (source rank, output pixels)
        complex(dp), intent(out) :: y(..) !! warped image
        complex(dp), intent(out) :: gradient(:,:) !! (source rank, output pixels)
        complex(dp), intent(out) :: hessian(:,:,:) !! (source rank, source rank, pixels)
        integer, intent(out) :: status !! zero on success
        real(dp), optional, intent(in) :: matrix(:,:) !! zero-based affine base
        real(dp), optional, intent(in) :: offset(:) !! zero-based affine offset
        integer, optional, intent(in) :: batch_size !! maximum batch size
        complex(dp), allocatable :: v(:)
        integer :: nr,np

        status=spline_status_invalid_input
        gradient=(0.0_dp,0.0_dp)
        hessian=(0.0_dp,0.0_dp)
        if (.not.allocated(plan%coefficients)) return
        if (.not.plan%complex_input) return
        nr=size(plan%original_shape)
        np=size(y)
        if (rank(y)<1 .or. rank(y)>15) return
        if (int(np,int64)>max_padded_elements/int(nr*nr,int64)) then
            status=spline_nd_status_too_expensive
            return
        end if
        if (any(shape(gradient)/=[nr,np]) .or. any(shape(hessian)/=[nr,nr,np])) return
        allocate(v(np))
        call evaluate_plan_transform_core(plan,shape(y),v,status, &
            matrix=matrix,offset=offset,batch_size=batch_size,displacements=displacements, &
            gradient=gradient,hessian=hessian)
        if (status==spline_status_success) call scatter_complex(v,y)
    end subroutine evaluate_plan_displacement_hessian_complex

    ! Grid-shaped displacement Hessians require two leading source axes and
    ! hence a destination of rank at most 13 under Fortran's rank-15 limit.
    subroutine evaluate_plan_displacement_grid_hessian_real(plan,field,y,gradient,hessian,status, &
                                                            matrix,offset,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! cached real image coefficients
        real(dp), intent(in) :: field(..) !! shape [source rank, shape(y)]
        real(dp), intent(out) :: y(..) !! warped grid, rank 1 through 13
        real(dp), intent(out) :: gradient(..) !! shape [source rank, shape(y)]
        real(dp), intent(out) :: hessian(..) !! shape [source rank, source rank, shape(y)]
        integer, intent(out) :: status !! zero on success
        real(dp), optional, intent(in) :: matrix(:,:) !! affine base map
        real(dp), optional, intent(in) :: offset(:) !! affine base offset
        integer, optional, intent(in) :: batch_size !! interpolation block size
        real(dp), allocatable :: flat_field(:,:),g(:,:),h(:,:,:)
        integer :: nr

        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (plan%complex_input) return
        if (rank(y)<1 .or. rank(y)>13) return
        nr=size(plan%original_shape)
        if (rank(gradient)/=rank(y)+1 .or. rank(hessian)/=rank(y)+2) return
        if (any(shape(gradient)/=[nr,shape(y)])) return
        if (any(shape(hessian)/=[nr,nr,shape(y)])) return
        if (size(y,kind=int64)>max_padded_elements/int(nr*nr,int64)) then
            status=spline_nd_status_too_expensive
            return
        end if
        call flatten_plan_displacements(field,nr,shape(y),flat_field,status)
        if (status/=spline_status_success) return
        allocate(g(nr,size(y)),h(nr,nr,size(y)))
        call evaluate_plan_displacement_hessian_real(plan,flat_field,y,g,h,status,matrix,offset,batch_size)
        if (status/=spline_status_success) return
        call scatter_real(reshape(g,[size(g)]),gradient)
        call scatter_real(reshape(h,[size(h)]),hessian)
    end subroutine evaluate_plan_displacement_grid_hessian_real

    subroutine evaluate_plan_displacement_grid_hessian_complex(plan,field,y,gradient,hessian,status, &
                                                               matrix,offset,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! cached complex image coefficients
        real(dp), intent(in) :: field(..) !! shape [source rank, shape(y)]
        complex(dp), intent(out) :: y(..) !! warped grid, rank 1 through 13
        complex(dp), intent(out) :: gradient(..) !! shape [source rank, shape(y)]
        complex(dp), intent(out) :: hessian(..) !! shape [source rank, source rank, shape(y)]
        integer, intent(out) :: status !! zero on success
        real(dp), optional, intent(in) :: matrix(:,:) !! affine base map
        real(dp), optional, intent(in) :: offset(:) !! affine base offset
        integer, optional, intent(in) :: batch_size !! interpolation block size
        real(dp), allocatable :: flat_field(:,:)
        complex(dp), allocatable :: g(:,:),h(:,:,:)
        integer :: nr

        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (.not.plan%complex_input) return
        if (rank(y)<1 .or. rank(y)>13) return
        nr=size(plan%original_shape)
        if (rank(gradient)/=rank(y)+1 .or. rank(hessian)/=rank(y)+2) return
        if (any(shape(gradient)/=[nr,shape(y)])) return
        if (any(shape(hessian)/=[nr,nr,shape(y)])) return
        if (size(y,kind=int64)>max_padded_elements/int(nr*nr,int64)) then
            status=spline_nd_status_too_expensive
            return
        end if
        call flatten_plan_displacements(field,nr,shape(y),flat_field,status)
        if (status/=spline_status_success) return
        allocate(g(nr,size(y)),h(nr,nr,size(y)))
        call evaluate_plan_displacement_hessian_complex(plan,flat_field,y,g,h,status,matrix,offset,batch_size)
        if (status/=spline_status_success) return
        call scatter_complex(reshape(g,[size(g)]),gradient)
        call scatter_complex(reshape(h,[size(h)]),hessian)
    end subroutine evaluate_plan_displacement_grid_hessian_complex

    ! Analytic derivatives of an affine warp with respect to the matrix and
    ! offset parameters. Matrix and offset use SciPy's zero-based convention:
    ! x_source=1+offset+matrix*(x_output-1). Chain rule gives
    ! dI/d(matrix(i,j))=(x_output(j)-1)*dI/d(x_source(i)),
    ! dI/d(offset(i))=dI/d(x_source(i)). No finite differences are used.
    ! Parameter derivatives are returned in flattened column-major pixel order.
    subroutine evaluate_plan_affine_gradient_real(plan,matrix,offset,y,dmatrix,doffset,status,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! prepared real-valued plan
        real(dp), intent(in) :: matrix(:,:) !! input-rank by output-rank affine coefficients
        real(dp), intent(in) :: offset(:) !! zero-based source offsets
        real(dp), intent(out) :: y(..) !! transformed image
        real(dp), intent(out) :: dmatrix(:,:,:) !! (source rank, output rank, pixels)
        real(dp), intent(out) :: doffset(:,:) !! (source rank, pixels)
        integer, intent(out) :: status !! success or validation/resource error
        integer, optional, intent(in) :: batch_size !! samples per interpolation batch
        complex(dp), allocatable :: v(:),dm(:,:,:),db(:,:)
        integer :: nr,nd,np

        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (plan%complex_input) return
        nr=size(plan%original_shape)
        nd=rank(y)
        if (nd<1 .or. nd>15) return
        if (any(shape(y)<1)) return
        if (size(y,kind=int64)>max_padded_elements) then
            status=spline_nd_status_too_expensive
            return
        end if
        np=size(y)
        if (any(shape(dmatrix)/=[nr,nd,np])) return
        if (any(shape(doffset)/=[nr,np])) return
        if (int(np,int64)>max_padded_elements/int(nr*nd,int64)) then
            status=spline_nd_status_too_expensive
            return
        end if
        allocate(v(np),dm(nr,nd,np),db(nr,np))
        call evaluate_plan_transform_core(plan,shape(y),v,status, &
            matrix=matrix,offset=offset,batch_size=batch_size,dmatrix=dm,doffset=db)
        if (status/=spline_status_success) return
        call scatter_real(real(v,dp),y)
        dmatrix=real(dm,dp)
        doffset=real(db,dp)
    end subroutine evaluate_plan_affine_gradient_real

    subroutine evaluate_plan_affine_gradient_complex(plan,matrix,offset,y,dmatrix,doffset,status,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! prepared complex-valued plan
        real(dp), intent(in) :: matrix(:,:) !! input-rank by output-rank affine coefficients
        real(dp), intent(in) :: offset(:) !! zero-based source offsets
        complex(dp), intent(out) :: y(..) !! transformed image
        complex(dp), intent(out) :: dmatrix(:,:,:) !! (source rank, output rank, pixels)
        complex(dp), intent(out) :: doffset(:,:) !! (source rank, pixels)
        integer, intent(out) :: status !! success or validation/resource error
        integer, optional, intent(in) :: batch_size !! samples per interpolation batch
        complex(dp), allocatable :: v(:)
        integer :: nr,nd,np

        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (.not.plan%complex_input) return
        nr=size(plan%original_shape)
        nd=rank(y)
        if (nd<1 .or. nd>15) return
        if (any(shape(y)<1)) return
        if (size(y,kind=int64)>max_padded_elements) then
            status=spline_nd_status_too_expensive
            return
        end if
        np=size(y)
        if (any(shape(dmatrix)/=[nr,nd,np])) return
        if (any(shape(doffset)/=[nr,np])) return
        if (int(np,int64)>max_padded_elements/int(nr*nd,int64)) then
            status=spline_nd_status_too_expensive
            return
        end if
        allocate(v(np))
        call evaluate_plan_transform_core(plan,shape(y),v,status, &
            matrix=matrix,offset=offset,batch_size=batch_size,dmatrix=dmatrix,doffset=doffset)
        if (status/=spline_status_success) return
        call scatter_complex(v,y)
    end subroutine evaluate_plan_affine_gradient_complex

    subroutine evaluate_plan_affine_real(plan,matrix,offset,y,status,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! real-valued prepared plan
        real(dp), intent(in) :: matrix(:,:) !! shape (input rank, output rank)
        real(dp), intent(in) :: offset(:) !! source offset, in zero-based units
        real(dp), intent(out) :: y(..) !! caller-allocated output, rank 1 through 15
        integer, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: batch_size !! maximum coordinates per batch
        complex(dp), allocatable :: flat(:)
        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (plan%complex_input) return
        if (rank(y)<1 .or. rank(y)>15) return
        if (size(y,kind=int64)>max_padded_elements) then
            status=spline_nd_status_too_expensive
            return
        end if
        allocate(flat(size(y)))
        call evaluate_plan_transform_core(plan,shape(y),flat,status, &
            matrix=matrix,offset=offset,batch_size=batch_size)
        if (status==spline_status_success) call scatter_real(real(flat,dp),y)
    end subroutine evaluate_plan_affine_real

    subroutine evaluate_plan_affine_complex(plan,matrix,offset,y,status,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! complex-valued prepared plan
        real(dp), intent(in) :: matrix(:,:) !! shape (input rank, output rank)
        real(dp), intent(in) :: offset(:) !! source offset, in zero-based units
        complex(dp), intent(out) :: y(..) !! caller-allocated output, rank 1 through 15
        integer, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: batch_size !! maximum coordinates per batch
        complex(dp), allocatable :: flat(:)
        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (.not.plan%complex_input) return
        if (rank(y)<1 .or. rank(y)>15) return
        if (size(y,kind=int64)>max_padded_elements) then
            status=spline_nd_status_too_expensive
            return
        end if
        allocate(flat(size(y)))
        call evaluate_plan_transform_core(plan,shape(y),flat,status, &
            matrix=matrix,offset=offset,batch_size=batch_size)
        if (status==spline_status_success) call scatter_complex(flat,y)
    end subroutine evaluate_plan_affine_complex

    ! SciFort extension: the callback transforms an entire block of points.
    ! This saves one procedure call per pixel and permits array-based maps.
    ! A final block may contain fewer columns than batch_size.
    subroutine evaluate_plan_geometric_batch_real(plan,mapping,y,status,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! prepared real source
        procedure(ndimage_coordinate_batch_map) :: mapping !! batched output-to-input mapping
        real(dp), intent(out) :: y(..) !! caller-allocated real output, rank 1 through 15
        integer, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: batch_size !! maximum columns per callback, 1..8192
        complex(dp), allocatable :: flat(:)

        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (plan%complex_input) return
        if (rank(y)<1 .or. rank(y)>15) return
        if (size(y,kind=int64)>max_padded_elements) then
            status=spline_nd_status_too_expensive
            return
        end if
        allocate(flat(size(y)))
        call evaluate_plan_transform_core(plan,shape(y),flat,status, &
            batch_mapping=mapping,batch_size=batch_size)
        if (status==spline_status_success) call scatter_real(real(flat,dp),y)
    end subroutine evaluate_plan_geometric_batch_real

    subroutine evaluate_plan_geometric_batch_complex(plan,mapping,y,status,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! prepared complex source
        procedure(ndimage_coordinate_batch_map) :: mapping !! batched output-to-input mapping
        complex(dp), intent(out) :: y(..) !! caller-allocated complex output, rank 1 through 15
        integer, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: batch_size !! maximum columns per callback, 1..8192
        complex(dp), allocatable :: flat(:)

        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (.not.plan%complex_input) return
        if (rank(y)<1 .or. rank(y)>15) return
        if (size(y,kind=int64)>max_padded_elements) then
            status=spline_nd_status_too_expensive
            return
        end if
        allocate(flat(size(y)))
        call evaluate_plan_transform_core(plan,shape(y),flat,status, &
            batch_mapping=mapping,batch_size=batch_size)
        if (status==spline_status_success) call scatter_complex(flat,y)
    end subroutine evaluate_plan_geometric_batch_complex

    subroutine evaluate_plan_geometric_real(plan,mapping,y,status,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! real-valued prepared plan
        procedure(ndimage_coordinate_map) :: mapping !! one-based output-to-source mapping
        real(dp), intent(out) :: y(..) !! caller-allocated output, rank 1 through 15
        integer, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: batch_size !! maximum coordinates per batch
        complex(dp), allocatable :: flat(:)
        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (plan%complex_input) return
        if (rank(y)<1 .or. rank(y)>15) return
        if (size(y,kind=int64)>max_padded_elements) then
            status=spline_nd_status_too_expensive
            return
        end if
        allocate(flat(size(y)))
        call evaluate_plan_transform_core(plan,shape(y),flat,status, &
            mapping=mapping,batch_size=batch_size)
        if (status==spline_status_success) call scatter_real(real(flat,dp),y)
    end subroutine evaluate_plan_geometric_real

    subroutine evaluate_plan_geometric_complex(plan,mapping,y,status,batch_size)
        type(spline_interpolation_plan), intent(in) :: plan !! complex-valued prepared plan
        procedure(ndimage_coordinate_map) :: mapping !! one-based output-to-source mapping
        complex(dp), intent(out) :: y(..) !! caller-allocated output, rank 1 through 15
        integer, intent(out) :: status !! zero on success
        integer, optional, intent(in) :: batch_size !! maximum coordinates per batch
        complex(dp), allocatable :: flat(:)
        status=spline_status_invalid_input
        if (.not.allocated(plan%coefficients)) return
        if (.not.plan%complex_input) return
        if (rank(y)<1 .or. rank(y)>15) return
        if (size(y,kind=int64)>max_padded_elements) then
            status=spline_nd_status_too_expensive
            return
        end if
        allocate(flat(size(y)))
        call evaluate_plan_transform_core(plan,shape(y),flat,status, &
            mapping=mapping,batch_size=batch_size)
        if (status==spline_status_success) call scatter_complex(flat,y)
    end subroutine evaluate_plan_geometric_complex

    subroutine coordinates_real(coeff,coords,values,order,status,mode,prefilter,cval)
        real(dp), intent(in) :: coeff(..) !! coefficients, rank 1..15
        real(dp), intent(in) :: coords(:,:) !! coordinates: (rank, npoints), ONE-BASED
        real(dp), intent(out) :: values(:) !! output values, length npoints
        integer, optional, intent(in) :: order !! order zero through five; default three
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! ndimage boundary extension
        logical, optional, intent(in) :: prefilter !! apply spline filtering first; default false
        real(dp), optional, intent(in) :: cval !! constant value beyond boundaries
        real(dp), allocatable :: input(:)
        complex(dp), allocatable :: work(:), result(:)
        complex(dp) :: fill
        character(len=:), allocatable :: boundary
        logical :: filter_first
        integer :: degree, code
        integer, allocatable :: dims(:)

        values=0.0_dp
        code=spline_status_invalid_input
        if (present(status)) status=code
        if (rank(coeff)<1 .or. rank(coeff)>15) return
        if (any(shape(coeff)<1)) return
        if (size(coeff,kind=int64)>int(huge(1),int64)) return
        if (size(coords,1)/=rank(coeff)) return
        if (size(coords,2)/=size(values)) return
        degree=3
        if (present(order)) degree=order
        if (degree<0 .or. degree>5) return
        boundary='mirror'
        if (present(mode)) boundary=trim(mode)
        if (.not. valid_boundary(boundary)) return
        fill=0.0_dp
        if (present(cval)) fill=cmplx(cval,0.0_dp,dp)
        filter_first=.false.
        if (present(prefilter)) filter_first=prefilter
        call pack_real(coeff,input)
        work=cmplx(input,0.0_dp,dp)
        if (any(.not. ieee_is_finite(real(work,dp))) .or. &
            any(.not. ieee_is_finite(aimag(work)))) return
        dims=shape(coeff)
        allocate(result(size(values)))
        if (filter_first .and. degree>1) then
            block
                complex(dp), allocatable :: coeff(:)
                integer, allocatable :: expanded(:), halo(:)
                call prepare_interpolation_coeff(work,dims,degree,boundary,fill, &
                                                 coeff,expanded,halo,code)
                if (code/=spline_status_success) then
                    if (present(status)) status=code
                    return
                end if
                call evaluate_interpolation(coeff,expanded,halo,dims,coords, &
                                            degree,boundary,fill,result,code)
            end block
        else
            call interpolate(work,dims,coords,result,degree,boundary,fill,code)
        end if
        if (code/=spline_status_success) return
        values=real(result,dp)
        if (present(status)) status=code
    end subroutine coordinates_real

    subroutine coordinates_complex(coeff,coords,values,order,status,mode,prefilter,cval)
        complex(dp), intent(in) :: coeff(..) !! coefficients, rank 1..15
        real(dp), intent(in) :: coords(:,:) !! coordinates: (rank, npoints), ONE-BASED
        complex(dp), intent(out) :: values(:) !! output values, length npoints
        integer, optional, intent(in) :: order !! order zero through five; default three
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! ndimage boundary extension
        logical, optional, intent(in) :: prefilter !! apply spline filtering first; default false
        complex(dp), optional, intent(in) :: cval !! constant value beyond boundaries
        complex(dp), allocatable :: input(:)
        complex(dp), allocatable :: work(:), result(:)
        complex(dp) :: fill
        character(len=:), allocatable :: boundary
        logical :: filter_first
        integer :: degree, code
        integer, allocatable :: dims(:)

        values=(0.0_dp,0.0_dp)
        code=spline_status_invalid_input
        if (present(status)) status=code
        if (rank(coeff)<1 .or. rank(coeff)>15) return
        if (any(shape(coeff)<1)) return
        if (size(coeff,kind=int64)>int(huge(1),int64)) return
        if (size(coords,1)/=rank(coeff)) return
        if (size(coords,2)/=size(values)) return
        degree=3
        if (present(order)) degree=order
        if (degree<0 .or. degree>5) return
        boundary='mirror'
        if (present(mode)) boundary=trim(mode)
        if (.not. valid_boundary(boundary)) return
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cval
        filter_first=.false.
        if (present(prefilter)) filter_first=prefilter
        call pack_complex(coeff,input)
        work=input
        if (any(.not. ieee_is_finite(real(work,dp))) .or. &
            any(.not. ieee_is_finite(aimag(work)))) return
        dims=shape(coeff)
        allocate(result(size(values)))
        if (filter_first .and. degree>1) then
            block
                complex(dp), allocatable :: coeff(:)
                integer, allocatable :: expanded(:), halo(:)
                call prepare_interpolation_coeff(work,dims,degree,boundary,fill, &
                                                 coeff,expanded,halo,code)
                if (code/=spline_status_success) then
                    if (present(status)) status=code
                    return
                end if
                call evaluate_interpolation(coeff,expanded,halo,dims,coords, &
                                            degree,boundary,fill,result,code)
            end block
        else
            call interpolate(work,dims,coords,result,degree,boundary,fill,code)
        end if
        if (code/=spline_status_success) return
        values=result
        if (present(status)) status=code
    end subroutine coordinates_complex

    ! Shape-preserving map_coordinates variant. The first coordinate axis
    ! selects input dimensions and every remaining axis shapes the output.
    ! Unlike map_coordinates_nd, this new interface defaults to prefiltering.
    subroutine coordinate_grid_real(x,coords,y,order,status,mode,prefilter,cval)
        real(dp), intent(in) :: x(..) !! real input, rank 1 through 15
        real(dp), intent(in) :: coords(..) !! shape (rank(x), shape(y))
        real(dp), intent(out) :: y(..) !! output, rank 1 through 14
        integer, optional, intent(in) :: order !! interpolation degree, default 3
        integer, optional, intent(out) :: status !! spline status
        character(len=*), optional, intent(in) :: mode !! boundary mode
        logical, optional, intent(in) :: prefilter !! default true
        real(dp), optional, intent(in) :: cval !! exterior value
        real(dp), allocatable :: packed(:), locations(:,:), values(:)
        integer, allocatable :: dims(:)
        integer :: degree, code
        logical :: do_prefilter

        code=spline_status_invalid_input
        if (present(status)) status=code
        if (rank(x)<1 .or. rank(x)>15) return
        if (rank(coords)<2 .or. rank(coords)>15) return
        if (rank(y)/=rank(coords)-1) return
        dims=shape(coords)
        if (dims(1)/=rank(x)) return
        if (any(dims(2:)/=shape(y))) return
        if (any(shape(x)<1) .or. any(shape(y)<1)) return
        if (size(y,kind=int64)>int(huge(1),int64)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        degree=3
        if (present(order)) degree=order
        do_prefilter=.true.
        if (present(prefilter)) do_prefilter=prefilter
        call pack_real(coords,packed)
        locations=reshape(packed,[rank(x),size(y)])
        allocate(values(size(y)))
        call map_coordinates_nd(x,locations,values,order=degree,status=code, &
            mode=mode,prefilter=do_prefilter,cval=cval)
        if (code==spline_status_success) call scatter_real(values,y)
        if (present(status)) status=code
    end subroutine coordinate_grid_real

    subroutine coordinate_grid_complex(x,coords,y,order,status,mode,prefilter,cval)
        complex(dp), intent(in) :: x(..) !! complex input, rank 1 through 15
        real(dp), intent(in) :: coords(..) !! shape (rank(x), shape(y))
        complex(dp), intent(out) :: y(..) !! output, rank 1 through 14
        integer, optional, intent(in) :: order !! interpolation degree, default 3
        integer, optional, intent(out) :: status !! spline status
        character(len=*), optional, intent(in) :: mode !! boundary mode
        logical, optional, intent(in) :: prefilter !! default true
        complex(dp), optional, intent(in) :: cval !! exterior value
        real(dp), allocatable :: packed(:), locations(:,:)
        complex(dp), allocatable :: values(:)
        integer, allocatable :: dims(:)
        integer :: degree, code
        logical :: do_prefilter

        code=spline_status_invalid_input
        if (present(status)) status=code
        if (rank(x)<1 .or. rank(x)>15) return
        if (rank(coords)<2 .or. rank(coords)>15) return
        if (rank(y)/=rank(coords)-1) return
        dims=shape(coords)
        if (dims(1)/=rank(x)) return
        if (any(dims(2:)/=shape(y))) return
        if (any(shape(x)<1) .or. any(shape(y)<1)) return
        if (size(y,kind=int64)>int(huge(1),int64)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        degree=3
        if (present(order)) degree=order
        do_prefilter=.true.
        if (present(prefilter)) do_prefilter=prefilter
        call pack_real(coords,packed)
        locations=reshape(packed,[rank(x),size(y)])
        allocate(values(size(y)))
        call map_coordinates_nd(x,locations,values,order=degree,status=code, &
            mode=mode,prefilter=do_prefilter,cval=cval)
        if (code==spline_status_success) call scatter_complex(values,y)
        if (present(status)) status=code
    end subroutine coordinate_grid_complex

    subroutine affine_samples(work,in_dims,out_dims,mat,offset,out_flat,degree,boundary,filter_first,fill,status)
        complex(dp), intent(in) :: work(:) !! input in column-major order
        integer, intent(in) :: in_dims(:) !! dimensions of input
        integer, intent(in) :: out_dims(:) !! dimensions of output
        real(dp), intent(in) :: mat(:,:) !! affine matrix
        real(dp), intent(in) :: offset(:) !! offset in 0-based input units
        complex(dp), intent(out) :: out_flat(:) !! column-major destination
        integer, intent(in) :: degree !! 0 through 5
        character(len=*), intent(in) :: boundary !! boundary mode
        logical, intent(in) :: filter_first !! optionally prefilter
        complex(dp), intent(in) :: fill !! constant extension value
        integer, intent(out) :: status !! zero if successful
        complex(dp), allocatable :: coeff(:), blockout(:)
        real(dp), allocatable :: coordinates(:,:)
        integer, allocatable :: coeff_dims(:), halo(:)
        integer, parameter :: batch_size=128
        integer(int64) :: p,base,remain,stride
        integer :: n,d,i,code

        status=spline_status_invalid_input
        n=size(in_dims)
        if (n<1 .or. n>15 .or. size(out_dims)/=n) return
        if (any(in_dims<1) .or. any(out_dims<1)) return
        if (size(mat,1)/=n .or. size(mat,2)/=n) return
        if (size(offset)/=n) return
        if (degree<0 .or. degree>5) return
        if (.not.valid_boundary(boundary)) return
        if (.not.all(ieee_is_finite(mat)) .or. .not.all(ieee_is_finite(offset))) return
        if (.not.ieee_is_finite(real(fill,dp)) .or. .not.ieee_is_finite(aimag(fill))) return
        if (filter_first .and. degree>1) then
            call prepare_interpolation_coeff(work,in_dims,degree,boundary,fill, &
                                             coeff,coeff_dims,halo,code)
            if (code/=spline_status_success) then
                status=code
                return
            end if
        else
            coeff=work
            coeff_dims=in_dims
            halo=spread(0,1,size(in_dims))
        end if
        allocate(coordinates(n,batch_size),blockout(batch_size))
        base=1_int64
        do while(base<=size(out_flat,kind=int64))
            remain=min(int(batch_size,int64),size(out_flat,kind=int64)-base+1_int64)
            do i=1,int(remain)
                p=base+int(i-2,int64)
                do d=1,n
                    stride=int(out_dims(d),int64)
                    coordinates(d,i)=real(modulo(p,stride),dp)
                    p=p/stride
                end do
            end do
            do i=1,int(remain)
                coordinates(:,i)=1.0_dp+offset+matmul(mat,coordinates(:,i))
            end do
            if (filter_first .and. degree>1) then
                call evaluate_interpolation(coeff,coeff_dims,halo,in_dims, &
                    coordinates(:,:int(remain)),degree,boundary,fill, &
                    blockout(:int(remain)),code)
            else
                call interpolate(coeff,in_dims,coordinates(:,:int(remain)), &
                    blockout(:int(remain)),degree,boundary,fill,code)
            end if
            if (code/=spline_status_success) then
                status=code
                return
            end if
            out_flat(base:base+remain-1_int64)=blockout(:int(remain))
            base=base+remain
        end do
        status=spline_status_success
    end subroutine affine_samples

    subroutine affine_real(x,matrix,offset,y,order,status,mode,prefilter,cval)
        real(dp), intent(in) :: x(..) !! rank 1..15 input
        real(dp), intent(in) :: matrix(:,:) !! square affine transform
        real(dp), intent(in) :: offset(:) !! 0-based offset into input
        real(dp), intent(out) :: y(..) !! caller-allocated destination
        integer, optional, intent(in) :: order !! interpolation order, default 3
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary extension
        logical, optional, intent(in) :: prefilter !! default false
        real(dp), optional, intent(in) :: cval !! outside value
        real(dp), allocatable :: flat(:), yflat(:)
        complex(dp), allocatable :: work(:), res(:)
        complex(dp) :: fill
        character(len=:), allocatable :: boundary
        logical :: do_filter
        integer :: degree,code

        code=spline_status_invalid_input
        if (present(status)) status=code
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x)) return
        if (any(shape(x)<1) .or. any(shape(y)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64) .or. &
            size(y,kind=int64)>int(huge(1),int64)) return
        degree=3
        if (present(order)) degree=order
        boundary='mirror'
        if (present(mode)) boundary=trim(mode)
        do_filter=.false.
        if (present(prefilter)) do_filter=prefilter
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cmplx(cval,0.0_dp,dp)
        call pack_real(x,flat)
        if (.not.all(ieee_is_finite(flat))) return
        work=cmplx(flat,0.0_dp,dp)
        allocate(res(size(y)))
        call affine_samples(work,shape(x),shape(y),matrix,offset,res,degree,boundary, &
                            do_filter,fill,code)
        if (code==spline_status_success) then
            yflat=real(res,dp)
            call scatter_real(yflat,y)
        end if
        if (present(status)) status=code
    end subroutine affine_real

    subroutine affine_complex(x,matrix,offset,y,order,status,mode,prefilter,cval)
        complex(dp), intent(in) :: x(..) !! rank 1..15 input
        real(dp), intent(in) :: matrix(:,:) !! square affine transform
        real(dp), intent(in) :: offset(:) !! 0-based offset into input
        complex(dp), intent(out) :: y(..) !! caller-allocated destination
        integer, optional, intent(in) :: order !! interpolation order, default 3
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary extension
        logical, optional, intent(in) :: prefilter !! default false
        complex(dp), optional, intent(in) :: cval !! outside value
        complex(dp), allocatable :: work(:), res(:)
        complex(dp) :: fill
        character(len=:), allocatable :: boundary
        logical :: do_filter
        integer :: degree,code

        code=spline_status_invalid_input
        if (present(status)) status=code
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)/=rank(x)) return
        if (any(shape(x)<1) .or. any(shape(y)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64) .or. &
            size(y,kind=int64)>int(huge(1),int64)) return
        degree=3
        if (present(order)) degree=order
        boundary='mirror'
        if (present(mode)) boundary=trim(mode)
        do_filter=.false.
        if (present(prefilter)) do_filter=prefilter
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cval
        call pack_complex(x,work)
        if (any(.not.ieee_is_finite(real(work,dp))) .or. &
            any(.not.ieee_is_finite(aimag(work)))) return
        allocate(res(size(y)))
        call affine_samples(work,shape(x),shape(y),matrix,offset,res,degree,boundary, &
                            do_filter,fill,code)
        if (code==spline_status_success) call scatter_complex(res,y)
        if (present(status)) status=code
    end subroutine affine_complex

    subroutine shift_real(x,displacement,y,order,status,mode,prefilter,cval)
        real(dp), intent(in) :: x(..) !! sampled real array
        real(dp), intent(in) :: displacement(:) !! shifts in samples along each axis
        real(dp), intent(out) :: y(..) !! shifted real array
        integer, optional, intent(in) :: order !! order 0 through 5
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary extension
        logical, optional, intent(in) :: prefilter !! default false
        real(dp), optional, intent(in) :: cval !! constant extension
        real(dp), allocatable :: identity(:,:)
        integer :: n,k

        n=rank(x)
        allocate(identity(n,n))
        identity=0.0_dp
        do k=1,n
            identity(k,k)=1.0_dp
        end do
        call affine_real(x,identity,-displacement,y,order,status,mode,prefilter,cval)
    end subroutine shift_real

    subroutine shift_complex(x,displacement,y,order,status,mode,prefilter,cval)
        complex(dp), intent(in) :: x(..) !! sampled complex array
        real(dp), intent(in) :: displacement(:) !! shifts in samples along each axis
        complex(dp), intent(out) :: y(..) !! shifted complex array
        integer, optional, intent(in) :: order !! order 0 through 5
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary extension
        logical, optional, intent(in) :: prefilter !! default false
        complex(dp), optional, intent(in) :: cval !! constant extension
        real(dp), allocatable :: identity(:,:)
        integer :: n,k

        n=rank(x)
        allocate(identity(n,n))
        identity=0.0_dp
        do k=1,n
            identity(k,k)=1.0_dp
        end do
        call affine_complex(x,identity,-displacement,y,order,status,mode,prefilter,cval)
    end subroutine shift_complex


    ! Shape helpers follow NumPy/Python's ties-to-even rounding, unlike NINT.
    pure integer function nearest_even_positive(x) result(y)
        real(dp), intent(in) :: x !! nonnegative finite extent
        integer :: lower
        real(dp) :: delta
        lower=int(x)
        delta=x-real(lower,dp)
        y=lower
        if (delta>0.5_dp) y=lower+1
        if (delta==0.5_dp .and. modulo(lower,2)/=0) y=lower+1
    end function nearest_even_positive

    function zoom_shape_vector(input_shape,factors,status) result(out_shape)
        integer, intent(in) :: input_shape(:) !! extents, 1 to 15 axes
        real(dp), intent(in) :: factors(:) !! positive zoom factors, per axis
        integer, optional, intent(out) :: status !! success or invalid input
        integer, allocatable :: out_shape(:)
        real(dp) :: scaled
        integer :: d,code

        allocate(out_shape(size(input_shape)))
        out_shape=0
        code=spline_status_invalid_input
        if (present(status)) status=code
        if (size(input_shape)<1 .or. size(input_shape)>15) return
        if (size(factors)/=size(input_shape)) return
        if (any(input_shape<1) .or. .not.all(ieee_is_finite(factors))) return
        if (any(factors<=0.0_dp)) return
        do d=1,size(input_shape)
            scaled=real(input_shape(d),dp)*factors(d)
            if (.not.ieee_is_finite(scaled)) return
            if (scaled<0.5_dp .or. scaled>real(huge(1)-1,dp)) return
            out_shape(d)=nearest_even_positive(scaled)
            if (out_shape(d)<1) then
                out_shape=0
                return
            end if
        end do
        if (present(status)) status=spline_status_success
    end function zoom_shape_vector

    function zoom_shape_scalar(input_shape,factor,status) result(out_shape)
        integer, intent(in) :: input_shape(:) !! extents, 1 to 15 axes
        real(dp), intent(in) :: factor !! common zoom factor
        integer, optional, intent(out) :: status !! zero on success
        integer, allocatable :: out_shape(:)
        out_shape=zoom_shape_vector(input_shape,spread(factor,1,size(input_shape)),status)
    end function zoom_shape_scalar

    subroutine zoom_matrix(input_shape,output_shape,grid,matrix,offset)
        integer, intent(in) :: input_shape(:) !! input array shape
        integer, intent(in) :: output_shape(:) !! output array shape
        logical, intent(in) :: grid !! count full pixel extent, not centers
        real(dp), intent(out) :: matrix(:,:) !! output-to-input mapping
        real(dp), intent(out) :: offset(:) !! zero-based mapping origin
        real(dp) :: spacing
        integer :: d,n

        n=size(input_shape)
        matrix=0.0_dp
        offset=0.0_dp
        do d=1,n
            if (grid) then
                spacing=real(input_shape(d),dp)/real(output_shape(d),dp)
                offset(d)=(spacing-1.0_dp)*0.5_dp
            else if (output_shape(d)>1) then
                spacing=real(input_shape(d)-1,dp)/real(output_shape(d)-1,dp)
            else
                spacing=1.0_dp
            end if
            matrix(d,d)=spacing
        end do
    end subroutine zoom_matrix

    subroutine zoom_real(x,factors,y,order,status,mode,prefilter,cval,grid_mode)
        real(dp), intent(in) :: x(..) !! real input, rank 1..15
        real(dp), intent(in) :: factors(:) !! positive factors for each axis
        real(dp), intent(out) :: y(..) !! caller-allocated zoomed array
        integer, optional, intent(in) :: order !! spline order, 0..5
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary mode
        logical, optional, intent(in) :: prefilter !! default false
        real(dp), optional, intent(in) :: cval !! constant boundary value
        logical, optional, intent(in) :: grid_mode !! default false
        real(dp), allocatable :: matrix(:,:),offset(:)
        integer, allocatable :: wanted(:)
        integer :: n,code
        logical :: grid

        code=spline_status_invalid_input
        if (present(status)) status=code
        n=rank(x)
        if (n<1 .or. n>15 .or. rank(y)/=n) return
        wanted=zoom_shape_vector(shape(x),factors,code)
        if (code/=spline_status_success) return
        if (any(shape(y)/=wanted)) return
        grid=.false.
        if (present(grid_mode)) grid=grid_mode
        allocate(matrix(n,n),offset(n))
        call zoom_matrix(shape(x),wanted,grid,matrix,offset)
        call affine_real(x,matrix,offset,y,order,code,mode,prefilter,cval)
        if (present(status)) status=code
    end subroutine zoom_real

    subroutine zoom_complex(x,factors,y,order,status,mode,prefilter,cval,grid_mode)
        complex(dp), intent(in) :: x(..) !! complex input, rank 1..15
        real(dp), intent(in) :: factors(:) !! positive factors for each axis
        complex(dp), intent(out) :: y(..) !! caller-allocated zoomed array
        integer, optional, intent(in) :: order !! spline order, 0..5
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary mode
        logical, optional, intent(in) :: prefilter !! default false
        complex(dp), optional, intent(in) :: cval !! constant boundary value
        logical, optional, intent(in) :: grid_mode !! default false
        real(dp), allocatable :: matrix(:,:),offset(:)
        integer, allocatable :: wanted(:)
        integer :: n,code
        logical :: grid

        code=spline_status_invalid_input
        if (present(status)) status=code
        n=rank(x)
        if (n<1 .or. n>15 .or. rank(y)/=n) return
        wanted=zoom_shape_vector(shape(x),factors,code)
        if (code/=spline_status_success) return
        if (any(shape(y)/=wanted)) return
        grid=.false.
        if (present(grid_mode)) grid=grid_mode
        allocate(matrix(n,n),offset(n))
        call zoom_matrix(shape(x),wanted,grid,matrix,offset)
        call affine_complex(x,matrix,offset,y,order,code,mode,prefilter,cval)
        if (present(status)) status=code
    end subroutine zoom_complex

    subroutine zoom_real_scalar(x,factor,y,order,status,mode,prefilter,cval,grid_mode)
        real(dp), intent(in) :: x(..) !! real input
        real(dp), intent(in) :: factor !! common factor for all axes
        real(dp), intent(out) :: y(..) !! caller-allocated output
        integer, optional, intent(in) :: order !! spline order 0..5
        integer, optional, intent(out) :: status !! success status
        character(len=*), optional, intent(in) :: mode !! extension mode
        logical, optional, intent(in) :: prefilter !! default false
        real(dp), optional, intent(in) :: cval !! constant extension
        logical, optional, intent(in) :: grid_mode !! pixel extent mapping
        real(dp), allocatable :: factors(:)
        factors=spread(factor,1,rank(x))
        call zoom_real(x,factors,y,order,status,mode,prefilter,cval,grid_mode)
    end subroutine zoom_real_scalar

    subroutine zoom_complex_scalar(x,factor,y,order,status,mode,prefilter,cval,grid_mode)
        complex(dp), intent(in) :: x(..) !! complex input
        real(dp), intent(in) :: factor !! common factor for all axes
        complex(dp), intent(out) :: y(..) !! caller-allocated output
        integer, optional, intent(in) :: order !! spline order 0..5
        integer, optional, intent(out) :: status !! success status
        character(len=*), optional, intent(in) :: mode !! extension mode
        logical, optional, intent(in) :: prefilter !! default false
        complex(dp), optional, intent(in) :: cval !! constant extension
        logical, optional, intent(in) :: grid_mode !! pixel extent mapping
        real(dp), allocatable :: factors(:)
        factors=spread(factor,1,rank(x))
        call zoom_complex(x,factors,y,order,status,mode,prefilter,cval,grid_mode)
    end subroutine zoom_complex_scalar

    pure subroutine rotation_trig(angle,c,s)
        real(dp), intent(in) :: angle !! degrees
        real(dp), intent(out) :: c,s !! cosine and sine
        real(dp) :: a
        integer :: q
        a=modulo(angle,360.0_dp)
        ! SciPy uses degree-aware trig to avoid almost-zero entries
        ! at right angles, which otherwise affect shape rounding.
        if (modulo(a,90.0_dp)==0.0_dp) then
            q=modulo(nint(a/90.0_dp),4)
            select case (q)
            case (0)
                c=1.0_dp
                s=0.0_dp
            case (1)
                c=0.0_dp
                s=1.0_dp
            case (2)
                c=-1.0_dp
                s=0.0_dp
            case default
                c=0.0_dp
                s=-1.0_dp
            end select
        else
            c=cos(a*acos(-1.0_dp)/180.0_dp)
            s=sin(a*acos(-1.0_dp)/180.0_dp)
        end if
    end subroutine rotation_trig

    function rotate_shape_nd(input_shape,angle,axes,reshape,status) result(out_shape)
        integer, intent(in) :: input_shape(:) !! extents of rank >= 2
        real(dp), intent(in) :: angle !! rotation angle in degrees
        integer, optional, intent(in) :: axes(:) !! two one-based axes
        logical, optional, intent(in) :: reshape !! default true
        integer, optional, intent(out) :: status !! zero on success
        integer, allocatable :: out_shape(:)
        integer :: a,b,n,code
        real(dp) :: c,s,extent_a,extent_b
        logical :: grow

        n=size(input_shape)
        allocate(out_shape(n))
        out_shape=0
        code=spline_status_invalid_input
        if (present(status)) status=code
        if (n<2 .or. n>15) return
        if (any(input_shape<1) .or. .not.ieee_is_finite(angle)) return
        a=1
        b=2
        if (present(axes)) then
            if (size(axes)/=2) return
            a=axes(1)
            b=axes(2)
        end if
        if (a<1 .or. b<1 .or. a>n .or. b>n .or. a==b) return
        out_shape=input_shape
        grow=.true.
        if (present(reshape)) grow=reshape
        if (grow) then
            call rotation_trig(angle,c,s)
            extent_a=abs(c)*real(input_shape(a),dp)+abs(s)*real(input_shape(b),dp)
            extent_b=abs(s)*real(input_shape(a),dp)+abs(c)*real(input_shape(b),dp)
            if (max(extent_a,extent_b)>real(huge(1)-1,dp)) then
                out_shape=0
                return
            end if
            out_shape(a)=int(extent_a+0.5_dp)
            out_shape(b)=int(extent_b+0.5_dp)
            if (any(out_shape<1)) then
                out_shape=0
                return
            end if
        end if
        if (present(status)) status=spline_status_success
    end function rotate_shape_nd

    subroutine rotate_transform(input_shape,output_shape,angle,axes,matrix,offset)
        integer, intent(in) :: input_shape(:) !! original shape
        integer, intent(in) :: output_shape(:) !! rotated shape
        real(dp), intent(in) :: angle !! degrees
        integer, intent(in) :: axes(2) !! sorted one-based axes
        real(dp), intent(out) :: matrix(:,:) !! inverse rotation matrix
        real(dp), intent(out) :: offset(:) !! zero-based input offset
        real(dp) :: c,s,cen_a,cen_b
        integer :: d,a,b

        a=minval(axes)
        b=maxval(axes)
        call rotation_trig(angle,c,s)
        matrix=0.0_dp
        do d=1,size(input_shape)
            matrix(d,d)=1.0_dp
        end do
        matrix(a,a)=c
        matrix(a,b)=s
        matrix(b,a)=-s
        matrix(b,b)=c
        offset=0.0_dp
        cen_a=0.5_dp*real(output_shape(a)-1,dp)
        cen_b=0.5_dp*real(output_shape(b)-1,dp)
        offset(a)=0.5_dp*real(input_shape(a)-1,dp)-c*cen_a-s*cen_b
        offset(b)=0.5_dp*real(input_shape(b)-1,dp)+s*cen_a-c*cen_b
    end subroutine rotate_transform

    subroutine rotate_real(x,angle,y,axes,reshape,order,status,mode,prefilter,cval)
        real(dp), intent(in) :: x(..) !! real input rank 2..15
        real(dp), intent(in) :: angle !! rotation in degrees
        real(dp), intent(out) :: y(..) !! rotated output of requested shape
        integer, optional, intent(in) :: axes(:) !! two one-based axes, default [1,2]
        logical, optional, intent(in) :: reshape !! default true
        integer, optional, intent(in) :: order !! spline order 0..5
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary extension
        logical, optional, intent(in) :: prefilter !! default false
        real(dp), optional, intent(in) :: cval !! constant fill
        integer, allocatable :: wanted(:)
        real(dp), allocatable :: matrix(:,:),offset(:)
        integer :: n,ab(2),code

        code=spline_status_invalid_input
        if (present(status)) status=code
        n=rank(x)
        if (n<2 .or. n>15 .or. rank(y)/=n) return
        wanted=rotate_shape_nd(shape(x),angle,axes,reshape,code)
        if (code/=spline_status_success) return
        if (any(shape(y)/=wanted)) return
        ab=[1,2]
        if (present(axes)) ab=axes
        allocate(matrix(n,n),offset(n))
        call rotate_transform(shape(x),wanted,angle,ab,matrix,offset)
        call affine_real(x,matrix,offset,y,order,code,mode,prefilter,cval)
        if (present(status)) status=code
    end subroutine rotate_real

    subroutine rotate_complex(x,angle,y,axes,reshape,order,status,mode,prefilter,cval)
        complex(dp), intent(in) :: x(..) !! complex input rank 2..15
        real(dp), intent(in) :: angle !! rotation in degrees
        complex(dp), intent(out) :: y(..) !! rotated output of requested shape
        integer, optional, intent(in) :: axes(:) !! two one-based axes, default [1,2]
        logical, optional, intent(in) :: reshape !! default true
        integer, optional, intent(in) :: order !! spline order 0..5
        integer, optional, intent(out) :: status !! zero on success
        character(len=*), optional, intent(in) :: mode !! boundary extension
        logical, optional, intent(in) :: prefilter !! default false
        complex(dp), optional, intent(in) :: cval !! constant fill
        integer, allocatable :: wanted(:)
        real(dp), allocatable :: matrix(:,:),offset(:)
        integer :: n,ab(2),code

        code=spline_status_invalid_input
        if (present(status)) status=code
        n=rank(x)
        if (n<2 .or. n>15 .or. rank(y)/=n) return
        wanted=rotate_shape_nd(shape(x),angle,axes,reshape,code)
        if (code/=spline_status_success) return
        if (any(shape(y)/=wanted)) return
        ab=[1,2]
        if (present(axes)) ab=axes
        allocate(matrix(n,n),offset(n))
        call rotate_transform(shape(x),wanted,angle,ab,matrix,offset)
        call affine_complex(x,matrix,offset,y,order,code,mode,prefilter,cval)
        if (present(status)) status=code
    end subroutine rotate_complex

    ! Build the full coordinate map in C-order, as SciPy invokes mapping.
    ! Coordinates remain one-based at the Fortran callback boundary.
    subroutine geometric_coordinate_table(input_rank,out_dims,mapping,coords,status)
        integer, intent(in) :: input_rank !! rank of source image, 1 through 15
        integer, intent(in) :: out_dims(:) !! positive extents of destination array
        procedure(ndimage_coordinate_map) :: mapping !! coordinate mapping callback
        real(dp), allocatable, intent(out) :: coords(:,:) !! source locations for output pixels
        integer, intent(out) :: status !! 0 success, 1 invalid, 2 resource limit
        integer(int64) :: n,cpos,rem,flat,stride(size(out_dims))
        integer :: nd,d,component,output_coords(size(out_dims))

        status=spline_status_invalid_input
        nd=size(out_dims)
        if (input_rank<1 .or. input_rank>15 .or. nd<1 .or. nd>15) return
        if (any(out_dims<1)) return
        n=1_int64
        do d=1,nd
            ! Check before multiplication, avoiding overflow for huge shapes.
            if (n>max_geometric_coordinates/int(input_rank,int64) / &
                int(out_dims(d),int64)) then
                status=spline_nd_status_too_expensive
                return
            end if
            n=n*int(out_dims(d),int64)
        end do
        allocate(coords(input_rank,int(n)))
        stride(1)=1_int64
        do d=2,nd
            stride(d)=stride(d-1)*int(out_dims(d-1),int64)
        end do
        do cpos=0_int64,n-1_int64
            rem=cpos
            flat=1_int64
            do d=nd,1,-1
                component=int(modulo(rem,int(out_dims(d),int64)))
                rem=rem/int(out_dims(d),int64)
                output_coords(d)=component+1
                flat=flat+int(component,int64)*stride(d)
            end do
            call mapping(output_coords,coords(:,flat))
        end do
        if (.not.all(ieee_is_finite(coords))) return
        status=spline_status_success
    end subroutine geometric_coordinate_table

    ! SciPy geometric_transform equivalent, with one-based callback coordinates.
    ! Destination and source ranks may differ, matching SciPy's output_shape.
    subroutine geometric_real(x,mapping,y,order,status,mode,prefilter,cval)
        real(dp), intent(in) :: x(..) !! real source image, ranks 1 through 15
        procedure(ndimage_coordinate_map) :: mapping !! output-to-input coordinate callback
        real(dp), intent(out) :: y(..) !! caller-shaped destination, ranks 1 through 15
        integer, optional, intent(in) :: order !! interpolation order 0..5, default 3
        integer, optional, intent(out) :: status !! 0 success, 1 invalid, 2 resource limit
        character(len=*), optional, intent(in) :: mode !! boundary mode, default constant
        logical, optional, intent(in) :: prefilter !! apply spline prefilter; default true
        real(dp), optional, intent(in) :: cval !! outside constant value, default zero
        real(dp), allocatable :: coords(:,:),values(:)
        character(len=:), allocatable :: boundary
        integer :: degree,code
        real(dp) :: fill
        logical :: do_prefilter

        code=spline_status_invalid_input
        if (present(status)) status=code
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)<1 .or. rank(y)>15) return
        if (any(shape(x)<1) .or. any(shape(y)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        degree=3
        if (present(order)) degree=order
        if (degree<0 .or. degree>5) return
        boundary='constant'
        if (present(mode)) boundary=trim(mode)
        if (.not.valid_boundary(boundary)) return
        do_prefilter=.true.
        if (present(prefilter)) do_prefilter=prefilter
        fill=0.0_dp
        if (present(cval)) fill=cval
        if (.not.ieee_is_finite(fill)) return
        call geometric_coordinate_table(rank(x),shape(y),mapping,coords,code)
        if (code/=spline_status_success) then
            if (present(status)) status=code
            return
        end if
        allocate(values(size(y)))
        call map_coordinates_nd(x,coords,values,order=degree,status=code, &
            mode=boundary,prefilter=do_prefilter,cval=fill)
        if (code==spline_status_success) call scatter_real(values,y)
        if (present(status)) status=code
    end subroutine geometric_real

    subroutine geometric_complex(x,mapping,y,order,status,mode,prefilter,cval)
        complex(dp), intent(in) :: x(..) !! complex source image, ranks 1 through 15
        procedure(ndimage_coordinate_map) :: mapping !! output-to-input coordinate callback
        complex(dp), intent(out) :: y(..) !! caller-shaped complex destination
        integer, optional, intent(in) :: order !! interpolation order 0..5, default 3
        integer, optional, intent(out) :: status !! 0 success, 1 invalid, 2 resource limit
        character(len=*), optional, intent(in) :: mode !! boundary mode, default constant
        logical, optional, intent(in) :: prefilter !! apply spline prefilter; default true
        complex(dp), optional, intent(in) :: cval !! exterior fill, default zero
        real(dp), allocatable :: coords(:,:)
        complex(dp), allocatable :: values(:)
        character(len=:), allocatable :: boundary
        integer :: degree,code
        complex(dp) :: fill
        logical :: do_prefilter

        code=spline_status_invalid_input
        if (present(status)) status=code
        if (rank(x)<1 .or. rank(x)>15 .or. rank(y)<1 .or. rank(y)>15) return
        if (any(shape(x)<1) .or. any(shape(y)<1)) return
        if (size(x,kind=int64)>int(huge(1),int64)) return
        degree=3
        if (present(order)) degree=order
        if (degree<0 .or. degree>5) return
        boundary='constant'
        if (present(mode)) boundary=trim(mode)
        if (.not.valid_boundary(boundary)) return
        do_prefilter=.true.
        if (present(prefilter)) do_prefilter=prefilter
        fill=(0.0_dp,0.0_dp)
        if (present(cval)) fill=cval
        if (.not.ieee_is_finite(real(fill,dp)) .or. .not.ieee_is_finite(aimag(fill))) return
        call geometric_coordinate_table(rank(x),shape(y),mapping,coords,code)
        if (code/=spline_status_success) then
            if (present(status)) status=code
            return
        end if
        allocate(values(size(y)))
        call map_coordinates_nd(x,coords,values,order=degree,status=code, &
            mode=boundary,prefilter=do_prefilter,cval=fill)
        if (code==spline_status_success) call scatter_complex(values,y)
        if (present(status)) status=code
    end subroutine geometric_complex

end module scifort_signal_spline_nd
