! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Machine-readable driver for independent randomized SciPy comparisons.
program tensor_grid_reference_driver
    use scifort_kinds, only : dp
    use scifort_interpolate, only : grid_axis, regular_grid_interpolator_type, &
        complex_regular_grid_interpolator_type, regular_grid_interpolator, complex_regular_grid_interpolator
    implicit none
    type(grid_axis), allocatable :: axes(:)
    type(regular_grid_interpolator_type) :: grid
    type(complex_regular_grid_interpolator_type) :: cgrid
    real(dp), allocatable :: r(:), im(:), x(:,:), yr(:,:)
    complex(dp), allocatable :: yc(:,:), z(:)
    integer, allocatable :: sizes(:), nu(:)
    integer :: nc, nd, nq, method_code, bounds_code, fill_code, complex_code
    integer :: case_id, cases, i, j, k, np, code, setup_code
    character(len=7) :: method
    logical :: bounds

    read(*,*) cases
    do case_id = 1, cases
        read(*,*) nd, nc, nq, method_code, bounds_code, fill_code, complex_code
        allocate(axes(nd),sizes(nd),nu(nd),x(nq,nd))
        np = 1
        do k = 1, nd
            read(*,*) sizes(k)
            allocate(axes(k)%points(sizes(k)))
            read(*,*) axes(k)%points
            np = np*sizes(k)
        end do
        read(*,*) nu
        allocate(r(np*nc),im(np*nc))
        read(*,*) r
        read(*,*) im
        do i = 1, nq
            read(*,*) x(i,:)
        end do
        method = 'linear'
        select case(method_code)
        case (1)
            method = 'nearest'
        case (2)
            method = 'slinear'
        case (3)
            method = 'cubic'
        case (4)
            method = 'quintic'
        case (5)
            method = 'pchip'
        end select
        bounds = bounds_code /= 0
        if (complex_code == 0) then
            select case(nd)
            case(1)
                grid = construct_real(axes,reshape(r,[sizes(1),nc]))
            case(2)
                grid = construct_real(axes,reshape(r,[sizes(1),sizes(2),nc]))
            case(3)
                grid = construct_real(axes,reshape(r,[sizes(1),sizes(2),sizes(3),nc]))
            case(4)
                grid = construct_real(axes,reshape(r,[sizes(1),sizes(2),sizes(3),sizes(4),nc]))
            case(5)
                grid = construct_real(axes,reshape(r,[sizes(1),sizes(2),sizes(3),sizes(4),sizes(5),nc]))
            case(6)
                grid = construct_real(axes,reshape(r,[sizes(1),sizes(2),sizes(3),sizes(4),sizes(5),sizes(6),nc]))
            end select
            if (nu(1) >= 0) then
                yr = grid%evaluate(x,status=code,nu=nu)
            else
                yr = grid%evaluate(x,status=code)
            end if
            write(*,'(a,1x,i0,1x,i0,1x,i0)') 'CASE',case_id,setup_code,code
            do i = 1, nq
                write(*,'(*(es25.16e3,1x))') (yr(i,j),0.0_dp,j=1,nc)
            end do
        else
            z = cmplx(r,im,kind=dp)
            select case(nd)
            case(1)
                cgrid = construct_complex(axes,reshape(z,[sizes(1),nc]))
            case(2)
                cgrid = construct_complex(axes,reshape(z,[sizes(1),sizes(2),nc]))
            case(3)
                cgrid = construct_complex(axes,reshape(z,[sizes(1),sizes(2),sizes(3),nc]))
            case(4)
                cgrid = construct_complex(axes,reshape(z,[sizes(1),sizes(2),sizes(3),sizes(4),nc]))
            case(5)
                cgrid = construct_complex(axes,reshape(z,[sizes(1),sizes(2),sizes(3),sizes(4),sizes(5),nc]))
            case(6)
                cgrid = construct_complex(axes,reshape(z,[sizes(1),sizes(2),sizes(3),sizes(4),sizes(5),sizes(6),nc]))
            end select
            if (nu(1) >= 0) then
                yc = cgrid%evaluate(x,status=code,nu=nu)
            else
                yc = cgrid%evaluate(x,status=code)
            end if
            write(*,'(a,1x,i0,1x,i0,1x,i0)') 'CASE',case_id,setup_code,code
            do i = 1, nq
                write(*,'(*(es25.16e3,1x))') (real(yc(i,j),dp),aimag(yc(i,j)),j=1,nc)
            end do
        end if
        deallocate(axes,sizes,nu,x,r,im)
    end do
contains
    function construct_real(coords,data) result(obj)
        type(grid_axis), intent(in) :: coords(:) !! input grid axes
        real(dp), intent(in) :: data(..) !! gridpoint-major sampled values
        type(regular_grid_interpolator_type) :: obj
        select case(fill_code)
        case(0)
            obj = regular_grid_interpolator(coords,data,method,bounds_error=bounds,status=setup_code)
        case(1)
            obj = regular_grid_interpolator(coords,data,method,bounds_error=bounds, &
                fill_value=-17.375_dp,status=setup_code)
        case(2)
            obj = regular_grid_interpolator(coords,data,method,bounds_error=bounds, &
                extrapolate=.true.,status=setup_code)
        end select
    end function construct_real

    function construct_complex(coords,data) result(obj)
        type(grid_axis), intent(in) :: coords(:) !! input grid axes
        complex(dp), intent(in) :: data(..) !! complex sampled values, optional trailing channels
        type(complex_regular_grid_interpolator_type) :: obj
        select case(fill_code)
        case(0)
            obj = complex_regular_grid_interpolator(coords,data,method,bounds_error=bounds,status=setup_code)
        case(1)
            obj = complex_regular_grid_interpolator(coords,data,method,bounds_error=bounds, &
                fill_value=cmplx(-17.375_dp,3.25_dp,kind=dp),status=setup_code)
        case(2)
            obj = complex_regular_grid_interpolator(coords,data,method,bounds_error=bounds, &
                extrapolate=.true.,status=setup_code)
        end select
    end function construct_complex
end program tensor_grid_reference_driver
