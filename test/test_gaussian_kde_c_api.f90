! SPDX-License-Identifier: MIT
program test_gaussian_kde_c_api
    use, intrinsic :: iso_c_binding, only : c_double, c_int, c_ptr, c_null_ptr, c_associated
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
    use scifort_gaussian_kde, only : gaussian_kde, gaussian_kde_init
    use scifort_c_api_generated_gaussian_kde
    implicit none
    type(gaussian_kde) :: native
    type(c_ptr) :: handle
    real(c_double) :: data(2,6), weights(6), points(2,3), actual(3), empty(0)
    real(c_double) :: neff, factor, scott, silverman
    real(c_double), allocatable :: expected(:)
    integer(c_int) :: status, d, n, mode
    integer :: info

    data(1,:) = [-2.0_c_double,-0.5_c_double,0.25_c_double,1.5_c_double,3.0_c_double,0.8_c_double]
    data(2,:) = [0.5_c_double,2.0_c_double,-1.0_c_double,0.2_c_double,2.8_c_double,-0.4_c_double]
    weights = [1.0_c_double,2.0_c_double,0.0_c_double,3.0_c_double,1.0_c_double,2.0_c_double]
    points(:,1) = [0.0_c_double,0.5_c_double]
    points(:,2) = [-1.0_c_double,1.5_c_double]
    points(:,3) = [2.0_c_double,-0.5_c_double]
    do mode = 0, 2
        select case (mode)
        case (0)
            call gaussian_kde_init(native,data,info,weights=weights,bw_method='scott')
        case (1)
            call gaussian_kde_init(native,data,info,weights=weights,bw_method='silverman')
        case (2)
            call gaussian_kde_init(native,data,info,weights=weights,bw_factor=0.7_c_double)
        end select
        if (info /= 0) stop 1
        handle = c_null_ptr
        call scifort_gaussian_kde_create_f64(2_c_int,6_c_int,data,weights,mode,0.7_c_double,handle,status)
        if (status /= 0 .or. .not. c_associated(handle)) stop 2
        call scifort_gaussian_kde_info_f64(handle,d,n,neff,factor,scott,silverman,status)
        if (status /= 0 .or. d /= 2 .or. n /= 6) stop 3
        if (neff /= native%neff .or. factor /= native%factor) stop 4
        if (scott /= native%scotts_factor() .or. silverman /= native%silverman_factor()) stop 5
        expected = native%pdf(points,info)
        call scifort_gaussian_kde_pdf_f64(2_c_int,3_c_int,handle,points,actual,status)
        if (status /= 0 .or. any(actual /= expected)) stop 6
        expected = native%logpdf(points,info)
        call scifort_gaussian_kde_logpdf_f64(2_c_int,3_c_int,handle,points,actual,status)
        if (status /= 0 .or. any(actual /= expected)) stop 7
        call scifort_gaussian_kde_pdf_f64(2_c_int,0_c_int,handle,points(:,1:0),empty,status)
        if (status /= 0) stop 8
        call scifort_gaussian_kde_pdf_f64(1_c_int,3_c_int,handle,points(1:1,:),actual,status)
        if (status /= 1 .or. .not. all(ieee_is_nan(actual))) stop 9
        call scifort_gaussian_kde_destroy_f64(handle,status)
        if (status /= 0 .or. c_associated(handle)) stop 10
        call scifort_gaussian_kde_destroy_f64(handle,status)
        if (status /= 0 .or. c_associated(handle)) stop 11
    end do
    print '(a)', 'Gaussian KDE native/C handle checks passed.'
end program test_gaussian_kde_c_api
