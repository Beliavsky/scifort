! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program scipy_speed_driver
    use, intrinsic :: iso_fortran_env, only : int64, output_unit
    use scifort_kinds, only : dp
    use scifort_normal, only : &
        normal_pdf, normal_logpdf, normal_cdf, &
        normal_sf, normal_logcdf, normal_logsf, &
        normal_ppf, normal_isf
    use scifort_exponential, only : &
        exponential_pdf, exponential_logpdf, exponential_cdf, &
        exponential_sf, exponential_logcdf, exponential_logsf, &
        exponential_ppf, exponential_isf
    use scifort_logistic, only : &
        logistic_pdf, logistic_logpdf, logistic_cdf, &
        logistic_sf, logistic_logcdf, logistic_logsf, &
        logistic_ppf, logistic_isf
    use scifort_laplace, only : &
        laplace_pdf, laplace_logpdf, laplace_cdf, &
        laplace_sf, laplace_logcdf, laplace_logsf, &
        laplace_ppf, laplace_isf
    use scifort_cauchy, only : &
        cauchy_pdf, cauchy_logpdf, cauchy_cdf, &
        cauchy_sf, cauchy_logcdf, cauchy_logsf, &
        cauchy_ppf, cauchy_isf
    use scifort_weibull, only : &
        weibull_pdf, weibull_logpdf, weibull_cdf, &
        weibull_sf, weibull_logcdf, weibull_logsf, &
        weibull_ppf, weibull_isf
    use scifort_lognormal, only : &
        lognormal_pdf, lognormal_logpdf, lognormal_cdf, &
        lognormal_sf, lognormal_logcdf, lognormal_logsf, &
        lognormal_ppf, lognormal_isf
    use scifort_signal_spectral, only : fft_plan
    use scifort_linalg_dense, only : lu_factor, lu_solve, cho_factor, cho_solve, &
        solve_triangular
    use scifort_interpolate_barycentric, only : barycentric_interpolator_type, &
        barycentric_interpolator
    implicit none
    character(40) :: mode, argument, family
    logical :: is_distribution, is_interpolation
    character(4096) :: result_path
    real(dp), allocatable :: inputs(:,:), values(:), matrix(:,:), lu(:,:), rhs(:,:)
    real(dp), allocatable :: nodes(:), samples(:), interpolated(:,:)
    complex(dp), allocatable :: signals(:,:), spectrum(:), output(:,:)
    integer, allocatable :: pivots(:)
    type(fft_plan) :: plan
    type(barycentric_interpolator_type) :: polynomial
    integer :: n, i, j, k, repetitions, repeat, status, unit, ios
    integer(int64) :: first, last, rate
    real(dp) :: t, elapsed, checksum

    call get_command_argument(1, mode)
    call get_command_argument(2, argument)
    read(argument,*,iostat=ios) n
    if (ios /= 0) error stop 'invalid benchmark size'
    if (n < 2) error stop 'benchmark size must be at least two'
    call get_command_argument(3, result_path)
    call system_clock(count_rate=rate)
    if (rate <= 0_int64) error stop 'no wall clock'
    allocate(values(n), output(n,2))

    i=index(trim(mode),'_')
    family=mode(:max(1,i-1))
    is_distribution=index(trim(mode),'_array') > 0 .or. trim(mode) == 'normal_cdf_scalar'
    is_interpolation=index(trim(mode),'barycentric') == 1
    if (is_distribution .or. is_interpolation) then
        allocate(inputs(n,2))
        do k=1,2
            do i=1,n
                t=(real(i-1,dp)+0.25_dp*real(k-1,dp))/real(n,dp)
                if (index(trim(mode),'_ppf_') > 0 .or. index(trim(mode),'_isf_') > 0) then
                    inputs(i,k)=1e-6_dp+(1.0_dp-2e-6_dp)*t
                else if (is_interpolation) then
                    inputs(i,k)=-0.99_dp+1.98_dp*t
                else if (family == 'exponential' .or. family == 'weibull' .or. &
                    family == 'lognormal') then
                    inputs(i,k)=0.01_dp+8.0_dp*t
                else
                    inputs(i,k)=-8.0_dp+16.0_dp*t
                end if
            end do
        end do
        if (is_interpolation) then
            allocate(nodes(33),samples(33))
            do i=1,33
                nodes(i)=cos(acos(-1.0_dp)*real(i-1,dp)/32.0_dp)
            end do
            samples=1.0_dp/(1.0_dp+25.0_dp*nodes**2)
            polynomial=barycentric_interpolator(nodes,samples,status)
            if (status /= 0) error stop 'interpolator setup failed'
        end if
    else
    select case (trim(mode))
    case ('fft_reuse', 'fft_setup', 'ifft_reuse', 'ifft_setup')
        allocate(signals(n,2),spectrum(n))
        do k=1,2
            do i=1,n
                t=real(i-1,dp)/real(n,dp)
                signals(i,k)=cmplx(sin(17.0_dp*t)+0.1_dp*real(k-1,dp), &
                    cos(11.0_dp*t),dp)
            end do
        end do
        call plan%initialize(n,status)
        if (status /= 0) error stop 'FFT setup failed'
    case ('lu_factor_solve', 'lu_reuse', 'lu_transpose', 'cho_factor_solve', &
        'cho_reuse', 'triangular_lower', 'triangular_upper', 'triangular_transpose')
        allocate(matrix(n,n),lu(n,n),rhs(n,2),pivots(n))
        do j=1,n
            do i=1,n
                matrix(i,j)=1.0_dp/real(i+j,dp)
            end do
            matrix(j,j)=matrix(j,j)+real(n,dp)
        end do
        if (index(trim(mode),'triangular') == 1) then
            do j=1,n
                do i=1,n
                    if (trim(mode) == 'triangular_lower') then
                        if (i < j) matrix(i,j)=0.0_dp
                    else
                        if (i > j) matrix(i,j)=0.0_dp
                    end if
                end do
            end do
        end if
        do k=1,2
            do i=1,n
                values(i)=sin(real(i,dp))+0.1_dp*real(k-1,dp)
            end do
            if (trim(mode) == 'lu_transpose' .or. trim(mode) == 'triangular_transpose') then
                rhs(:,k)=matmul(transpose(matrix),values)
            else
                rhs(:,k)=matmul(matrix,values)
            end if
        end do
        if (index(trim(mode),'cho') == 1) then
            call cho_factor(matrix,lu,status,lower=.true.)
        else
            call lu_factor(matrix,lu,pivots,status)
        end if
        if (status /= 0) error stop 'LU setup failed'
    case default
        error stop 'unknown benchmark'
    end select
    end if

    ! The parent starts timing only after setup and full-result validation.
    write(*,'(a,1x,i0)') 'READY',rate
    flush(output_unit)
    do
        read(*,*,iostat=ios) repetitions
        if (ios /= 0) exit
        if (repetitions < 0) exit
        if (repetitions == 0) then
            do k=1,2
                call calculate(k)
                if (allocated(spectrum)) then
                    output(:,k)=spectrum
                else
                    output(:,k)=cmplx(values,0.0_dp,dp)
                end if
            end do
            open(newunit=unit,file=trim(result_path),access='stream',form='unformatted', &
                status='replace',action='write',iostat=ios)
            if (ios /= 0) error stop 'cannot write validation results'
            write(unit) output
            close(unit)
            write(*,'(a)') 'VALIDATED'
        else
            checksum=0.0_dp
            call system_clock(first)
            do repeat=1,repetitions
                k=1+mod(repeat-1,2)
                call calculate(k)
                if (allocated(spectrum)) then
                    checksum=checksum+real(spectrum(1)+spectrum(n),dp)+aimag(spectrum(1)+spectrum(n))
                else
                    checksum=checksum+values(1)+values(n)
                end if
            end do
            call system_clock(last)
            elapsed=real(last-first,dp)/real(rate,dp)
            write(*,'(es25.16e3,1x,es25.16e3)') elapsed,checksum
        end if
        flush(output_unit)
    end do
contains
    subroutine calculate(which)
        implicit none
        integer, intent(in) :: which !! index of one of two prepared input sets, 1 or 2
        integer :: idx, code
        select case (trim(mode))
        case ('normal_cdf_scalar')
            do idx=1,n
                values(idx)=normal_cdf(inputs(idx,which))
            end do
        case ('normal_pdf_array')
            values=normal_pdf(inputs(:,which))
        case ('normal_logpdf_array')
            values=normal_logpdf(inputs(:,which))
        case ('normal_cdf_array')
            values=normal_cdf(inputs(:,which))
        case ('normal_sf_array')
            values=normal_sf(inputs(:,which))
        case ('normal_logcdf_array')
            values=normal_logcdf(inputs(:,which))
        case ('normal_logsf_array')
            values=normal_logsf(inputs(:,which))
        case ('normal_ppf_array')
            values=normal_ppf(inputs(:,which))
        case ('normal_isf_array')
            values=normal_isf(inputs(:,which))
        case ('exponential_pdf_array')
            values=exponential_pdf(inputs(:,which))
        case ('exponential_logpdf_array')
            values=exponential_logpdf(inputs(:,which))
        case ('exponential_cdf_array')
            values=exponential_cdf(inputs(:,which))
        case ('exponential_sf_array')
            values=exponential_sf(inputs(:,which))
        case ('exponential_logcdf_array')
            values=exponential_logcdf(inputs(:,which))
        case ('exponential_logsf_array')
            values=exponential_logsf(inputs(:,which))
        case ('exponential_ppf_array')
            values=exponential_ppf(inputs(:,which))
        case ('exponential_isf_array')
            values=exponential_isf(inputs(:,which))
        case ('logistic_pdf_array')
            values=logistic_pdf(inputs(:,which))
        case ('logistic_logpdf_array')
            values=logistic_logpdf(inputs(:,which))
        case ('logistic_cdf_array')
            values=logistic_cdf(inputs(:,which))
        case ('logistic_sf_array')
            values=logistic_sf(inputs(:,which))
        case ('logistic_logcdf_array')
            values=logistic_logcdf(inputs(:,which))
        case ('logistic_logsf_array')
            values=logistic_logsf(inputs(:,which))
        case ('logistic_ppf_array')
            values=logistic_ppf(inputs(:,which))
        case ('logistic_isf_array')
            values=logistic_isf(inputs(:,which))
        case ('laplace_pdf_array')
            values=laplace_pdf(inputs(:,which))
        case ('laplace_logpdf_array')
            values=laplace_logpdf(inputs(:,which))
        case ('laplace_cdf_array')
            values=laplace_cdf(inputs(:,which))
        case ('laplace_sf_array')
            values=laplace_sf(inputs(:,which))
        case ('laplace_logcdf_array')
            values=laplace_logcdf(inputs(:,which))
        case ('laplace_logsf_array')
            values=laplace_logsf(inputs(:,which))
        case ('laplace_ppf_array')
            values=laplace_ppf(inputs(:,which))
        case ('laplace_isf_array')
            values=laplace_isf(inputs(:,which))
        case ('cauchy_pdf_array')
            values=cauchy_pdf(inputs(:,which))
        case ('cauchy_logpdf_array')
            values=cauchy_logpdf(inputs(:,which))
        case ('cauchy_cdf_array')
            values=cauchy_cdf(inputs(:,which))
        case ('cauchy_sf_array')
            values=cauchy_sf(inputs(:,which))
        case ('cauchy_logcdf_array')
            values=cauchy_logcdf(inputs(:,which))
        case ('cauchy_logsf_array')
            values=cauchy_logsf(inputs(:,which))
        case ('cauchy_ppf_array')
            values=cauchy_ppf(inputs(:,which))
        case ('cauchy_isf_array')
            values=cauchy_isf(inputs(:,which))
        case ('weibull_pdf_array')
            values=weibull_pdf(inputs(:,which),1.7_dp)
        case ('weibull_logpdf_array')
            values=weibull_logpdf(inputs(:,which),1.7_dp)
        case ('weibull_cdf_array')
            values=weibull_cdf(inputs(:,which),1.7_dp)
        case ('weibull_sf_array')
            values=weibull_sf(inputs(:,which),1.7_dp)
        case ('weibull_logcdf_array')
            values=weibull_logcdf(inputs(:,which),1.7_dp)
        case ('weibull_logsf_array')
            values=weibull_logsf(inputs(:,which),1.7_dp)
        case ('weibull_ppf_array')
            values=weibull_ppf(inputs(:,which),1.7_dp)
        case ('weibull_isf_array')
            values=weibull_isf(inputs(:,which),1.7_dp)
        case ('lognormal_pdf_array')
            values=lognormal_pdf(inputs(:,which),0.8_dp)
        case ('lognormal_logpdf_array')
            values=lognormal_logpdf(inputs(:,which),0.8_dp)
        case ('lognormal_cdf_array')
            values=lognormal_cdf(inputs(:,which),0.8_dp)
        case ('lognormal_sf_array')
            values=lognormal_sf(inputs(:,which),0.8_dp)
        case ('lognormal_logcdf_array')
            values=lognormal_logcdf(inputs(:,which),0.8_dp)
        case ('lognormal_logsf_array')
            values=lognormal_logsf(inputs(:,which),0.8_dp)
        case ('lognormal_ppf_array')
            values=lognormal_ppf(inputs(:,which),0.8_dp)
        case ('lognormal_isf_array')
            values=lognormal_isf(inputs(:,which),0.8_dp)
        case ('fft_reuse','fft_setup','ifft_reuse','ifft_setup')
            if (index(trim(mode),'setup') > 0) then
                call plan%initialize(n,code)
                if (code /= 0) error stop 'FFT setup failed'
            end if
            spectrum=signals(:,which)
            if (index(trim(mode),'ifft') == 1) then
                call plan%inverse(spectrum,code)
            else
                call plan%execute(spectrum,code)
            end if
            if (code /= 0) error stop 'FFT execution failed'
        case ('lu_factor_solve','lu_reuse','lu_transpose')
            if (trim(mode) == 'lu_factor_solve') then
                call lu_factor(matrix,lu,pivots,code)
                if (code /= 0) error stop 'LU factorization failed'
            end if
            if (trim(mode) == 'lu_transpose') then
                call lu_solve(lu,pivots,rhs(:,which),values,code,trans='T')
            else
                call lu_solve(lu,pivots,rhs(:,which),values,code)
            end if
            if (code /= 0) error stop 'LU solve failed'
        case ('cho_factor_solve','cho_reuse')
            if (trim(mode) == 'cho_factor_solve') then
                call cho_factor(matrix,lu,code,lower=.true.)
                if (code /= 0) error stop 'Cholesky factorization failed'
            end if
            call cho_solve(lu,rhs(:,which),values,code,lower=.true.)
            if (code /= 0) error stop 'Cholesky solve failed'
        case ('triangular_lower','triangular_upper','triangular_transpose')
            call solve_triangular(matrix,rhs(:,which),values,code, &
                lower=trim(mode) == 'triangular_lower', &
                trans=merge('T','N',trim(mode) == 'triangular_transpose'))
            if (code /= 0) error stop 'triangular solve failed'
        case ('barycentric','barycentric_setup','barycentric_derivative')
            if (trim(mode) == 'barycentric_setup') then
                polynomial=barycentric_interpolator(nodes,samples,code)
                if (code /= 0) error stop 'interpolator construction failed'
            end if
            interpolated=polynomial%values(inputs(:,which), &
                nu=merge(1,0,trim(mode) == 'barycentric_derivative'))
            values=interpolated(:,1)
        end select
    end subroutine calculate
end program scipy_speed_driver
