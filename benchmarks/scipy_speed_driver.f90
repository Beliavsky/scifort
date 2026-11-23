! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program scipy_speed_driver
    use, intrinsic :: iso_fortran_env, only : int64, output_unit
    use scifort_kinds, only : dp
    use scifort_normal, only : normal_cdf, normal_ppf
    use scifort_signal_spectral, only : fft_plan
    use scifort_linalg_dense, only : lu_factor, lu_solve
    use scifort_interpolate_barycentric, only : barycentric_interpolator_type, &
        barycentric_interpolator
    implicit none
    character(40) :: mode, argument
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

    select case (trim(mode))
    case ('normal_cdf_array', 'normal_cdf_scalar', 'normal_ppf_array', 'barycentric')
        allocate(inputs(n,2))
        do k=1,2
            do i=1,n
                t=(real(i-1,dp)+0.25_dp*real(k-1,dp))/real(n,dp)
                select case (trim(mode))
                case ('normal_ppf_array')
                    inputs(i,k)=1e-6_dp+(1.0_dp-2e-6_dp)*t
                case ('barycentric')
                    inputs(i,k)=-0.99_dp+1.98_dp*t
                case default
                    inputs(i,k)=-8.0_dp+16.0_dp*t
                end select
            end do
        end do
        if (trim(mode) == 'barycentric') then
            allocate(nodes(33),samples(33))
            do i=1,33
                nodes(i)=cos(acos(-1.0_dp)*real(i-1,dp)/32.0_dp)
            end do
            samples=1.0_dp/(1.0_dp+25.0_dp*nodes**2)
            polynomial=barycentric_interpolator(nodes,samples,status)
            if (status /= 0) error stop 'interpolator setup failed'
        end if
    case ('fft_reuse', 'fft_setup')
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
    case ('lu_factor_solve', 'lu_reuse')
        allocate(matrix(n,n),lu(n,n),rhs(n,2),pivots(n))
        do j=1,n
            do i=1,n
                matrix(i,j)=1.0_dp/real(i+j,dp)
            end do
            matrix(j,j)=matrix(j,j)+real(n,dp)
        end do
        do k=1,2
            do i=1,n
                values(i)=sin(real(i,dp))+0.1_dp*real(k-1,dp)
            end do
            rhs(:,k)=matmul(matrix,values)
        end do
        call lu_factor(matrix,lu,pivots,status)
        if (status /= 0) error stop 'LU setup failed'
    case default
        error stop 'unknown benchmark'
    end select

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
                    checksum=checksum+real(spectrum(1)+spectrum(n),dp)
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
        integer :: index, code
        select case (trim(mode))
        case ('normal_cdf_array')
            values=normal_cdf(inputs(:,which))
        case ('normal_cdf_scalar')
            do index=1,n
                values(index)=normal_cdf(inputs(index,which))
            end do
        case ('normal_ppf_array')
            values=normal_ppf(inputs(:,which))
        case ('fft_reuse','fft_setup')
            if (trim(mode) == 'fft_setup') then
                call plan%initialize(n,code)
                if (code /= 0) error stop 'FFT setup failed'
            end if
            spectrum=signals(:,which)
            call plan%execute(spectrum,code)
            if (code /= 0) error stop 'FFT execution failed'
        case ('lu_factor_solve','lu_reuse')
            if (trim(mode) == 'lu_factor_solve') then
                call lu_factor(matrix,lu,pivots,code)
                if (code /= 0) error stop 'LU factorization failed'
            end if
            call lu_solve(lu,pivots,rhs(:,which),values,code)
            if (code /= 0) error stop 'LU solve failed'
        case ('barycentric')
            interpolated=polynomial%values(inputs(:,which))
            values=interpolated(:,1)
        end select
    end subroutine calculate
end program scipy_speed_driver
