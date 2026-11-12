! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Independent structured-matrix reference driver for the Python/SciPy oracle.
program linalg_structure_diagsvd_driver
    use, intrinsic :: iso_fortran_env, only : int32, int64
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success
    use scifort_linalg_structure, only : bandwidth, issymmetric, ishermitian
    use scifort_linalg_diagsvd, only : diagsvd
    implicit none
    integer :: t, ncases, op, dtype, ndim, batch, m, n, ns
    real(dp) :: atol, rtol

    interface call_structure_2d
        procedure call_structure_2d_i32
        procedure call_structure_2d_i64
        procedure call_structure_2d_r
        procedure call_structure_2d_c
    end interface
    interface call_structure_3d
        procedure call_structure_3d_i32
        procedure call_structure_3d_i64
        procedure call_structure_3d_r
        procedure call_structure_3d_c
    end interface

    read(*,*) ncases
    do t = 1, ncases
        read(*,*) op, dtype, ndim, batch, m, n, ns, atol, rtol
        select case(dtype)
        case(1)
            block
                integer(int32), allocatable :: a(:, :, :), singular_values(:, :), sigma(:, :, :)
                integer, allocatable :: low(:), high(:)
                logical, allocatable :: flags(:)
                integer :: b, j, i, status

                allocate(a(batch, m, n), singular_values(batch, ns))
                allocate(sigma(batch, m, n), low(batch), high(batch), flags(batch))
                if (op == 4) then
                    do b = 1, batch
                        do i = 1, ns
                            read(*,*) singular_values(b,i)
                        end do
                    end do
                else
                    do b = 1, batch
                        do j = 1, n
                            do i = 1, m
                                read(*,*) a(b,i,j)
                            end do
                        end do
                    end do
                end if
                select case(op)
                case(1)
                    if (ndim == 2) then
                        call bandwidth(a(1,:,:), low(1), high(1), status)
                    else
                        call bandwidth(a, low, high, status)
                    end if
                case(2,3)
                    if (ndim == 2) then
                        call call_structure_2d(a(1,:,:), op, flags(1), status, atol, rtol)
                    else
                        call call_structure_3d(a, op, flags, status, atol, rtol)
                    end if
                case(4)
                    if (ndim == 2) then
                        call diagsvd(singular_values(1,:), m, n, sigma(1,:,:), status)
                    else
                        call diagsvd(singular_values, m, n, sigma, status)
                    end if
                end select
                write(*,'(i0)') status
                if (status == linalg_status_success) then
                    select case(op)
                    case(1)
                        do b = 1, batch
                            write(*,'(i0,1x,i0)') low(b), high(b)
                        end do
                    case(2,3)
                        do b = 1, batch
                            write(*,'(i0)') merge(1,0,flags(b))
                        end do
                    case(4)
                        do b = 1, batch
                            do j = 1, n
                                do i = 1, m
                                    write(*,'(i0)') sigma(b,i,j)
                                end do
                            end do
                        end do
                    end select
                end if
            end block
        case(2)
            block
                integer(int64), allocatable :: a(:, :, :), singular_values(:, :), sigma(:, :, :)
                integer, allocatable :: low(:), high(:)
                logical, allocatable :: flags(:)
                integer :: b, j, i, status

                allocate(a(batch, m, n), singular_values(batch, ns))
                allocate(sigma(batch, m, n), low(batch), high(batch), flags(batch))
                if (op == 4) then
                    do b = 1, batch
                        do i = 1, ns
                            read(*,*) singular_values(b,i)
                        end do
                    end do
                else
                    do b = 1, batch
                        do j = 1, n
                            do i = 1, m
                                read(*,*) a(b,i,j)
                            end do
                        end do
                    end do
                end if
                select case(op)
                case(1)
                    if (ndim == 2) then
                        call bandwidth(a(1,:,:), low(1), high(1), status)
                    else
                        call bandwidth(a, low, high, status)
                    end if
                case(2,3)
                    if (ndim == 2) then
                        call call_structure_2d(a(1,:,:), op, flags(1), status, atol, rtol)
                    else
                        call call_structure_3d(a, op, flags, status, atol, rtol)
                    end if
                case(4)
                    if (ndim == 2) then
                        call diagsvd(singular_values(1,:), m, n, sigma(1,:,:), status)
                    else
                        call diagsvd(singular_values, m, n, sigma, status)
                    end if
                end select
                write(*,'(i0)') status
                if (status == linalg_status_success) then
                    select case(op)
                    case(1)
                        do b = 1, batch
                            write(*,'(i0,1x,i0)') low(b), high(b)
                        end do
                    case(2,3)
                        do b = 1, batch
                            write(*,'(i0)') merge(1,0,flags(b))
                        end do
                    case(4)
                        do b = 1, batch
                            do j = 1, n
                                do i = 1, m
                                    write(*,'(i0)') sigma(b,i,j)
                                end do
                            end do
                        end do
                    end select
                end if
            end block
        case(3)
            block
                real(dp), allocatable :: a(:, :, :), singular_values(:, :), sigma(:, :, :)
                integer, allocatable :: low(:), high(:)
                logical, allocatable :: flags(:)
                integer :: b, j, i, status

                allocate(a(batch, m, n), singular_values(batch, ns))
                allocate(sigma(batch, m, n), low(batch), high(batch), flags(batch))
                if (op == 4) then
                    do b = 1, batch
                        do i = 1, ns
                            read(*,*) singular_values(b,i)
                        end do
                    end do
                else
                    do b = 1, batch
                        do j = 1, n
                            do i = 1, m
                                read(*,*) a(b,i,j)
                            end do
                        end do
                    end do
                end if
                select case(op)
                case(1)
                    if (ndim == 2) then
                        call bandwidth(a(1,:,:), low(1), high(1), status)
                    else
                        call bandwidth(a, low, high, status)
                    end if
                case(2,3)
                    if (ndim == 2) then
                        call call_structure_2d(a(1,:,:), op, flags(1), status, atol, rtol)
                    else
                        call call_structure_3d(a, op, flags, status, atol, rtol)
                    end if
                case(4)
                    if (ndim == 2) then
                        call diagsvd(singular_values(1,:), m, n, sigma(1,:,:), status)
                    else
                        call diagsvd(singular_values, m, n, sigma, status)
                    end if
                end select
                write(*,'(i0)') status
                if (status == linalg_status_success) then
                    select case(op)
                    case(1)
                        do b = 1, batch
                            write(*,'(i0,1x,i0)') low(b), high(b)
                        end do
                    case(2,3)
                        do b = 1, batch
                            write(*,'(i0)') merge(1,0,flags(b))
                        end do
                    case(4)
                        do b = 1, batch
                            do j = 1, n
                                do i = 1, m
                                    write(*,'(es25.16e3)') sigma(b,i,j)
                                end do
                            end do
                        end do
                    end select
                end if
            end block
        case(4)
            block
                complex(dp), allocatable :: a(:, :, :), singular_values(:, :), sigma(:, :, :)
                integer, allocatable :: low(:), high(:)
                logical, allocatable :: flags(:)
                integer :: b, j, i, status
                real(dp) :: re, im

                allocate(a(batch, m, n), singular_values(batch, ns))
                allocate(sigma(batch, m, n), low(batch), high(batch), flags(batch))
                if (op == 4) then
                    do b = 1, batch
                        do i = 1, ns
                            read(*,*) re, im
                            singular_values(b,i) = cmplx(re, im, dp)
                        end do
                    end do
                else
                    do b = 1, batch
                        do j = 1, n
                            do i = 1, m
                                read(*,*) re, im
                                a(b,i,j) = cmplx(re, im, dp)
                            end do
                        end do
                    end do
                end if
                select case(op)
                case(1)
                    if (ndim == 2) then
                        call bandwidth(a(1,:,:), low(1), high(1), status)
                    else
                        call bandwidth(a, low, high, status)
                    end if
                case(2,3)
                    if (ndim == 2) then
                        call call_structure_2d(a(1,:,:), op, flags(1), status, atol, rtol)
                    else
                        call call_structure_3d(a, op, flags, status, atol, rtol)
                    end if
                case(4)
                    if (ndim == 2) then
                        call diagsvd(singular_values(1,:), m, n, sigma(1,:,:), status)
                    else
                        call diagsvd(singular_values, m, n, sigma, status)
                    end if
                end select
                write(*,'(i0)') status
                if (status == linalg_status_success) then
                    select case(op)
                    case(1)
                        do b = 1, batch
                            write(*,'(i0,1x,i0)') low(b), high(b)
                        end do
                    case(2,3)
                        do b = 1, batch
                            write(*,'(i0)') merge(1,0,flags(b))
                        end do
                    case(4)
                        do b = 1, batch
                            do j = 1, n
                                do i = 1, m
                                    write(*,'(es25.16e3,1x,es25.16e3)') &
                                        real(sigma(b,i,j),dp), aimag(sigma(b,i,j))
                                end do
                            end do
                        end do
                    end select
                end if
            end block
        end select
    end do
contains
    subroutine call_structure_2d_i32(a, op, flags, status, at, rt)
        integer(int32), intent(in) :: a(:,:) !! input matrix
        integer, intent(in) :: op !! 2 for symmetric, 3 for Hermitian
        logical, intent(out) :: flags !! results
        integer, intent(out) :: status !! status
        real(dp), intent(in) :: at !! absolute tolerance or negative for omitted
        real(dp), intent(in) :: rt !! relative tolerance or negative for omitted

        if (op == 2) then
            if (at >= 0.0_dp .and. rt >= 0.0_dp) then
                call issymmetric(a, flags, status, atol=at, rtol=rt)
            else if (at >= 0.0_dp) then
                call issymmetric(a, flags, status, atol=at)
            else if (rt >= 0.0_dp) then
                call issymmetric(a, flags, status, rtol=rt)
            else
                call issymmetric(a, flags, status)
            end if
        else if (op == 3) then
            if (at >= 0.0_dp .and. rt >= 0.0_dp) then
                call ishermitian(a, flags, status, atol=at, rtol=rt)
            else if (at >= 0.0_dp) then
                call ishermitian(a, flags, status, atol=at)
            else if (rt >= 0.0_dp) then
                call ishermitian(a, flags, status, rtol=rt)
            else
                call ishermitian(a, flags, status)
            end if
        end if
    end subroutine call_structure_2d_i32

    subroutine call_structure_3d_i32(a, op, flags, status, at, rt)
        integer(int32), intent(in) :: a(:,:,:) !! input matrix
        integer, intent(in) :: op !! 2 for symmetric, 3 for Hermitian
        logical, intent(out) :: flags(:) !! results
        integer, intent(out) :: status !! status
        real(dp), intent(in) :: at !! absolute tolerance or negative for omitted
        real(dp), intent(in) :: rt !! relative tolerance or negative for omitted

        if (op == 2) then
            if (at >= 0.0_dp .and. rt >= 0.0_dp) then
                call issymmetric(a, flags, status, atol=at, rtol=rt)
            else if (at >= 0.0_dp) then
                call issymmetric(a, flags, status, atol=at)
            else if (rt >= 0.0_dp) then
                call issymmetric(a, flags, status, rtol=rt)
            else
                call issymmetric(a, flags, status)
            end if
        else if (op == 3) then
            if (at >= 0.0_dp .and. rt >= 0.0_dp) then
                call ishermitian(a, flags, status, atol=at, rtol=rt)
            else if (at >= 0.0_dp) then
                call ishermitian(a, flags, status, atol=at)
            else if (rt >= 0.0_dp) then
                call ishermitian(a, flags, status, rtol=rt)
            else
                call ishermitian(a, flags, status)
            end if
        end if
    end subroutine call_structure_3d_i32

    subroutine call_structure_2d_i64(a, op, flags, status, at, rt)
        integer(int64), intent(in) :: a(:,:) !! input matrix
        integer, intent(in) :: op !! 2 for symmetric, 3 for Hermitian
        logical, intent(out) :: flags !! results
        integer, intent(out) :: status !! status
        real(dp), intent(in) :: at !! absolute tolerance or negative for omitted
        real(dp), intent(in) :: rt !! relative tolerance or negative for omitted

        if (op == 2) then
            if (at >= 0.0_dp .and. rt >= 0.0_dp) then
                call issymmetric(a, flags, status, atol=at, rtol=rt)
            else if (at >= 0.0_dp) then
                call issymmetric(a, flags, status, atol=at)
            else if (rt >= 0.0_dp) then
                call issymmetric(a, flags, status, rtol=rt)
            else
                call issymmetric(a, flags, status)
            end if
        else if (op == 3) then
            if (at >= 0.0_dp .and. rt >= 0.0_dp) then
                call ishermitian(a, flags, status, atol=at, rtol=rt)
            else if (at >= 0.0_dp) then
                call ishermitian(a, flags, status, atol=at)
            else if (rt >= 0.0_dp) then
                call ishermitian(a, flags, status, rtol=rt)
            else
                call ishermitian(a, flags, status)
            end if
        end if
    end subroutine call_structure_2d_i64

    subroutine call_structure_3d_i64(a, op, flags, status, at, rt)
        integer(int64), intent(in) :: a(:,:,:) !! input matrix
        integer, intent(in) :: op !! 2 for symmetric, 3 for Hermitian
        logical, intent(out) :: flags(:) !! results
        integer, intent(out) :: status !! status
        real(dp), intent(in) :: at !! absolute tolerance or negative for omitted
        real(dp), intent(in) :: rt !! relative tolerance or negative for omitted

        if (op == 2) then
            if (at >= 0.0_dp .and. rt >= 0.0_dp) then
                call issymmetric(a, flags, status, atol=at, rtol=rt)
            else if (at >= 0.0_dp) then
                call issymmetric(a, flags, status, atol=at)
            else if (rt >= 0.0_dp) then
                call issymmetric(a, flags, status, rtol=rt)
            else
                call issymmetric(a, flags, status)
            end if
        else if (op == 3) then
            if (at >= 0.0_dp .and. rt >= 0.0_dp) then
                call ishermitian(a, flags, status, atol=at, rtol=rt)
            else if (at >= 0.0_dp) then
                call ishermitian(a, flags, status, atol=at)
            else if (rt >= 0.0_dp) then
                call ishermitian(a, flags, status, rtol=rt)
            else
                call ishermitian(a, flags, status)
            end if
        end if
    end subroutine call_structure_3d_i64

    subroutine call_structure_2d_r(a, op, flags, status, at, rt)
        real(dp), intent(in) :: a(:,:) !! input matrix
        integer, intent(in) :: op !! 2 for symmetric, 3 for Hermitian
        logical, intent(out) :: flags !! results
        integer, intent(out) :: status !! status
        real(dp), intent(in) :: at !! absolute tolerance or negative for omitted
        real(dp), intent(in) :: rt !! relative tolerance or negative for omitted

        if (op == 2) then
            if (at >= 0.0_dp .and. rt >= 0.0_dp) then
                call issymmetric(a, flags, status, atol=at, rtol=rt)
            else if (at >= 0.0_dp) then
                call issymmetric(a, flags, status, atol=at)
            else if (rt >= 0.0_dp) then
                call issymmetric(a, flags, status, rtol=rt)
            else
                call issymmetric(a, flags, status)
            end if
        else if (op == 3) then
            if (at >= 0.0_dp .and. rt >= 0.0_dp) then
                call ishermitian(a, flags, status, atol=at, rtol=rt)
            else if (at >= 0.0_dp) then
                call ishermitian(a, flags, status, atol=at)
            else if (rt >= 0.0_dp) then
                call ishermitian(a, flags, status, rtol=rt)
            else
                call ishermitian(a, flags, status)
            end if
        end if
    end subroutine call_structure_2d_r

    subroutine call_structure_3d_r(a, op, flags, status, at, rt)
        real(dp), intent(in) :: a(:,:,:) !! input matrix
        integer, intent(in) :: op !! 2 for symmetric, 3 for Hermitian
        logical, intent(out) :: flags(:) !! results
        integer, intent(out) :: status !! status
        real(dp), intent(in) :: at !! absolute tolerance or negative for omitted
        real(dp), intent(in) :: rt !! relative tolerance or negative for omitted

        if (op == 2) then
            if (at >= 0.0_dp .and. rt >= 0.0_dp) then
                call issymmetric(a, flags, status, atol=at, rtol=rt)
            else if (at >= 0.0_dp) then
                call issymmetric(a, flags, status, atol=at)
            else if (rt >= 0.0_dp) then
                call issymmetric(a, flags, status, rtol=rt)
            else
                call issymmetric(a, flags, status)
            end if
        else if (op == 3) then
            if (at >= 0.0_dp .and. rt >= 0.0_dp) then
                call ishermitian(a, flags, status, atol=at, rtol=rt)
            else if (at >= 0.0_dp) then
                call ishermitian(a, flags, status, atol=at)
            else if (rt >= 0.0_dp) then
                call ishermitian(a, flags, status, rtol=rt)
            else
                call ishermitian(a, flags, status)
            end if
        end if
    end subroutine call_structure_3d_r

    subroutine call_structure_2d_c(a, op, flags, status, at, rt)
        complex(dp), intent(in) :: a(:,:) !! input matrix
        integer, intent(in) :: op !! 2 for symmetric, 3 for Hermitian
        logical, intent(out) :: flags !! results
        integer, intent(out) :: status !! status
        real(dp), intent(in) :: at !! absolute tolerance or negative for omitted
        real(dp), intent(in) :: rt !! relative tolerance or negative for omitted

        if (op == 2) then
            if (at >= 0.0_dp .and. rt >= 0.0_dp) then
                call issymmetric(a, flags, status, atol=at, rtol=rt)
            else if (at >= 0.0_dp) then
                call issymmetric(a, flags, status, atol=at)
            else if (rt >= 0.0_dp) then
                call issymmetric(a, flags, status, rtol=rt)
            else
                call issymmetric(a, flags, status)
            end if
        else if (op == 3) then
            if (at >= 0.0_dp .and. rt >= 0.0_dp) then
                call ishermitian(a, flags, status, atol=at, rtol=rt)
            else if (at >= 0.0_dp) then
                call ishermitian(a, flags, status, atol=at)
            else if (rt >= 0.0_dp) then
                call ishermitian(a, flags, status, rtol=rt)
            else
                call ishermitian(a, flags, status)
            end if
        end if
    end subroutine call_structure_2d_c

    subroutine call_structure_3d_c(a, op, flags, status, at, rt)
        complex(dp), intent(in) :: a(:,:,:) !! input matrix
        integer, intent(in) :: op !! 2 for symmetric, 3 for Hermitian
        logical, intent(out) :: flags(:) !! results
        integer, intent(out) :: status !! status
        real(dp), intent(in) :: at !! absolute tolerance or negative for omitted
        real(dp), intent(in) :: rt !! relative tolerance or negative for omitted

        if (op == 2) then
            if (at >= 0.0_dp .and. rt >= 0.0_dp) then
                call issymmetric(a, flags, status, atol=at, rtol=rt)
            else if (at >= 0.0_dp) then
                call issymmetric(a, flags, status, atol=at)
            else if (rt >= 0.0_dp) then
                call issymmetric(a, flags, status, rtol=rt)
            else
                call issymmetric(a, flags, status)
            end if
        else if (op == 3) then
            if (at >= 0.0_dp .and. rt >= 0.0_dp) then
                call ishermitian(a, flags, status, atol=at, rtol=rt)
            else if (at >= 0.0_dp) then
                call ishermitian(a, flags, status, atol=at)
            else if (rt >= 0.0_dp) then
                call ishermitian(a, flags, status, rtol=rt)
            else
                call ishermitian(a, flags, status)
            end if
        end if
    end subroutine call_structure_3d_c

end program linalg_structure_diagsvd_driver
