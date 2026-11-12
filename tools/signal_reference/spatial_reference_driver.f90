! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Batch reference driver: see compare_signal_spatial_scipy.py for input format.
program spatial_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_spatial, only : order_filter, sepfir2d
    implicit none
    integer :: ncases, caseid, method, xtype, dtype, ndim, rank, code
    integer :: dims(3), fdims(3), n, m
    real(dp), allocatable :: xr(:), xi(:), dr(:), rr(:), ri(:), cr(:), ci(:)
    real(dp), allocatable :: o_real_1(:)
    real(dp), allocatable :: o_real_2(:,:)
    real(dp), allocatable :: o_real_3(:,:,:)
    integer, allocatable :: o_int_1(:)
    integer, allocatable :: o_int_2(:,:)
    integer, allocatable :: o_int_3(:,:,:)
    real(dp), allocatable :: sr(:,:)
    complex(dp), allocatable :: sc(:,:)

    read(*,*) ncases
    do caseid = 1, ncases
        read(*,*) method, xtype, dtype, ndim, rank
        read(*,*) dims
        read(*,*) fdims
        n = product(dims(:ndim))
        allocate(xr(n), xi(n))
        read(*,*) xr
        read(*,*) xi
        if (method == 0) then
            m = product(fdims(:ndim))
            allocate(dr(m))
            read(*,*) dr
            select case (xtype)
            case (0)
                select case (ndim)
                case (1)
                    select case (dtype)
                    case (0)
                        o_real_1 = order_filter(reshape(xr,dims(:1)), &
                            nint(reshape(dr, fdims(:1))), rank, code)
                    case (1)
                        o_real_1 = order_filter(reshape(xr,dims(:1)), &
                            reshape(dr, fdims(:1)), rank, code)
                    case (2)
                        o_real_1 = order_filter(reshape(xr,dims(:1)), &
                            reshape(dr /= 0.0_dp, fdims(:1)), rank, code)
                    end select
                    call emit(code, reshape(o_real_1, [size(o_real_1)]))
                case (2)
                    select case (dtype)
                    case (0)
                        o_real_2 = order_filter(reshape(xr,dims(:2)), &
                            nint(reshape(dr, fdims(:2))), rank, code)
                    case (1)
                        o_real_2 = order_filter(reshape(xr,dims(:2)), &
                            reshape(dr, fdims(:2)), rank, code)
                    case (2)
                        o_real_2 = order_filter(reshape(xr,dims(:2)), &
                            reshape(dr /= 0.0_dp, fdims(:2)), rank, code)
                    end select
                    call emit(code, reshape(o_real_2, [size(o_real_2)]))
                case (3)
                    select case (dtype)
                    case (0)
                        o_real_3 = order_filter(reshape(xr,dims(:3)), &
                            nint(reshape(dr, fdims(:3))), rank, code)
                    case (1)
                        o_real_3 = order_filter(reshape(xr,dims(:3)), &
                            reshape(dr, fdims(:3)), rank, code)
                    case (2)
                        o_real_3 = order_filter(reshape(xr,dims(:3)), &
                            reshape(dr /= 0.0_dp, fdims(:3)), rank, code)
                    end select
                    call emit(code, reshape(o_real_3, [size(o_real_3)]))
                end select
            case (1)
                select case (ndim)
                case (1)
                    select case (dtype)
                    case (0)
                        o_int_1 = order_filter(reshape(nint(xr),dims(:1)), &
                            nint(reshape(dr, fdims(:1))), rank, code)
                    case (1)
                        o_int_1 = order_filter(reshape(nint(xr),dims(:1)), &
                            reshape(dr, fdims(:1)), rank, code)
                    case (2)
                        o_int_1 = order_filter(reshape(nint(xr),dims(:1)), &
                            reshape(dr /= 0.0_dp, fdims(:1)), rank, code)
                    end select
                    call emit(code, real(reshape(o_int_1, [size(o_int_1)]), dp))
                case (2)
                    select case (dtype)
                    case (0)
                        o_int_2 = order_filter(reshape(nint(xr),dims(:2)), &
                            nint(reshape(dr, fdims(:2))), rank, code)
                    case (1)
                        o_int_2 = order_filter(reshape(nint(xr),dims(:2)), &
                            reshape(dr, fdims(:2)), rank, code)
                    case (2)
                        o_int_2 = order_filter(reshape(nint(xr),dims(:2)), &
                            reshape(dr /= 0.0_dp, fdims(:2)), rank, code)
                    end select
                    call emit(code, real(reshape(o_int_2, [size(o_int_2)]), dp))
                case (3)
                    select case (dtype)
                    case (0)
                        o_int_3 = order_filter(reshape(nint(xr),dims(:3)), &
                            nint(reshape(dr, fdims(:3))), rank, code)
                    case (1)
                        o_int_3 = order_filter(reshape(nint(xr),dims(:3)), &
                            reshape(dr, fdims(:3)), rank, code)
                    case (2)
                        o_int_3 = order_filter(reshape(nint(xr),dims(:3)), &
                            reshape(dr /= 0.0_dp, fdims(:3)), rank, code)
                    end select
                    call emit(code, real(reshape(o_int_3, [size(o_int_3)]), dp))
                end select
            end select
            deallocate(dr)
        else
            allocate(rr(fdims(1)), ri(fdims(1)), cr(fdims(2)), ci(fdims(2)))
            read(*,*) rr
            read(*,*) ri
            read(*,*) cr
            read(*,*) ci
            ! xtype = 0 (real), 1 (complex), or 2 (integer)
            ! dtype bits: 0/1 for hrow real/complex, 0/2 for hcol real/complex.
            select case (xtype)
            case (0)
                select case (dtype)
                case (0)
                    sr = sepfir2d(reshape(xr,dims(:2)), &
                        rr, cr, code)
                    call emit(code, reshape(sr, [size(sr)]))
                case (2)
                    sc = sepfir2d(reshape(xr,dims(:2)), &
                        rr, cmplx(cr,ci,dp), code)
                    call emit(code, real(reshape(sc,[size(sc)]),dp), &
                        aimag(reshape(sc,[size(sc)])))
                case (1)
                    sc = sepfir2d(reshape(xr,dims(:2)), &
                        cmplx(rr,ri,dp), cr, code)
                    call emit(code, real(reshape(sc,[size(sc)]),dp), &
                        aimag(reshape(sc,[size(sc)])))
                case (3)
                    sc = sepfir2d(reshape(xr,dims(:2)), &
                        cmplx(rr,ri,dp), cmplx(cr,ci,dp), code)
                    call emit(code, real(reshape(sc,[size(sc)]),dp), &
                        aimag(reshape(sc,[size(sc)])))
                end select
            case (1)
                select case (dtype)
                case (0)
                    sc = sepfir2d(reshape(cmplx(xr,xi,dp),dims(:2)), &
                        rr, cr, code)
                    call emit(code, real(reshape(sc,[size(sc)]),dp), &
                        aimag(reshape(sc,[size(sc)])))
                case (2)
                    sc = sepfir2d(reshape(cmplx(xr,xi,dp),dims(:2)), &
                        rr, cmplx(cr,ci,dp), code)
                    call emit(code, real(reshape(sc,[size(sc)]),dp), &
                        aimag(reshape(sc,[size(sc)])))
                case (1)
                    sc = sepfir2d(reshape(cmplx(xr,xi,dp),dims(:2)), &
                        cmplx(rr,ri,dp), cr, code)
                    call emit(code, real(reshape(sc,[size(sc)]),dp), &
                        aimag(reshape(sc,[size(sc)])))
                case (3)
                    sc = sepfir2d(reshape(cmplx(xr,xi,dp),dims(:2)), &
                        cmplx(rr,ri,dp), cmplx(cr,ci,dp), code)
                    call emit(code, real(reshape(sc,[size(sc)]),dp), &
                        aimag(reshape(sc,[size(sc)])))
                end select
            case (2)
                select case (dtype)
                case (0)
                    sr = sepfir2d(reshape(nint(xr),dims(:2)), &
                        rr, cr, code)
                    call emit(code, reshape(sr, [size(sr)]))
                case (2)
                    sc = sepfir2d(reshape(nint(xr),dims(:2)), &
                        rr, cmplx(cr,ci,dp), code)
                    call emit(code, real(reshape(sc,[size(sc)]),dp), &
                        aimag(reshape(sc,[size(sc)])))
                case (1)
                    sc = sepfir2d(reshape(nint(xr),dims(:2)), &
                        cmplx(rr,ri,dp), cr, code)
                    call emit(code, real(reshape(sc,[size(sc)]),dp), &
                        aimag(reshape(sc,[size(sc)])))
                case (3)
                    sc = sepfir2d(reshape(nint(xr),dims(:2)), &
                        cmplx(rr,ri,dp), cmplx(cr,ci,dp), code)
                    call emit(code, real(reshape(sc,[size(sc)]),dp), &
                        aimag(reshape(sc,[size(sc)])))
                end select
            end select
            deallocate(rr,ri,cr,ci)
        end if
        deallocate(xr,xi)
    end do

contains

    subroutine emit(status, re, im)
        integer, intent(in) :: status !! status from SciFort routine
        real(dp), intent(in) :: re(:) !! flattened real results
        real(dp), intent(in), optional :: im(:) !! optional imaginary parts
        integer :: i

        write(*,'(i0)', advance='no') status
        do i = 1, size(re)
            if (present(im)) then
                write(*,'(1x,es25.16e3,1x,es25.16e3)', advance='no') re(i), im(i)
            else
                write(*,'(1x,es25.16e3,1x,es25.16e3)', advance='no') re(i), 0.0_dp
            end if
        end do
        write(*,*)
    end subroutine emit
end program spatial_reference_driver
