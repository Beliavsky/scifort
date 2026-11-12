! SPDX-License-Identifier: MIT
! Binary stream driver for independent scipy.ndimage spline comparisons.
program spline_nd_reference_driver
    use, intrinsic :: iso_fortran_env, only : int32
    use scifort_kinds, only : dp
    use scifort_signal_spline_nd, only : spline_filter_nd, map_coordinates_nd
    implicit none
    integer(int32) :: rank_value, degree, is_complex, mode_value, npts
    integer, allocatable :: dims(:)
    real(dp), allocatable :: rv(:), coordinates(:,:)
    complex(dp), allocatable :: zv(:)
    integer :: ios
    character(len=*), parameter :: infile='spline_nd_reference_input.bin'
    character(len=*), parameter :: outfile='spline_nd_reference_output.bin'
    open(unit=10,file=infile,status='old',access='stream',form='unformatted')
    open(unit=11,file=outfile,status='replace',access='stream',form='unformatted')
    do
        read(10,iostat=ios) rank_value, degree, is_complex, mode_value, npts
        if (ios/=0) exit
        allocate(dims(rank_value),coordinates(rank_value,npts))
        read(10) dims
        if (is_complex==0) then
            allocate(rv(product(dims)))
            read(10) rv
        else
            allocate(zv(product(dims)))
            read(10) zv
        end if
        read(10) coordinates
        if (is_complex==0) then
            call dispatch_real(rv,dims,coordinates,degree,mode_value)
            deallocate(rv)
        else
            call dispatch_complex(zv,dims,coordinates,degree,mode_value)
            deallocate(zv)
        end if
        deallocate(dims,coordinates)
    end do
    close(10)
    close(11)
contains
    subroutine dispatch_real(rv,dims,coords,degree,mode_value)
        real(dp), intent(in) :: rv(:) !! packed samples
        integer, intent(in) :: dims(:), degree, mode_value !! shape and filter setup
        real(dp), intent(in) :: coords(:,:) !! query coordinates
        character(len=8) :: mode
        integer(int32) :: status, status2, status3
        real(dp) :: yy(size(coords,2)), zz(size(coords,2))
        mode='mirror'
        if (mode_value==1) mode='constant'
        select case (size(dims))
        case (1)
            block
                real(dp) :: a(dims(1))
                real(dp) :: b(dims(1))
                a=reshape(rv,shape(a))
                call spline_filter_nd(a,b,degree,status)
                call map_coordinates_nd(b,coords,yy,degree,status2,mode=trim(mode),cval=-2.25_dp)
                call map_coordinates_nd(a,coords,zz,degree,status3,mode=trim(mode), &
                    prefilter=.true.,cval=-2.25_dp)
                write(11) status, status2, status3
                write(11) b
            end block
        case (2)
            block
                real(dp) :: a(dims(1),dims(2))
                real(dp) :: b(dims(1),dims(2))
                a=reshape(rv,shape(a))
                call spline_filter_nd(a,b,degree,status)
                call map_coordinates_nd(b,coords,yy,degree,status2,mode=trim(mode),cval=-2.25_dp)
                call map_coordinates_nd(a,coords,zz,degree,status3,mode=trim(mode), &
                    prefilter=.true.,cval=-2.25_dp)
                write(11) status, status2, status3
                write(11) b
            end block
        case (3)
            block
                real(dp) :: a(dims(1),dims(2),dims(3))
                real(dp) :: b(dims(1),dims(2),dims(3))
                a=reshape(rv,shape(a))
                call spline_filter_nd(a,b,degree,status)
                call map_coordinates_nd(b,coords,yy,degree,status2,mode=trim(mode),cval=-2.25_dp)
                call map_coordinates_nd(a,coords,zz,degree,status3,mode=trim(mode), &
                    prefilter=.true.,cval=-2.25_dp)
                write(11) status, status2, status3
                write(11) b
            end block
        case (4)
            block
                real(dp) :: a(dims(1),dims(2),dims(3),dims(4))
                real(dp) :: b(dims(1),dims(2),dims(3),dims(4))
                a=reshape(rv,shape(a))
                call spline_filter_nd(a,b,degree,status)
                call map_coordinates_nd(b,coords,yy,degree,status2,mode=trim(mode),cval=-2.25_dp)
                call map_coordinates_nd(a,coords,zz,degree,status3,mode=trim(mode), &
                    prefilter=.true.,cval=-2.25_dp)
                write(11) status, status2, status3
                write(11) b
            end block
        case (5)
            block
                real(dp) :: a(dims(1),dims(2),dims(3),dims(4),dims(5))
                real(dp) :: b(dims(1),dims(2),dims(3),dims(4),dims(5))
                a=reshape(rv,shape(a))
                call spline_filter_nd(a,b,degree,status)
                call map_coordinates_nd(b,coords,yy,degree,status2,mode=trim(mode),cval=-2.25_dp)
                call map_coordinates_nd(a,coords,zz,degree,status3,mode=trim(mode), &
                    prefilter=.true.,cval=-2.25_dp)
                write(11) status, status2, status3
                write(11) b
            end block
        case (15)
            block
                real(dp) :: a(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9),dims(10),dims(11),dims(12),dims(13),dims(14),dims(15))
                real(dp) :: b(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9),dims(10),dims(11),dims(12),dims(13),dims(14),dims(15))
                a=reshape(rv,shape(a))
                call spline_filter_nd(a,b,degree,status)
                call map_coordinates_nd(b,coords,yy,degree,status2,mode=trim(mode),cval=-2.25_dp)
                call map_coordinates_nd(a,coords,zz,degree,status3,mode=trim(mode), &
                    prefilter=.true.,cval=-2.25_dp)
                write(11) status, status2, status3
                write(11) b
            end block
        end select
        write(11) yy,zz
    end subroutine dispatch_real
    subroutine dispatch_complex(zv,dims,coords,degree,mode_value)
        complex(dp), intent(in) :: zv(:) !! packed complex samples
        integer, intent(in) :: dims(:), degree, mode_value !! shape and filter setup
        real(dp), intent(in) :: coords(:,:) !! query coordinates
        character(len=8) :: mode
        integer(int32) :: status, status2, status3
        complex(dp) :: yy(size(coords,2)), zz(size(coords,2))
        mode='mirror'
        if (mode_value==1) mode='constant'
        select case (size(dims))
        case (1)
            block
                complex(dp) :: a(dims(1))
                complex(dp) :: b(dims(1))
                a=reshape(zv,shape(a))
                call spline_filter_nd(a,b,degree,status)
                call map_coordinates_nd(b,coords,yy,degree,status2,mode=trim(mode), &
                    cval=cmplx(-2.25_dp,0.125_dp,dp))
                call map_coordinates_nd(a,coords,zz,degree,status3,mode=trim(mode), &
                    prefilter=.true.,cval=cmplx(-2.25_dp,0.125_dp,dp))
                write(11) status, status2, status3
                write(11) b
            end block
        case (2)
            block
                complex(dp) :: a(dims(1),dims(2))
                complex(dp) :: b(dims(1),dims(2))
                a=reshape(zv,shape(a))
                call spline_filter_nd(a,b,degree,status)
                call map_coordinates_nd(b,coords,yy,degree,status2,mode=trim(mode), &
                    cval=cmplx(-2.25_dp,0.125_dp,dp))
                call map_coordinates_nd(a,coords,zz,degree,status3,mode=trim(mode), &
                    prefilter=.true.,cval=cmplx(-2.25_dp,0.125_dp,dp))
                write(11) status, status2, status3
                write(11) b
            end block
        case (3)
            block
                complex(dp) :: a(dims(1),dims(2),dims(3))
                complex(dp) :: b(dims(1),dims(2),dims(3))
                a=reshape(zv,shape(a))
                call spline_filter_nd(a,b,degree,status)
                call map_coordinates_nd(b,coords,yy,degree,status2,mode=trim(mode), &
                    cval=cmplx(-2.25_dp,0.125_dp,dp))
                call map_coordinates_nd(a,coords,zz,degree,status3,mode=trim(mode), &
                    prefilter=.true.,cval=cmplx(-2.25_dp,0.125_dp,dp))
                write(11) status, status2, status3
                write(11) b
            end block
        case (4)
            block
                complex(dp) :: a(dims(1),dims(2),dims(3),dims(4))
                complex(dp) :: b(dims(1),dims(2),dims(3),dims(4))
                a=reshape(zv,shape(a))
                call spline_filter_nd(a,b,degree,status)
                call map_coordinates_nd(b,coords,yy,degree,status2,mode=trim(mode), &
                    cval=cmplx(-2.25_dp,0.125_dp,dp))
                call map_coordinates_nd(a,coords,zz,degree,status3,mode=trim(mode), &
                    prefilter=.true.,cval=cmplx(-2.25_dp,0.125_dp,dp))
                write(11) status, status2, status3
                write(11) b
            end block
        case (5)
            block
                complex(dp) :: a(dims(1),dims(2),dims(3),dims(4),dims(5))
                complex(dp) :: b(dims(1),dims(2),dims(3),dims(4),dims(5))
                a=reshape(zv,shape(a))
                call spline_filter_nd(a,b,degree,status)
                call map_coordinates_nd(b,coords,yy,degree,status2,mode=trim(mode), &
                    cval=cmplx(-2.25_dp,0.125_dp,dp))
                call map_coordinates_nd(a,coords,zz,degree,status3,mode=trim(mode), &
                    prefilter=.true.,cval=cmplx(-2.25_dp,0.125_dp,dp))
                write(11) status, status2, status3
                write(11) b
            end block
        case (15)
            block
                complex(dp) :: a(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9),dims(10),dims(11),dims(12),dims(13),dims(14),dims(15))
                complex(dp) :: b(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9),dims(10),dims(11),dims(12),dims(13),dims(14),dims(15))
                a=reshape(zv,shape(a))
                call spline_filter_nd(a,b,degree,status)
                call map_coordinates_nd(b,coords,yy,degree,status2,mode=trim(mode), &
                    cval=cmplx(-2.25_dp,0.125_dp,dp))
                call map_coordinates_nd(a,coords,zz,degree,status3,mode=trim(mode), &
                    prefilter=.true.,cval=cmplx(-2.25_dp,0.125_dp,dp))
                write(11) status, status2, status3
                write(11) b
            end block
        end select
        write(11) yy,zz
    end subroutine dispatch_complex
end program spline_nd_reference_driver
