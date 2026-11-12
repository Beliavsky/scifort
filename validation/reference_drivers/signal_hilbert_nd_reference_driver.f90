! SPDX-License-Identifier: MIT
! Pipe/file-based SciPy independent comparison driver.
program signal_hilbert_nd_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_hilbert_nd, only : hilbert_into, hilbert2_into
    implicit none
    integer :: input, output, cases, icase, ndim, mode, a1, a2, d, code, j
    integer :: dims(15), outdims(15), ax(2), lengths(2)
    character(1024) :: inpath, outpath
    real(dp), allocatable :: v(:)
    complex(dp), allocatable :: result(:)

    real(dp), allocatable :: x1(:)
    complex(dp), allocatable :: y1(:)
    real(dp), allocatable :: x2(:,:)
    complex(dp), allocatable :: y2(:,:)
    real(dp), allocatable :: x3(:,:,:)
    complex(dp), allocatable :: y3(:,:,:)
    real(dp), allocatable :: x5(:,:,:,:,:)
    complex(dp), allocatable :: y5(:,:,:,:,:)
    real(dp), allocatable :: x15(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
    complex(dp), allocatable :: y15(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)

    call get_command_argument(1,inpath)
    call get_command_argument(2,outpath)
    open(newunit=input,file=trim(inpath),status='old',action='read')
    open(newunit=output,file=trim(outpath),status='replace',action='write')
    read(input,*) cases
    do icase=1,cases
        read(input,*) ndim, mode, a1, a2, dims, outdims
        if (ndim < 1 .or. ndim > 15) error stop 'invalid rank in driver'
        allocate(v(product(dims(1:ndim))))
        read(input,*) v
        ax = [a1,a2]
        if (ax(1) < 0) ax(1) = ax(1)+ndim+1
        if (mode == 2 .and. ax(2) < 0) ax(2) = ax(2)+ndim+1
        lengths = 1
        lengths(1) = outdims(ax(1))
        if (mode == 2) lengths(2) = outdims(ax(2))
        select case(ndim)

        case (1)
            allocate(x1(dims(1)))
            allocate(y1(outdims(1)))
            x1 = reshape(v,shape(x1))
            y1 = cmplx(-999.0_dp,999.0_dp,dp)
            if (mode == 1) then
                call hilbert_into(x1,y1,axis=a1,n=lengths(1),status=code)
            else
                call hilbert2_into(x1,y1,axes=[a1,a2],n=lengths,status=code)
            end if
            result = reshape(y1,[size(y1)])
            deallocate(x1,y1)
        case (2)
            allocate(x2(dims(1),dims(2)))
            allocate(y2(outdims(1),outdims(2)))
            x2 = reshape(v,shape(x2))
            y2 = cmplx(-999.0_dp,999.0_dp,dp)
            if (mode == 1) then
                call hilbert_into(x2,y2,axis=a1,n=lengths(1),status=code)
            else
                call hilbert2_into(x2,y2,axes=[a1,a2],n=lengths,status=code)
            end if
            result = reshape(y2,[size(y2)])
            deallocate(x2,y2)
        case (3)
            allocate(x3(dims(1),dims(2),dims(3)))
            allocate(y3(outdims(1),outdims(2),outdims(3)))
            x3 = reshape(v,shape(x3))
            y3 = cmplx(-999.0_dp,999.0_dp,dp)
            if (mode == 1) then
                call hilbert_into(x3,y3,axis=a1,n=lengths(1),status=code)
            else
                call hilbert2_into(x3,y3,axes=[a1,a2],n=lengths,status=code)
            end if
            result = reshape(y3,[size(y3)])
            deallocate(x3,y3)
        case (5)
            allocate(x5(dims(1),dims(2),dims(3),dims(4),dims(5)))
            allocate(y5(outdims(1),outdims(2),outdims(3),outdims(4),outdims(5)))
            x5 = reshape(v,shape(x5))
            y5 = cmplx(-999.0_dp,999.0_dp,dp)
            if (mode == 1) then
                call hilbert_into(x5,y5,axis=a1,n=lengths(1),status=code)
            else
                call hilbert2_into(x5,y5,axes=[a1,a2],n=lengths,status=code)
            end if
            result = reshape(y5,[size(y5)])
            deallocate(x5,y5)
        case (15)
            allocate(x15(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8), &
                dims(9),dims(10),dims(11),dims(12),dims(13),dims(14),dims(15)))
            allocate(y15(outdims(1),outdims(2),outdims(3),outdims(4),outdims(5),outdims(6),outdims(7),outdims(8), &
                outdims(9),outdims(10),outdims(11),outdims(12),outdims(13),outdims(14),outdims(15)))
            x15 = reshape(v,shape(x15))
            y15 = cmplx(-999.0_dp,999.0_dp,dp)
            if (mode == 1) then
                call hilbert_into(x15,y15,axis=a1,n=lengths(1),status=code)
            else
                call hilbert2_into(x15,y15,axes=[a1,a2],n=lengths,status=code)
            end if
            result = reshape(y15,[size(y15)])
            deallocate(x15,y15)
        case default
            error stop 'unsupported rank'
        end select
        write(output,'(i0)') code
        write(output,'(*(es25.16,1x))') (real(result(j),dp),aimag(result(j)),j=1,size(result))
        deallocate(v,result)
    end do
    close(input)
    close(output)
end program signal_hilbert_nd_reference_driver
