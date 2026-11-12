! SPDX-License-Identifier: MIT
! Reference driver for assumed-rank multirate signal APIs.
program multirate_assumed_rank_driver
    use scifort_signal, only : dp, resample_into, resample_poly_into, upfirdn_into, decimate_into
    implicit none
    character(32) :: op, dtype, htype, mode
    integer :: nrank, ax, p1, p2, st, i, j
    integer :: dims(15), yd(15)
    real(dp), allocatable :: re(:), im(:)
    complex(dp), allocatable :: output(:)
    real(dp) :: h(3)
    complex(dp) :: hc(3)
    read(*,*) op,dtype,htype,mode,ax,p1,p2,nrank
    read(*,*) (dims(i),i=1,nrank)
    read(*,*) (yd(i),i=1,nrank)
    allocate(re(product(dims(:nrank))),im(product(dims(:nrank))))
    read(*,*) (re(i),i=1,size(re))
    read(*,*) (im(i),i=1,size(im))
    h=[0.2_dp,0.6_dp,0.2_dp]
    hc=cmplx(h,[0.07_dp,-0.02_dp,0.11_dp],dp)
    st=-1
    select case(nrank)
    case (4)
        block
            real(dp), allocatable :: xr(:,:,:,:), yr(:,:,:,:)
            complex(dp), allocatable :: xc(:,:,:,:), yc(:,:,:,:)
            allocate(xr(dims(1),dims(2),dims(3),dims(4)))
            allocate(yr(yd(1),yd(2),yd(3),yd(4)))
            allocate(xc(dims(1),dims(2),dims(3),dims(4)))
            allocate(yc(yd(1),yd(2),yd(3),yd(4)))
            xr=reshape(re,shape(xr))
            xc=reshape(cmplx(re,im,dp),shape(xc))
            yr=0.0_dp
            yc=cmplx(0.0_dp,0.0_dp,dp)
            select case(trim(op))
            case ('resample')
                if (trim(dtype)=='r') then
                    if (trim(mode)=='hann') then
                        call resample_into(xr,yr,p1,window='hann',axis=ax,status=st)
                    else
                        call resample_into(xr,yr,p1,axis=ax,status=st)
                    end if
                else
                    if (trim(mode)=='hann') then
                        call resample_into(xc,yc,p1,window='hann',axis=ax,status=st)
                    else
                        call resample_into(xc,yc,p1,domain=trim(mode),axis=ax,status=st)
                    end if
                end if
            case ('poly')
                if (trim(dtype)=='r') then
                    if (trim(htype)=='r') then
                        call resample_poly_into(xr,yr,p1,p2,taps=h,padtype=trim(mode),axis=ax,status=st)
                    else
                        call resample_poly_into(xr,yc,p1,p2,taps=hc,padtype=trim(mode),axis=ax,status=st)
                    end if
                else
                    if (trim(htype)=='r') then
                        call resample_poly_into(xc,yc,p1,p2,taps=h,padtype=trim(mode),axis=ax,status=st)
                    else
                        call resample_poly_into(xc,yc,p1,p2,taps=hc,padtype=trim(mode),axis=ax,status=st)
                    end if
                end if
            case ('upfirdn')
                if (trim(dtype)=='r') then
                    if (trim(htype)=='r') then
                        call upfirdn_into(h,xr,yr,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    else
                        call upfirdn_into(hc,xr,yc,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    end if
                else
                    if (trim(htype)=='r') then
                        call upfirdn_into(h,xc,yc,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    else
                        call upfirdn_into(hc,xc,yc,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    end if
                end if
            case ('decimate')
                if (trim(dtype)=='r') then
                    call decimate_into(xr,yr,p1,n=8,ftype=trim(mode),zero_phase=.false.,axis=ax,status=st)
                else
                    call decimate_into(xc,yc,p1,n=8,ftype=trim(mode),zero_phase=.false.,axis=ax,status=st)
                end if
            end select
            if (trim(dtype)=='r' .and. &
                (trim(htype)=='r' .or. trim(op)=='resample' .or. trim(op)=='decimate')) then
                output=cmplx(reshape(yr,[size(yr)]),0.0_dp,dp)
            else
                output=reshape(yc,[size(yc)])
            end if
        end block
    case (5)
        block
            real(dp), allocatable :: xr(:,:,:,:,:), yr(:,:,:,:,:)
            complex(dp), allocatable :: xc(:,:,:,:,:), yc(:,:,:,:,:)
            allocate(xr(dims(1),dims(2),dims(3),dims(4),dims(5)))
            allocate(yr(yd(1),yd(2),yd(3),yd(4),yd(5)))
            allocate(xc(dims(1),dims(2),dims(3),dims(4),dims(5)))
            allocate(yc(yd(1),yd(2),yd(3),yd(4),yd(5)))
            xr=reshape(re,shape(xr))
            xc=reshape(cmplx(re,im,dp),shape(xc))
            yr=0.0_dp
            yc=cmplx(0.0_dp,0.0_dp,dp)
            select case(trim(op))
            case ('resample')
                if (trim(dtype)=='r') then
                    if (trim(mode)=='hann') then
                        call resample_into(xr,yr,p1,window='hann',axis=ax,status=st)
                    else
                        call resample_into(xr,yr,p1,axis=ax,status=st)
                    end if
                else
                    if (trim(mode)=='hann') then
                        call resample_into(xc,yc,p1,window='hann',axis=ax,status=st)
                    else
                        call resample_into(xc,yc,p1,domain=trim(mode),axis=ax,status=st)
                    end if
                end if
            case ('poly')
                if (trim(dtype)=='r') then
                    if (trim(htype)=='r') then
                        call resample_poly_into(xr,yr,p1,p2,taps=h,padtype=trim(mode),axis=ax,status=st)
                    else
                        call resample_poly_into(xr,yc,p1,p2,taps=hc,padtype=trim(mode),axis=ax,status=st)
                    end if
                else
                    if (trim(htype)=='r') then
                        call resample_poly_into(xc,yc,p1,p2,taps=h,padtype=trim(mode),axis=ax,status=st)
                    else
                        call resample_poly_into(xc,yc,p1,p2,taps=hc,padtype=trim(mode),axis=ax,status=st)
                    end if
                end if
            case ('upfirdn')
                if (trim(dtype)=='r') then
                    if (trim(htype)=='r') then
                        call upfirdn_into(h,xr,yr,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    else
                        call upfirdn_into(hc,xr,yc,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    end if
                else
                    if (trim(htype)=='r') then
                        call upfirdn_into(h,xc,yc,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    else
                        call upfirdn_into(hc,xc,yc,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    end if
                end if
            case ('decimate')
                if (trim(dtype)=='r') then
                    call decimate_into(xr,yr,p1,n=8,ftype=trim(mode),zero_phase=.false.,axis=ax,status=st)
                else
                    call decimate_into(xc,yc,p1,n=8,ftype=trim(mode),zero_phase=.false.,axis=ax,status=st)
                end if
            end select
            if (trim(dtype)=='r' .and. &
                (trim(htype)=='r' .or. trim(op)=='resample' .or. trim(op)=='decimate')) then
                output=cmplx(reshape(yr,[size(yr)]),0.0_dp,dp)
            else
                output=reshape(yc,[size(yc)])
            end if
        end block
    case (6)
        block
            real(dp), allocatable :: xr(:,:,:,:,:,:), yr(:,:,:,:,:,:)
            complex(dp), allocatable :: xc(:,:,:,:,:,:), yc(:,:,:,:,:,:)
            allocate(xr(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)))
            allocate(yr(yd(1),yd(2),yd(3),yd(4),yd(5),yd(6)))
            allocate(xc(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)))
            allocate(yc(yd(1),yd(2),yd(3),yd(4),yd(5),yd(6)))
            xr=reshape(re,shape(xr))
            xc=reshape(cmplx(re,im,dp),shape(xc))
            yr=0.0_dp
            yc=cmplx(0.0_dp,0.0_dp,dp)
            select case(trim(op))
            case ('resample')
                if (trim(dtype)=='r') then
                    if (trim(mode)=='hann') then
                        call resample_into(xr,yr,p1,window='hann',axis=ax,status=st)
                    else
                        call resample_into(xr,yr,p1,axis=ax,status=st)
                    end if
                else
                    if (trim(mode)=='hann') then
                        call resample_into(xc,yc,p1,window='hann',axis=ax,status=st)
                    else
                        call resample_into(xc,yc,p1,domain=trim(mode),axis=ax,status=st)
                    end if
                end if
            case ('poly')
                if (trim(dtype)=='r') then
                    if (trim(htype)=='r') then
                        call resample_poly_into(xr,yr,p1,p2,taps=h,padtype=trim(mode),axis=ax,status=st)
                    else
                        call resample_poly_into(xr,yc,p1,p2,taps=hc,padtype=trim(mode),axis=ax,status=st)
                    end if
                else
                    if (trim(htype)=='r') then
                        call resample_poly_into(xc,yc,p1,p2,taps=h,padtype=trim(mode),axis=ax,status=st)
                    else
                        call resample_poly_into(xc,yc,p1,p2,taps=hc,padtype=trim(mode),axis=ax,status=st)
                    end if
                end if
            case ('upfirdn')
                if (trim(dtype)=='r') then
                    if (trim(htype)=='r') then
                        call upfirdn_into(h,xr,yr,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    else
                        call upfirdn_into(hc,xr,yc,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    end if
                else
                    if (trim(htype)=='r') then
                        call upfirdn_into(h,xc,yc,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    else
                        call upfirdn_into(hc,xc,yc,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    end if
                end if
            case ('decimate')
                if (trim(dtype)=='r') then
                    call decimate_into(xr,yr,p1,n=8,ftype=trim(mode),zero_phase=.false.,axis=ax,status=st)
                else
                    call decimate_into(xc,yc,p1,n=8,ftype=trim(mode),zero_phase=.false.,axis=ax,status=st)
                end if
            end select
            if (trim(dtype)=='r' .and. &
                (trim(htype)=='r' .or. trim(op)=='resample' .or. trim(op)=='decimate')) then
                output=cmplx(reshape(yr,[size(yr)]),0.0_dp,dp)
            else
                output=reshape(yc,[size(yc)])
            end if
        end block
    case (8)
        block
            real(dp), allocatable :: xr(:,:,:,:,:,:,:,:), yr(:,:,:,:,:,:,:,:)
            complex(dp), allocatable :: xc(:,:,:,:,:,:,:,:), yc(:,:,:,:,:,:,:,:)
            allocate(xr(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8)))
            allocate(yr(yd(1),yd(2),yd(3),yd(4),yd(5),yd(6),yd(7),yd(8)))
            allocate(xc(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8)))
            allocate(yc(yd(1),yd(2),yd(3),yd(4),yd(5),yd(6),yd(7),yd(8)))
            xr=reshape(re,shape(xr))
            xc=reshape(cmplx(re,im,dp),shape(xc))
            yr=0.0_dp
            yc=cmplx(0.0_dp,0.0_dp,dp)
            select case(trim(op))
            case ('resample')
                if (trim(dtype)=='r') then
                    if (trim(mode)=='hann') then
                        call resample_into(xr,yr,p1,window='hann',axis=ax,status=st)
                    else
                        call resample_into(xr,yr,p1,axis=ax,status=st)
                    end if
                else
                    if (trim(mode)=='hann') then
                        call resample_into(xc,yc,p1,window='hann',axis=ax,status=st)
                    else
                        call resample_into(xc,yc,p1,domain=trim(mode),axis=ax,status=st)
                    end if
                end if
            case ('poly')
                if (trim(dtype)=='r') then
                    if (trim(htype)=='r') then
                        call resample_poly_into(xr,yr,p1,p2,taps=h,padtype=trim(mode),axis=ax,status=st)
                    else
                        call resample_poly_into(xr,yc,p1,p2,taps=hc,padtype=trim(mode),axis=ax,status=st)
                    end if
                else
                    if (trim(htype)=='r') then
                        call resample_poly_into(xc,yc,p1,p2,taps=h,padtype=trim(mode),axis=ax,status=st)
                    else
                        call resample_poly_into(xc,yc,p1,p2,taps=hc,padtype=trim(mode),axis=ax,status=st)
                    end if
                end if
            case ('upfirdn')
                if (trim(dtype)=='r') then
                    if (trim(htype)=='r') then
                        call upfirdn_into(h,xr,yr,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    else
                        call upfirdn_into(hc,xr,yc,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    end if
                else
                    if (trim(htype)=='r') then
                        call upfirdn_into(h,xc,yc,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    else
                        call upfirdn_into(hc,xc,yc,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    end if
                end if
            case ('decimate')
                if (trim(dtype)=='r') then
                    call decimate_into(xr,yr,p1,n=8,ftype=trim(mode),zero_phase=.false.,axis=ax,status=st)
                else
                    call decimate_into(xc,yc,p1,n=8,ftype=trim(mode),zero_phase=.false.,axis=ax,status=st)
                end if
            end select
            if (trim(dtype)=='r' .and. &
                (trim(htype)=='r' .or. trim(op)=='resample' .or. trim(op)=='decimate')) then
                output=cmplx(reshape(yr,[size(yr)]),0.0_dp,dp)
            else
                output=reshape(yc,[size(yc)])
            end if
        end block
    case (14)
        block
            real(dp), allocatable :: xr(:,:,:,:,:,:,:,:,:,:,:,:,:,:), yr(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
            complex(dp), allocatable :: xc(:,:,:,:,:,:,:,:,:,:,:,:,:,:), yc(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
            allocate(xr(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9), &
                dims(10),dims(11),dims(12),dims(13),dims(14)))
            allocate(yr(yd(1),yd(2),yd(3),yd(4),yd(5),yd(6),yd(7),yd(8),yd(9),yd(10),yd(11),yd(12),yd(13),yd(14)))
            allocate(xc(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9), &
                dims(10),dims(11),dims(12),dims(13),dims(14)))
            allocate(yc(yd(1),yd(2),yd(3),yd(4),yd(5),yd(6),yd(7),yd(8),yd(9),yd(10),yd(11),yd(12),yd(13),yd(14)))
            xr=reshape(re,shape(xr))
            xc=reshape(cmplx(re,im,dp),shape(xc))
            yr=0.0_dp
            yc=cmplx(0.0_dp,0.0_dp,dp)
            select case(trim(op))
            case ('resample')
                if (trim(dtype)=='r') then
                    if (trim(mode)=='hann') then
                        call resample_into(xr,yr,p1,window='hann',axis=ax,status=st)
                    else
                        call resample_into(xr,yr,p1,axis=ax,status=st)
                    end if
                else
                    if (trim(mode)=='hann') then
                        call resample_into(xc,yc,p1,window='hann',axis=ax,status=st)
                    else
                        call resample_into(xc,yc,p1,domain=trim(mode),axis=ax,status=st)
                    end if
                end if
            case ('poly')
                if (trim(dtype)=='r') then
                    if (trim(htype)=='r') then
                        call resample_poly_into(xr,yr,p1,p2,taps=h,padtype=trim(mode),axis=ax,status=st)
                    else
                        call resample_poly_into(xr,yc,p1,p2,taps=hc,padtype=trim(mode),axis=ax,status=st)
                    end if
                else
                    if (trim(htype)=='r') then
                        call resample_poly_into(xc,yc,p1,p2,taps=h,padtype=trim(mode),axis=ax,status=st)
                    else
                        call resample_poly_into(xc,yc,p1,p2,taps=hc,padtype=trim(mode),axis=ax,status=st)
                    end if
                end if
            case ('upfirdn')
                if (trim(dtype)=='r') then
                    if (trim(htype)=='r') then
                        call upfirdn_into(h,xr,yr,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    else
                        call upfirdn_into(hc,xr,yc,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    end if
                else
                    if (trim(htype)=='r') then
                        call upfirdn_into(h,xc,yc,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    else
                        call upfirdn_into(hc,xc,yc,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    end if
                end if
            case ('decimate')
                if (trim(dtype)=='r') then
                    call decimate_into(xr,yr,p1,n=8,ftype=trim(mode),zero_phase=.false.,axis=ax,status=st)
                else
                    call decimate_into(xc,yc,p1,n=8,ftype=trim(mode),zero_phase=.false.,axis=ax,status=st)
                end if
            end select
            if (trim(dtype)=='r' .and. &
                (trim(htype)=='r' .or. trim(op)=='resample' .or. trim(op)=='decimate')) then
                output=cmplx(reshape(yr,[size(yr)]),0.0_dp,dp)
            else
                output=reshape(yc,[size(yc)])
            end if
        end block
    case (15)
        block
            real(dp), allocatable :: xr(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:), yr(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
            complex(dp), allocatable :: xc(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:), yc(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
            allocate(xr(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9), &
                dims(10),dims(11),dims(12),dims(13),dims(14),dims(15)))
            allocate(yr(yd(1),yd(2),yd(3),yd(4),yd(5),yd(6),yd(7),yd(8),yd(9),yd(10),yd(11),yd(12),yd(13),yd(14),yd(15)))
            allocate(xc(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7),dims(8),dims(9), &
                dims(10),dims(11),dims(12),dims(13),dims(14),dims(15)))
            allocate(yc(yd(1),yd(2),yd(3),yd(4),yd(5),yd(6),yd(7),yd(8),yd(9),yd(10),yd(11),yd(12),yd(13),yd(14),yd(15)))
            xr=reshape(re,shape(xr))
            xc=reshape(cmplx(re,im,dp),shape(xc))
            yr=0.0_dp
            yc=cmplx(0.0_dp,0.0_dp,dp)
            select case(trim(op))
            case ('resample')
                if (trim(dtype)=='r') then
                    if (trim(mode)=='hann') then
                        call resample_into(xr,yr,p1,window='hann',axis=ax,status=st)
                    else
                        call resample_into(xr,yr,p1,axis=ax,status=st)
                    end if
                else
                    if (trim(mode)=='hann') then
                        call resample_into(xc,yc,p1,window='hann',axis=ax,status=st)
                    else
                        call resample_into(xc,yc,p1,domain=trim(mode),axis=ax,status=st)
                    end if
                end if
            case ('poly')
                if (trim(dtype)=='r') then
                    if (trim(htype)=='r') then
                        call resample_poly_into(xr,yr,p1,p2,taps=h,padtype=trim(mode),axis=ax,status=st)
                    else
                        call resample_poly_into(xr,yc,p1,p2,taps=hc,padtype=trim(mode),axis=ax,status=st)
                    end if
                else
                    if (trim(htype)=='r') then
                        call resample_poly_into(xc,yc,p1,p2,taps=h,padtype=trim(mode),axis=ax,status=st)
                    else
                        call resample_poly_into(xc,yc,p1,p2,taps=hc,padtype=trim(mode),axis=ax,status=st)
                    end if
                end if
            case ('upfirdn')
                if (trim(dtype)=='r') then
                    if (trim(htype)=='r') then
                        call upfirdn_into(h,xr,yr,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    else
                        call upfirdn_into(hc,xr,yc,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    end if
                else
                    if (trim(htype)=='r') then
                        call upfirdn_into(h,xc,yc,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    else
                        call upfirdn_into(hc,xc,yc,up=p1,down=p2,mode=trim(mode),axis=ax,status=st)
                    end if
                end if
            case ('decimate')
                if (trim(dtype)=='r') then
                    call decimate_into(xr,yr,p1,n=8,ftype=trim(mode),zero_phase=.false.,axis=ax,status=st)
                else
                    call decimate_into(xc,yc,p1,n=8,ftype=trim(mode),zero_phase=.false.,axis=ax,status=st)
                end if
            end select
            if (trim(dtype)=='r' .and. &
                (trim(htype)=='r' .or. trim(op)=='resample' .or. trim(op)=='decimate')) then
                output=cmplx(reshape(yr,[size(yr)]),0.0_dp,dp)
            else
                output=reshape(yc,[size(yc)])
            end if
        end block
    case default
        error stop 'invalid rank'
    end select
    write(*,'(a,1x,i0)') 'STATUS',st
    do j=1,size(output)
        write(*,'(2(es26.17,1x))') real(output(j),dp),aimag(output(j))
    end do
end program multirate_assumed_rank_driver
