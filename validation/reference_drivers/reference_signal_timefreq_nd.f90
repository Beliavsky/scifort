! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Generated with scripts/gen_signal_timefreq_nd_reference.py.
program reference_signal_timefreq_nd
    use scifort_signal, only : dp, stft_into, spectrogram_into
    implicit none
    integer :: ui,uo,ios,cid,rk,method,xt,ax,wi,di,si,bi,pa
    character(len=256) :: arg1,arg2
    call get_command_argument(1,arg1)
    call get_command_argument(2,arg2)
    if (len_trim(arg1)==0 .or. len_trim(arg2)==0) error stop 'usage: cases.txt values.bin'
    open(newunit=ui,file=trim(arg1),status='old',action='read')
    open(newunit=uo,file=trim(arg2),status='replace',access='stream',form='unformatted')
    do
        read(ui,*,iostat=ios) cid,rk,method,xt,ax,wi,di,si,bi,pa
        if (ios<0) exit
        if (ios/=0) error stop 'malformed cases'
        select case(rk)
        case(2)
            call process2(cid,method,xt,ax,wi,di,si,bi,pa,uo)
        case(3)
            call process3(cid,method,xt,ax,wi,di,si,bi,pa,uo)
        case(6)
            call process6(cid,method,xt,ax,wi,di,si,bi,pa,uo)
        case(14)
            call process14(cid,method,xt,ax,wi,di,si,bi,pa,uo)
        case default
            error stop 'unknown rank'
        end select
    end do
    close(ui)
    close(uo)
contains
    subroutine samples(cid,r,c)
        integer, intent(in) :: cid !! case id
        real(dp), intent(out) :: r(:) !! real samples
        complex(dp), intent(out) :: c(:) !! complex samples
        integer :: i
        do i=1,size(r)
            r(i)=sin(0.071_dp*i+0.013_dp*cid)+0.2_dp*cos(0.047_dp*i-0.011_dp*cid)
            c(i)=cmplx(r(i),0.22_dp*cos(0.081_dp*i+0.017_dp*cid),dp)
        end do
    end subroutine samples

    subroutine names(wi,di,si,bi,window,detrend,scaling,boundary)
        integer, intent(in) :: wi,di,si,bi !! options
        character(:), allocatable, intent(out) :: window,detrend,scaling,boundary
        select case(wi)
        case(0)
            window='boxcar'
        case(1)
            window='hann'
        case default
            window='tukey'
        end select
        select case(di)
        case(0)
            detrend='constant'
        case(1)
            detrend='linear'
        case default
            detrend='none'
        end select
        scaling='density'
        if (si==1) scaling='spectrum'
        boundary='none'
        if (bi==1) boundary='zeros'
    end subroutine names

    subroutine process2(cid,method,xt,ax,wi,di,si,bi,pa,uo)
        integer, intent(in) :: cid,method,xt,ax,wi,di,si,bi,pa,uo !! test configuration
        integer, parameter :: nd=2
        integer :: d(nd),e(nd+1),nf,nt,nall,nw,no,hop,st
        real(dp), allocatable :: r(:),f(:),t(:)
        complex(dp), allocatable :: c(:)
        real(dp), allocatable :: xr(:,:)
        complex(dp), allocatable :: xc(:,:)
        real(dp), allocatable :: power(:,:,:)
        complex(dp), allocatable :: spec(:,:,:)
        character(:), allocatable :: window,detrend,scaling,boundary
        d=[12,8]
        nf=9
        if (xt==1) nf=16
        nw=6
        no=2
        hop=nw-no
        nall=d(ax)
        if (method==1 .and. bi==1) nall=nall+nw
        nt=1+(nall-nw)/hop
        if (method==1 .and. pa==1) nt=1+(nall-nw+hop-1)/hop
        e=[d(:ax-1),nf,d(ax+1:),nt]
        allocate(r(product(d)),c(product(d)))
        call samples(cid,r,c)
        allocate(xr(d(1),d(2)))
        allocate(xc(d(1),d(2)))
        xr=reshape(r,d)
        xc=reshape(c,d)
        allocate(power(e(1),e(2),e(3)))
        allocate(spec(e(1),e(2),e(3)))
        power=-900.0_dp
        spec=(-900.0_dp,20.0_dp)
        call names(wi,di,si,bi,window,detrend,scaling,boundary)
        select case(method)
        case(1)
            if (si==0) scaling='psd'
            if (xt==0) then
                call stft_into(xr,spec,f,t,fs=7.0_dp,window=window,nperseg=nw, &
                    noverlap=no,nfft=16,detrend=detrend,boundary=boundary, &
                    padded=(pa==1),scaling=scaling,axis=ax,status=st)
            else
                call stft_into(xc,spec,f,t,fs=7.0_dp,window=window,nperseg=nw, &
                    noverlap=no,nfft=16,detrend=detrend,boundary=boundary, &
                    padded=(pa==1),scaling=scaling,axis=ax,status=st)
            end if
        case(2:6)
            if (method==3) then
                if (xt==0) then
                    call spectrogram_into(xr,spec,f,t,fs=7.0_dp,window=window, &
                        nperseg=nw,noverlap=no,nfft=16,detrend=detrend,scaling=scaling, &
                        mode='complex',axis=ax,status=st)
                else
                    call spectrogram_into(xc,spec,f,t,fs=7.0_dp,window=window, &
                        nperseg=nw,noverlap=no,nfft=16,detrend=detrend,scaling=scaling, &
                        mode='complex',axis=ax,status=st)
                end if
            else
                select case(method)
                case(2)
                    boundary='psd'
                case(4)
                    boundary='magnitude'
                case(5)
                    boundary='angle'
                case(6)
                    boundary='phase'
                end select
                if (xt==0) then
                    call spectrogram_into(xr,power,f,t,fs=7.0_dp,window=window, &
                        nperseg=nw,noverlap=no,nfft=16,detrend=detrend,scaling=scaling, &
                        mode=boundary,axis=ax,status=st)
                else
                    call spectrogram_into(xc,power,f,t,fs=7.0_dp,window=window, &
                        nperseg=nw,noverlap=no,nfft=16,detrend=detrend,scaling=scaling, &
                        mode=boundary,axis=ax,status=st)
                end if
            end if
        case default
            error stop 'bad method'
        end select
        if (st/=0 .or. size(f)/=nf .or. size(t)/=nt) then
            print *, 'failed case',cid,nd,method,xt,ax,st,size(f),size(t),nf,nt
            error stop 'time-frequency reference failure'
        end if
        write(uo) st
        write(uo) f
        write(uo) t
        if (method==1 .or. method==3) then
            write(uo) reshape(spec,[size(spec)])
        else
            write(uo) cmplx(reshape(power,[size(power)]),0.0_dp,dp)
        end if
    end subroutine process2

    subroutine process3(cid,method,xt,ax,wi,di,si,bi,pa,uo)
        integer, intent(in) :: cid,method,xt,ax,wi,di,si,bi,pa,uo !! test configuration
        integer, parameter :: nd=3
        integer :: d(nd),e(nd+1),nf,nt,nall,nw,no,hop,st
        real(dp), allocatable :: r(:),f(:),t(:)
        complex(dp), allocatable :: c(:)
        real(dp), allocatable :: xr(:,:,:)
        complex(dp), allocatable :: xc(:,:,:)
        real(dp), allocatable :: power(:,:,:,:)
        complex(dp), allocatable :: spec(:,:,:,:)
        character(:), allocatable :: window,detrend,scaling,boundary
        d=[12,8,10]
        nf=9
        if (xt==1) nf=16
        nw=6
        no=2
        hop=nw-no
        nall=d(ax)
        if (method==1 .and. bi==1) nall=nall+nw
        nt=1+(nall-nw)/hop
        if (method==1 .and. pa==1) nt=1+(nall-nw+hop-1)/hop
        e=[d(:ax-1),nf,d(ax+1:),nt]
        allocate(r(product(d)),c(product(d)))
        call samples(cid,r,c)
        allocate(xr(d(1),d(2),d(3)))
        allocate(xc(d(1),d(2),d(3)))
        xr=reshape(r,d)
        xc=reshape(c,d)
        allocate(power(e(1),e(2),e(3),e(4)))
        allocate(spec(e(1),e(2),e(3),e(4)))
        power=-900.0_dp
        spec=(-900.0_dp,20.0_dp)
        call names(wi,di,si,bi,window,detrend,scaling,boundary)
        select case(method)
        case(1)
            if (si==0) scaling='psd'
            if (xt==0) then
                call stft_into(xr,spec,f,t,fs=7.0_dp,window=window,nperseg=nw, &
                    noverlap=no,nfft=16,detrend=detrend,boundary=boundary, &
                    padded=(pa==1),scaling=scaling,axis=ax,status=st)
            else
                call stft_into(xc,spec,f,t,fs=7.0_dp,window=window,nperseg=nw, &
                    noverlap=no,nfft=16,detrend=detrend,boundary=boundary, &
                    padded=(pa==1),scaling=scaling,axis=ax,status=st)
            end if
        case(2:6)
            if (method==3) then
                if (xt==0) then
                    call spectrogram_into(xr,spec,f,t,fs=7.0_dp,window=window, &
                        nperseg=nw,noverlap=no,nfft=16,detrend=detrend,scaling=scaling, &
                        mode='complex',axis=ax,status=st)
                else
                    call spectrogram_into(xc,spec,f,t,fs=7.0_dp,window=window, &
                        nperseg=nw,noverlap=no,nfft=16,detrend=detrend,scaling=scaling, &
                        mode='complex',axis=ax,status=st)
                end if
            else
                select case(method)
                case(2)
                    boundary='psd'
                case(4)
                    boundary='magnitude'
                case(5)
                    boundary='angle'
                case(6)
                    boundary='phase'
                end select
                if (xt==0) then
                    call spectrogram_into(xr,power,f,t,fs=7.0_dp,window=window, &
                        nperseg=nw,noverlap=no,nfft=16,detrend=detrend,scaling=scaling, &
                        mode=boundary,axis=ax,status=st)
                else
                    call spectrogram_into(xc,power,f,t,fs=7.0_dp,window=window, &
                        nperseg=nw,noverlap=no,nfft=16,detrend=detrend,scaling=scaling, &
                        mode=boundary,axis=ax,status=st)
                end if
            end if
        case default
            error stop 'bad method'
        end select
        if (st/=0 .or. size(f)/=nf .or. size(t)/=nt) then
            print *, 'failed case',cid,nd,method,xt,ax,st,size(f),size(t),nf,nt
            error stop 'time-frequency reference failure'
        end if
        write(uo) st
        write(uo) f
        write(uo) t
        if (method==1 .or. method==3) then
            write(uo) reshape(spec,[size(spec)])
        else
            write(uo) cmplx(reshape(power,[size(power)]),0.0_dp,dp)
        end if
    end subroutine process3

    subroutine process6(cid,method,xt,ax,wi,di,si,bi,pa,uo)
        integer, intent(in) :: cid,method,xt,ax,wi,di,si,bi,pa,uo !! test configuration
        integer, parameter :: nd=6
        integer :: d(nd),e(nd+1),nf,nt,nall,nw,no,hop,st
        real(dp), allocatable :: r(:),f(:),t(:)
        complex(dp), allocatable :: c(:)
        real(dp), allocatable :: xr(:,:,:,:,:,:)
        complex(dp), allocatable :: xc(:,:,:,:,:,:)
        real(dp), allocatable :: power(:,:,:,:,:,:,:)
        complex(dp), allocatable :: spec(:,:,:,:,:,:,:)
        character(:), allocatable :: window,detrend,scaling,boundary
        d=[1,12,1,8,1,10]
        nf=9
        if (xt==1) nf=16
        nw=6
        no=2
        hop=nw-no
        nall=d(ax)
        if (method==1 .and. bi==1) nall=nall+nw
        nt=1+(nall-nw)/hop
        if (method==1 .and. pa==1) nt=1+(nall-nw+hop-1)/hop
        e=[d(:ax-1),nf,d(ax+1:),nt]
        allocate(r(product(d)),c(product(d)))
        call samples(cid,r,c)
        allocate(xr(d(1),d(2),d(3),d(4),d(5),d(6)))
        allocate(xc(d(1),d(2),d(3),d(4),d(5),d(6)))
        xr=reshape(r,d)
        xc=reshape(c,d)
        allocate(power(e(1),e(2),e(3),e(4),e(5),e(6),e(7)))
        allocate(spec(e(1),e(2),e(3),e(4),e(5),e(6),e(7)))
        power=-900.0_dp
        spec=(-900.0_dp,20.0_dp)
        call names(wi,di,si,bi,window,detrend,scaling,boundary)
        select case(method)
        case(1)
            if (si==0) scaling='psd'
            if (xt==0) then
                call stft_into(xr,spec,f,t,fs=7.0_dp,window=window,nperseg=nw, &
                    noverlap=no,nfft=16,detrend=detrend,boundary=boundary, &
                    padded=(pa==1),scaling=scaling,axis=ax,status=st)
            else
                call stft_into(xc,spec,f,t,fs=7.0_dp,window=window,nperseg=nw, &
                    noverlap=no,nfft=16,detrend=detrend,boundary=boundary, &
                    padded=(pa==1),scaling=scaling,axis=ax,status=st)
            end if
        case(2:6)
            if (method==3) then
                if (xt==0) then
                    call spectrogram_into(xr,spec,f,t,fs=7.0_dp,window=window, &
                        nperseg=nw,noverlap=no,nfft=16,detrend=detrend,scaling=scaling, &
                        mode='complex',axis=ax,status=st)
                else
                    call spectrogram_into(xc,spec,f,t,fs=7.0_dp,window=window, &
                        nperseg=nw,noverlap=no,nfft=16,detrend=detrend,scaling=scaling, &
                        mode='complex',axis=ax,status=st)
                end if
            else
                select case(method)
                case(2)
                    boundary='psd'
                case(4)
                    boundary='magnitude'
                case(5)
                    boundary='angle'
                case(6)
                    boundary='phase'
                end select
                if (xt==0) then
                    call spectrogram_into(xr,power,f,t,fs=7.0_dp,window=window, &
                        nperseg=nw,noverlap=no,nfft=16,detrend=detrend,scaling=scaling, &
                        mode=boundary,axis=ax,status=st)
                else
                    call spectrogram_into(xc,power,f,t,fs=7.0_dp,window=window, &
                        nperseg=nw,noverlap=no,nfft=16,detrend=detrend,scaling=scaling, &
                        mode=boundary,axis=ax,status=st)
                end if
            end if
        case default
            error stop 'bad method'
        end select
        if (st/=0 .or. size(f)/=nf .or. size(t)/=nt) then
            print *, 'failed case',cid,nd,method,xt,ax,st,size(f),size(t),nf,nt
            error stop 'time-frequency reference failure'
        end if
        write(uo) st
        write(uo) f
        write(uo) t
        if (method==1 .or. method==3) then
            write(uo) reshape(spec,[size(spec)])
        else
            write(uo) cmplx(reshape(power,[size(power)]),0.0_dp,dp)
        end if
    end subroutine process6

    subroutine process14(cid,method,xt,ax,wi,di,si,bi,pa,uo)
        integer, intent(in) :: cid,method,xt,ax,wi,di,si,bi,pa,uo !! test configuration
        integer, parameter :: nd=14
        integer :: d(nd),e(nd+1),nf,nt,nall,nw,no,hop,st
        real(dp), allocatable :: r(:),f(:),t(:)
        complex(dp), allocatable :: c(:)
        real(dp), allocatable :: xr(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: xc(:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        real(dp), allocatable :: power(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        complex(dp), allocatable :: spec(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
        character(:), allocatable :: window,detrend,scaling,boundary
        d=[1,1,1,1,1,1,1,1,1,1,1,12,8,10]
        nf=9
        if (xt==1) nf=16
        nw=6
        no=2
        hop=nw-no
        nall=d(ax)
        if (method==1 .and. bi==1) nall=nall+nw
        nt=1+(nall-nw)/hop
        if (method==1 .and. pa==1) nt=1+(nall-nw+hop-1)/hop
        e=[d(:ax-1),nf,d(ax+1:),nt]
        allocate(r(product(d)),c(product(d)))
        call samples(cid,r,c)
        allocate(xr(d(1),d(2),d(3),d(4),d(5),d(6),d(7),d(8),d(9),d(10),d(11),d(12),d(13),d(14)))
        allocate(xc(d(1),d(2),d(3),d(4),d(5),d(6),d(7),d(8),d(9),d(10),d(11),d(12),d(13),d(14)))
        xr=reshape(r,d)
        xc=reshape(c,d)
        allocate(power(e(1),e(2),e(3),e(4),e(5),e(6),e(7),e(8),e(9),e(10),e(11),e(12),e(13),e(14),e(15)))
        allocate(spec(e(1),e(2),e(3),e(4),e(5),e(6),e(7),e(8),e(9),e(10),e(11),e(12),e(13),e(14),e(15)))
        power=-900.0_dp
        spec=(-900.0_dp,20.0_dp)
        call names(wi,di,si,bi,window,detrend,scaling,boundary)
        select case(method)
        case(1)
            if (si==0) scaling='psd'
            if (xt==0) then
                call stft_into(xr,spec,f,t,fs=7.0_dp,window=window,nperseg=nw, &
                    noverlap=no,nfft=16,detrend=detrend,boundary=boundary, &
                    padded=(pa==1),scaling=scaling,axis=ax,status=st)
            else
                call stft_into(xc,spec,f,t,fs=7.0_dp,window=window,nperseg=nw, &
                    noverlap=no,nfft=16,detrend=detrend,boundary=boundary, &
                    padded=(pa==1),scaling=scaling,axis=ax,status=st)
            end if
        case(2:6)
            if (method==3) then
                if (xt==0) then
                    call spectrogram_into(xr,spec,f,t,fs=7.0_dp,window=window, &
                        nperseg=nw,noverlap=no,nfft=16,detrend=detrend,scaling=scaling, &
                        mode='complex',axis=ax,status=st)
                else
                    call spectrogram_into(xc,spec,f,t,fs=7.0_dp,window=window, &
                        nperseg=nw,noverlap=no,nfft=16,detrend=detrend,scaling=scaling, &
                        mode='complex',axis=ax,status=st)
                end if
            else
                select case(method)
                case(2)
                    boundary='psd'
                case(4)
                    boundary='magnitude'
                case(5)
                    boundary='angle'
                case(6)
                    boundary='phase'
                end select
                if (xt==0) then
                    call spectrogram_into(xr,power,f,t,fs=7.0_dp,window=window, &
                        nperseg=nw,noverlap=no,nfft=16,detrend=detrend,scaling=scaling, &
                        mode=boundary,axis=ax,status=st)
                else
                    call spectrogram_into(xc,power,f,t,fs=7.0_dp,window=window, &
                        nperseg=nw,noverlap=no,nfft=16,detrend=detrend,scaling=scaling, &
                        mode=boundary,axis=ax,status=st)
                end if
            end if
        case default
            error stop 'bad method'
        end select
        if (st/=0 .or. size(f)/=nf .or. size(t)/=nt) then
            print *, 'failed case',cid,nd,method,xt,ax,st,size(f),size(t),nf,nt
            error stop 'time-frequency reference failure'
        end if
        write(uo) st
        write(uo) f
        write(uo) t
        if (method==1 .or. method==3) then
            write(uo) reshape(spec,[size(spec)])
        else
            write(uo) cmplx(reshape(power,[size(power)]),0.0_dp,dp)
        end if
    end subroutine process14

end program reference_signal_timefreq_nd
