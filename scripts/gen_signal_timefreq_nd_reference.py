#!/usr/bin/env python3
"""Generate the deterministic Fortran time-frequency reference driver."""
from pathlib import Path

root = Path(__file__).resolve().parents[1]
ranks = {2:[12,8], 3:[12,8,10], 6:[1,12,1,8,1,10], 14:[1]*11+[12,8,10]}
txt = '''! SPDX-License-Identifier: MIT
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
'''
for rank in ranks:
    txt += f'        case({rank})\n            call process{rank}(cid,method,xt,ax,wi,di,si,bi,pa,uo)\n'
txt += '''        case default
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

'''
for rk,dims in ranks.items():
    rankp=rk+1
    col=','.join(':' for _ in range(rk))
    outcol=','.join(':' for _ in range(rankp))
    txt+=f'''    subroutine process{rk}(cid,method,xt,ax,wi,di,si,bi,pa,uo)
        integer, intent(in) :: cid,method,xt,ax,wi,di,si,bi,pa,uo !! test configuration
        integer, parameter :: nd={rk}
        integer :: d(nd),e(nd+1),nf,nt,nall,nw,no,hop,st
        real(dp), allocatable :: r(:),f(:),t(:)
        complex(dp), allocatable :: c(:)
        real(dp), allocatable :: xr({col})
        complex(dp), allocatable :: xc({col})
        real(dp), allocatable :: power({outcol})
        complex(dp), allocatable :: spec({outcol})
        character(:), allocatable :: window,detrend,scaling,boundary
        d=[{','.join(map(str,dims))}]
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
        allocate(xr({','.join(f'd({i})' for i in range(1,rk+1))}))
        allocate(xc({','.join(f'd({i})' for i in range(1,rk+1))}))
        xr=reshape(r,d)
        xc=reshape(c,d)
        allocate(power({','.join(f'e({i})' for i in range(1,rankp+1))}))
        allocate(spec({','.join(f'e({i})' for i in range(1,rankp+1))}))
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
    end subroutine process{rk}

'''
txt += 'end program reference_signal_timefreq_nd\n'
p = root/'validation/reference_drivers/reference_signal_timefreq_nd.f90'
p.write_text(txt)
print(p,'lines',len(txt.splitlines()),'maxlen',max(map(len,txt.splitlines())))
