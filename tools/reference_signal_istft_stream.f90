! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! One-case stdin/stdout reference for SciPy comparison of streaming iSTFT.
program reference_signal_istft_stream
    use scifort_signal, only : dp, istft_stream, istft_stream_init, &
        istft_stream_push, istft_stream_finish, istft_status_invalid_input
    implicit none
    type(istft_stream) :: state
    integer :: nw, nf, no, nframes, onesided_flag, boundary_flag, real_flag
    integer :: nfreq, st, j, k
    real(dp) :: fs, re, im
    character(40) :: win, scale
    complex(dp), allocatable :: z(:,:), total_complex(:), complex_part(:)
    real(dp), allocatable :: total_real(:), real_part(:), times(:), part_times(:)

    read(*,*) nw, nf, no, nframes, onesided_flag, boundary_flag, real_flag, fs
    read(*,*) win, scale
    nfreq = nf
    if (onesided_flag==1) nfreq=nf/2+1
    allocate(z(nfreq,nframes), total_complex(0), total_real(0),times(0))
    do j=1,nframes
        do k=1,nfreq
            read(*,*) re,im
            z(k,j)=cmplx(re,im,dp)
        end do
    end do
    call istft_stream_init(state,nfreq,st,fs=fs,window=trim(win), &
        nperseg=nw,noverlap=no,nfft=nf,boundary=boundary_flag==1, &
        scaling=trim(scale),input_onesided=onesided_flag==1)
    if (st == istft_status_invalid_input) then
        print '(i0,1x,i0)', st, 0
        stop
    end if
    do j=1,nframes
        if (real_flag==1) then
            call istft_stream_push(state,z(:,j),real_part,st,part_times)
            total_real=[total_real,real_part]
        else
            call istft_stream_push(state,z(:,j),complex_part,st,part_times)
            total_complex=[total_complex,complex_part]
        end if
        if (st == istft_status_invalid_input) error stop 'invalid stream frame'
        times=[times,part_times]
    end do
    if (real_flag==1) then
        call istft_stream_finish(state,real_part,st,part_times)
        total_real=[total_real,real_part]
    else
        call istft_stream_finish(state,complex_part,st,part_times)
        total_complex=[total_complex,complex_part]
    end if
    times=[times,part_times]
    if (real_flag==1) then
        print '(i0,1x,i0)',st,size(total_real)
        do j=1,size(total_real)
            print '(3(es26.17,1x))',times(j),total_real(j),0.0_dp
        end do
    else
        print '(i0,1x,i0)',st,size(total_complex)
        do j=1,size(total_complex)
            print '(3(es26.17,1x))',times(j),real(total_complex(j),dp),aimag(total_complex(j))
        end do
    end if
end program reference_signal_istft_stream
