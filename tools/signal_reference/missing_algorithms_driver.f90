! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Deterministic input families for the SciPy cross-language signal tests.
program missing_algorithms_driver
    use scifort_kinds, only : dp
    use scifort_signal_stft_dual, only : closest_stft_dual_window
    use scifort_signal_czt, only : czt_points
    use scifort_signal_fir2d, only : firwin_2d
    use scifort_signal_unique_roots, only : unique_roots
    use scifort_signal_analog, only : findfreqs
    use scifort_signal_symiir, only : symiirorder1, symiirorder2
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: sizes(4) = [8,11,16,23]
    integer, parameter :: hops(6) = [1,2,3,4,5,7]
    integer, parameter :: dims(2,5) = reshape([5,5,5,7,7,5,9,11,11,13],[2,5])
    character(len=9), parameter :: window_names(4) = ['hamming  ','hann     ','blackman ','kaiser   ']
    integer :: i, j, n, hop, dflag, cflag, sflag, st, t, ifilt, imod, ir, ic
    integer :: hsize(2)
    complex(dp), allocatable :: input(:), target(:), dual(:), points(:), roots(:), unique(:)
    complex(dp), allocatable :: response_complex(:), response_complex2(:,:)
    real(dp), allocatable :: rinput(:), rtarget(:), rdual(:), taps(:,:), omega(:), response(:), response2(:,:)
    complex(dp) :: alpha, zp(8), w, a
    real(dp) :: alpha_real, par1, par2
    complex(dp), allocatable :: signals(:), signals2(:,:)
    integer, allocatable :: multiplicities(:)
    character(len=8) :: method
    real(dp), allocatable :: cutoffs(:)

    do t = 1,size(sizes)
        n = sizes(t)
        allocate(input(n),target(n),rinput(n),rtarget(n))
        do i = 1,n
            input(i) = cmplx(0.1_dp+sin(0.17_dp*real(i,dp))**2, &
                0.2_dp*cos(0.19_dp*real(i,dp)),dp)
            target(i) = cmplx(0.7_dp+0.25_dp*cos(0.27_dp*real(i,dp)), &
                0.4_dp*sin(0.23_dp*real(i,dp)),dp)
        end do
        rinput = real(input,dp)
        rtarget = real(target,dp)
        do j = 1,size(hops)
            hop = hops(j)
            do cflag = 0,1
                do dflag = 0,1
                    do sflag = 0,1
                        if (cflag == 0) then
                            if (dflag == 1) then
                                call closest_stft_dual_window(rinput,hop,rdual,alpha_real, &
                                    desired_dual=rtarget,scaled=(sflag==1),status=st)
                            else
                                call closest_stft_dual_window(rinput,hop,rdual,alpha_real, &
                                    scaled=(sflag==1),status=st)
                            end if
                            dual = cmplx(rdual,0.0_dp,dp)
                            alpha = cmplx(alpha_real,0.0_dp,dp)
                        else
                            if (dflag == 1) then
                                call closest_stft_dual_window(input,hop,dual,alpha, &
                                    desired_dual=target,scaled=(sflag==1),status=st)
                            else
                                call closest_stft_dual_window(input,hop,dual,alpha, &
                                    scaled=(sflag==1),status=st)
                            end if
                        end if
                        write(*,'(a,7(1x,i0),2(1x,es24.16))',advance='no') &
                            'D',n,hop,cflag,dflag,sflag,st,size(dual),real(alpha,dp),aimag(alpha)
                        call show_complex(dual)
                    end do
                end do
            end do
        end do
        deallocate(input,target,rinput,rtarget)
    end do

    do t = 1,size(dims,2)
        hsize = dims(:,t)
        do ir = 1,2
            do ifilt = 1,size(window_names)
                do imod = 0,1
                    if (ir == 2 .and. imod == 1) cycle
                    if (imod == 0) then
                        cutoffs = [0.10_dp]
                        method = 'lowpass'
                    else
                        cutoffs = [0.17_dp,0.32_dp]
                        method = 'bandpass'
                    end if
                    ! For even lengths SciPy's symmetric firwin refuses some
                    ! passbands; keep this reference strictly on odd lengths.
                    if (ifilt == 4) then
                        taps = firwin_2d(hsize,trim(window_names(ifilt)),cutoffs, &
                            window_parameter=5.0_dp,circular=(ir==2),pass_zero=trim(method),status=st)
                    else
                        taps = firwin_2d(hsize,trim(window_names(ifilt)),cutoffs, &
                            circular=(ir==2),pass_zero=trim(method),status=st)
                    end if
                    write(*,'(a,8(1x,i0))',advance='no') 'F',hsize,ir,ifilt,imod,st, &
                        size(taps,1),size(taps,2)
                    call show_real(reshape(taps,[size(taps)]))
                end do
            end do
        end do
    end do

    do t = 1,2
        do ir = 1,3
            do imod = 1,3
                if (imod == 1) method = 'min'
                if (imod == 2) method = 'max'
                if (imod == 3) method = 'avg'
                do i = 1,8
                    zp(i) = cmplx(real(i/4,dp)+0.007_dp*real(mod(i,4),dp), &
                        0.007_dp*real(mod(i+1,4),dp),dp)
                end do
                if (t == 1) zp = cmplx(real(zp,dp),0.0_dp,dp)
                call unique_roots(zp,unique,multiplicities, &
                    tol=0.002_dp*real(1+5*(ir-1),dp),rtype=trim(method),status=st)
                write(*,'(a,5(1x,i0))',advance='no') 'U',t,ir,imod,st,size(unique)
                do i=1,size(unique)
                    write(*,'(2(1x,es24.16),1x,i0)',advance='no') &
                        real(unique(i),dp),aimag(unique(i)),multiplicities(i)
                end do
                write(*,*)
            end do
        end do
    end do

    do t = 1,5
        select case(t)
        case (1)
            omega = findfreqs([1.0_dp,0.0_dp],[1.0_dp,8.0_dp,25.0_dp],13,status=st)
        case (2)
            omega = findfreqs([1.0_dp,2.0_dp],[1.0_dp,2.0_dp,3.0_dp],17,status=st)
        case (3)
            omega = findfreqs([cmplx(0.0_dp,0.0_dp,dp)], &
                [cmplx(-1.0_dp,2.0_dp,dp),cmplx(-1.0_dp,-2.0_dp,dp)],11,kind='zp',status=st)
        case (4)
            omega = findfreqs([cmplx(-2.0_dp,1.0_dp,dp)], &
                [cmplx(-2.0_dp,2.0_dp,dp),cmplx(-2.0_dp,-2.0_dp,dp)],12,kind='zp',status=st)
        case (5)
            omega = findfreqs([cmplx(0.0_dp,0.0_dp,dp)], &
                [cmplx(-0.1_dp,0.0_dp,dp)],7,kind='zp',status=st)
        end select
        write(*,'(a,3(1x,i0))',advance='no') 'W',t,st,size(omega)
        call show_real(omega)
    end do

    do t = 1,4
        n = 1+5*t
        do ir = 0,1
            w = exp(cmplx(0.0_dp,-0.19_dp,dp))
            a = cmplx(1.1_dp,0.4_dp,dp)
            if (ir == 0) then
                points = czt_points(n,status=st)
            else
                points = czt_points(n,w=w,a=a,status=st)
            end if
            write(*,'(a,4(1x,i0))',advance='no') 'C',n,ir,st,size(points)
            call show_complex(points)
        end do
    end do


    allocate(signals(120),signals2(2,120))
    do i = 1,120
        signals(i) = cmplx(sin(0.1_dp*real(i,dp))+0.1_dp*cos(0.023_dp*real(i,dp)), &
            0.2_dp*cos(0.07_dp*real(i,dp))-0.08_dp*sin(0.03_dp*real(i,dp)),dp)
    end do
    do i = 1,2
        signals2(i,:) = signals*cmplx(1.0_dp+0.25_dp*real(i,dp),0.03_dp*real(i,dp),dp)
    end do
    do ifilt = 1,2
        do t = 1,3
            if (ifilt == 1) then
                select case (t)
                case (1)
                    par1 = 0.64_dp
                    par2 = 0.2_dp
                case (2)
                    par1 = 0.25_dp
                    par2 = 0.5_dp
                case (3)
                    par1 = 1.69_dp
                    par2 = -0.3_dp
                end select
            else
                select case (t)
                case (1)
                    par1 = 0.2_dp
                    par2 = 0.4_dp
                case (2)
                    par1 = 0.5_dp
                    par2 = 0.8_dp
                case (3)
                    par1 = 0.7_dp
                    par2 = 1.2_dp
                end select
            end if
            do cflag = 0,1
                do dflag = 1,2
                    if (ifilt == 1) then
                        if (cflag == 1) then
                            if (dflag == 1) then
                                response_complex = symiirorder1(signals,par1,par2, &
                                    precision=1.0e-12_dp,status=st)
                            else
                                response_complex2 = symiirorder1(signals2,par1,par2, &
                                    precision=1.0e-12_dp,status=st)
                            end if
                        else
                            if (dflag == 1) then
                                response = symiirorder1(real(signals,dp),par1,par2, &
                                    precision=1.0e-12_dp,status=st)
                            else
                                response2 = symiirorder1(real(signals2,dp),par1,par2, &
                                    precision=1.0e-12_dp,status=st)
                            end if
                        end if
                    else
                        if (cflag == 1) then
                            if (dflag == 1) then
                                response_complex = symiirorder2(signals,par1,par2, &
                                    precision=1.0e-12_dp,status=st)
                            else
                                response_complex2 = symiirorder2(signals2,par1,par2, &
                                    precision=1.0e-12_dp,status=st)
                            end if
                        else
                            if (dflag == 1) then
                                response = symiirorder2(real(signals,dp),par1,par2, &
                                    precision=1.0e-12_dp,status=st)
                            else
                                response2 = symiirorder2(real(signals2,dp),par1,par2, &
                                    precision=1.0e-12_dp,status=st)
                            end if
                        end if
                    end if
                    write(*,'(a,5(1x,i0))',advance='no') 'S',ifilt,t,cflag,dflag,st
                    if (cflag == 0) then
                        if (dflag == 1) then
                            write(*,'(2(1x,i0))',advance='no') 1,size(response)
                            call show_complex(cmplx(response,0.0_dp,dp))
                        else
                            write(*,'(2(1x,i0))',advance='no') size(response2,1),size(response2,2)
                            call show_complex(cmplx(reshape(response2,[size(response2)]),0.0_dp,dp))
                        end if
                    else
                        if (dflag == 1) then
                            write(*,'(2(1x,i0))',advance='no') 1,size(response_complex)
                            call show_complex(response_complex)
                        else
                            write(*,'(2(1x,i0))',advance='no') size(response_complex2,1),size(response_complex2,2)
                            call show_complex(reshape(response_complex2,[size(response_complex2)]))
                        end if
                    end if
                end do
            end do
        end do
    end do

contains

    subroutine show_real(x)
        real(dp), intent(in) :: x(:) !! data to serialize
        integer :: j
        do j = 1,size(x)
            write(*,'(1x,es24.16)',advance='no') x(j)
        end do
        write(*,*)
    end subroutine show_real

    subroutine show_complex(x)
        complex(dp), intent(in) :: x(:) !! data to serialize
        integer :: j
        do j = 1,size(x)
            write(*,'(2(1x,es24.16))',advance='no') real(x(j),dp),aimag(x(j))
        end do
        write(*,*)
    end subroutine show_complex
end program missing_algorithms_driver
