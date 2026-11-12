! SPDX-License-Identifier: MIT
! Text protocol for independent scipy.signal.sosfilt comparisons.
program sos_simulation_driver
    use scifort_kinds, only : dp
    use scifort_signal_discrete_simulation, only : dlsim_sos, discrete_simulation_result
    implicit none
    integer :: ncases, icase, variant, nsec, nin, time_flag, state_flag
    integer :: i, j, nt
    real(dp) :: dt, re, im
    real(dp), allocatable :: r_sos(:,:), r_u(:), timestamps(:), r_zi(:,:)
    complex(dp), allocatable :: c_sos(:,:), c_u(:), c_zi(:,:)
    type(discrete_simulation_result) :: out

    read(*,*) ncases
    do icase=1,ncases
        read(*,*) variant,nsec,nin,time_flag,state_flag,dt
        allocate(r_sos(nsec,6),c_sos(nsec,6),r_u(nin),c_u(nin), &
            timestamps(nin),r_zi(nsec,2),c_zi(nsec,2))
        do i = 1,nsec
            do j = 1,6
                read(*,*) re,im
                c_sos(i,j) = cmplx(re,im,dp)
                r_sos(i,j) = re
            end do
        end do
        do i = 1,nin
            read(*,*) re,im
            c_u(i) = cmplx(re,im,dp)
            r_u(i) = re
        end do
        if (time_flag/=0) then
            do i = 1,nin
                read(*,*) timestamps(i)
            end do
        end if
        if (state_flag/=0) then
            do i = 1,nsec
                do j = 1,2
                    read(*,*) re,im
                    c_zi(i,j) = cmplx(re,im,dp)
                    r_zi(i,j) = re
                end do
            end do
        end if
        select case (variant)
        case (0)
            if (time_flag/=0 .and. state_flag/=0) then
                out = dlsim_sos(r_sos,r_u,dt,t=timestamps,zi=r_zi)
            else if (time_flag/=0) then
                out = dlsim_sos(r_sos,r_u,dt,t=timestamps)
            else if (state_flag/=0) then
                out = dlsim_sos(r_sos,r_u,dt,zi=r_zi)
            else
                out = dlsim_sos(r_sos,r_u,dt)
            end if
        case (1)
            if (time_flag/=0 .and. state_flag/=0) then
                out = dlsim_sos(r_sos,c_u,dt,t=timestamps,zi=c_zi)
            else if (time_flag/=0) then
                out = dlsim_sos(r_sos,c_u,dt,t=timestamps)
            else if (state_flag/=0) then
                out = dlsim_sos(r_sos,c_u,dt,zi=c_zi)
            else
                out = dlsim_sos(r_sos,c_u,dt)
            end if
        case (2)
            if (time_flag/=0 .and. state_flag/=0) then
                out = dlsim_sos(c_sos,r_u,dt,t=timestamps,zi=c_zi)
            else if (time_flag/=0) then
                out = dlsim_sos(c_sos,r_u,dt,t=timestamps)
            else if (state_flag/=0) then
                out = dlsim_sos(c_sos,r_u,dt,zi=c_zi)
            else
                out = dlsim_sos(c_sos,r_u,dt)
            end if
        case (3)
            if (time_flag/=0 .and. state_flag/=0) then
                out = dlsim_sos(c_sos,c_u,dt,t=timestamps,zi=c_zi)
            else if (time_flag/=0) then
                out = dlsim_sos(c_sos,c_u,dt,t=timestamps)
            else if (state_flag/=0) then
                out = dlsim_sos(c_sos,c_u,dt,zi=c_zi)
            else
                out = dlsim_sos(c_sos,c_u,dt)
            end if
        end select
        write(*,'(a,1x,i0,1x,i0)') 'CASE',icase,out%status
        if (out%status==0) then
            nt = size(out%t)
            write(*,'(a,1x,i0)') 'LEN',nt
            do i = 1,nt
                write(*,'(es26.17e3)') out%t(i)
            end do
            do i = 1,nt
                write(*,'(2es26.17e3)') real(out%y(i,1),dp),aimag(out%y(i,1))
            end do
            do i = 1,nt
                do j = 1,2*nsec
                    write(*,'(2es26.17e3)') real(out%x(i,j),dp),aimag(out%x(i,j))
                end do
            end do
            do i = 1,nsec
                do j = 1,2
                    write(*,'(2es26.17e3)') real(out%zf(i,j),dp),aimag(out%zf(i,j))
                end do
            end do
        end if
        deallocate(r_sos,c_sos,r_u,c_u,timestamps,r_zi,c_zi)
    end do
end program sos_simulation_driver
