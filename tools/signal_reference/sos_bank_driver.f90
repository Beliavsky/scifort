! SPDX-License-Identifier: MIT
! Text protocol for SciPy comparisons of parallel SOS filter banks.
program sos_bank_driver
    use scifort_kinds, only : dp
    use scifort_signal_discrete_simulation, only : dlsim_sos_bank, discrete_sos_bank_result
    implicit none
    integer :: ncases, icase, variant, nsec, nin, nch, time_flag, state_flag
    integer :: i, j, k, nt
    real(dp) :: dt, re, im
    real(dp), allocatable :: rsos(:,:,:), ru(:,:), times(:), rzi(:,:,:)
    complex(dp), allocatable :: csos(:,:,:), cu(:,:), czi(:,:,:)
    type(discrete_sos_bank_result) :: out

    read(*,*) ncases
    do icase = 1, ncases
        read(*,*) variant,nsec,nin,nch,time_flag,state_flag,dt
        allocate(rsos(nsec,6,nch),ru(nin,nch),rzi(nsec,2,nch), &
            csos(nsec,6,nch),cu(nin,nch),czi(nsec,2,nch))
        do k = 1, nch
            do i = 1, nsec
                do j = 1, 6
                    read(*,*) re, im
                    rsos(i,j,k) = re
                    csos(i,j,k) = cmplx(re,im,dp)
                end do
            end do
        end do
        do k = 1, nch
            do i = 1, nin
                read(*,*) re, im
                ru(i,k) = re
                cu(i,k) = cmplx(re,im,dp)
            end do
        end do
        if (time_flag /= 0) then
            allocate(times(nin))
            do i = 1, nin
                read(*,*) times(i)
            end do
        end if
        if (state_flag /= 0) then
            do k = 1, nch
                do i = 1, nsec
                    do j = 1, 2
                        read(*,*) re, im
                        rzi(i,j,k) = re
                        czi(i,j,k) = cmplx(re,im,dp)
                    end do
                end do
            end do
        end if
        select case (variant)
        case (0:3)
            select case (variant)
            case (0)
                if (state_flag /= 0) then
                    out = dlsim_sos_bank(rsos(:,:,1),ru,dt,t=times,zi=rzi)
                else
                    out = dlsim_sos_bank(rsos(:,:,1),ru,dt,t=times)
                end if
            case (1)
                if (state_flag /= 0) then
                    out = dlsim_sos_bank(rsos(:,:,1),cu,dt,t=times,zi=czi)
                else
                    out = dlsim_sos_bank(rsos(:,:,1),cu,dt,t=times)
                end if
            case (2)
                if (state_flag /= 0) then
                    out = dlsim_sos_bank(csos(:,:,1),ru,dt,t=times,zi=czi)
                else
                    out = dlsim_sos_bank(csos(:,:,1),ru,dt,t=times)
                end if
            case (3)
                if (state_flag /= 0) then
                    out = dlsim_sos_bank(csos(:,:,1),cu,dt,t=times,zi=czi)
                else
                    out = dlsim_sos_bank(csos(:,:,1),cu,dt,t=times)
                end if
            end select
        case (4:7)
            select case (variant)
            case (4)
                if (state_flag /= 0) then
                    out = dlsim_sos_bank(rsos,ru,dt,t=times,zi=rzi)
                else
                    out = dlsim_sos_bank(rsos,ru,dt,t=times)
                end if
            case (5)
                if (state_flag /= 0) then
                    out = dlsim_sos_bank(rsos,cu,dt,t=times,zi=czi)
                else
                    out = dlsim_sos_bank(rsos,cu,dt,t=times)
                end if
            case (6)
                if (state_flag /= 0) then
                    out = dlsim_sos_bank(csos,ru,dt,t=times,zi=czi)
                else
                    out = dlsim_sos_bank(csos,ru,dt,t=times)
                end if
            case (7)
                if (state_flag /= 0) then
                    out = dlsim_sos_bank(csos,cu,dt,t=times,zi=czi)
                else
                    out = dlsim_sos_bank(csos,cu,dt,t=times)
                end if
            end select
        end select
        write(*,'(a,1x,i0,1x,i0)') 'CASE',icase,out%status
        if (out%status == 0) then
            nt = size(out%t)
            write(*,'(a,1x,i0)') 'LEN',nt
            do i = 1, nt
                write(*,'(es26.17e3)') out%t(i)
            end do
            do k = 1, nch
                do i = 1, nt
                    write(*,'(2es26.17e3)') real(out%y(i,k),dp),aimag(out%y(i,k))
                end do
            end do
            do k = 1, nch
                do i = 1, nt
                    do j = 1, 2*nsec
                        write(*,'(2es26.17e3)') real(out%x(i,j,k),dp),aimag(out%x(i,j,k))
                    end do
                end do
            end do
            do k = 1, nch
                do i = 1, nsec
                    do j = 1, 2
                        write(*,'(2es26.17e3)') real(out%zf(i,j,k),dp),aimag(out%zf(i,j,k))
                    end do
                end do
            end do
        end if
        deallocate(rsos,ru,rzi,csos,cu,czi)
        if (allocated(times)) deallocate(times)
    end do
end program sos_bank_driver
