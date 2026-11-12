! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program filter_axes_driver
    use scifort_kinds, only : dp
    use scifort_signal_sos, only : sosfilt
    use scifort_signal_filtfilt, only : filtfilt, sosfiltfilt
    implicit none
    integer :: cases, it, op, kind, rank, n1,n2,n3, along, ns, haszi, mode_id, padlen, st
    integer :: i,j,k, sh(3), zs(4)
    real(dp) :: rv, iv
    complex(dp), allocatable :: sos(:,:), x(:,:,:), zi(:,:,:,:), values(:,:,:), delays(:,:,:,:)
    complex(dp), allocatable :: b(:),a(:)
    real(dp), allocatable :: ry2(:,:),ry3(:,:,:),rz3(:,:,:),rz4(:,:,:,:)
    complex(dp), allocatable :: cy2(:,:),cy3(:,:,:),cz3(:,:,:),cz4(:,:,:,:)
    character(len=8) :: mode

    read(*,*) cases
    do it=1,cases
        read(*,*) op,kind,rank,n1,n2,n3,along,ns,haszi,mode_id,padlen
        if (allocated(sos)) deallocate(sos,x,b,a,zi,values)
        if (allocated(delays)) deallocate(delays)
        allocate(sos(ns,6),x(n1,n2,n3),b(3),a(3))
        do j=1,6
            do i=1,ns
                read(*,*) rv,iv
                sos(i,j)=cmplx(rv,iv,dp)
            end do
        end do
        do i=1,3
            read(*,*) rv,iv
            b(i)=cmplx(rv,iv,dp)
        end do
        do i=1,3
            read(*,*) rv,iv
            a(i)=cmplx(rv,iv,dp)
        end do
        do k=1,n3
            do j=1,n2
                do i=1,n1
                    read(*,*) rv,iv
                    x(i,j,k)=cmplx(rv,iv,dp)
                end do
            end do
        end do
        sh=[n1,n2,n3]
        sh(along)=2
        allocate(zi(ns,sh(1),sh(2),sh(3)))
        if(haszi==1)then
            do k=1,sh(3)
                do j=1,sh(2)
                    do i=1,sh(1)
                        block
                            integer :: q
                            do q=1,ns
                                read(*,*) rv,iv
                                zi(q,i,j,k)=cmplx(rv,iv,dp)
                            end do
                        end block
                    end do
                end do
            end do
        else
            zi=cmplx(0.0_dp,0.0_dp,dp)
        end if
        select case(mode_id)
        case(1)
            mode='odd'
        case(2)
            mode='even'
        case(3)
            mode='constant'
        case default
            mode='none'
        end select
        zs=0
        allocate(values(n1,n2,n3))
        if(op==1)then
            if(rank==2)then
                if(kind==0)then
                    if(haszi==1)then
                        ry2=sosfilt(real(sos,dp),real(x(:,:,1),dp),zi=real(zi(:,:,:,1),dp), &
                            zf=rz3,axis=along,status=st)
                    else
                        ry2=sosfilt(real(sos,dp),real(x(:,:,1),dp),zf=rz3,axis=along,status=st)
                    end if
                    values(:,:,1)=cmplx(ry2,0.0_dp,dp)
                    zs=[shape(rz3),1]
                    allocate(delays(ns,sh(1),sh(2),1))
                    delays(:,:,:,1)=cmplx(rz3,0.0_dp,dp)
                else
                    if(haszi==1)then
                        cy2=sosfilt(sos,x(:,:,1),zi=zi(:,:,:,1),zf=cz3,axis=along,status=st)
                    else
                        cy2=sosfilt(sos,x(:,:,1),zf=cz3,axis=along,status=st)
                    end if
                    values(:,:,1)=cy2
                    zs=[shape(cz3),1]
                    allocate(delays(ns,sh(1),sh(2),1))
                    delays(:,:,:,1)=cz3
                end if
            else
                if(kind==0)then
                    if(haszi==1)then
                        ry3=sosfilt(real(sos,dp),real(x,dp),zi=real(zi,dp), &
                            zf=rz4,axis=along,status=st)
                    else
                        ry3=sosfilt(real(sos,dp),real(x,dp),zf=rz4,axis=along,status=st)
                    end if
                    values=cmplx(ry3,0.0_dp,dp)
                    zs=shape(rz4)
                    allocate(delays(ns,sh(1),sh(2),sh(3)))
                    delays=cmplx(rz4,0.0_dp,dp)
                else
                    if(haszi==1)then
                        cy3=sosfilt(sos,x,zi=zi,zf=cz4,axis=along,status=st)
                    else
                        cy3=sosfilt(sos,x,zf=cz4,axis=along,status=st)
                    end if
                    values=cy3
                    zs=shape(cz4)
                    allocate(delays(ns,sh(1),sh(2),sh(3)))
                    delays=cz4
                end if
            end if
        else if(op==2)then
            if(rank==2)then
                if(kind==0)then
                    ry2=filtfilt(real(b,dp),real(a,dp),real(x(:,:,1),dp), &
                        axis=along,padtype=trim(mode),padlen=padlen,status=st)
                    values(:,:,1)=cmplx(ry2,0.0_dp,dp)
                else
                    cy2=filtfilt(b,a,x(:,:,1),axis=along,padtype=trim(mode), &
                        padlen=padlen,status=st)
                    values(:,:,1)=cy2
                end if
            else
                if(kind==0)then
                    ry3=filtfilt(real(b,dp),real(a,dp),real(x,dp), &
                        axis=along,padtype=trim(mode),padlen=padlen,status=st)
                    values=cmplx(ry3,0.0_dp,dp)
                else
                    cy3=filtfilt(b,a,x,axis=along,padtype=trim(mode),padlen=padlen,status=st)
                    values=cy3
                end if
            end if
        else
            if(rank==2)then
                if(kind==0)then
                    ry2=sosfiltfilt(real(sos,dp),real(x(:,:,1),dp),axis=along, &
                        padtype=trim(mode),padlen=padlen,status=st)
                    values(:,:,1)=cmplx(ry2,0.0_dp,dp)
                else
                    cy2=sosfiltfilt(sos,x(:,:,1),axis=along, &
                        padtype=trim(mode),padlen=padlen,status=st)
                    values(:,:,1)=cy2
                end if
            else
                if(kind==0)then
                    ry3=sosfiltfilt(real(sos,dp),real(x,dp),axis=along, &
                        padtype=trim(mode),padlen=padlen,status=st)
                    values=cmplx(ry3,0.0_dp,dp)
                else
                    cy3=sosfiltfilt(sos,x,axis=along, &
                        padtype=trim(mode),padlen=padlen,status=st)
                    values=cy3
                end if
            end if
        end if
        write(*,'(8(i0,1x))') st,n1,n2,n3,zs
        do k=1,n3
            do j=1,n2
                do i=1,n1
                    write(*,'(2(es25.16,1x))') real(values(i,j,k),dp),aimag(values(i,j,k))
                end do
            end do
        end do
        if(op==1)then
            do k=1,sh(3)
                do j=1,sh(2)
                    do i=1,sh(1)
                        block
                            integer :: q
                            do q=1,ns
                                write(*,'(2(es25.16,1x))') real(delays(q,i,j,k),dp),aimag(delays(q,i,j,k))
                            end do
                        end block
                    end do
                end do
            end do
        end if
    end do
end program
