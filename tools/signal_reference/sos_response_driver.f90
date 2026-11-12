! SPDX-License-Identifier: MIT
! Input/output protocol for randomized SOS response comparisons with SciPy.
program sos_response_driver
    use scifort_kinds, only : dp
    use scifort_signal_discrete_simulation, only : dstep_sos, dimpulse_sos, &
        discrete_response_result
    implicit none
    integer :: cases, case_id, nsections, nsamples, iscomplex, isimpulse, has_zi, &
        stride, i, j
    real(dp) :: dt, br(6), bi(6), zr(2), zimag(2)
    complex(dp), allocatable :: sections(:,:), states(:,:)
    real(dp), allocatable :: real_sections(:,:), real_states(:,:)
    type(discrete_response_result) :: out

    read(*,*) cases
    do case_id=1,cases
        read(*,*) nsections, nsamples, dt, iscomplex, isimpulse, has_zi, stride
        allocate(sections(stride*nsections,6),states(stride*nsections,2))
        sections = (0.0_dp,0.0_dp)
        states = (0.0_dp,0.0_dp)
        do i=1,nsections
            read(*,*) (br(j),j=1,6),(bi(j),j=1,6)
            sections(1+(i-1)*stride,:) = cmplx(br,bi,dp)
        end do
        if (has_zi /= 0) then
            do i=1,nsections
                read(*,*) zr(1),zr(2),zimag(1),zimag(2)
                states(1+(i-1)*stride,:) = cmplx(zr,zimag,dp)
            end do
        end if
        if (iscomplex /= 0) then
            if (isimpulse /= 0) then
                if (has_zi /= 0) then
                    out = dimpulse_sos(sections(1::stride,:),dt,n=nsamples,zi=states(1::stride,:))
                else
                    out = dimpulse_sos(sections(1::stride,:),dt,n=nsamples)
                end if
            else
                if (has_zi /= 0) then
                    out = dstep_sos(sections(1::stride,:),dt,n=nsamples,zi=states(1::stride,:))
                else
                    out = dstep_sos(sections(1::stride,:),dt,n=nsamples)
                end if
            end if
        else
            real_sections = real(sections,dp)
            real_states = real(states,dp)
            if (isimpulse /= 0) then
                if (has_zi /= 0) then
                    out = dimpulse_sos(real_sections(1::stride,:),dt,n=nsamples,zi=real_states(1::stride,:))
                else
                    out = dimpulse_sos(real_sections(1::stride,:),dt,n=nsamples)
                end if
            else
                if (has_zi /= 0) then
                    out = dstep_sos(real_sections(1::stride,:),dt,n=nsamples,zi=real_states(1::stride,:))
                else
                    out = dstep_sos(real_sections(1::stride,:),dt,n=nsamples)
                end if
            end if
        end if
        write(*,'(i0,1x,i0)') out%status, nsamples
        if (out%status == 0) then
            do i=1,nsamples
                write(*,'(3(es25.16e3,1x))') out%t(i),real(out%y(i,1,1),dp),aimag(out%y(i,1,1))
            end do
        end if
        deallocate(sections,states)
        if (allocated(real_sections)) deallocate(real_sections)
        if (allocated(real_states)) deallocate(real_states)
    end do
end program sos_response_driver
