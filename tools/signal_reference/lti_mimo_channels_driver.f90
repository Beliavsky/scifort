! SPDX-License-Identifier: MIT
! Stream probe for comparing selected prepared MIMO responses against SciPy.
program lti_mimo_channels_driver
    use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, iostat_end
    use scifort_signal, only : dp, prepare_mimo_frequency, mimo_frequency_plan, &
        freqresp_mimo, dfreqresp_mimo, mimo_frequency_response
    implicit none
    complex(dp), allocatable :: a(:,:), b(:,:), c(:,:), d(:,:)
    real(dp), allocatable :: w(:)
    integer, allocatable :: inputs(:), outputs(:)
    integer :: mode, ns, ni, no, nw, nsi, nso, ios, i, j, k
    real(dp) :: xr, xi
    type(mimo_frequency_plan) :: plan
    type(mimo_frequency_response) :: reduced, direct

    do
        read(input_unit,*,iostat=ios) mode,ns,ni,no,nw,nsi,nso
        if (ios == iostat_end) exit
        if (ios /= 0) error stop 'Invalid driver header'
        allocate(a(ns,ns),b(ns,ni),c(no,ns),d(no,ni),w(nw), &
            inputs(nsi),outputs(nso))
        do j=1,ns
            do i=1,ns
                read(input_unit,*) xr,xi
                a(i,j) = cmplx(xr,xi,dp)
            end do
        end do
        do j=1,ni
            do i=1,ns
                read(input_unit,*) xr,xi
                b(i,j) = cmplx(xr,xi,dp)
            end do
        end do
        do j=1,ns
            do i=1,no
                read(input_unit,*) xr,xi
                c(i,j) = cmplx(xr,xi,dp)
            end do
        end do
        do j=1,ni
            do i=1,no
                read(input_unit,*) xr,xi
                d(i,j) = cmplx(xr,xi,dp)
            end do
        end do
        do i=1,nw
            read(input_unit,*) w(i)
        end do
        do i=1,nsi
            read(input_unit,*) inputs(i)
        end do
        do i=1,nso
            read(input_unit,*) outputs(i)
        end do
        plan = prepare_mimo_frequency(a,b,c,d)
        if (mode == 1) then
            reduced = freqresp_mimo(plan,w=w,inputs=inputs,outputs=outputs)
            direct = freqresp_mimo(a,b,c,d,w=w,inputs=inputs,outputs=outputs,solver='direct')
        else
            reduced = dfreqresp_mimo(plan,w=w,inputs=inputs,outputs=outputs)
            direct = dfreqresp_mimo(a,b,c,d,w=w,inputs=inputs,outputs=outputs,solver='direct')
        end if
        write(output_unit,'(4(i0,1x))') plan%status,reduced%status,direct%status,size(w)
        do k=1,nw
            do j=1,nsi
                do i=1,nso
                    write(output_unit,'(4(es26.17e3,1x))') &
                        real(reduced%response(i,j,k),dp),aimag(reduced%response(i,j,k)), &
                        real(direct%response(i,j,k),dp),aimag(direct%response(i,j,k))
                end do
            end do
        end do
        deallocate(a,b,c,d,w,inputs,outputs)
    end do
end program lti_mimo_channels_driver
