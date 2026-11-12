! SPDX-License-Identifier: MIT
! Standard-input reference driver for optional complex Schur MIMO plans.
program lti_schur_driver
    use scifort_signal, only : dp, prepare_mimo_frequency, mimo_frequency_plan, &
        freqresp_mimo, dfreqresp_mimo, mimo_frequency_response
    implicit none
    complex(dp), allocatable :: a(:,:), b(:,:), c(:,:), d(:,:)
    real(dp), allocatable :: w(:)
    real(dp) :: re, im
    integer :: ncases, ns, ni, no, nw, discrete, k, i, j, h, icase
    type(mimo_frequency_plan) :: plan
    type(mimo_frequency_response) :: out

    read(*,*) ncases
    do icase = 1,ncases
        read(*,*) ns,ni,no,nw,discrete
        allocate(a(ns,ns),b(ns,ni),c(no,ns),d(no,ni),w(nw))
        do j = 1,ns
            do i = 1,ns
                read(*,*) re, im
                a(i,j) = cmplx(re,im,dp)
            end do
        end do
        do j = 1,ni
            do i = 1,ns
                read(*,*) re, im
                b(i,j) = cmplx(re,im,dp)
            end do
        end do
        do j = 1,ns
            do i = 1,no
                read(*,*) re, im
                c(i,j) = cmplx(re,im,dp)
            end do
        end do
        do j = 1,ni
            do i = 1,no
                read(*,*) re, im
                d(i,j) = cmplx(re,im,dp)
            end do
        end do
        read(*,*) w
        plan = prepare_mimo_frequency(a,b,c,d,solver='schur')
        if (discrete == 0) then
            out = freqresp_mimo(plan,w=w)
        else
            out = dfreqresp_mimo(plan,w=w)
        end if
        write(*,'(a,1x,i0)') trim(plan%method),out%status
        do k = 1,nw
            do j = 1,ni
                do i = 1,no
                    write(*,'(2(es25.16,1x))') real(out%response(i,j,k),dp), &
                        aimag(out%response(i,j,k))
                end do
            end do
        end do
        deallocate(a,b,c,d,w)
    end do
end program lti_schur_driver
