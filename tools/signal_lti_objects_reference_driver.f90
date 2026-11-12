! SPDX-License-Identifier: MIT
! Reference driver: simple formatted matrix I/O for Python/SciPy comparisons.
program signal_lti_objects_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_dlti, only : dlti_system,dlti_ss,dlti_balanced_result,dlti_success
    use scifort_signal_clti, only : clti_system,clti_ss,clti_success
    use scifort_signal_continuous_simulation, only : continuous_simulation_result
    use scifort_signal_lti_frequency, only : mimo_frequency_response
    implicit none
    complex(dp), allocatable :: a(:,:),b(:,:),c(:,:),d(:,:),u(:,:)
    real(dp) :: dt,w(5),tt(10)
    integer :: cases,k,mode,n,nu,ny,keep,i,j
    type(dlti_system) :: disc,dsys
    type(dlti_balanced_result) :: balanced
    type(clti_system) :: continuous
    type(mimo_frequency_response) :: frequency
    type(continuous_simulation_result) :: simulation

    read(*,*) cases
    w = [0.0_dp,0.3_dp,0.9_dp,2.4_dp,5.0_dp]
    do k=1,cases
        read(*,*) mode,n,nu,ny,keep,dt
        allocate(a(n,n),b(n,nu),c(ny,n),d(ny,nu))
        call read_mat(a)
        call read_mat(b)
        call read_mat(c)
        call read_mat(d)
        if (mode == 1) then
            dsys = dlti_ss(a,b,c,d,dt)
            balanced = dsys%balanced_truncate(keep)
            write(*,'(i0)') balanced%status
            if (balanced%status == dlti_success) then
                do i=1,n
                    write(*,'(es25.16)') balanced%hsv(i)
                end do
                write(*,'(es25.16)') balanced%error_bound
                call write_mat(balanced%model%a)
                call write_mat(balanced%model%b)
                call write_mat(balanced%model%c)
                call write_mat(balanced%model%d)
            end if
        else
            continuous = clti_ss(a,b,c,d)
            write(*,'(i0)') continuous%status()
            if (continuous%status() == clti_success) then
                frequency = continuous%freqresp_mimo(w=w,solver='direct')
                write(*,'(i0)') frequency%status
                if (frequency%status == 0) then
                    do i=1,size(w)
                        call write_mat(frequency%response(:,:,i))
                    end do
                end if
                do i=1,size(tt)
                    tt(i) = dt*real(i-1,dp)
                end do
                allocate(u(size(tt),nu))
                do i=1,size(tt)
                    do j=1,nu
                        u(i,j) = 0.5_dp*cmplx(cos(0.3_dp*i+0.2_dp*j), &
                            sin(0.2_dp*i-0.12_dp*j),dp)
                    end do
                end do
                simulation = continuous%output(u,tt,interp=.false.)
                write(*,'(i0)') simulation%status
                if (simulation%status == 0) call write_mat(simulation%y)
                disc = continuous%to_discrete(dt)
                write(*,'(i0)') disc%status
                if (disc%status == dlti_success) then
                    call write_mat(disc%a)
                    call write_mat(disc%b)
                    call write_mat(disc%c)
                    call write_mat(disc%d)
                end if
                deallocate(u)
            end if
        end if
        deallocate(a,b,c,d)
    end do
contains
    subroutine read_mat(x)
        complex(dp), intent(out) :: x(:,:) !! row-wise real/imaginary input
        real(dp) :: re,im
        integer :: p,q
        do p=1,size(x,1)
            do q=1,size(x,2)
                read(*,*) re,im
                x(p,q) = cmplx(re,im,dp)
            end do
        end do
    end subroutine read_mat
    subroutine write_mat(x)
        complex(dp), intent(in) :: x(:,:) !! row-wise complex values
        integer :: p,q
        do p=1,size(x,1)
            do q=1,size(x,2)
                write(*,'(es25.16,1x,es25.16)') real(x(p,q),dp),aimag(x(p,q))
            end do
        end do
    end subroutine write_mat
end program signal_lti_objects_reference_driver
