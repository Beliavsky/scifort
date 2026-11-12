! SPDX-License-Identifier: MIT
! Read a batch of discrete-time LTI systems and report their numerical results.
program dlti_object_driver
    use scifort_kinds, only : dp
    use scifort_signal_dlti, only : dlti_tf, dlti_zpk, dlti_ss, dlti_system
    use scifort_signal_discrete_simulation, only : discrete_simulation_result, discrete_response_result
    use scifort_signal_lti_frequency, only : mimo_frequency_response, bode_result, mimo_bode_result
    use scifort_signal_design, only : frequency_response
    implicit none
    integer :: ncase, icase, mode, nt, nw, ns, ni, no, nb, na, nz, np, k, j, i
    real(dp) :: dt
    complex(dp), allocatable :: num(:), den(:), z(:), p(:), a(:,:), b(:,:), c(:,:), d(:,:), u(:,:)
    complex(dp) :: gain
    real(dp), allocatable :: w(:)
    type(dlti_system) :: system, converted
    type(discrete_simulation_result) :: response
    type(discrete_response_result) :: step_result, impulse_result
    type(frequency_response) :: fr
    type(mimo_frequency_response) :: fr_mimo
    type(bode_result) :: br
    type(mimo_bode_result) :: br_mimo

    read(*,*) ncase
    do icase = 1, ncase
        read(*,*) mode, nt, nw, ns, ni, no, nb, na, nz, np, dt
        allocate(w(nw),u(nt,ni))
        read(*,*) w
        select case (mode)
        case (1)
            allocate(num(nb),den(na))
            call get_vector(num)
            call get_vector(den)
            system = dlti_tf(num,den,dt=dt)
        case (2)
            allocate(z(nz),p(np))
            call get_vector(z)
            call get_vector(p)
            read(*,*) gain
            system = dlti_zpk(z,p,gain,dt=dt)
        case (3)
            allocate(a(ns,ns),b(ns,ni),c(no,ns),d(no,ni))
            call get_matrix(a)
            call get_matrix(b)
            call get_matrix(c)
            call get_matrix(d)
            system = dlti_ss(a,b,c,d,dt=dt)
        end select
        call get_matrix(u)
        response = system%output(u)
        step_result = system%step(n=nt)
        impulse_result = system%impulse(n=nt)
        converted = system%to_ss()
        write(*,'(7(i0,1x))') icase, system%status, response%status, &
            step_result%status, impulse_result%status, converted%status, system%nstates()
        if (system%status == 0 .and. response%status == 0 .and. &
            step_result%status == 0 .and. impulse_result%status == 0) then
            if (mode == 3 .and. (ni /= 1 .or. no /= 1)) then
                fr_mimo = system%freqresp_mimo(w=w,solver='direct')
                br_mimo = system%bode_mimo(w=w,solver='direct')
                write(*,'(2(i0,1x))') fr_mimo%status, br_mimo%status
                if (fr_mimo%status == 0 .and. br_mimo%status == 0) then
                    do k=1,nw
                        do j=1,ni
                            do i=1,no
                                call put_complex(fr_mimo%response(i,j,k))
                                write(*,'(2(es25.16e3,1x))') br_mimo%magnitude(i,j,k), br_mimo%phase(i,j,k)
                            end do
                        end do
                    end do
                end if
            else
                fr = system%freqresp(w=w)
                br = system%bode(w=w)
                write(*,'(2(i0,1x))') fr%status, br%status
                if (fr%status == 0 .and. br%status == 0) then
                    do k=1,nw
                        call put_complex(fr%response(k))
                        write(*,'(2(es25.16e3,1x))') br%magnitude(k), br%phase(k)
                    end do
                end if
            end if
            do k=1,nt
                do i=1,no
                    call put_complex(response%y(k,i))
                end do
            end do
            do k=1,nt
                do j=1,ni
                    do i=1,no
                        call put_complex(step_result%y(k,i,j))
                        call put_complex(impulse_result%y(k,i,j))
                    end do
                end do
            end do
        end if
        if (allocated(num)) deallocate(num,den)
        if (allocated(z)) deallocate(z,p)
        if (allocated(a)) deallocate(a,b,c,d)
        deallocate(w,u)
    end do
contains
    subroutine get_vector(v)
        complex(dp), intent(out) :: v(:) !! destination
        integer :: idx
        do idx=1,size(v)
            read(*,*) v(idx)
        end do
    end subroutine get_vector

    subroutine get_matrix(m)
        complex(dp), intent(out) :: m(:,:) !! matrix in column-major order
        integer :: row, column
        do column=1,size(m,2)
            do row=1,size(m,1)
                read(*,*) m(row,column)
            end do
        end do
    end subroutine get_matrix

    subroutine put_complex(v)
        complex(dp), intent(in) :: v !! output value
        write(*,'(2(es25.16e3,1x))') real(v,dp),aimag(v)
    end subroutine put_complex
end program dlti_object_driver
