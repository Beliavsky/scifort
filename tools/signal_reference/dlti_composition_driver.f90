! SPDX-License-Identifier: MIT
! Batch driver for SciFort discrete-time series and parallel connections.
program dlti_composition_driver
    use scifort_kinds, only : dp
    use scifort_signal_dlti, only : dlti_system, dlti_ss
    use scifort_signal_discrete_simulation, only : discrete_simulation_result, discrete_response_result
    use scifort_signal_lti_frequency, only : mimo_frequency_response
    implicit none
    complex(dp), allocatable :: a1(:,:),b1(:,:),c1(:,:),d1(:,:), &
        a2(:,:),b2(:,:),c2(:,:),d2(:,:),u(:,:)
    real(dp), allocatable :: w(:)
    type(dlti_system) :: first, second, combined
    type(discrete_simulation_result) :: output
    type(discrete_response_result) :: sr, ir
    type(mimo_frequency_response) :: fr
    integer :: ncase, cas, mode, n1, n2, nu1, ny1, nu2, ny2, nt, nw, nu, ny, i, j, k
    real(dp) :: dt
    read(*,*) ncase
    do cas=1,ncase
        read(*,*) mode,n1,n2,nu1,ny1,nu2,ny2,nt,nw,dt
        allocate(a1(n1,n1),b1(n1,nu1),c1(ny1,n1),d1(ny1,nu1), &
            a2(n2,n2),b2(n2,nu2),c2(ny2,n2),d2(ny2,nu2),w(nw))
        read(*,*) w
        call input_matrix(a1)
        call input_matrix(b1)
        call input_matrix(c1)
        call input_matrix(d1)
        call input_matrix(a2)
        call input_matrix(b2)
        call input_matrix(c2)
        call input_matrix(d2)
        first = dlti_ss(a1,b1,c1,d1,dt)
        second = dlti_ss(a2,b2,c2,d2,dt)
        if (mode == 0) then
            combined = first%series(second)
        else
            combined = first%parallel(second)
        end if
        nu = nu2
        ny = ny1
        allocate(u(nt,nu))
        call input_matrix(u)
        output = combined%output(u)
        sr = combined%step(n=nt)
        ir = combined%impulse(n=nt)
        fr = combined%freqresp_mimo(w=w,solver='direct')
        write(*,'(6(i0,1x))') cas,combined%status,output%status,sr%status,ir%status,fr%status
        if (combined%status == 0 .and. output%status == 0 .and. fr%status == 0 &
            .and. sr%status == 0 .and. ir%status == 0) then
            do k=1,nw
                do j=1,nu
                    do i=1,ny
                        call put(fr%response(i,j,k))
                    end do
                end do
            end do
            do k=1,nt
                do i=1,ny
                    call put(output%y(k,i))
                end do
            end do
            do k=1,nt
                do j=1,nu
                    do i=1,ny
                        call put(sr%y(k,i,j))
                        call put(ir%y(k,i,j))
                    end do
                end do
            end do
        end if
        deallocate(a1,b1,c1,d1,a2,b2,c2,d2,u,w)
    end do
contains
    subroutine input_matrix(a)
        complex(dp), intent(out) :: a(:,:) !! matrix in column-major order
        integer :: row,column
        do column=1,size(a,2)
            do row=1,size(a,1)
                read(*,*) a(row,column)
            end do
        end do
    end subroutine input_matrix
    subroutine put(z)
        complex(dp), intent(in) :: z !! complex value
        write(*,'(2(es25.16e3,1x))') real(z,dp),aimag(z)
    end subroutine put
end program dlti_composition_driver
