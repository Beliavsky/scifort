! SPDX-License-Identifier: MIT
! Read MIMO systems and report feedback, append, selection, Gramians and Kalman matrices.
program dlti_interconnections_driver
    use scifort_kinds, only : dp
    use scifort_signal_dlti, only : dlti_system, dlti_ss, dlti_gramian_result
    implicit none
    type(dlti_system) :: plant,controller,feedback,appended,selected
    type(dlti_gramian_result) :: gr
    complex(dp), allocatable :: ap(:,:),bp(:,:),cp(:,:),dp1(:,:), &
        ak(:,:),bk(:,:),ck(:,:),dk(:,:),ctr(:,:),obs(:,:)
    integer :: ncase, cas, np, nk, nu, ny, sign, i
    real(dp) :: dt

    read(*,*) ncase
    do cas=1,ncase
        read(*,*) np,nk,nu,ny,sign,dt
        allocate(ap(np,np),bp(np,nu),cp(ny,np),dp1(ny,nu), &
            ak(nk,nk),bk(nk,ny),ck(nu,nk),dk(nu,ny))
        call input_matrix(ap)
        call input_matrix(bp)
        call input_matrix(cp)
        call input_matrix(dp1)
        call input_matrix(ak)
        call input_matrix(bk)
        call input_matrix(ck)
        call input_matrix(dk)
        plant = dlti_ss(ap,bp,cp,dp1,dt)
        controller = dlti_ss(ak,bk,ck,dk,dt)
        feedback = plant%feedback(controller,sign=sign)
        appended = plant%append(controller)
        selected = appended%select_io([nu+1,1,nu+1],[ny+1,1])
        gr = feedback%gramians()
        write(*,'(5(i0,1x))') cas,feedback%status,appended%status,selected%status,gr%status
        if (feedback%status == 0 .and. appended%status == 0 .and. &
            selected%status == 0 .and. gr%status == 0) then
            call output_matrix(feedback%a)
            call output_matrix(feedback%b)
            call output_matrix(feedback%c)
            call output_matrix(feedback%d)
            call output_matrix(appended%d)
            call output_matrix(selected%d)
            ctr = feedback%controllability()
            obs = feedback%observability()
            call output_matrix(ctr)
            call output_matrix(obs)
            call output_matrix(gr%wc)
            call output_matrix(gr%wo)
            write(*,'(es25.16e3,1x,i0)') gr%h2_norm,gr%iterations
        end if
        deallocate(ap,bp,cp,dp1,ak,bk,ck,dk)
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
    subroutine output_matrix(a)
        complex(dp), intent(in) :: a(:,:) !! matrix in column-major order
        integer :: row,column
        do column=1,size(a,2)
            do row=1,size(a,1)
                write(*,'(2(es25.16e3,1x))') real(a(row,column),dp),aimag(a(row,column))
            end do
        end do
    end subroutine output_matrix
end program dlti_interconnections_driver
