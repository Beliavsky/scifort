! SPDX-License-Identifier: MIT
! Read complex state-space systems and report reduced realizations.
program dlti_minreal_driver
    use scifort_kinds, only : dp
    use scifort_signal_dlti, only : dlti_system, dlti_ss
    implicit none
    type(dlti_system) :: original,reduced
    complex(dp), allocatable :: a(:,:),b(:,:),c(:,:),d(:,:)
    integer :: ncase,cas,n,nu,ny
    real(dp) :: dt
    read(*,*) ncase
    do cas=1,ncase
        read(*,*) n,nu,ny,dt
        allocate(a(n,n),b(n,nu),c(ny,n),d(ny,nu))
        call input_matrix(a)
        call input_matrix(b)
        call input_matrix(c)
        call input_matrix(d)
        original = dlti_ss(a,b,c,d,dt)
        reduced = original%minreal()
        write(*,'(3(i0,1x))') cas,reduced%status,reduced%nstates()
        if (reduced%status == 0) then
            call output_matrix(reduced%a)
            call output_matrix(reduced%b)
            call output_matrix(reduced%c)
            call output_matrix(reduced%d)
        end if
        deallocate(a,b,c,d)
    end do
contains
    subroutine input_matrix(matrix)
        complex(dp), intent(out) :: matrix(:,:) !! column-major input
        integer :: i,j
        do j=1,size(matrix,2)
            do i=1,size(matrix,1)
                read(*,*) matrix(i,j)
            end do
        end do
    end subroutine input_matrix
    subroutine output_matrix(matrix)
        complex(dp), intent(in) :: matrix(:,:) !! column-major output
        integer :: i,j
        do j=1,size(matrix,2)
            do i=1,size(matrix,1)
                write(*,'(2(es25.16e3,1x))') real(matrix(i,j),dp),aimag(matrix(i,j))
            end do
        end do
    end subroutine output_matrix
end program dlti_minreal_driver
