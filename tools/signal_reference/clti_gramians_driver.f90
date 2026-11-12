! SPDX-License-Identifier: MIT
! Test driver for continuous-time Gramians and balanced model reduction.
program clti_gramians_driver
    use scifort_kinds, only : dp
    use scifort_signal_clti, only : clti_system,clti_ss,clti_gramian_result, &
        clti_balanced_result
    use scifort_signal_lti_frequency, only : mimo_frequency_response
    implicit none
    type(clti_system) :: original
    type(clti_gramian_result) :: gr
    type(clti_balanced_result) :: red
    type(mimo_frequency_response) :: freq
    complex(dp), allocatable :: a(:,:),b(:,:),c(:,:),d(:,:)
    real(dp) :: w(8)
    integer :: case_count,cas,n,nu,ny,order,j,k,i

    w = [0.0_dp,0.10_dp,0.30_dp,0.60_dp,1.0_dp,1.6_dp,2.5_dp,4.0_dp]
    read(*,*) case_count
    do cas=1,case_count
        read(*,*) n,nu,ny,order
        allocate(a(n,n),b(n,nu),c(ny,n),d(ny,nu))
        call input_matrix(a)
        call input_matrix(b)
        call input_matrix(c)
        call input_matrix(d)
        original = clti_ss(a,b,c,d)
        gr = original%gramians()
        red = original%balanced_truncate(order)
        write(*,'(3(i0,1x))') cas,gr%status,red%status
        if (gr%status == 0 .and. red%status == 0) then
            call output_matrix(gr%wc)
            call output_matrix(gr%wo)
            write(*,'(es26.17e3)') gr%h2_norm
            do i=1,size(red%hsv)
                write(*,'(es26.17e3)') red%hsv(i)
            end do
            write(*,'(es26.17e3)') red%error_bound
            freq = red%model%freqresp_mimo(w=w,solver='direct')
            if (freq%status /= 0) error stop 'continuous frequency response failed'
            do k=1,size(w)
                call output_matrix(freq%response(:,:,k))
            end do
        end if
        deallocate(a,b,c,d)
    end do
contains
    subroutine input_matrix(mat)
        complex(dp), intent(out) :: mat(:,:) !! matrix to read by columns
        integer :: row,col
        do col=1,size(mat,2)
            do row=1,size(mat,1)
                read(*,*) mat(row,col)
            end do
        end do
    end subroutine input_matrix

    subroutine output_matrix(mat)
        complex(dp), intent(in) :: mat(:,:) !! matrix to serialize by columns
        integer :: row,col
        do col=1,size(mat,2)
            do row=1,size(mat,1)
                write(*,'(2(es26.17e3,1x))') real(mat(row,col),dp),aimag(mat(row,col))
            end do
        end do
    end subroutine output_matrix
end program clti_gramians_driver
