! SPDX-License-Identifier: MIT
! Ensure arbitrary-rank analytic transforms are available through scifort_signal.
program test_signal_hilbert_nd_public
    use scifort_signal, only : dp, hilbert, hilbert2, hilbert_into, hilbert2_into, &
        hilbert_status_success
    implicit none
    real(dp) :: x(6,5,4)
    complex(dp) :: y(6,5,7), z(8,5,7)
    complex(dp), allocatable :: ref(:), ref2(:,:)
    integer :: j, k, st

    do k=1,4
        do j=1,5
            x(:,j,k) = cos(real([1,2,3,4,5,6]*j+3*k,dp))
        end do
    end do
    call hilbert_into(x,y,axis=3,n=7,status=st)
    if (st /= hilbert_status_success) error stop 'public Hilbert status'
    do j=1,5
        do k=1,6
            ref = hilbert(x(k,j,:),n=7)
            if (maxval(abs(y(k,j,:)-ref)) > 1.0e-11_dp) &
                error stop 'public Hilbert values'
        end do
    end do
    call hilbert2_into(x,z,axes=[1,3],n=[8,7],status=st)
    if (st /= hilbert_status_success) error stop 'public Hilbert2 status'
    do j=1,5
        ref2 = hilbert2(x(:,j,:),n=[8,7])
        if (maxval(abs(z(:,j,:)-ref2)) > 1.0e-11_dp) &
            error stop 'public Hilbert2 values'
    end do
    print '(a)', 'test_signal_hilbert_nd_public passed'
end program test_signal_hilbert_nd_public
