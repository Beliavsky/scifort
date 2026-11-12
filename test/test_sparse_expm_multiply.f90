! SPDX-License-Identifier: MIT
! Native regression for the action of a sparse real matrix exponential.
program test_sparse_expm_multiply
    use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix, csr_from_coo, csr_from_dense
    use scifort_sparse_linalg, only : expm_multiply, expm_multiply_interval, &
        sparse_linalg_invalid_input, sparse_linalg_breakdown
    implicit none
    type(csr_matrix) :: a
    real(dp), allocatable :: y(:), ym(:,:), yg(:,:), ymg(:,:,:)
    real(dp) :: d(2,2), rhs(2,3), zero(0,0), nan_value, t, pi, big(3,3)
    integer :: info, i, nsteps, nprod

    pi = acos(-1.0_dp)
    d = reshape([0.0_dp,1.0_dp,-1.0_dp,0.0_dp],[2,2])
    a = csr_from_dense(d)
    call expm_multiply(a,[1.0_dp,0.0_dp],y,info,t=pi/2.0_dp, &
        scaling_steps=nsteps,products=nprod)
    call require(info == 0 .and. nsteps >= 1 .and. nprod > 0,'rotation status')
    call close(y,[0.0_dp,1.0_dp],1.0e-13_dp,'rotation')
    rhs = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp,2.0_dp,-3.0_dp],[2,3])
    call expm_multiply(a,rhs,ym,info,t=pi/2.0_dp)
    call require(info == 0 .and. all(shape(ym) == [2,3]), 'block shape')
    call close(ym(:,1),[0.0_dp,1.0_dp],1.0e-13_dp,'block first')
    call close(ym(:,2),[-1.0_dp,0.0_dp],1.0e-13_dp,'block second')
    call close(ym(:,3),[3.0_dp,2.0_dp],1.0e-13_dp,'block third')
    call expm_multiply(a,[1.0_dp,0.0_dp],y,info,t=0.0_dp,products=nprod)
    call require(info == 0 .and. nprod == 0, 'time zero no products')
    call close(y,[1.0_dp,0.0_dp],0.0_dp,'identity action')

    ! Normal and reverse grids, excluding the endpoint, and single-sample grid.
    call expm_multiply_interval(a,[1.0_dp,0.0_dp],-1.0_dp,1.0_dp,5,yg,info)
    call require(info == 0 .and. all(shape(yg) == [2,5]), 'interval shape')
    do i = 1,5
        t = -1.0_dp+real(i-1,dp)/2.0_dp
        call close(yg(:,i),[cos(t),sin(t)],2.0e-13_dp,'grid rotation')
    end do
    call expm_multiply_interval(a,[1.0_dp,0.0_dp],1.0_dp,-1.0_dp,4,yg,info,endpoint=.false.)
    do i = 1,4
        t = 1.0_dp-real(i-1,dp)/2.0_dp
        call close(yg(:,i),[cos(t),sin(t)],2.0e-13_dp,'reverse grid')
    end do
    call expm_multiply_interval(a,rhs,2.0_dp,4.0_dp,1,ymg,info)
    call require(info == 0 .and. all(shape(ymg) == [2,3,1]),'one-time block shape')
    call close(ymg(:,1,1),[cos(2.0_dp),sin(2.0_dp)],1.0e-12_dp,'single grid point')

    ! For an upper triangular nilpotent matrix the series ends exactly.
    big = 0.0_dp
    big(1,2) = 20.0_dp
    big(2,3) = 20.0_dp
    a = csr_from_dense(big)
    call expm_multiply(a,[0.0_dp,0.0_dp,1.0_dp],y,info,t=1.0_dp)
    call require(info == 0, 'nilpotent status')
    call close(y,[200.0_dp,20.0_dp,1.0_dp],1.0e-12_dp,'nilpotent result')

    ! The diagonal trace shift should avoid excessive Taylor scaling.
    a = csr_from_dense(reshape([2.0_dp,0.0_dp,0.0_dp,2.0_dp],[2,2]))
    call expm_multiply(a,[1.0_dp,-1.0_dp],y,info,scaling_steps=nsteps,products=nprod)
    call require(info == 0 .and. nprod == 0, 'pure shift products')
    call close(y,exp(2.0_dp)*[1.0_dp,-1.0_dp],1.0e-13_dp,'diagonal exponential')
    call expm_multiply(a,[1.0_dp,-1.0_dp],y,info,tracea=0.0_dp)
    call require(info == 0, 'supplied trace status')
    call close(y,exp(2.0_dp)*[1.0_dp,-1.0_dp],1.0e-13_dp,'alternate trace')
    call expm_multiply(a,1.0e-300_dp*[1.0_dp,-1.0_dp],y,info)
    call require(info == 0, 'tiny b success')
    call require(maxval(abs(y/1.0e-300_dp-exp(2.0_dp)*[1.0_dp,-1.0_dp])) < 1.0e-12_dp, &
        'tiny b relative accuracy')

    ! The first grid time can underflow even though subsequent times are
    ! representable: anchor at the time closest to zero, not the first time.
    a = csr_from_dense(reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2]))
    call expm_multiply_interval(a,[1.0_dp,2.0_dp],-1000.0_dp,0.0_dp,2,yg,info)
    call require(info == 0, 'underflow interval status')
    call close(yg(:,2),[1.0_dp,2.0_dp],1.0e-14_dp,'underflow interval recovery')
    call close(yg(:,1),[0.0_dp,0.0_dp],0.0_dp,'underflow remote point')

    ! A large sparse example must not densify the operator.
    block
        integer, parameter :: n = 18000
        integer, allocatable :: rows(:), cols(:)
        real(dp), allocatable :: values(:), vec(:)
        allocate(rows(n),cols(n),values(n),vec(n))
        do i = 1,n
            rows(i) = i
            cols(i) = i
        end do
        values = 0.2_dp
        a = csr_from_coo(n,n,rows,cols,values)
        vec = 0.0_dp
        vec(1) = 1.0_dp
        vec(n) = -2.0_dp
        call expm_multiply(a,vec,y,info,products=nprod)
        call require(info == 0, 'large sparse status')
        call require(nprod == 0, 'large constant diagonal trace optimization')
        call require(size(y) == n, 'large sparse result shape')
        call require(abs(y(1)-exp(0.2_dp)) < 1.0e-12_dp, 'large first entry')
        call require(abs(y(n)+2.0_dp*exp(0.2_dp)) < 1.0e-12_dp, 'large last entry')
    end block

    a = csr_from_dense(zero)
    call expm_multiply(a,[real(dp)::],y,info)
    call require(info == 0 .and. size(y) == 0, 'empty zero-by-zero')
    call expm_multiply_interval(a,[real(dp)::],0.0_dp,1.0_dp,3,yg,info)
    call require(info == 0 .and. all(shape(yg) == [0,3]), 'empty interval')

    a = csr_from_dense(reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp],[2,3]))
    call expm_multiply(a,[1.0_dp,2.0_dp],y,info)
    call require(info == sparse_linalg_invalid_input, 'rectangular operator rejected')
    a = csr_from_dense(d)
    nan_value = ieee_value(0.0_dp,ieee_quiet_nan)
    call expm_multiply(a,[1.0_dp,nan_value],y,info)
    call require(info == sparse_linalg_invalid_input, 'nonfinite RHS rejected')
    call expm_multiply(a,[1.0_dp,0.0_dp],y,info,t=nan_value)
    call require(info == sparse_linalg_invalid_input, 'nonfinite time rejected')
    call expm_multiply(a,[1.0_dp,0.0_dp],y,info,tol=-1.0_dp)
    call require(info == sparse_linalg_invalid_input, 'negative tolerance rejected')
    call expm_multiply(a,[1.0_dp,0.0_dp],y,info,tracea=nan_value)
    call require(info == sparse_linalg_invalid_input, 'invalid trace rejected')
    call expm_multiply_interval(a,[1.0_dp,0.0_dp],0.0_dp,1.0_dp,0,yg,info)
    call require(info == sparse_linalg_invalid_input, 'invalid grid length')
    call expm_multiply_interval(a,[1.0_dp,0.0_dp],nan_value,1.0_dp,3,yg,info)
    call require(info == sparse_linalg_invalid_input, 'invalid start time')
    a = csr_from_dense(reshape([1000.0_dp,0.0_dp,0.0_dp,1000.0_dp],[2,2]))
    call expm_multiply(a,[1.0_dp,0.0_dp],y,info)
    call require(info == sparse_linalg_breakdown, 'overflow reported')

    print '(a)', 'test_sparse_expm_multiply: PASS'
contains
    subroutine require(condition, message)
        logical, intent(in) :: condition !! Test assertion.
        character(len=*), intent(in) :: message !! Failure message.
        if (.not. condition) then
            print '(2a)', 'FAIL: ', message
            error stop 1
        end if
    end subroutine require

    subroutine close(actual, expected, tolerance, message)
        real(dp), intent(in) :: actual(:) !! Actual values.
        real(dp), intent(in) :: expected(:) !! Reference values.
        real(dp), intent(in) :: tolerance !! Maximum absolute deviation.
        character(len=*), intent(in) :: message !! Failure message.
        call require(size(actual) == size(expected), message//' size')
        call require(all(abs(actual-expected) <= tolerance), message)
    end subroutine close
end program test_sparse_expm_multiply
