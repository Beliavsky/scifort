! SPDX-License-Identifier: MIT
program test_sparse_constructors
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use scifort_kinds, only : dp
    use scifort_random, only : rng_state, rng_seed, rng_get_state
    use scifort_sparse, only : csr_matrix, csr_to_dense, valid_csr, sparse_success, sparse_invalid_input
    use scifort_sparse_constructors, only : sparse_diagonal, csr_diags_array, csr_random, csr_random_array, csr_rand
    implicit none

    type(csr_matrix) :: a, b
    type(rng_state) :: rng, rng2
    type(sparse_diagonal) :: unset(1)
    real(dp), allocatable :: dense(:,:)
    real(dp) :: nan, mean, variance
    integer(int64) :: state_before(6), state_after(6)
    integer :: status, i

    a = csr_diags_array([sparse_diagonal([1.0_dp,2.0_dp,3.0_dp,4.0_dp]), &
                         sparse_diagonal([5.0_dp,6.0_dp,7.0_dp]), &
                         sparse_diagonal([8.0_dp,9.0_dp])], [0,-1,2],status=status)
    call require(status == sparse_success .and. valid_csr(a),'diags inferred shape')
    call require(a%n_rows == 4 .and. a%n_cols == 4,'diags inferred dimensions')
    dense = csr_to_dense(a)
    call require(all(dense(1,:) == [1.0_dp,0.0_dp,8.0_dp,0.0_dp]),'diags first row')
    call require(all(dense(3,:) == [0.0_dp,6.0_dp,3.0_dp,0.0_dp]),'diags middle row')
    call require(all(dense(4,:) == [0.0_dp,0.0_dp,7.0_dp,4.0_dp]),'diags last row')

    a = csr_diags_array([sparse_diagonal([1.0_dp]),sparse_diagonal([-2.0_dp]), &
                         sparse_diagonal([1.0_dp])], [-1,0,1],8,8,status)
    call require(status == sparse_success .and. size(a%data) == 22,'scalar broadcast')
    dense = csr_to_dense(a)
    do i = 1, 8
        call require(dense(i,i) == -2.0_dp,'scalar main diagonal')
    end do
    call require(dense(1,2) == 1.0_dp .and. dense(8,7) == 1.0_dp,'scalar off diagonal')

    a = csr_diags_array([sparse_diagonal([1.0_dp,2.0_dp,3.0_dp,4.0_dp])], &
                        [1],2,3,status)
    call require(status == sparse_success .and. size(a%data) == 2,'truncate long diagonal')
    dense = csr_to_dense(a)
    call require(dense(1,2) == 1.0_dp .and. dense(2,3) == 2.0_dp,'truncated content')
    a = csr_diags_array([sparse_diagonal([1.0_dp,2.0_dp,3.0_dp])],[2],status=status)
    call require(status == sparse_success .and. a%n_rows == 5 .and. a%n_cols == 5, 'inferred offset')
    a = csr_diags_array([sparse_diagonal([0.0_dp,2.0_dp])],[0],2,2,status)
    call require(status == sparse_success .and. size(a%data) == 2,'explicit diagonal zero')
    a = csr_diags_array([sparse_diagonal([1.0_dp])],[0],0,6,status)
    call require(status == sparse_success .and. valid_csr(a),'empty row dimensions')
    a = csr_diags_array([sparse_diagonal([1.0_dp])],[0],5,0,status)
    call require(status == sparse_success .and. valid_csr(a),'empty column dimensions')

    a = csr_diags_array([sparse_diagonal([1.0_dp])],[0],4,4,status)
    call require(status == sparse_success,'scalar diagonal allowed with shape')
    a = csr_diags_array([sparse_diagonal([1.0_dp,2.0_dp])],[0],4,4,status)
    call require(status == sparse_invalid_input,'short diagonal rejected')
    a = csr_diags_array([sparse_diagonal([1.0_dp]),sparse_diagonal([2.0_dp])], &
                        [0,0],4,4,status)
    call require(status == sparse_invalid_input,'duplicate offsets rejected')
    a = csr_diags_array([sparse_diagonal([1.0_dp])],[0],4,status=status)
    call require(status == sparse_invalid_input,'partial explicit shape rejected')
    a = csr_diags_array(unset,[0],2,2,status)
    call require(status == sparse_invalid_input,'unallocated diagonal rejected')
    nan = ieee_value(0.0_dp,ieee_quiet_nan)
    a = csr_diags_array([sparse_diagonal([nan])],[0],2,2,status)
    call require(status == sparse_invalid_input,'nonfinite diagonal rejected')
    a = csr_diags_array([sparse_diagonal([1.0_dp])],[8],2,2,status)
    call require(status == sparse_invalid_input,'offset outside bounds rejected')

    call rng_seed(rng,77)
    rng2 = rng
    a = csr_random(13,17,0.2_dp,rng,status)
    b = csr_random(13,17,0.2_dp,rng2,status)
    call require(status == sparse_success .and. valid_csr(a),'random CSR canonical')
    call require(size(a%data) == 44,'random exact nonzero count')
    call require(all(a%indptr == b%indptr),'random reproducible pointers')
    call require(all(a%indices == b%indices) .and. all(a%data == b%data),'random reproducible values')
    call require(all(a%data > 0.0_dp) .and. all(a%data < 1.0_dp),'uniform range')
    call rng_get_state(rng,state_before)
    a = csr_random(13,17,-0.1_dp,rng,status)
    call rng_get_state(rng,state_after)
    call require(status == sparse_invalid_input .and. all(state_before == state_after),'invalid density keeps state')
    a = csr_random(13,17,0.2_dp,rng,status,distribution='unknown')
    call rng_get_state(rng,state_after)
    call require(status == sparse_invalid_input .and. all(state_before == state_after),'invalid distribution keeps state')
    a = csr_random(13,17,0.2_dp,rng,status,data_values=[0.0_dp])
    call rng_get_state(rng,state_after)
    call require(status == sparse_invalid_input .and. all(state_before == state_after),'invalid values size keeps state')

    call rng_seed(rng,192)
    a = csr_random(3,5,1.0_dp,rng,status)
    call require(status == sparse_success .and. size(a%data) == 15,'density one')
    a = csr_random_array(3,5,0.0_dp,rng,status)
    call require(status == sparse_success .and. size(a%data) == 0,'density zero')
    a = csr_random(0,5,0.9_dp,rng,status)
    call require(status == sparse_success .and. valid_csr(a),'zero-row random')
    a = csr_random(5,0,0.9_dp,rng,status)
    call require(status == sparse_success .and. valid_csr(a),'zero-column random')
    a = csr_random(2,2,0.125_dp,rng,status)
    call require(status == sparse_success .and. size(a%data) == 0,'bankers round 0.5')
    a = csr_random(2,2,0.375_dp,rng,status)
    call require(status == sparse_success .and. size(a%data) == 2,'bankers round 1.5')
    a = csr_random(2,2,0.625_dp,rng,status)
    call require(status == sparse_success .and. size(a%data) == 2,'bankers round 2.5')

    a = csr_random(3,1000000000,3.0e-9_dp,rng,status)
    call require(status == sparse_success .and. valid_csr(a),'billion-column random CSR')
    call require(size(a%data) == 9,'billion-column exact count')
    call require(all(a%indices >= 1) .and. all(a%indices <= 1000000000),'large coordinates')

    a = csr_random(2,3,0.5_dp,rng,status,data_values=[-2.0_dp,0.0_dp,8.0_dp])
    call require(status == sparse_success .and. size(a%data) == 3,'custom values')
    call require(any(a%data == 0.0_dp) .and. any(a%data == -2.0_dp),'explicit custom zero')
    a = csr_rand(400,400,0.2_dp,rng,status)
    call require(status == sparse_success .and. size(a%data) == 32000,'rand alias')
    mean = sum(a%data)/real(size(a%data),dp)
    call require(abs(mean-0.5_dp) < 0.02_dp,'uniform mean sanity')

    a = csr_random(400,400,0.2_dp,rng,status,'normal')
    call require(status == sparse_success .and. size(a%data) == 32000,'normal values count')
    mean = sum(a%data)/real(size(a%data),dp)
    variance = sum((a%data-mean)**2)/real(size(a%data),dp)
    call require(abs(mean) < 0.04_dp,'normal mean sanity')
    call require(abs(variance-1.0_dp) < 0.06_dp,'normal variance sanity')

    a = csr_random(-1,12,0.1_dp,rng,status)
    call require(status == sparse_invalid_input,'invalid shape')
    a = csr_random(100000,100000,1.0_dp,rng,status)
    call require(status == sparse_invalid_input,'unrepresentable full CSR')
    a = csr_random(3,3,nan,rng,status)
    call require(status == sparse_invalid_input,'nonfinite density')
    print '(a)', 'test_sparse_constructors: PASS'
contains
    subroutine require(condition,message)
        logical, intent(in) :: condition !! Assertion condition.
        character(len=*), intent(in) :: message !! Error description.
        if (.not. condition) then
            print '(a)', 'FAILED: '//trim(message)
            error stop 1
        end if
    end subroutine require
end program test_sparse_constructors
