! SPDX-License-Identifier: MIT
program test_linalg_banded_sturm
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape, &
        linalg_status_no_convergence
    use scifort_linalg_banded_eigen, only : eigvals_banded, eig_banded
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, &
        ieee_positive_inf, ieee_negative_inf
    implicit none

    real(dp) :: lo(2,2), up(2,2), w(2), short(1), w3(3), zero(3,3)
    complex(dp) :: bc(2,2)
    real(dp), allocatable :: large(:, :)
    real(dp) :: small(3,3), svalues(3), trial(5)
    integer :: status, nfound, j, offset, iblock, n, workcount

    lo = 0.0_dp
    lo(1,:) = 2.0_dp
    lo(2,1) = 1.0_dp
    call eigvals_banded(lo,w,status,lower=.true.,method='sturm')
    call expect(status==linalg_status_success,'real lower, full spectrum')
    call expect(maxval(abs(w-[1.0_dp,3.0_dp]))<2.0e-12_dp,'real eigenvalues')

    up = 0.0_dp
    up(2,:) = 2.0_dp
    up(1,2) = 1.0_dp
    call eigvals_banded(up,w,status,subset_by_index=[2,2],num_selected=nfound)
    call expect(status==linalg_status_success .and. nfound==1,'upper index selection')
    call expect(abs(w(1)-3.0_dp)<2.0e-12_dp,'upper index eigenvalue')
    call expect(w(2)==0.0_dp,'unused output zeroed')

    call eigvals_banded(up,w,status,subset_by_value=[1.0_dp+1.0e-6_dp,4.0_dp], &
        num_selected=nfound)
    call expect(status==linalg_status_success .and. nfound==1,'half-open value selection')
    call expect(abs(w(1)-3.0_dp)<2.0e-12_dp,'value selection result')
    call eigvals_banded(up,w,status,subset_by_value=[ &
        ieee_value(0.0_dp,ieee_negative_inf),ieee_value(0.0_dp,ieee_positive_inf)], &
        num_selected=nfound)
    call expect(status==linalg_status_success .and. nfound==2,'infinite interval')
    call eigvals_banded(up,w,status,subset_by_value=[4.0_dp,5.0_dp],num_selected=nfound)
    call expect(status==linalg_status_success .and. nfound==0,'empty selection')
    call expect(all(w==0.0_dp),'empty-selection output')

    bc = cmplx(0.0_dp,0.0_dp,dp)
    bc(1,:) = cmplx(2.0_dp,0.0_dp,dp)
    bc(2,1) = cmplx(0.0_dp,1.0_dp,dp)
    call eigvals_banded(bc,w,status,lower=.true.,method='sturm')
    call expect(status==linalg_status_success,'complex Hermitian full')
    call expect(maxval(abs(w-[1.0_dp,3.0_dp]))<2.0e-12_dp,'complex Hermitian eigenvalues')

    call eigvals_banded(lo,short,status,lower=.true.,subset_by_index=[1,2],num_selected=nfound)
    call expect(status==linalg_status_invalid_shape .and. nfound==2,'small output capacity')
    call expect(short(1)==0.0_dp,'small buffer zeroed')
    call eigvals_banded(lo,w,status,subset_by_index=[1,1],subset_by_value=[0.0_dp,4.0_dp])
    call expect(status==linalg_status_invalid_shape,'conflicting selection options')
    call eigvals_banded(lo,w,status,subset_by_index=[0,1])
    call expect(status==linalg_status_invalid_shape,'invalid index')
    call eigvals_banded(lo,w,status,subset_by_value=[2.0_dp,2.0_dp])
    call expect(status==linalg_status_invalid_shape,'empty interval bounds invalid')
    call eigvals_banded(lo,w,status,method='unknown')
    call expect(status==linalg_status_invalid_shape,'invalid method')
    call eigvals_banded(lo,w,status,method='sturm',tol=-1.0_dp)
    call expect(status==linalg_status_invalid_shape,'invalid tolerance')
    call eigvals_banded(lo,w,status,method='sturm',max_iterations=0)
    call expect(status==linalg_status_invalid_shape,'invalid iteration limit')
    up(1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    call eigvals_banded(up,w,status,subset_by_index=[1,2])
    call expect(status==linalg_status_success,'unused padding ignored')
    up(1,1) = 0.0_dp
    call eigvals_banded(up,w,status,method='dense',subset_by_index=[1,2])
    call expect(status==linalg_status_success .and. maxval(abs(w-[1.0_dp,3.0_dp]))<1e-12_dp, &
        'dense selected-method compatibility')

    zero = 0.0_dp
    zero(1,1) = -1.0e150_dp
    zero(1,2) = 2.0e-150_dp
    zero(1,3) = 1.0e150_dp
    call eigvals_banded(zero,w3,status,lower=.true.,method='sturm')
    call expect(status==linalg_status_success,'extreme diagonal scales')
    call expect(all(w3==[-1.0e150_dp,2.0e-150_dp,1.0e150_dp]),'diagonal eigenvalues exact')
    zero(2,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    call eigvals_banded(zero,w3,status,lower=.true.,method='sturm')
    call expect(status==linalg_status_no_convergence,'nonfinite stored coefficient')

    ! An order-10000, bandwidth-2 block-diagonal Hermitian problem.
    ! The 3x3 block's eigenvalues are independently obtained with the
    ! previously-tested full dense band solver, then repeated in each block.
    small = 0.0_dp
    small(1,:) = [2.0_dp,3.0_dp,4.0_dp]
    small(2,1:2) = [0.5_dp,-0.75_dp]
    small(3,1) = 0.25_dp
    call eig_banded(small,svalues,status,lower=.true.)
    call expect(status==linalg_status_success,'small block reference')
    n = 10000
    allocate(large(3,n))
    large = 0.0_dp
    do iblock = 1, n/3
        offset = 3*(iblock-1)
        large(1,offset+1:offset+3) = small(1,1:3)
        large(2,offset+1:offset+2) = small(2,1:2)
        large(3,offset+1) = small(3,1)
    end do
    large(1,n) = 10.0_dp
    call eigvals_banded(large,trial,status,lower=.true., &
        subset_by_index=[5000,5004],num_selected=nfound)
    call expect(status==linalg_status_success .and. nfound==5,'order-10000 selected spectrum')
    call expect(maxval(abs(trial-svalues(2)))<2.0e-10_dp,'order-10000 repeated eigenvalues')

    write(*,'(A)') 'banded Sturm native tests passed (including n=10000,kd=2)'
contains
    subroutine expect(test,message)
        logical, intent(in) :: test
        character(len=*), intent(in) :: message
        if (.not. test) then
            print '(A)', 'FAILED: '//message
            error stop 1
        end if
    end subroutine expect
end program test_linalg_banded_sturm
