! SPDX-License-Identifier: MIT
program test_linalg_companion
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    use scifort_linalg_companion, only : companion, fiedler_companion
    use scifort_linalg_fiedler, only : fiedler
    implicit none
    integer :: s, n, i, j, k, kind_case
    real(dp) :: a(8), c(7,7), fc(7,7), f(8,8), expected(4,4)
    complex(dp) :: z(8), zc(7,7), zfc(7,7)

    ! SciPy's published companion example.
    block
        real(dp) :: m(3,3)
        call companion([1.0_dp,-10.0_dp,31.0_dp,-30.0_dp],m,s)
        call assert_success(s,'real companion example')
        if (any(m(1,:) /= [10.0_dp,-31.0_dp,30.0_dp]) .or. &
            any(m(2,:) /= [1.0_dp,0.0_dp,0.0_dp]) .or. &
            any(m(3,:) /= [0.0_dp,1.0_dp,0.0_dp])) error stop 'companion example'
    end block

    ! SciPy's published Fiedler companion example, degree four.
    block
        real(dp) :: m(4,4)
        expected(1,:) = [16.0_dp,-86.0_dp,1.0_dp,0.0_dp]
        expected(2,:) = [1.0_dp,0.0_dp,0.0_dp,0.0_dp]
        expected(3,:) = [0.0_dp,176.0_dp,0.0_dp,-105.0_dp]
        expected(4,:) = [0.0_dp,1.0_dp,0.0_dp,0.0_dp]
        call fiedler_companion([1.0_dp,-16.0_dp,86.0_dp,-176.0_dp,105.0_dp],m,s)
        call assert_success(s,'Fiedler companion example')
        if (any(m /= expected)) error stop 'Fiedler companion example mismatch'
    end block

    ! A Fiedler distance matrix, including integer and complex overloads.
    block
        integer :: fi(5,5)
        real(dp) :: fr(2,2)
        call fiedler([1,4,12,45,77],fi,s)
        call assert_success(s,'integer Fiedler example')
        if (any(fi(:,1) /= [0,3,11,44,76]) .or. &
            any(fi(2,:) /= [3,0,8,41,73])) error stop 'integer Fiedler example'
        call fiedler([cmplx(0.0_dp,1.0_dp,dp),cmplx(3.0_dp,5.0_dp,dp)],fr,s)
        call assert_success(s,'complex Fiedler')
        if (abs(fr(1,2)-5.0_dp) > 1.0e-14_dp .or. &
            abs(fr(2,1)-5.0_dp) > 1.0e-14_dp) error stop 'complex Fiedler modulus'
    end block

    ! Exercise all degrees, parity branches, coefficient scaling and dtype.
    do kind_case=1,5
        do n=2,8
            do i=1,n
                select case (kind_case)
                case (1)
                    a(i) = real(i,dp)
                case (2)
                    a(i) = real(i*i-2*i-6,dp)
                case (3)
                    a(i) = (-1.0_dp)**i * (1.0_dp+real(i,dp)/10.0_dp)
                case (4)
                    a(i) = (real(i,dp)/10.0_dp)**3
                case (5)
                    a(i) = 1.0e-25_dp*real(i,dp)
                end select
                z(i) = cmplx(a(i),a(i)*0.3_dp-0.14_dp,dp)
            end do
            call companion(a(:n),c(:n-1,:n-1),s)
            call assert_success(s,'real companion loop')
            call fiedler_companion(a(:n),fc(:n-1,:n-1),s)
            call assert_success(s,'real Fiedler companion loop')
            call companion(z(:n),zc(:n-1,:n-1),s)
            call assert_success(s,'complex companion loop')
            call fiedler_companion(z(:n),zfc(:n-1,:n-1),s)
            call assert_success(s,'complex Fiedler companion loop')
            call fiedler(a(:n),f(:n,:n),s)
            call assert_success(s,'real Fiedler loop')
            if (maxval(abs(c(1,:n-1)+a(2:n)/a(1))) > 1.0e-13_dp) &
                error stop 'companion row relation'
            do i=2,n-1
                if (c(i,i-1) /= 1.0_dp) error stop 'companion subdiagonal'
            end do
            if (maxval(abs(f(:n,:n)-transpose(f(:n,:n)))) > 1.0e-13_dp) &
                error stop 'Fiedler symmetry'
            do i=1,n
                do j=1,n
                    if (abs(f(i,j)-abs(a(i)-a(j))) > 1.0e-13_dp) &
                        error stop 'Fiedler distance'
                end do
            end do
            if (abs(c(1,1)-fc(1,1)) > 1.0e-13_dp) &
                error stop 'companion trace entry'
            if (abs(zc(1,1)-zfc(1,1)) > 1.0e-13_dp) &
                error stop 'complex companion trace entry'
        end do
    end do

    ! Degree zero and empty Fiedler matrices.
    block
        real(dp) :: one(1,1), zero(0,0), not_square(2,1), bad_comp(3,3)
        integer :: izero(0,0), i_fiedler(3,3)
        call companion([2.0_dp,6.0_dp],one,s)
        call assert_success(s,'degree one companion')
        if (one(1,1) /= -3.0_dp) error stop 'degree one companion value'
        call fiedler_companion([2.0_dp,6.0_dp],one,s)
        call assert_success(s,'degree one Fiedler companion')
        if (one(1,1) /= -3.0_dp) error stop 'degree one Fiedler value'
        call fiedler_companion([0.0_dp],zero,s)
        call assert_success(s,'degree zero Fiedler companion')
        call fiedler_companion([integer::],zero,s)
        call assert_success(s,'empty Fiedler companion')
        call fiedler([integer::],izero,s)
        call assert_success(s,'empty Fiedler')
        call companion([1.0_dp],zero,s)
        if (s /= linalg_status_invalid_shape) error stop 'companion accepted degree zero'
        call companion([0.0_dp,2.0_dp],one,s)
        if (s /= linalg_status_invalid_shape .or. any(one /= 0.0_dp)) &
            error stop 'companion accepted leading zero'
        call fiedler_companion([0.0_dp,2.0_dp],one,s)
        if (s /= linalg_status_invalid_shape .or. any(one /= 0.0_dp)) &
            error stop 'Fiedler companion accepted leading zero'
        call companion([1.0_dp,2.0_dp,3.0_dp],bad_comp,s)
        if (s /= linalg_status_invalid_shape .or. any(bad_comp /= 0.0_dp)) &
            error stop 'companion wrong output shape'
        call fiedler_companion([1.0_dp,2.0_dp,3.0_dp],bad_comp,s)
        if (s /= linalg_status_invalid_shape .or. any(bad_comp /= 0.0_dp)) &
            error stop 'Fiedler companion wrong output shape'
        call fiedler([1.0_dp,2.0_dp],not_square,s)
        if (s /= linalg_status_invalid_shape .or. any(not_square /= 0.0_dp)) &
            error stop 'Fiedler wrong output shape'
        call fiedler([1,4,7],i_fiedler,s)
        call assert_success(s,'integer Fiedler')
        if (any(i_fiedler(:,1) /= [0,3,6])) error stop 'integer Fiedler output'
    end block

    ! Batching must equal independent per-slice operations and fail atomically.
    block
        real(dp) :: batch(3,6), batch_c(3,5,5), one_c(5,5), batch_f(3,6,6)
        real(dp) :: invalid(3,4,4), batch_bad(3,6)
        real(dp) :: int_c(2,2,2), z_fiedler(3,6,6), one_fiedler(6,6)
        complex(dp) :: batch_z(3,6), batch_zc(3,5,5), one_zc(5,5)
        integer :: int_a(2,3), fi(2,3,3)
        do j=1,3
            do i=1,6
                batch(j,i) = 1.0_dp + real(i+j,dp)*0.1_dp
                batch_z(j,i) = cmplx(batch(j,i),real(i-j,dp)*0.3_dp,dp)
            end do
        end do
        call companion(batch,batch_c,s)
        call assert_success(s,'real companion batch')
        do k=1,3
            call companion(batch(k,:),one_c,s)
            if (any(batch_c(k,:,:) /= one_c)) error stop 'real companion batch mismatch'
        end do
        call fiedler_companion(batch,batch_c,s)
        call assert_success(s,'real Fiedler companion batch')
        do k=1,3
            call fiedler_companion(batch(k,:),one_c,s)
            if (any(batch_c(k,:,:) /= one_c)) error stop 'real Fiedler companion batch mismatch'
        end do
        call fiedler(batch,batch_f,s)
        call assert_success(s,'real Fiedler batch')
        call companion(batch_z,batch_zc,s)
        call assert_success(s,'complex companion batch')
        do k=1,3
            call companion(batch_z(k,:),one_zc,s)
            if (any(batch_zc(k,:,:) /= one_zc)) error stop 'complex companion batch mismatch'
        end do
        call fiedler_companion(batch_z,batch_zc,s)
        call assert_success(s,'complex Fiedler companion batch')
        do k=1,3
            call fiedler_companion(batch_z(k,:),one_zc,s)
            if (any(batch_zc(k,:,:) /= one_zc)) error stop 'complex Fiedler companion batch mismatch'
        end do
        call fiedler(batch_z,z_fiedler,s)
        call assert_success(s,'complex Fiedler batch')
        do k=1,3
            call fiedler(batch_z(k,:),one_fiedler,s)
            if (any(z_fiedler(k,:,:) /= one_fiedler)) &
                error stop 'complex Fiedler batch mismatch'
        end do
        batch_bad = batch
        batch_bad(2,1) = 0.0_dp
        call companion(batch_bad,batch_c,s)
        if (s /= linalg_status_invalid_shape .or. any(batch_c /= 0.0_dp)) &
            error stop 'batch companion leading zero handling'
        call fiedler_companion(batch_bad,batch_c,s)
        if (s /= linalg_status_invalid_shape .or. any(batch_c /= 0.0_dp)) &
            error stop 'batch Fiedler leading zero handling'
        call companion(batch,invalid,s)
        if (s /= linalg_status_invalid_shape .or. any(invalid /= 0.0_dp)) &
            error stop 'invalid batch shape'
        call fiedler(batch,invalid,s)
        if (s /= linalg_status_invalid_shape .or. any(invalid /= 0.0_dp)) &
            error stop 'invalid Fiedler batch shape'
        batch_z(2,1) = (0.0_dp,0.0_dp)
        call companion(batch_z,batch_zc,s)
        if (s /= linalg_status_invalid_shape .or. &
            any(batch_zc /= (0.0_dp,0.0_dp))) &
            error stop 'complex batch companion leading zero handling'
        call fiedler_companion(batch_z,batch_zc,s)
        if (s /= linalg_status_invalid_shape .or. &
            any(batch_zc /= (0.0_dp,0.0_dp))) &
            error stop 'complex batch Fiedler leading zero handling'
        int_a(1,:) = [1,-3,2]
        int_a(2,:) = [2,6,4]
        call companion(int_a,int_c,s)
        call assert_success(s,'integer companion batch')
        if (any(int_c(1,1,:) /= [3.0_dp,-2.0_dp])) error stop 'integer companion batch values'
        call fiedler_companion(int_a,int_c,s)
        call assert_success(s,'integer Fiedler batch')
        call fiedler(int_a,fi,s)
        call assert_success(s,'integer distance batch')
        if (any(fi(1,1,:) /= [0,4,1])) error stop 'integer distance values'
    end block

    print '(a)', 'test_linalg_companion: PASS'
contains
    subroutine assert_success(status, tag)
        integer, intent(in) :: status !! actual return code
        character(*), intent(in) :: tag !! name of failing assertion
        if (status /= linalg_status_success) then
            print '(a,1x,a,1x,i0)', 'FAILED:',tag,status
            error stop 'nonzero status'
        end if
    end subroutine assert_success
end program test_linalg_companion
