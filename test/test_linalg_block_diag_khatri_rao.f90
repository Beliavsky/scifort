! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_linalg_block_diag_khatri_rao
    use, intrinsic :: iso_fortran_env, only : int32, int64
    use scifort_kinds, only : dp
    use scifort_linalg, only : linalg_status_success, linalg_status_invalid_shape
    use scifort_linalg_block_diag, only : linalg_block, as_block, block_diag
    use scifort_linalg_khatri_rao, only : khatri_rao
    implicit none
    integer :: status, i, j, k, q
    type(linalg_block) :: no_blocks(0)
    real(dp) :: a(2,2), r(4,5), v(2), zero3(3,2)
    real(dp) :: batch1(1,2,2), batch2(3,1,1), bout(3,3,3), badbat(2,2,2)
    real(dp) :: zero_batch(0,2,2), zero_result(0,4,4), bad_shape(3,5)
    real(dp) :: noout(1,0), evec(0), emat(0,3), ebout(1,3)
    complex(dp) :: zc(2,3), zout(4,5)
    integer(int64) :: li(2,2), il(2,2)
    integer(int32) :: i32out(2,2)
    real(dp) :: ka(2,3), kb(3,3), kc(6,3), kcwrong(4,3)
    real(dp) :: kba(2,2,3), kbb(1,3,3), kbc(2,6,3), kbc2(2,6,3)
    real(dp) :: kzero(0,2,3), krzero(0,6,3), kempty(0,3), kemptyout(0,3)
    real(dp) :: knocols(2,0), knocolsb(3,0), knocolsout(6,0)
    real(dp) :: wrongk(3,4)
    complex(dp) :: kcomplex(2,3), kcomplexout(6,3)
    integer(int32) :: ka32(2,2), kb32(2,2), kc32(4,2)
    integer(int64) :: ka64(2,2), kb64(2,2), kc64(4,2)

    a(1,:) = [1.0_dp,2.0_dp]
    a(2,:) = [3.0_dp,4.0_dp]
    v = [6.0_dp,7.0_dp]
    call block_diag([as_block(a), as_block(v), as_block(8.0_dp)], r, status)
    call check(status == linalg_status_success, 'mixed rank real block status')
    call check(all(r(1:2,1:2) == a), 'first real block')
    call check(all(r(3,3:4) == v) .and. r(4,5) == 8.0_dp, 'vector and scalar blocks')
    call check(all(r(1:2,3:5) == 0.0_dp) .and. all(r(3:4,1:2) == 0.0_dp), &
        'off-diagonal zeros')

    li = reshape([9007199254740993_int64, -11_int64, 77_int64, 13_int64], [2,2])
    call block_diag([as_block(li)], il, status)
    call check(status == 0 .and. all(il == li), 'exact int64 block')
    call block_diag([as_block(li)], i32out, status)
    call check(status == linalg_status_invalid_shape .and. all(i32out == 0), &
        'int64 to int32 overflow rejected')

    zc(1,:) = [cmplx(1.0_dp,2.0_dp,dp),cmplx(2.0_dp,-1.0_dp,dp),cmplx(3.0_dp,0.0_dp,dp)]
    zc(2,:) = [cmplx(4.0_dp,0.0_dp,dp),cmplx(5.0_dp,1.0_dp,dp),cmplx(6.0_dp,7.0_dp,dp)]
    call block_diag([as_block(a), as_block(zc)], zout, status)
    call check(status == 0 .and. all(zout(1:2,1:2) == cmplx(a,kind=dp)), &
        'real and complex diagonal blocks')
    call check(all(zout(3:4,3:5) == zc), 'complex block values')
    call block_diag([as_block(zc)], bad_shape, status)
    call check(status == linalg_status_invalid_shape .and. all(bad_shape == 0.0_dp), &
        'complex to real disallowed')

    batch1(1,:,:) = a
    do q=1,3
        batch2(q,1,1) = real(q,dp)*11.0_dp
    end do
    call block_diag([as_block(batch1),as_block(batch2)], bout, status)
    call check(status == 0, 'broadcast blocks')
    do q=1,3
        call check(all(bout(q,1:2,1:2) == a), 'broadcast block contents')
        call check(bout(q,3,3) == batch2(q,1,1), 'batch-specific scalar matrix')
    end do
    badbat = 7.0_dp
    call block_diag([as_block(badbat), as_block(batch2)], bout, status)
    call check(status == linalg_status_invalid_shape .and. all(bout == 0.0_dp), &
        'incompatible batch sizes')
    call block_diag([as_block(batch1)], a, status)
    call check(status == linalg_status_invalid_shape .and. all(a == 0.0_dp), &
        'rank3 input not accepted in rank2 output')
    a(1,:) = [1.0_dp,2.0_dp]
    a(2,:) = [3.0_dp,4.0_dp]

    call block_diag([as_block(zero_batch),as_block(a)], zero_result, status)
    call check(status == linalg_status_success, 'empty batch allowed in block diag')
    call block_diag(no_blocks,noout,status)
    call check(status == 0, 'zero blocks returns a 1-by-0 output')
    call block_diag([as_block(evec),as_block(a)], ebout, status)
    call check(status == linalg_status_invalid_shape .and. all(ebout == 0.0_dp), &
        'empty vector still occupies a row')
    call block_diag([as_block(emat),as_block(v)], zero3, status)
    call check(status == linalg_status_invalid_shape .and. all(zero3 == 0.0_dp), &
        'empty matrix retains row and column count')

    ka(1,:) = [1.0_dp,2.0_dp,3.0_dp]
    ka(2,:) = [4.0_dp,5.0_dp,6.0_dp]
    kb(1,:) = [3.0_dp,4.0_dp,5.0_dp]
    kb(2,:) = [6.0_dp,7.0_dp,8.0_dp]
    kb(3,:) = [2.0_dp,3.0_dp,9.0_dp]
    call khatri_rao(ka,kb,kc,status)
    call check(status == 0, 'Khatri-Rao real status')
    do k=1,3
        do i=1,2
            do j=1,3
                call check(kc((i-1)*3+j,k) == ka(i,k)*kb(j,k), 'Khatri-Rao row order')
            end do
        end do
    end do
    call khatri_rao(ka,kb,kcwrong,status)
    call check(status == linalg_status_invalid_shape .and. all(kcwrong == 0.0_dp), &
        'bad Khatri-Rao output shape')

    kcomplex = cmplx(ka,ka*0.5_dp,dp)
    call khatri_rao(kcomplex,kb,kcomplexout,status)
    call check(status == 0, 'complex-real Khatri-Rao status')
    call check(all(abs(kcomplexout(:,2)-cmplx(kc(:,2),0.5_dp*kc(:,2),dp)) &
        < 1.0e-13_dp), 'complex-real Khatri-Rao values')
    call khatri_rao(ka,kcomplex,kcomplexout,status)
    call check(status == linalg_status_invalid_shape .and. all(kcomplexout == 0.0_dp), &
        'wrong mixed output shape rejected')

    ka32 = reshape([1_int32,2_int32,3_int32,4_int32],[2,2])
    kb32 = reshape([5_int32,6_int32,7_int32,8_int32],[2,2])
    call khatri_rao(ka32,kb32,kc32,status)
    call check(status == 0 .and. kc32(4,2) == 32_int32, 'int32 Khatri-Rao')
    ka64 = int(ka32,int64)
    kb64 = int(kb32,int64)
    call khatri_rao(ka64,kb64,kc64,status)
    call check(status == 0 .and. all(kc64 == int(kc32,int64)), 'int64 Khatri-Rao')

    do q=1,2
        kba(q,:,:) = ka*real(q,dp)
    end do
    kbb(1,:,:) = kb
    call khatri_rao(kba,kbb,kbc,status)
    call check(status == 0 .and. all(kbc(2,:,:) == 2.0_dp*kc), &
        '3D/3D singleton broadcasting')
    call khatri_rao(ka,kba,kbc2,status)
    call check(status == linalg_status_invalid_shape .and. all(kbc2 == 0.0_dp), &
        '2D/3D mismatched dimensions')
    call khatri_rao(ka,kbb,kbc,status)
    call check(status == linalg_status_invalid_shape, '2D/3D output batch mismatch')
    call khatri_rao(kb,kba,kbc,status)
    call check(status == 0, '2D/3D valid product')
    call check(kbc(2,1,1) == kb(1,1)*kba(2,1,1), '2D/3D row order')
    call khatri_rao(kba,kb,kbc,status)
    call check(status == 0, '3D/2D valid product')
    call check(kbc(2,1,1) == kba(2,1,1)*kb(1,1), '3D/2D row order')
    call khatri_rao(kzero,kbb,krzero,status)
    call check(status == linalg_status_invalid_shape, 'zero-sized Khatri-Rao batch rejected')
    call khatri_rao(kempty,kb,kemptyout,status)
    call check(status == 0, 'zero rows accepted')
    call khatri_rao(knocols,knocolsb,knocolsout,status)
    call check(status == 0, 'zero columns accepted')
    wrongk = 5.0_dp
    call khatri_rao(ka,wrongk,kc,status)
    call check(status == linalg_status_invalid_shape .and. all(kc == 0.0_dp), &
        'incompatible column counts rejected')

    print *, 'SciFort block_diag and khatri_rao native tests passed'
contains
    subroutine check(ok, message)
        logical, intent(in) :: ok !! assertion result
        character(*), intent(in) :: message !! failure label
        if (.not. ok) then
            print *, 'FAILED: ', message
            error stop 1
        end if
    end subroutine check
end program test_linalg_block_diag_khatri_rao
