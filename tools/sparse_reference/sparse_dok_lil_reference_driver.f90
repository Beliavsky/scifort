! SPDX-License-Identifier: MIT
! Driven by compare_sparse_dok_lil_scipy.py for randomized independent tests.
program sparse_dok_lil_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_matrix
    use scifort_sparse_coo, only : coo_matrix, coo_from_triplets
    use scifort_sparse_dok_lil
    implicit none
    type(dok_matrix) :: dok, dok2
    type(lil_matrix) :: lil, lil2
    type(coo_matrix) :: coo
    type(csr_matrix) :: csr_d, csr_l
    integer, allocatable :: rr(:), cc(:)
    real(dp), allocatable :: vv(:), x(:), yd(:), yl(:), dd(:,:), dl(:,:)
    real(dp) :: v, result_d, result_l
    integer :: ncases, icase, m, n, ntrip, nops, k, i, j, info, op

    read(*,*) ncases
    do icase = 1, ncases
        read(*,*) m, n, ntrip, nops
        allocate(rr(ntrip),cc(ntrip),vv(ntrip),x(n))
        do k = 1, ntrip
            read(*,*) rr(k),cc(k),vv(k)
        end do
        if (n > 0) read(*,*) x
        coo = coo_from_triplets(m,n,rr,cc,vv,info)
        call require(info == 0, 'COO case input')
        dok = dok_from_coo(coo,info)
        call require(info == 0, 'DOK from COO')
        lil = lil_from_coo(coo,info)
        call require(info == 0, 'LIL from COO')
        do k = 1, nops
            read(*,*) op, i, j, v
            select case (op)
            case (1)
                call dok_set(dok,i,j,v,info)
                call require(info == 0, 'DOK set')
                call lil_set(lil,i,j,v,info)
                call require(info == 0, 'LIL set')
            case (2)
                call dok_add(dok,i,j,v,info)
                call require(info == 0, 'DOK add')
                call lil_add(lil,i,j,v,info)
                call require(info == 0, 'LIL add')
            case (3)
                call dok_delete(dok,i,j,info)
                call require(info == 0, 'DOK delete')
                call lil_delete(lil,i,j,info)
                call require(info == 0, 'LIL delete')
            case (4)
                ! Read-only lookup, including missing entries.
            case default
                error stop 'unknown reference operation'
            end select
            result_d = dok_get(dok,i,j,info)
            call require(info == 0, 'DOK get')
            result_l = lil_get(lil,i,j,info)
            call require(info == 0, 'LIL get')
            write(*,'(A,1X,I0,1X,I0,1X,ES26.17E3,1X,ES26.17E3,1X,I0,1X,I0)') &
                'EDIT',icase,k,result_d,result_l,dok_nnz(dok),lil_nnz(lil)
        end do
        call require(valid_dok(dok) .and. valid_lil(lil), 'final format invariants')
        dok2 = dok_transpose(dok,info)
        call require(info == 0 .and. valid_dok(dok2), 'DOK transpose')
        lil2 = lil_transpose(lil,info)
        call require(info == 0 .and. valid_lil(lil2), 'LIL transpose')
        dd = dok_to_dense(dok,info)
        call require(info == 0, 'DOK dense')
        dl = lil_to_dense(lil,info)
        call require(info == 0, 'LIL dense')
        call require(all(dd == dl), 'dense agreement')
        call require(all(dok_to_dense(dok2) == transpose(dd)), &
            'DOK transpose values')
        call require(all(lil_to_dense(lil2) == transpose(dl)), &
            'LIL transpose values')
        call require(dok_nnz(dok2) == dok_nnz(dok), 'DOK transpose structure')
        call require(lil_nnz(lil2) == lil_nnz(lil), 'LIL transpose structure')
        csr_d = dok_to_csr(dok,info)
        call require(info == 0, 'DOK CSR')
        csr_l = lil_to_csr(lil,info)
        call require(info == 0, 'LIL CSR')
        yd = dok_matvec(dok,x,info)
        call require(info == 0, 'DOK matvec')
        yl = lil_matvec(lil,x,info)
        call require(info == 0, 'LIL matvec')
        write(*,'(A,1X,I0,1X,I0,1X,I0,1X,I0,1X,I0)') &
            'CASE',icase,m,n,dok_nnz(dok),lil_nnz(lil)
        do i = 1,m
            do j = 1,n
                write(*,'(A,1X,I0,1X,I0,1X,ES26.17E3,1X,ES26.17E3)') &
                    'D',i,j,dd(i,j),dl(i,j)
            end do
            write(*,'(A,1X,I0,1X,ES26.17E3,1X,ES26.17E3)') 'Y',i,yd(i),yl(i)
        end do
        do k = 1, size(csr_d%data)
            write(*,'(A,1X,I0,1X,ES26.17E3)') 'CSRD',csr_d%indices(k),csr_d%data(k)
        end do
        do i = 1, size(csr_d%indptr)
            write(*,'(A,1X,I0)') 'PTRD',csr_d%indptr(i)
        end do
        do k = 1, size(csr_l%data)
            write(*,'(A,1X,I0,1X,ES26.17E3)') 'CSRL',csr_l%indices(k),csr_l%data(k)
        end do
        do i = 1, size(csr_l%indptr)
            write(*,'(A,1X,I0)') 'PTRL',csr_l%indptr(i)
        end do
        deallocate(rr,cc,vv,x,yd,yl,dd,dl)
    end do
contains
    subroutine require(ok,message)
        logical, intent(in) :: ok !! Test success condition.
        character(len=*), intent(in) :: message !! Test failure description.
        if (.not. ok) then
            print *, 'FAIL: ', message
            error stop 1
        end if
    end subroutine require
end program sparse_dok_lil_reference_driver
