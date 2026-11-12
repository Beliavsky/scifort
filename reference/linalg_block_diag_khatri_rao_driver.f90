! SPDX-License-Identifier: MIT
! Original MIT test driver for tools/compare_linalg_block_diag_khatri_rao.py
! Protocol: number of cases; mode B or K; metadata and complex input entries.
program linalg_block_diag_khatri_rao_driver
    use, intrinsic :: iso_fortran_env, only : int32, int64
    use scifort_kinds, only : dp
    use scifort_linalg_block_diag, only : linalg_block, as_block, block_diag
    use scifort_linalg_khatri_rao, only : khatri_rao
    implicit none
    integer :: count, case_id
    character(1) :: mode

    read(*,*) count
    do case_id=1,count
        read(*,*) mode
        select case(mode)
        case('B')
            call run_block()
        case('K')
            call run_khatri()
        case default
            error stop 'unknown mode'
        end select
    end do
contains

    subroutine emit(z)
        complex(dp), intent(in) :: z !! value to write
        write(*,'(es26.17,1x,es26.17)') real(z,dp),aimag(z)
    end subroutine emit

    subroutine run_block()
        integer :: num, rnk, r, c, b, i, j, q, rnk_out, bo, no, ko, status
        character(1) :: typ, output_type
        type(linalg_block), allocatable :: blocks(:)
        complex(dp), allocatable :: temp(:,:,:)
        read(*,*) num,output_type,rnk_out,bo,no,ko
        allocate(blocks(num))
        do b=1,num
            read(*,*) typ,rnk,q,r,c
            allocate(temp(q,r,c))
            do i=1,q
                do j=1,r
                    block
                        integer :: k
                        do k=1,c
                            read(*,*) temp(i,j,k)
                        end do
                    end block
                end do
            end do
            select case(typ)
            case('J')
                select case(rnk)
                case(0)
                    blocks(b) = as_block(int(real(temp(1,1,1),dp),int32))
                case(1)
                    blocks(b) = as_block(int(real(temp(1,1,:),dp),int32))
                case(2)
                    blocks(b) = as_block(int(real(temp(1,:,:),dp),int32))
                case(3)
                    blocks(b) = as_block(int(real(temp,dp),int32))
                end select
            case('I')
                select case(rnk)
                case(0)
                    blocks(b) = as_block(int(real(temp(1,1,1),dp),int64))
                case(1)
                    blocks(b) = as_block(int(real(temp(1,1,:),dp),int64))
                case(2)
                    blocks(b) = as_block(int(real(temp(1,:,:),dp),int64))
                case(3)
                    blocks(b) = as_block(int(real(temp,dp),int64))
                end select
            case('R')
                select case(rnk)
                case(0)
                    blocks(b) = as_block(real(temp(1,1,1),dp))
                case(1)
                    blocks(b) = as_block(real(temp(1,1,:),dp))
                case(2)
                    blocks(b) = as_block(real(temp(1,:,:),dp))
                case(3)
                    blocks(b) = as_block(real(temp,dp))
                end select
            case('C')
                select case(rnk)
                case(0)
                    blocks(b) = as_block(temp(1,1,1))
                case(1)
                    blocks(b) = as_block(temp(1,1,:))
                case(2)
                    blocks(b) = as_block(temp(1,:,:))
                case(3)
                    blocks(b) = as_block(temp)
                end select
            end select
            deallocate(temp)
        end do
        select case(output_type)
        case('J')
            if(rnk_out == 2) then
                block
                    integer(int32), allocatable :: out(:,:)
                    allocate(out(no,ko))
                    call block_diag(blocks,out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do i=1,no
                            do j=1,ko
                                call emit(cmplx(real(out(i,j),dp),0.0_dp,dp))
                            end do
                        end do
                    end if
                end block
            else
                block
                    integer(int32), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call block_diag(blocks,out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                    call emit(cmplx(real(out(q,i,j),dp),0.0_dp,dp))
                                end do
                            end do
                        end do
                    end if
                end block
            end if
        case('I')
            if(rnk_out == 2) then
                block
                    integer(int64), allocatable :: out(:,:)
                    allocate(out(no,ko))
                    call block_diag(blocks,out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do i=1,no
                            do j=1,ko
                                call emit(cmplx(real(out(i,j),dp),0.0_dp,dp))
                            end do
                        end do
                    end if
                end block
            else
                block
                    integer(int64), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call block_diag(blocks,out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                    call emit(cmplx(real(out(q,i,j),dp),0.0_dp,dp))
                                end do
                            end do
                        end do
                    end if
                end block
            end if
        case('R')
            if(rnk_out == 2) then
                block
                    real(dp), allocatable :: out(:,:)
                    allocate(out(no,ko))
                    call block_diag(blocks,out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do i=1,no
                            do j=1,ko
                                call emit(cmplx(out(i,j),kind=dp))
                            end do
                        end do
                    end if
                end block
            else
                block
                    real(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call block_diag(blocks,out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                    call emit(cmplx(out(q,i,j),kind=dp))
                                end do
                            end do
                        end do
                    end if
                end block
            end if
        case('C')
            if(rnk_out == 2) then
                block
                    complex(dp), allocatable :: out(:,:)
                    allocate(out(no,ko))
                    call block_diag(blocks,out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do i=1,no
                            do j=1,ko
                                call emit(out(i,j))
                            end do
                        end do
                    end if
                end block
            else
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call block_diag(blocks,out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                    call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            end if
        end select
    end subroutine run_block

    subroutine run_khatri()
        integer :: ra, rb, ba, bb, n, m, k, bo, no, ko, i, j, q, status
        character(1) :: ta,tb
        complex(dp), allocatable :: aa(:,:,:), xb(:,:,:)
        read(*,*) ta,tb,ra,rb,ba,bb,n,m,k,bo,no,ko
        allocate(aa(ba,n,k),xb(bb,m,k))
        do q=1,ba
            do i=1,n
                do j=1,k
                    read(*,*) aa(q,i,j)
                end do
            end do
        end do
        do q=1,bb
            do i=1,m
                do j=1,k
                    read(*,*) xb(q,i,j)
                end do
            end do
        end do
        select case(ta//tb)
        case('JJ')
            if (ra == 2 .and. rb == 2) then
                block
                    integer(int32), allocatable :: out(:,:)
                    allocate(out(no,ko))
                    call khatri_rao(int(real(aa(1,:,:),dp),int32),int(real(xb(1,:,:),dp),int32),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do i=1,no
                            do j=1,ko
                                call emit(cmplx(real(out(i,j),dp),0.0_dp,dp))
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 3) then
                block
                    integer(int32), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa,dp),int32),int(real(xb,dp),int32),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(real(out(q,i,j),dp),0.0_dp,dp))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 2 .and. rb == 3) then
                block
                    integer(int32), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa(1,:,:),dp),int32),int(real(xb,dp),int32),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(real(out(q,i,j),dp),0.0_dp,dp))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 2) then
                block
                    integer(int32), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa,dp),int32),int(real(xb(1,:,:),dp),int32),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(real(out(q,i,j),dp),0.0_dp,dp))
                                end do
                            end do
                        end do
                    end if
                end block
            end if
        case('JI')
            if (ra == 2 .and. rb == 2) then
                block
                    integer(int64), allocatable :: out(:,:)
                    allocate(out(no,ko))
                    call khatri_rao(int(real(aa(1,:,:),dp),int32),int(real(xb(1,:,:),dp),int64),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do i=1,no
                            do j=1,ko
                                call emit(cmplx(real(out(i,j),dp),0.0_dp,dp))
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 3) then
                block
                    integer(int64), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa,dp),int32),int(real(xb,dp),int64),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(real(out(q,i,j),dp),0.0_dp,dp))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 2 .and. rb == 3) then
                block
                    integer(int64), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa(1,:,:),dp),int32),int(real(xb,dp),int64),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(real(out(q,i,j),dp),0.0_dp,dp))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 2) then
                block
                    integer(int64), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa,dp),int32),int(real(xb(1,:,:),dp),int64),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(real(out(q,i,j),dp),0.0_dp,dp))
                                end do
                            end do
                        end do
                    end if
                end block
            end if
        case('JR')
            if (ra == 2 .and. rb == 2) then
                block
                    real(dp), allocatable :: out(:,:)
                    allocate(out(no,ko))
                    call khatri_rao(int(real(aa(1,:,:),dp),int32),real(xb(1,:,:),dp),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do i=1,no
                            do j=1,ko
                                call emit(cmplx(out(i,j),kind=dp))
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 3) then
                block
                    real(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa,dp),int32),real(xb,dp),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(out(q,i,j),kind=dp))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 2 .and. rb == 3) then
                block
                    real(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa(1,:,:),dp),int32),real(xb,dp),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(out(q,i,j),kind=dp))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 2) then
                block
                    real(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa,dp),int32),real(xb(1,:,:),dp),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(out(q,i,j),kind=dp))
                                end do
                            end do
                        end do
                    end if
                end block
            end if
        case('JC')
            if (ra == 2 .and. rb == 2) then
                block
                    complex(dp), allocatable :: out(:,:)
                    allocate(out(no,ko))
                    call khatri_rao(int(real(aa(1,:,:),dp),int32),xb(1,:,:),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do i=1,no
                            do j=1,ko
                                call emit(out(i,j))
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 3) then
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa,dp),int32),xb,out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 2 .and. rb == 3) then
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa(1,:,:),dp),int32),xb,out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 2) then
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa,dp),int32),xb(1,:,:),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            end if
        case('IJ')
            if (ra == 2 .and. rb == 2) then
                block
                    integer(int64), allocatable :: out(:,:)
                    allocate(out(no,ko))
                    call khatri_rao(int(real(aa(1,:,:),dp),int64),int(real(xb(1,:,:),dp),int32),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do i=1,no
                            do j=1,ko
                                call emit(cmplx(real(out(i,j),dp),0.0_dp,dp))
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 3) then
                block
                    integer(int64), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa,dp),int64),int(real(xb,dp),int32),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(real(out(q,i,j),dp),0.0_dp,dp))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 2 .and. rb == 3) then
                block
                    integer(int64), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa(1,:,:),dp),int64),int(real(xb,dp),int32),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(real(out(q,i,j),dp),0.0_dp,dp))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 2) then
                block
                    integer(int64), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa,dp),int64),int(real(xb(1,:,:),dp),int32),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(real(out(q,i,j),dp),0.0_dp,dp))
                                end do
                            end do
                        end do
                    end if
                end block
            end if
        case('II')
            if (ra == 2 .and. rb == 2) then
                block
                    integer(int64), allocatable :: out(:,:)
                    allocate(out(no,ko))
                    call khatri_rao(int(real(aa(1,:,:),dp),int64),int(real(xb(1,:,:),dp),int64),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do i=1,no
                            do j=1,ko
                                call emit(cmplx(real(out(i,j),dp),0.0_dp,dp))
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 3) then
                block
                    integer(int64), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa,dp),int64),int(real(xb,dp),int64),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(real(out(q,i,j),dp),0.0_dp,dp))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 2 .and. rb == 3) then
                block
                    integer(int64), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa(1,:,:),dp),int64),int(real(xb,dp),int64),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(real(out(q,i,j),dp),0.0_dp,dp))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 2) then
                block
                    integer(int64), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa,dp),int64),int(real(xb(1,:,:),dp),int64),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(real(out(q,i,j),dp),0.0_dp,dp))
                                end do
                            end do
                        end do
                    end if
                end block
            end if
        case('IR')
            if (ra == 2 .and. rb == 2) then
                block
                    real(dp), allocatable :: out(:,:)
                    allocate(out(no,ko))
                    call khatri_rao(int(real(aa(1,:,:),dp),int64),real(xb(1,:,:),dp),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do i=1,no
                            do j=1,ko
                                call emit(cmplx(out(i,j),kind=dp))
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 3) then
                block
                    real(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa,dp),int64),real(xb,dp),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(out(q,i,j),kind=dp))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 2 .and. rb == 3) then
                block
                    real(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa(1,:,:),dp),int64),real(xb,dp),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(out(q,i,j),kind=dp))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 2) then
                block
                    real(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa,dp),int64),real(xb(1,:,:),dp),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(out(q,i,j),kind=dp))
                                end do
                            end do
                        end do
                    end if
                end block
            end if
        case('IC')
            if (ra == 2 .and. rb == 2) then
                block
                    complex(dp), allocatable :: out(:,:)
                    allocate(out(no,ko))
                    call khatri_rao(int(real(aa(1,:,:),dp),int64),xb(1,:,:),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do i=1,no
                            do j=1,ko
                                call emit(out(i,j))
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 3) then
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa,dp),int64),xb,out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 2 .and. rb == 3) then
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa(1,:,:),dp),int64),xb,out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 2) then
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(int(real(aa,dp),int64),xb(1,:,:),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            end if
        case('RJ')
            if (ra == 2 .and. rb == 2) then
                block
                    real(dp), allocatable :: out(:,:)
                    allocate(out(no,ko))
                    call khatri_rao(real(aa(1,:,:),dp),int(real(xb(1,:,:),dp),int32),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do i=1,no
                            do j=1,ko
                                call emit(cmplx(out(i,j),kind=dp))
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 3) then
                block
                    real(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(real(aa,dp),int(real(xb,dp),int32),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(out(q,i,j),kind=dp))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 2 .and. rb == 3) then
                block
                    real(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(real(aa(1,:,:),dp),int(real(xb,dp),int32),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(out(q,i,j),kind=dp))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 2) then
                block
                    real(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(real(aa,dp),int(real(xb(1,:,:),dp),int32),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(out(q,i,j),kind=dp))
                                end do
                            end do
                        end do
                    end if
                end block
            end if
        case('RI')
            if (ra == 2 .and. rb == 2) then
                block
                    real(dp), allocatable :: out(:,:)
                    allocate(out(no,ko))
                    call khatri_rao(real(aa(1,:,:),dp),int(real(xb(1,:,:),dp),int64),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do i=1,no
                            do j=1,ko
                                call emit(cmplx(out(i,j),kind=dp))
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 3) then
                block
                    real(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(real(aa,dp),int(real(xb,dp),int64),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(out(q,i,j),kind=dp))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 2 .and. rb == 3) then
                block
                    real(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(real(aa(1,:,:),dp),int(real(xb,dp),int64),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(out(q,i,j),kind=dp))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 2) then
                block
                    real(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(real(aa,dp),int(real(xb(1,:,:),dp),int64),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(out(q,i,j),kind=dp))
                                end do
                            end do
                        end do
                    end if
                end block
            end if
        case('RR')
            if (ra == 2 .and. rb == 2) then
                block
                    real(dp), allocatable :: out(:,:)
                    allocate(out(no,ko))
                    call khatri_rao(real(aa(1,:,:),dp),real(xb(1,:,:),dp),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do i=1,no
                            do j=1,ko
                                call emit(cmplx(out(i,j),kind=dp))
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 3) then
                block
                    real(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(real(aa,dp),real(xb,dp),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(out(q,i,j),kind=dp))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 2 .and. rb == 3) then
                block
                    real(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(real(aa(1,:,:),dp),real(xb,dp),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(out(q,i,j),kind=dp))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 2) then
                block
                    real(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(real(aa,dp),real(xb(1,:,:),dp),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(cmplx(out(q,i,j),kind=dp))
                                end do
                            end do
                        end do
                    end if
                end block
            end if
        case('RC')
            if (ra == 2 .and. rb == 2) then
                block
                    complex(dp), allocatable :: out(:,:)
                    allocate(out(no,ko))
                    call khatri_rao(real(aa(1,:,:),dp),xb(1,:,:),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do i=1,no
                            do j=1,ko
                                call emit(out(i,j))
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 3) then
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(real(aa,dp),xb,out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 2 .and. rb == 3) then
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(real(aa(1,:,:),dp),xb,out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 2) then
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(real(aa,dp),xb(1,:,:),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            end if
        case('CJ')
            if (ra == 2 .and. rb == 2) then
                block
                    complex(dp), allocatable :: out(:,:)
                    allocate(out(no,ko))
                    call khatri_rao(aa(1,:,:),int(real(xb(1,:,:),dp),int32),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do i=1,no
                            do j=1,ko
                                call emit(out(i,j))
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 3) then
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(aa,int(real(xb,dp),int32),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 2 .and. rb == 3) then
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(aa(1,:,:),int(real(xb,dp),int32),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 2) then
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(aa,int(real(xb(1,:,:),dp),int32),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            end if
        case('CI')
            if (ra == 2 .and. rb == 2) then
                block
                    complex(dp), allocatable :: out(:,:)
                    allocate(out(no,ko))
                    call khatri_rao(aa(1,:,:),int(real(xb(1,:,:),dp),int64),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do i=1,no
                            do j=1,ko
                                call emit(out(i,j))
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 3) then
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(aa,int(real(xb,dp),int64),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 2 .and. rb == 3) then
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(aa(1,:,:),int(real(xb,dp),int64),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 2) then
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(aa,int(real(xb(1,:,:),dp),int64),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            end if
        case('CR')
            if (ra == 2 .and. rb == 2) then
                block
                    complex(dp), allocatable :: out(:,:)
                    allocate(out(no,ko))
                    call khatri_rao(aa(1,:,:),real(xb(1,:,:),dp),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do i=1,no
                            do j=1,ko
                                call emit(out(i,j))
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 3) then
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(aa,real(xb,dp),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 2 .and. rb == 3) then
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(aa(1,:,:),real(xb,dp),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 2) then
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(aa,real(xb(1,:,:),dp),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            end if
        case('CC')
            if (ra == 2 .and. rb == 2) then
                block
                    complex(dp), allocatable :: out(:,:)
                    allocate(out(no,ko))
                    call khatri_rao(aa(1,:,:),xb(1,:,:),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do i=1,no
                            do j=1,ko
                                call emit(out(i,j))
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 3) then
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(aa,xb,out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 2 .and. rb == 3) then
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(aa(1,:,:),xb,out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            else if (ra == 3 .and. rb == 2) then
                block
                    complex(dp), allocatable :: out(:,:,:)
                    allocate(out(bo,no,ko))
                    call khatri_rao(aa,xb(1,:,:),out,status)
                    write(*,'(i0)') status
                    if (status == 0) then
                        do q=1,bo
                            do i=1,no
                                do j=1,ko
                                call emit(out(q,i,j))
                                end do
                            end do
                        end do
                    end if
                end block
            end if
        case default
            error stop 'unknown types in reference driver'
        end select
    end subroutine run_khatri
end program linalg_block_diag_khatri_rao_driver
