! SPDX-License-Identifier: MIT
! Streaming reference driver: dense input is used only to assemble a CSR/CSC
! matrix and supply an independent matrix-free matvec callback.
program sparse_complex_expm_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        csr_complex_from_dense, csr_complex_to_csc
    use scifort_sparse_linalg_complex_expm, only : expm_multiply_complex, &
        expm_multiply_complex_operator, expm_multiply_complex_interval, &
        expm_multiply_complex_operator_interval
    implicit none
    complex(dp), allocatable :: dense(:,:), b(:,:), zcsr(:,:), zcsc(:,:), zop(:,:)
    complex(dp), allocatable :: series_csr(:,:,:), series_csc(:,:,:), series_op(:,:,:)
    type(csr_complex_matrix) :: a
    type(csc_complex_matrix) :: c
    integer :: cases, n, nrhs, case_no, i, j, method, info, num, include_stop, sample
    real(dp) :: re, im, time, bound, start_time, stop_time
    complex(dp) :: tracea

    read(*,*) cases
    do case_no = 1,cases
        read(*,*) n, nrhs, time, start_time, stop_time, num, include_stop
        allocate(dense(n,n),b(n,nrhs))
        do i = 1,n
            do j = 1,n
                read(*,*) re, im
                dense(i,j) = cmplx(re,im,dp)
            end do
        end do
        do i = 1,n
            do j = 1,nrhs
                read(*,*) re, im
                b(i,j) = cmplx(re,im,dp)
            end do
        end do
        a = csr_complex_from_dense(dense,info)
        if (info /= 0) error stop 'invalid randomized CSR input'
        c = csr_complex_to_csc(a,info)
        if (info /= 0) error stop 'invalid randomized CSC input'
        tracea = cmplx(0.0_dp,0.0_dp,dp)
        do i = 1,n
            tracea = tracea+dense(i,i)
        end do
        bound = maxval(sum(abs(dense),dim=2))
        call expm_multiply_complex(a,b,zcsr,info,t=time)
        if (info /= 0) then
            print '(a,i0)', 'CSR error: ', info
            error stop 1
        end if
        call expm_multiply_complex(c,b,zcsc,info,t=time,tracea=tracea)
        if (info /= 0) then
            print '(a,i0)', 'CSC error: ', info
            error stop 1
        end if
        call expm_multiply_complex_operator(n,b,apply_dense,bound,zop,info,t=time,tracea=tracea)
        if (info /= 0) then
            print '(a,i0)', 'operator error: ', info
            error stop 1
        end if
        call expm_multiply_complex_interval(a,b,start_time,stop_time,num,series_csr,info, &
            endpoint=include_stop == 1)
        if (info /= 0) error stop 'CSR interval error'
        call expm_multiply_complex_interval(c,b,start_time,stop_time,num,series_csc,info, &
            endpoint=include_stop == 1)
        if (info /= 0) error stop 'CSC interval error'
        call expm_multiply_complex_operator_interval(n,b,apply_dense,bound,start_time,stop_time,num, &
            series_op,info,endpoint=include_stop == 1,tracea=tracea)
        if (info /= 0) error stop 'operator interval error'
        do method = 1,3
            do i = 1,n
                do j = 1,nrhs
                    select case (method)
                    case (1)
                        print '(2(es26.17e3,1x))', real(zcsr(i,j),dp), aimag(zcsr(i,j))
                    case (2)
                        print '(2(es26.17e3,1x))', real(zcsc(i,j),dp), aimag(zcsc(i,j))
                    case (3)
                        print '(2(es26.17e3,1x))', real(zop(i,j),dp), aimag(zop(i,j))
                    end select
                end do
            end do
        end do
        do method = 1,3
            do sample = 1,num
                do i = 1,n
                    do j = 1,nrhs
                        select case (method)
                        case (1)
                            print '(2(es26.17e3,1x))', real(series_csr(i,j,sample),dp), aimag(series_csr(i,j,sample))
                        case (2)
                            print '(2(es26.17e3,1x))', real(series_csc(i,j,sample),dp), aimag(series_csc(i,j,sample))
                        case (3)
                            print '(2(es26.17e3,1x))', real(series_op(i,j,sample),dp), aimag(series_op(i,j,sample))
                        end select
                    end do
                end do
            end do
        end do
        deallocate(dense,b)
    end do
contains
    subroutine apply_dense(x, y, status)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: status
        y = matmul(dense,x)
        status = 0
    end subroutine apply_dense
end program sparse_complex_expm_reference_driver
