! SPDX-License-Identifier: MIT
! Line-oriented interface used by compare_linalg_spectral_norm_scipy.py.
! Each case specifies dtype, operation, input matrices and initial vector.
program driver_linalg_spectral_norm
    use scifort_kinds, only : dp
    use scifort_linalg_spectral_norm, only : estimate_spectral_norm, &
        estimate_spectral_norm_diff, estimate_spectral_norm_operator, &
        estimate_spectral_norm_operator_complex, estimate_spectral_norm_diff_operator, &
        estimate_spectral_norm_diff_operator_complex
    implicit none
    real(dp), allocatable :: a(:, :), b(:, :), start(:)
    complex(dp), allocatable :: c(:, :), d(:, :), start_c(:)
    real(dp) :: re, im, first, second
    integer :: ncases, testcase, kind, m, n, its, i, j, status1, status2
    read(*,*) ncases
    do testcase = 1, ncases
        read(*,*) kind, m, n, its
        allocate(a(m,n),b(m,n),c(m,n),d(m,n),start(n),start_c(n))
        do i = 1, m
            do j = 1, n
                read(*,*) re, im
                a(i,j) = re
                c(i,j) = cmplx(re,im,dp)
            end do
        end do
        do i = 1, m
            do j = 1, n
                read(*,*) re, im
                b(i,j) = re
                d(i,j) = cmplx(re,im,dp)
            end do
        end do
        do j = 1, n
            read(*,*) re, im
            start(j) = re
            start_c(j) = cmplx(re,im,dp)
        end do
        select case (kind)
        case (0)
            call estimate_spectral_norm(a,first,status1,its=its,start=start)
            call estimate_spectral_norm_operator(m,n,mv_a,rmv_a, &
                second,status2,its=its,start=start)
        case (1)
            call estimate_spectral_norm_diff(a,b,first,status1,its=its,start=start)
            call estimate_spectral_norm_diff_operator(m,n,mv_a,rmv_a,mv_b, &
                rmv_b,second,status2,its=its,start=start)
        case (2)
            call estimate_spectral_norm(c,first,status1,its=its,start=start_c)
            call estimate_spectral_norm_operator_complex(m,n,mv_c,rmv_c, &
                second,status2,its=its,start=start_c)
        case (3)
            call estimate_spectral_norm_diff(c,d,first,status1,its=its,start=start_c)
            call estimate_spectral_norm_diff_operator_complex(m,n,mv_c,rmv_c,mv_d, &
                rmv_d,second,status2,its=its,start=start_c)
        case default
            error stop 'unknown spectral norm test case'
        end select
        write(*,'(a,1x,i0,1x,i0,1x,i0,1x,es26.17,1x,es26.17)') &
            'CASE',testcase,status1,status2,first,second
        deallocate(a,b,c,d,start,start_c)
    end do
contains
    subroutine mv_a(v,y,status)
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: y(:)
        integer, intent(out) :: status
        y = matmul(a,v)
        status = 0
    end subroutine mv_a
    subroutine rmv_a(v,y,status)
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: y(:)
        integer, intent(out) :: status
        y = matmul(transpose(a),v)
        status = 0
    end subroutine rmv_a
    subroutine mv_b(v,y,status)
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: y(:)
        integer, intent(out) :: status
        y = matmul(b,v)
        status = 0
    end subroutine mv_b
    subroutine rmv_b(v,y,status)
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: y(:)
        integer, intent(out) :: status
        y = matmul(transpose(b),v)
        status = 0
    end subroutine rmv_b
    subroutine mv_c(v,y,status)
        complex(dp), intent(in) :: v(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: status
        y = matmul(c,v)
        status = 0
    end subroutine mv_c
    subroutine rmv_c(v,y,status)
        complex(dp), intent(in) :: v(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: status
        y = matmul(transpose(conjg(c)),v)
        status = 0
    end subroutine rmv_c
    subroutine mv_d(v,y,status)
        complex(dp), intent(in) :: v(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: status
        y = matmul(d,v)
        status = 0
    end subroutine mv_d
    subroutine rmv_d(v,y,status)
        complex(dp), intent(in) :: v(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: status
        y = matmul(transpose(conjg(d)),v)
        status = 0
    end subroutine rmv_d
end program driver_linalg_spectral_norm
