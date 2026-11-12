! SPDX-License-Identifier: MIT
! Reference driver for randomized complex Krylov function-action comparisons.
program sparse_complex_funm_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        csr_complex_from_dense, csr_complex_to_csc
    use scifort_sparse_linalg_complex_funm, only : &
        funm_multiply_krylov_complex, funm_multiply_krylov_complex_operator, &
        complex_funm_dense_function
    implicit none
    type(csr_complex_matrix) :: csr
    type(csc_complex_matrix) :: csc
    procedure(complex_funm_dense_function), pointer :: f
    complex(dp), allocatable :: dense(:,:), b(:,:), yvec(:), yblock(:,:)
    complex(dp) :: t
    integer :: ncases, i, j, k, testcase, n, nrhs, which, herm, m, cap
    integer :: statuses(6), status
    real(dp) :: r, im

    read(*,*) ncases
    do testcase = 1, ncases
        read(*,*) n, nrhs, which, herm, m, cap, r, im
        t = cmplx(r,im,dp)
        allocate(dense(n,n),b(n,nrhs))
        do i = 1,n
            do j = 1,n
                read(*,*) r, im
                dense(i,j) = cmplx(r,im,dp)
            end do
        end do
        do i = 1,n
            do j = 1,nrhs
                read(*,*) r, im
                b(i,j) = cmplx(r,im,dp)
            end do
        end do
        csr = csr_complex_from_dense(dense,status)
        if (status /= 0) error stop 'bad CSR'
        csc = csr_complex_to_csc(csr,status)
        if (status /= 0) error stop 'bad CSC'
        select case(which)
        case(0)
            f => poly_f
        case(1)
            f => exp_f
        case(2)
            f => sin_f
        case default
            error stop 'bad function id'
        end select
        call funm_multiply_krylov_complex(csr,b(:,1),f,yvec,status,t=t, &
            restart_every_m=m,max_restarts=cap,assume_hermitian=herm==1,rtol=1.0e-9_dp)
        statuses(1) = status
        call output_vector(yvec,n)
        call funm_multiply_krylov_complex(csc,b(:,1),f,yvec,status,t=t, &
            restart_every_m=m,max_restarts=cap,assume_hermitian=herm==1,rtol=1.0e-9_dp)
        statuses(2) = status
        call output_vector(yvec,n)
        call funm_multiply_krylov_complex_operator(n,b(:,1),dense_product,f,yvec,status,t=t, &
            restart_every_m=m,max_restarts=cap,assume_hermitian=herm==1,rtol=1.0e-9_dp)
        statuses(3) = status
        call output_vector(yvec,n)
        call funm_multiply_krylov_complex(csr,b,f,yblock,status,t=t, &
            restart_every_m=m,max_restarts=cap,assume_hermitian=herm==1,rtol=1.0e-9_dp)
        statuses(4) = status
        call output_block(yblock,n,nrhs)
        call funm_multiply_krylov_complex(csc,b,f,yblock,status,t=t, &
            restart_every_m=m,max_restarts=cap,assume_hermitian=herm==1,rtol=1.0e-9_dp)
        statuses(5) = status
        call output_block(yblock,n,nrhs)
        call funm_multiply_krylov_complex_operator(n,b,dense_product,f,yblock,status,t=t, &
            restart_every_m=m,max_restarts=cap,assume_hermitian=herm==1,rtol=1.0e-9_dp)
        statuses(6) = status
        call output_block(yblock,n,nrhs)
        write(*,'(6(i0,1x))') statuses
        deallocate(dense,b)
    end do
contains
    subroutine output_vector(x,nexpected)
        complex(dp), intent(in) :: x(:) !! Result vector (possibly empty).
        integer, intent(in) :: nexpected !! Number of entries in reference.
        integer :: row
        do row=1,nexpected
            if (size(x) == nexpected) then
                write(*,'(2(es25.17,1x))') real(x(row),dp),aimag(x(row))
            else
                write(*,'(2(es25.17,1x))') 0.0_dp,0.0_dp
            end if
        end do
    end subroutine output_vector

    subroutine output_block(x,nexpected,kexpected)
        complex(dp), intent(in) :: x(:,:) !! Result block.
        integer, intent(in) :: nexpected !! Number of rows.
        integer, intent(in) :: kexpected !! Number of columns.
        integer :: row, col
        do row=1,nexpected
            do col=1,kexpected
                if (size(x,1) == nexpected .and. size(x,2) == kexpected) then
                    write(*,'(2(es25.17,1x))') real(x(row,col),dp),aimag(x(row,col))
                else
                    write(*,'(2(es25.17,1x))') 0.0_dp,0.0_dp
                end if
            end do
        end do
    end subroutine output_block

    subroutine dense_product(x,ax,stat)
        complex(dp), intent(in) :: x(:) !! Matrix-free input.
        complex(dp), intent(out) :: ax(:) !! Dense reference matvec.
        integer, intent(out) :: stat !! Always zero.
        ax = matmul(dense,x)
        stat = 0
    end subroutine dense_product

    subroutine poly_f(h,fh,stat)
        complex(dp), intent(in) :: h(:,:) !! Small input matrix.
        complex(dp), intent(out) :: fh(:,:) !! Polynomial evaluation.
        integer, intent(out) :: stat !! Status.
        integer :: j
        fh = matmul(h,h) + cmplx(0.2_dp,-0.3_dp,dp)*h
        do j=1,size(h,1)
            fh(j,j) = fh(j,j) + cmplx(0.75_dp,0.1_dp,dp)
        end do
        stat=0
    end subroutine poly_f

    subroutine exp_f(h,fh,stat)
        complex(dp), intent(in) :: h(:,:) !! Small matrix.
        complex(dp), intent(out) :: fh(:,:) !! exp(h).
        integer, intent(out) :: stat !! Status.
        call small_exp(h,fh)
        stat=0
    end subroutine exp_f

    subroutine sin_f(h,fh,stat)
        complex(dp), intent(in) :: h(:,:) !! Small matrix.
        complex(dp), intent(out) :: fh(:,:) !! sin(h).
        integer, intent(out) :: stat !! Status.
        complex(dp) :: eplus(size(h,1),size(h,2)), eminus(size(h,1),size(h,2))
        call small_exp(cmplx(0.0_dp,1.0_dp,dp)*h,eplus)
        call small_exp(cmplx(0.0_dp,-1.0_dp,dp)*h,eminus)
        fh = (eplus-eminus)/cmplx(0.0_dp,2.0_dp,dp)
        stat=0
    end subroutine sin_f

    subroutine small_exp(h,fh)
        complex(dp), intent(in) :: h(:,:) !! Dense projected matrix.
        complex(dp), intent(out) :: fh(:,:) !! Exponential by scaling and Taylor series.
        complex(dp) :: term(size(h,1),size(h,2)), scaled(size(h,1),size(h,2))
        real(dp) :: norm_bound
        integer :: n, row, k, squarings
        n=size(h,1)
        norm_bound=0.0_dp
        do row=1,n
            norm_bound = max(norm_bound,sum(abs(h(row,:))))
        end do
        squarings=0
        if (norm_bound > 0.5_dp) squarings=max(0,ceiling(log(norm_bound/0.5_dp)/log(2.0_dp)))
        squarings=min(squarings,30)
        scaled=h*2.0_dp**(-squarings)
        fh=cmplx(0.0_dp,0.0_dp,dp)
        term=cmplx(0.0_dp,0.0_dp,dp)
        do row=1,n
            fh(row,row)=cmplx(1.0_dp,0.0_dp,dp)
            term(row,row)=cmplx(1.0_dp,0.0_dp,dp)
        end do
        do k=1,100
            term=matmul(term,scaled)/real(k,dp)
            fh=fh+term
            if (maxval(abs(term))<=2.0e-16_dp*max(1.0_dp,maxval(abs(fh)))) exit
        end do
        do row=1,squarings
            fh=matmul(fh,fh)
        end do
    end subroutine small_exp
end program sparse_complex_funm_reference_driver
