! SPDX-License-Identifier: MIT
! Helper for independent SciPy comparisons; build with the sparse core objects.
program funm_krylov_driver
    use scifort_kinds, only : dp
    use scifort_sparse, only : csr_from_dense, csr_matrix
    use scifort_sparse_linalg, only : funm_multiply_krylov
    implicit none
    type(csr_matrix) :: a
    real(dp), allocatable :: dense(:,:), b(:), y(:), mf(:)
    real(dp) :: t, discrepancy
    integer :: n, m, cap, symmetry, which, status, mf_status, ios, cycles, products

    do
        read(*,*,iostat=ios) n,m,cap,symmetry,which,t
        if (ios < 0) exit
        if (ios /= 0) error stop 'bad case header'
        allocate(dense(n,n),b(n))
        read(*,*) dense
        read(*,*) b
        a = csr_from_dense(dense)
        call funm_multiply_krylov(a,b,matrix_function,y,status,t=t, &
            restart_every_m=m,max_restarts=cap,assume_hermitian=symmetry == 1, &
            cycles=cycles, products=products)
        call funm_multiply_krylov(n,b,dense_product,matrix_function,mf,mf_status,t=t, &
            restart_every_m=m,max_restarts=cap,assume_hermitian=symmetry == 1)
        discrepancy = 0.0_dp
        if (size(y) /= size(mf)) error stop 'matrix-free shape failure'
        if (size(y) > 0) discrepancy = maxval(abs(y-mf))
        if (status /= mf_status) error stop 'matrix-free status failure'
        write(*,'(3(i0,1x),es24.15)') status,cycles,products,discrepancy
        write(*,'(*(es25.16,1x))') y
        deallocate(dense,b,y,mf)
    end do
contains
    subroutine dense_product(x, z)
        real(dp), intent(in) :: x(:) !! Input vector for current matrix.
        real(dp), intent(out) :: z(:) !! Matrix-free multiplication.
        z = matmul(dense,x)
    end subroutine dense_product

    subroutine matrix_function(h, result, info)
        real(dp), intent(in) :: h(:,:) !! Krylov projected matrix.
        real(dp), intent(out) :: result(:,:) !! Matrix function on projection.
        integer, intent(out) :: info !! Zero on success.
        real(dp) :: hm(size(h,1),size(h,1)), term(size(h,1),size(h,1))
        integer :: i,j,scales
        info = 0
        select case(which)
        case(0)  ! exp(H), scaling/squaring Taylor for modest cases
            result = 0.0_dp
            do i=1,size(h,1)
                result(i,i) = 1.0_dp
            end do
            scales = 0
            if (size(h,1) > 0) scales = max(0,ceiling(log(max(1.0_dp,maxval(sum(abs(h),dim=1))))/log(2.0_dp)))
            hm = h/real(2**scales,dp)
            term = result
            do j = 1,90
                term = matmul(term,hm)/real(j,dp)
                result = result+term
                if (maxval(abs(term)) < 1.0e-16_dp) exit
            end do
            do j=1,scales
                result = matmul(result,result)
            end do
        case(1)  ! H^3-H
            result = matmul(h,matmul(h,h))-h
        case(2)  ! I+H+0.1 H^2
            result = h+0.1_dp*matmul(h,h)
            do i = 1,size(h,1)
                result(i,i) = result(i,i)+1.0_dp
            end do
        case(3)  ! sin(H), Taylor for bounded arguments
            result = h
            term = h
            hm = matmul(h,h)
            do j=1,60
                term = -matmul(term,hm)/(real(2*j,dp)*real(2*j+1,dp))
                result = result+term
                if (maxval(abs(term)) < 1.0e-16_dp) exit
            end do
        case default
            info = 1
        end select
    end subroutine matrix_function
end program funm_krylov_driver
