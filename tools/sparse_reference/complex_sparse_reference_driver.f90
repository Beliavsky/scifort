! SPDX-License-Identifier: MIT
! Reference runner: multiple cases, one per input block. Each matrix is
! provided as possibly duplicated, unsorted complex COO triplets.
program complex_sparse_reference_driver
    use scifort_kinds, only : dp
    use scifort_sparse_complex
    use scifort_sparse_linalg_complex
    implicit none
    integer :: n, m, nz, solver, format_id, use_jacobi, restart, cap, ios, k, info, steps
    integer, allocatable :: rows(:), cols(:)
    complex(dp), allocatable :: data(:), b(:), x0(:), x(:)
    type(csr_complex_matrix) :: csr
    type(csc_complex_matrix) :: csc
    real(dp) :: re, im, residual

    do
        read(*,*,iostat=ios) n,m,nz,solver,format_id,use_jacobi,restart,cap
        if (ios < 0) exit
        if (ios /= 0) error stop 'bad input header'
        allocate(rows(nz),cols(nz),data(nz),b(n),x0(m))
        do k = 1, nz
            read(*,*) rows(k),cols(k),re,im
            data(k) = cmplx(re,im,dp)
        end do
        do k = 1, n
            read(*,*) re,im
            b(k) = cmplx(re,im,dp)
        end do
        do k = 1, m
            read(*,*) re,im
            x0(k) = cmplx(re,im,dp)
        end do
        csr = csr_complex_from_coo(n,m,rows,cols,data,info)
        if (info /= 0) error stop 'CSR construction failed'
        csc = csr_complex_to_csc(csr,info)
        if (info /= 0) error stop 'CSC conversion failed'
        select case (solver)
        case (0) ! Matrix-vector product.
            if (format_id == 1) then
                x = csr_complex_matvec(csr,x0,info)
            else
                x = csc_complex_matvec(csc,x0,info)
            end if
            residual = 0.0_dp
            steps = 0
        case (1) ! Hermitian positive definite CG.
            if (format_id == 1) then
                call cg_complex(csr,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                    maxiter=cap,x0=x0,jacobi=(use_jacobi/=0),iterations=steps,residual_norm=residual)
            else
                call cg_complex(csc,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                    maxiter=cap,x0=x0,jacobi=(use_jacobi/=0),iterations=steps,residual_norm=residual)
            end if
        case (2) ! BiCGSTAB.
            if (format_id == 1) then
                call bicgstab_complex(csr,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                    maxiter=cap,x0=x0,jacobi=(use_jacobi/=0),iterations=steps,residual_norm=residual)
            else
                call bicgstab_complex(csc,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                    maxiter=cap,x0=x0,jacobi=(use_jacobi/=0),iterations=steps,residual_norm=residual)
            end if
        case (3) ! Right-preconditioned restarted GMRES.
            if (format_id == 1) then
                call gmres_complex(csr,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                    maxiter=cap,restart=restart,x0=x0,jacobi=(use_jacobi/=0), &
                    iterations=steps,residual_norm=residual)
            else
                call gmres_complex(csc,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                    maxiter=cap,restart=restart,x0=x0,jacobi=(use_jacobi/=0), &
                    iterations=steps,residual_norm=residual)
            end if
        case default
            error stop 'invalid solver'
        end select
        write(*,'(i0,1x,i0,1x,es24.16)') info,steps,residual
        do k = 1, size(x)
            write(*,'(es24.16,1x,es24.16)') real(x(k),dp),aimag(x(k))
        end do
        deallocate(rows,cols,data,b,x0,x)
    end do
end program complex_sparse_reference_driver
