! SPDX-License-Identifier: MIT
! Streaming reference driver for complex ILU(0), BiCG and QMR.
program complex_sparse_ilu_qmr_driver
    use scifort_kinds, only : dp
    use scifort_sparse_complex, only : csr_complex_matrix, csc_complex_matrix, &
        csr_complex_from_coo, csr_complex_to_csc
    use scifort_sparse_complex_ilu, only : complex_ilu0_factor, ilu0_complex, ilu0_complex_solve
    use scifort_sparse_linalg_complex, only : bicg_complex, qmr_complex, &
        bicgstab_complex, gmres_complex
    implicit none
    integer :: n, nnz, solver, format_id, mode, cap, ios, k, info, steps
    integer, allocatable :: rows(:), cols(:)
    complex(dp), allocatable :: values(:), b(:), x0(:), x(:)
    type(csr_complex_matrix) :: csr
    type(csc_complex_matrix) :: csc
    type(complex_ilu0_factor) :: factor
    real(dp) :: re, im, residual

    do
        read(*,*,iostat=ios) n,nnz,solver,format_id,mode,cap
        if (ios < 0) exit
        if (ios /= 0) error stop 'invalid header'
        allocate(rows(nnz),cols(nnz),values(nnz),b(n),x0(n))
        do k = 1, nnz
            read(*,*) rows(k),cols(k),re,im
            values(k) = cmplx(re,im,dp)
        end do
        do k = 1, n
            read(*,*) re,im
            b(k) = cmplx(re,im,dp)
        end do
        do k = 1, n
            read(*,*) re,im
            x0(k) = cmplx(re,im,dp)
        end do
        csr = csr_complex_from_coo(n,n,rows,cols,values,info)
        if (info /= 0) error stop 'bad CSR input'
        csc = csr_complex_to_csc(csr,info)
        if (info /= 0) error stop 'bad CSC conversion'
        steps = 0
        residual = 0.0_dp
        if (mode == 2) then
            call ilu0_complex(csr,factor,info)
            if (info /= 0) error stop 'bad ILU input'
        end if
        if (solver == 5) then
            if (mode /= 2) error stop 'ILU solve needs factor'
            call ilu0_complex_solve(factor,b,x,info,adjoint=(format_id==2))
        else
            if (format_id == 1) then
                if (mode == 2) then
                    select case (solver)
                    case (1)
                        call bicg_complex(csr,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                            maxiter=cap,x0=x0,ilu=factor,iterations=steps,residual_norm=residual)
                    case (2)
                        call qmr_complex(csr,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                            maxiter=cap,x0=x0,ilu=factor,iterations=steps,residual_norm=residual)
                    case (3)
                        call bicgstab_complex(csr,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                            maxiter=cap,x0=x0,ilu=factor,iterations=steps,residual_norm=residual)
                    case (4)
                        call gmres_complex(csr,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                            maxiter=cap,restart=min(20,n),x0=x0,ilu=factor, &
                            iterations=steps,residual_norm=residual)
                    end select
                else
                    select case (solver)
                    case (1)
                        call bicg_complex(csr,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                            maxiter=cap,x0=x0,jacobi=(mode==1),iterations=steps,residual_norm=residual)
                    case (2)
                        call qmr_complex(csr,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                            maxiter=cap,x0=x0,jacobi=(mode==1),iterations=steps,residual_norm=residual)
                    case (3)
                        call bicgstab_complex(csr,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                            maxiter=cap,x0=x0,jacobi=(mode==1),iterations=steps,residual_norm=residual)
                    case (4)
                        call gmres_complex(csr,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                            maxiter=cap,restart=min(20,n),x0=x0,jacobi=(mode==1), &
                            iterations=steps,residual_norm=residual)
                    end select
                end if
            else
                if (mode == 2) then
                    select case (solver)
                    case (1)
                        call bicg_complex(csc,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                            maxiter=cap,x0=x0,ilu=factor,iterations=steps,residual_norm=residual)
                    case (2)
                        call qmr_complex(csc,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                            maxiter=cap,x0=x0,ilu=factor,iterations=steps,residual_norm=residual)
                    case (3)
                        call bicgstab_complex(csc,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                            maxiter=cap,x0=x0,ilu=factor,iterations=steps,residual_norm=residual)
                    case (4)
                        call gmres_complex(csc,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                            maxiter=cap,restart=min(20,n),x0=x0,ilu=factor, &
                            iterations=steps,residual_norm=residual)
                    end select
                else
                    select case (solver)
                    case (1)
                        call bicg_complex(csc,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                            maxiter=cap,x0=x0,jacobi=(mode==1),iterations=steps,residual_norm=residual)
                    case (2)
                        call qmr_complex(csc,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                            maxiter=cap,x0=x0,jacobi=(mode==1),iterations=steps,residual_norm=residual)
                    case (3)
                        call bicgstab_complex(csc,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                            maxiter=cap,x0=x0,jacobi=(mode==1),iterations=steps,residual_norm=residual)
                    case (4)
                        call gmres_complex(csc,b,x,info,rtol=1.0e-10_dp,atol=1.0e-12_dp, &
                            maxiter=cap,restart=min(20,n),x0=x0,jacobi=(mode==1), &
                            iterations=steps,residual_norm=residual)
                    end select
                end if
            end if
        end if
        write(*,'(i0,1x,i0,1x,es24.16)') info,steps,residual
        do k = 1, size(x)
            write(*,'(es24.16,1x,es24.16)') real(x(k),dp),aimag(x(k))
        end do
        deallocate(rows,cols,values,b,x0,x)
    end do
end program complex_sparse_ilu_qmr_driver
