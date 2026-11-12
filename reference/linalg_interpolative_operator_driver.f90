! SPDX-License-Identifier: MIT
! Numerical driver for test/compare_linalg_interpolative_operator_scipy.py.
! stdin: m n type(0=real/1=complex) mode(0=fixed/1=tol) parameter seed
! followed by m lines of n real values or n (real,imag) complex tuples.
program driver_linalg_interpolative_operator
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg_interpolative_operator, only : randomized_interp_decomp_operator, &
        randomized_interp_decomp_operator_complex
    implicit none
    real(dp), allocatable :: a(:, :), proj(:, :)
    complex(dp), allocatable :: z(:, :), projz(:, :)
    integer, allocatable :: idx(:)
    integer :: m, n, typ, mode, i, j, k, status, n_forward, n_adjoint
    real(dp) :: param, est
    integer(int64) :: seed
    read(*,*) m,n,typ,mode,param,seed
    allocate(a(m,n),z(m,n))
    if (typ == 0) then
        do i=1,m
            if (n > 0) read(*,*) a(i,:)
        end do
    else
        do i=1,m
            if (n > 0) read(*,*) z(i,:)
        end do
    end if
    n_forward = 0
    n_adjoint = 0
    est = -1.0_dp
    if (typ == 0) then
        if (mode == 0) then
            k = int(param)
            call randomized_interp_decomp_operator(m,n,k,forward_real,adjoint_real, &
                idx,proj,status,seed=seed)
        else
            call randomized_interp_decomp_operator(m,n,param,forward_real,adjoint_real, &
                k,idx,proj,status,seed=seed,estimated_error=est)
        end if
    else
        if (mode == 0) then
            k = int(param)
            call randomized_interp_decomp_operator_complex(m,n,k,forward_complex,adjoint_complex, &
                idx,projz,status,seed=seed)
        else
            call randomized_interp_decomp_operator_complex(m,n,param,forward_complex,adjoint_complex, &
                k,idx,projz,status,seed=seed,estimated_error=est)
        end if
    end if
    write(*,'(a,1x,i0,1x,i0,1x,es24.16e3)') 'RESULT',status,k,est
    write(*,'(a,1x,i0,1x,i0)') 'CALLS',n_forward,n_adjoint
    write(*,'(a)') 'INDICES'
    write(*,'(*(i0,1x))') idx
    write(*,'(a)') 'PROJ'
    if (typ == 0) then
        do j=1,size(proj,2)
            do i=1,size(proj,1)
                write(*,'(2(i0,1x),es25.16e3)') i,j,proj(i,j)
            end do
        end do
    else
        do j=1,size(projz,2)
            do i=1,size(projz,1)
                write(*,'(2(i0,1x),2(es25.16e3,1x))') i,j,real(projz(i,j),dp),aimag(projz(i,j))
            end do
        end do
    end if
contains
    subroutine forward_real(x,y,istat)
        real(dp), intent(in) :: x(:) !! input
        real(dp), intent(out) :: y(:) !! output
        integer, intent(out) :: istat !! callback status
        n_forward = n_forward + 1
        y = matmul(a,x)
        istat = 0
    end subroutine forward_real
    subroutine adjoint_real(x,y,istat)
        real(dp), intent(in) :: x(:) !! input
        real(dp), intent(out) :: y(:) !! output
        integer, intent(out) :: istat !! callback status
        n_adjoint = n_adjoint + 1
        y = matmul(transpose(a),x)
        istat = 0
    end subroutine adjoint_real
    subroutine forward_complex(x,y,istat)
        complex(dp), intent(in) :: x(:) !! input
        complex(dp), intent(out) :: y(:) !! output
        integer, intent(out) :: istat !! callback status
        n_forward = n_forward + 1
        y = matmul(z,x)
        istat = 0
    end subroutine forward_complex
    subroutine adjoint_complex(x,y,istat)
        complex(dp), intent(in) :: x(:) !! input
        complex(dp), intent(out) :: y(:) !! output
        integer, intent(out) :: istat !! callback status
        n_adjoint = n_adjoint + 1
        y = matmul(transpose(conjg(z)),x)
        istat = 0
    end subroutine adjoint_complex
end program driver_linalg_interpolative_operator
