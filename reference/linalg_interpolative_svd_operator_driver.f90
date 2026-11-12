! SPDX-License-Identifier: MIT
! Driver for test/compare_linalg_interpolative_svd_operator_scipy.py.
! stdin: m n type(0=real,1=complex) mode(0=fixed,1=tol) parameter seed
! followed by m rows containing n real values or n (real,imag) pairs.
program driver_linalg_interpolative_svd_operator
    use, intrinsic :: iso_fortran_env, only : int64
    use scifort_kinds, only : dp
    use scifort_linalg_interpolative_svd_operator, only : randomized_interpolative_svd_operator
    implicit none
    real(dp), allocatable :: a(:, :), u(:, :), s(:), v(:, :)
    complex(dp), allocatable :: z(:, :), uc(:, :), vc(:, :)
    integer :: m, n, typ, mode, i, j, k, status, nforward, nadjoint
    real(dp) :: param, est
    integer(int64) :: seed

    read(*,*) m,n,typ,mode,param,seed
    if (m < 0 .or. n < 0) error stop 'negative dimension'
    allocate(a(m,n),z(m,n))
    if (typ == 0) then
        do i = 1,m
            if (n > 0) read(*,*) a(i,:)
        end do
    else
        do i = 1,m
            if (n > 0) read(*,*) z(i,:)
        end do
    end if
    nforward = 0
    nadjoint = 0
    est = -1.0_dp
    if (typ == 0) then
        if (mode == 0) then
            k = int(param)
            allocate(u(m,k),s(k),v(n,k))
            call randomized_interpolative_svd_operator(m,n,k,mv,rmv,u,s,v,status,seed=seed)
        else
            call randomized_interpolative_svd_operator(m,n,param,mv,rmv,k,u,s,v, &
                status,seed=seed,estimated_error=est)
        end if
    else
        if (mode == 0) then
            k = int(param)
            allocate(uc(m,k),s(k),vc(n,k))
            call randomized_interpolative_svd_operator(m,n,k,mvc,rmvc,uc,s,vc,status,seed=seed)
        else
            call randomized_interpolative_svd_operator(m,n,param,mvc,rmvc,k,uc,s,vc, &
                status,seed=seed,estimated_error=est)
        end if
    end if
    write(*,'(a,1x,i0,1x,i0,1x,es25.16e3)') 'RESULT',status,k,est
    write(*,'(a,1x,i0,1x,i0)') 'CALLS',nforward,nadjoint
    do j = 1,k
        write(*,'(a,1x,i0,1x,es25.16e3)') 'S',j,s(j)
    end do
    if (typ == 0) then
        do j = 1,k
            do i = 1,m
                write(*,'(a,1x,i0,1x,i0,1x,es25.16e3)') 'U',i,j,u(i,j)
            end do
            do i = 1,n
                write(*,'(a,1x,i0,1x,i0,1x,es25.16e3)') 'V',i,j,v(i,j)
            end do
        end do
    else
        do j = 1,k
            do i = 1,m
                write(*,'(a,1x,i0,1x,i0,1x,es25.16e3,1x,es25.16e3)') &
                    'U',i,j,real(uc(i,j),dp),aimag(uc(i,j))
            end do
            do i = 1,n
                write(*,'(a,1x,i0,1x,i0,1x,es25.16e3,1x,es25.16e3)') &
                    'V',i,j,real(vc(i,j),dp),aimag(vc(i,j))
            end do
        end do
    end if
contains
    subroutine mv(x,y,istat)
        real(dp), intent(in) :: x(:) !! input length n
        real(dp), intent(out) :: y(:) !! output length m
        integer, intent(out) :: istat !! 0 on success
        nforward = nforward+1
        y = matmul(a,x)
        istat = 0
    end subroutine mv
    subroutine rmv(x,y,istat)
        real(dp), intent(in) :: x(:) !! input length m
        real(dp), intent(out) :: y(:) !! output length n
        integer, intent(out) :: istat !! 0 on success
        nadjoint = nadjoint+1
        y = matmul(transpose(a),x)
        istat = 0
    end subroutine rmv
    subroutine mvc(x,y,istat)
        complex(dp), intent(in) :: x(:) !! input length n
        complex(dp), intent(out) :: y(:) !! output length m
        integer, intent(out) :: istat !! 0 on success
        nforward = nforward+1
        y = matmul(z,x)
        istat = 0
    end subroutine mvc
    subroutine rmvc(x,y,istat)
        complex(dp), intent(in) :: x(:) !! input length m
        complex(dp), intent(out) :: y(:) !! output length n
        integer, intent(out) :: istat !! 0 on success
        nadjoint = nadjoint+1
        y = matmul(transpose(conjg(z)),x)
        istat = 0
    end subroutine rmvc
end program driver_linalg_interpolative_svd_operator
