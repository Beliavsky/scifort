! SPDX-License-Identifier: MIT
program bernstein_hermite_driver
    use scifort_kinds, only : dp
    use scifort_interpolate, only : bpoly_type, bpoly_from_derivatives
    implicit none
    type(bpoly_type) :: bp
    real(dp), allocatable :: x(:), derivatives(:,:,:), q(:), c(:,:,:)
    integer, allocatable :: counts(:)
    integer :: nt, icase, n, nd, nch, nq, iex, status, i, j
    read(*,*) nt
    do icase=1,nt
        read(*,*) n,nd,nch,nq,iex
        allocate(x(n),counts(n),derivatives(n,nd,nch),q(nq))
        read(*,*) x
        read(*,*) counts
        read(*,*) derivatives
        read(*,*) q
        bp=bpoly_from_derivatives(x,derivatives,counts,extrapolate=iex==1,status=status)
        if (status/=0 .or. .not. bp%valid()) error stop 'Hermite constructor failed'
        c=bp%coefficients()
        write(*,'(a,1x,i0,1x,i0)') 'CASE',icase,size(c,1)
        write(*,'(*(1x,es24.16))') c
        do j=0,3
            do i=1,nq
                write(*,'(*(1x,es24.16))') bp%evaluate(q(i),nu=j)
            end do
        end do
        write(*,'(*(1x,es24.16))') bp%integrate(x(1),x(n))
        write(*,'(*(1x,es24.16))') bp%integrate(x(n),x(1))
        deallocate(x,counts,derivatives,q,c)
    end do
end program bernstein_hermite_driver
