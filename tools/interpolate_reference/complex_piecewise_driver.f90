! SPDX-License-Identifier: MIT
! Input/output reference driver for complex PPoly/BPoly and basis conversion.
program complex_piecewise_driver
    use scifort_kinds, only : dp
    use scifort_interpolate, only : ppoly,bpoly,complex_ppoly_type,complex_bpoly_type, &
        complex_ppoly_from_bernstein,complex_bpoly_from_power
    implicit none
    type(complex_ppoly_type) :: p, p_from_b, p_deriv, p_anti
    type(complex_bpoly_type) :: b, b_from_p, b_deriv, b_anti
    real(dp), allocatable :: x(:),re(:,:,:),im(:,:,:),q(:)
    complex(dp), allocatable :: c(:,:,:), vals(:,:),integ(:)
    real(dp) :: a,bound
    integer :: i, cases, k, nseg, nch, nq,npairs, periodic, ext, j, nu, obj, idx, code

    read(*,*) cases
    do i=1,cases
        read(*,*) k,nseg,nch,nq,npairs,periodic,ext
        allocate(x(nseg+1),re(k+1,nseg,nch),im(k+1,nseg,nch),c(k+1,nseg,nch),q(nq))
        read(*,*) x
        read(*,*) re
        read(*,*) im
        c=cmplx(re,im,dp)
        read(*,*) q
        p=ppoly(c,x,periodic=periodic/=0,extrapolate=ext/=0,status=code)
        if (code/=0) error stop 'complex PPoly invalid'
        b=bpoly(c,x,periodic=periodic/=0,extrapolate=ext/=0,status=code)
        if (code/=0) error stop 'complex BPoly invalid'
        b_from_p=complex_bpoly_from_power(p)
        p_from_b=complex_ppoly_from_bernstein(b)
        p_deriv=p%derivative(1)
        b_deriv=b%derivative(1)
        p_anti=p%antiderivative(1)
        b_anti=b%antiderivative(1)
        write(*,'(a,1x,i0)') 'CASE',i
        do obj=1,4
            do nu=0,3
                select case(obj)
                case(1)
                    vals=p%evaluate(q,nu)
                case(2)
                    vals=b%evaluate(q,nu)
                case(3)
                    vals=b_from_p%evaluate(q,nu)
                case(4)
                    vals=p_from_b%evaluate(q,nu)
                end select
                call print_complex(vals)
            end do
        end do
        do obj=1,4
            do j=1,npairs
                read(*,*) a,bound
                select case(obj)
                case(1)
                    integ=p%integrate(a,bound)
                case(2)
                    integ=b%integrate(a,bound)
                case(3)
                    integ=b_from_p%integrate(a,bound)
                case(4)
                    integ=p_from_b%integrate(a,bound)
                end select
                call print_complex(reshape(integ,[1,nch]))
            end do
        end do
        call print_complex(p_deriv%evaluate(q))
        call print_complex(b_deriv%evaluate(q))
        call print_complex(p_anti%evaluate(q,extrapolate=.true.))
        call print_complex(b_anti%evaluate(q,extrapolate=.true.))
        deallocate(x,re,im,c,q)
    end do
contains
    subroutine print_complex(v)
        complex(dp), intent(in) :: v(:,:) !! complex array written in Fortran storage order
        write(*,'(*(1x,es26.17e3))') real(v,dp)
        write(*,'(*(1x,es26.17e3))') aimag(v)
    end subroutine print_complex
end program complex_piecewise_driver
