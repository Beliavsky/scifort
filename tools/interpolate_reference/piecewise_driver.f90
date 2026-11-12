! SPDX-License-Identifier: MIT
! Batch numerical reference driver for PPoly/BPoly and basis conversions.
program piecewise_driver
    use scifort_kinds, only : dp
    use scifort_interpolate, only : ppoly_type, bpoly_type, ppoly, bpoly, &
        ppoly_from_bernstein, bpoly_from_power
    implicit none
    type(ppoly_type) :: p, converted_p, pd, ap
    type(bpoly_type) :: b, converted_b, db, ab
    real(dp), allocatable :: x(:), cp(:,:,:), cb(:,:,:), query(:)
    real(dp), allocatable :: y(:)
    real(dp) :: a, z
    integer :: ncase, icase, k, m, nch, nq, nints, iper, iex, it, nu, id
    read(*,*) ncase
    do icase = 1,ncase
        read(*,*) k,m,nch,nq,nints,iper,iex
        allocate(x(m+1),cp(k+1,m,nch),cb(k+1,m,nch),query(nq))
        read(*,*) x
        read(*,*) cp
        read(*,*) cb
        read(*,*) query
        p=ppoly(cp,x,periodic=iper==1,extrapolate=iex==1)
        b=bpoly(cb,x,periodic=iper==1,extrapolate=iex==1)
        converted_p=ppoly_from_bernstein(b)
        converted_b=bpoly_from_power(p)
        if (.not. p%valid() .or. .not. b%valid() .or. &
            .not. converted_p%valid() .or. .not. converted_b%valid()) error stop 'constructor failed'
        write(*,'(a,1x,i0)') 'CASE',icase
        do id=1,4
            do nu=0,4
                do it=1,nq
                    select case(id)
                    case(1)
                        y=p%evaluate(query(it),nu=nu)
                    case(2)
                        y=b%evaluate(query(it),nu=nu)
                    case(3)
                        y=converted_p%evaluate(query(it),nu=nu)
                    case(4)
                        y=converted_b%evaluate(query(it),nu=nu)
                    end select
                    write(*,'(*(1x,es24.16))') y
                end do
            end do
            do it=1,nints
                read(*,*) a,z
                select case(id)
                case(1)
                    y=p%integrate(a,z)
                case(2)
                    y=b%integrate(a,z)
                case(3)
                    y=converted_p%integrate(a,z)
                case(4)
                    y=converted_b%integrate(a,z)
                end select
                write(*,'(*(1x,es24.16))') y
            end do
            select case(id)
            case(1)
                pd=p%derivative()
                ap=p%antiderivative()
                do it=1,nq
                    write(*,'(*(1x,es24.16))') pd%evaluate(query(it))
                    write(*,'(*(1x,es24.16))') ap%evaluate(query(it),extrapolate=.true.)
                end do
            case(2)
                db=b%derivative()
                ab=b%antiderivative()
                do it=1,nq
                    write(*,'(*(1x,es24.16))') db%evaluate(query(it))
                    write(*,'(*(1x,es24.16))') ab%evaluate(query(it),extrapolate=.true.)
                end do
            case(3)
                pd=converted_p%derivative()
                ap=converted_p%antiderivative()
                do it=1,nq
                    write(*,'(*(1x,es24.16))') pd%evaluate(query(it))
                    write(*,'(*(1x,es24.16))') ap%evaluate(query(it),extrapolate=.true.)
                end do
            case(4)
                db=converted_b%derivative()
                ab=converted_b%antiderivative()
                do it=1,nq
                    write(*,'(*(1x,es24.16))') db%evaluate(query(it))
                    write(*,'(*(1x,es24.16))') ab%evaluate(query(it),extrapolate=.true.)
                end do
            end select
        end do
        deallocate(x,cp,cb,query)
    end do
end program piecewise_driver
