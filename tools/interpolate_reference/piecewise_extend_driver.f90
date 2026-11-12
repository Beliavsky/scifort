! SPDX-License-Identifier: MIT
! Batched public API regression driver for SciPy PPoly/BPoly.extend.
program piecewise_extend_driver
    use scifort_kinds, only : dp
    use scifort_interpolate, only : ppoly, bpoly, ppoly_type, bpoly_type, &
        complex_ppoly_type, complex_bpoly_type
    implicit none
    type(ppoly_type) :: pp
    type(bpoly_type) :: bp
    type(complex_ppoly_type) :: cpp
    type(complex_bpoly_type) :: cbp
    real(dp), allocatable :: x(:), xadded(:), oldr(:,:,:), oldi(:,:,:)
    real(dp), allocatable :: newr(:,:,:), newi(:,:,:), coeffr(:,:,:), knots(:), q(:)
    complex(dp), allocatable :: oldz(:,:,:), newz(:,:,:), coeffz(:,:,:)
    integer :: nc, icase, kind, deg0, seg0, deg1, seg1, nch, nq, status, i, der
    read(*,*) nc
    do icase=1,nc
        read(*,*) kind, deg0, seg0, deg1, seg1, nch, nq
        allocate(x(seg0+1),xadded(seg1),oldr(deg0+1,seg0,nch), &
            oldi(deg0+1,seg0,nch),newr(deg1+1,seg1,nch),newi(deg1+1,seg1,nch),q(nq))
        read(*,*) x
        read(*,*) oldr
        read(*,*) oldi
        read(*,*) xadded
        read(*,*) newr
        read(*,*) newi
        read(*,*) q
        if (kind<=2) then
            if (kind==1) then
                pp=ppoly(oldr,x,status=status)
                if (status/=0) error stop 'PPoly constructor failed'
                call pp%extend(newr,xadded,status)
                if (status/=0) error stop 'PPoly extend failed'
                coeffr=pp%coefficients()
                knots=pp%knots()
            else
                bp=bpoly(oldr,x,status=status)
                if (status/=0) error stop 'BPoly constructor failed'
                call bp%extend(newr,xadded,status)
                if (status/=0) error stop 'BPoly extend failed'
                coeffr=bp%coefficients()
                knots=bp%knots()
            end if
            write(*,'(a,1x,i0,1x,i0,1x,i0)') 'CASE',icase,size(coeffr,1),size(coeffr,2)
            write(*,'(*(1x,es24.16))') knots
            write(*,'(*(1x,es24.16))') reshape(coeffr,[size(coeffr)])
            do der=0,2
                do i=1,nq
                    if (kind==1) then
                        write(*,'(*(1x,es24.16))') pp%evaluate(q(i),nu=der)
                    else
                        write(*,'(*(1x,es24.16))') bp%evaluate(q(i),nu=der)
                    end if
                end do
            end do
            deallocate(coeffr)
        else
            oldz=cmplx(oldr,oldi,dp)
            newz=cmplx(newr,newi,dp)
            if (kind==3) then
                cpp=ppoly(oldz,x,status=status)
                if (status/=0) error stop 'complex PPoly constructor failed'
                call cpp%extend(newz,xadded,status)
                if (status/=0) error stop 'complex PPoly extend failed'
                coeffz=cpp%coefficients()
                knots=cpp%knots()
            else
                cbp=bpoly(oldz,x,status=status)
                if (status/=0) error stop 'complex BPoly constructor failed'
                call cbp%extend(newz,xadded,status)
                if (status/=0) error stop 'complex BPoly extend failed'
                coeffz=cbp%coefficients()
                knots=cbp%knots()
            end if
            write(*,'(a,1x,i0,1x,i0,1x,i0)') 'CASE',icase,size(coeffz,1),size(coeffz,2)
            write(*,'(*(1x,es24.16))') knots
            write(*,'(*(1x,es24.16))') reshape(real(coeffz,dp),[size(coeffz)])
            write(*,'(*(1x,es24.16))') reshape(aimag(coeffz),[size(coeffz)])
            do der=0,2
                do i=1,nq
                    if (kind==3) then
                        write(*,'(*(1x,es24.16))') real(cpp%evaluate(q(i),nu=der),dp)
                        write(*,'(*(1x,es24.16))') aimag(cpp%evaluate(q(i),nu=der))
                    else
                        write(*,'(*(1x,es24.16))') real(cbp%evaluate(q(i),nu=der),dp)
                        write(*,'(*(1x,es24.16))') aimag(cbp%evaluate(q(i),nu=der))
                    end if
                end do
            end do
            deallocate(coeffz,oldz,newz)
        end if
        deallocate(x,xadded,oldr,oldi,newr,newi,q,knots)
    end do
end program piecewise_extend_driver
