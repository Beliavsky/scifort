! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program aaa_large_driver
    use scifort_kinds, only : dp
    use scifort_interpolate_aaa, only : aaa,aaa_type,complex_aaa_type
    implicit none
    integer :: n,mode,clean,io
    do
        read(*,*,iostat=io) n,mode,clean
        if (io /= 0) exit
        block
            type(aaa_type) :: real_fit
            type(complex_aaa_type) :: complex_fit
            complex(dp), allocatable :: x(:),f(:),s(:),v(:),w(:),p(:),z(:),r(:)
            real(dp) :: xr,xi,fr,fi,ctol
            integer :: i,status,ps,zs,rs,removed
            allocate(x(n),f(n))
            read(*,*) ctol
            do i = 1, n
                read(*,*) xr,xi,fr,fi
                x(i) = cmplx(xr,xi,dp)
                f(i) = cmplx(fr,fi,dp)
            end do
            removed = 0
            if (mode == 0) then
                real_fit = aaa(real(x,dp),real(f,dp),rtol=0.0_dp,max_terms=n, &
                    clean_up=clean==1,status=status)
                if (clean == 2) removed = real_fit%clean_up(ctol,status)
                s = cmplx(real_fit%support_points(),0.0_dp,dp)
                v = cmplx(real_fit%support_values(),0.0_dp,dp)
                w = cmplx(real_fit%weights(),0.0_dp,dp)
                p = real_fit%poles(ps)
                z = real_fit%roots(zs)
                r = real_fit%residues(rs)
            else
                complex_fit = aaa(x,f,rtol=0.0_dp,max_terms=n,clean_up=clean==1,status=status)
                if (clean == 2) removed = complex_fit%clean_up(ctol,status)
                s = complex_fit%support_points()
                v = complex_fit%support_values()
                w = complex_fit%weights()
                p = complex_fit%poles(ps)
                z = complex_fit%roots(zs)
                r = complex_fit%residues(rs)
            end if
            print '(9(i0,1x))',status,removed,size(s),ps,zs,rs,size(p),size(z),size(r)
            do i = 1, size(s)
                print '(6(es26.17e3,1x))',real(s(i),dp),aimag(s(i)),real(v(i),dp), &
                    aimag(v(i)),real(w(i),dp),aimag(w(i))
            end do
            do i = 1, size(p)
                print '(2(es26.17e3,1x))',real(p(i),dp),aimag(p(i))
            end do
            do i = 1, size(z)
                print '(2(es26.17e3,1x))',real(z(i),dp),aimag(z(i))
            end do
            do i = 1, size(r)
                print '(2(es26.17e3,1x))',real(r(i),dp),aimag(r(i))
            end do
        end block
    end do
end program aaa_large_driver
