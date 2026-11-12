! SPDX-License-Identifier: MIT
! Stream-mode driver for SciPy generalized eig/eigvals comparisons.
program driver_generalized_eigen
    use scifort_kinds, only : dp
    use scifort_linalg_generalized_eigen, only : eig, eigvals
    implicit none
    real(dp), allocatable :: ar(:, :), br(:, :)
    complex(dp), allocatable :: ac(:, :), bc(:, :)
    complex(dp), allocatable :: w(:), wval(:), alpha(:), beta(:), av(:), bv(:)
    complex(dp), allocatable :: vl(:, :), vr(:, :)
    real(dp) :: ra, ia, rb, ib
    integer :: nc, k, n, kind, i, j, status, statusval

    read (*,*) nc
    do k=1,nc
        read (*,*) kind, n
        allocate(ar(n,n),br(n,n),ac(n,n),bc(n,n))
        allocate(w(n),wval(n),alpha(n),beta(n),av(n),bv(n),vl(n,n),vr(n,n))
        do j=1,n
            do i=1,n
                read (*,*) ra,ia,rb,ib
                ac(i,j) = cmplx(ra,ia,dp)
                bc(i,j) = cmplx(rb,ib,dp)
            end do
        end do
        ar = real(ac,dp)
        br = real(bc,dp)
        select case(kind)
        case(1)
            call eig(ar,br,w,status,vr=vr,vl=vl,alpha=alpha,beta=beta)
            call eigvals(ar,br,wval,statusval,alpha=av,beta=bv)
        case(2)
            call eig(ar,bc,w,status,vr=vr,vl=vl,alpha=alpha,beta=beta)
            call eigvals(ar,bc,wval,statusval,alpha=av,beta=bv)
        case(3)
            call eig(ac,br,w,status,vr=vr,vl=vl,alpha=alpha,beta=beta)
            call eigvals(ac,br,wval,statusval,alpha=av,beta=bv)
        case(4)
            call eig(ac,bc,w,status,vr=vr,vl=vl,alpha=alpha,beta=beta)
            call eigvals(ac,bc,wval,statusval,alpha=av,beta=bv)
        case default
            error stop 'invalid data kind'
        end select
        print '(i0,1x,i0)',status,statusval
        do i=1,n
            print '(12(es26.17e3,1x))',real(w(i),dp),aimag(w(i)), &
                real(alpha(i),dp),aimag(alpha(i)),real(beta(i),dp),aimag(beta(i)), &
                real(wval(i),dp),aimag(wval(i)),real(av(i),dp),aimag(av(i)), &
                real(bv(i),dp),aimag(bv(i))
        end do
        do j=1,n
            do i=1,n
                print '(4(es26.17e3,1x))',real(vr(i,j),dp),aimag(vr(i,j)), &
                    real(vl(i,j),dp),aimag(vl(i,j))
            end do
        end do
        deallocate(ar,br,ac,bc,w,wval,alpha,beta,av,bv,vl,vr)
    end do
end program driver_generalized_eigen
