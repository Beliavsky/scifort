! SPDX-License-Identifier: MIT
! Compare ordered generalized complex Schur forms with independent SciPy data.
program driver_ordqz
    use scifort_kinds, only : dp
    use scifort_linalg_ordqz, only : ordqz, ordqz_factors
    implicit none
    complex(dp), allocatable :: a(:,:), b(:,:), aa(:,:), bb(:,:), q(:,:), z(:,:)
    complex(dp), allocatable :: alpha(:), beta(:)
    real(dp), allocatable :: ar(:,:), br(:,:)
    character(len=4), parameter :: sorts(5) = ['lhp ', 'rhp ', 'iuc ', 'ouc ', 'none']
    integer :: n, family, supplied, sort_index, status, sdim, i, j, ios
    real(dp) :: xr, xi, yr, yi

    do
        read(*,*,iostat=ios) n, family, supplied, sort_index
        if (ios /= 0) exit
        allocate(a(n,n),b(n,n),aa(n,n),bb(n,n),q(n,n),z(n,n))
        allocate(ar(n,n),br(n,n),alpha(n),beta(n))
        do j=1,n
            do i=1,n
                read(*,*) xr,xi,yr,yi
                a(i,j) = cmplx(xr,xi,dp)
                b(i,j) = cmplx(yr,yi,dp)
                ar(i,j) = xr
                br(i,j) = yr
            end do
        end do
        if (supplied == 1) then
            do j=1,n
                do i=1,n
                    read(*,*) xr,xi,yr,yi
                    q(i,j) = cmplx(xr,xi,dp)
                    z(i,j) = cmplx(yr,yi,dp)
                end do
            end do
            aa = a
            bb = b
            call ordqz_factors(aa,bb,q,z,status,sort=trim(sorts(sort_index)),sdim=sdim)
        else
            select case(family)
            case(1)
                call ordqz(ar,br,aa,bb,q,z,status,alpha,beta, &
                    sort=trim(sorts(sort_index)),sdim=sdim)
            case(2)
                call ordqz(ar,b,aa,bb,q,z,status,alpha,beta, &
                    sort=trim(sorts(sort_index)),sdim=sdim)
            case(3)
                call ordqz(a,br,aa,bb,q,z,status,alpha,beta, &
                    sort=trim(sorts(sort_index)),sdim=sdim)
            case default
                call ordqz(a,b,aa,bb,q,z,status,alpha,beta, &
                    sort=trim(sorts(sort_index)),sdim=sdim)
            end select
        end if
        if (supplied == 1) then
            do i=1,n
                alpha(i) = aa(i,i)
                beta(i) = bb(i,i)
            end do
        end if
        print '(2(i0,1x))', status,sdim
        do j=1,n
            do i=1,n
                print '(8(es26.17e3,1x))', real(aa(i,j),dp),aimag(aa(i,j)), &
                    real(bb(i,j),dp),aimag(bb(i,j)), &
                    real(q(i,j),dp),aimag(q(i,j)),real(z(i,j),dp),aimag(z(i,j))
            end do
        end do
        do i=1,n
            print '(4(es26.17e3,1x))', real(alpha(i),dp),aimag(alpha(i)), &
                real(beta(i),dp),aimag(beta(i))
        end do
        deallocate(a,b,aa,bb,q,z,ar,br,alpha,beta)
    end do
end program driver_ordqz
