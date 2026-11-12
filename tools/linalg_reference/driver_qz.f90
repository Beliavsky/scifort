! SPDX-License-Identifier: MIT
! Text protocol for independent Python/SciPy generalized Schur comparisons.
program driver_qz
    use scifort_kinds, only : dp
    use scifort_linalg_qz, only : qz
    implicit none
    real(dp), allocatable :: ar(:,:), br(:,:)
    complex(dp), allocatable :: a(:,:), b(:,:), aa(:,:), bb(:,:), q(:,:), z(:,:)
    complex(dp), allocatable :: alpha(:), beta(:)
    real(dp) :: x, y, u, v
    integer :: total, k, n, family, i, j, status
    character(len=20) :: method

    method = 'auto'
    if (command_argument_count() >= 1) call get_command_argument(1,method)
    read(*,*) total
    do k=1,total
        read(*,*) family, n
        allocate(ar(n,n),br(n,n),a(n,n),b(n,n),aa(n,n),bb(n,n), &
            q(n,n),z(n,n),alpha(n),beta(n))
        do j=1,n
            do i=1,n
                read(*,*) x,y,u,v
                a(i,j) = cmplx(x,y,dp)
                b(i,j) = cmplx(u,v,dp)
                ar(i,j) = x
                br(i,j) = u
            end do
        end do
        select case(family)
        case(1)
            call qz(ar,br,aa,bb,q,z,status,alpha,beta,method=trim(method))
        case(2)
            call qz(ar,b,aa,bb,q,z,status,alpha,beta,method=trim(method))
        case(3)
            call qz(a,br,aa,bb,q,z,status,alpha,beta,method=trim(method))
        case default
            call qz(a,b,aa,bb,q,z,status,alpha,beta,method=trim(method))
        end select
        print '(i0)', status
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
        deallocate(ar,br,a,b,aa,bb,q,z,alpha,beta)
    end do
end program driver_qz
