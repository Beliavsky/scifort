! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program real_qz_driver
    use scifort_kinds, only: dp
    use scifort_linalg_qz, only: qz
    use scifort_linalg_ordqz, only: ordqz,ordqz_factors
    implicit none
    integer :: n,mode,io,i,j,status,sdim
    character(4), parameter :: sorts(5)=[character(4)::'none','lhp','rhp','iuc','ouc']
    do
        read(*,*,iostat=io) n,mode
        if (io/=0) exit
        block
            real(dp), allocatable :: a(:,:),b(:,:),aa(:,:),bb(:,:),q(:,:),z(:,:),beta(:)
            complex(dp), allocatable :: alpha(:)
            allocate(a(n,n),b(n,n),aa(n,n),bb(n,n),q(n,n),z(n,n),alpha(n),beta(n))
            do i=1,n
                read(*,*) a(i,:)
            end do
            do i=1,n
                read(*,*) b(i,:)
            end do
            sdim=0
            if (mode==0) then
                call qz(a,b,aa,bb,q,z,status,alpha,beta)
            else if (mode==5) then
                call ordqz(a,b,aa,bb,q,z,status,alpha,beta,sdim=sdim,selector=positive_imag)
            else if (mode==6) then
                call qz(a,b,aa,bb,q,z,status)
                if (status==0) call ordqz_factors(aa,bb,q,z,status,sort='lhp',sdim=sdim)
            else
                call ordqz(a,b,aa,bb,q,z,status,alpha,beta,sort=trim(sorts(mode+1)),sdim=sdim)
            end if
            print '(2(i0,1x))',status,sdim
            if (status/=0) cycle
            do j=1,4
                do i=1,n
                    select case(j)
                    case(1)
                        print '(*(es26.17e3,1x))',aa(i,:)
                    case(2)
                        print '(*(es26.17e3,1x))',bb(i,:)
                    case(3)
                        print '(*(es26.17e3,1x))',q(i,:)
                    case(4)
                        print '(*(es26.17e3,1x))',z(i,:)
                    end select
                end do
            end do
            if (mode/=6) then
                print '(*(es26.17e3,1x))',real(alpha,dp)
                print '(*(es26.17e3,1x))',aimag(alpha)
                print '(*(es26.17e3,1x))',beta
            end if
        end block
    end do
contains
    logical function positive_imag(alpha,beta) result(selected)
        implicit none
        complex(dp), intent(in) :: alpha !! Homogeneous numerator supplied by the decomposition.
        complex(dp), intent(in) :: beta !! Homogeneous denominator supplied by the decomposition.
        selected=aimag(alpha)>0.0_dp
    end function positive_imag
end program real_qz_driver
