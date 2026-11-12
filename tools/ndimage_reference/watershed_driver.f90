! SPDX-License-Identifier: MIT
program watershed_driver
    use scifort_ndimage, only : watershed_ift_nd
    implicit none
    integer :: a(6,7),m(6,7),res(6,7)
    integer :: a3(4,3,3),m3(4,3,3),r3(4,3,3)
    integer :: i,j,k,t,stat,unit
    logical :: st(3,3)
    open(newunit=unit,file='watershed_results.bin',form='unformatted', &
        access='stream',status='replace')
    do t=1,24
        do j=1,7
            do i=1,6
                a(i,j)=modulo(11*i*i+17*i*j+13*j*j+9*t*i+2*t*j,13+t*3)
                m(i,j)=0
            end do
        end do
        m(2,2)=1
        m(5,6)=2
        if (modulo(t,3)==0) m(1,7)=-3
        if (modulo(t,4)==0) m(4,3)=4
        select case(modulo(t,4))
        case(0)
            st=.true.
        case(1)
            st=.false.
            st(2,:)=.true.
            st(:,2)=.true.
        case(2)
            st=.false.
            st(1,1)=.true.
            st(2,2)=.true.
            st(3,3)=.true.
        case(3)
            st=.false.
            st(1,2)=.true.
            st(2,2)=.true.
            st(3,2)=.true.
        end select
        call watershed_ift_nd(a,m,res,stat,structure=st)
        if (stat/=0) error stop 'watershed 2d'
        write(unit) res
        call watershed_ift_nd(a,m,res,stat,structure=st,legacy_linear_bounds=.true.)
        if (stat/=0) error stop 'watershed 2d compat'
        write(unit) res
    end do
    do t=1,12
        do k=1,3
            do j=1,3
                do i=1,4
                    a3(i,j,k)=modulo(5*i*i+9*i*j+7*k+11*j*k+t*i,19+t)
                    m3(i,j,k)=0
                end do
            end do
        end do
        m3(1,1,1)=1
        m3(4,3,3)=-2
        if (modulo(t,2)==0) m3(2,2,2)=3
        call watershed_ift_nd(a3,m3,r3,stat)
        if (stat/=0) error stop 'watershed 3d'
        write(unit) r3
        call watershed_ift_nd(a3,m3,r3,stat,legacy_linear_bounds=.true.)
        if (stat/=0) error stop 'watershed 3d compat'
        write(unit) r3
    end do
    close(unit)
    print '(a)', 'PASS: watershed reference driver'
end program watershed_driver
