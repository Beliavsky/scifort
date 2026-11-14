! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program generate_knots_driver
    use scifort_kinds, only: dp
    use scifort_interpolate, only: generate_knots,knot_generator_type
    implicit none
    type(knot_generator_type) :: gen
    real(dp), allocatable :: x(:),y(:,:),w(:),t(:)
    real(dp) :: s,xb,xe
    integer :: m,ndim,k,periodic,nest,i,code,weighted
    logical :: done
    character(10) :: bc
    read(*,*) m,ndim,k,periodic,nest,weighted,s,xb,xe
    allocate(x(m),y(m,ndim),w(m))
    read(*,*) x
    do i=1,m
        read(*,*) y(i,:)
    end do
    read(*,*) w
    bc='not-a-knot'
    if (periodic==1) bc='periodic'
    if (s==0.0_dp) then
        gen=generate_knots(x,y,k=k,xb=xb,xe=xe,bc_type=bc,status=code)
    else if (nest==0) then
        if (weighted==1) then
            gen=generate_knots(x,y,k=k,w=w,s=s,xb=xb,xe=xe,bc_type=bc,status=code)
        else
            gen=generate_knots(x,y,k=k,s=s,xb=xb,xe=xe,bc_type=bc,status=code)
        end if
    else
        gen=generate_knots(x,y,k=k,w=w,s=s,xb=xb,xe=xe,nest=nest,bc_type=bc,status=code)
    end if
    write(*,'(i0)') code
    if (code/=0) stop
    do
        call gen%next(t,done,code)
        write(*,'(i0,1x,i0)') code,size(t)
        if (code/=0 .or. done) exit
        write(*,'(*(es25.16e3,1x))') t
    end do
end program generate_knots_driver
