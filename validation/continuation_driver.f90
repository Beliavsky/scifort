! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
! Validation only: link the exact, separately cached upstream module for comparison.
program continuation_driver
    use scifort_kinds, only: dp
    use scifort_interpolate, only: splrep,splprep,fitpack_tck_type, &
        fitpack_parametric_tck_type,fitpack_continuation_type
    use fitpack_core, only: curfit,percur,parcur,clocur
    implicit none
    character(20) :: argument
    logical :: reference,periodic
    integer :: m,ndim,k,per,steps,io,i,j,step,mode,nest,n,lwrk,ier,original,nc,code,active
    real(dp) :: target,xscale,wscale,ub,ue,fp
    real(dp), allocatable :: u(:),x(:,:),w(:),targets(:),t(:),c(:),wrk(:),xx(:,:),ww(:),cc(:,:)
    integer, allocatable :: iw(:)
    type(fitpack_continuation_type) :: state
    type(fitpack_tck_type) :: scalar
    type(fitpack_parametric_tck_type) :: curve
    call get_command_argument(1,argument)
    reference=trim(argument)=='reference'
    do
        read(*,*,iostat=io) m,ndim,k,per,steps
        if (io/=0) exit
        block
            allocate(u(m),x(m,ndim),w(m),targets(steps))
            read(*,*) u
            do i=1,m
                read(*,*) x(i,:)
            end do
            read(*,*) w
            read(*,*) targets
            periodic=per/=0
            active=m-per
            nest=m+k+1
            if (periodic) nest=m+2*k
            lwrk=m*(k+1)+nest*(20+2*ndim+5*k)
            allocate(t(nest),c(nest*ndim),wrk(lwrk),iw(nest),xx(ndim,m),ww(m))
            xscale=maxval(abs(x(:active,:)))
            if (xscale==0.0_dp) xscale=1.0_dp
            wscale=maxval(w(:active))
            xx=transpose(x/xscale)
            ww=w/wscale
            if (periodic) then
                xx(:,m)=xx(:,1)
                ww(m)=ww(1)
            end if
            t=0.0_dp
            c=0.0_dp
            wrk=0.0_dp
            iw=0
            n=2*(k+1)
            call state%reset()
            do step=1,steps
                mode=1
                if (step==1) mode=0
                if (reference) then
                    target=0.0_dp
                    if (targets(step)>0.0_dp) target=exp(log(targets(step))-2*log(xscale)-2*log(wscale))
                    call upstream(mode)
                    original=ier
                    ! Record the unmodified upstream failure before independent fresh restart.
                    if (ier>0 .and. periodic .and. ndim>1 .and. step>1) then
                        if (targets(step-1)==0.0_dp .and. targets(step)>0.0_dp) then
                            n=2*(k+1)
                            call upstream(0)
                        end if
                    end if
                    code=0
                    if (ier>0) code=3
                    print '(2(i0,1x))',code,original
                    if (code/=0) exit
                    nc=n-k-1
                    allocate(cc(nc,ndim))
                    do j=1,ndim
                        cc(:,j)=c((j-1)*n+1:(j-1)*n+nc)*xscale
                    end do
                else
                    if (ndim==1) then
                        scalar=splrep(u,x(:,1),w=w,k=k,s=targets(step),per=periodic, &
                            task=mode,state=state,status=code)
                        print '(2(i0,1x))',code,code
                        if (code/=0) exit
                        n=size(scalar%t)
                        t(:n)=scalar%t
                        cc=reshape(scalar%c(:n-k-1),[n-k-1,1])
                    else
                        curve=splprep(x,u=u,w=w,k=k,s=targets(step),per=periodic, &
                            task=mode,state=state,status=code)
                        print '(2(i0,1x))',code,code
                        if (code/=0) exit
                        n=size(curve%t)
                        t(:n)=curve%t
                        cc=curve%c
                    end if
                end if
                print '(i0)',n
                print '(*(es26.17e3,1x))',t(:n)
                do j=1,ndim
                    print '(*(es26.17e3,1x))',cc(:,j)
                end do
                deallocate(cc)
            end do
            deallocate(u,x,w,targets,t,c,wrk,iw,xx,ww)
        end block
    end do
contains
    subroutine upstream(task)
        implicit none
        integer, intent(in) :: task !! Zero starts fresh; one uses retained upstream work arrays.
        ub=0.0_dp
        ue=1.0_dp
        fp=0.0_dp
        ier=0
        if (ndim==1) then
            if (periodic) then
                call percur(task,m,u,xx(1,:),ww,k,target,nest,n,t,c,fp,wrk,lwrk,iw,ier)
            else
                call curfit(task,m,u,xx(1,:),ww,ub,ue,k,target,nest,n,t,c,fp,wrk,lwrk,iw,ier)
            end if
        else
            if (periodic) then
                call clocur(task,1,ndim,m,u,m*ndim,xx,ww,k,target,nest,n,t, &
                    nest*ndim,c,fp,wrk,lwrk,iw,ier)
            else
                call parcur(task,1,ndim,m,u,m*ndim,xx,ww,ub,ue,k,target,nest,n,t, &
                    nest*ndim,c,fp,wrk,lwrk,iw,ier)
            end if
        end if
    end subroutine upstream
end program continuation_driver
