! Read batches of dense linear programs from stdin for independent SciPy comparisons.
! Compile manually; this is not an FPM auto-test executable.
program linprog_reference_driver
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
    use scifort_kinds, only : dp
    use scifort_optimize_linprog, only : linprog, linprog_result
    implicit none
    type(linprog_result) :: res
    real(dp), allocatable :: c(:), lo(:), hi(:), aub(:,:), bub(:), aeq(:,:), beq(:)
    integer :: n, nub, neq, i, status
    do
        read(*,*,iostat=status) n,nub,neq
        if (status/=0 .or. n<0) exit
        allocate(c(n),lo(n),hi(n),aub(nub,n),bub(nub),aeq(neq,n),beq(neq))
        read(*,*) c
        read(*,*) lo
        read(*,*) hi
        do i=1,nub
            read(*,*) aub(i,:)
        end do
        if (nub>0) read(*,*) bub
        do i=1,neq
            read(*,*) aeq(i,:)
        end do
        if (neq>0) read(*,*) beq
        res=linprog(c,a_ub=aub,b_ub=bub,a_eq=aeq,b_eq=beq,lower=lo,upper=hi)
        write(*,'(i0,1x,es26.17e3,1x,i0)') res%status,res%fun,res%nit
        if (allocated(res%x)) then
            write(*,'(*(es26.17e3,1x))') res%x
        else
            write(*,'(a)') 'no-solution'
        end if
        deallocate(c,lo,hi,aub,bub,aeq,beq)
    end do
end program linprog_reference_driver
