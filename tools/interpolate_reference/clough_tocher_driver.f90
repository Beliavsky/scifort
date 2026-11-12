! SPDX-License-Identifier: MIT
! Batched SciPy numerical comparison driver. Standard input contains cases:
! n sites, m channels, q queries, rescale (0/1); then n sites, n complex
! multichannel value rows, and q query coordinates. Terminate with n=0.
program clough_tocher_reference_driver
    use scifort_kinds, only: dp
    use scifort_interpolate_clough_tocher, only: clough_tocher_2d_interpolator, &
        complex_clough_tocher_2d_interpolator_type
    implicit none
    real(dp), allocatable :: sites(:,:), queries(:,:)
    complex(dp), allocatable :: values(:,:), predicted(:,:)
    type(complex_clough_tocher_2d_interpolator_type) :: interp
    real(dp) :: u,v
    integer :: n,m,q,scaled,i,j,stat
    do
        read(*,*,iostat=stat) n,m,q,scaled
        if (stat/=0 .or. n==0) exit
        allocate(sites(n,2),queries(q,2),values(n,m))
        do i=1,n
            read(*,*) sites(i,:)
        end do
        do i=1,n
            do j=1,m
                read(*,*) u,v
                values(i,j)=cmplx(u,v,dp)
            end do
        end do
        do i=1,q
            read(*,*) queries(i,:)
        end do
        interp=clough_tocher_2d_interpolator(sites,values,rescale=scaled==1, &
            tol=1.0e-11_dp,maxiter=1600,status=stat)
        write(*,'(i0,1x,l1)') stat,interp%converged()
        if (stat==0) then
            predicted=interp%evaluate(queries,status=stat)
            do i=1,q
                do j=1,m
                    write(*,'(2(es25.16e3,1x))') real(predicted(i,j),dp),aimag(predicted(i,j))
                end do
            end do
        end if
        deallocate(sites,queries,values)
    end do
end program clough_tocher_reference_driver
