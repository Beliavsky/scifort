! SPDX-License-Identifier: MIT
! Text protocol for independent NumPy comparison of prepared harmonic fits.
program reference_signal_harmonic_plan
    use scifort_kinds, only : dp
    use scifort_signal_irregular, only : harmonic_fit_plan, harmonic_fit_result, &
        lombscargle_harmonic_prepare, lombscargle_harmonic_fit_batch
    implicit none

    integer :: ncases, icase, n, m, nch, nh, weighted, floating, center, covariance
    integer :: j, k, h, a, b, p
    real(dp), allocatable :: x(:), w(:), f(:), y(:,:)
    type(harmonic_fit_plan) :: plan
    type(harmonic_fit_result), allocatable :: fits(:)

    read (*,*) ncases
    do icase = 1, ncases
        read (*,*) n, m, nch, nh, weighted, floating, center, covariance
        allocate(x(n),w(n),f(m),y(n,nch))
        read (*,*) x
        read (*,*) f
        read (*,*) w
        do j = 1, nch
            read (*,*) y(:,j)
        end do
        if (weighted == 1) then
            plan = lombscargle_harmonic_prepare(x,f,nh,weights=w,floating_mean=floating == 1)
        else
            plan = lombscargle_harmonic_prepare(x,f,nh,floating_mean=floating == 1)
        end if
        fits = lombscargle_harmonic_fit_batch(plan,y,precenter=center == 1, &
            estimate_covariance=covariance == 1)
        write (*,'(i0,1x,i0)') plan%status,size(fits)
        p = 1+2*nh
        do j = 1, nch
            write (*,'(i0,1x,i0)') fits(j)%status,fits(j)%degrees_of_freedom
            do k = 1, m
                write (*,'(i0)') merge(1,0,fits(j)%valid(k))
                if (.not. fits(j)%valid(k)) cycle
                write (*,'(es26.17e3)') fits(j)%offset(k)
                do h = 1, nh
                    write (*,'(es26.17e3)') fits(j)%cosine(h,k)
                    write (*,'(es26.17e3)') fits(j)%sine(h,k)
                end do
                write (*,'(es26.17e3)') fits(j)%power(k)
                write (*,'(es26.17e3)') fits(j)%normalized_power(k)
                write (*,'(es26.17e3)') fits(j)%rms(k)
                if (covariance /= 1) cycle
                write (*,'(es26.17e3)') fits(j)%residual_variance(k)
                do a = 1, p
                    write (*,'(es26.17e3)') fits(j)%standard_error(a,k)
                end do
                do b = 1, p
                    do a = 1, p
                        write (*,'(es26.17e3)') fits(j)%covariance(a,b,k)
                    end do
                end do
            end do
        end do
        deallocate(x,w,f,y)
    end do
end program reference_signal_harmonic_plan
