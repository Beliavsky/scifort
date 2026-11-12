! SPDX-License-Identifier: MIT
! Independent deterministic data for SciPy 1.17 mixed-type convolution checks.
! CASE family (1=conv,2=corr,3=fftconv,4=oaconv) order mode algorithm n m nout
! mode = full/same/valid; algorithm = direct/fft/auto for families 1,2.
program signal_mixed_convolution_driver
    use scifort_kinds, only : dp
    use scifort_signal_convolve, only : convolve, correlate, convolution_status_success
    use scifort_signal_fftconvolve, only : fftconvolve, oaconvolve
    implicit none
    integer, parameter :: nlens(5) = [1,2,11,179,701], mlens(4) = [1,2,19,101]
    character(len=6), parameter :: methods(3) = ['direct','fft   ','auto  ']
    character(len=5), parameter :: modes(3) = ['full ','same ','valid']
    real(dp), allocatable :: a(:), b(:)
    complex(dp), allocatable :: za(:), zb(:), output(:)
    integer :: i, j, k, l, family, order, alg, method_count, state, t

    do i=1,size(nlens)
        allocate(a(nlens(i)),za(nlens(i)))
        do t=1,nlens(i)
            a(t) = sin(0.137_dp*real(t,dp)) + 0.25_dp*cos(0.031_dp*real(t,dp))
            za(t) = cmplx(a(t),0.23_dp*cos(0.191_dp*real(t,dp)),dp)
        end do
        do j=1,size(mlens)
            allocate(b(mlens(j)),zb(mlens(j)))
            do t=1,mlens(j)
                b(t) = cos(0.221_dp*real(t,dp))-0.17_dp*sin(0.052_dp*real(t,dp))
                zb(t) = cmplx(b(t),0.29_dp*sin(0.109_dp*real(t,dp)),dp)
            end do
            do family=1,4
                method_count = 1
                if (family<=2) method_count = 3
                do order=1,2
                    do k=1,3
                        do alg=1,method_count
                            select case(family)
                            case(1)
                                if (order==1) then
                                    output=convolve(a,zb,trim(modes(k)),trim(methods(alg)),state)
                                else
                                    output=convolve(za,b,trim(modes(k)),trim(methods(alg)),state)
                                end if
                            case(2)
                                if (order==1) then
                                    output=correlate(a,zb,trim(modes(k)),trim(methods(alg)),state)
                                else
                                    output=correlate(za,b,trim(modes(k)),trim(methods(alg)),state)
                                end if
                            case(3)
                                if (order==1) then
                                    output=fftconvolve(a,zb,trim(modes(k)),state)
                                else
                                    output=fftconvolve(za,b,trim(modes(k)),state)
                                end if
                            case(4)
                                if (order==1) then
                                    output=oaconvolve(a,zb,trim(modes(k)),state)
                                else
                                    output=oaconvolve(za,b,trim(modes(k)),state)
                                end if
                            end select
                            if (state /= convolution_status_success) error stop 'mixed convolution reference status'
                            write(*,'(a,7(1x,i0))') 'CASE',family,order,k,alg,nlens(i),mlens(j),size(output)
                            do t=1,size(output)
                                write(*,'(es25.16e3,1x,es25.16e3)') real(output(t),dp),aimag(output(t))
                            end do
                        end do
                    end do
                end do
            end do
            deallocate(b,zb)
        end do
        deallocate(a,za)
    end do
end program signal_mixed_convolution_driver
