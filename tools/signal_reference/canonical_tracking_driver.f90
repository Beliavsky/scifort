! SPDX-License-Identifier: MIT
! Read independent spectral canonical vectors for matching/phase reference tests.
program canonical_tracking_driver
    use scifort_kinds, only : dp
    use scifort_signal_canonical_coherence, only : canonical_coherence_result
    use scifort_signal_canonical_tracking, only : canonical_mode_tracks_result, track_canonical_modes
    implicit none
    type(canonical_coherence_result) :: input
    type(canonical_mode_tracks_result) :: output
    integer :: ncase, c, f, i, k, nf, na, nb, nm
    real(dp) :: re, im, minimum

    read(*,*) ncase
    do c = 1, ncase
        read(*,*) nf, na, nb, nm, minimum
        allocate(input%frequencies(nf),input%valid(nf),input%directions_stable(nf,nm), &
            input%squared_modes(nf,nm),input%coefficients_a(nf,na,nm), &
            input%coefficients_b(nf,nb,nm))
        input%status = 0
        do f = 1, nf
            read(*,*) input%frequencies(f),input%valid(f)
            read(*,*) input%squared_modes(f,:)
            read(*,*) input%directions_stable(f,:)
            do k = 1, nm
                do i = 1, na
                    read(*,*) re,im
                    input%coefficients_a(f,i,k) = cmplx(re,im,dp)
                end do
                do i = 1, nb
                    read(*,*) re,im
                    input%coefficients_b(f,i,k) = cmplx(re,im,dp)
                end do
            end do
        end do
        output = track_canonical_modes(input,min_overlap=minimum)
        write(*,'(a,1x,i0,1x,i0)') 'CASE',c,output%status
        if (output%status == 0) then
            do f = 1, nf
                do k = 1, nm
                    write(*,'(i0,1x,i0,1x,i0,1x,l1,1x,es25.16e3,1x,es25.16e3)', &
                        advance='no') f,k,output%source_mode(f,k), &
                        output%matched(f,k),output%overlap(f,k),output%squared_modes(f,k)
                    do i = 1, na
                        write(*,'(1x,es25.16e3,1x,es25.16e3)',advance='no') &
                            real(output%coefficients_a(f,i,k),dp),aimag(output%coefficients_a(f,i,k))
                    end do
                    do i = 1, nb
                        write(*,'(1x,es25.16e3,1x,es25.16e3)',advance='no') &
                            real(output%coefficients_b(f,i,k),dp),aimag(output%coefficients_b(f,i,k))
                    end do
                    write(*,*)
                end do
            end do
        end if
        deallocate(input%frequencies,input%valid,input%directions_stable, &
            input%squared_modes,input%coefficients_a,input%coefficients_b)
    end do
end program canonical_tracking_driver
