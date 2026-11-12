! SPDX-License-Identifier: MIT
! Driven by test/compare_signal_cosine_dpss.py, output one line per case.
program signal_cosine_dpss_reference_driver
    use scifort_kinds, only : dp
    use scifort_signal_spectral, only : get_window, general_cosine
    implicit none
    integer :: kind, n, sy, k, code, ios, i
    real(dp) :: nw
    real(dp), allocatable :: coeff(:), w(:)

    do
        read(*,*,iostat=ios) kind, n, sy, k, nw
        if (ios /= 0) exit
        allocate(coeff(k))
        if (k > 0) read(*,*,iostat=ios) coeff
        if (ios /= 0) error stop 'missing coefficients'
        select case(kind)
        case(1)
            w = general_cosine(n,coeff,sym=sy==1,status=code)
        case(2)
            w = get_window('general_cosine',n,coefficients=coeff,fftbins=sy==0,status=code)
        case(3)
            w = get_window('dpss',n,parameter=nw,fftbins=sy==0,status=code)
        case(4)
            w = get_window('slepian',n,parameter=nw,fftbins=sy==0,status=code)
        case default
            error stop 'unknown family'
        end select
        write(*,'(i0,1x,i0,*(1x,es26.17e3))') code, size(w), w
        deallocate(coeff)
    end do
end program signal_cosine_dpss_reference_driver
