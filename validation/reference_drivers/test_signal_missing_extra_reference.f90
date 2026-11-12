! SPDX-License-Identifier: MIT
! Data-driven SciPy reference comparison driver; see compare_signal_missing_extra_scipy.py.
program test_signal_missing_extra_reference
    use scifort_kinds, only : dp
    use scifort_signal, only : convolve, correlate, normalize, normalized_tf_result, &
        whittaker_henderson, whittaker_result
    implicit none
    character(10) :: tag, mode, method
    integer :: n, m, p, nr, nc, na, status, ios, i
    real(dp) :: lam
    real(dp), allocatable :: x(:), y(:), w(:), result(:), den(:), num(:), rows(:,:)
    complex(dp), allocatable :: cx(:), cy(:), cz(:)
    type(normalized_tf_result) :: tf
    type(whittaker_result) :: smoothed

    do
        read(*,*,iostat=ios) tag
        if (ios /= 0) exit
        select case(trim(tag))
        case ('WH')
            read(*,*) n,p,lam
            allocate(x(n),w(n))
            read(*,*) x
            read(*,*) w
            smoothed = whittaker_henderson(x,lam,order=p,weights=w)
            write(*,'(i0,1x,i0)') smoothed%status,size(smoothed%x)
            if (smoothed%status==0) write(*,'(*(es24.16,1x))') smoothed%x
            deallocate(x,w)
        case ('NOR')
            read(*,*) n,m
            allocate(num(n),den(m))
            read(*,*) num
            read(*,*) den
            tf = normalize(num,den)
            write(*,'(4(i0,1x))') tf%status,tf%trimmed,size(tf%b,2),size(tf%a)
            if (tf%status==0) then
                write(*,'(*(es24.16,1x))') real(tf%b(1,:),dp)
                write(*,'(*(es24.16,1x))') real(tf%a,dp)
            end if
            deallocate(num,den)
        case ('NOM')
            read(*,*) nr,nc,na
            allocate(rows(nr,nc),den(na))
            do i = 1,nr
                read(*,*) rows(i,:)
            end do
            read(*,*) den
            tf = normalize(rows,den)
            write(*,'(5(i0,1x))') tf%status,tf%trimmed,size(tf%b,1),size(tf%b,2),size(tf%a)
            if (tf%status==0) then
                do i = 1,size(tf%b,1)
                    write(*,'(*(es24.16,1x))') real(tf%b(i,:),dp)
                end do
                write(*,'(*(es24.16,1x))') real(tf%a,dp)
            end if
            deallocate(rows,den)
        case ('CONR','CORR')
            read(*,*) n,m,mode,method
            allocate(x(n),y(m))
            read(*,*) x
            read(*,*) y
            if (tag=='CONR') then
                result = convolve(x,y,trim(mode),trim(method),status)
            else
                result = correlate(x,y,trim(mode),trim(method),status)
            end if
            write(*,'(2(i0,1x))') status,size(result)
            if (status==0) write(*,'(*(es24.16,1x))') result
            deallocate(x,y,result)
        case ('CONC','CORC')
            read(*,*) n,m,mode,method
            allocate(cx(n),cy(m))
            ! Complex list-directed input accepts pairs (real,imag).
            read(*,*) cx
            read(*,*) cy
            if (tag=='CONC') then
                cz = convolve(cx,cy,trim(mode),trim(method),status)
            else
                cz = correlate(cx,cy,trim(mode),trim(method),status)
            end if
            write(*,'(2(i0,1x))') status,size(cz)
            if (status==0) then
                write(*,'(*(es24.16,1x))') real(cz,dp)
                write(*,'(*(es24.16,1x))') aimag(cz)
            end if
            deallocate(cx,cy,cz)
        case default
            error stop 'invalid reference input record'
        end select
    end do
end program test_signal_missing_extra_reference
