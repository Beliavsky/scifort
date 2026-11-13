! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
program test_fft_nd
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_nan
    use scifort_fft
    implicit none
    complex(dp) :: x(3,5), y(3,5), z(3,5), expected(3,5), half(3,3), padded(4,7)
    complex(dp) :: reordered(2,5), gap(6,5), single(1,1), line(3), line_out(3), padded_back(4,7)
    complex(dp) :: volume(2,3,4), vout(2,3,4), vback(2,3,4)
    complex(dp) :: high(2,1,1,1,1,1,1,1,1,1,1,1,1,1,2)
    complex(dp) :: hout(2,1,1,1,1,1,1,1,1,1,1,1,1,1,2)
    complex(dp) :: hback(2,1,1,1,1,1,1,1,1,1,1,1,1,1,2)
    complex(dp) :: empty_result(0,5), empty_half(0,3)
    real(dp) :: real_x(3,5), real_y(3,5), real_z(3,5), small(1,1), small_out(1,1), empty(0,5)
    real(dp) :: nan, inf, angle, factor
    real(dp), parameter :: pi = acos(-1.0_dp)
    character(8), parameter :: norms(3) = [character(8) :: 'backward','forward','ortho']
    integer :: i, j, k, l, mode, kind, code, checks, failures
    checks = 0
    failures = 0
    nan = ieee_value(0.0_dp,ieee_quiet_nan)
    inf = ieee_value(0.0_dp,ieee_positive_inf)
    do j = 1, 5
        do i = 1, 3
            x(i,j) = cmplx(sin(real(i+2*j,dp)),cos(real(i-j,dp)),dp)
        end do
    end do
    real_x = real(x,dp)
    do l = 1, 5
        do k = 1, 3
            expected(k,l) = (0.0_dp,0.0_dp)
            do j = 1, 5
                do i = 1, 3
                    angle = -2*pi*(real((i-1)*(k-1),dp)/3+real((j-1)*(l-1),dp)/5)
                    expected(k,l) = expected(k,l)+x(i,j)*cmplx(cos(angle),sin(angle),dp)
                end do
            end do
        end do
    end do
    do mode = 1, 3
        factor = 1.0_dp
        if (mode == 2) factor = 1.0_dp/15
        if (mode == 3) factor = 1.0_dp/sqrt(15.0_dp)
        call fftn(x,y,norm=trim(norms(mode)),status=code)
        call check(code == 0,'forward status')
        call close(y,expected*factor,'independent two-dimensional DFT')
        call ifftn(y,z,norm=trim(norms(mode)),status=code)
        call check(code == 0,'inverse status')
        call close(z,x,'complex round trip')
        call fftn(real_x,y,norm=trim(norms(mode)))
        call fftn(cmplx(real_x,0.0_dp,dp),z,norm=trim(norms(mode)))
        call close(y,z,'real forward overload')
        call ifftn(real_x,y,norm=trim(norms(mode)))
        call ifftn(cmplx(real_x,0.0_dp,dp),z,norm=trim(norms(mode)))
        call close(y,z,'real inverse overload')
        call rfftn(real_x,half,norm=trim(norms(mode)),status=code)
        call check(code == 0,'rfftn status')
        call irfftn(half,real_y,s=[3,5],norm=trim(norms(mode)),status=code)
        call check(code == 0,'irfftn status')
        call check(maxval(abs(real_y-real_x)) < 2e-12_dp,'odd real round trip')
        call ihfftn(real_x,half,norm=trim(norms(mode)),status=code)
        call check(code == 0,'ihfftn status')
        call hfftn(half,real_y,s=[3,5],norm=trim(norms(mode)),status=code)
        call check(code == 0,'hfftn status')
        call check(maxval(abs(real_y-real_x)) < 2e-12_dp,'Hermitian round trip')
        call rfftn(real_x,reordered,axes=[2,1],norm=trim(norms(mode)),status=code)
        call check(code == 0,'half dimension follows axes order')
        call irfftn(reordered,real_y,s=[5,3],axes=[2,1],norm=trim(norms(mode)))
        call check(maxval(abs(real_y-real_x)) < 2e-12_dp,'reordered axes round trip')
        do kind = 1, 4
            call dctn(x,y,type=kind,norm=trim(norms(mode)),status=code)
            call check(code == 0,'complex DCT status')
            call idctn(y,z,type=kind,norm=trim(norms(mode)))
            call close(z,x,'complex DCT round trip')
            call dstn(x,y,type=kind,norm=trim(norms(mode)))
            call idstn(y,z,type=kind,norm=trim(norms(mode)))
            call close(z,x,'complex DST round trip')
            call dctn(real_x,real_y,type=kind,norm=trim(norms(mode)))
            call idctn(real_y,real_z,type=kind,norm=trim(norms(mode)))
            call check(maxval(abs(real_z-real_x)) < 2e-12_dp,'real DCT round trip')
            call dstn(real_x,real_y,type=kind,norm=trim(norms(mode)))
            call idstn(real_y,real_z,type=kind,norm=trim(norms(mode)))
            call check(maxval(abs(real_z-real_x)) < 2e-12_dp,'real DST round trip')
        end do
    end do
    call fftn(x,y,axes=[1])
    do j = 1, 5
        call check(maxval(abs(y(:,j)-fft(x(:,j)))) < 2e-12_dp,'axis one matches 1D')
    end do
    call fftn(x,y,axes=[2])
    do i = 1, 3
        call check(maxval(abs(y(i,:)-fft(x(i,:)))) < 2e-12_dp,'axis two matches 1D')
    end do
    call fftn(x,padded,s=[4,7])
    call ifftn(padded,padded_back,status=code)
    call check(code == 0,'padded inverse status')
    call close(padded_back(:3,:5),x,'zero padding preserves source')
    call check(maxval(abs(padded_back(4,:))) < 2e-12_dp,'row padding')
    call check(maxval(abs(padded_back(:,6:))) < 2e-12_dp,'column padding')
    call fftn(x,y,axes=[integer ::],status=code)
    call check(code == 0,'empty axes identity status')
    call close(y,x,'empty axes identity')
    gap = (-123.0_dp,7.0_dp)
    gap(1:6:2,:) = x
    call fftn(gap(1:6:2,:),gap(2:6:2,:),status=code)
    call check(code == 0,'strided input and output')
    call close(gap(2:6:2,:),expected,'strided output values')
    call close(gap(1:6:2,:),x,'input preserved')
    line = x(:,1)
    call fftn(line,line_out)
    call check(maxval(abs(line_out-fft(line))) < 2e-12_dp,'rank one')
    volume = reshape([(cmplx(sin(real(i,dp)),cos(real(i,dp)),dp),i=1,24)],shape(volume))
    call fftn(volume,vout,axes=[3,1])
    call ifftn(vout,vback,axes=[3,1])
    call check(maxval(abs(vback-volume)) < 2e-12_dp,'rank three and subset axes')
    high = reshape([x(1,1),x(2,1),x(1,2),x(2,2)],shape(high))
    call fftn(high,hout,axes=[15,1])
    call ifftn(hout,hback,axes=[15,1])
    call check(maxval(abs(hback-high)) < 2e-12_dp,'rank fifteen')
    call fftn(empty,y,s=[3,5],status=code)
    call check(code == 0 .and. maxval(abs(y)) == 0,'empty input zero padded')
    call fftn(empty,empty_result,axes=[2],status=code)
    call check(code == 0,'empty unselected batch complex FFT')
    call rfftn(empty,empty_half,axes=[2],status=code)
    call check(code == 0,'empty unselected batch real FFT')
    call fftn(x,y,axes=[1,1],status=code)
    call invalid(code,y,'duplicate axes')
    call fftn(x,y,axes=[0],status=code)
    call invalid(code,y,'zero axis')
    call fftn(x,y,axes=[3],status=code)
    call invalid(code,y,'out of range axis')
    call fftn(x,y,s=[3],axes=[1,2],status=code)
    call invalid(code,y,'length and axes mismatch')
    call fftn(x,y,s=[0,5],status=code)
    call invalid(code,y,'zero transform length')
    call fftn(x,y,norm='invalid',status=code)
    call invalid(code,y,'invalid normalization')
    call fftn(x,single,status=code)
    call invalid(code,single,'wrong output shape')
    call rfftn(real_x,half,axes=[integer ::],status=code)
    call invalid(code,half,'empty real transform axes')
    small = 1
    call dctn(small,small_out,type=1,status=code)
    call check(code == 1,'DCT I length one rejected')
    call dstn(real_x,real_y,type=5,status=code)
    call check(code == 1 .and. all(ieee_is_nan(real_y)),'invalid DST type')
    x(1,1) = cmplx(nan,0.0_dp,dp)
    call fftn(x,y,status=code)
    call invalid(code,y,'NaN input')
    real_x(1,1) = inf
    call rfftn(real_x,half,status=code)
    call invalid(code,half,'infinite real input')
    call irfftn(real_x,real_y,s=[3,5],status=code)
    call check(code == 1,'real half overload finite validation')
    call hfftn(real_x,real_y,s=[3,5],status=code)
    call check(code == 1,'real Hermitian overload finite validation')
    if (failures /= 0) error stop 'multidimensional FFT checks failed'
    print '(a,i0,a)', 'Multidimensional FFT: ',checks,' checks passed'
contains
    subroutine check(ok,label)
        implicit none
        logical, intent(in) :: ok !! condition required to pass
        character(*), intent(in) :: label !! diagnostic description
        checks = checks+1
        if (ok) return
        failures = failures+1
        print '(a)', 'FAIL: '//label
    end subroutine check
    subroutine close(a,b,label)
        implicit none
        complex(dp), intent(in) :: a(:,:), b(:,:) !! actual and expected values
        character(*), intent(in) :: label !! diagnostic description
        call check(maxval(abs(a-b)) < 2e-12_dp*max(1.0_dp,maxval(abs(b))),label)
    end subroutine close
    subroutine invalid(status,y,label)
        implicit none
        integer, intent(in) :: status !! observed completion code
        complex(dp), intent(in) :: y(:,:) !! invalid result to inspect
        character(*), intent(in) :: label !! diagnostic description
        call check(status == 1 .and. all(ieee_is_nan(real(y,dp))),label)
    end subroutine invalid
end program test_fft_nd
