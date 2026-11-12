! SPDX-License-Identifier: MIT
! Rank-independent SOS filtering, streaming states and mixed real/complex inputs.
program test_signal_sos_nd
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan
    use scifort_signal, only : dp, sosfilt, sosfilt_nd, &
        sos_status_success, sos_status_invalid_input
    implicit none

    real(dp) :: sos(2,6), bad_sos(2,6)
    complex(dp) :: csos(2,6)
    real(dp) :: x1(11), y1(11), z1(2,2), end1(2,2)
    real(dp) :: x5(6,4,2,3,5), y5(6,4,2,3,5), split5(6,4,2,3,5)
    real(dp) :: zi5(2,6,2,2,3,5), zf5(2,6,2,2,3,5)
    real(dp) :: block5(2,6,2,2,3,5), final5(2,6,2,2,3,5)
    real(dp) :: shortzi(2,6,2,2,3,4), wrongshape(6,4,2,3,4)
    real(dp) :: x14(2,2,2,1,1,1,1,1,1,1,1,1,1,5)
    real(dp) :: y14(2,2,2,1,1,1,1,1,1,1,1,1,1,5)
    real(dp) :: y14p(2,2,2,1,1,1,1,1,1,1,1,1,1,5)
    real(dp) :: zi14(2,2,2,2,1,1,1,1,1,1,1,1,1,1,2)
    real(dp) :: zf14(2,2,2,2,1,1,1,1,1,1,1,1,1,1,2)
    real(dp) :: packed14(32), final14(32)
    real(dp) :: x15(4,2,3,1,1,1,1,1,1,1,1,1,1,1,6)
    real(dp) :: y15(4,2,3,1,1,1,1,1,1,1,1,1,1,1,6)
    real(dp) :: zi15(96), zf15(96), st(2,2)
    real(dp) :: xe(0,3), ye(0,3), ze(2,2,3), end_e(2,2,3)
    complex(dp) :: cx5(6,4,2,3,5), cy5(6,4,2,3,5)
    complex(dp) :: czi5(2,6,2,2,3,5), czf5(2,6,2,2,3,5)
    complex(dp) :: cx15(4,2,3,1,1,1,1,1,1,1,1,1,1,1,6)
    complex(dp) :: cy15(4,2,3,1,1,1,1,1,1,1,1,1,1,1,6)
    complex(dp) :: czi15(96), czf15(96), cst(2,2)
    real(dp), allocatable :: reference(:), last_st(:,:)
    complex(dp), allocatable :: creference(:), clast(:,:)
    integer :: i, code, offset

    sos(1,:) = [0.2_dp,0.3_dp,0.1_dp,1.0_dp,-0.25_dp,0.08_dp]
    sos(2,:) = [1.0_dp,-0.1_dp,0.2_dp,1.0_dp,-0.3_dp,0.1_dp]
    csos = cmplx(sos,0.0_dp,dp)
    csos(1,2) = csos(1,2)+(0.0_dp,0.04_dp)
    csos(2,5) = csos(2,5)-(0.0_dp,0.07_dp)

    x1 = [(sin(real(i,dp)),i=1,11)]
    z1 = reshape([(0.03_dp*i,i=1,4)],[2,2])
    call sosfilt_nd(sos,x1,y1,zi=z1,zf=end1,status=code)
    if (code /= sos_status_success) error stop 'rank1 status'
    reference = sosfilt(sos,x1,zi=z1,zf=last_st)
    if (any(abs(end1-last_st)>1.e-13_dp)) error stop 'rank1 final state'
    if (any(abs(y1-reference)>1.e-13_dp)) error stop 'rank1 result'

    x5 = reshape([(0.2_dp*sin(0.07_dp*i)+0.003_dp*i,i=1,size(x5))],shape(x5))
    zi5 = reshape([(0.05_dp*cos(0.11_dp*i),i=1,size(zi5))],shape(zi5))
    cx5 = cmplx(x5,-0.2_dp*x5+0.15_dp,dp)
    czi5 = cmplx(zi5,-0.1_dp*zi5,dp)

    call sosfilt_nd(sos,x5,y5,zi=zi5,zf=zf5,axis=2,status=code)
    if (code /= sos_status_success) error stop 'rank5 real status'
    reference = sosfilt(sos,x5(3,:,1,2,4),zi=zi5(:,3,:,1,2,4),zf=last_st)
    if (any(abs(y5(3,:,1,2,4)-reference)>1.e-13_dp)) error stop 'rank5 real lane'
    if (any(abs(zf5(:,3,:,1,2,4)-last_st)>1.e-13_dp)) error stop 'rank5 real state'
    call sosfilt_nd(sos,x5(:,1:2,:,:,:),split5(:,1:2,:,:,:), &
        zi=zi5,zf=block5,axis=2,status=code)
    if (code /= sos_status_success) error stop 'first block status'
    call sosfilt_nd(sos,x5(:,3:4,:,:,:),split5(:,3:4,:,:,:), &
        zi=block5,zf=final5,axis=2,status=code)
    if (code /= sos_status_success) error stop 'second block status'
    if (any(abs(y5-split5)>1.e-13_dp)) error stop 'rank5 split result'
    if (any(abs(zf5-final5)>1.e-13_dp)) error stop 'rank5 split state'

    ! The three complex-output combinations must match 1D SOS references.
    call sosfilt_nd(csos,cx5,cy5,zi=czi5,zf=czf5,axis=2,status=code)
    if (code /= sos_status_success) error stop 'complex SOS and samples'
    creference = sosfilt(csos,cx5(3,:,1,2,4),zi=czi5(:,3,:,1,2,4),zf=clast)
    if (any(abs(cy5(3,:,1,2,4)-creference)>1.e-13_dp)) error stop 'rank5 cc output'
    if (any(abs(czf5(:,3,:,1,2,4)-clast)>1.e-13_dp)) error stop 'rank5 cc final state'
    call sosfilt_nd(sos,cx5,cy5,zi=czi5,zf=czf5,axis=2,status=code)
    if (code /= sos_status_success) error stop 'real SOS complex samples'
    creference = sosfilt(cmplx(sos,0.0_dp,dp),cx5(3,:,1,2,4),zi=czi5(:,3,:,1,2,4))
    if (any(abs(cy5(3,:,1,2,4)-creference)>1.e-13_dp)) error stop 'rank5 rc output'
    call sosfilt_nd(csos,x5,cy5,zi=czi5,zf=czf5,axis=2,status=code)
    if (code /= sos_status_success) error stop 'complex SOS real samples'
    creference = sosfilt(csos,cmplx(x5(3,:,1,2,4),0.0_dp,dp),zi=czi5(:,3,:,1,2,4))
    if (any(abs(cy5(3,:,1,2,4)-creference)>1.e-13_dp)) error stop 'rank5 cr output'

    ! The CONTIGUOUS dummy makes strided actual arguments copy in/out.
    y5 = -999.0_dp
    zf5 = -888.0_dp
    call sosfilt_nd(sos,x5(1:6:2,:,:,:,:),y5(1:6:2,:,:,:,:), &
        zi=zi5(:,1:6:2,:,:,:,:),zf=zf5(:,1:6:2,:,:,:,:),axis=2,status=code)
    if (code /= sos_status_success) error stop 'strided status'
    if (any(y5(2,:,:,:,:) /= -999.0_dp)) error stop 'strided output clobber'
    if (any(zf5(:,2,:,:,:,:) /= -888.0_dp)) error stop 'strided states clobber'
    reference = sosfilt(sos,x5(3,:,1,2,4),zi=zi5(:,3,:,1,2,4))
    if (any(abs(y5(3,:,1,2,4)-reference)>1.e-13_dp)) error stop 'strided lane'

    ! Rank 14 can use rank-15 shaped state arrays or packed 1D buffers.
    x14 = reshape([(0.2_dp*sin(0.12_dp*i),i=1,size(x14))],shape(x14))
    zi14 = reshape([(0.004_dp*i,i=1,size(zi14))],shape(zi14))
    packed14 = reshape(zi14,[size(zi14)])
    call sosfilt_nd(sos,x14,y14,zi=zi14,zf=zf14,axis=14,status=code)
    if (code /= sos_status_success) error stop 'rank14 shaped state'
    call sosfilt_nd(sos,x14,y14p,zi=packed14,zf=final14,axis=14,status=code)
    if (code /= sos_status_success) error stop 'rank14 packed state'
    if (any(abs(y14-y14p)>1.e-13_dp)) error stop 'rank14 packed output'
    if (any(abs(reshape(zf14,[size(zf14)])-final14)>1.e-13_dp)) &
        error stop 'rank14 packed state comparison'

    ! Rank 15 requires packed rank-one SOS state buffers: 2 states per section per lane.
    x15 = reshape([(cos(0.043_dp*i)+0.002_dp*i,i=1,size(x15))],shape(x15))
    zi15 = [(0.005_dp*i,i=1,size(zi15))]
    call sosfilt_nd(sos,x15,y15,zi=zi15,zf=zf15,axis=15,status=code)
    if (code /= sos_status_success) error stop 'rank15 packed states'
    offset = (2-1)+4*(1-1)+8*(2-1)
    st(:,1) = zi15(2*offset+1:2*offset+2)
    st(:,2) = zi15(2*offset+49:2*offset+50)
    reference = sosfilt(sos,x15(2,1,2,1,1,1,1,1,1,1,1,1,1,1,:),zi=st,zf=last_st)
    if (any(abs(y15(2,1,2,1,1,1,1,1,1,1,1,1,1,1,:)-reference)>1.e-13_dp)) &
        error stop 'rank15 real lane'
    if (any(abs(zf15(2*offset+1:2*offset+2)-last_st(:,1))>1.e-13_dp)) &
        error stop 'rank15 real state 1'
    if (any(abs(zf15(2*offset+49:2*offset+50)-last_st(:,2))>1.e-13_dp)) &
        error stop 'rank15 real state 2'
    call sosfilt_nd(sos,x15,y15,axis=1,status=code)
    if (code /= sos_status_success) error stop 'rank15 axis 1'
    reference = sosfilt(sos,x15(:,1,2,1,1,1,1,1,1,1,1,1,1,1,3))
    if (any(abs(y15(:,1,2,1,1,1,1,1,1,1,1,1,1,1,3)-reference)>1.e-13_dp)) &
        error stop 'rank15 axis1 output'

    cx15 = cmplx(x15,0.15_dp*x15-0.05_dp,dp)
    czi15 = cmplx(zi15,-0.03_dp*zi15,dp)
    call sosfilt_nd(csos,cx15,cy15,zi=czi15,zf=czf15,axis=15,status=code)
    if (code /= sos_status_success) error stop 'rank15 complex status'
    cst(:,1) = czi15(2*offset+1:2*offset+2)
    cst(:,2) = czi15(2*offset+49:2*offset+50)
    creference = sosfilt(csos,cx15(2,1,2,1,1,1,1,1,1,1,1,1,1,1,:),zi=cst,zf=clast)
    if (any(abs(cy15(2,1,2,1,1,1,1,1,1,1,1,1,1,1,:)-creference)>1.e-13_dp)) &
        error stop 'rank15 complex lane'
    call sosfilt_nd(csos,x15,cy15,zi=czi15,zf=czf15,axis=15,status=code)
    if (code /= sos_status_success) error stop 'rank15 mixed status'

    ! Empty filtered dimension must preserve initial states for all other lanes.
    ze = reshape([(0.07_dp*i,i=1,size(ze))],shape(ze))
    call sosfilt_nd(sos,xe,ye,zi=ze,zf=end_e,axis=1,status=code)
    if (code /= sos_status_success) error stop 'empty filtered axis status'
    if (any(abs(end_e-ze)>0.0_dp)) error stop 'empty filtered axis state'
    call sosfilt_nd(sos,x5,y5,axis=0,status=code)
    if (code /= sos_status_invalid_input .or. .not. all(ieee_is_nan(y5))) &
        error stop 'bad axis status or NaN result'
    call sosfilt_nd(sos,x5,wrongshape,status=code)
    if (code /= sos_status_invalid_input) error stop 'bad output shape'
    call sosfilt_nd(sos,x5,y5,zi=shortzi,axis=2,status=code)
    if (code /= sos_status_invalid_input) error stop 'bad state shape'
    bad_sos = sos
    bad_sos(1,4) = 2.0_dp
    call sosfilt_nd(bad_sos,x5,y5,zf=zf5,axis=2,status=code)
    if (code /= sos_status_invalid_input .or. .not. all(ieee_is_nan(y5))) &
        error stop 'bad SOS matrix result'
    if (.not. all(ieee_is_nan(zf5))) error stop 'bad SOS matrix state'
    bad_sos = sos
    bad_sos(2,3) = ieee_value(0.0_dp,ieee_quiet_nan)
    call sosfilt_nd(bad_sos,x5,y5,axis=2,status=code)
    if (code /= sos_status_invalid_input) error stop 'NaN SOS coeff'
    zi5(1,1,1,1,1,1) = ieee_value(0.0_dp,ieee_quiet_nan)
    call sosfilt_nd(sos,x5,y5,zi=zi5,axis=2,status=code)
    if (code /= sos_status_invalid_input) error stop 'NaN state'
    print '(a)', 'test_signal_sos_nd passed'
end program test_signal_sos_nd
