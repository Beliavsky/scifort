! SPDX-License-Identifier: MIT
! Streamed reference probe for real, complex and mixed assumed-rank IIR filters.
program lfilter_nd_driver
    use scifort_kinds, only : dp
    use scifort_signal_lfilter, only : lfilter_nd
    implicit none
    integer :: case_count, c, rank_id, ax, combo, with_state, code, k
    integer :: dims5(5), dims15(15)
    real(dp) :: b(3), a(2), r5(4,3,2,3,5), ry5(4,3,2,3,5)
    complex(dp) :: cb(3), ca(2), cx5(4,3,2,3,5), cy5(4,3,2,3,5)
    real(dp) :: r15(3,2,1,2,1,1,1,1,1,1,1,1,1,1,4)
    real(dp) :: ry15(3,2,1,2,1,1,1,1,1,1,1,1,1,1,4)
    complex(dp) :: cx15(3,2,1,2,1,1,1,1,1,1,1,1,1,1,4)
    complex(dp) :: cy15(3,2,1,2,1,1,1,1,1,1,1,1,1,1,4)
    real(dp), allocatable :: rzi5(:,:,:,:,:), rzf5(:,:,:,:,:)
    complex(dp), allocatable :: czi5(:,:,:,:,:), czf5(:,:,:,:,:)
    real(dp), allocatable :: rzi15(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
    real(dp), allocatable :: rzf15(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
    complex(dp), allocatable :: czi15(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
    complex(dp), allocatable :: czf15(:,:,:,:,:,:,:,:,:,:,:,:,:,:,:)
    complex(dp), allocatable :: flat(:)

    read(*,*) case_count
    do c = 1, case_count
        read(*,*) rank_id,ax,combo,with_state
        read(*,*) cb
        read(*,*) ca
        b = real(cb,dp)
        a = real(ca,dp)
        select case (rank_id)
        case (5)
            read(*,*) cx5
            r5 = real(cx5,dp)
            dims5 = shape(cx5)
            dims5(ax) = max(size(a),size(b))-1
            allocate(czi5(dims5(1),dims5(2),dims5(3),dims5(4),dims5(5)))
            allocate(czf5(dims5(1),dims5(2),dims5(3),dims5(4),dims5(5)))
            allocate(rzi5(dims5(1),dims5(2),dims5(3),dims5(4),dims5(5)))
            allocate(rzf5(dims5(1),dims5(2),dims5(3),dims5(4),dims5(5)))
            czi5 = (0.0_dp,0.0_dp)
            if (with_state /= 0) read(*,*) czi5
            rzi5 = real(czi5,dp)
            select case (combo)
            case (0)
                call lfilter_nd(b,a,r5,ry5,zi=rzi5,zf=rzf5,status=code,axis=ax)
                cy5 = cmplx(ry5,0.0_dp,dp)
                czf5 = cmplx(rzf5,0.0_dp,dp)
            case (1)
                call lfilter_nd(cb,a,r5,cy5,zi=czi5,zf=czf5,status=code,axis=ax)
            case (2)
                call lfilter_nd(b,ca,r5,cy5,zi=czi5,zf=czf5,status=code,axis=ax)
            case (3)
                call lfilter_nd(cb,ca,r5,cy5,zi=czi5,zf=czf5,status=code,axis=ax)
            case (4)
                call lfilter_nd(b,a,cx5,cy5,zi=czi5,zf=czf5,status=code,axis=ax)
            case (5)
                call lfilter_nd(cb,a,cx5,cy5,zi=czi5,zf=czf5,status=code,axis=ax)
            case (6)
                call lfilter_nd(b,ca,cx5,cy5,zi=czi5,zf=czf5,status=code,axis=ax)
            case (7)
                call lfilter_nd(cb,ca,cx5,cy5,zi=czi5,zf=czf5,status=code,axis=ax)
            end select
            flat = reshape(cy5,[size(cy5)])
            write(*,'(i0)') code
            do k = 1, size(flat)
                write(*,'(es26.17,1x,es26.17)') real(flat(k),dp),aimag(flat(k))
            end do
            flat = reshape(czf5,[size(czf5)])
            do k = 1, size(flat)
                write(*,'(es26.17,1x,es26.17)') real(flat(k),dp),aimag(flat(k))
            end do
            deallocate(czi5,czf5,rzi5,rzf5)
        case (15)
            read(*,*) cx15
            r15 = real(cx15,dp)
            dims15 = shape(cx15)
            dims15(ax) = max(size(a),size(b))-1
            allocate(czi15(dims15(1),dims15(2),dims15(3),dims15(4),dims15(5),dims15(6),dims15(7), &
                dims15(8),dims15(9),dims15(10),dims15(11),dims15(12),dims15(13),dims15(14),dims15(15)))
            allocate(czf15(dims15(1),dims15(2),dims15(3),dims15(4),dims15(5),dims15(6),dims15(7), &
                dims15(8),dims15(9),dims15(10),dims15(11),dims15(12),dims15(13),dims15(14),dims15(15)))
            allocate(rzi15(dims15(1),dims15(2),dims15(3),dims15(4),dims15(5),dims15(6),dims15(7), &
                dims15(8),dims15(9),dims15(10),dims15(11),dims15(12),dims15(13),dims15(14),dims15(15)))
            allocate(rzf15(dims15(1),dims15(2),dims15(3),dims15(4),dims15(5),dims15(6),dims15(7), &
                dims15(8),dims15(9),dims15(10),dims15(11),dims15(12),dims15(13),dims15(14),dims15(15)))
            czi15 = (0.0_dp,0.0_dp)
            if (with_state /= 0) read(*,*) czi15
            rzi15 = real(czi15,dp)
            select case (combo)
            case (0)
                call lfilter_nd(b,a,r15,ry15,zi=rzi15,zf=rzf15,status=code,axis=ax)
                cy15 = cmplx(ry15,0.0_dp,dp)
                czf15 = cmplx(rzf15,0.0_dp,dp)
            case (1)
                call lfilter_nd(cb,a,r15,cy15,zi=czi15,zf=czf15,status=code,axis=ax)
            case (2)
                call lfilter_nd(b,ca,r15,cy15,zi=czi15,zf=czf15,status=code,axis=ax)
            case (3)
                call lfilter_nd(cb,ca,r15,cy15,zi=czi15,zf=czf15,status=code,axis=ax)
            case (4)
                call lfilter_nd(b,a,cx15,cy15,zi=czi15,zf=czf15,status=code,axis=ax)
            case (5)
                call lfilter_nd(cb,a,cx15,cy15,zi=czi15,zf=czf15,status=code,axis=ax)
            case (6)
                call lfilter_nd(b,ca,cx15,cy15,zi=czi15,zf=czf15,status=code,axis=ax)
            case (7)
                call lfilter_nd(cb,ca,cx15,cy15,zi=czi15,zf=czf15,status=code,axis=ax)
            end select
            flat = reshape(cy15,[size(cy15)])
            write(*,'(i0)') code
            do k = 1, size(flat)
                write(*,'(es26.17,1x,es26.17)') real(flat(k),dp),aimag(flat(k))
            end do
            flat = reshape(czf15,[size(czf15)])
            do k = 1, size(flat)
                write(*,'(es26.17,1x,es26.17)') real(flat(k),dp),aimag(flat(k))
            end do
            deallocate(czi15,czf15,rzi15,rzf15)
        end select
    end do
end program lfilter_nd_driver
