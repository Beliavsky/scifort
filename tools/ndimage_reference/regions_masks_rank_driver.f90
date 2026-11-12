! SPDX-License-Identifier: MIT
! Reference data for labeled extrema, masked morphology and rank footprints.
program regions_masks_rank_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : binary_dilation_nd, binary_erosion_nd, &
        binary_opening_nd, binary_closing_nd, median_filter_nd, rank_filter_nd, &
        percentile_filter_nd, extrema_labels_nd, bounding_boxes_nd, find_objects_nd
    implicit none
    integer, parameter :: nids=7
    logical :: mask(5,6),editable(5,6),output(5,6),binary_foot(2,3)
    logical :: foot_a(2,3),foot_b(3,3)
    real(dp) :: image(5,6),filtered(5,6),lo(nids),hi(nids)
    integer :: labels(5,6),ids(nids),pmin(2,nids),pmax(2,nids)
    integer :: lower(2,nids),upper(2,nids),objects_low(2,5),objects_high(2,5)
    integer :: origins(2,3),origin(2),i,j,op,it,org,m,ft,unit,status
    integer :: passes(4)
    character(len=14) :: modes(8)
    logical :: mask3(4,3,2),editable3(4,3,2),result3(4,3,2),foot3(3,2,2)
    real(dp) :: image3(4,3,2),filter3(4,3,2),mins3(6),maxs3(6)
    integer :: labels3(4,3,2),ids3(6),pmin3(3,6),pmax3(3,6),low3(3,6),up3(3,6)
    integer :: origin3(3,3),shifts3(3),k

    ids=[-2,0,1,2,3,9,1]
    passes=[1,2,0,-1]
    origins(:,1)=[-1,0]
    origins(:,2)=[0,-1]
    origins(:,3)=[0,1]
    modes=[character(len=14) :: 'reflect','grid-mirror','mirror','nearest', &
        'wrap','grid-wrap','constant','grid-constant']
    binary_foot=reshape([.true.,.false.,.true.,.true.,.false.,.true.],[2,3])
    foot_a=reshape([.true.,.false.,.true.,.true.,.true.,.false.],[2,3])
    foot_b=reshape([.true.,.false.,.true.,.false.,.true.,.true.,.true.,.true.,.false.],[3,3])
    do j=1,6
        do i=1,5
            mask(i,j)=modulo(3*i+7*j,5)<2
            editable(i,j)=modulo(i+j,4)/=0
            image(i,j)=sin(real(i+2*j,dp))+real(modulo(i*i+3*j,17),dp)
            labels(i,j)=modulo(i+2*j,6)-2
        end do
    end do
    open(newunit=unit,file='regions_masks_rank_results.bin',form='unformatted', &
        access='stream',status='replace')
    do op=1,4
        do it=1,size(passes)
            do org=1,3
                origin=origins(:,org)
                select case(op)
                case(1)
                    call binary_dilation_nd(mask,output,status,iterations=passes(it), &
                        mask=editable,structure=binary_foot,origin=origin)
                case(2)
                    call binary_erosion_nd(mask,output,status,iterations=passes(it), &
                        mask=editable,structure=binary_foot,origin=origin)
                case(3)
                    call binary_opening_nd(mask,output,status,iterations=passes(it), &
                        mask=editable,structure=binary_foot,origin=origin)
                case(4)
                    call binary_closing_nd(mask,output,status,iterations=passes(it), &
                        mask=editable,structure=binary_foot,origin=origin)
                end select
                if (status/=0) error stop 'masked morphology'
                write(unit) merge(1.0_dp,0.0_dp,output)
            end do
        end do
    end do
    do ft=1,2
        do m=1,size(modes)
            do org=1,3
                origin=origins(:,org)
                ! The 3x3 footprint allows only origins from -1 to +1.
                if (ft==1) then
                    call median_filter_nd(image,[2,3],filtered,status,trim(modes(m)),cval=-3.0_dp, &
                        footprint=foot_a,origin=origin)
                else
                    call median_filter_nd(image,[3,3],filtered,status,trim(modes(m)),cval=-3.0_dp, &
                        footprint=foot_b,origin=origin)
                end if
                if (status/=0) error stop 'median footprint'
                write(unit) filtered
                if (ft==1) then
                    call rank_filter_nd(image,-2,[2,3],filtered,status,trim(modes(m)),cval=-3.0_dp, &
                        footprint=foot_a,origin=origin)
                else
                    call rank_filter_nd(image,-2,[3,3],filtered,status,trim(modes(m)),cval=-3.0_dp, &
                        footprint=foot_b,origin=origin)
                end if
                if (status/=0) error stop 'rank footprint'
                write(unit) filtered
                if (ft==1) then
                    call percentile_filter_nd(image,63.0_dp,[2,3],filtered,status,trim(modes(m)), &
                        cval=-3.0_dp,footprint=foot_a,origin=origin)
                else
                    call percentile_filter_nd(image,63.0_dp,[3,3],filtered,status,trim(modes(m)), &
                        cval=-3.0_dp,footprint=foot_b,origin=origin)
                end if
                if (status/=0) error stop 'percentile footprint'
                write(unit) filtered
            end do
        end do
    end do
    call extrema_labels_nd(image,labels,ids,lo,hi,pmin,pmax,status)
    if (status/=0) error stop 'label extrema'
    write(unit) lo,hi,real(pmin,dp),real(pmax,dp)
    call bounding_boxes_nd(labels,ids,lower,upper,status)
    if (status/=0) error stop 'bounding boxes'
    write(unit) real(lower,dp),real(upper,dp)
    call find_objects_nd(labels,objects_low,objects_high,status)
    if (status/=0) error stop 'find objects'
    write(unit) real(objects_low,dp),real(objects_high,dp)
    ids3=[-2,0,1,2,5,9]
    origin3(:,1)=[0,0,0]
    origin3(:,2)=[-1,-1,0]
    origin3(:,3)=[1,0,-1]
    foot3=reshape([.true.,.false.,.true.,.false.,.true.,.true., &
        .false.,.true.,.true.,.false.,.true.,.true.],[3,2,2])
    do k=1,2
        do j=1,3
            do i=1,4
                mask3(i,j,k)=modulo(7*i+3*j+2*k,6)<3
                editable3(i,j,k)=modulo(2*i+j+4*k,5)/=1
                image3(i,j,k)=sin(real(3*i-j+2*k,dp)) + &
                    real(modulo(i*i+2*j+5*k,13),dp)
                labels3(i,j,k)=modulo(i+2*j+3*k,6)-2
            end do
        end do
    end do
    do op=1,4
        do it=1,3
            do org=1,3
                shifts3=origin3(:,org)
                select case(op)
                case(1)
                    call binary_dilation_nd(mask3,result3,status,mask=editable3, &
                        iterations=passes(it),structure=foot3,origin=shifts3)
                case(2)
                    call binary_erosion_nd(mask3,result3,status,mask=editable3, &
                        iterations=passes(it),structure=foot3,origin=shifts3)
                case(3)
                    call binary_opening_nd(mask3,result3,status,mask=editable3, &
                        iterations=passes(it),structure=foot3,origin=shifts3)
                case(4)
                    call binary_closing_nd(mask3,result3,status,mask=editable3, &
                        iterations=passes(it),structure=foot3,origin=shifts3)
                end select
                if (status/=0) error stop '3D masked morphology'
                write(unit) merge(1.0_dp,0.0_dp,result3)
            end do
        end do
    end do
    do m=1,size(modes)
        do org=1,3
            shifts3=origin3(:,org)
            call median_filter_nd(image3,y=filter3,status=status,mode=trim(modes(m)), &
                cval=-3.0_dp,footprint=foot3,origin=shifts3)
            if (status/=0) error stop '3D median footprint'
            write(unit) filter3
            call rank_filter_nd(image3,2,y=filter3,status=status,mode=trim(modes(m)), &
                cval=-3.0_dp,footprint=foot3,origin=shifts3)
            if (status/=0) error stop '3D rank footprint'
            write(unit) filter3
            call percentile_filter_nd(image3,75.0_dp,y=filter3,status=status, &
                mode=trim(modes(m)),cval=-3.0_dp,footprint=foot3,origin=shifts3)
            if (status/=0) error stop '3D percentile footprint'
            write(unit) filter3
        end do
    end do
    call extrema_labels_nd(image3,labels3,ids3,mins3,maxs3,pmin3,pmax3,status)
    if (status/=0) error stop '3D extrema'
    write(unit) mins3,maxs3,real(pmin3,dp),real(pmax3,dp)
    call bounding_boxes_nd(labels3,ids3,low3,up3,status)
    if (status/=0) error stop '3D bounding boxes'
    write(unit) real(low3,dp),real(up3,dp)
    close(unit)
    print '(a)', 'PASS: region, mask and rank filter driver'
end program regions_masks_rank_driver
