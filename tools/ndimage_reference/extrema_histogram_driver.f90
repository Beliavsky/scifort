! SPDX-License-Identifier: MIT
! Reference results for custom extrema footprints and selected-region statistics.
program extrema_histogram_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : minimum_filter_nd, maximum_filter_nd, &
        standard_deviation_labels_nd, count_labels_nd, histogram_labels_nd
    implicit none
    real(dp) :: image2(5,6),res2(5,6),image3(4,3,2),res3(4,3,2)
    real(dp) :: stds(6)
    integer :: lab2(5,6),lab3(4,3,2),ids(6),counts(6),hist(7,6)
    logical :: foot_a(2,3),foot_b(3,3),foot3(3,2,2)
    integer :: origins2(2,3),origins3(3,3)
    integer :: i,j,k,ft,m,org,op,status,unit
    character(len=14) :: modes(8)

    ids=[-2,0,1,2,9,1]
    modes=[character(len=14) :: 'reflect','grid-mirror','mirror','nearest', &
        'wrap','grid-wrap','constant','grid-constant']
    foot_a=reshape([.true.,.false.,.true.,.true.,.true.,.false.],[2,3])
    foot_b=reshape([.true.,.false.,.true.,.false.,.true.,.true.,.true.,.true.,.false.],[3,3])
    foot3=reshape([.true.,.false.,.true.,.false.,.true.,.true., &
        .false.,.true.,.true.,.false.,.true.,.true.],[3,2,2])
    origins2(:,1)=[0,0]
    origins2(:,2)=[-1,-1]
    origins2(:,3)=[0,1]
    origins3(:,1)=[0,0,0]
    origins3(:,2)=[-1,-1,0]
    origins3(:,3)=[1,0,-1]
    do j=1,6
        do i=1,5
            image2(i,j)=sin(real(i*2+j,dp))+real(modulo(i*i+3*j,13),dp)
            lab2(i,j)=modulo(i+2*j,6)-2
        end do
    end do
    do k=1,2
        do j=1,3
            do i=1,4
                image3(i,j,k)=cos(real(2*i-j+3*k,dp))+real(modulo(i*i+2*j+4*k,13),dp)
                lab3(i,j,k)=modulo(i+2*j+3*k,6)-2
            end do
        end do
    end do
    open(newunit=unit,file='extrema_histogram_results.bin',form='unformatted', &
        access='stream',status='replace')
    do ft=1,2
        do m=1,size(modes)
            do org=1,3
                do op=1,2
                    if (ft==1) then
                        if (op==1) then
                            call minimum_filter_nd(image2,y=res2,status=status,footprint=foot_a, &
                                origin=origins2(:,org),mode=trim(modes(m)),cval=-3.25_dp)
                        else
                            call maximum_filter_nd(image2,y=res2,status=status,footprint=foot_a, &
                                origin=origins2(:,org),mode=trim(modes(m)),cval=-3.25_dp)
                        end if
                    else
                        if (op==1) then
                            call minimum_filter_nd(image2,y=res2,status=status,footprint=foot_b, &
                                origin=origins2(:,org),mode=trim(modes(m)),cval=-3.25_dp)
                        else
                            call maximum_filter_nd(image2,y=res2,status=status,footprint=foot_b, &
                                origin=origins2(:,org),mode=trim(modes(m)),cval=-3.25_dp)
                        end if
                    end if
                    if (status/=0) error stop '2D extrema footprint'
                    write(unit) res2
                end do
            end do
        end do
    end do
    do m=1,size(modes)
        do org=1,3
            do op=1,2
                if (op==1) then
                    call minimum_filter_nd(image3,y=res3,status=status,footprint=foot3, &
                        origin=origins3(:,org),mode=trim(modes(m)),cval=-3.25_dp)
                else
                    call maximum_filter_nd(image3,y=res3,status=status,footprint=foot3, &
                        origin=origins3(:,org),mode=trim(modes(m)),cval=-3.25_dp)
                end if
                if (status/=0) error stop '3D extrema footprint'
                write(unit) res3
            end do
        end do
    end do
    do m=1,size(modes)
        do op=1,2
            if (op==1) then
                call minimum_filter_nd(image2,[2,3],res2,status,mode=trim(modes(m)), &
                    origin=[-1,1],cval=-3.25_dp)
            else
                call maximum_filter_nd(image2,[2,3],res2,status,mode=trim(modes(m)), &
                    origin=[-1,1],cval=-3.25_dp)
            end if
            if (status/=0) error stop 'rectangular extrema origin'
            write(unit) res2
        end do
    end do
    call count_labels_nd(lab2,ids,counts,status)
    if (status/=0) error stop '2D counts'
    write(unit) real(counts,dp)
    call histogram_labels_nd(image2,lab2,ids,-2.0_dp,14.0_dp,hist,status)
    if (status/=0) error stop '2D histogram'
    write(unit) real(hist,dp)
    call standard_deviation_labels_nd(image2,lab2,ids,stds,status)
    if (status/=0) error stop '2D standard deviation'
    write(unit) stds
    call count_labels_nd(lab3,ids,counts,status)
    if (status/=0) error stop '3D counts'
    write(unit) real(counts,dp)
    call histogram_labels_nd(image3,lab3,ids,-2.0_dp,14.0_dp,hist,status)
    if (status/=0) error stop '3D histogram'
    write(unit) real(hist,dp)
    call standard_deviation_labels_nd(image3,lab3,ids,stds,status)
    if (status/=0) error stop '3D standard deviation'
    write(unit) stds
    close(unit)
end program extrema_histogram_driver
