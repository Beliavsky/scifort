! SPDX-License-Identifier: MIT
! Emit deterministic ndimage footprint, grayscale and feature-transform results.
program footprint_feature_driver
    use scifort_kinds, only : dp
    use scifort_ndimage
    implicit none
    real(dp) :: a(6,7),out(6,7),dist(6,7),a3(4,3,3),out3(4,3,3)
    logical :: m(6,7),b(6,7),seed(6,7),allowed(6,7),fp3(3,3),fp23(2,3),sym3(3,3)
    logical :: full(6,7),h(6,7),r(6,7),m3(4,3,3),r3(4,3,3),fp3d(2,1,3)
    integer :: ind(2,42),labels(6,7),nfeat,stat,i,j,k,l,unit
    character(len=15), parameter :: modes(8)=[character(len=15) :: 'reflect','mirror','nearest', &
        'wrap','constant','grid-constant','grid-mirror','grid-wrap']

    do j=1,7
        do i=1,6
            a(i,j)=sin(0.38_dp*i+0.67_dp*j)+real(modulo(7*i+11*j+i*j,13),dp)
            m(i,j)=modulo(i*i+3*j+i*j,7)<3
            allowed(i,j)=modulo(i+2*j,5)/=0
            seed(i,j)=modulo(2*i+j,11)==0
        end do
    end do
    fp3=.false.
    fp3(1,2)=.true.
    fp3(2,2)=.true.
    fp3(3,1)=.true.
    fp3(2,3)=.true.
    fp23=.false.
    fp23(1,1)=.true.
    fp23(2,2)=.true.
    fp23(2,3)=.true.
    sym3=.false.
    sym3(2,2)=.true.
    sym3(1,1)=.true.
    sym3(3,3)=.true.

    open(newunit=unit,file='footprint_feature_results.bin',form='unformatted',access='stream',status='replace')
    do k=1,2
        do i=1,2
            if (k==1) then
                call binary_dilation_nd(m,b,stat,iterations=i,border_value=.false.,structure=fp3)
                call binary_erosion_nd(m,r,stat,iterations=i,border_value=.false.,structure=fp3)
                call binary_opening_nd(m,h,stat,iterations=i,border_value=.false.,structure=fp3)
                call binary_closing_nd(m,full,stat,iterations=i,border_value=.false.,structure=fp3)
            else
                call binary_dilation_nd(m,b,stat,iterations=i,border_value=.true.,structure=fp23)
                call binary_erosion_nd(m,r,stat,iterations=i,border_value=.true.,structure=fp23)
                call binary_opening_nd(m,h,stat,iterations=i,border_value=.true.,structure=fp23)
                call binary_closing_nd(m,full,stat,iterations=i,border_value=.true.,structure=fp23)
            end if
            if (stat/=0) error stop 'binary footprint operation status'
            write(unit) merge(1,0,reshape(b,[42]))
            write(unit) merge(1,0,reshape(r,[42]))
            write(unit) merge(1,0,reshape(h,[42]))
            write(unit) merge(1,0,reshape(full,[42]))
        end do
    end do
    call binary_propagation_nd(seed,allowed,b,stat,structure=fp3)
    if (stat/=0) error stop 'binary propagation footprint'
    write(unit) merge(1,0,reshape(b,[42]))
    call binary_propagation_nd(seed,allowed,b,stat,structure=fp23)
    if (stat/=0) error stop 'binary propagation footprint even'
    write(unit) merge(1,0,reshape(b,[42]))
    call binary_fill_holes_nd(m,b,stat,structure=fp3)
    if (stat/=0) error stop 'binary fill holes footprint'
    write(unit) merge(1,0,reshape(b,[42]))
    call label_nd(m,labels,nfeat,stat,structure=sym3)
    if (stat/=0) error stop 'label custom symmetric footprint'
    write(unit) reshape(labels,[42])
    write(unit) nfeat
    call label_nd(m,labels,nfeat,stat,structure=fp3)
    if (stat/=1) error stop 'label asym structure should be rejected'
    call binary_hit_or_miss_nd(m,b,stat)
    if (stat/=0) error stop 'hit-or-miss defaults'
    write(unit) merge(1,0,reshape(b,[42]))
    call binary_hit_or_miss_nd(m,b,stat,structure1=fp3)
    if (stat/=0) error stop 'hit-or-miss asymmetric pattern'
    write(unit) merge(1,0,reshape(b,[42]))
    call binary_hit_or_miss_nd(m,b,stat,structure1=fp23)
    if (stat/=0) error stop 'hit-or-miss even pattern'
    write(unit) merge(1,0,reshape(b,[42]))
    call binary_hit_or_miss_nd(m,b,stat,structure1=fp3,structure2=sym3)
    if (stat/=0) error stop 'hit-or-miss two patterns'
    write(unit) merge(1,0,reshape(b,[42]))
    call binary_hit_or_miss_nd(m,b,stat,structure2=sym3)
    if (stat/=0) error stop 'hit-or-miss explicit background'
    write(unit) merge(1,0,reshape(b,[42]))

    do k=1,8
        call grey_dilation_nd(a,out,stat,footprint=fp3,mode=trim(modes(k)),cval=-2.7_dp)
        if (stat/=0) error stop 'grey dilation status'
        write(unit) reshape(out,[42])
        call grey_erosion_nd(a,out,stat,footprint=fp3,mode=trim(modes(k)),cval=-2.7_dp)
        if (stat/=0) error stop 'grey erosion status'
        write(unit) reshape(out,[42])
        call grey_opening_nd(a,out,stat,footprint=fp3,mode=trim(modes(k)),cval=-2.7_dp)
        if (stat/=0) error stop 'grey opening status'
        write(unit) reshape(out,[42])
        call grey_closing_nd(a,out,stat,footprint=fp3,mode=trim(modes(k)),cval=-2.7_dp)
        if (stat/=0) error stop 'grey closing status'
        write(unit) reshape(out,[42])
        call morphological_gradient_nd(a,out,stat,footprint=fp3,mode=trim(modes(k)),cval=-2.7_dp)
        if (stat/=0) error stop 'morph gradient status'
        write(unit) reshape(out,[42])
        call morphological_laplace_nd(a,out,stat,footprint=fp3,mode=trim(modes(k)),cval=-2.7_dp)
        if (stat/=0) error stop 'morph laplace status'
        write(unit) reshape(out,[42])
        call white_tophat_nd(a,out,stat,footprint=fp3,mode=trim(modes(k)),cval=-2.7_dp)
        if (stat/=0) error stop 'white tophat status'
        write(unit) reshape(out,[42])
        call black_tophat_nd(a,out,stat,footprint=fp3,mode=trim(modes(k)),cval=-2.7_dp)
        if (stat/=0) error stop 'black tophat status'
        write(unit) reshape(out,[42])
        call grey_dilation_nd(a,out,stat,footprint=fp23,mode=trim(modes(k)),cval=-2.7_dp)
        if (stat/=0) error stop 'even grey dilation status'
        write(unit) reshape(out,[42])
        call grey_erosion_nd(a,out,stat,footprint=fp23,mode=trim(modes(k)),cval=-2.7_dp)
        if (stat/=0) error stop 'even grey erosion status'
        write(unit) reshape(out,[42])
        call grey_opening_nd(a,out,stat,footprint=fp23,mode=trim(modes(k)),cval=-2.7_dp)
        if (stat/=0) error stop 'even grey opening status'
        write(unit) reshape(out,[42])
        call grey_closing_nd(a,out,stat,footprint=fp23,mode=trim(modes(k)),cval=-2.7_dp)
        if (stat/=0) error stop 'even grey closing status'
        write(unit) reshape(out,[42])
        call grey_dilation_nd(a,out,stat,sizes=[2,3],mode=trim(modes(k)),cval=-2.7_dp)
        if (stat/=0) error stop 'sizes grey dilation status'
        write(unit) reshape(out,[42])
        call grey_erosion_nd(a,out,stat,sizes=[2,3],mode=trim(modes(k)),cval=-2.7_dp)
        if (stat/=0) error stop 'sizes grey erosion status'
        write(unit) reshape(out,[42])
    end do
    call distance_transform_edt_nd(m,dist,stat,sampling=[1.7_dp,0.6_dp],feature_indices=ind)
    if (stat/=0) error stop 'edt features status'
    write(unit) reshape(dist,[42])
    write(unit) reshape(ind,[2*42])
    full=.true.
    call distance_transform_edt_nd(full,dist,stat,feature_indices=ind)
    if (stat/=0) error stop 'edt all-foreground features status'
    write(unit) reshape(dist,[42])
    write(unit) reshape(ind,[2*42])
    do l=1,3
        do j=1,3
            do i=1,4
                a3(i,j,l)=cos(.24_dp*i+.49_dp*j-.31_dp*l)+real(modulo(11*i+7*j+3*l,9),dp)
                m3(i,j,l)=modulo(i*i+3*j+5*l,6)<3
            end do
        end do
    end do
    fp3d=.false.
    fp3d(1,1,1)=.true.
    fp3d(1,1,3)=.true.
    fp3d(2,1,2)=.true.
    call binary_dilation_nd(m3,r3,stat,structure=fp3d)
    if (stat/=0) error stop '3d footprint dilation'
    write(unit) merge(1,0,reshape(r3,[36]))
    call binary_erosion_nd(m3,r3,stat,structure=fp3d)
    if (stat/=0) error stop '3d footprint erosion'
    write(unit) merge(1,0,reshape(r3,[36]))
    call binary_hit_or_miss_nd(m3,r3,stat,structure1=fp3d)
    if (stat/=0) error stop '3d hit or miss'
    write(unit) merge(1,0,reshape(r3,[36]))
    do k=1,8
        call grey_dilation_nd(a3,out3,stat,footprint=fp3d,mode=trim(modes(k)),cval=-1.5_dp)
        if (stat/=0) error stop '3d grayscale dilation'
        write(unit) reshape(out3,[36])
        call grey_erosion_nd(a3,out3,stat,footprint=fp3d,mode=trim(modes(k)),cval=-1.5_dp)
        if (stat/=0) error stop '3d grayscale erosion'
        write(unit) reshape(out3,[36])
        call morphological_gradient_nd(a3,out3,stat,footprint=fp3d,mode=trim(modes(k)),cval=-1.5_dp)
        if (stat/=0) error stop '3d grayscale gradient'
        write(unit) reshape(out3,[36])
    end do
    close(unit)
    print '(a)', 'PASS: footprint/feature Fortran comparison driver'
end program footprint_feature_driver
