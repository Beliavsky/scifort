! SPDX-License-Identifier: MIT
program measurements_driver
    use scifort_kinds, only : dp
    use scifort_ndimage, only : sum_labels_nd, mean_labels_nd, &
        variance_labels_nd, center_of_mass_nd
    implicit none
    integer, parameter :: ns=7
    real(dp) :: image2(5,6),image3(3,2,4),total(ns),means(ns),vars(ns)
    real(dp) :: com2(2,ns),com3(3,ns)
    real(dp) :: zero(2),cm1(1,1)
    integer :: labels2(5,6),labels3(3,2,4),idx(ns),one(2)
    integer :: i,j,k,unit,code
    idx=[-1,0,1,2,3,99,1]
    do j=1,6
        do i=1,5
            image2(i,j)=sin(real(i+2*j,dp)) + real(modulo(i*i+3*j,17),dp)
            labels2(i,j)=modulo(i+2*j,5)-1
        end do
    end do
    do k=1,4
        do j=1,2
            do i=1,3
                image3(i,j,k)=cos(real(i-j*3+k,dp)) + real(modulo(i*j+5*k,13),dp)
                labels3(i,j,k)=modulo(i+3*j+2*k,5)-1
            end do
        end do
    end do
    open(newunit=unit,file='measurements_results.bin',form='unformatted', &
        access='stream',status='replace')
    call sum_labels_nd(image2,labels2,idx,total,code)
    if (code/=0) error stop 'sum 2D'
    write(unit) total
    call mean_labels_nd(image2,labels2,idx,means,code)
    if (code/=0) error stop 'mean 2D'
    write(unit) means
    call variance_labels_nd(image2,labels2,idx,vars,code)
    if (code/=0) error stop 'variance 2D'
    write(unit) vars
    call center_of_mass_nd(image2,labels2,idx,com2,code)
    if (code/=0) error stop 'centers 2D'
    write(unit) com2
    call sum_labels_nd(image3,labels3,idx,total,code)
    if (code/=0) error stop 'sum 3D'
    write(unit) total
    call mean_labels_nd(image3,labels3,idx,means,code)
    if (code/=0) error stop 'mean 3D'
    write(unit) means
    call variance_labels_nd(image3,labels3,idx,vars,code)
    if (code/=0) error stop 'variance 3D'
    write(unit) vars
    call center_of_mass_nd(image3,labels3,idx,com3,code)
    if (code/=0) error stop 'centers 3D'
    write(unit) com3
    close(unit)
    zero=[-1.0_dp,1.0_dp]
    one=[7,7]
    call center_of_mass_nd(zero,one,[7],cm1,code)
    if (code/=0 .or. cm1(1,1)<huge(1.0_dp)) error stop 'zero-mass positive infinity'
    print '(a)', 'PASS: labeled measurements reference driver'
end program measurements_driver
