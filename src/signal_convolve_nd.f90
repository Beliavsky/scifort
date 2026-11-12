! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! SciPy-style N-dimensional direct/FFT convolution and correlation.
! The public caller-output interfaces accept ranks 1..15 and all real/complex
! operand combinations. Axis subsets and singleton broadcasting mirror the
! SciFort fftconvolve_nd interface (one-based axes).
module scifort_signal_convolve_nd
    use, intrinsic :: iso_c_binding, only : c_loc, c_f_pointer
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    use scifort_signal_fftconvolve, only : fftconvolve_flat
    implicit none
    private

    integer, parameter, public :: convolve_nd_success = 0
    integer, parameter, public :: convolve_nd_invalid_input = 1

    interface convolve_nd
        module procedure convolve_nd_rr, convolve_nd_cc, convolve_nd_rc, convolve_nd_cr
    end interface convolve_nd
    interface correlate_nd
        module procedure correlate_nd_rr, correlate_nd_cc, correlate_nd_rc, correlate_nd_cr
    end interface correlate_nd
    public :: convolve_nd, correlate_nd, convolution_shape_nd

contains

    pure function checked_product(dims) result(total)
        integer, intent(in) :: dims(:) !! array extents
        integer(int64) :: total
        integer :: j

        total = 1_int64
        do j=1,size(dims)
            if (dims(j) < 1) then
                total = 0_int64
                return
            end if
            if (total > huge(total)/int(dims(j),int64)) then
                total = 0_int64
                return
            end if
            total = total*int(dims(j),int64)
        end do
    end function checked_product

    pure subroutine selected_axes(rank, axes, selected, ok)
        integer, intent(in) :: rank !! number of array dimensions
        integer, intent(in), optional :: axes(:) !! selected 1-based axes, default all
        logical, intent(out) :: selected(rank) !! true for convolutions axes
        logical, intent(out) :: ok !! false for repeated, empty, invalid axes
        integer :: i

        selected = .true.
        ok = rank >= 1 .and. rank <= 15
        if (.not. ok) return
        if (.not. present(axes)) return
        selected = .false.
        ok = size(axes)>0
        if (.not. ok) return
        do i=1,size(axes)
            if (axes(i)<1 .or. axes(i)>rank) then
                ok = .false.
                return
            end if
            if (selected(axes(i))) then
                ok = .false.
                return
            end if
            selected(axes(i)) = .true.
        end do
    end subroutine selected_axes

    function convolution_shape_nd(n,m,mode,axes,status) result(count)
        integer, intent(in) :: n(:) !! shape of first operand
        integer, intent(in) :: m(:) !! shape of second operand
        character(*), intent(in), optional :: mode !! full (default), same or valid
        integer, intent(in), optional :: axes(:) !! one-based convolution axes
        integer, intent(out), optional :: status !! zero on success
        integer, allocatable :: count(:)
        integer :: j, r
        logical :: ok
        logical, allocatable :: selected(:)
        character(:), allocatable :: which

        if (present(status)) status = convolve_nd_invalid_input
        allocate(count(0))
        r = size(n)
        if (r<1 .or. r>15 .or. size(m)/=r) return
        if (any(n<1) .or. any(m<1)) return
        allocate(selected(r))
        call selected_axes(r,axes,selected,ok)
        if (.not. ok) return
        which = 'full'
        if (present(mode)) which = trim(mode)
        select case(which)
        case ('full','same','valid')
        case default
            return
        end select
        if (which == 'valid') then
            if (.not. (all(pack(n>=m,selected)) .or. &
                       all(pack(n<=m,selected)))) return
        end if
        deallocate(count)
        allocate(count(r))
        do j=1,r
            if (selected(j)) then
                if (n(j)>huge(j)-m(j)+1) exit
                select case(which)
                case ('full')
                    count(j) = n(j)+m(j)-1
                case ('same')
                    count(j) = n(j)
                case ('valid')
                    count(j) = abs(n(j)-m(j))+1
                end select
            else
                if (n(j)/=m(j) .and. n(j)/=1 .and. m(j)/=1) exit
                count(j) = max(n(j),m(j))
                if (which == 'same') count(j) = n(j)
            end if
        end do
        if (j <= r) then
            deallocate(count)
            allocate(count(0))
            return
        end if
        if (checked_product(count)==0_int64) then
            deallocate(count)
            allocate(count(0))
            return
        end if
        if (present(status)) status = convolve_nd_success
    end function convolution_shape_nd

    pure subroutine reverse_conjugate_kernel(b,m,selected,flipped)
        complex(dp), intent(in) :: b(:) !! original column-major kernel
        integer, intent(in) :: m(:) !! kernel dimensions
        logical, intent(in) :: selected(:) !! axes along which to reverse
        complex(dp), intent(out) :: flipped(:) !! reversed and conjugated kernel
        integer(int64) :: k,q,pos,stride,coord
        integer :: j

        do k=0_int64,size(b,kind=int64)-1_int64
            q = k
            pos = 0_int64
            stride = 1_int64
            do j=1,size(m)
                coord = mod(q,int(m(j),int64))
                q = q/int(m(j),int64)
                if (selected(j)) coord = int(m(j)-1,int64)-coord
                pos = pos+coord*stride
                stride = stride*int(m(j),int64)
            end do
            flipped(k+1_int64) = conjg(b(pos+1_int64))
        end do
    end subroutine reverse_conjugate_kernel

    subroutine direct_flat(a,b,n,m,count,selected,which,result)
        complex(dp), intent(in) :: a(:) !! packed first operand
        complex(dp), intent(in) :: b(:) !! packed second operand
        integer, intent(in) :: n(:) !! shape of a
        integer, intent(in) :: m(:) !! shape of b
        integer, intent(in) :: count(:) !! shape of the output
        logical, intent(in) :: selected(:) !! convolution axes
        character(*), intent(in) :: which !! full/same/valid
        complex(dp), intent(out) :: result(:) !! packed output
        integer(int64) :: k, v, q, qb, ai, bi, sa, sb, bc, ac, fc
        integer(int64) :: kernels, stride_a(size(n)), stride_b(size(n))
        integer :: j, first(size(n)), fullcoord(size(n))
        logical :: inside

        stride_a(1) = 1_int64
        stride_b(1) = 1_int64
        do j=2,size(n)
            stride_a(j) = stride_a(j-1)*int(n(j-1),int64)
            stride_b(j) = stride_b(j-1)*int(m(j-1),int64)
        end do
        kernels = 1_int64
        first = 1
        do j=1,size(n)
            if (.not. selected(j)) cycle
            kernels = kernels*int(m(j),int64)
            select case(which)
            case ('same')
                first(j) = 1+(m(j)-1)/2
            case ('valid')
                first(j) = min(n(j),m(j))
            end select
        end do
        result = (0.0_dp,0.0_dp)
        do k=0_int64,size(result,kind=int64)-1_int64
            q = k
            ai = 1_int64
            bi = 1_int64
            do j=1,size(n)
                fc = mod(q,int(count(j),int64))
                q = q/int(count(j),int64)
                if (selected(j)) then
                    fullcoord(j) = int(fc)+first(j)-1
                else
                    ! For same mode, SciPy centers the larger untransformed
                    ! broadcast dimension before extracting n(j) samples.
                    if (which == 'same') fc = fc+int((max(n(j),m(j))-n(j))/2,int64)
                    ai = ai+min(fc,int(n(j)-1,int64))*stride_a(j)
                    bi = bi+min(fc,int(m(j)-1,int64))*stride_b(j)
                end if
            end do
            do v=0_int64,kernels-1_int64
                qb = v
                sa = ai
                sb = bi
                inside = .true.
                do j=1,size(n)
                    if (.not. selected(j)) cycle
                    bc = mod(qb,int(m(j),int64))
                    qb = qb/int(m(j),int64)
                    ac = int(fullcoord(j),int64)-bc
                    if (ac<0_int64 .or. ac>=int(n(j),int64)) then
                        inside = .false.
                        exit
                    end if
                    sa = sa+ac*stride_a(j)
                    sb = sb+bc*stride_b(j)
                end do
                if (inside) result(k+1_int64) = result(k+1_int64)+a(sa)*b(sb)
            end do
        end do
    end subroutine direct_flat

    subroutine convolve_core(a,b,n,m,ydims,result,mode,method,axes,correlation,status)
        complex(dp), intent(in) :: a(:) !! packed first operand
        complex(dp), intent(in) :: b(:) !! packed second operand
        integer, intent(in) :: n(:) !! shape of a
        integer, intent(in) :: m(:) !! shape of b
        integer, intent(in) :: ydims(:) !! shape of preallocated output
        complex(dp), allocatable, intent(out) :: result(:) !! computed output
        character(*), intent(in), optional :: mode !! full/same/valid
        character(*), intent(in), optional :: method !! auto/direct/fft
        integer, intent(in), optional :: axes(:) !! one-based axes
        logical, intent(in) :: correlation !! conjugate and reverse kernel
        integer, intent(out) :: status !! zero on success
        integer, allocatable :: count(:)
        complex(dp), allocatable :: kernel(:)
        integer(int64) :: work, ny
        logical :: ok
        logical, allocatable :: selected(:)
        character(:), allocatable :: which, how

        status = convolve_nd_invalid_input
        if (size(n)<1 .or. size(n)>15 .or. size(n)/=size(m)) return
        if (size(ydims)/=size(n)) return
        if (checked_product(n)/=size(a,kind=int64)) return
        if (checked_product(m)/=size(b,kind=int64)) return
        if (.not. all(ieee_is_finite(real(a,dp))) .or. &
            .not. all(ieee_is_finite(aimag(a)))) return
        if (.not. all(ieee_is_finite(real(b,dp))) .or. &
            .not. all(ieee_is_finite(aimag(b)))) return
        which = 'full'
        if (present(mode)) which = trim(mode)
        how = 'auto'
        if (present(method)) how = trim(method)
        if (how/='auto' .and. how/='direct' .and. how/='fft') return
        count = convolution_shape_nd(n,m,which,axes)
        if (size(count)/=size(n)) return
        if (any(ydims/=count)) return
        ny = checked_product(count)
        if (ny==0_int64) return
        allocate(selected(size(n)))
        call selected_axes(size(n),axes,selected,ok)
        if (.not. ok) return
        if (correlation) then
            allocate(kernel(size(b)))
            call reverse_conjugate_kernel(b,m,selected,kernel)
        else
            kernel = b
        end if
        if (how=='auto') then
            work = 1_int64
            ! Worst-case per-output kernel operations, guarded for overflow.
            block
                integer :: j
                do j=1,size(m)
                    if (.not. selected(j)) cycle
                    if (work>15000_int64/int(m(j),int64)) then
                        work = 15001_int64
                        exit
                    end if
                    work = work*int(m(j),int64)
                end do
            end block
            if (work>15000_int64/max(1_int64,ny)) then
                how = 'fft'
            else
                how = 'direct'
            end if
        end if
        select case(how)
        case ('direct')
            allocate(result(ny))
            call direct_flat(a,kernel,n,m,count,selected,which,result)
            status = convolve_nd_success
        case ('fft')
            call fftconvolve_flat(a,kernel,n,m,count,result,which,axes,status)
        end select
    end subroutine convolve_core

    subroutine convolve_nd_rr(in1,in2,y,mode,method,axes,status)
        real(dp), contiguous, target, intent(in) :: in1(..) !! first operand, rank 1..15
        real(dp), contiguous, target, intent(in) :: in2(..) !! second operand, same rank
        real(dp), contiguous, target, intent(inout) :: y(..) !! caller-owned output, unchanged on error
        character(*), intent(in), optional :: mode !! full (default), same, valid
        character(*), intent(in), optional :: method !! auto (default), direct, fft
        integer, intent(in), optional :: axes(:) !! selected 1-based axes; default all
        integer, intent(out), optional :: status !! zero on success
        real(dp), pointer :: a(:) => null()
        real(dp), pointer :: b(:) => null()
        real(dp), pointer :: out(:) => null()
        complex(dp), allocatable :: tmp(:)
        integer :: code

        if (present(status)) status = convolve_nd_invalid_input
        if (rank(in1)<1 .or. rank(in1)>15) return
        if (rank(in2)/=rank(in1) .or. rank(y)/=rank(in1)) return
        if (size(in1,kind=int64)==0_int64 .or. size(in2,kind=int64)==0_int64) return
        if (size(y,kind=int64)==0_int64) return
        call c_f_pointer(c_loc(in1),a,[size(in1,kind=int64)])
        call c_f_pointer(c_loc(in2),b,[size(in2,kind=int64)])
        call convolve_core(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp),shape(in1),shape(in2),shape(y), &
            tmp,mode,method,axes,.false.,code)
        if (code/=convolve_nd_success) return
        call c_f_pointer(c_loc(y),out,[size(y,kind=int64)])
        out = real(tmp,dp)
        if (present(status)) status = convolve_nd_success
    end subroutine convolve_nd_rr

    subroutine convolve_nd_cc(in1,in2,y,mode,method,axes,status)
        complex(dp), contiguous, target, intent(in) :: in1(..) !! first operand, rank 1..15
        complex(dp), contiguous, target, intent(in) :: in2(..) !! second operand, same rank
        complex(dp), contiguous, target, intent(inout) :: y(..) !! caller-owned output, unchanged on error
        character(*), intent(in), optional :: mode !! full (default), same, valid
        character(*), intent(in), optional :: method !! auto (default), direct, fft
        integer, intent(in), optional :: axes(:) !! selected 1-based axes; default all
        integer, intent(out), optional :: status !! zero on success
        complex(dp), pointer :: a(:) => null()
        complex(dp), pointer :: b(:) => null()
        complex(dp), pointer :: out(:) => null()
        complex(dp), allocatable :: tmp(:)
        integer :: code

        if (present(status)) status = convolve_nd_invalid_input
        if (rank(in1)<1 .or. rank(in1)>15) return
        if (rank(in2)/=rank(in1) .or. rank(y)/=rank(in1)) return
        if (size(in1,kind=int64)==0_int64 .or. size(in2,kind=int64)==0_int64) return
        if (size(y,kind=int64)==0_int64) return
        call c_f_pointer(c_loc(in1),a,[size(in1,kind=int64)])
        call c_f_pointer(c_loc(in2),b,[size(in2,kind=int64)])
        call convolve_core(a,b,shape(in1),shape(in2),shape(y), &
            tmp,mode,method,axes,.false.,code)
        if (code/=convolve_nd_success) return
        call c_f_pointer(c_loc(y),out,[size(y,kind=int64)])
        out = tmp
        if (present(status)) status = convolve_nd_success
    end subroutine convolve_nd_cc

    subroutine convolve_nd_rc(in1,in2,y,mode,method,axes,status)
        real(dp), contiguous, target, intent(in) :: in1(..) !! first operand, rank 1..15
        complex(dp), contiguous, target, intent(in) :: in2(..) !! second operand, same rank
        complex(dp), contiguous, target, intent(inout) :: y(..) !! caller-owned output, unchanged on error
        character(*), intent(in), optional :: mode !! full (default), same, valid
        character(*), intent(in), optional :: method !! auto (default), direct, fft
        integer, intent(in), optional :: axes(:) !! selected 1-based axes; default all
        integer, intent(out), optional :: status !! zero on success
        real(dp), pointer :: a(:) => null()
        complex(dp), pointer :: b(:) => null()
        complex(dp), pointer :: out(:) => null()
        complex(dp), allocatable :: tmp(:)
        integer :: code

        if (present(status)) status = convolve_nd_invalid_input
        if (rank(in1)<1 .or. rank(in1)>15) return
        if (rank(in2)/=rank(in1) .or. rank(y)/=rank(in1)) return
        if (size(in1,kind=int64)==0_int64 .or. size(in2,kind=int64)==0_int64) return
        if (size(y,kind=int64)==0_int64) return
        call c_f_pointer(c_loc(in1),a,[size(in1,kind=int64)])
        call c_f_pointer(c_loc(in2),b,[size(in2,kind=int64)])
        call convolve_core(cmplx(a,0.0_dp,dp),b,shape(in1),shape(in2),shape(y), &
            tmp,mode,method,axes,.false.,code)
        if (code/=convolve_nd_success) return
        call c_f_pointer(c_loc(y),out,[size(y,kind=int64)])
        out = tmp
        if (present(status)) status = convolve_nd_success
    end subroutine convolve_nd_rc

    subroutine convolve_nd_cr(in1,in2,y,mode,method,axes,status)
        complex(dp), contiguous, target, intent(in) :: in1(..) !! first operand, rank 1..15
        real(dp), contiguous, target, intent(in) :: in2(..) !! second operand, same rank
        complex(dp), contiguous, target, intent(inout) :: y(..) !! caller-owned output, unchanged on error
        character(*), intent(in), optional :: mode !! full (default), same, valid
        character(*), intent(in), optional :: method !! auto (default), direct, fft
        integer, intent(in), optional :: axes(:) !! selected 1-based axes; default all
        integer, intent(out), optional :: status !! zero on success
        complex(dp), pointer :: a(:) => null()
        real(dp), pointer :: b(:) => null()
        complex(dp), pointer :: out(:) => null()
        complex(dp), allocatable :: tmp(:)
        integer :: code

        if (present(status)) status = convolve_nd_invalid_input
        if (rank(in1)<1 .or. rank(in1)>15) return
        if (rank(in2)/=rank(in1) .or. rank(y)/=rank(in1)) return
        if (size(in1,kind=int64)==0_int64 .or. size(in2,kind=int64)==0_int64) return
        if (size(y,kind=int64)==0_int64) return
        call c_f_pointer(c_loc(in1),a,[size(in1,kind=int64)])
        call c_f_pointer(c_loc(in2),b,[size(in2,kind=int64)])
        call convolve_core(a,cmplx(b,0.0_dp,dp),shape(in1),shape(in2),shape(y), &
            tmp,mode,method,axes,.false.,code)
        if (code/=convolve_nd_success) return
        call c_f_pointer(c_loc(y),out,[size(y,kind=int64)])
        out = tmp
        if (present(status)) status = convolve_nd_success
    end subroutine convolve_nd_cr

    subroutine correlate_nd_rr(in1,in2,y,mode,method,axes,status)
        real(dp), contiguous, target, intent(in) :: in1(..) !! first operand, rank 1..15
        real(dp), contiguous, target, intent(in) :: in2(..) !! second operand, same rank
        real(dp), contiguous, target, intent(inout) :: y(..) !! caller-owned output, unchanged on error
        character(*), intent(in), optional :: mode !! full (default), same, valid
        character(*), intent(in), optional :: method !! auto (default), direct, fft
        integer, intent(in), optional :: axes(:) !! selected 1-based axes; default all
        integer, intent(out), optional :: status !! zero on success
        real(dp), pointer :: a(:) => null()
        real(dp), pointer :: b(:) => null()
        real(dp), pointer :: out(:) => null()
        complex(dp), allocatable :: tmp(:)
        integer :: code

        if (present(status)) status = convolve_nd_invalid_input
        if (rank(in1)<1 .or. rank(in1)>15) return
        if (rank(in2)/=rank(in1) .or. rank(y)/=rank(in1)) return
        if (size(in1,kind=int64)==0_int64 .or. size(in2,kind=int64)==0_int64) return
        if (size(y,kind=int64)==0_int64) return
        call c_f_pointer(c_loc(in1),a,[size(in1,kind=int64)])
        call c_f_pointer(c_loc(in2),b,[size(in2,kind=int64)])
        call convolve_core(cmplx(a,0.0_dp,dp),cmplx(b,0.0_dp,dp),shape(in1),shape(in2),shape(y), &
            tmp,mode,method,axes,.true.,code)
        if (code/=convolve_nd_success) return
        call c_f_pointer(c_loc(y),out,[size(y,kind=int64)])
        out = real(tmp,dp)
        if (present(status)) status = convolve_nd_success
    end subroutine correlate_nd_rr

    subroutine correlate_nd_cc(in1,in2,y,mode,method,axes,status)
        complex(dp), contiguous, target, intent(in) :: in1(..) !! first operand, rank 1..15
        complex(dp), contiguous, target, intent(in) :: in2(..) !! second operand, same rank
        complex(dp), contiguous, target, intent(inout) :: y(..) !! caller-owned output, unchanged on error
        character(*), intent(in), optional :: mode !! full (default), same, valid
        character(*), intent(in), optional :: method !! auto (default), direct, fft
        integer, intent(in), optional :: axes(:) !! selected 1-based axes; default all
        integer, intent(out), optional :: status !! zero on success
        complex(dp), pointer :: a(:) => null()
        complex(dp), pointer :: b(:) => null()
        complex(dp), pointer :: out(:) => null()
        complex(dp), allocatable :: tmp(:)
        integer :: code

        if (present(status)) status = convolve_nd_invalid_input
        if (rank(in1)<1 .or. rank(in1)>15) return
        if (rank(in2)/=rank(in1) .or. rank(y)/=rank(in1)) return
        if (size(in1,kind=int64)==0_int64 .or. size(in2,kind=int64)==0_int64) return
        if (size(y,kind=int64)==0_int64) return
        call c_f_pointer(c_loc(in1),a,[size(in1,kind=int64)])
        call c_f_pointer(c_loc(in2),b,[size(in2,kind=int64)])
        call convolve_core(a,b,shape(in1),shape(in2),shape(y), &
            tmp,mode,method,axes,.true.,code)
        if (code/=convolve_nd_success) return
        call c_f_pointer(c_loc(y),out,[size(y,kind=int64)])
        out = tmp
        if (present(status)) status = convolve_nd_success
    end subroutine correlate_nd_cc

    subroutine correlate_nd_rc(in1,in2,y,mode,method,axes,status)
        real(dp), contiguous, target, intent(in) :: in1(..) !! first operand, rank 1..15
        complex(dp), contiguous, target, intent(in) :: in2(..) !! second operand, same rank
        complex(dp), contiguous, target, intent(inout) :: y(..) !! caller-owned output, unchanged on error
        character(*), intent(in), optional :: mode !! full (default), same, valid
        character(*), intent(in), optional :: method !! auto (default), direct, fft
        integer, intent(in), optional :: axes(:) !! selected 1-based axes; default all
        integer, intent(out), optional :: status !! zero on success
        real(dp), pointer :: a(:) => null()
        complex(dp), pointer :: b(:) => null()
        complex(dp), pointer :: out(:) => null()
        complex(dp), allocatable :: tmp(:)
        integer :: code

        if (present(status)) status = convolve_nd_invalid_input
        if (rank(in1)<1 .or. rank(in1)>15) return
        if (rank(in2)/=rank(in1) .or. rank(y)/=rank(in1)) return
        if (size(in1,kind=int64)==0_int64 .or. size(in2,kind=int64)==0_int64) return
        if (size(y,kind=int64)==0_int64) return
        call c_f_pointer(c_loc(in1),a,[size(in1,kind=int64)])
        call c_f_pointer(c_loc(in2),b,[size(in2,kind=int64)])
        call convolve_core(cmplx(a,0.0_dp,dp),b,shape(in1),shape(in2),shape(y), &
            tmp,mode,method,axes,.true.,code)
        if (code/=convolve_nd_success) return
        call c_f_pointer(c_loc(y),out,[size(y,kind=int64)])
        out = tmp
        if (present(status)) status = convolve_nd_success
    end subroutine correlate_nd_rc

    subroutine correlate_nd_cr(in1,in2,y,mode,method,axes,status)
        complex(dp), contiguous, target, intent(in) :: in1(..) !! first operand, rank 1..15
        real(dp), contiguous, target, intent(in) :: in2(..) !! second operand, same rank
        complex(dp), contiguous, target, intent(inout) :: y(..) !! caller-owned output, unchanged on error
        character(*), intent(in), optional :: mode !! full (default), same, valid
        character(*), intent(in), optional :: method !! auto (default), direct, fft
        integer, intent(in), optional :: axes(:) !! selected 1-based axes; default all
        integer, intent(out), optional :: status !! zero on success
        complex(dp), pointer :: a(:) => null()
        real(dp), pointer :: b(:) => null()
        complex(dp), pointer :: out(:) => null()
        complex(dp), allocatable :: tmp(:)
        integer :: code

        if (present(status)) status = convolve_nd_invalid_input
        if (rank(in1)<1 .or. rank(in1)>15) return
        if (rank(in2)/=rank(in1) .or. rank(y)/=rank(in1)) return
        if (size(in1,kind=int64)==0_int64 .or. size(in2,kind=int64)==0_int64) return
        if (size(y,kind=int64)==0_int64) return
        call c_f_pointer(c_loc(in1),a,[size(in1,kind=int64)])
        call c_f_pointer(c_loc(in2),b,[size(in2,kind=int64)])
        call convolve_core(a,cmplx(b,0.0_dp,dp),shape(in1),shape(in2),shape(y), &
            tmp,mode,method,axes,.true.,code)
        if (code/=convolve_nd_success) return
        call c_f_pointer(c_loc(y),out,[size(y,kind=int64)])
        out = tmp
        if (present(status)) status = convolve_nd_success
    end subroutine correlate_nd_cr

end module scifort_signal_convolve_nd
