! SPDX-License-Identifier: MIT
! Copyright (c) 2026 SciFort contributors
!
! Stable barycentric polynomial interpolation on distinct real nodes.
! The weights are normalized using logarithms to avoid products of large
! or small node distances. No source from SciPy or other projects is used.
module scifort_interpolate_barycentric
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use scifort_kinds, only : dp
    implicit none
    private

    integer, parameter, public :: barycentric_success = 0
    integer, parameter, public :: barycentric_invalid_input = 1

    type, public :: barycentric_interpolator_type
        private
        real(dp), allocatable :: x(:), w(:), y(:,:)
    contains
        procedure :: valid => real_valid
        procedure :: nodes => real_nodes
        procedure :: weights => real_weights
        procedure :: channels => real_channels
        procedure :: value => real_value
        procedure :: values => real_values
        procedure :: derivative => real_derivative
        procedure :: set_values => real_set_values
        procedure :: add_nodes => real_add_nodes
    end type barycentric_interpolator_type

    type, public :: complex_barycentric_interpolator_type
        private
        real(dp), allocatable :: x(:), w(:)
        complex(dp), allocatable :: y(:,:)
    contains
        procedure :: valid => complex_valid
        procedure :: nodes => complex_nodes
        procedure :: weights => complex_weights
        procedure :: channels => complex_channels
        procedure :: value => complex_value
        procedure :: values => complex_values
        procedure :: derivative => complex_derivative
        procedure :: set_values => complex_set_values
        procedure :: add_nodes => complex_add_nodes
    end type complex_barycentric_interpolator_type

    public :: barycentric_interpolator, barycentric_interpolate
    interface barycentric_interpolator
        module procedure create_real_vector
        module procedure create_real_matrix
        module procedure create_complex_vector
        module procedure create_complex_matrix
    end interface
    interface barycentric_interpolate
        module procedure interpolate_real_vector
        module procedure interpolate_complex_vector
    end interface
contains

    subroutine calculate_weights(x, w, status)
        real(dp), intent(in) :: x(:) !! interpolation nodes in arbitrary order
        real(dp), allocatable, intent(out) :: w(:) !! normalized nonzero barycentric weights
        integer, intent(out) :: status !! zero on success, nonzero on invalid input
        real(dp), allocatable :: logs(:)
        real(dp) :: d, xmax, scale, sign_w
        integer :: i, j, n
        status = barycentric_invalid_input
        n = size(x)
        if (n < 1) return
        if (.not. all(ieee_is_finite(x))) return
        allocate(logs(n), w(n))
        w = 0.0_dp
        if (n == 1) then
            w = 1.0_dp
            status = barycentric_success
            return
        end if
        scale = maxval(x)-minval(x)
        if (scale <= 0.0_dp .or. .not. ieee_is_finite(scale)) return
        do i = 1, n
            sign_w = 1.0_dp
            logs(i) = 0.0_dp
            do j = 1, n
                if (i == j) cycle
                d = (x(i)-x(j))/scale
                if (d == 0.0_dp) return
                if (d < 0.0_dp) sign_w = -sign_w
                logs(i) = logs(i)-log(abs(d))
            end do
            w(i) = sign_w
        end do
        xmax = maxval(logs)
        w = w*exp(logs-xmax)
        if (any(w == 0.0_dp)) return
        status = barycentric_success
    end subroutine calculate_weights

    function create_real_vector(x, y, status) result(p)
        real(dp), intent(in) :: x(:) !! distinct interpolation nodes
        real(dp), intent(in) :: y(:) !! ordinates at each node
        integer, intent(out), optional :: status !! zero on valid interpolator
        type(barycentric_interpolator_type) :: p
        integer :: code
        code = barycentric_invalid_input
        if (size(x) == size(y) .and. all(ieee_is_finite(y))) then
            call calculate_weights(x,p%w,code)
            if (code == barycentric_success) then
                p%x = x
                p%y = reshape(y,[size(y),1])
            end if
        end if
        if (present(status)) status = code
    end function create_real_vector

    function create_real_matrix(x, y, status) result(p)
        real(dp), intent(in) :: x(:) !! interpolation nodes
        real(dp), intent(in) :: y(:,:) !! ordinates (nodes,channels)
        integer, intent(out), optional :: status !! zero on valid interpolator
        type(barycentric_interpolator_type) :: p
        integer :: code
        code = barycentric_invalid_input
        if (size(x) == size(y,1) .and. size(y,2)>0 .and. all(ieee_is_finite(y))) then
            call calculate_weights(x,p%w,code)
            if (code == barycentric_success) then
                p%x = x
                p%y = y
            end if
        end if
        if (present(status)) status = code
    end function create_real_matrix

    function create_complex_vector(x, y, status) result(p)
        real(dp), intent(in) :: x(:) !! distinct real interpolation nodes
        complex(dp), intent(in) :: y(:) !! complex ordinates
        integer, intent(out), optional :: status !! zero on success
        type(complex_barycentric_interpolator_type) :: p
        integer :: code
        code = barycentric_invalid_input
        if (size(x) == size(y) .and. finite_complex(y)) then
            call calculate_weights(x,p%w,code)
            if (code == barycentric_success) then
                p%x = x
                p%y = reshape(y,[size(y),1])
            end if
        end if
        if (present(status)) status = code
    end function create_complex_vector

    function create_complex_matrix(x, y, status) result(p)
        real(dp), intent(in) :: x(:) !! distinct real interpolation nodes
        complex(dp), intent(in) :: y(:,:) !! complex ordinates (nodes,channels)
        integer, intent(out), optional :: status !! zero on success
        type(complex_barycentric_interpolator_type) :: p
        integer :: code
        code = barycentric_invalid_input
        if (size(x) == size(y,1) .and. size(y,2)>0 .and. finite_complex(reshape(y,[size(y)]))) then
            call calculate_weights(x,p%w,code)
            if (code == barycentric_success) then
                p%x = x
                p%y = y
            end if
        end if
        if (present(status)) status = code
    end function create_complex_matrix

    pure logical function finite_complex(y)
        complex(dp), intent(in) :: y(:) !! values to validate
        finite_complex = all(ieee_is_finite(real(y,dp))) .and. all(ieee_is_finite(aimag(y)))
    end function finite_complex

    pure logical function real_valid(self)
        class(barycentric_interpolator_type), intent(in) :: self !! interpolator to validate
        real_valid = allocated(self%x) .and. allocated(self%w) .and. allocated(self%y)
    end function real_valid

    pure logical function complex_valid(self)
        class(complex_barycentric_interpolator_type), intent(in) :: self !! interpolator to validate
        complex_valid = allocated(self%x) .and. allocated(self%w) .and. allocated(self%y)
    end function complex_valid

    function real_nodes(self) result(x)
        class(barycentric_interpolator_type), intent(in) :: self !! source interpolator
        real(dp), allocatable :: x(:)
        if (self%valid()) x = self%x
    end function real_nodes

    function complex_nodes(self) result(x)
        class(complex_barycentric_interpolator_type), intent(in) :: self !! source interpolator
        real(dp), allocatable :: x(:)
        if (self%valid()) x = self%x
    end function complex_nodes

    function real_weights(self) result(w)
        class(barycentric_interpolator_type), intent(in) :: self !! source interpolator
        real(dp), allocatable :: w(:)
        if (self%valid()) w = self%w
    end function real_weights

    function complex_weights(self) result(w)
        class(complex_barycentric_interpolator_type), intent(in) :: self !! source interpolator
        real(dp), allocatable :: w(:)
        if (self%valid()) w = self%w
    end function complex_weights

    pure integer function real_channels(self)
        class(barycentric_interpolator_type), intent(in) :: self !! source interpolator
        real_channels = 0
        if (self%valid()) real_channels = size(self%y,2)
    end function real_channels

    pure integer function complex_channels(self)
        class(complex_barycentric_interpolator_type), intent(in) :: self !! source interpolator
        complex_channels = 0
        if (self%valid()) complex_channels = size(self%y,2)
    end function complex_channels

    ! When a query matches a node, the exact stored ordinate is returned.
    ! Otherwise normalize ratios to avoid large intermediate 1/(q-xi).
    function interpolate_real_at(x,w,y,q) result(v)
        real(dp), intent(in) :: x(:) !! interpolation nodes
        real(dp), intent(in) :: w(:) !! normalized barycentric weights
        real(dp), intent(in) :: y(:,:) !! ordinates (nodes,channels)
        real(dp), intent(in) :: q !! evaluation coordinate
        real(dp) :: v(size(y,2))
        real(dp) :: terms(size(x)), scale, den
        integer :: i
        v = 0.0_dp
        if (size(x) == 0) return
        do i=1,size(x)
            if (q == x(i)) then
                v = y(i,:)
                return
            end if
        end do
        terms = w/((q-x)/max(1.0_dp,abs(q),maxval(abs(x))))
        scale = maxval(abs(terms))
        if (scale == 0.0_dp) return
        terms = terms/scale
        den = sum(terms)
        if (den == 0.0_dp) return
        do i=1,size(x)
            v = v + terms(i)*y(i,:)
        end do
        v = v/den
    end function interpolate_real_at

    function interpolate_complex_at(x,w,y,q) result(v)
        real(dp), intent(in) :: x(:) !! interpolation nodes
        real(dp), intent(in) :: w(:) !! normalized barycentric weights
        complex(dp), intent(in) :: y(:,:) !! complex ordinates (nodes,channels)
        real(dp), intent(in) :: q !! evaluation coordinate
        complex(dp) :: v(size(y,2))
        real(dp) :: terms(size(x)), scale, den
        integer :: i
        v = (0.0_dp,0.0_dp)
        if (size(x) == 0) return
        do i=1,size(x)
            if (q == x(i)) then
                v = y(i,:)
                return
            end if
        end do
        terms = w/((q-x)/max(1.0_dp,abs(q),maxval(abs(x))))
        scale = maxval(abs(terms))
        if (scale == 0.0_dp) return
        terms = terms/scale
        den = sum(terms)
        if (den == 0.0_dp) return
        do i=1,size(x)
            v = v + terms(i)*y(i,:)
        end do
        v = v/den
    end function interpolate_complex_at

    ! Differentiate the interpolation polynomial at the interpolation nodes
    ! using the barycentric spectral differentiation matrix. Successive
    ! multiplication yields derivatives of any order below the node count.
    function real_node_derivatives(x,w,y,nu) result(d)
        real(dp), intent(in) :: x(:) !! distinct nodes
        real(dp), intent(in) :: w(:) !! nonzero scaled weights
        real(dp), intent(in) :: y(:,:) !! ordinates
        integer, intent(in) :: nu !! nonnegative derivative order
        real(dp) :: d(size(y,1),size(y,2))
        real(dp), allocatable :: mat(:,:), next(:,:)
        integer :: i,j,k,n
        n = size(x)
        d = y
        if (nu == 0) return
        if (nu >= n) then
            d = 0.0_dp
            return
        end if
        allocate(mat(n,n))
        mat = 0.0_dp
        do i=1,n
            do j=1,n
                if (i==j) cycle
                mat(i,j) = (w(j)/w(i))/(x(i)-x(j))
                mat(i,i) = mat(i,i)-mat(i,j)
            end do
        end do
        do k=1,nu
            next = matmul(mat,d)
            d = next
        end do
    end function real_node_derivatives

    function complex_node_derivatives(x,w,y,nu) result(d)
        real(dp), intent(in) :: x(:) !! distinct nodes
        real(dp), intent(in) :: w(:) !! nonzero scaled weights
        complex(dp), intent(in) :: y(:,:) !! complex ordinates
        integer, intent(in) :: nu !! nonnegative derivative order
        complex(dp) :: d(size(y,1),size(y,2))
        real(dp), allocatable :: mat(:,:)
        complex(dp), allocatable :: next(:,:)
        integer :: i,j,k,n
        n = size(x)
        d = y
        if (nu == 0) return
        if (nu >= n) then
            d = (0.0_dp,0.0_dp)
            return
        end if
        allocate(mat(n,n))
        mat = 0.0_dp
        do i=1,n
            do j=1,n
                if (i==j) cycle
                mat(i,j) = (w(j)/w(i))/(x(i)-x(j))
                mat(i,i) = mat(i,i)-mat(i,j)
            end do
        end do
        do k=1,nu
            next = matmul(mat,d)
            d = next
        end do
    end function complex_node_derivatives

    function real_values(self,q,nu) result(v)
        class(barycentric_interpolator_type), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: q(:) !! query coordinates
        integer, intent(in), optional :: nu !! derivative order, default zero
        real(dp), allocatable :: v(:,:)
        real(dp), allocatable :: d(:,:)
        integer :: j,order
        if (.not. self%valid()) then
            allocate(v(0,0))
            return
        end if
        allocate(v(size(q),self%channels()))
        v = 0.0_dp
        order = 0
        if (present(nu)) order = nu
        if (order<0 .or. order>=size(self%x)) return
        if (order==0) then
            do j=1,size(q)
                v(j,:) = interpolate_real_at(self%x,self%w,self%y,q(j))
            end do
        else
            d = real_node_derivatives(self%x,self%w,self%y,order)
            do j=1,size(q)
                v(j,:) = interpolate_real_at(self%x,self%w,d,q(j))
            end do
        end if
    end function real_values

    function complex_values(self,q,nu) result(v)
        class(complex_barycentric_interpolator_type), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: q(:) !! query coordinates
        integer, intent(in), optional :: nu !! derivative order, default zero
        complex(dp), allocatable :: v(:,:)
        complex(dp), allocatable :: d(:,:)
        integer :: j,order
        if (.not. self%valid()) then
            allocate(v(0,0))
            return
        end if
        allocate(v(size(q),self%channels()))
        v = (0.0_dp,0.0_dp)
        order = 0
        if (present(nu)) order = nu
        if (order<0 .or. order>=size(self%x)) return
        if (order==0) then
            do j=1,size(q)
                v(j,:) = interpolate_complex_at(self%x,self%w,self%y,q(j))
            end do
        else
            d = complex_node_derivatives(self%x,self%w,self%y,order)
            do j=1,size(q)
                v(j,:) = interpolate_complex_at(self%x,self%w,d,q(j))
            end do
        end if
    end function complex_values

    real(dp) function real_value(self,q,channel,nu) result(v)
        class(barycentric_interpolator_type), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: q !! query coordinate
        integer, intent(in), optional :: channel !! selected channel, 1-based
        integer, intent(in), optional :: nu !! derivative order, default zero
        real(dp), allocatable :: a(:,:)
        integer :: ch
        v = 0.0_dp
        if (.not. self%valid()) return
        ch=1
        if (present(channel)) ch=channel
        if (ch<1 .or. ch>self%channels()) return
        a=self%values([q],nu)
        v=a(1,ch)
    end function real_value

    complex(dp) function complex_value(self,q,channel,nu) result(v)
        class(complex_barycentric_interpolator_type), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: q !! query coordinate
        integer, intent(in), optional :: channel !! selected channel, 1-based
        integer, intent(in), optional :: nu !! derivative order, default zero
        complex(dp), allocatable :: a(:,:)
        integer :: ch
        v = (0.0_dp,0.0_dp)
        if (.not. self%valid()) return
        ch=1
        if (present(channel)) ch=channel
        if (ch<1 .or. ch>self%channels()) return
        a=self%values([q],nu)
        v=a(1,ch)
    end function complex_value

    real(dp) function real_derivative(self,q,nu,channel) result(v)
        class(barycentric_interpolator_type), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: q !! evaluation coordinate
        integer, intent(in), optional :: nu !! derivative order, default first
        integer, intent(in), optional :: channel !! output channel, 1-based
        integer :: order
        order=1
        if (present(nu)) order=nu
        v=self%value(q,channel,order)
    end function real_derivative

    complex(dp) function complex_derivative(self,q,nu,channel) result(v)
        class(complex_barycentric_interpolator_type), intent(in) :: self !! fitted interpolator
        real(dp), intent(in) :: q !! evaluation coordinate
        integer, intent(in), optional :: nu !! derivative order, default first
        integer, intent(in), optional :: channel !! output channel, 1-based
        integer :: order
        order=1
        if (present(nu)) order=nu
        v=self%value(q,channel,order)
    end function complex_derivative

    subroutine real_set_values(self,y,status)
        class(barycentric_interpolator_type), intent(inout) :: self !! interpolator to update
        real(dp), intent(in) :: y(:,:) !! new values (nodes,channels)
        integer, intent(out), optional :: status !! zero if updated
        integer :: code
        code=barycentric_invalid_input
        if (allocated(self%x)) then
            if (size(y,1)==size(self%x) .and. size(y,2)>0 .and. all(ieee_is_finite(y))) then
                self%y=y
                code=barycentric_success
            end if
        end if
        if (present(status)) status=code
    end subroutine real_set_values

    subroutine complex_set_values(self,y,status)
        class(complex_barycentric_interpolator_type), intent(inout) :: self !! interpolator to update
        complex(dp), intent(in) :: y(:,:) !! new complex values (nodes,channels)
        integer, intent(out), optional :: status !! zero if updated
        integer :: code
        code=barycentric_invalid_input
        if (allocated(self%x)) then
            if (size(y,1)==size(self%x) .and. size(y,2)>0 .and. finite_complex(reshape(y,[size(y)]))) then
                self%y=y
                code=barycentric_success
            end if
        end if
        if (present(status)) status=code
    end subroutine complex_set_values

    subroutine real_add_nodes(self,xnew,ynew,status)
        class(barycentric_interpolator_type), intent(inout) :: self !! interpolator to update
        real(dp), intent(in) :: xnew(:) !! distinct additional nodes
        real(dp), intent(in) :: ynew(:,:) !! values (additional nodes,channels)
        integer, intent(out), optional :: status !! zero on successful transactional insertion
        real(dp), allocatable :: combined(:), weights(:), new_y(:,:) 
        integer :: code
        code=barycentric_invalid_input
        if (self%valid() .and. size(xnew)>0) then
            if (size(ynew,1)==size(xnew) .and. size(ynew,2)==self%channels()) then
                if (all(ieee_is_finite(ynew))) then
                    combined=[self%x,xnew]
                    call calculate_weights(combined,weights,code)
                    if (code==barycentric_success) then
                        self%x=combined
                        self%w=weights
                        allocate(new_y(size(combined),self%channels()))
                        new_y(:size(self%y,1),:) = self%y
                        new_y(size(self%y,1)+1:,:) = ynew
                        call move_alloc(new_y,self%y)
                    end if
                end if
            end if
        end if
        if (present(status)) status=code
    end subroutine real_add_nodes

    subroutine complex_add_nodes(self,xnew,ynew,status)
        class(complex_barycentric_interpolator_type), intent(inout) :: self !! interpolator to update
        real(dp), intent(in) :: xnew(:) !! distinct additional nodes
        complex(dp), intent(in) :: ynew(:,:) !! new complex values (additional nodes,channels)
        integer, intent(out), optional :: status !! zero on successful transactional insertion
        real(dp), allocatable :: combined(:), weights(:)
        complex(dp), allocatable :: new_y(:,:)
        integer :: code
        code=barycentric_invalid_input
        if (self%valid() .and. size(xnew)>0) then
            if (size(ynew,1)==size(xnew) .and. size(ynew,2)==self%channels()) then
                if (finite_complex(reshape(ynew,[size(ynew)]))) then
                    combined=[self%x,xnew]
                    call calculate_weights(combined,weights,code)
                    if (code==barycentric_success) then
                        self%x=combined
                        self%w=weights
                        allocate(new_y(size(combined),self%channels()))
                        new_y(:size(self%y,1),:) = self%y
                        new_y(size(self%y,1)+1:,:) = ynew
                        call move_alloc(new_y,self%y)
                    end if
                end if
            end if
        end if
        if (present(status)) status=code
    end subroutine complex_add_nodes

    function interpolate_real_vector(x,y,q,nu,status) result(v)
        real(dp), intent(in) :: x(:) !! distinct interpolation nodes
        real(dp), intent(in) :: y(:) !! nodal ordinates
        real(dp), intent(in) :: q(:) !! evaluation points
        integer, intent(in), optional :: nu !! derivative order (default zero)
        integer, intent(out), optional :: status !! zero if construction succeeds
        real(dp), allocatable :: v(:)
        type(barycentric_interpolator_type) :: p
        real(dp), allocatable :: tmp(:,:)
        integer :: code
        p=barycentric_interpolator(x,y,status=code)
        if (p%valid()) then
            tmp=p%values(q,nu)
            v=tmp(:,1)
        else
            allocate(v(0))
        end if
        if (present(status)) status=code
    end function interpolate_real_vector

    function interpolate_complex_vector(x,y,q,nu,status) result(v)
        real(dp), intent(in) :: x(:) !! distinct interpolation nodes
        complex(dp), intent(in) :: y(:) !! complex nodal ordinates
        real(dp), intent(in) :: q(:) !! evaluation points
        integer, intent(in), optional :: nu !! derivative order (default zero)
        integer, intent(out), optional :: status !! zero if construction succeeds
        complex(dp), allocatable :: v(:)
        type(complex_barycentric_interpolator_type) :: p
        complex(dp), allocatable :: tmp(:,:)
        integer :: code
        p=barycentric_interpolator(x,y,status=code)
        if (p%valid()) then
            tmp=p%values(q,nu)
            v=tmp(:,1)
        else
            allocate(v(0))
        end if
        if (present(status)) status=code
    end function interpolate_complex_vector
end module scifort_interpolate_barycentric
