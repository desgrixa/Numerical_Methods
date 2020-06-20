module simpleinteg
    contains
        subroutine integ_polynomials(coeff, lowrange, highrange, integval)
            use types
            real(kind=rkind), dimension(:), intent(in) :: coeff
            real(kind=rkind), intent(in) :: lowrange, highrange
            real(kind=rkind) :: integval
            integer :: i

            !Series which will evaluate
            real(kind = rkind) :: F_low, F_high !values of primitive function for low and high ranges

            F_low = 0

            do i=1, ubound(coeff, 1) ! we don;t exactly know the upper boundary
                !it's always coeff to the
                F_low = F_low + coeff(i)/i * lowrange**i
            end do

            F_high = 0

            do i=1, ubound(coeff, 1)
                F_high = F_high + coeff(i)/i * highrange**i
            end do

            integval = F_high - F_low

        end subroutine integ_polynomials

    function collocation(f,alpha,beta, n) result(integval)
        use types
        interface
            function f(x) result(y)
                use types
                real(kind = rkind), intent(in) :: x
                real(kind = rkind) :: y
            end function f
        end interface

        real(kind = rkind), intent(in) :: alpha, beta
        real(kind = rkind) :: integval
        ! number of integration nodes
        integer :: n
        real(kind = rkind), dimension(:,:), allocatable :: mtx
        real(kind = rkind), dimension(:), allocatable :: bside, xvect, coeff

        !most important step to have coeff
        call gen_collocation(f, alpha, beta, n, mtx, bside, xvect, coeff)
        !
        !all integ_polynomials(coeff, alpha, beta, integval)

        !print *, "function with known exact solution: ", integval

    end function collocation

    subroutine gen_collocation(f, alpha, beta, n, mtx,bside, xvect, coeff)
    !xvect - (vector of x)
        use types
        use linalg
        interface
            function f(x) result(y)
                use types
                real(kind = rkind), intent(in) :: x
                real(kind = rkind) :: y
            end function f
        end interface

        real(kind = rkind), intent(in) :: alpha, beta
        real(kind = rkind), dimension(:,:), allocatable, intent(out) :: mtx
        real(kind = rkind), dimension(:), allocatable, intent(out) :: bside, xvect, coeff
        ! degree of polynomial: (n-1); n - (number of internal points)
        integer, intent(in) :: n
        integer :: i, j

        allocate(xvect(n-2)) !n-2 is a dimension
         !xvect and bside have a dimension given by number of n
        allocate(bside(n)) !dimension refers to number of coeff
        allocate(mtx(n,n)) !matrix for the given dimension

        do i=1, ubound(xvect,1) ! loop until dimension of xvect (to the range of xvect (ubound x, n-2) )
            xvect(i) = (beta-alpha)/(n-1) * i + alpha
        end do

        ! u(x) = α exp(2x) + β exp(−2x) + 1/13*sin 3x,

        print *, xvect !
        ! now the matrix is allocated

        mtx(1:n,1) = 1.0_rkind ! (1:n, 1) -> means row number 1 till n, and first col
    !rkind - just to be accurate enough
        mtx(3:n,1) = 4.0_rkind
        do i=2, n ! from the second col till n
            mtx(1,i) = alpha**(i-1) ! 1 row till i, equal to
        end do



        do i=2, n ! for the second row in the matrix
            ! mtx(2,i) = exp(-2.0) + sin(3.0)/13
            mtx(2,i) = beta**(i-1)
        end do

        do j=3, n ! rest of the matrix (row 3, 4, 5)
            do i=2, n
                mtx(j,i) = -((i-2)*(i-1)*(xvect(j-2))**(i-3)) + 4*(xvect(j-2)**(i-1))
            end do
        end do
        ! procedure to bside
        bside(1) = f(alpha)
        bside(2) = f(beta)
        ! print *, bside(2) , "2"
        do i=1, n-2 ! beacuse first to are reserved
            bside(i+2) = f(xvect(i))
        end do

        ! print *, bside
        call printmatrix(mtx)
        call gem(mtx, coeff, bside, xvect) ! when we call gausiion elemination the
        !coeff will be a vector of coeff of our polynomial approximation
        !so then we will have our procedure ready.
        ! this gives us ability to integrate functions
    ! matrix generated
    end subroutine gen_collocation

    function finite_method(f, lowrange, highrange, deltax) result(integval)
        use types
        interface
            function f(x) result (y)
                use types
                real(kind =rkind), intent(in) :: x
                real(kind =rkind)  :: y
            end function f
        end interface

        real(kind =rkind), intent(in) :: lowrange, highrange, deltax
        real(kind =rkind) :: integval
        real(kind =rkind), dimension(:,:), allocatable:: mtx
        real(kind =rkind), dimension(:),  allocatable :: bside, xvect, coeff


        call FDM(f, lowrange, highrange, deltax, mtx,  bside, xvect, coeff)
        ! call gen_coeffmat(f, lowrange, highrange, n, mtx,  bside, xvect, coeff)

    end function finite_method

    subroutine FDM(f, lowrange, highrange, deltax, mtx, bside, xvect, coeff)
        use types
        use linalg
        interface
            function f(x) result (y)
                use types
                real(kind =rkind), intent(in) :: x
                real(kind =rkind)  :: y
            end function f
        end interface

        real(kind =rkind), intent(in) :: lowrange, highrange
        real(kind =rkind), dimension(:,:), allocatable, intent (out) :: mtx
        real(kind =rkind), dimension(:), allocatable, intent (out) :: bside, xvect, coeff
        real(kind =rkind) :: deltax
        integer :: i
        integer :: n

        n = ((highrange-lowrange)/deltax)+1

        allocate(xvect(n))
        allocate(bside(n))
        allocate(mtx(n,n))

        xvect(1) = lowrange

        do i=2, n-1
            xvect(i) = xvect(i-1) + deltax
        end do

        xvect(n) = highrange
        mtx = 0

        do i = 1, n
        mtx(i,i) = 2/deltax**2 + 4
        end do

        do i = 1, n-1
        mtx(i,i+1) = -1/deltax**2
        end do

        mtx(n,n) = -1/deltax**2

        do i = 2, n
        mtx(i,i-1)= -1/deltax**2
        end do

        do i=1, n
            bside(i) = f(xvect(i))
        end do

        ! print *, bside
        !call printmatrix(mtx)
        call gem_FDM(mtx, coeff, bside, xvect)

    end subroutine FDM

    subroutine gen_coeffmat(f,lowrange,highrange,n,mtx,bside,xvec,coeff)
        use types
        use linalg
        interface
            function f(x) result(y)
                use types
                real(kind=rkind), intent(in):: x
                real(kind=rkind):: y
            end function f
        end interface

        real(kind=rkind), intent(in) :: lowrange,highrange
        real(kind=rkind), dimension(:,:), intent(out), allocatable :: mtx
        real(kind=rkind), dimension(:),intent(out), allocatable :: bside,xvec,coeff
        integer, intent(in) :: n
        integer :: i,j
        real(kind=rkind):: h

        allocate(xvec(n-2))

        allocate(bside(n))

        allocate(mtx(n,n))

        h= (highrange - lowrange)/2

        do i=1, ubound(xvec,1)
            xvec(i) = (highrange - lowrange)/(n+1)*i +lowrange
        end do

        print*,"Nodes: ", xvec

        ! do i =1, n
        !     mtx(1,i) = 2/h**(2) +4
        ! end do

        mtx(1:n,i) = 2/h**(2) +4



        ! do i =2, n
        !     mtx(1,i) = -1/h**(2)
        ! end do
        mtx(1:n,i+1) = -1/h**(2)

        ! do i =2, n
        !     mtx(2,i-1) = -1/h**(2)
        ! end do

        mtx(2,i-1) = -1/h**(2)

        do j =3, n
            do i=2,n
                mtx(j,i) = xvec(j-2)**(i-1)
            end do
        end do

        bside(1) = f(lowrange)
        bside(2) = f(highrange)

        do i=2, n-2
            bside(i+2) = f(xvec(i))
        end do
        call printmatrix(mtx)

        call gem(mtx, coeff, bside, xvec)


    end subroutine gen_coeffmat



end module simpleinteg
