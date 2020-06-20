module simpleinteg
    ! make different calls for different values of n.
    !public :: newton_cottas ! for closed interval
    private :: gen_coefmatrix_open, gen_coefmatrix_closed
    private :: integ_polynomials
    contains
        function newton_cottas_func(f,lowrange,highrange, n) result(integval)
            use types
            interface
                function f(x) result(y)
                    use types
                    real(kind = rkind), intent(in) :: x
                    real(kind = rkind) :: y
                end function f
            end interface

            real(kind = rkind), intent(in) :: lowrange, highrange
            real(kind = rkind) :: integval
            ! number of integration nodes
            integer :: n
            real(kind = rkind), dimension(:,:), allocatable :: mtx
            real(kind = rkind), dimension(:), allocatable :: bside, xvect, coeff

            !most important step to have coeff
            call gen_coefmatrix_open(f, lowrange, highrange, n, mtx, bside, xvect, coeff)
            !
            call integ_polynomials(coeff, lowrange, highrange, integval)

            !print *, "function with known exact solution: ", integval

        end function newton_cottas_func

        subroutine newton_cottas_subr(f,lowrange,highrange, n, integval)
            use types
            interface
                function f(x) result(y)
                    use types
                    real(kind = rkind), intent(in) :: x
                    real(kind = rkind) :: y
                end function f
            end interface

            real(kind = rkind), intent(in) :: lowrange, highrange, integval
            !real(kind = rkind) :: integval
            ! number of integration nodes
            integer :: n
            real(kind = rkind), dimension(:,:), allocatable :: mtx
            real(kind = rkind), dimension(:), allocatable :: bside, xvect, coeff

            !most important step to have coeff
            call gen_coefmatrix_closed(f, lowrange, highrange, n, mtx, bside, xvect, coeff)
            !
            call integ_polynomials(coeff, lowrange, highrange, integval)

            print *, "function with known exact solution closed interval: ", integval

        end subroutine newton_cottas_subr


        subroutine gen_coefmatrix_closed(f,lowrange,highrange,n,mtx,bside, xvect, coeff)
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

            real(kind = rkind), intent(in) :: lowrange, highrange
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
                xvect(i) = (highrange-lowrange)/(n-1) * i + lowrange
            end do

            print *, xvect
            ! now the matrix is allocated

            mtx(1:n,1) = 1.0_rkind ! (1:n, 1) -> means row number 1 till n, and first col
        !rkind - just to be accurate enough
            do i=2, n  ! from the second col till n
                mtx(1,i) = lowrange**(i-1) ! 1 row till i, equal to
            end do

            do i=2, n ! for the second row in the matrix
                mtx(2,i) = highrange**(i-1)
            end do

            do j=3, n ! rest of the matrix (row 3, 4, 5)
                do i=2, n
                    mtx(j,i) = xvect(j-2)**(i-1)
                end do
            end do
            ! procedure to bside
            bside(1) = f(lowrange)
            bside(2) = f(highrange)

            do i=1, n-2 ! beacuse first to are reserved
                bside(i+2) = f(xvect(i))
            end do

            call gem(mtx, coeff, bside) ! when we call gausiion elemination the
            !coeff will be a vector of coeff of our polynomial approximation
            !so then we will have our procedure ready.
            ! this gives us ability to integrate functions
        ! matrix generated
        end subroutine gen_coefmatrix_closed
        !
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

        subroutine gen_coefmatrix_open(f,lowrange,highrange,n,mtx,bside, xvect, coeff)
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

            real(kind = rkind), intent(in) :: lowrange, highrange
            real(kind = rkind), dimension(:,:), allocatable, intent(out) :: mtx
            real(kind = rkind), dimension(:), allocatable, intent(out) :: bside, xvect, coeff
            ! degree of polynomial: (n-1); n - (number of internal points)
            integer, intent(in) :: n
            integer :: i, j

            allocate(xvect(n)) !n-2 is a dimension (internal nodes)
             !xvect and bide have a dimension given by number of n
            allocate(bside(n)) !dimension refers to number of coeff
            allocate(mtx(n,n)) !matrix for the given dimension

            do i=1, ubound(xvect,1) ! loop until dimension of xvect (to the range of xvect (ubound x, n-2) )
                xvect(i) = (highrange-lowrange)/(n+1) * i + lowrange
            end do

            !print *, "Open interval values", xvect
            ! now the matrix is allocated

            mtx(1:n,1) = 1.0_rkind ! (1:n, 1) -> means row number 1 till n, and first col
        !rkind - just to be accurate enough
            ! do i=2, n-2  ! from the second col till n
            !     mtx(1,i) = lowrange**(i-1) ! 1 row till i, equal to
            ! end do
            !
            ! do i=2, n ! for the second row in the matrix
            !     mtx(2,i) = highrange**(i-1)
            ! end do

            do j=1, n! rows (2, 3, 4)
                do i=2, n
                    mtx(j,i) = xvect(j)**(i-1)
                end do
            end do
            ! procedure to bside
            ! bside(1) = f(xvect(1))
            ! bside(2) = f(xvect(n-2))

            do i=1, n
                bside(i) = f(xvect(i))
            end do

            call gem(mtx, coeff, bside) ! when we call gausiion elemination the
            !coeff will be a vector of coeff of our polynomial approximation
            !so then we will have our procedure ready.
            ! this gives us ability to integrate functions
        ! matrix generated
    end subroutine gen_coefmatrix_open


end module simpleinteg
