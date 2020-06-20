module linalg
    public :: genmatrix
    public :: gem
    contains

        subroutine genmatrix(A, matdim)
            use types

            real(kind = rkind), dimension(:, :), intent(out), allocatable :: A
            integer, intent(in) :: matdim

            integer :: i, j
            allocate(A(matdim, matdim))
            do i =1, ubound(A, 1)
                do j = 1, ubound(A,1)
                    A(i,j) = 1.0_rkind/(i +j)
                end do
            end do

        end subroutine genmatrix

        subroutine readmatrix(A)
            use types
            real(kind = rkind), dimension(:,:), intent(in out), allocatable :: A

            !real(kind = rkind), dimension(:), intent(out) :: x
            !real(kind = rkind), dimension(:), intent(out) :: b

            integer :: filematrix, Adim, i

            open(newunit = filematrix, file = "matrix.alex")

            read(filematrix, *) Adim

            !print *, "matrix dimension is", Adim

            allocate(A(Adim, Adim))

            do i = 1, ubound(A,2)
                read(filematrix, *) A(i, :)
            end do

            close(filematrix)

        end subroutine readmatrix

        subroutine printmatrix(A)
            use types
            real(kind = rkind), dimension(:,:), intent(in) :: A

            integer :: i

            do i = 1, ubound(A, 1)
                print *, "row no.", i, "is", A(i, :)
            end do

        end subroutine printmatrix

        subroutine gem(A, x, b) !gassian elemenation method
            use types
            real(kind = rkind), dimension(:,:), intent(in) :: A
            real(kind = rkind), dimension(:), intent(out), allocatable :: x
            real(kind = rkind), dimension(:), intent(in) :: b

            real(kind = rkind), dimension(:,:), allocatable :: aside
            !real(kind = rkind), dimension(:,:) :: aside
            integer :: i, j, n

            n = ubound(A, 1)

            allocate(aside(n, n+1))


            if (.not. allocated(x)) then
                allocate(x(n))
            end if

            aside(:,1:n) = A
            aside(:,n+1) = b


            do j = 2, n
                do i = j, n
                    aside(i,:) = aside(i,:) + aside(j-1,:)/aside(j-1, j-1)*(-aside(i,j-1))
                    !b(i) = b(i) + b(j-1) / A(j-1, j-1) * (-A(i,j-1))
                end do
            end do
            !print *, "the updated matrix is:"

            !call printmatrix(aside)

            x(n) = aside(n,n+1)/aside(n, n)

            do i = n-1, 1, -1
                x(i) = (aside(i,n+1) - dot_product(x(i+1:n), aside(i, i+1:n))) / aside(i,i)
            end do

            !print *, x

        end subroutine gem

end module linalg
! open(newunit = fileid, file = 'data2.dat', status = 'old')
! write(fileid, *) x(n), b(n)
! do i = n-1, 1, -1
!     x(i) = (aside(i,n+1) - dot_product(x(i+1:n), aside(i, i+1:n))) / aside(i,i)
!     write(fileid, *) x(i), b(i)
! end do
