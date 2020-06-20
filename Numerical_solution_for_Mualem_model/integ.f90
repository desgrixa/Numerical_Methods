module integ
    contains
        function trapinteg(a, b, f, dx) result(integ)
            use types

            real(kind = rkind), intent(in) :: a, b, dx

            interface
                function f(x) result(y)
                    use types
                    real(kind = rkind), intent(in) :: x
                    real(kind = rkind) :: y

                end function f
            end interface

            real(kind = rkind) :: integ

            real(kind = rkind) :: length, trap_area, x1, x2, dxupdated
            integer :: notraps, i

            length = abs(b - a)
            notraps = int(length / dx) + 1 ! sequance to compute number of traps

            dxupdated = length / notraps

            !print *, length/dx, int(length/dx), huge(1) ; stop
            !print *, "domain lenght is:", length
            !print *, "we have updated original dx, which was:", dx, "to value", dxupdated

            x1 = a
            x2 = x1 + dx
            integ = 0 ! sum of our trapezoids is eqal to integral

            do i =1, notraps
                trap_area = (f(x1) + f(x2)) / 2 * dxupdated
                integ = integ + trap_area
                x1 = x2
                x2 = x1 + dx
            end do

        end function trapinteg

end module integ
