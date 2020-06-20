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


        function Gauss_quad(f, lowrange, highrange) result(ss)
            use types
            interface
                function f(x) result(y)
                    use types
                    real(kind = rkind), intent(in) :: x
                    real(kind = rkind) :: y
                end function f
            end interface

            real(kind = rkind) :: lowrange, highrange, dx,xm,xr
            real(kind = rkind) :: w(6),n(6)
            real(kind = rkind) :: ss
            integer :: i

            !The Gauss nodes and weights for a 6th-degree integration for ξ = -1 to ξ = +1 are

            n(1) = 0.93246951420315261
            n(2) = 0.66120938646626592
            n(3) = 0.23861918608319749
            n(4) = -n(3)
            n(5) = -n(2)
            n(6) = -n(1)

            w(1) = 0.17132449237916891
            w(2) = 0.36076157304813894
            w(3) = 0.46791393457269215
            w(4) = w(3)
            w(5) = w(2)
            w(6) = w(1)

            xm = (highrange+lowrange)/2
            xr = (highrange-lowrange)/2
            ss = 0

            do i = 1, 6
                dx = xr * n(i) ! point on actual range from a to b
                !ss = ss + w(i) * (f(xm+dx)+f(f(xm-dx)))
                ss = ss + w(i) * (f(xm - dx))

            end do

            ss = ss * xr



        end function Gauss_quad


end module integ
