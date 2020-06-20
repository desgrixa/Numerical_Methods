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
            print *, "No traps:", notraps
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

        function Gauss_quad_10(f, lowrange, highrange) result(ss)
            use types
            interface
                function f(x) result(y)
                    use types
                    real(kind = rkind), intent(in) :: x
                    real(kind = rkind) :: y
                end function f
            end interface

            real(kind = rkind) :: lowrange, highrange, dx,xm,xr
            real(kind = rkind) :: w(5),n(5)
            real(kind = rkind) :: ss
            integer :: i

            !The Gauss nodes and weights for a 6th-degree integration for ξ = -1 to ξ = +1 are

            n(1) = 0.97390652851715043
            n(2) = 0.86506336668888839
            n(3) = 0.67940956829889076
            n(4) = 0.43339539412919476
            n(5) = 0.14887433898164026


            w(1) = 6.6671344308743441E-002
            w(2) = 0.14945134915065394
            w(3) = 0.21908636251596847
            w(4) = 0.26926671930987656
            w(5) = 0.29552422471475753

            xm = (highrange+lowrange)/2
            xr = (highrange-lowrange)/2
            ss = 0

            do i = 1, 5
                dx = xr * n(i) ! point on actual range from a to b
                ss = ss + w(i) * (f(xm+dx)+f(xm-dx))
                !ss = ss + w(i) * (f(xm - dx))

            end do

            ss = ss * xr

        end function Gauss_quad_10

        function Gauss_quad_8(f, lowrange, highrange) result(ss)
            use types
            interface
                function f(x) result(y)
                    use types
                    real(kind = rkind), intent(in) :: x
                    real(kind = rkind) :: y
                end function f
            end interface

            real(kind = rkind) :: lowrange, highrange, dx,xm,xr
            real(kind = rkind) :: w(8),n(8)
            real(kind = rkind) :: ss
            integer :: i

            !The Gauss nodes and weights for a 6th-degree integration for ξ = -1 to ξ = +1 are

            n(1) = 0.96028985649754117
            n(2) = 0.79666647741365249
            n(3) = 0.52553240991637584
            n(4) = 0.18343464249567387
            n(5) = -n(4)
            n(6) = -n(3)
            n(7) = -n(2)
            n(8) = -n(1)





            w(1) = 0.10122853629036332
            w(2) = 0.22238103445334745
            w(3) = 0.31370664587788299
            w(4) = 0.36268378337840629
            w(5) = w(4)
            w(6) = w(3)
            w(7) = w(2)
            w(8) = w(1)


            xm = (highrange+lowrange)/2
            xr = (highrange-lowrange)/2
            ss = 0

            do i = 1, 8
                dx = xr * n(i) ! point on actual range from a to b
                !ss = ss + w(i) * (f(xm+dx)+f(f(xm-dx)))
                ss = ss + w(i) * (f(xm - dx))

            end do

            ss = ss * xr


        end function Gauss_quad_8

        function Gauss_quad_4(f, lowrange, highrange) result(ss)
            use types
            interface
                function f(x) result(y)
                    use types
                    real(kind = rkind), intent(in) :: x
                    real(kind = rkind) :: y
                end function f
            end interface

            real(kind = rkind) :: lowrange, highrange, dx,xm,xr
            real(kind = rkind) :: w(2),n(2)
            real(kind = rkind) :: ss
            integer :: i

            !The Gauss nodes and weights for a 6th-degree integration for ξ = -1 to ξ = +1 are

            n(1) = 0.86113631159405246
            n(2) = 0.33998104358485604

            w(1) = 0.34785484513745413
            w(2) = 0.65214515486254587

            xm = (highrange+lowrange)/2
            xr = (highrange-lowrange)/2
            ss = 0

            do i = 1, 2
                dx = xr * n(i) ! point on actual range from a to b
                ss = ss + w(i) * (f(xm+dx)+f(xm-dx))
                !ss = ss + w(i) * (f(xm - dx))

            end do

            ss = ss * xr


        end function Gauss_quad_4

        function Gauss_quad_5(f, lowrange, highrange) result(ss)
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

            n(1) = 0.90617984593866419
            n(2) = 0.53846931010568377
            n(3) = 0.0000000000000000
            n(4) = -n(3)
            n(5) = -n(2)
            n(6) = -n(1)

            w(1) = 0.23692688505618875
            w(2) = 0.47862867049936592
            w(3) = 0.56888888888889078
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

        end function Gauss_quad_5


        function Gauss_quad_14(f, lowrange, highrange) result(ss)
            use types
            interface
                function f(x) result(y)
                    use types
                    real(kind = rkind), intent(in) :: x
                    real(kind = rkind) :: y
                end function f
            end interface

            real(kind = rkind) :: lowrange, highrange, dx,xm,xr
            real(kind = rkind) :: w(7),n(7)
            real(kind = rkind) :: ss
            integer :: i

            !The Gauss nodes and weights for a 6th-degree integration for ξ = -1 to ξ = +1 are

            n(1) = 0.98628380869598375
            n(2) = 0.92843488365859739
            n(3) = 0.82720131505416383
            n(4) = 0.68729290477313043
            n(5) = 0.51524863628658224
            n(6) = 0.31911236884487243
            n(7) = 0.10805494866898753


            w(1) = 3.5119460333944197E-002
            w(2) = 8.0158087166334990E-002
            w(3) = 0.12151857070369276
            w(4) = 0.15720316718850080
            w(5) = 0.18553839750797962
            w(6) = 0.20519846370689887
            w(7) = 0.21526385339264870

            xm = (highrange+lowrange)/2
            xr = (highrange-lowrange)/2
            ss = 0

            do i = 1, 7
                dx = xr * n(i) ! point on actual range from a to b
                ss = ss + w(i) * (f(xm+dx)+f(xm-dx))
                !ss = ss + w(i) * (f(xm - dx))

            end do

            ss = ss * xr


        end function Gauss_quad_14

        function Gauss_quad_16(f, lowrange, highrange) result(ss)
            use types
            interface
                function f(x) result(y)
                    use types
                    real(kind = rkind), intent(in) :: x
                    real(kind = rkind) :: y
                end function f
            end interface

            real(kind = rkind) :: lowrange, highrange, dx,xm,xr
            real(kind = rkind) :: w(8),n(8)
            real(kind = rkind) :: ss
            integer :: i

            !The Gauss nodes and weights for a 6th-degree integration for ξ = -1 to ξ = +1 are

            n(1)	=	0.98940093490327152
            n(2)	=	0.94457502263193316
            n(3)	=	0.86563120140699679
            n(4)	=	0.75540440678842946
            n(5)	=	0.61787624235454719
            n(6)	=	0.45801677539872154
            n(7)	=	0.28160354887422495
            n(8)	=	9.5012509068704981E-002


            w(1)	=	2.7152459635771731E-002
            w(2)	=	6.2253524405201821E-002
            w(3)	=	9.5158512270877793E-002
            w(4)	=	0.12462897181308819
            w(5)	=	0.14959598919746009
            w(6)	=	0.16915651938365350
            w(7)	=	0.18260341430016042
            w(8)	=	0.18945060899378635

            xm = (highrange+lowrange)/2
            xr = (highrange-lowrange)/2
            ss = 0

            do i = 1, 8
                dx = xr * n(i) ! point on actual range from a to b
                ss = ss + w(i) * (f(xm+dx)+f(xm-dx))
                !ss = ss + w(i) * (f(xm - dx))

            end do

            ss = ss * xr


        end function Gauss_quad_16

        function Gauss_quad_20(f, lowrange, highrange) result(ss)
            use types
            interface
                function f(x) result(y)
                    use types
                    real(kind = rkind), intent(in) :: x
                    real(kind = rkind) :: y

                end function f
            end interface

            real(kind = rkind) :: lowrange, highrange, dx,xm,xr
            real(kind = rkind) :: w(10),n(10)
            real(kind = rkind) :: ss
            integer :: i

            !The Gauss nodes and weights for a 6th-degree integration for ξ = -1 to ξ = +1 are

            n(1)	=	0.99312841421418818
            n(2)	=	0.96397097178340663
            n(3)	=	0.91223216858868539
            n(4)	=	0.83911303088165534
            n(5)	=	0.74632616583371703
            n(6)	=	0.63604640451484151
            n(7)	=	0.51085895347401455
            n(8)	=	0.37369854969503713
            n(9)	=	0.22778040108574302
            n(10)	=	7.6524522568086201E-002


            w(1)	=	1.7614479691862694E-002
            w(2)	=	4.0602485937562438E-002
            w(3)	=	6.2673574830429679E-002
            w(4)	=	8.3278534933615589E-002
            w(5)	=	0.10193186212156746
            w(6)	=	0.11819577618974186
            w(7)	=	0.13168884817770543
            w(8)	=	0.14209482541707841
            w(9)	=	0.14917012651161810
            w(10)	=	0.15274948618881828

            xm = (highrange+lowrange)/2
            xr = (highrange-lowrange)/2
            ss = 0

            do i = 1, 10
                dx = xr * n(i) ! point on actual range from a to b
                ss = ss + w(i) * (f(xm+dx)+f(xm-dx))
                !ss = ss + w(i) * (f(xm - dx))

            end do

            ss = ss * xr


        end function Gauss_quad_20
end module integ
