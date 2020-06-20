program main
  use types
  use integ
  use functions

  real(kind = rkind) :: lowrange, highrange, dx, output
  real(kind = rkind) :: lowrange_1, highrange_1, dx_1, output_1
  real(kind = rkind) :: lowrange_2, highrange_2, dx_2, output_2
  real(kind = rkind) :: ThetaR, ThetaS, alpha, n
  real(kind = rkind) :: Kh
  real(kind = rkind) :: h1, h2, h3, h, Theta, Ks
  integer :: fileid, i

      ! 1. Use function with known exact solution to validate your numerical method.
      lowrange = 0
      highrange = 1
      dx = 1e-8

      output = trapinteg(lowrange, highrange, benchmark, dx)
      print *, "function with known exact solution: ", output

      !I have Average Silt Loam and Clay where:
      !log scale y
      Ks = 9.03e-7
      ThetaR = 0.0675
      ThetaS = 0.415
      alpha = 1.19
      n = 1.25
      !l = 0.5!!
      h1 = -0.1_rkind
      h2 = -12.0_rkind
      h3 = -1995.0_rkind
      Theta = retention_curve(h1,ThetaR,ThetaS, alpha, n) ! retention curve
      h = inverse(Theta,ThetaR,ThetaS, alpha, n)
      print *, "Value of h: ", h2
      print *, "Value of retention_curve: ", Theta
      print *, "Value of inverse retention_curve: ", h
      ! 10e-5
    ! First integral, where lower boundary is ThetaR, upper - Theta(h)
      lowrange_1 = ThetaR
      highrange_1 = Theta
      dx_1 = 1e-7

    ! Second integral, where lb - ThetaR, ub - ThetaS
      lowrange_2 = ThetaR + 1e-7
      highrange_2 = ThetaS - 1e-7
      dx_2 = 1e-7

      ! inv = inverse_1(Theta)
      ! print *, "1/h(Theta):", inv

      output_1 = trapinteg(lowrange_1, highrange_1, inverse_1, dx_1)
      output_2 = trapinteg(lowrange_2, highrange_2, inverse_1, dx_2)

      Kh = Ks * sqrt((Theta - ThetaR)/(ThetaS - ThetaR))*(output_1/output_2)**2

      print *, "Value of the first integral:", output_1
      print *, "Value of the second integral: ", output_2
      print *, "The result is: ", Kh

      !output data into a file
      open(fileid, file = "data1.dat", status="old", position="append", action="write", form='formatted')
      !-1000
      do i = 1, 1000
          Theta = retention_curve(h1,ThetaR,ThetaS, alpha, n)
          output_1 = trapinteg(lowrange_1, highrange_1, inverse_1, dx_1)
          output_2 = trapinteg(lowrange_2, highrange_2, inverse_1, dx_2)
          Kh = Ks * sqrt((Theta - ThetaR)/(ThetaS - ThetaR))*(output_1/output_2)**2
          write(fileid, *) h1, kh
          h1 = h1 - 0.01
      end do

end program main
