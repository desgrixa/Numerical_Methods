program main
  use types
  use integ
  use functions
  use linalg
  use simpleinteg

    real(kind= rkind) ::  output_nc, output_trap, div_check
    real(kind = rkind) :: lowrange_1, highrange_1, output_1, dx_1, lowrange_3
    real(kind = rkind) :: lowrange_2, highrange_2, output_2, dx_2, highrange_3
    real(kind = rkind) :: output_3, output_4, output_5, output_6
    real(kind = rkind) :: ThetaR, ThetaS, alpha, n
    real(kind = rkind) :: Kh, Kh_nc, Kh_qauss
    real(kind = rkind) :: h1, h2, h3, h, Theta, Ks
    ! integer :: i

  !call newton_cottas_subr(benchmarknc, 0.5_rkind, 1.5_rkind, 4, integval)

  output_nc = Gauss_quad_8(benchmarknc, 0.5_rkind, 1.5_rkind)
  !print *, "function with known exact solution open interval: ", output_nc
  print *, "Gauss_quad: ", output_nc



  output_trap = trapinteg(0.5_rkind, 1.5_rkind, benchmarknc, 1e-8_rkind)
  print *, "function with known exact solution(trapinteg): ", output_trap


    div_check = newton_cottas_func(benchmarknc, 0.5_rkind, 1.5_rkind, 10)
    print *, "p=10", "function with known exact solution open interval: ", div_check



! For any iterative numerical technique, each successive iteration results in a
! solution that moves progressively closer to the true solution. This is
! known as convergence. A numerical method is not always guaranteed to produce
! converging results. Convergence is subject to satisfying certain conditions.
! If these conditions are not met, each successive iteration produces a result
! that progressively moves away from the true solution. This is known as divergence.

!I have Average Silt Loam and Clay where:
   !log scale y
    Ks = 9.03e-7
    ThetaR = 0.0675
    ThetaS = 0.415
    alpha = 1.19
    n = 1.25
    !l = 0.5!!
    h1 = -0.5_rkind
    h2 = -12.0_rkind
    h3 = -1995.0_rkind
    Theta = retention_curve(h1,ThetaR,ThetaS, alpha, n) ! retention curve
    h = inverse(Theta,ThetaR,ThetaS, alpha, n)

    lowrange_1 = ThetaR
    highrange_1 = Theta
    dx_1 = 1e-7

    ! Second integral, where lb - ThetaR, ub - ThetaS
    lowrange_2 = ThetaR + 1e-7
    highrange_2 = ThetaS - 1e-7
    dx_2 = 1e-7

    lowrange_3 = ThetaR
    highrange_3 = ThetaS


    output_1 = trapinteg(lowrange_1, highrange_1, inverse_1, dx_1)
    output_2 = trapinteg(lowrange_2, highrange_2, inverse_1, dx_2)
    Kh = Ks * sqrt((Theta - ThetaR)/(ThetaS - ThetaR))*(output_1/output_2)**2


    output_3 = newton_cottas_func(inverse_1, lowrange_1, highrange_1, 46)
    output_4 = newton_cottas_func(inverse_1, lowrange_2, highrange_2, 46)
    Kh_nc = Ks * sqrt((Theta - ThetaR)/(ThetaS - ThetaR))*(output_3/output_4)**2
    ! print *, "j", j
    ! print *, "Value of the first integral (newton_cottas):", output_3
    ! print *, "Value of the second integral (newton_cottas): ", output_4
    ! print *, "The result is (newton_cottas): ", Kh_nc

    output_5 = Gauss_quad_20(inverse_1, lowrange_1, highrange_1)
    output_6 = Gauss_quad_20(inverse_1, lowrange_2, highrange_2)
    Kh_qauss = Ks * sqrt((Theta - ThetaR)/(ThetaS - ThetaR))*(output_5/output_6)**2

    ! print *, "Value of the first integral (trap):", output_1
    ! print *, "Value of the second integral (trap): ", output_2
    ! print *, "Value of the first integral (newton_cottas):", output_3
    ! print *, "Value of the second integral (newton_cottas): ", output_4
    print *, "The result is (trap): ", Kh
    print *, "The result is (newton_cottas): ", Kh_nc
    print *, "The result is (Gauss_quad): ", Kh_qauss


end program main
