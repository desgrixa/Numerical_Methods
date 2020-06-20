module functions
  contains
      function retention_curve(h, ThetaR, ThetaS, alpha, n) result(Thetah)
          use types
          real(kind = rkind), intent(in) :: h, ThetaR, ThetaS, alpha, n
          real(kind = rkind) :: Thetah

              Thetah = ThetaR + ((ThetaS-ThetaR)/((1+(-alpha*h)**n)**(1-1/n)))
             ! print *, "Theta: ", Thetah

      end function retention_curve

      function inverse(Thetah, ThetaR, ThetaS, alpha, n) result(h)
          use types
          real(kind = rkind), intent(in) :: Thetah, ThetaR, ThetaS, alpha, n
          real(kind = rkind) :: h

              h = -(1/alpha)* ((((ThetaS - ThetaR)/(Thetah - ThetaR))**(1/(1-(1/n)))-1)**(1/n))
              !h = 1/h
              !print *, "inverse:", h

      end function inverse

      function inverse_1(Thetah) result(h)
          use types
          real(kind = rkind), intent(in) :: Thetah
          real(kind = rkind) :: ThetaR, ThetaS, alpha, n
          real(kind = rkind) :: h

              ThetaR = 0.0675
              ThetaS = 0.415
              alpha = 1.19
              n = 1.25
              
              h = -(1/alpha)* ((((ThetaS - ThetaR)/(Thetah - ThetaR))**(1/(1-(1/n)))-1)**(1/n))
              h = 1/h
              !print *, "inverse:", h

      end function inverse_1
      ! function mualem() result(Kh)
      !     use types
      !     use integ
      !
      !     real(kind = rkind) :: ThetaR, ThetaS, Ks, Thetah
      !     real(kind = rkind) :: lowrange_1, highrange_2, dx_1, output_1
      !     real(kind = rkind) :: lowrange_2, highrange_2, dx_2, output_2
      !     real(kind = rkind) :: Kh
      !
      !     output_1 = trapinteg(lowrange_1, highrange_1, inverse, dx_1)
      !     output_2 = trapinteg(lowrange_2, highrange_2, inverse, dx_2)
      !     Kh = Ks * sqer((Thetah - ThetaR)/(ThetaS - ThetaR))*(output_1/output_2)**2
      !
      ! end function mualem

      function benchmark(x) result(y)
          use types
          real(kind = rkind), intent(in) :: x
          real(kind = rkind) :: y

            y = x**3

      end function benchmark


end module functions
