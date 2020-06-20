module functions
    public :: benchmark
    public :: benchmarknc

    contains

      function benchmark(x) result(y)
          use types
          real(kind = rkind), intent(in) :: x
          real(kind = rkind) :: y

            y = x**3

      end function benchmark

      function benchmarknc(x) result(y)
          use types
          real(kind = rkind), intent(in) :: x
          real(kind = rkind) :: y

            y = x*sin(x*x)

      end function benchmarknc


      function ODE2(x) result(ux)
      		use types
      		real(kind=rkind), intent(in):: x
      		real(kind=rkind):: alpha=0.95, beta=1.95, ux

      		ux= (alpha*exp(2*x)) + (beta*exp(-2*x)) + (1/13 * sin(3*x))

      end function ODE2

      function ODE2_sin(x) result(y)
              use types
              real(kind=rkind), intent(in):: x
              real(kind=rkind):: y

                y = sin(3*x)

      end function ODE2_sin

      function difficult(x,coeff) result (y)

          use types
          real(kind=rkind), intent(in):: x
          real(kind = rkind), dimension(:), intent(in) :: coeff
          real(kind=rkind):: y
          real(kind=rkind)::double_def_u,u
          integer:: n,i

          n = ubound(coeff,1)
          double_def_u = 0.0_rkind

          do i=2, n

              double_def_u= double_def_u + i*(i-1)*coeff(i)*x**(i-2)

          end do

          u = 1.0_rkind

          do i = 1,n

              u= u+coeff(i)*x**i
          end do

          y= -double_def_u+4*u

          !y = -(2*coeff(2)+6*coeff(3)*x+12*coeff(4)*x**2)+4+4*(coeff(1)*x+coeff(2)*x**2+coeff(3)*x**3+coeff(4)*x**4)

      end function


end module functions
