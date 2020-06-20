program main
  use types
  use integ
  use functions
  use linalg
  use simpleinteg

    real(kind= rkind) ::  output_c, output_f
    real(kind = rkind) ::  x1, exact
    integer :: fileid, i

  output_c = collocation(ODE2_sin, 0.0_rkind, 1.0_rkind, 5)
  print *, "collocation: ", output_c


  output_f =  finite_method(ODE2_sin, 0.0_rkind, 1.0_rkind, 0.02_rkind)
  !print *, "FDM: ", output_f

  open(newunit = fileid, file = 'data1.dat', status = 'old')
  x1 = 0.000001

  do i =1 , 100
     exact = ODE2(x1)
     write(fileid, *) x1, exact
     x1 =  x1 + 0.01
 end do


end program main
