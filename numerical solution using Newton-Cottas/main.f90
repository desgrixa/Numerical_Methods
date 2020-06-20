program main
  use types
  use integ
  use functions
  use linalg
  use simpleinteg

    ! real(kind= rkind) ::  output_nc, output_trap, integval
    ! integer :: i

    call newton_cottas(benchmark, 0_rkind, 1_rkind, 4, integval)
  ! output_nc = newton_cottas(benchmarknc, 0.5_rkind, 1.5_rkind, 10)
  ! print *, "function with known exact solution open interval: ", output_nc
  ! do i = 2, 25
  !     output_nc = newton_cottas(benchmarknc, 0.5_rkind, 1.5_rkind, i)
  !     print *, "function with known exact solution open interval: ", output_nc
  ! end do
  ! output_trap = trapinteg(0.5_rkind, 1.5_rkind, benchmarknc, 1e-8_rkind)
  ! print *, "function with known exact solution: ", output_trap



end program main
