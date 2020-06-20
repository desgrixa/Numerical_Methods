! here we specifying real and integer numbers
module types

  integer, parameter, public :: rkind = selected_real_kind(15,99)

  integer, parameter, public :: ikind = selected_int_kind(8)
end module types
