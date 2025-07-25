program example_is_diagonal
  use stdlib_linalg, only: is_diagonal
  implicit none
  real :: A(2, 2), B(2, 2)
  logical :: res
  A = reshape([1., 0., 0., 4.], shape(A))
  B = reshape([1., 0., 3., 4.], shape(B))
  res = is_diagonal(A) ! returns .true.
  res = is_diagonal(B) ! returns .false.
end program example_is_diagonal
