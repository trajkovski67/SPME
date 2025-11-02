program test
  use, intrinsic :: iso_c_binding
  implicit none

  interface
    function fftw_plan_dft_1d(n, in, out, sign, flags) bind(C, name="fftw_plan_dft_1d")
      use iso_c_binding
      integer(c_int), value :: n, sign, flags
      type(c_ptr), value :: in, out
      type(c_ptr) :: fftw_plan_dft_1d
    end function

    subroutine fftw_execute_dft(plan, in, out) bind(C, name="fftw_execute_dft")
      use iso_c_binding
      type(c_ptr), value :: plan, in, out
    end subroutine

    subroutine fftw_destroy_plan(plan) bind(C, name="fftw_destroy_plan")
      use iso_c_binding
      type(c_ptr), value :: plan
    end subroutine
  end interface

  integer, parameter :: n = 8
  complex(c_double_complex), target :: in(n), out(n)
  type(c_ptr) :: plan

  integer :: i

  ! Initialize input
  do i = 1, n
    in(i) = cmplx(dble(i), 0.0d0, kind=c_double)
  end do

  ! Create FFT plan
  plan = fftw_plan_dft_1d(n, c_loc(in), c_loc(out), -1, 64)

  ! Execute FFT
  call fftw_execute_dft(plan, c_loc(in), c_loc(out))

  ! Print output
  print *, 'FFT output:'
  do i = 1, n
    print '(I2, 2F12.6)', i, real(out(i)), aimag(out(i))
  end do

  ! Destroy plan
  call fftw_destroy_plan(plan)

end program test

