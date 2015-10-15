! Copyright (c) 2015 Florian Wende (flwende@gmail.com)
!
! Distributed under the BSD 2-clause Software License 
! (See accompanying file LICENSE)

#if defined(__MIC__)
#define VL 8
#elif defined(__AVX__) || defined(__AVX2__)
#define VL 4
#elif defined(__SSE2__) || defined(__SSE3__) || defined(__SSE4_1__) || defined(__SSE4_2__)
#define VL 2
#else
#define VL 1
#endif

#define DEPTH 20

#if defined MANUALVEC
module simd
      use, intrinsic :: iso_c_binding
      type, public :: simd_real8
      real(c_double) :: x(0:VL-1)
      end type simd_real8
      type, public :: simd_int4
      real(c_int) :: x(0:VL-1)
      end type simd_int4
      type, public :: simd_mask8
      logical :: x(0:VL-1)
      end type simd_mask8
end module simd
#endif

program main
      use, intrinsic :: iso_c_binding
#if defined MANUALVEC
      use simd
#endif
      implicit none
      ! srand48
      interface 
      subroutine srand48(s) bind(C, name="srand48")
      use, intrinsic :: iso_c_binding
      integer(c_long), value :: s
      end subroutine srand48
      end interface
      ! drand48
      interface 
      function drand48() bind(C, name="drand48")
      use, intrinsic :: iso_c_binding
      real(c_double) :: drand48
      end function drand48
      end interface
      ! openmp timer
      interface 
      function omp_get_wtime() bind(C, name="omp_get_wtime")
      use, intrinsic :: iso_c_binding
      real(c_double) :: omp_get_wtime
      end function omp_get_wtime
      end interface
      ! variables
      integer :: i, idx
      integer, parameter :: n = 8 * 1024 * 1024
      real*8, allocatable :: x1(:), x2(:), y(:), yref(:)
      real*8 :: dev
      real*8 :: t0, t1
#if defined MANUALVEC
      integer :: ii
      type(simd_real8) :: buffer_x1, buffer_x2, buffer_y
      type(simd_mask8) :: m
#endif
      allocate(x1(n), x2(n), y(n), yref(n))
      ! init random number generator
      call srand48(1_8)
      ! initialize x1 x2 array
      do i=1, n
         x1(i) = 2.0 * drand48()
         x2(i) = drand48()
      enddo
      ! call function foo
      t0 = omp_get_wtime()
#if defined MANUALVEC
      do i=1, n, VL
!$omp simd
         do ii=0, VL-1
            m%x(ii) = .false.
            if ((i + ii) .le.  n) then 
               m%x(ii) = .true.
               buffer_x1%x(ii) = x1(i + ii)
               buffer_x2%x(ii) = x2(i + ii)
            endif
         enddo
!dir$ noinline
         call foo_simd(buffer_x1, buffer_x2, DEPTH, buffer_y, m)
!$omp simd
         do ii=0, VL-1
            if (m%x(ii)) y(i + ii) = buffer_y%x(ii)
         enddo
      enddo
#else
#if defined AUTOVEC
!$omp simd
#endif
      do i=1, n
!dir$ noinline
         call foo(x1(i), x2(i), DEPTH, y(i))
      enddo
#endif
      t1 = omp_get_wtime()
      ! check results
      call srand48(1_8)
      do i=1, n
         x1(i) = 2.0 * drand48()
         x2(i) = drand48()
      enddo
!dir$ novector
      do i=1, n
!dir$ noinline
         call foo(x1(i), x2(i), DEPTH, yref(i))
      enddo
      dev = 0.0
      do i=1, n
         dev = dev + (y(i) - yref(i))**2
      enddo
      dev = sqrt(dev)
      ! print output
      idx = 1 + int((n - 16) * drand48())
      print *, y(idx:idx + 16)
      print *, "elapsed time=", (t1 - t0) * 1.0E3, "ms"
      print *, "deviation=", dev
      deallocate(x1, x2, y, yref)
end program main

subroutine foo(x1, x2, d, y)
#if defined(AUTOVEC)
!$omp declare simd(foo) simdlen(VL)
#endif
      use, intrinsic :: iso_c_binding
      real*8, intent(in) :: x1, x2
      integer, intent(in) :: d
      real*8, intent(out) :: y
      integer :: i, i_max
      i = 1
      i_max = int(d * x2)
      y = 0.0
      do while (i .le. i_max)
         y = sqrt(x1 + y)
         if (y .gt. 1.0) y = log(y)
         i = i + 1
      end do
end subroutine foo

#if defined MANUALVEC
subroutine foo_simd(x1, x2, d, y, m0)
      use, intrinsic :: iso_c_binding
      use simd
      type(simd_real8), intent(in) :: x1, x2
      integer, intent(in) :: d
      type(simd_real8), intent(out) :: y
      type(simd_mask8), intent(in) :: m0
      type(simd_real8) :: a
      type(simd_mask8) :: m1
      type(simd_int4) :: i, i_max
      real*8 :: temp
      integer(c_int) :: itemp
      integer :: ii
      logical :: true_for_any
      true_for_any = .false.
#if defined __INTEL_COMPILER
!$omp simd reduction(.or. : true_for_any) aligned(x2,y)
#else
!$omp simd reduction(.or. : true_for_any)
#endif
      do ii=0, VL-1
         i%x(ii) = 1
         i_max%x(ii) = int(d * x2%x(ii))
         if (m0%x(ii) .and. (i%x(ii) .le. i_max%x(ii))) then
            m1%x(ii) = .true.
            true_for_any = .true.
         else
            m1%x(ii) = .false.
         endif
         y%x(ii) = 0.0
      end do
      do while (true_for_any)
         true_for_any = .false.
#if defined __INTEL_COMPILER
!$omp simd reduction(.or. : true_for_any) aligned(x1,x2,y)
#else
!$omp simd reduction(.or. : true_for_any)
#endif
         do ii=0, VL-1
            temp = sqrt(x1%x(ii) + y%x(ii))
            if (temp .gt. 1.0) temp = log(temp)
            i%x(ii) = i%x(ii) + 1
            if (m1%x(ii)) then
               y%x(ii) = temp
               if (i%x(ii) .le. i_max%x(ii)) then
                  true_for_any = .true.
               else
                  m1%x(ii) = .false.
               endif
            endif
         end do
      end do
end subroutine foo_simd
#endif
