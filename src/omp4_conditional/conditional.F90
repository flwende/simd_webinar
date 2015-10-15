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

module simd
      use, intrinsic :: iso_c_binding
      type, public :: simd_real8
      real(c_double) :: x(0:VL-1)
      end type simd_real8
      type, public :: simd_mask8
      logical :: x(0:VL-1)
      end type simd_mask8
end module simd

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
      real(c_double), allocatable :: x1(:), x2(:), y(:), yref(:)
      real(c_double) :: dev
      real(c_double) :: t0, t1
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
         call foo_simd(buffer_x1, buffer_x2, buffer_y, m)
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
         call foo(x1(i), x2(i), y(i))
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
         call foo(x1(i), x2(i), yref(i))
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

subroutine foo(x1, x2, y)
#if defined AUTOVEC
!$omp declare simd(foo) simdlen(VL)
#endif
      use, intrinsic :: iso_c_binding
      real(c_double), intent(in) :: x1, x2
      real(c_double), intent(out) :: y
      if (x2 .gt. 0.5) then
         y = sqrt(x1)
         if (y .gt. 1.0) y = log(y)
      else  
         y = 0.0
      endif
end subroutine foo

subroutine foo_simd(x1, x2, y, m0)
      use simd
      type(simd_real8), intent(in) :: x1, x2
      type(simd_real8), intent(out) :: y
      type(simd_mask8), intent(in) :: m0
      real(c_double) :: temp
#if defined __INTEL_COMPILER
!$omp simd aligned(x1,x2,y)
#else
!$omp simd
#endif
      do ii=0, VL-1
         if (x2%x(ii) .gt. 0.5) then
            temp = sqrt(x1%x(ii))
            if (temp .gt. 1.0) temp = log(temp)
         else  
            temp = 0.0
         endif
         if (m0%x(ii)) y%x(ii) = temp
      enddo
end subroutine foo_simd
