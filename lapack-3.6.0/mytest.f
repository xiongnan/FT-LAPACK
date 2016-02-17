      program hello
      
      real, dimension (:,:), allocatable :: upper
      real, dimension (:,:), allocatable :: lower
      integer :: N
      integer :: i, j
      real :: temp
      
      external DPOTRF
      
      N = 16
      
      allocate ( upper(N, N) )
      allocate ( lower(N, N) )
      
      do i = 1, N
         do j = 1, N
            if (j <= i) then
               temp = int(rand() * 10)
            else
               temp = 0.0;
            end if
            lower(i, j) = temp
            upper(j, i) = temp   
         end do
      end do
      
      do i=1, N
         Print *, ( lower(i,j), j=1,N )
      end do

      do i = 1, N
         do j = 1, N
!            write(*, advance='no') lower(i, j), "  "
          end do
!         Print *, "\n"
      end do
 100   format (1x, 16(1x,f5.1)
      Print *, "Hello World!"
      end program Hello
