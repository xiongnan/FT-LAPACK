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
               temp = rand()
            else
               temp = 0.0;
            end if
            lower(i, j) = temp
            upper(j, i) = temp   
         end do
      end do

      do i = 1, N
         do j = 1, N
            write(*, advance='no') lower(i, j), "  "
         end do
         Print *, "\n"
      end do

      Print *, "Hello World!"
      end program Hello
