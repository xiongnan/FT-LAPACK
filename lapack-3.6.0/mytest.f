      program hello
      
      real, dimension (:,:), allocatable :: upper
      real, dimension (:,:), allocatable :: lower
      real, dimension (:,:), allocatable :: matrix
      integer :: N
      integer :: i, j
      real :: temp
      
      external DGEMM, DPOTRF
      
      N = 16
      
      allocate ( upper(N, N) )
      allocate ( lower(N, N) )
      allocate ( matrix(N, N))
      do i = 1, N
         do j = 1, N
            if (j <= i) then
               temp = int(rand() * 9 + 1)
            else
               temp = 0.0;
            end if
            lower(i, j) = temp
            upper(j, i) = temp   
         end do
      end do
      Print *, "Lower:"
      do i=1, N
         Print 100, ( lower(i,j), j=1,N )
      end do
      
      Print *, "Upper"
      do i=1, N
         Print 100, ( upper(i,j), j=1,N )
      end do

      call DGEMM('N', 'N', N, N, N, 1.0, lower, N, upper, N, 0.0, matrix
     +, N)
      Print *, "Matrix"
      do i=1, N
         Print 100, ( matrix(i,j), j=1,N )
      end do
      
 100  format (1x, 16(1x,f5.1))
      Print *, "Hello World!"
      end program Hello
