      program hello

      integer,parameter :: N = 8      
      real :: upper (N, N)
      real :: lower (N, N)
      real :: matrix (N, N)

      integer :: i, j
      real :: temp
      real :: one, zero
      external DGEMM, DPOTRF
      
      one = 1.0
      zero = 0.0

      !allocate ( upper(N, N) )
      !allocate ( lower(N, N) )
      !allocate ( matrix(N, N))
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

      call DGEMM('N', 'N', N, N, N, one, lower, N, upper, N, zero,
     +matrix, N)

      Print *, "Matrix"
c      do i=1, N
c         Print 100, ( matrix(i,j), j=1,N )
c      end do
      
 100  format (1x, 16(1x,f5.1))
      Print *, "Hello World!"
      end program Hello
