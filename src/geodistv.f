      subroutine geodistv(lat1,lon1,lat2,lon2,dist, n)
      double precision lat1(n), lon1(n), lat2(n), lon2(n), dist(n)
      integer n

      double precision faz, baz
      double precision a,f
      integer i

      a = 6378137.00 
      f = 1./298.257223563
      do 200 i=1,n
      call geoddist(lat1(i), lon1(i), lat2(i), lon2(i), a, f, 
     1 faz, baz, dist(i))
 200  continue
      end


