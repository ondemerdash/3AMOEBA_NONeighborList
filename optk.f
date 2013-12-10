c      subroutine optkvec(a,optvec)
      program optk
      real*8 a
      real*8 k,optvec,optvec2
      real*8 eps,box,cut
      eps=1.0d-8
      k=15.0d0
      optvec=0.0d0
c      a=0.5d0
      box=20.0d0
      cut=20.0d0
      call ewaldcof(a,cut)
      print*,"aewald=",a      
c      ratio = eps + 1.0d0
c      x = 0.5d0
c      i = 0
      optvec=k-( exp(-k*k/4/a/a)-eps*k*k)/
     & (exp(-k*k/4/a/a)*(-2*k/4/a/a)-2*eps*k)
      diff=optvec-k
      do while (abs(diff) .ge. eps)
         
c         i = i + 1
c         x = 2.0d0 * x
c         y = x * cutoff
c         ratio = erfc(y) / cutoff
      optvec=k-( exp(-k*k/4/a/a)-eps*k*k)/
     & (exp(-k*k/4/a/a)*(-2*k/4/a/a)-2*eps*k)

        diff=optvec-k
        k=optvec
      end do
      optvec2=optvec*box/2/3.14
      print*,"Optvec=",optvec 
      print*,"Optvec2=",optvec2
      optvec_frenkel=cut*a*box*a/3.14
      print*,"Optvec_frenkel=",optvec_frenkel
      end

      subroutine ewaldcof (alpha,cutoff)
      implicit none
      integer i,k
      real*8 alpha,cutoff,eps
      real*8 x,xlo,xhi,y
      real*8 ratio,erfc
      external erfc
c
c
c     set the tolerance value; use of 1.0d-8 results in large
c     Ewald coefficients that ensure continuity in the gradient
c
      eps = 1.0d-8
c
c     get approximate value from cutoff and tolerance
c
      ratio = eps + 1.0d0
      x = 0.5d0
      i = 0
      do while (ratio .ge. eps)
         i = i + 1
         x = 2.0d0 * x
         y = x * cutoff
         ratio = erfc(y) / cutoff
      end do
c
c     use a binary search to refine the coefficient
c
      k = i + 60
      xlo = 0.0d0
      xhi = x
      do i = 1, k
         x = (xlo+xhi) / 2.0d0
         y = x * cutoff
         ratio = erfc(y) / cutoff
         if (ratio .ge. eps) then
            xlo = x
         else
            xhi = x
         end if
      end do
      alpha = x
      return
      end

      function erfc (x)
      implicit none
      integer mode
      real*8 erfc,x
      real*8 result
c
c
c     get the complementary error function via Chebyshev fitting
c
      mode = 1
      call erfcore (x,result,mode)
      erfc = result
      return
      end

      subroutine erfcore (arg,result,mode)
      implicit none
      integer i,mode
      real*8 arg,result
      real*8 x,y,ysq
      real*8 del,sqrpi
      real*8 xnum,xden
      real*8 xtiny,xbig
      real*8 a(5),b(4)
      real*8 c(9),d(8)
      real*8 p(6),q(5)
c
c     mathematical and machine-dependent constants
c
      data sqrpi  / 5.6418958354775628695d-1 /
      data xtiny  / 1.11d-16 /
      data xbig   / 26.543d0 /
c
c     coefficients for approximation to erf in first interval
c
      data a  / 3.16112374387056560d0,  1.13864154151050156d2,
     &          3.77485237685302021d2,  3.20937758913846947d3,
     &          1.85777706184603153d-1 /
      data b  / 2.36012909523441209d1,  2.44024637934444173d2,
     &          1.28261652607737228d3,  2.84423683343917062d3 /
c
c     coefficients for approximation to erfc in second interval
c
      data c  / 5.64188496988670089d-1, 8.88314979438837594d0,
     &          6.61191906371416295d1,  2.98635138197400131d2,
     &          8.81952221241769090d2,  1.71204761263407058d3,
     &          2.05107837782607147d3,  1.23033935479799725d3,
     &          2.15311535474403846d-8 /
      data d  / 1.57449261107098347d1,  1.17693950891312499d2,
     &          5.37181101862009858d2,  1.62138957456669019d3,
     &          3.29079923573345963d3,  4.36261909014324716d3,
     &          3.43936767414372164d3,  1.23033935480374942d3 /
c
c     coefficients for approximation to erfc in third interval
c
      data p  / 3.05326634961232344d-1, 3.60344899949804439d-1,
     &          1.25781726111229246d-1, 1.60837851487422766d-2,
     &          6.58749161529837803d-4, 1.63153871373020978d-2 /
      data q  / 2.56852019228982242d0,  1.87295284992346047d0,
     &          5.27905102951428412d-1, 6.05183413124413191d-2,
     &          2.33520497626869185d-3 /
c
c
c     store the argument and its absolute value
c
      x = arg
      y = abs(x)
c
c     evaluate error function for |x| less than 0.46875
c
      if (y .le. 0.46875d0) then
         ysq = 0.0d0
         if (y .gt. xtiny)  ysq = y * y
         xnum = a(5) * ysq
         xden = ysq
         do i = 1, 3
            xnum = (xnum + a(i)) * ysq
            xden = (xden + b(i)) * ysq
         end do
         result = x * (xnum + a(4)) / (xden + b(4))
         if (mode .ne. 0)  result = 1.0d0 - result
c
c     get complementary error function for 0.46875 <= |x| <= 4.0
c
      else if (y .le. 4.0d0) then
         xnum = c(9) * y
         xden = y
         do i = 1, 7
            xnum = (xnum + c(i)) * y
            xden = (xden + d(i)) * y
         end do
         result = (xnum + c(8)) / (xden + d(8))
         ysq = aint(16.0d0*y) / 16.0d0
         del = (y-ysq) * (y+ysq)
c        result = exp(-ysq*ysq) * exp(-del) * result
         result = exp(-ysq*ysq-del) * result
         if (mode .eq. 0) then
            result = 1.0d0 - result
            if (x .lt. 0.0d0)  result = -result
         else
            if (x .lt. 0.0d0)  result = 2.0d0 - result
         end if
c
c     get complementary error function for |x| greater than 4.0
c
      else
         result = 0.0d0
         if (y .lt. xbig) then
            ysq = 1.0d0 / (y * y)
            xnum = p(6) * ysq
            xden = ysq
            do i = 1, 4
               xnum = (xnum + p(i)) * ysq
               xden = (xden + q(i)) * ysq
            end do
            result = ysq * (xnum + p(5)) / (xden + q(5))
            result = (sqrpi - result) / y
            ysq = aint(16.0d0*y) / 16.0d0
            del = (y-ysq) * (y+ysq)
c           result = exp(-ysq*ysq) * exp(-del) * result
            result = exp(-ysq*ysq-del) * result
         end if
         if (mode .eq. 0) then
            result = 1.0d0 - result
            if (x .lt. 0.0d0)  result = -result
         else
            if (x .lt. 0.0d0)  result = 2.0d0 - result
         end if
      end if
      return
      end

