c
c     ##################################################################
c     ##                                                              ##
c     ##   subroutine grid_matrix3b -- Evaluate grid Matrix for PME   ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "grid_matrix3b" evaluates the matrix elements of involved in the
c     Matrix for the calculating the induced dipole field with the PME
c     method
c
      subroutine grid_matrix3b
      implicit none
      include 'pme.i'

      integer i,j,k
      integer isite,iatm
      real*8 u0,u1,u2
      real*8 u02,u12,u22
      real*8 v0,v1,v2
      real*8 v02,v12,v22
      real*8 a(3,3)
c
c     convert Cartesian dipoles to fractional coordinates
c
      do i = 1, 3
         a(1,i) = dble(nfft1) * recip(i,1)
         a(2,i) = dble(nfft2) * recip(i,2)
         a(3,i) = dble(nfft3) * recip(i,3)
      end do
     
c  Set up Matrix 

      do isite = 1, npole
         iatm = ipole(isite)
         igrd0 = igrid(1,iatm)
         jgrd0 = igrid(2,iatm)
         kgrd0 = igrid(3,iatm)
         k0 = kgrd0
         do it3 = 1, bsorder
            k0 = k0 + 1
            k = k0 + 1 + (nfft3-isign(nfft3,k0))/2
            v0 = thetai3(1,it3,iatm)
            v02 = v0*v0
            v1 = thetai3(2,it3,iatm)
            v12 = v1*v1
            v2 = thetai3(3,it3,iatm)
            v22 = v2*v2
            v3 = thetai3(4,it3,iatm)
            j0 = jgrd0
            do it2 = 1, bsorder
               j0 = j0 + 1
               j = j0 + 1 + (nfft2-isign(nfft2,j0))/2
               u0 = thetai2(1,it2,iatm)
               u02 = u0*u0)
               u1 = thetai2(2,it2,iatm)
               u12 = u1*u1
               u2 = thetai2(3,it2,iatm)
               u22 = u2*u2
               u3 = thetai2(4,it2,iatm)
               t3 = 0.0d0
               i0 = igrd0
               do it1 = 1, bsorder
                  w0 = thetai2(1,it1,iatm)
                  w1 = thetai2(2,it1,iatm)
                  w2 = thetai2(3,it1,iatm)

                  i0 = i0 + 1
                  i = i0 + 1 + (nfft1-isign(nfft1,i0))/2

                  M11 =  w1*v0*u0
                  M21 =  w1*v0*u0
                  M31 =  w1*v0*u0

                  M1_1(1) = M1_1 + a(1,1)*M11*
                  M1_1(2) = M1_1 + a(1,1)*M21*
                  M1_1(3) = M1_1 + a(1,1)*M31*

                  M1_2(1) = M1_2 + M12*
                  M1_2(2) = M1_2 + M22*
                  M1_2(3) = M1_2 + M32*
                  M1_2(4) = M1_2 + M12*
                  M1_2(5) = M1_2 + M22*
                  M1_2(6) = M1_2 + M32*
                  M1_2(7) = M1_2 + M12*
                  M1_2(8) = M1_2 + M22*
                  M1_2(9) = M1_2 + M32*

                  M1_3(1) = M1_3 + 
                  M1_3(2) = M1_3 + 
                  M1_3(3) = M1_3 + 
                  M1_3(4) = M1_3 + 
                  M1_3(5) = M1_3 + 
                  M1_3(6) = M1_3 + 
                  M1_3(7) = M1_3 + 
                  M1_3(8) = M1_3 + 
                  M1_3(9) = M1_3 + 

                  M2_1(1) = M2_1 +
                  M2_1(2) = M2_1 +
                  M2_1(3) = M2_1 +
                  M2_1(4) = M2_1 + 
                  M2_1(5) = M2_1 + 
                  M2_1(6) = M2_1 + 
                  M2_1(7) = M2_1 + 
                  M2_1(8) = M2_1 + 
                  M2_1(9) = M2_1 + 

                  M2_2(1) = M2_2 + 
                  M2_2(2) = M2_2 + 
                  M2_2(3) = M2_2 + 
                  M2_2(4) = M2_2 + 
                  M2_2(5) = M2_2 + 
                  M2_2(6) = M2_2 + 
                  M2_2(7) = M2_2 + 
                  M2_2(8) = M2_2 + 
                  M2_2(9) = M2_2 + 

                  M2_3(1) = M2_3 + 
                  M2_3(2) = M2_3 + 
                  M2_3(3) = M2_3 + 
                  M2_3(4) = M2_3 + 
                  M2_3(5) = M2_3 + 
                  M2_3(6) = M2_3 + 
                  M2_3(7) = M2_3 + 
                  M2_3(8) = M2_3 + 
                  M2_3(9) = M2_3 + 

                  M3_1(1) = M3_1 + 
                  M3_1(2) = M3_1 + 
                  M3_1(3) = M3_1 + 
                  M3_1(4) = M3_1 + 
                  M3_1(5) = M3_1 + 
                  M3_1(6) = M3_1 + 
                  M3_1(7) = M3_1 + 
                  M3_1(8) = M3_1 + 
                  M3_1(9) = M3_1 + 

                  M3_2(1) = M3_2 + 
                  M3_2(2) = M3_2 + 
                  M3_2(3) = M3_2 + 
                  M3_2(4) = M3_2 + 
                  M3_2(5) = M3_2 + 
                  M3_2(6) = M3_2 + 
                  M3_2(7) = M3_2 + 
                  M3_2(8) = M3_2 + 
                  M3_2(9) = M3_2 + 

                  M3_3(1) = M3_3 + 
                  M3_3(2) = M3_3 + 
                  M3_3(3) = M3_3 + 
                  M3_3(4) = M3_3 + 
                  M3_3(5) = M3_3 + 
                  M3_3(6) = M3_3 + 
                  M3_3(7) = M3_3 + 
                  M3_3(8) = M3_3 + 
                  M3_3(9) = M3_3 + 

               end do
            end do
         end do
      end do
 
