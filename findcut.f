c
c
c     ###################################################
c     ## COPYRIGHT(C) 2012 Liam Dennis O'Suilleabhain  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine findr2 -- Assess 2 and 3 body for cutoff r ##
c     ##                                                         ##
c     #############################################################
c
c
c     "analyze" computes and displays the total potential energy;
c     options are provided to display system and force field info,
c     partition the energy by atom or by potential function type,
c     show force field parameters by atom; output the large energy
c     interactions and find electrostatic and inertial properties
c
c
      subroutine findr2
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'molcul.i'
      include 'combo.i'
      integer i,j,epvr,count1,count2
      real*8 xr,yr,zr,r2
 
      moli = 1
      count1 = 1
      count2 = 4
      do moli1 = 1, nmol 
         do moli2 = moli1+1, nmol
            i = count1
            j = count2
            xr = x(i) - x(j)
            yr = y(i) - y(j)
            zr = z(i) - z(j)
            r2 = xr*xr + yr*yr + zr*zr
            r2b(moli)=sqrt(r2)
            moli = moli + 1
            count2 = count2 + 3
            print*,"liamo:findcut",i,j,xr
         end do
         count1 = count1 + 3
         count2 = count1 + 3
      end do
      return
      end

c
c
c     ###################################################
c     ## COPYRIGHT(C) 2012 Liam Dennis O'Suilleabhain  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine findr3 -- Find max radius between 3 bodies  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "analyze" computes and displays the total potential energy;
c     options are provided to display system and force field info,
c     partition the energy by atom or by potential function type,
c     show force field parameters by atom; output the large energy
c     interactions and find electrostatic and inertial properties
c
c
      subroutine findr2
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'molcul.i'
      include 'combo.i'
      integer i,j,k
      real*8 xr1,yr1,zr1,r2_1
      real*8 xr2,yr2,zr2,r2_2
      real*8 xr3,yr3,zr3,r2_3
 
      moli = 1
      do moli1 = 1, nmol 
         do moli2 = moli1+1, nmol
            do moli3 = moli2+1, nmol
               i = imol(1,moli1)
               j = imol(1,moli2)
               k = imol(1,moli3)
               xr1 = x(i) - x(j)
               xr2 = x(i) - x(k)
               xr3 = x(j) - x(k)
               yr1 = y(i) - y(j)
               yr2 = y(i) - y(k)
               yr3 = y(j) - y(k)
               zr1 = z(i) - z(j)
               zr2 = z(i) - z(k)
               zr3 = z(j) - z(k)
               r2_1 = xr1*xr1 + yr1*yr1 + zr1*zr1
               r2_2 = xr2*xr2 + yr2*yr2 + zr2*zr2
               r2_3 = xr3*xr3 + yr3*yr3 + zr3*zr3
               r3b_1(moli) = sqrt(r2_1)
               r3b_2(moli) = sqrt(r2_2)
               r3b_3(moli) = sqrt(r2_3)
               moli = moli + 1
            end do
         end do
      end do
      return
      end



