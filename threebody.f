c
c
c     #############################################################
c     ##  COPYRIGHT (C) 2012 Liam Denis O'Suilleabhain           ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  routines below find the respective numbers of atoms in    ##
c     ##  for the interactions between combinations of molecules    ##
c     ##                                                            ##
c     ################################################################
c
c     #############################################################
c     ##                                                         ##
c     ##       subroutine combo1  --  combinations of atoms      ##
c     ##                                                         ##
c     #############################################################
c
c
c     "combo1" finds the index of atoms for which the interactions
c     are being evaluated in empole3
c
c
      subroutine combo1
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      include 'molcul.i'
      include 'combo.i'
      integer i,ii,index
      integer hatom1,latom1
      integer hpole1,lpole1
      logical goon
c
c Find the pole numbers and atom indices
c for each of the 1 body interaction
c
      goon = .true.
      latom1 = imol(1,moli1)
      hatom1 = imol(2,moli1)
      do i = 1, npole
         ii = ipole(i)
         if (ii .ge. latom1 .and. goon) then
            lpole1 = i
            hpole1 = lpole1
            goon = .false.
          else 
            if (ii .le. hatom1 .and.
     &       ii .gt. latom1) hpole1 = hpole1+1
          end if
      end do
      np1 = (hpole1 - lpole1) + 1
      index = 0
      do i = 1, np1
         pnum(i) = lpole1 + index
         index = index + 1
      end do
      npole3b = np1
      return
      end

c     #############################################################
c     ##                                                         ##
c     ##       subroutine combo2  --  combinations of atoms      ##
c     ##                                                         ##
c     #############################################################
c
c
c     "combo2" finds the index of atoms for which the interactions
c     are being evaluated in empole3
c
      subroutine combo2
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      include 'molcul.i'
      include 'combo.i'
      integer i,ii,index
      integer hatom1,hatom2
      integer hpole1,hpole2
      integer latom1,latom2
      integer lpole1,lpole2
      logical goon

      goon = .true.
      latom2 = imol(1,moli2)
      hatom2 = imol(2,moli2)
      do i = 1, npole
         ii = ipole(i)
         if (ii .ge. latom2 .and. goon) then
            lpole2 = i
            hpole2 = lpole2
            goon = .false.
         else
            if (ii .le. hatom2 .and.
     &          ii .gt. latom2) hpole2 = hpole2+1
         end if
      end do
 
      np2 = np1 + (hpole2 - lpole2) + 1
      index = 0
      do i = (np1+1), np2
          pnum(i) = lpole2 + index
          index = index + 1
      end do
      npole3b = np2
      return
      end



c     #############################################################
c     ##                                                         ##
c     ##       subroutine combo2new   combinations of atoms      ##
c     ##                                                         ##
c     #############################################################
c
c
c
c     "combo2new" finds the index of atoms for which the interactions
c     are being evaluated in empole3
c
      subroutine combo2new
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      include 'molcul.i'
      include 'combo.i'
      integer i,ii,index
      integer hatom1,hatom2
      integer hpole1,hpole2
      integer latom1,latom2
      integer lpole1,lpole2
      logical goon

c
c Find the pole numbers and atom indices
c for each of the 1 body interaction
c
      goon = .true.
      latom1 = imol(1,moli1)
      hatom1 = imol(2,moli1)
      do i = 1, npole
         ii = ipole(i)
         if (ii .ge. latom1 .and. goon) then
            lpole1 = i
            hpole1 = lpole1
            goon = .false.
          else 
            if (ii .le. hatom1 .and.
     &       ii .gt. latom1) hpole1 = hpole1+1
          end if
      end do
      np1 = (hpole1 - lpole1) + 1
      index = 0
      do i = 1, np1
         pnum(i) = lpole1 + index
         index = index + 1
      end do

      goon = .true.
      latom2 = imol(1,moli2)
      hatom2 = imol(2,moli2)
      do i = 1, npole
         ii = ipole(i)
         if (ii .ge. latom2 .and. goon) then
            lpole2 = i
            hpole2 = lpole2
            goon = .false.
         else
            if (ii .le. hatom2 .and.
     &          ii .gt. latom2) hpole2 = hpole2+1
         end if
      end do

      np2 = np1 + (hpole2 - lpole2) + 1
      index = 0
      do i = (np1+1), np2
          pnum(i) = lpole2 + index
          index = index + 1
      end do
      npole3b = np2
      return
      end



c     #############################################################
c     ##                                                         ##
c     ##       subroutine combo3  --  combinations of atoms      ##
c     ##                                                         ##
c     #############################################################
c
c
c     "combo3" finds the index of atoms for which the interactions
c     are being evaluated in empole3
c
      subroutine combo3
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      include 'molcul.i'
      include 'combo.i'
      integer i,ii,index
      integer hatom3, latom3
      integer hpole3, lpole3
      logical goon
  
c
c Get atomic indices of molecule 3
c
      goon = .true.
      latom3 = imol(1,moli3)
      hatom3 = imol(2,moli3)
      do i = 1, npole
         ii = ipole(i)
         if (ii .ge. latom3 .and. goon) then
            lpole3 = i
            hpole3 = lpole3
            goon = .false.
         else
            if (ii .le. hatom3 .and.
     &          ii .gt. latom3) hpole3 = hpole3+1
         end if
      end do
 
c      print*,"lpole3=",lpole3,"hpole3=",hpole3
      np3 = np2 + (hpole3 - lpole3) + 1
      index = 0
      do i = (np2+1), np3
          pnum(i) = lpole3 + index
          index = index + 1
      end do
      npole3b = np3
      return
      end


c     #############################################################
c     ##                                                         ##
c     ##       subroutine combo3new   combinations of atoms      ##
c     ##                                                         ##
c     #############################################################
c
c
c     "combo3new" finds the index of atoms for which the interactions
c     are being evaluated in empole3
c
      subroutine combo3new
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      include 'molcul.i'
      include 'combo.i'
      integer i,ii,index
      integer hatom3, latom3
      integer hpole3, lpole3
      integer hatom1,hatom2
      integer hpole1,hpole2
      integer latom1,latom2
      integer lpole1,lpole2
      logical goon

c
c Get atomic indices of molecule 3
c
      goon = .true.
      latom1 = imol(1,moli1)
      hatom1 = imol(2,moli1)
      do i = 1, npole
         ii = ipole(i)
         if (ii .ge. latom1 .and. goon) then
            lpole1 = i
            hpole1 = lpole1
            goon = .false.
          else 
            if (ii .le. hatom1 .and.
     &       ii .gt. latom1) hpole1 = hpole1+1
          end if
      end do
      np1 = (hpole1 - lpole1) + 1
      index = 0
      do i = 1, np1
         pnum(i) = lpole1 + index
         index = index + 1
      end do

      goon = .true.
      latom2 = imol(1,moli2)
      hatom2 = imol(2,moli2)
      do i = 1, npole
         ii = ipole(i)
         if (ii .ge. latom2 .and. goon) then
            lpole2 = i
            hpole2 = lpole2
            goon = .false.
         else
            if (ii .le. hatom2 .and.
     &          ii .gt. latom2) hpole2 = hpole2+1
         end if
      end do

      np2 = np1 + (hpole2 - lpole2) + 1
      index = 0
      do i = (np1+1), np2
          pnum(i) = lpole2 + index
          index = index + 1
      end do

      goon = .true.
      latom3 = imol(1,moli3)
      hatom3 = imol(2,moli3)
      do i = 1, npole
         ii = ipole(i)
         if (ii .ge. latom3 .and. goon) then
            lpole3 = i
            hpole3 = lpole3
            goon = .false.
         else
            if (ii .le. hatom3 .and.
     &          ii .gt. latom3) hpole3 = hpole3+1
         end if
      end do

c      print*,"lpole3=",lpole3,"hpole3=",hpole3
      np3 = np2 + (hpole3 - lpole3) + 1
      index = 0
      do i = (np2+1), np3
          pnum(i) = lpole3 + index
          index = index + 1
      end do
      npole3b = np3
      return
      end

c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine findr2 -- Assess 2 and 3 body for cutoff r ##
c     ##                                                         ##
c     #############################################################
c
c
c
      subroutine findr2 (r_123)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'molcul.i'
      include 'combo.i'
      include 'atmtyp.i'
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2
      real*8 xcm,ycm,zcm,r_cm,r2,r_123
      real*8 xr,yr,zr
      integer i,j,l1,mol1,mol2,M1,M2,M

      xcm1 = 0.0d0
      ycm1 = 0.0d0
      zcm1 = 0.0d0
      xcm2 = 0.0d0
      ycm2 = 0.0d0
      zcm2 = 0.0d0
      xcm = 0.0d0
      ycm = 0.0d0
      zcm = 0.0d0
      M1 = 0.0d0
      M2 = 0.0d0
      M = 0.0d0
c
c Find center of mass of each molecule
c
      do l1 = 1, np1
         i = pnum(l1)
         xcm1 = xcm1 + (mass(i)*x(i))
         ycm1 = ycm1 + (mass(i)*y(i))
         zcm1 = zcm1 + (mass(i)*z(i))
         M1 = M1 + mass(i)
         xcm = xcm + (mass(i)*x(i))
         ycm = ycm + (mass(i)*y(i))
         zcm = zcm + (mass(i)*z(i))
         M = M + mass(i)
      end do

      do l1 = np1+1, np2
         i = pnum(l1)
         xcm2 = xcm2 + (mass(i)*x(i))
         ycm2 = ycm2 + (mass(i)*y(i))
         zcm2 = zcm2 + (mass(i)*z(i))
         M2 = M2 + mass(i)
         xcm = xcm + (mass(i)*x(i))
         ycm = ycm + (mass(i)*y(i))
         zcm = zcm + (mass(i)*z(i))
         M = M + mass(i)
      end do

      xcm1 = xcm1/M1
      ycm1 = ycm1/M1
      zcm1 = zcm1/M1

      xcm2 = xcm2/M2
      ycm2 = ycm2/M2
      zcm2 = zcm2/M2

c
c Find center of mass of system
c
c      do l1 = 1, npole3b
c         i = pnum(l1)
c         xcm = xcm + (mass(i)*x(i))
c         ycm = ycm + (mass(i)*y(i))
c         zcm = zcm + (mass(i)*z(i))
c         M = M + mass(i)
c      end do

      xcm = xcm/M
      ycm = ycm/M
      zcm = zcm/M
c
c Get distance between the two molecular
c COMs in reference to the system COM
c
      r_123 = ((xcm1 - xcm)*(xcm1 - xcm)
     &            + (ycm1 - ycm)*(ycm1 - ycm)
     &            + (zcm1 - zcm)*(zcm1 - zcm))
     &        + ((xcm2 - xcm)*(xcm2 - xcm)
     &            + (ycm2 - ycm)*(ycm2 - ycm)
     &            + (zcm2 - zcm)*(zcm2 - zcm))
      r_123 = r_123/2.0d0

      return
      end


c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine findr2 -- Assess 2 and 3 body for cutoff r ##
c     ##                                                         ##
c     #############################################################
c
c
c
      subroutine findr2_minimage (r_123)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'molcul.i'
      include 'combo.i'
      include 'atmtyp.i'
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2
      real*8 xcm,ycm,zcm,r_cm,r2,r_123
      real*8 xr,yr,zr,xcm_tmp1,ycm_tmp1
      real*8 zcm_tmp1,xr1,yr1,zr1,xr2,yr2,zr2
      integer i,j,l1,mol1,mol2,M1,M2,M

      xcm1 = 0.0d0
      ycm1 = 0.0d0
      zcm1 = 0.0d0
      xcm2 = 0.0d0
      ycm2 = 0.0d0
      zcm2 = 0.0d0
      xcm = 0.0d0
      ycm = 0.0d0
      zcm = 0.0d0
      M1 = 0.0d0
      M2 = 0.0d0
      M = 0.0d0
c
c Find center of mass of each molecule
c
      do l1 = 1, np1
         i = pnum(l1)
         xcm1 = xcm1 + (mass(i)*x(i))
         ycm1 = ycm1 + (mass(i)*y(i))
         zcm1 = zcm1 + (mass(i)*z(i))
         M1 = M1 + mass(i)
         xcm = xcm + (mass(i)*x(i))
         ycm = ycm + (mass(i)*y(i))
         zcm = zcm + (mass(i)*z(i))
         M = M + mass(i)
      end do

      do l1 = np1+1, np2
         i = pnum(l1)
         xcm2 = xcm2 + (mass(i)*x(i))
         ycm2 = ycm2 + (mass(i)*y(i))
         zcm2 = zcm2 + (mass(i)*z(i))
         M2 = M2 + mass(i)
         xcm = xcm + (mass(i)*x(i))
         ycm = ycm + (mass(i)*y(i))
         zcm = zcm + (mass(i)*z(i))
         M = M + mass(i)
      end do

      xcm1 = xcm1/M1
      ycm1 = ycm1/M1
      zcm1 = zcm1/M1

      xcm2 = xcm2/M2
      ycm2 = ycm2/M2
      zcm2 = zcm2/M2

      xcm = xcm/M
      ycm = ycm/M
      zcm = zcm/M

      xr = xcm1 - xcm
      yr = ycm1 - ycm
      zr = zcm1 - zcm

      call image(xr,yr,zr)
      
      xcm_tmp1 = xcm1-xr
      ycm_tmp1 = ycm1-yr
      zcm_tmp1 = zcm1-zr

      xr1=xr
      yr1=yr
      zr1=zr

      xr = xcm2 - xcm_tmp1
      yr = ycm2 - ycm_tmp1
      zr = zcm2 - zcm_tmp1

      call image(xr,yr,zr)

      xr2=xr
      yr2=yr
      zr2=zr


c
c Get distance between the two molecular
c COMs in reference to the system COM
c
c      r_123 = ((xcm1 - xcm)*(xcm1 - xcm)
c     &            + (ycm1 - ycm)*(ycm1 - ycm)
c     &            + (zcm1 - zcm)*(zcm1 - zcm))
c     &        + ((xcm2 - xcm)*(xcm2 - xcm)
c     &            + (ycm2 - ycm)*(ycm2 - ycm)
c     &            + (zcm2 - zcm)*(zcm2 - zcm))
      r_123=(xr1*xr1+yr1*yr1+zr1*zr1+xr2*xr2+
     &      +yr2*yr2+zr2*zr2)/2.0d0

c      r_123 = r_123/2.0d0

      return
      end

c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine findr2 -- Assess 2 and 3 body for cutoff r ##
c     ##                                                         ##
c     #############################################################
c
c
c
      subroutine findr2_minimage2 (r_123)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
c      include 'molcul.i'
      include 'combo.i'

      include 'atmtyp.i'
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2
      real*8 xcm,ycm,zcm,r_cm,r2,r_123
      real*8 xr,yr,zr,xcm_tmp1,ycm_tmp1
      real*8 zcm_tmp1,xr1,yr1,zr1,xr2,yr2,zr2
      integer i,j,l1,mol1,mol2,M1,M2,M

      xcm1 = 0.0d0
      ycm1 = 0.0d0
      zcm1 = 0.0d0
      xcm2 = 0.0d0
      ycm2 = 0.0d0
      zcm2 = 0.0d0
      xcm = 0.0d0
      ycm = 0.0d0
      zcm = 0.0d0
      M1 = 0.0d0
      M2 = 0.0d0
      M = 0.0d0
c
c Find center of mass of each molecule
c
      do l1 = 1, np1
         i = pnum(l1)
         xcm1 = xcm1 + (mass(i)*x(i))
         ycm1 = ycm1 + (mass(i)*y(i))
         zcm1 = zcm1 + (mass(i)*z(i))
         M1 = M1 + mass(i)
c         xcm = xcm + (mass(i)*x(i))
c         ycm = ycm + (mass(i)*y(i))
c         zcm = zcm + (mass(i)*z(i))
c         M = M + mass(i)
      end do

      do l1 = np1+1, np2
         i = pnum(l1)
         xcm2 = xcm2 + (mass(i)*x(i))
         ycm2 = ycm2 + (mass(i)*y(i))
         zcm2 = zcm2 + (mass(i)*z(i))
         M2 = M2 + mass(i)
c         xcm = xcm + (mass(i)*x(i))
c         ycm = ycm + (mass(i)*y(i))
c         zcm = zcm + (mass(i)*z(i))
c         M = M + mass(i)
      end do

      xcm1 = xcm1/M1
      ycm1 = ycm1/M1
      zcm1 = zcm1/M1

      xcm2 = xcm2/M2
      ycm2 = ycm2/M2
      zcm2 = zcm2/M2

c      xcm = xcm/M
c      ycm = ycm/M
c      zcm = zcm/M

      xr = xcm1 - xcm2
      yr = ycm1 - ycm2
      zr = zcm1 - zcm2

      call image(xr,yr,zr)


c
c Get distance between the two molecular
c COMs in reference to the system COM
c
c      r_123 = ((xcm1 - xcm)*(xcm1 - xcm)
c     &            + (ycm1 - ycm)*(ycm1 - ycm)
c     &            + (zcm1 - zcm)*(zcm1 - zcm))
c     &        + ((xcm2 - xcm)*(xcm2 - xcm)
c     &            + (ycm2 - ycm)*(ycm2 - ycm)
c     &            + (zcm2 - zcm)*(zcm2 - zcm))
c      r_123=(xr1*xr1+yr1*yr1+zr1*zr1+xr2*xr2+
c     &      +yr2*yr2+zr2*zr2)/2.0d0
      r_123=xr*xr + yr*yr + zr*zr

c      r_123 = r_123/2.0d0

      return
      end

c
      function findr2_minimage2func (pnum,np1,np2)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
c      include 'molcul.i'
c      include 'combo.i'
      include 'atmtyp.i'
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2
      real*8 xcm,ycm,zcm,r_cm,r2,r_123
      real*8 xr,yr,zr,xcm_tmp1,ycm_tmp1
      real*8 zcm_tmp1,xr1,yr1,zr1,xr2,yr2,zr2
      real*8 findr2_minimage2func
      integer i,j,l1,mol1,mol2,M1,M2,M
      integer np1,np2
      integer pnum(*)
      xcm1 = 0.0d0
      ycm1 = 0.0d0
      zcm1 = 0.0d0
      xcm2 = 0.0d0
      ycm2 = 0.0d0
      zcm2 = 0.0d0
      xcm = 0.0d0
      ycm = 0.0d0
      zcm = 0.0d0
      M1 = 0.0d0
      M2 = 0.0d0
      M = 0.0d0
c
c Find center of mass of each molecule
c
      do l1 = 1, np1
         i = pnum(l1)
         xcm1 = xcm1 + (mass(i)*x(i))
         ycm1 = ycm1 + (mass(i)*y(i))
         zcm1 = zcm1 + (mass(i)*z(i))
         M1 = M1 + mass(i)
c         xcm = xcm + (mass(i)*x(i))
c         ycm = ycm + (mass(i)*y(i))
c         zcm = zcm + (mass(i)*z(i))
c         M = M + mass(i)
      end do

      do l1 = np1+1, np2
         i = pnum(l1)
         xcm2 = xcm2 + (mass(i)*x(i))
         ycm2 = ycm2 + (mass(i)*y(i))
         zcm2 = zcm2 + (mass(i)*z(i))
         M2 = M2 + mass(i)
c         xcm = xcm + (mass(i)*x(i))
c         ycm = ycm + (mass(i)*y(i))
c         zcm = zcm + (mass(i)*z(i))
c         M = M + mass(i)
      end do

      xcm1 = xcm1/M1
      ycm1 = ycm1/M1
      zcm1 = zcm1/M1

      xcm2 = xcm2/M2
      ycm2 = ycm2/M2
      zcm2 = zcm2/M2

c      xcm = xcm/M
c      ycm = ycm/M
c      zcm = zcm/M

      xr = xcm1 - xcm2
      yr = ycm1 - ycm2
      zr = zcm1 - zcm2

      call image(xr,yr,zr)


c
c Get distance between the two molecular
c COMs in reference to the system COM
c
c      r_123 = ((xcm1 - xcm)*(xcm1 - xcm)
c     &            + (ycm1 - ycm)*(ycm1 - ycm)
c     &            + (zcm1 - zcm)*(zcm1 - zcm))
c     &        + ((xcm2 - xcm)*(xcm2 - xcm)
c     &            + (ycm2 - ycm)*(ycm2 - ycm)
c     &            + (zcm2 - zcm)*(zcm2 - zcm))
c      r_123=(xr1*xr1+yr1*yr1+zr1*zr1+xr2*xr2+
c     &      +yr2*yr2+zr2*zr2)/2.0d0
      r_123=xr*xr + yr*yr + zr*zr

c      r_123 = r_123/2.0d0
      findr2_minimage2func=r_123
      return
      end



c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine findr2 -- Assess 2 and 3 body for cutoff r ##
c     ##                                                         ##
c     #############################################################
c
c
c
      subroutine findr2new (r_123)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'molcul.i'
      include 'combo.i'
      include 'atmtyp.i'
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2
      real*8 xcm,ycm,zcm,r_cm,r2,r_123
      real*8 xr,yr,zr
      integer i,j,l1,mol1,mol2,M1,M2,M

      xcm1 = 0.0d0
      ycm1 = 0.0d0
      zcm1 = 0.0d0
      xcm2 = 0.0d0
      ycm2 = 0.0d0
      zcm2 = 0.0d0
      xcm = 0.0d0
      ycm = 0.0d0
      zcm = 0.0d0
      M1 = 0.0d0
      M2 = 0.0d0
      M = 0.0d0
c
c Find center of mass of each molecule
c
      do l1 = 1, np1
         i = pnum(l1)
         xcm1 = xcm1 + (mass(i)*x(i))
         ycm1 = ycm1 + (mass(i)*y(i))
         zcm1 = zcm1 + (mass(i)*z(i))
         M1 = M1 + mass(i)
         xcm = xcm + (mass(i)*x(i))
         ycm = ycm + (mass(i)*y(i))
         zcm = zcm + (mass(i)*z(i))
         M = M + mass(i)
      end do

      do l1 = np1+1, np2
         i = pnum(l1)
         xcm2 = xcm2 + (mass(i)*x(i))
         ycm2 = ycm2 + (mass(i)*y(i))
         zcm2 = zcm2 + (mass(i)*z(i))
         M2 = M2 + mass(i)
         xcm = xcm + (mass(i)*x(i))
         ycm = ycm + (mass(i)*y(i))
         zcm = zcm + (mass(i)*z(i))
         M = M + mass(i)
      end do

      xcm1 = xcm1/M1
      ycm1 = ycm1/M1
      zcm1 = zcm1/M1

      xcm2 = xcm2/M2
      ycm2 = ycm2/M2
      zcm2 = zcm2/M2

c
c Find center of mass of system
c
c      do l1 = 1, npole3b
c         i = pnum(l1)
c         xcm = xcm + (mass(i)*x(i))
c         ycm = ycm + (mass(i)*y(i))
c         zcm = zcm + (mass(i)*z(i))
c         M = M + mass(i)
c      end do

      xcm = xcm/M
      ycm = ycm/M
      zcm = zcm/M
c
c Get distance between the two molecular
c COMs in reference to the system COM
c
      r_123 = ((xcm1 - xcm)*(xcm1 - xcm)
     &            + (ycm1 - ycm)*(ycm1 - ycm)
     &            + (zcm1 - zcm)*(zcm1 - zcm))
     &        + ((xcm2 - xcm)*(xcm2 - xcm)
     &            + (ycm2 - ycm)*(ycm2 - ycm)
     &            + (zcm2 - zcm)*(zcm2 - zcm))
      r_123 = r_123/2.0d0

      return
      end




c     ##                                                         ##
c     #############################################################
c
c
c
      subroutine findr2geom (r_123)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'molcul.i'
      include 'combo.i'
      include 'atmtyp.i'
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2
      real*8 xcm,ycm,zcm,r_cm,r2,r_123
      real*8 xr,yr,zr
      integer i,j,l1,mol1,mol2,M1,M2,M

      xcm1 = 0.0d0
      ycm1 = 0.0d0
      zcm1 = 0.0d0
      xcm2 = 0.0d0
      ycm2 = 0.0d0
      zcm2 = 0.0d0
      xcm = 0.0d0
      ycm = 0.0d0
      zcm = 0.0d0
      M1 = 0.0d0
      M2 = 0.0d0
      M = 0.0d0
c
c Find center of mass of each molecule
c
      do l1 = 1, np1
         i = pnum(l1)
          xcm1 = xcm1 + x(i)
          ycm1 = ycm1 + y(i)
          zcm1 = zcm1 + z(i)
          xcm = xcm +x(i)
          ycm = ycm +y(i)
          zcm = zcm +z(i)
      end do

      do l1 = np1+1, np2
         i = pnum(l1)
         xcm2 = xcm2 + x(i)
         ycm2 = ycm2 + y(i)
         zcm2 = zcm2 + z(i)
          xcm = xcm +x(i)
          ycm = ycm +y(i)
          zcm = zcm +z(i)

      end do

      xcm1 = xcm1/3.0d0
      ycm1 = ycm1/3.0d0
      zcm1 = zcm1/3.0d0

      xcm2 = xcm2/3.0d0
      ycm2 = ycm2/3.0d0
      zcm2 = zcm2/3.0d0

      xcm = xcm/6.0d0
      ycm = ycm/6.0d0
      zcm = zcm/6.0d0

c
c Find center of mass of system
c

c
c Get distance between the two molecular
c COMs in reference to the system COM
c
      r_123 = ((xcm1 - xcm)*(xcm1 - xcm)
     &            + (ycm1 - ycm)*(ycm1 - ycm)
     &            + (zcm1 - zcm)*(zcm1 - zcm))
     &        + ((xcm2 - xcm)*(xcm2 - xcm)
     &            + (ycm2 - ycm)*(ycm2 - ycm)
     &            + (zcm2 - zcm)*(zcm2 - zcm))
      r_123 = r_123/2.0d0

      return
      end

c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine findr2 -- Assess 2 and 3 body for cutoff r ##
c     ##                                                         ##
c     #############################################################
c
c
c
      subroutine findr2test (r_123)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'molcul.i'
      include 'combo.i'
      include 'atmtyp.i'
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2
      real*8 xcm,ycm,zcm,r_cm,r2,r_123
      real*8 xr,yr,zr
      integer i,j,l1,mol1,mol2,M1,M2,M

      xcm1 = 0.0d0
      ycm1 = 0.0d0
      zcm1 = 0.0d0
      xcm2 = 0.0d0
      ycm2 = 0.0d0
      zcm2 = 0.0d0
      xcm = 0.0d0
      ycm = 0.0d0
      zcm = 0.0d0
      M1 = 0.0d0
      M2 = 0.0d0
      M = 0.0d0
c
c Find center of mass of each molecule
c
c      do l1 = 1, np1
c         i = pnum(l1)
c         xcm1 = xcm1 + (mass(i)*x(i))
c         ycm1 = ycm1 + (mass(i)*y(i))
c         zcm1 = zcm1 + (mass(i)*z(i))
c         M1 = M1 + mass(i)
c      end do

c      do l1 = np1+1, np2
c         i = pnum(l1)
c         xcm2 = xcm2 + (mass(i)*x(i))
c         ycm2 = ycm2 + (mass(i)*y(i))
c         zcm2 = zcm2 + (mass(i)*z(i))
c         M2 = M2 + mass(i)
c      end do

      i=pnum(1)
      xcm1 = x(i)
      ycm1 = y(i)
      zcm1 = z(i)

      i=pnum(2)
      xcm2 = x(i)
      ycm2 = y(i)
      zcm2 = z(i)

c
c Find center of mass of system
c
      do l1 = 1, npole3b
         i = pnum(l1)
c         xcm = xcm + (mass(i)*x(i))
c         ycm = ycm + (mass(i)*y(i))
c         zcm = zcm + (mass(i)*z(i))
c         M = M + mass(i)
         xcm = xcm + x(i)
         ycm = ycm + y(i)
         zcm = zcm + z(i)
      end do

      xcm = xcm/2.0d0
      ycm = ycm/2.0d0
      zcm = zcm/2.0d0
c
c Get distance between the two molecular
c COMs in reference to the system COM
c
      r_123 = ((xcm1 - xcm)*(xcm1 - xcm)
     &            + (ycm1 - ycm)*(ycm1 - ycm)
     &            + (zcm1 - zcm)*(zcm1 - zcm))
     &        + ((xcm2 - xcm)*(xcm2 - xcm)
     &            + (ycm2 - ycm)*(ycm2 - ycm)
     &            + (zcm2 - zcm)*(zcm2 - zcm))
      r_123 = r_123/2.0d0

      return
      end



c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine findr3 -- Find max radius between 3 bodies  ##
c     ##                                                         ##
c     #############################################################
c
c     findr3 finds the distance of each particle from the center of
c     mass in each direction. The radius of interaction of three 
c     bodies is defined as the sum of the distances of each particle
c     from the center of mass.
c
      subroutine findr3 (r_123)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'molcul.i'
      include 'combo.i'
      include 'atmtyp.i'
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2
      real*8 xcm3,ycm3,zcm3,xr,yr,zr
      real*8 xcm,ycm,zcm,r_cm,r2,r_123
      integer i,j,l1,mol1,mol2,M1,M2,M3,M


      xcm1 = 0.0d0
      ycm1 = 0.0d0
      zcm1 = 0.0d0
      xcm2 = 0.0d0
      ycm2 = 0.0d0
      zcm2 = 0.0d0
      xcm3 = 0.0d0
      ycm3 = 0.0d0
      zcm3 = 0.0d0
      xcm = 0.0d0
      ycm = 0.0d0
      zcm = 0.0d0
      M1 = 0.0d0
      M2 = 0.0d0
      M3 = 0.0d0
      M = 0.0d0
c
c Find center of mass of each molecule
c
      do l1 = 1, np1
         i = pnum(l1)
         xcm1 = xcm1 + (mass(i)*x(i))
         ycm1 = ycm1 + (mass(i)*y(i))
         zcm1 = zcm1 + (mass(i)*z(i))
         M1 = M1 + mass(i)
         xcm = xcm + (mass(i)*x(i))
         ycm = ycm + (mass(i)*y(i))
         zcm = zcm + (mass(i)*z(i))
         M = M + mass(i)
      end do

      do l1 = np1+1, np2
         i = pnum(l1)
         xcm2 = xcm2 + (mass(i)*x(i))
         ycm2 = ycm2 + (mass(i)*y(i))
         zcm2 = zcm2 + (mass(i)*z(i))
         M2 = M2 + mass(i)
         xcm = xcm + (mass(i)*x(i))
         ycm = ycm + (mass(i)*y(i))
         zcm = zcm + (mass(i)*z(i))
         M = M + mass(i)
      end do

      do l1 = np2+1, np3
         i = pnum(l1)
         xcm3 = xcm3 + (mass(i)*x(i))
         ycm3 = ycm3 + (mass(i)*y(i))
         zcm3 = zcm3 + (mass(i)*z(i))
         M3 = M3 + mass(i)
         xcm = xcm + (mass(i)*x(i))
         ycm = ycm + (mass(i)*y(i))
         zcm = zcm + (mass(i)*z(i))
         M = M + mass(i)
      end do

      xcm1 = xcm1/M1
      ycm1 = ycm1/M1
      zcm1 = zcm1/M1

      xcm2 = xcm2/M2
      ycm2 = ycm2/M2
      zcm2 = zcm2/M2

      xcm3 = xcm3/M3
      ycm3 = ycm3/M3
      zcm3 = zcm3/M3

c
c Find center of mass of system
c
c      do l1 = 1, npole3b
c         i = pnum(l1)
c         xcm = xcm + (mass(i)*x(i))
c         ycm = ycm + (mass(i)*y(i))
c         zcm = zcm + (mass(i)*z(i))
c         M = M + mass(i)
c      end do

      xcm = xcm/M
      ycm = ycm/M
      zcm = zcm/M

c
c Get distance between the two molecular
c COMs in reference to the system COM
c
       r_123 = ((xcm1 - xcm)*(xcm1 - xcm)
     &         + (ycm1 - ycm)*(ycm1 - ycm)
     &         + (zcm1 - zcm)*(zcm1 - zcm))
     &       + ((xcm2 - xcm)*(xcm2 - xcm)
     &         + (ycm2 - ycm)*(ycm2 - ycm)
     &         + (zcm2 - zcm)*(zcm2 - zcm))
     &       + ((xcm3 - xcm)*(xcm3 - xcm)
     &         + (ycm3 - ycm)*(ycm3 - ycm)
     &         + (zcm3 - zcm)*(zcm3 - zcm))
      r_123 = r_123/3.0d0
      return
      end


c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine findr3 -- Find max radius between 3 bodies  ##
c     ##                                                         ##
c     #############################################################
c
c     findr3 finds the distance of each particle from the center of
c     mass in each direction. The radius of interaction of three 
c     bodies is defined as the sum of the distances of each particle
c     from the center of mass.
c
      subroutine findr3_minimage (r_123)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'molcul.i'
      include 'combo.i'
      include 'atmtyp.i'
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2
      real*8 xcm3,ycm3,zcm3,xr,yr,zr,xcm_tmp1,ycm_tmp1
      real*8 zcm_tmp1,xr1,yr1,zr1,xr2,yr2,zr2,xr3,yr3,zr3
      real*8 xcm_tmp2,ycm_tmp2,zcm_tmp2
      real*8 xcm,ycm,zcm,r_cm,r2,r_123
      integer i,j,l1,mol1,mol2,M1,M2,M3,M

      xcm1 = 0.0d0
      ycm1 = 0.0d0
      zcm1 = 0.0d0
      xcm2 = 0.0d0
      ycm2 = 0.0d0
      zcm2 = 0.0d0
      xcm3 = 0.0d0
      ycm3 = 0.0d0
      zcm3 = 0.0d0
      xcm = 0.0d0
      ycm = 0.0d0
      zcm = 0.0d0
      M1 = 0.0d0
      M2 = 0.0d0
      M3 = 0.0d0
      M = 0.0d0
c
c Find center of mass of each molecule
c
      do l1 = 1, np1
         i = pnum(l1)
         xcm1 = xcm1 + (mass(i)*x(i))
         ycm1 = ycm1 + (mass(i)*y(i))
         zcm1 = zcm1 + (mass(i)*z(i))
         M1 = M1 + mass(i)
         xcm = xcm + (mass(i)*x(i))
         ycm = ycm + (mass(i)*y(i))
         zcm = zcm + (mass(i)*z(i))
         M = M + mass(i)
      end do

      do l1 = np1+1, np2
         i = pnum(l1)
         xcm2 = xcm2 + (mass(i)*x(i))
         ycm2 = ycm2 + (mass(i)*y(i))
         zcm2 = zcm2 + (mass(i)*z(i))
         M2 = M2 + mass(i)
         xcm = xcm + (mass(i)*x(i))
         ycm = ycm + (mass(i)*y(i))
         zcm = zcm + (mass(i)*z(i))
         M = M + mass(i)
      end do

      do l1 = np2+1, np3
         i = pnum(l1)
         xcm3 = xcm3 + (mass(i)*x(i))
         ycm3 = ycm3 + (mass(i)*y(i))
         zcm3 = zcm3 + (mass(i)*z(i))
         M3 = M3 + mass(i)
         xcm = xcm + (mass(i)*x(i))
         ycm = ycm + (mass(i)*y(i))
         zcm = zcm + (mass(i)*z(i))
         M = M + mass(i)
      end do

      xcm1 = xcm1/M1
      ycm1 = ycm1/M1
      zcm1 = zcm1/M1

      xcm2 = xcm2/M2
      ycm2 = ycm2/M2
      zcm2 = zcm2/M2

      xcm3 = xcm3/M3
      ycm3 = ycm3/M3
      zcm3 = zcm3/M3

      xcm = xcm/M
      ycm = ycm/M
      zcm = zcm/M

      xr = xcm1 - xcm
      yr = ycm1 - ycm
      zr = zcm1 - zcm

      call image(xr,yr,zr)
      
      xcm_tmp1 = xcm1-xr
      ycm_tmp1 = ycm1-yr
      zcm_tmp1 = zcm1-zr

      xr1=xr
      yr1=yr
      zr1=zr

      xr = xcm2 - xcm_tmp1
      yr = ycm2 - ycm_tmp1
      zr = zcm2 - zcm_tmp1

      call image(xr,yr,zr)

      xcm_tmp2 = xcm2 - xr
      ycm_tmp2 = ycm2 - yr
      zcm_tmp2 = zcm2 - zr

      xr2=xr
      yr2=yr
      zr2=zr

      xcm = (xcm_tmp1 + xcm_tmp2)/2
      ycm = (ycm_tmp1 + ycm_tmp2)/2
      zcm = (zcm_tmp1 + zcm_tmp2)/2

      xr = xcm3 - xcm
      yr = ycm3 - ycm
      zr = zcm3 - zcm

      call image(xr,yr,zr)

      xr3=xr
      yr3=yr
      zr3=zr
c
c Get distance between the two molecular
c COMs in reference to the system COM
c
c       r_123 = ((xcm1 - xcm)*(xcm1 - xcm)
c     &         + (ycm1 - ycm)*(ycm1 - ycm)
c     &         + (zcm1 - zcm)*(zcm1 - zcm))
c     &       + ((xcm2 - xcm)*(xcm2 - xcm)
c     &         + (ycm2 - ycm)*(ycm2 - ycm)
c     &         + (zcm2 - zcm)*(zcm2 - zcm))
c     &       + ((xcm3 - xcm)*(xcm3 - xcm)
c     &         + (ycm3 - ycm)*(ycm3 - ycm)
c     &         + (zcm3 - zcm)*(zcm3 - zcm))
      r_123 = (xr1*xr1 + yr1*yr1 + zr1*zr1+
     &        xr2*xr2 + yr2*yr2 + yr2*yr2 +
     &        xr3*xr3 + yr3*yr3 + zr3*zr3)/3.0d0
c      r_123 = r_123/3.0d0
      return
      end


c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine findr3 -- Find max radius between 3 bodies  ##
c     ##                                                         ##
c     #############################################################
c
c     findr3 finds the distance of each particle from the center of
c     mass in each direction. The radius of interaction of three 
c     bodies is defined as the sum of the distances of each particle
c     from the center of mass.
c
      subroutine findr3_minimage2 (r_123)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'molcul.i'
      include 'combo.i'
      include 'atmtyp.i'
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2
      real*8 xcm3,ycm3,zcm3,xr,yr,zr,xcm_tmp1,ycm_tmp1
      real*8 zcm_tmp1,xr1,yr1,zr1,xr2,yr2,zr2,xr3,yr3,zr3
      real*8 xcm_tmp2,ycm_tmp2,zcm_tmp2
      real*8 xcm,ycm,zcm,r_cm,r2,r_123
      integer i,j,l1,mol1,mol2,M1,M2,M3,M

      xcm1 = 0.0d0
      ycm1 = 0.0d0
      zcm1 = 0.0d0
      xcm2 = 0.0d0
      ycm2 = 0.0d0
      zcm2 = 0.0d0
      xcm3 = 0.0d0
      ycm3 = 0.0d0
      zcm3 = 0.0d0
      xcm = 0.0d0
      ycm = 0.0d0
      zcm = 0.0d0
      M1 = 0.0d0
      M2 = 0.0d0
      M3 = 0.0d0
      M = 0.0d0
c
c Find center of mass of each molecule
c
      do l1 = 1, np1
         i = pnum(l1)
         xcm1 = xcm1 + (mass(i)*x(i))
         ycm1 = ycm1 + (mass(i)*y(i))
         zcm1 = zcm1 + (mass(i)*z(i))
         M1 = M1 + mass(i)
c         xcm = xcm + (mass(i)*x(i))
c         ycm = ycm + (mass(i)*y(i))
c         zcm = zcm + (mass(i)*z(i))
c         M = M + mass(i)
      end do

      do l1 = np1+1, np2
         i = pnum(l1)
         xcm2 = xcm2 + (mass(i)*x(i))
         ycm2 = ycm2 + (mass(i)*y(i))
         zcm2 = zcm2 + (mass(i)*z(i))
         M2 = M2 + mass(i)
c         xcm = xcm + (mass(i)*x(i))
c         ycm = ycm + (mass(i)*y(i))
c         zcm = zcm + (mass(i)*z(i))
c         M = M + mass(i)
      end do

      do l1 = np2+1, np3
         i = pnum(l1)
         xcm3 = xcm3 + (mass(i)*x(i))
         ycm3 = ycm3 + (mass(i)*y(i))
         zcm3 = zcm3 + (mass(i)*z(i))
         M3 = M3 + mass(i)
c         xcm = xcm + (mass(i)*x(i))
c         ycm = ycm + (mass(i)*y(i))
c         zcm = zcm + (mass(i)*z(i))
c         M = M + mass(i)
      end do

      xcm1 = xcm1/M1
      ycm1 = ycm1/M1
      zcm1 = zcm1/M1

      xcm2 = xcm2/M2
      ycm2 = ycm2/M2
      zcm2 = zcm2/M2

      xcm3 = xcm3/M3
      ycm3 = ycm3/M3
      zcm3 = zcm3/M3

c      xcm = xcm/M
c      ycm = ycm/M
c      zcm = zcm/M

      xr = xcm1 - xcm2
      yr = ycm1 - ycm2
      zr = zcm1 - zcm2

      call image(xr,yr,zr)

      xr1=xr
      yr1=yr
      zr1=zr

      xr = xcm1 - xcm3
      yr = ycm1 - ycm3
      zr = zcm1 - zcm3

      call image(xr,yr,zr)

      xr2=xr
      yr2=yr
      zr2=zr

c      xcm = (xcm_tmp1 + xcm_tmp2)/2
c      ycm = (ycm_tmp1 + ycm_tmp2)/2
c      zcm = (zcm_tmp1 + zcm_tmp2)/2

      xr = xcm2 - xcm3
      yr = ycm2 - ycm3
      zr = zcm2 - zcm3

      call image(xr,yr,zr)

      xr3=xr
      yr3=yr
      zr3=zr
c
c Get distance between the two molecular
c COMs in reference to the system COM
c
c       r_123 = ((xcm1 - xcm)*(xcm1 - xcm)
c     &         + (ycm1 - ycm)*(ycm1 - ycm)
c     &         + (zcm1 - zcm)*(zcm1 - zcm))
c     &       + ((xcm2 - xcm)*(xcm2 - xcm)
c     &         + (ycm2 - ycm)*(ycm2 - ycm)
c     &         + (zcm2 - zcm)*(zcm2 - zcm))
c     &       + ((xcm3 - xcm)*(xcm3 - xcm)
c     &         + (ycm3 - ycm)*(ycm3 - ycm)
c     &         + (zcm3 - zcm)*(zcm3 - zcm))
      r_123 = (xr1*xr1 + yr1*yr1 + zr1*zr1+
     &        xr2*xr2 + yr2*yr2 + yr2*yr2 +
     &        xr3*xr3 + yr3*yr3 + zr3*zr3)/3.0d0
c      r_123 = r_123/3.0d0
      return
      end

c
      function findr3_minimage2func(pnum,np1,np2,np3)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
c      include 'molcul.i'
c      include 'combo.i'
      include 'atmtyp.i'
      real*8 findr3_minimage2func
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2
      real*8 xcm3,ycm3,zcm3,xr,yr,zr,xcm_tmp1,ycm_tmp1
      real*8 zcm_tmp1,xr1,yr1,zr1,xr2,yr2,zr2,xr3,yr3,zr3
      real*8 xcm_tmp2,ycm_tmp2,zcm_tmp2
      real*8 xcm,ycm,zcm,r_cm,r2,r_123
      integer i,j,l1,mol1,mol2,M1,M2,M3,M
      integer np1,np2,np3 
      integer pnum(*)
      xcm1 = 0.0d0
      ycm1 = 0.0d0
      zcm1 = 0.0d0
      xcm2 = 0.0d0
      ycm2 = 0.0d0
      zcm2 = 0.0d0
      xcm3 = 0.0d0
      ycm3 = 0.0d0
      zcm3 = 0.0d0
      xcm = 0.0d0
      ycm = 0.0d0
      zcm = 0.0d0
      M1 = 0.0d0
      M2 = 0.0d0
      M3 = 0.0d0
      M = 0.0d0
c
c Find center of mass of each molecule
c
      do l1 = 1, np1
         i = pnum(l1)
         xcm1 = xcm1 + (mass(i)*x(i))
         ycm1 = ycm1 + (mass(i)*y(i))
         zcm1 = zcm1 + (mass(i)*z(i))
         M1 = M1 + mass(i)
c         xcm = xcm + (mass(i)*x(i))
c         ycm = ycm + (mass(i)*y(i))
c         zcm = zcm + (mass(i)*z(i))
c         M = M + mass(i)
      end do

      do l1 = np1+1, np2
         i = pnum(l1)
         xcm2 = xcm2 + (mass(i)*x(i))
         ycm2 = ycm2 + (mass(i)*y(i))
         zcm2 = zcm2 + (mass(i)*z(i))
         M2 = M2 + mass(i)
c         xcm = xcm + (mass(i)*x(i))
c         ycm = ycm + (mass(i)*y(i))
c         zcm = zcm + (mass(i)*z(i))
c         M = M + mass(i)
      end do

      do l1 = np2+1, np3
         i = pnum(l1)
         xcm3 = xcm3 + (mass(i)*x(i))
         ycm3 = ycm3 + (mass(i)*y(i))
         zcm3 = zcm3 + (mass(i)*z(i))
         M3 = M3 + mass(i)
c         xcm = xcm + (mass(i)*x(i))
c         ycm = ycm + (mass(i)*y(i))
c         zcm = zcm + (mass(i)*z(i))
c         M = M + mass(i)
      end do
      xcm1 = xcm1/M1
      ycm1 = ycm1/M1
      zcm1 = zcm1/M1

      xcm2 = xcm2/M2
      ycm2 = ycm2/M2
      zcm2 = zcm2/M2

      xcm3 = xcm3/M3
      ycm3 = ycm3/M3
      zcm3 = zcm3/M3

c      xcm = xcm/M
c      ycm = ycm/M
c      zcm = zcm/M

      xr = xcm1 - xcm2
      yr = ycm1 - ycm2
      zr = zcm1 - zcm2

      call image(xr,yr,zr)

      xr1=xr
      yr1=yr
      zr1=zr

      xr = xcm1 - xcm3
      yr = ycm1 - ycm3
      zr = zcm1 - zcm3

      call image(xr,yr,zr)

      xr2=xr
      yr2=yr
      zr2=zr

c      xcm = (xcm_tmp1 + xcm_tmp2)/2
c      ycm = (ycm_tmp1 + ycm_tmp2)/2
c      zcm = (zcm_tmp1 + zcm_tmp2)/2

      xr = xcm2 - xcm3
      yr = ycm2 - ycm3
      zr = zcm2 - zcm3

      call image(xr,yr,zr)

      xr3=xr
      yr3=yr
      zr3=zr
c
c Get distance between the two molecular
c COMs in reference to the system COM
c
c       r_123 = ((xcm1 - xcm)*(xcm1 - xcm)
c     &         + (ycm1 - ycm)*(ycm1 - ycm)
c     &         + (zcm1 - zcm)*(zcm1 - zcm))
c     &       + ((xcm2 - xcm)*(xcm2 - xcm)
c     &         + (ycm2 - ycm)*(ycm2 - ycm)
c     &         + (zcm2 - zcm)*(zcm2 - zcm))
c     &       + ((xcm3 - xcm)*(xcm3 - xcm)
c     &         + (ycm3 - ycm)*(ycm3 - ycm)
c     &         + (zcm3 - zcm)*(zcm3 - zcm))
      r_123 = (xr1*xr1 + yr1*yr1 + zr1*zr1+
     &        xr2*xr2 + yr2*yr2 + yr2*yr2 +
     &        xr3*xr3 + yr3*yr3 + zr3*zr3)/3.0d0
c      r_123 = r_123/3.0d0
c      r_123=findr3_minimage2func()
      findr3_minimage2func=r_123
      return
      end




c     bodies is defined as the sum of the distances of each particle
c     from the center of mass.
c
      subroutine findr3geom (r_123)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'molcul.i'
      include 'combo.i'
      include 'atmtyp.i'
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2
      real*8 xcm3,ycm3,zcm3,xr,yr,zr
      real*8 xcm,ycm,zcm,r_cm,r2,r_123
      integer i,j,l1,mol1,mol2,M1,M2,M3,M


      xcm1 = 0.0d0
      ycm1 = 0.0d0
      zcm1 = 0.0d0
      xcm2 = 0.0d0
      ycm2 = 0.0d0
      zcm2 = 0.0d0
      xcm3 = 0.0d0
      ycm3 = 0.0d0
      zcm3 = 0.0d0
      xcm = 0.0d0
      ycm = 0.0d0
      zcm = 0.0d0
      M1 = 0.0d0
      M2 = 0.0d0
      M3 = 0.0d0
      M = 0.0d0
c
c Find center of mass of each molecule
c
      do l1 = 1, np1
         i = pnum(l1)
         xcm1 = xcm1 + x(i)
         ycm1 = ycm1 + y(i)
         zcm1 = zcm1 + z(i)
         xcm = xcm + x(i)
         ycm = ycm + y(i)
         zcm = zcm + z(i)
      end do

      do l1 = np1+1, np2
         i = pnum(l1)
         xcm2 = xcm2 + x(i)
         ycm2 = ycm2 + y(i)
         zcm2 = zcm2 + z(i)
         xcm = xcm + x(i)
         ycm = ycm + y(i)
         zcm = zcm + z(i)
      end do

      do l1 = np2+1, np3
         i = pnum(l1)
         xcm3 = xcm3 + x(i)
         ycm3 = ycm3 + y(i)
         zcm3 = zcm3 + z(i)
         xcm = xcm + x(i)
         ycm = ycm + y(i)
         zcm = zcm + z(i)
      end do

      xcm1 = xcm1/3.0d0
      ycm1 = ycm1/3.0d0
      zcm1 = zcm1/3.0d0

      xcm2 = xcm2/3.0d0
      ycm2 = ycm2/3.0d0
      zcm2 = zcm2/3.0d0

      xcm3 = xcm3/3.0d0
      ycm3 = ycm3/3.0d0
      zcm3 = zcm3/3.0d0

      xcm = xcm/9.0d0
      ycm = ycm/9.0d0
      zcm = zcm/9.0d0


c
c Find center of mass of system
c
c      do l1 = 1, npole3b
c         i = pnum(l1)
c         xcm = xcm + (mass(i)*x(i))
c         ycm = ycm + (mass(i)*y(i))
c         zcm = zcm + (mass(i)*z(i))
c         M = M + mass(i)
c      end do

c      xcm = xcm/M
c      ycm = ycm/M
c      zcm = zcm/M

c
c Get distance between the two molecular
c COMs in reference to the system COM
c
       r_123 = ((xcm1 - xcm)*(xcm1 - xcm)
     &         + (ycm1 - ycm)*(ycm1 - ycm)
     &         + (zcm1 - zcm)*(zcm1 - zcm))
     &       + ((xcm2 - xcm)*(xcm2 - xcm)
     &         + (ycm2 - ycm)*(ycm2 - ycm)
     &         + (zcm2 - zcm)*(zcm2 - zcm))
     &       + ((xcm3 - xcm)*(xcm3 - xcm)
     &         + (ycm3 - ycm)*(ycm3 - ycm)
     &         + (zcm3 - zcm)*(zcm3 - zcm))
      r_123 = r_123/3.0d0
      return
      end


c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine findr3 -- Find max radius between 3 bodies  ##
c     ##                                                         ##
c     #############################################################
c
c     findr3 finds the distance of each particle from the center of
c     mass in each direction. The radius of interaction of three 
c     bodies is defined as the sum of the distances of each particle
c     from the center of mass.
c
      subroutine findr3test (r_123)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'molcul.i'
      include 'combo.i'
      include 'atmtyp.i'
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2
      real*8 xcm3,ycm3,zcm3,xr,yr,zr
      real*8 xcm,ycm,zcm,r_cm,r2,r_123
      integer i,j,l1,mol1,mol2,M1,M2,M3,M


      xcm1 = 0.0d0
      ycm1 = 0.0d0
      zcm1 = 0.0d0
      xcm2 = 0.0d0
      ycm2 = 0.0d0
      zcm2 = 0.0d0
      xcm3 = 0.0d0
      ycm3 = 0.0d0
      zcm3 = 0.0d0
      xcm = 0.0d0
      ycm = 0.0d0
      zcm = 0.0d0
      M1 = 0.0d0
      M2 = 0.0d0
      M3 = 0.0d0
      M = 0.0d0
c
c Find center of mass of each molecule
c
c      do l1 = 1, np1
c         i = pnum(l1)
c         xcm1 = xcm1 + (mass(i)*x(i))
c         ycm1 = ycm1 + (mass(i)*y(i))
c         zcm1 = zcm1 + (mass(i)*z(i))
c         M1 = M1 + mass(i)
c      end do

c      do l1 = np1+1, np2
c         i = pnum(l1)
c         xcm2 = xcm2 + (mass(i)*x(i))
c         ycm2 = ycm2 + (mass(i)*y(i))
c         zcm2 = zcm2 + (mass(i)*z(i))
c         M2 = M2 + mass(i)
c      end do

c      do l1 = np2+1, np3
c         i = pnum(l1)
c         xcm3 = xcm3 + (mass(i)*x(i))
c         ycm3 = ycm3 + (mass(i)*y(i))
c         zcm3 = zcm3 + (mass(i)*z(i))
c         M3 = M3 + mass(i)
c      end do

c      xcm1 = xcm1/M1
c      ycm1 = ycm1/M1
c      zcm1 = zcm1/M1

c      xcm2 = xcm2/M2
c      ycm2 = ycm2/M2
c      zcm2 = zcm2/M2

c      xcm3 = xcm3/M3
c      ycm3 = ycm3/M3
c      zcm3 = zcm3/M3
      i=pnum(1)
      xcm1 = x(i)
      ycm1 = y(i)
      zcm1 = z(i)

      i=pnum(2)
      xcm2 = x(i)
      ycm2 = y(i)
      zcm2 = z(i)

      i=pnum(3)
      xcm3 = x(i)
      ycm3 = y(i)
      zcm3 = z(i)
         
c
c Find center of mass of system
c
      do l1 = 1, npole3b
         i = pnum(l1)
c         xcm = xcm + (mass(i)*x(i))
c         ycm = ycm + (mass(i)*y(i))
c         zcm = zcm + (mass(i)*z(i))
         xcm = xcm + x(i)
         ycm = ycm + y(i)
         zcm = zcm + z(i)
c         M = M + mass(i)
      end do

      xcm = xcm/3.0d0
      ycm = ycm/3.0d0
      zcm = zcm/3.0d0

c
c Get distance between the two molecular
c COMs in reference to the system COM
c
       r_123 = ((xcm1 - xcm)*(xcm1 - xcm)
     &         + (ycm1 - ycm)*(ycm1 - ycm)
     &         + (zcm1 - zcm)*(zcm1 - zcm))
     &       + ((xcm2 - xcm)*(xcm2 - xcm)
     &         + (ycm2 - ycm)*(ycm2 - ycm)
     &         + (zcm2 - zcm)*(zcm2 - zcm))
     &       + ((xcm3 - xcm)*(xcm3 - xcm)
     &         + (ycm3 - ycm)*(ycm3 - ycm)
     &         + (zcm3 - zcm)*(zcm3 - zcm))
      r_123 = r_123/3.0d0
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
c      subroutine findr2
c      implicit none
c      include 'sizes.i'
c      include 'atoms.i'
c     include 'molcul.i'
c      include 'combo.i'
c      integer i,j
c      real*8 xr,yr,zr,r2
c 
c      moli = 1
c      do moli1 = 1, nmol 
c         do moli2 = moli1+1, nmol
c            i = imol(1,moli1)
c            j = imol(1,moli2)
c            xr = x(i) - x(j)
c            yr = y(i) - y(j)
c            zr = z(i) - z(j)
c            r2 = xr*xr + yr*yr + zr*zr
c            r2b(moli)=sqrt(r2)
c            print*,"liamo:findr2",moli1,moli2,i,j,sqrt(r2)
c            moli = moli + 1
c         end do
c      end do
c      return
c      end

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
c      subroutine findr3
c      implicit none
c      include 'sizes.i'
c      include 'atoms.i'
c      include 'molcul.i'
c      include 'combo.i'
c      integer i,j,k
c      real*8 xr1,yr1,zr1,r2_1
c      real*8 xr2,yr2,zr2,r2_2
c      real*8 xr3,yr3,zr3,r2_3
 
c      moli = 1
c      do moli1 = 1, nmol 
c         do moli2 = moli1+1, nmol
c            do moli3 = moli2+1, nmol
c               i = imol(1,moli1)
c               j = imol(1,moli2)
c               k = imol(1,moli3)
c               xr1 = x(i) - x(j)
c               xr2 = x(i) - x(k)
c               xr3 = x(j) - x(k)
c               yr1 = y(i) - y(j)
c               yr2 = y(i) - y(k)
c               yr3 = y(j) - y(k)
c               zr1 = z(i) - z(j)
c               zr2 = z(i) - z(k)
c               zr3 = z(j) - z(k)
c               r2_1 = xr1*xr1 + yr1*yr1 + zr1*zr1
c               r2_2 = xr2*xr2 + yr2*yr2 + zr2*zr2
c               r2_3 = xr3*xr3 + yr3*yr3 + zr3*zr3
c               r3b_1(moli) = sqrt(r2_1)
c               r3b_2(moli) = sqrt(r2_2)
c               r3b_3(moli) = sqrt(r2_3)
c               moli = moli + 1
c            end do
c         end do
c      end do
c      return
c      end



