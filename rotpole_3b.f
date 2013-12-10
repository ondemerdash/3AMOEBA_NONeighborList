c
c
c     ############################################################
c     ##  COPYRIGHT (C) 1995 by Yong Kong & Jay William Ponder  ##
c     ##                  All Rights Reserved                   ##
c     ############################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine rotpole_3b  --  rotate multipoles to global frame  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "rotpole" constructs the set of atomic multipoles in the global
c     frame by applying the correct rotation matrix for each site
c
c
      subroutine rotpole_3b
      implicit none
      include 'sizes.i'
      include 'mpole.i'
      include 'molcul.i'
      include 'combo.i'
      integer l1
      integer i
      real*8 a(3,3)
c
c
c     rotate the atomic multipoles at each site in turn
c
      do l1 = 1, npole3b
         i = pnum(l1)
         call rotmat (i,a)
         call rotsite (i,a)
      end do
      return
      end
