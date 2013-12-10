c
c
c     #############################################################
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine empole  --  multipole & polarization energy  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "empole" calculates the electrostatic energy due to
c     atomic multipole and dipole polarizability interactions
c
c
      subroutine empole_3b
      implicit none
      include 'cutoff.i'
      include 'energi.i'
      include 'potent.i'
c
c Liam's variables and inclusions
c
      include 'sizes.i'
      include 'combo.i'
      include 'molcul.i'
      integer i,ii,b1,b2,b3
      real*8 em1(nmol),em2(nmol,nmol),em3(nmol,nmol,nmol)
      real*8 ep1(nmol),ep2(nmol,nmol),ep3(nmol,nmol,nmol)
      return
      end
