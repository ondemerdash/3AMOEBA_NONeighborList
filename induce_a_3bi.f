c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine induce1_3bi  --  double loop multipole analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "induce1_3bi" calculates the direct field from the interaction
c     between atoms in 1-body
c
c
      subroutine induce1_a_3bi (fdir,field,fieldp,dscale,pscale)
      implicit none
      include 'boxes.i'
      include 'sizes.i'
      include 'ewald.i'
      include 'math.i'
      include 'atoms.i'
      include 'mpole.i'
      include 'combo.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'potent.i'
      include 'shunt.i'
      include 'units.i'
      include 'uprior.i'
      include 'molcul.i'
c      real*8 term
c      real*8 ucell(3)
      real*8 fdir(3,*)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8 dscale(*)
      real*8 pscale(*)
      real*8 xr,yr,zr
      real*8 fgrp,r,r2
      real*8 rr3,rr5,rr7
      real*8 ci,dix,diy,diz
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,duir,puir
      real*8 dkr,dukr,pukr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 udsum,upsum
      real*8 eps,epsold
      real*8 epsd,epsp
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      logical proceed
      character*6 mode
      integer i,ii,j,l1,k,m,kk
      integer l2,l3,k1,k2,i1,i2
c
c     get the reciprical space part of the electrostatic field
c
c      call udirect1_3b (field)
c
c     get the real space portion of the electrostatic field
c
c      do l1 = 1, npole3b
c         i = pnum(l1)
c         do j = 1, 3
c            fieldp(j,i) = field(j,i)
c         end do
c      end do
c      call udirect2a_3b (field,fieldp)
c
c     get the self-energy portion of the electrostatic field
c
c      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi


c      do l1 = 1, npole3b
c         i = pnum(l1)
c         do j = 1, 3
c            field(j,i) = field(j,i) + term*rpole(j+1,i)
c            fieldp(j,i) = fieldp(j,i) + term*rpole(j+1,i)
c         end do
c      end do

c      if (boundary .eq. 'VACUUM') then
c         do i = 1, 3
c            ucell(i) = 0.0d0
c         end do
c         do l1 = 1, npole3b
c            i = pnum(l1)
c            ii = ipole(i)
c            ucell(1) = ucell(1) + rpole(2,i) + rpole(1,i)*x(ii)
c            ucell(2) = ucell(2) + rpole(3,i) + rpole(1,i)*y(ii)
c            ucell(3) = ucell(3) + rpole(4,i) + rpole(1,i)*z(ii)
c         end do
c         term = (4.0d0/3.0d0) * pi/volbox
c         do l1 = 1, npole3b
c            i = pnum(l1)
c            do j = 1, 3
c               field(j,i) = field(j,i) - term*ucell(j)
c               fieldp(j,i) = fieldp(j,i) - term*ucell(j)
c            end do
c         end do
      return
      end


c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine induce2_3bi  --  double loop multipole analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "induce2_3bi" calculates the direct field from the interaction
c     between atoms in 1-body
c
c
      subroutine induce2_a_3bi (fdir,field,fieldp,dscale,pscale)
      implicit none
      include 'boxes.i'
      include 'sizes.i'
      include 'ewald.i'
      include 'math.i'
      include 'atoms.i'
      include 'mpole.i'
      include 'combo.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'potent.i'
      include 'shunt.i'
      include 'units.i'
      include 'uprior.i'
      include 'molcul.i'
c      real*8 term
c      real*8 ucell(3)
      real*8 fdir(3,*)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      real*8 dscale(*)
      real*8 pscale(*)
      real*8 xr,yr,zr
      real*8 fgrp,r,r2
      real*8 rr3,rr5,rr7
      real*8 ci,dix,diy,diz
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,duir,puir
      real*8 dkr,dukr,pukr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 udsum,upsum
      real*8 eps,epsold
      real*8 epsd,epsp
      real*8 damp,expdamp
      real*8 scale3,scale5
      real*8 scale7
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      character*6 mode
      integer i,ii,j,l1,k,m,kk
      integer l2,l3,k1,k2,i1,i2
      logical proceed
c
c     get the reciprical space part of the electrostatic field
c
c      call udirect1_3b (field)
c
c     get the real space portion of the electrostatic field
c
c      do l1 = 1, npole3b
c         i = pnum(l1)
c         do j = 1, 3
c            fieldp(j,i) = field(j,i)
c         end do
c      end do
c      call udirect2a_3b (field,fieldp)
c      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
c      do l1 = 1, npole3b
c         i = pnum(l1)
c         do j = 1, 3
c            field(j,i) = field(j,i) + term*rpole(j+1,i)
c            fieldp(j,i) = fieldp(j,i) + term*rpole(j+1,i)
c         end do
c      end do

c      if (boundary .eq. 'VACUUM') then
c         do i = 1, 3
c            ucell(i) = 0.0d0
c         end do
c         do l1 = 1, npole3b
c            i = pnum(l1)
c            ii = ipole(i)
c            ucell(1) = ucell(1) + rpole(2,i) + rpole(1,i)*x(ii)
c            ucell(2) = ucell(2) + rpole(3,i) + rpole(1,i)*y(ii)
c            ucell(3) = ucell(3) + rpole(4,i) + rpole(1,i)*z(ii)
c         end do
c         term = (4.0d0/3.0d0) * pi/volbox
c         do l1 = 1, npole3b
c            i = pnum(l1)
c            do j = 1, 3
c               field(j,i) = field(j,i) - term*ucell(j)
c               fieldp(j,i) = fieldp(j,i) - term*ucell(j)
c            end do
c         end do
c      end if

      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine induce3_3bi  --  double loop multipole analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "induce3_3bi" calculates the direct field from the interaction
c     between atoms in 1-body
c
c
      subroutine induce3_a_3bi (fdir,field,fieldp)
      return
      end

