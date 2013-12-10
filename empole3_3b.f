
c     ##  COPYRIGHT (C) 1999 by Pengyu Ren & Jay William Ponder  ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine empole3_3b  --  mpole/polar energy & analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "empole3" calculates the electrostatic energy due to
c     atomic multipole and dipole polarizability interactions,
c     and partitions the energy among the atoms
c
c
      subroutine empole3_3b
      implicit none
      include 'sizes.i'
      include 'analyz.i'
      include 'atoms.i'
      include 'cutoff.i'
      include 'energi.i'
      include 'mpole.i'
      include 'potent.i'
      include 'molcul.i'
      include 'combo.i'
      integer i,ii,b1,b2,b3
      real*8 em1(nmol),em2(nmol,nmol),em3(nmol,nmol,nmol)
      real*8 ep1(nmol),ep2(nmol,nmol),ep3(nmol,nmol,nmol)
      real*8 emliam,epliam
      emliam = 0
      epliam = 0
      delta1 = 0
      delta2 = 0
      delta3 = 0
      b1 = 1
      b2 = 2
      b3 = 3
c
c     choose the method for summing over multipole interactions
c
      if (use_ewald) then
         if (use_mlist) then
            call empole3d_3bi
         else
            call empole3c_3bi
            print*,"liam",ep,em
            call empole3c
            print*,"regular tinker",ep,em
        end if
      else
        if (use_mlist) then
            call empole3b_3bi
         else
            call empole3a_3bi
            print*,"liam",ep,em
            call empole3a
           print*,"regular tinker",ep,em
         end if
      end if
c
c     zero out energy terms and analysis which are not in use
c
      if (.not. use_mpole) then
         em = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            aem(ii) = 0.0d0
         end do
      end if
      if (.not. use_polar) then
         ep = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            aep(ii) = 0.0d0
         end do
      end if

      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine empole3a_3b  --  double loop multipole analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "empole3a" calculates the atomic multipole and dipole
c     polarizability interaction energy using a double loop,
c     and partitions the energy among the atoms
c
c
      subroutine empole3a_3b
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
      include 'math.i'
      include 'molcul.i'
      include 'mplpot.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'potent.i'
      include 'shunt.i'
      include 'usage.i'
      include 'combo.i'
      integer l1,l2,l3
      integer i,j,k
      integer ii,kk
      integer ix,iy,iz
      integer kx,ky,kz
      real*8 e,ei,fgrp
      real*8 f,fm,fp
      real*8 r,r2,xr,yr,zr
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 ukx,uky,ukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qix,qiy,qiz
      real*8 qkx,qky,qkz
      real*8 scale3,scale5
      real*8 scale7
      real*8 sc(10),sci(8)
      real*8 gl(0:4),gli(3)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: pscale(:)
      logical proceed
      logical header,huge
      logical usei,usek
      logical muse,puse
      character*6 mode
c
c
c     zero out multipole and polarization energy and partitioning
c
      nem = 0
      nep = 0
      em = 0.0d0
      ep = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
         aep(i) = 0.0d0
      end do
      header = .true.
c
c     check the sign of multipole components at chiral sites
c
      call chkpole_3b
c
c     rotate the multipole components into the global frame
c
      call rotpole_3b
c
c     compute the induced dipoles at each polarizable atom
c
      call induce_3b
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (pscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      if (npole3b .eq. 0)  return
      do i = 1, n
         mscale(i) = 1.0d0
         pscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     calculate the multipole interaction energy term
c
      do l1 = 1, npole3b-1
         i = pnum(l1)
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
c
c     decide whether to compute the current interaction
c
         do l3 = l1+1, npole3b
            k = pnum(l3)
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
                  ukx = uind(1,k)
                  uky = uind(2,k)
                  ukz = uind(3,k)
c
c     construct some intermediate quadrupole values
c
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
                  sc(2) = dix*dkx + diy*dky + diz*dkz
                  sc(3) = dix*xr + diy*yr + diz*zr
                  sc(4) = dkx*xr + dky*yr + dkz*zr
                  sc(5) = qix*xr + qiy*yr + qiz*zr
                  sc(6) = qkx*xr + qky*yr + qkz*zr
                  sc(7) = qix*dkx + qiy*dky + qiz*dkz
                  sc(8) = qkx*dix + qky*diy + qkz*diz
                  sc(9) = qix*qkx + qiy*qky + qiz*qkz
                  sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                        + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
                  sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                        + diy*uky + uiz*dkz + diz*ukz
                  sci(3) = uix*xr + uiy*yr + uiz*zr
                  sci(4) = ukx*xr + uky*yr + ukz*zr
                  sci(7) = qix*ukx + qiy*uky + qiz*ukz
                  sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
                  gl(0) = ci*ck
                  gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
                  gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
     &                       + 2.0d0*(sc(7)-sc(8)+sc(10))
                  gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
                  gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
                  gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
                  gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                        - sc(3)*sci(4)
                  gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
                  rr1 = 1.0d0 / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = 1.0d0 - expdamp
                        scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                        scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                          *expdamp
                     end if
                  end if
                  e = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
     &                   + gl(3)*rr7 + gl(4)*rr9
                  ei = gli(1)*rr3*scale3 + gli(2)*rr5*scale5
     &                    + gli(3)*rr7*scale7
c
c     apply the energy adjustments for scaled interactions
c
                  fm = f * mscale(kk)
                  fp = f * pscale(kk)
                  e = fm * e
                  ei = 0.5d0 * fp * ei
c
c     scale the interaction based on its group membership;
c     polarization cannot be group scaled as it is not pairwise
c
                  if (use_group) then
                     e = e * fgrp
c                    ei = ei * fgrp
                  end if
c
c     increment the overall multipole and polarization energies
c
                  muse = (use_mpole .and. mscale(kk).ne.0.0d0)
                  puse = (use_polar .and. pscale(kk).ne.0.0d0)
                  if (muse)  nem = nem + 1
                  if (puse)  nep = nep + 1
                  em = em + e
                  ep = ep + ei
                  aem(ii) = aem(ii) + 0.5d0*e
                  aem(kk) = aem(kk) + 0.5d0*e
                  aep(ii) = aep(ii) + 0.5d0*ei
                  aep(kk) = aep(kk) + 0.5d0*ei
c
c     increment the total intermolecular energy
c
                  if (molcule(ii) .ne. molcule(kk)) then
                     einter = einter + e + ei
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (max(abs(e),abs(ei)) .gt. 100.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (muse .or. puse) then
                        if (header) then
                           header = .false.
                           write (iout,10)
   10                      format (/,' Individual Multipole and',
     &                                ' Polarization Interactions :',
     &                             //,' Type',14x,'Atom Names',
     &                                15x,'Distance',6x,'Energies',
     &                                ' (MPol,Polar)',/)
                        end if
                        write (iout,20)  ii,name(ii),kk,name(kk),r,e,ei
   20                   format (' M-Pole',4x,2(i7,'-',a3),9x,
     &                             f10.4,2x,2f12.4)
                     end if
                  end if
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      do l1 = 1, npole3b
         i = pnum(l1)
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(k)
         pdi = pdamp(i)
         pti = thole(i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
c
c     decide whether to compute the current interaction
c
         do l3 = l1, npole3b
            k = pnum(l3)
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               do j = 1, ncell
                  xr = x(kk) - x(ii)
                  yr = y(kk) - y(ii)
                  zr = z(kk) - z(ii)
                  call imager (xr,yr,zr,j)
                  r2 = xr*xr + yr* yr + zr*zr
                  if (r2 .le. off2) then
                     r = sqrt(r2)
                     ck = rpole(1,k)
                     dkx = rpole(2,k)
                     dky = rpole(3,k)
                     dkz = rpole(4,k)
                     qkxx = rpole(5,k)
                     qkxy = rpole(6,k)
                     qkxz = rpole(7,k)
                     qkyy = rpole(9,k)
                     qkyz = rpole(10,k)
                     qkzz = rpole(13,k)
                     ukx = uind(1,k)
                     uky = uind(2,k)
                     ukz = uind(3,k)
c
c     construct some intermediate quadrupole values
c
                     qix = qixx*xr + qixy*yr + qixz*zr
                     qiy = qixy*xr + qiyy*yr + qiyz*zr
                     qiz = qixz*xr + qiyz*yr + qizz*zr
                     qkx = qkxx*xr + qkxy*yr + qkxz*zr
                     qky = qkxy*xr + qkyy*yr + qkyz*zr
                     qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
                     sc(2) = dix*dkx + diy*dky + diz*dkz
                     sc(3) = dix*xr + diy*yr + diz*zr
                     sc(4) = dkx*xr + dky*yr + dkz*zr
                     sc(5) = qix*xr + qiy*yr + qiz*zr
                     sc(6) = qkx*xr + qky*yr + qkz*zr
                     sc(7) = qix*dkx + qiy*dky + qiz*dkz
                     sc(8) = qkx*dix + qky*diy + qkz*diz
                     sc(9) = qix*qkx + qiy*qky + qiz*qkz
                     sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                           + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
                     sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                           + diy*uky + uiz*dkz + diz*ukz
                     sci(3) = uix*xr + uiy*yr + uiz*zr
                     sci(4) = ukx*xr + uky*yr + ukz*zr
                     sci(7) = qix*ukx + qiy*uky + qiz*ukz
                     sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
                     gl(0) = ci*ck
                     gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
                     gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
     &                          + 2.0d0*(sc(7)-sc(8)+sc(10))
                     gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
                     gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
                     gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
                     gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                           - sc(3)*sci(4)
                     gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
                     rr1 = 1.0d0 / r
                     rr3 = rr1 / r2
                     rr5 = 3.0d0 * rr3 / r2
                     rr7 = 5.0d0 * rr5 / r2
                     rr9 = 7.0d0 * rr7 / r2
                     scale3 = 1.0d0
                     scale5 = 1.0d0
                     scale7 = 1.0d0
                     damp = pdi * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(pti,thole(k))
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           expdamp = exp(damp)
                           scale3 = 1.0d0 - expdamp
                           scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                           scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                             *expdamp
                        end if
                     end if
                     e = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
     &                      + gl(3)*rr7 + gl(4)*rr9
                     ei = gli(1)*rr3*scale3 + gli(2)*rr5*scale5
     &                       + gli(3)*rr7*scale7
c
c     apply the energy adjustments for scaled interactions
c
                     fm = f
                     fp = f
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           fm = fm * mscale(kk)
                           fp = fp * pscale(kk)
                        end if
                     end if
                     e = fm * e
                     ei = 0.5d0 * fp * ei
c
c     scale the interaction based on its group membership;
c     polarization cannot be group scaled as it is not pairwise
c
                     if (use_group) then
                        e = e * fgrp
c                       ei = ei * fgrp
                     end if
c
c     increment the overall multipole and polarization energies
c
                     if (ii .eq. kk) then
                        e = 0.5d0 * e
                        ei = 0.5d0 * ei
                     end if
                     nem = nem + 1
                     nep = nep + 1
                     em = em + e
                     ep = ep + ei
                     aem(ii) = aem(ii) + 0.5d0*e
                     aem(kk) = aem(kk) + 0.5d0*e
                     aep(ii) = aep(ii) + 0.5d0*ei
                     aep(kk) = aep(kk) + 0.5d0*ei
c
c     increment the total intermolecular energy
c
                     einter = einter + e + ei
c
c     print a message if the energy of this interaction is large
c
                     huge = (max(abs(e),abs(ei)) .gt. 100.0d0)
                     if (debug .or. (verbose.and.huge)) then
                        if (header) then
                           header = .false.
                           write (iout,30)
   30                      format (/,' Individual Multipole and',
     &                                ' Polarization Interactions :',
     &                             //,' Type',14x,'Atom Names',
     &                                15x,'Distance',6x,'Energies',
     &                                ' (MPol,Polar)',/)
                        end if
                        write (iout,40)  ii,name(ii),kk,name(kk),r,e,ei
   40                   format (' M-Pole',4x,2(i7,'-',a3),1x,
     &                             '(X)',5x,f10.4,2x,2f12.4)
                     end if
                  end if
               end do
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c

      deallocate (mscale)
      deallocate (pscale)
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine empole3b_3b  --  neighbor list multipole analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "empole3b" calculates the atomic multipole and dipole
c     polarizability interaction energy using a neighbor list,
c     and partitions the energy among the atoms
c
c
      subroutine empole3b_3b
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'group.i'
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
      include 'math.i'
      include 'molcul.i'
      include 'mplpot.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'potent.i'
      include 'shunt.i'
      include 'usage.i'
      integer i,j,k
      integer ii,kk,kkk
      integer ix,iy,iz
      integer kx,ky,kz
      real*8 e,ei,fgrp
      real*8 f,fm,fp
      real*8 r,r2,xr,yr,zr
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 ukx,uky,ukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qix,qiy,qiz
      real*8 qkx,qky,qkz
      real*8 scale3,scale5
      real*8 scale7
      real*8 sc(10),sci(8)
      real*8 gl(0:4),gli(3)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: pscale(:)
      logical proceed
      logical header,huge
      logical usei,usek
      logical muse,puse
      character*6 mode
c
c
c     zero out multipole and polarization energy and partitioning
c
      nem = 0
      nep = 0
      em = 0.0d0
      ep = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
         aep(i) = 0.0d0
      end do
      header = .true.
c
c     check the sign of multipole components at chiral sites
c
      call chkpole_3b
c
c     rotate the multipole components into the global frame
c
      call rotpole_3b
c
c     compute the induced dipoles at each polarizable atom
c
      call induce_3b
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (pscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      if (npole .eq. 0)  return
      do i = 1, n
         mscale(i) = 1.0d0
         pscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'MPOLE'
      call switch (mode)
c
c     calculate the multipole interaction energy term
c
      do i = 1, npole-1
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
c
c     decide whether to compute the current interaction
c
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
c
c     compute the energy contribution for this interaction
c
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
                  ukx = uind(1,k)
                  uky = uind(2,k)
                  ukz = uind(3,k)
c
c     construct some intermediate quadrupole values
c
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
                  sc(2) = dix*dkx + diy*dky + diz*dkz
                  sc(3) = dix*xr + diy*yr + diz*zr
                  sc(4) = dkx*xr + dky*yr + dkz*zr
                  sc(5) = qix*xr + qiy*yr + qiz*zr
                  sc(6) = qkx*xr + qky*yr + qkz*zr
                  sc(7) = qix*dkx + qiy*dky + qiz*dkz
                  sc(8) = qkx*dix + qky*diy + qkz*diz
                  sc(9) = qix*qkx + qiy*qky + qiz*qkz
                  sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                        + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
                  sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                        + diy*uky + uiz*dkz + diz*ukz
                  sci(3) = uix*xr + uiy*yr + uiz*zr
                  sci(4) = ukx*xr + uky*yr + ukz*zr
                  sci(7) = qix*ukx + qiy*uky + qiz*ukz
                  sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
                  gl(0) = ci*ck
                  gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
                  gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
     &                       + 2.0d0*(sc(7)-sc(8)+sc(10))
                  gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
                  gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
                  gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
                  gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                        - sc(3)*sci(4)
                  gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
                  rr1 = 1.0d0 / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = 1.0d0 - expdamp
                        scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                        scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                          *expdamp
                     end if
                  end if
                  e = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
     &                   + gl(3)*rr7 + gl(4)*rr9
                  ei = gli(1)*rr3*scale3 + gli(2)*rr5*scale5
     &                    + gli(3)*rr7*scale7
c
c     apply the energy adjustments for scaled interactions
c
                  fm = f * mscale(kk)
                  fp = f * pscale(kk)
                  e = fm * e
                  ei = 0.5d0 * fp * ei
c
c     scale the interaction based on its group membership;
c     polarization cannot be group scaled as it is not pairwise
c
                  if (use_group) then
                     e = e * fgrp
c                    ei = ei * fgrp
                  end if
c
c     increment the overall multipole and polarization energies
c
                  muse = (use_mpole .and. mscale(kk).ne.0.0d0)
                  puse = (use_polar .and. pscale(kk).ne.0.0d0)
                  if (muse)  nem = nem + 1
                  if (puse)  nep = nep + 1
                  em = em + e
                  ep = ep + ei
                  aem(ii) = aem(ii) + 0.5d0*e
                  aem(kk) = aem(kk) + 0.5d0*e
                  aep(ii) = aep(ii) + 0.5d0*ei
                  aep(kk) = aep(kk) + 0.5d0*ei
c
c     increment the total intermolecular energy
c
                  if (molcule(ii) .ne. molcule(kk)) then
                     einter = einter + e + ei
                  end if
c
c     print a message if the energy of this interaction is large
c
                  huge = (max(abs(e),abs(ei)) .gt. 100.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (muse .or. puse) then
                        if (header) then
                           header = .false.
                           write (iout,10)
   10                      format (/,' Individual Multipole and',
     &                                ' Polarization Interactions :',
     &                             //,' Type',14x,'Atom Names',
     &                                15x,'Distance',6x,'Energies',
     &                                ' (MPol,Polar)',/)
                        end if
                        write (iout,20)  ii,name(ii),kk,name(kk),r,e,ei
   20                   format (' M-Pole',4x,2(i7,'-',a3),9x,
     &                             f10.4,2x,2f12.4)
                     end if
                  end if
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (pscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine empole3c_3b  --  Ewald multipole analysis via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "empole3c" calculates the atomic multipole and dipole
c     polarizability interaction energy using a particle mesh
c     Ewald summation and double loop, and partitions the energy
c     among the atoms
c
c
      subroutine empole3c_3b
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'energi.i'
      include 'ewald.i'
      include 'inter.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'molcul.i'
      include 'combo.i'
      integer l1
      integer i,ii
      real*8 e,ei,eintra
      real*8 f,term,fterm
      real*8 cii,dii,qii,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
c
c
c     zero out the multipole and polarization energies
c
      nem = 0
      nep = 0
      em = 0.0d0
      ep = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
         aep(i) = 0.0d0
      end do
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
c
c     compute the induced dipoles at each polarizable atom
c
      call induce_3b
c
c     compute the reciprocal space part of the Ewald summation
c
      call erecip3_3b
c
c     compute the real space part of the Ewald summation
c
      call ereal3c_3b (eintra)
c
c     compute the self-energy part of the Ewald summation
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do l1 = 1, npole3b
         i = pnum(l1)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = qixx*qixx + qiyy*qiyy + qizz*qizz
     &            + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
         uii = dix*uix + diy*uiy + diz*uiz
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         ei = fterm * term * uii / 3.0d0
         nem = nem + 1
         nep = nep + 1
         em = em + e
         ep = ep + ei
         aem(i) = aem(i) + e
         aep(i) = aep(i) + ei
      end do

c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xu = 0.0d0
         yu = 0.0d0
         zu = 0.0d0
	 do l1 = 1, npole3b
            i = pnum(l1)
            ii = ipole(i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            uix = uind(1,i)
            uiy = uind(2,i)
            uiz = uind(3,i)
            xd = xd + dix + rpole(1,i)*x(ii)
            yd = yd + diy + rpole(1,i)*y(ii)
            zd = zd + diz + rpole(1,i)*z(ii)
            xu = xu + uix
            yu = yu + uiy
            zu = zu + uiz
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         nem = nem + 1
         nep = nep + 1
         em = em + term*(xd*xd+yd*yd+zd*zd)
         ep = ep + term*(xd*xu+yd*yu+zd*zu)
      end if
c
c     intermolecular energy is total minus intramolecular part
c
      einter = einter + em + ep - eintra
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ereal3c_3b  --  real space mpole analysis via loop  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ereal3c" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipole interactions and dipole
c     polarizability and partitions the energy among the atoms
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine ereal3c_3b (eintra)
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'ewald.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'molcul.i'
      include 'mplpot.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'potent.i'
      include 'shunt.i'
      include 'combo.i'
      integer l1,l3
      integer i,j,k,m
      integer ii,kk
      real*8 e,ei,eintra
      real*8 f,erfc
      real*8 r,r2,xr,yr,zr
      real*8 bfac,exp2a
      real*8 efix,eifix
      real*8 ralpha
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 alsq2,alsq2n
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 ukx,uky,ukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qix,qiy,qiz
      real*8 qkx,qky,qkz
      real*8 scale3,scale5
      real*8 scale7
      real*8 sc(10),sci(8)
      real*8 gl(0:4),gli(3)
      real*8 bn(0:4)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: pscale(:)
      logical header,huge
      logical muse,puse
      character*6 mode
      external erfc
c
c
c     zero out the intramolecular portion of the Ewald energy
c
      eintra = 0.0d0
      if (npole .eq. 0)  return
      header = .true.
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (pscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         mscale(i) = 1.0d0
         pscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     compute the real space portion of the Ewald summation
c
      do l1 = 1, npole3b-1
         i = pnum(l1)
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
         do l3 = l1+1, npole3b
            k = pnum(l3)
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
               ukx = uind(1,k)
               uky = uind(2,k)
               ukz = uind(3,k)
c
c     calculate the error function damping terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)
     &            alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do m = 1, 4
                  bfac = dble(m+m-1)
                  alsq2n = alsq2 * alsq2n
                  bn(m) = (bfac*bn(m-1)+alsq2n*exp2a) / r2
               end do
c
c     construct some intermediate quadrupole values
c
               qix = qixx*xr + qixy*yr + qixz*zr
               qiy = qixy*xr + qiyy*yr + qiyz*zr
               qiz = qixz*xr + qiyz*yr + qizz*zr
               qkx = qkxx*xr + qkxy*yr + qkxz*zr
               qky = qkxy*xr + qkyy*yr + qkyz*zr
               qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
               sc(2) = dix*dkx + diy*dky + diz*dkz
               sc(3) = dix*xr + diy*yr + diz*zr
               sc(4) = dkx*xr + dky*yr + dkz*zr
               sc(5) = qix*xr + qiy*yr + qiz*zr
               sc(6) = qkx*xr + qky*yr + qkz*zr
               sc(7) = qix*dkx + qiy*dky + qiz*dkz
               sc(8) = qkx*dix + qky*diy + qkz*diz
               sc(9) = qix*qkx + qiy*qky + qiz*qkz
               sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                     + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
               sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                     + diy*uky + uiz*dkz + diz*ukz
               sci(3) = uix*xr + uiy*yr + uiz*zr
               sci(4) = ukx*xr + uky*yr + ukz*zr
               sci(7) = qix*ukx + qiy*uky + qiz*ukz
               sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
               gl(0) = ci*ck
               gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
               gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
     &                    + 2.0d0*(sc(7)-sc(8)+sc(10))
               gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
               gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
               gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
               gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                     - sc(3)*sci(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
               e = gl(0)*bn(0) + gl(1)*bn(1) + gl(2)*bn(2)
     &                + gl(3)*bn(3) + gl(4)*bn(4)
               ei = gli(1)*bn(1) + gli(2)*bn(2) + gli(3)*bn(3)
c
c     full real space energies needed for scaled interactions
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               scale3 = pscale(kk)
               scale5 = pscale(kk)
               scale7 = pscale(kk)
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = scale3 * (1.0d0-expdamp)
                     scale5 = scale5 * (1.0d0-(1.0d0-damp)*expdamp)
                     scale7 = scale7 * (1.0d0-(1.0d0-damp
     &                                    +0.6d0*damp**2)*expdamp)
                  end if
               end if
               efix = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
     &                   + gl(3)*rr7 + gl(4)*rr9
               eifix = gli(1)*rr3*(1.0d0-scale3)
     &                    + gli(2)*rr5*(1.0d0-scale5)
     &                    + gli(3)*rr7*(1.0d0-scale7)
c
c     apply the energy adjustments for scaled interactions
c
               e = e - efix*(1.0d0-mscale(kk))
               ei = ei - eifix
               e = f * e
               ei = 0.5d0 * f * ei
c
c     increment the overall multipole and polarization energies
c
               muse = use_mpole
               puse = use_polar
               if (muse)  nem = nem + 1
               if (puse)  nep = nep + 1
               em = em + e
               ep = ep + ei
               aem(ii) = aem(ii) + 0.5d0*e
               aem(kk) = aem(kk) + 0.5d0*e
               aep(ii) = aep(ii) + 0.5d0*ei
               aep(kk) = aep(kk) + 0.5d0*ei
c
c     increment the total intramolecular energy
c
               efix = f * efix * mscale(kk)
               eifix = gli(1)*rr3*scale3 + gli(2)*rr5*scale5
     &                    + gli(3)*rr7*scale7
               eifix = 0.5d0 * f * eifix
               if (molcule(ii) .eq. molcule(kk)) then
                  eintra = eintra + efix + eifix
               end if
c
c     print a message if the energy of this interaction is large
c
               huge = (max(abs(efix),abs(eifix)) .gt. 100.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (muse .or. puse) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Real Space Multipole and',
     &                             ' Polarization Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',6x,'Energies',
     &                             ' (MPol,Polar)',/)
                     end if
                     write (iout,20)  ii,name(ii),kk,name(kk),r,
     &                                efix,eifix
   20                format (' M-Pole',4x,2(i7,'-',a3),9x,
     &                          f10.4,2x,2f12.4)
                  end if
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (.not. use_replica)  return
c
c     calculate interaction energy with other unit cells
c
      print *, "e = ", e
      do l1 = 1, npole3b
         i = pnum(l1)
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
         do l3 = 1, npole3b
            k = pnum(l3)
            kk = ipole(k)
            do j = 1, ncell
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               call imager (xr,yr,zr,j)
               r2 = xr*xr + yr* yr + zr*zr
               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
                  ukx = uind(1,k)
                  uky = uind(2,k)
                  ukz = uind(3,k)
c
c     calculate the error function damping terms
c
                  ralpha = aewald * r
                  bn(0) = erfc(ralpha) / r
                  alsq2 = 2.0d0 * aewald**2
                  alsq2n = 0.0d0
                  if (aewald .gt. 0.0d0)
     &               alsq2n = 1.0d0 / (sqrtpi*aewald)
                  exp2a = exp(-ralpha**2)
                  do m = 1, 4
                     bfac = dble(m+m-1)
                     alsq2n = alsq2 * alsq2n
                     bn(m) = (bfac*bn(m-1)+alsq2n*exp2a) / r2
                  end do
c
c     construct some intermediate quadrupole values
c
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
                  sc(2) = dix*dkx + diy*dky + diz*dkz
                  sc(3) = dix*xr + diy*yr + diz*zr
                  sc(4) = dkx*xr + dky*yr + dkz*zr
                  sc(5) = qix*xr + qiy*yr + qiz*zr
                  sc(6) = qkx*xr + qky*yr + qkz*zr
                  sc(7) = qix*dkx + qiy*dky + qiz*dkz
                  sc(8) = qkx*dix + qky*diy + qkz*diz
                  sc(9) = qix*qkx + qiy*qky + qiz*qkz
                  sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                        + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
                  sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                        + diy*uky + uiz*dkz + diz*ukz
                  sci(3) = uix*xr + uiy*yr + uiz*zr
                  sci(4) = ukx*xr + uky*yr + ukz*zr
                  sci(7) = qix*ukx + qiy*uky + qiz*ukz
                  sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
                  gl(0) = ci*ck
                  gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
                  gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
     &                       + 2.0d0*(sc(7)-sc(8)+sc(10))
                  gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
                  gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
                  gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
                  gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                        - sc(3)*sci(4)
                  gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
                  e = gl(0)*bn(0) + gl(1)*bn(1) + gl(2)*bn(2)
     &                   + gl(3)*bn(3) + gl(4)*bn(4)
                  ei = gli(1)*bn(1) + gli(2)*bn(2) + gli(3)*bn(3)
c
c     full real space energies needed for scaled interactions
c
                  rr1 = 1.0d0 / r
                  rr3 = rr1 / r2
                  rr5 = 3.0d0 * rr3 / r2
                  rr7 = 5.0d0 * rr5 / r2
                  rr9 = 7.0d0 * rr7 / r2
                  scale3 = 1.0d0
                  scale5 = 1.0d0
                  scale7 = 1.0d0
                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3 = 1.0d0 - expdamp
                        scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                        scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                          *expdamp
                        if (use_polymer .and. r2.le.polycut2) then
                           scale3 = scale3 * pscale(kk)
                           scale5 = scale5 * pscale(kk)
                           scale7 = scale7 * pscale(kk)
                        end if
                     end if
                  end if
                  efix = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
     &                      + gl(3)*rr7 + gl(4)*rr9
                  eifix = gli(1)*rr3*(1.0d0-scale3)
     &                       + gli(2)*rr5*(1.0d0-scale5)
     &                       + gli(3)*rr7*(1.0d0-scale7)
c
c     apply the energy adjustments for scaled interactions
c
                  if (use_polymer .and. r2.le.polycut2)
     &               e = e - efix*(1.0d0-mscale(kk))
                  ei = ei - eifix
c
c     increment the overall multipole and polarization energies
c
                  e = f * e
                  ei = 0.5d0 * f * ei
                  if (ii .eq. kk) then
                     e = 0.5d0 * e
                     ei = 0.5d0 * ei
                  end if
                  nem = nem + 1
                  nep = nep + 1
                  em = em + e
                  ep = ep + ei
                  aem(ii) = aem(ii) + 0.5d0*e
                  aem(kk) = aem(kk) + 0.5d0*e
                  aep(ii) = aep(ii) + 0.5d0*ei
                  aep(kk) = aep(kk) + 0.5d0*ei
c
c     print a message if the energy of this interaction is large
c
                  efix = f * efix * mscale(kk)
                  eifix = gli(1)*rr3*scale3 + gli(2)*rr5*scale5
     &                       + gli(3)*rr7*scale7
                  eifix = 0.5d0 * f * eifix
                  huge = (max(abs(efix),abs(eifix)) .gt. 100.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,30)
   30                   format (/,' Real Space Multipole and',
     &                             ' Polarization Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',6x,'Energies',
     &                             ' (MPol,Polar)',/)
                     end if
                     write (iout,40)  ii,name(ii),kk,name(kk),r,
     &                                efix,eifix
   40                format (' M-Pole',4x,2(i7,'-',a3),1x,
     &                          '(X)',5x,f10.4,2x,2f12.4)
                  end if
               end if
            end do
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (pscale)
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine empole3d_3b  --  Ewald multipole analysis via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "empole3d" calculates the atomic multipole and dipole
c     polarizability interaction energy using a particle mesh
c     Ewald summation and a neighbor list, and partitions the
c     energy among the atoms
c
c
      subroutine empole3d_3b
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'energi.i'
      include 'ewald.i'
      include 'inter.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      integer i,ii
      real*8 e,ei,eintra
      real*8 f,term,fterm
      real*8 cii,dii,qii,uii
      real*8 xd,yd,zd
      real*8 xu,yu,zu
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
c
c
c     zero out the multipole and polarization energies
c
      nem = 0
      nep = 0
      em = 0.0d0
      ep = 0.0d0
      do i = 1, n
         aem(i) = 0.0d0
         aep(i) = 0.0d0
      end do
c
c     set the energy unit conversion factor
c
      f = electric / dielec
c
c     check the sign of multipole components at chiral sites
c
      call chkpole_3b
c
c     rotate the multipole components into the global frame
c
      call rotpole_3b
c
c     compute the induced dipoles at each polarizable atom
c
      call induce_3b
c
c     compute the reciprocal space part of the Ewald summation
c
      call emrecip_3b
c
c     compute the real space part of the Ewald summation
c
      call ereal3d_3b (eintra)
c
c     compute the self-energy part of the Ewald summation
c
      term = 2.0d0 * aewald * aewald
      fterm = -f * aewald / sqrtpi
      do i = 1, npole
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
         cii = ci*ci
         dii = dix*dix + diy*diy + diz*diz
         qii = qixx*qixx + qiyy*qiyy + qizz*qizz
     &            + 2.0d0*(qixy*qixy+qixz*qixz+qiyz*qiyz)
         uii = dix*uix + diy*uiy + diz*uiz
         e = fterm * (cii + term*(dii/3.0d0+2.0d0*term*qii/5.0d0))
         ei = fterm * term * uii / 3.0d0
         nem = nem + 1
         nep = nep + 1
         em = em + e
         ep = ep + ei
         aem(i) = aem(i) + e
         aep(i) = aep(i) + ei
      end do
c
c     compute the cell dipole boundary correction term
c
      if (boundary .eq. 'VACUUM') then
         xd = 0.0d0
         yd = 0.0d0
         zd = 0.0d0
         xu = 0.0d0
         yu = 0.0d0
         zu = 0.0d0
         do i = 1, npole
            ii = ipole(i)
            dix = rpole(2,i)
            diy = rpole(3,i)
            diz = rpole(4,i)
            uix = uind(1,i)
            uiy = uind(2,i)
            uiz = uind(3,i)
            xd = xd + dix + rpole(1,i)*x(ii)
            yd = yd + diy + rpole(1,i)*y(ii)
            zd = zd + diz + rpole(1,i)*z(ii)
            xu = xu + uix
            yu = yu + uiy
            zu = zu + uiz
         end do
         term = (2.0d0/3.0d0) * f * (pi/volbox)
         nem = nem + 1
         nep = nep + 1
         em = em + term*(xd*xd+yd*yd+zd*zd)
         ep = ep + term*(xd*xu+yd*yu+zd*zu)
      end if
c
c     intermolecular energy is total minus intramolecular part
c
      einter = einter + em + ep - eintra
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine ereal3d_3b  --  real space mpole analysis via list  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "ereal3d" evaluates the real space portion of the Ewald sum
c     energy due to atomic multipole interactions and dipole
c     polarizability and partitions the energy among the atoms
c     using a pairwise neighbor list
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine ereal3d_3b (eintra)
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'analyz.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'energi.i'
      include 'ewald.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'molcul.i'
      include 'mplpot.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'potent.i'
      include 'shunt.i'
      integer i,j,k,m
      integer ii,kk,kkk
      real*8 e,ei,eintra
      real*8 f,erfc
      real*8 r,r2,xr,yr,zr
      real*8 bfac,exp2a
      real*8 efix,eifix
      real*8 ralpha
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 alsq2,alsq2n
      real*8 rr1,rr3,rr5
      real*8 rr7,rr9
      real*8 ci,dix,diy,diz
      real*8 uix,uiy,uiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 ukx,uky,ukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 qix,qiy,qiz
      real*8 qkx,qky,qkz
      real*8 scale3,scale5
      real*8 scale7
      real*8 sc(10),sci(8)
      real*8 gl(0:4),gli(3)
      real*8 bn(0:4)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: pscale(:)
      logical header,huge
      logical muse,puse
      character*6 mode
      external erfc
c
c
c     zero out the intramolecular portion of the Ewald energy
c
      eintra = 0.0d0
      if (npole .eq. 0)  return
      header = .true.
c
c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (pscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         mscale(i) = 1.0d0
         pscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
      mode = 'EWALD'
      call switch (mode)
c
c     compute the real space portion of the Ewald summation
c
      do i = 1, npole-1
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         uix = uind(1,i)
         uiy = uind(2,i)
         uiz = uind(3,i)
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
         do kkk = 1, nelst(i)
            k = elst(kkk,i)
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dkx = rpole(2,k)
               dky = rpole(3,k)
               dkz = rpole(4,k)
               qkxx = rpole(5,k)
               qkxy = rpole(6,k)
               qkxz = rpole(7,k)
               qkyy = rpole(9,k)
               qkyz = rpole(10,k)
               qkzz = rpole(13,k)
               ukx = uind(1,k)
               uky = uind(2,k)
               ukz = uind(3,k)
c
c     calculate the error function damping terms
c
               ralpha = aewald * r
               bn(0) = erfc(ralpha) / r
               alsq2 = 2.0d0 * aewald**2
               alsq2n = 0.0d0
               if (aewald .gt. 0.0d0)
     &            alsq2n = 1.0d0 / (sqrtpi*aewald)
               exp2a = exp(-ralpha**2)
               do m = 1, 4
                  bfac = dble(m+m-1)
                  alsq2n = alsq2 * alsq2n
                  bn(m) = (bfac*bn(m-1)+alsq2n*exp2a) / r2
               end do
c
c     construct some intermediate quadrupole values
c
               qix = qixx*xr + qixy*yr + qixz*zr
               qiy = qixy*xr + qiyy*yr + qiyz*zr
               qiz = qixz*xr + qiyz*yr + qizz*zr
               qkx = qkxx*xr + qkxy*yr + qkxz*zr
               qky = qkxy*xr + qkyy*yr + qkyz*zr
               qkz = qkxz*xr + qkyz*yr + qkzz*zr
c
c     calculate the scalar products for permanent multipoles
c
               sc(2) = dix*dkx + diy*dky + diz*dkz
               sc(3) = dix*xr + diy*yr + diz*zr
               sc(4) = dkx*xr + dky*yr + dkz*zr
               sc(5) = qix*xr + qiy*yr + qiz*zr
               sc(6) = qkx*xr + qky*yr + qkz*zr
               sc(7) = qix*dkx + qiy*dky + qiz*dkz
               sc(8) = qkx*dix + qky*diy + qkz*diz
               sc(9) = qix*qkx + qiy*qky + qiz*qkz
               sc(10) = 2.0d0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz)
     &                     + qixx*qkxx + qiyy*qkyy + qizz*qkzz
c
c     calculate the scalar products for polarization components
c
               sci(2) = uix*dkx + dix*ukx + uiy*dky
     &                     + diy*uky + uiz*dkz + diz*ukz
               sci(3) = uix*xr + uiy*yr + uiz*zr
               sci(4) = ukx*xr + uky*yr + ukz*zr
               sci(7) = qix*ukx + qiy*uky + qiz*ukz
               sci(8) = qkx*uix + qky*uiy + qkz*uiz
c
c     calculate the gl functions for permanent multipoles
c
               gl(0) = ci*ck
               gl(1) = ck*sc(3) - ci*sc(4) + sc(2)
               gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
     &                    + 2.0d0*(sc(7)-sc(8)+sc(10))
               gl(3) = sc(3)*sc(6) - sc(4)*sc(5) - 4.0d0*sc(9)
               gl(4) = sc(5)*sc(6)
c
c     calculate the gl functions for polarization components
c
               gli(1) = ck*sci(3) - ci*sci(4) + sci(2)
               gli(2) = 2.0d0*(sci(7)-sci(8)) - sci(3)*sc(4)
     &                     - sc(3)*sci(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
c
c     compute the energy contributions for this interaction
c
               e = gl(0)*bn(0) + gl(1)*bn(1) + gl(2)*bn(2)
     &                + gl(3)*bn(3) + gl(4)*bn(4)
               ei = gli(1)*bn(1) + gli(2)*bn(2) + gli(3)*bn(3)
c
c     full real space energies needed for scaled interactions
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               scale3 = pscale(kk)
               scale5 = pscale(kk)
               scale7 = pscale(kk)
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = scale3 * (1.0d0-expdamp)
                     scale5 = scale5 * (1.0d0-(1.0d0-damp)*expdamp)
                     scale7 = scale7 * (1.0d0-(1.0d0-damp+0.6d0*damp**2)
     &                                               *expdamp)
                  end if
               end if
               efix = gl(0)*rr1 + gl(1)*rr3 + gl(2)*rr5
     &                   + gl(3)*rr7 + gl(4)*rr9
               eifix = gli(1)*rr3*(1.0d0-scale3)
     &                    + gli(2)*rr5*(1.0d0-scale5)
     &                    + gli(3)*rr7*(1.0d0-scale7)
c
c     apply the energy adjustments for scaled interactions
c
               e = e - efix*(1.0d0-mscale(kk))
               ei = ei - eifix
c
c     increment the overall multipole and polarization energies
c
               e = f * e
               ei = 0.5d0 * f * ei
               muse = use_mpole
               puse = use_polar
               if (muse)  nem = nem + 1
               if (puse)  nep = nep + 1
               em = em + e
               ep = ep + ei
               aem(ii) = aem(ii) + 0.5d0*e
               aem(kk) = aem(kk) + 0.5d0*e
               aep(ii) = aep(ii) + 0.5d0*ei
               aep(kk) = aep(kk) + 0.5d0*ei
c
c     increment the total intramolecular energy
c
               efix = f * efix * mscale(kk)
               eifix = gli(1)*rr3*scale3 + gli(2)*rr5*scale5
     &                    + gli(3)*rr7*scale7
               eifix = 0.5d0 * f * eifix
               if (molcule(ii) .eq. molcule(kk)) then
                  eintra = eintra + efix + eifix
               end if
c
c     print a message if the energy of this interaction is large
c
               huge = (max(abs(efix),abs(eifix)) .gt. 100.0d0)
               if (debug .or. (verbose.and.huge)) then
                  if (muse .or. puse) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Real Space Multipole and',
     &                             ' Polarization Interactions :',
     &                          //,' Type',14x,'Atom Names',
     &                             15x,'Distance',6x,'Energies',
     &                             ' (MPol,Polar)',/)
                     end if
                     write (iout,20)  ii,name(ii),kk,name(kk),r,
     &                                efix,eifix
   20                format (' M-Pole',4x,2(i7,'-',a3),9x,
     &                          f10.4,2x,2f12.4)
                  end if
               end if
            end if
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (pscale)
      return
      end
c
c
c     ###################################################################
c     ##                                                               ##
c     ##  subroutine emrecip_3b  --  PME recip space multipole energy  ##
c     ##                                                               ##
c     ###################################################################
c
c
c     "emrecip" evaluates the reciprocal space portion of the particle
c     mesh Ewald energy due to atomic multipole interactions and
c     dipole polarizability
c
c     literature reference:
c
c     C. Sagui, L. G. Pedersen and T. A. Darden, "Towards an Accurate
c     Representation of Electrostatics in Classical Force Fields:
c     Efficient Implementation of Multipolar Interactions in
c     Biomolecular Simulations", Journal of Chemical Physics, 120,
c     73-87 (2004)
c
c     modifications for nonperiodic systems suggested by Tom Darden
c     during May 2007
c
c
      subroutine emrecip_3b
      implicit none
      include 'sizes.i'
      include 'bound.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'energi.i'
      include 'ewald.i'
      include 'math.i'
      include 'mpole.i'
      include 'pme.i'
      include 'polar.i'
      include 'potent.i'
      include 'molcul.i'
      include 'combo.i'
      integer l1
      integer i,j,k,ntot
      integer k1,k2,k3
      integer m1,m2,m3
      integer nff,nf1,nf2,nf3
      real*8 e,r1,r2,r3
      real*8 h1,h2,h3
      real*8 volterm,denom
      real*8 hsq,expterm
      real*8 term,pterm
      real*8 a(3,3)
      real*8, allocatable :: fuind(:,:)
      real*8, allocatable :: cmp(:,:)
      real*8, allocatable :: fmp(:,:)
      real*8, allocatable :: fphi(:,:)
c
c
c     return if the Ewald coefficient is zero
c
      if (aewald .lt. 1.0d-6)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (fuind(3,npole))
      allocate (cmp(10,npole))
      allocate (fmp(10,npole))
      allocate (fphi(20,npole))
c
c     copy the multipole moments into local storage areas
c
      do l1 = 1, npole3b
         i = pnum(l1)
         cmp(1,i) = rpole(1,i)
         cmp(2,i) = rpole(2,i)
         cmp(3,i) = rpole(3,i)
         cmp(4,i) = rpole(4,i)
         cmp(5,i) = rpole(5,i)
         cmp(6,i) = rpole(9,i)
         cmp(7,i) = rpole(13,i)
         cmp(8,i) = 2.0d0 * rpole(6,i)
         cmp(9,i) = 2.0d0 * rpole(7,i)
         cmp(10,i) = 2.0d0 * rpole(10,i)
      end do
c
c     compute B-spline coefficients and spatial decomposition
c
      if (.not. use_polar) then
         call bspline_fill
         call table_fill
      end if
c
c     convert Cartesian multipoles to fractional coordinates
c
      call cmp_to_fmp_3b (cmp,fmp)
c
c     assign PME grid and perform 3-D FFT forward transform
c
      call grid_mpole_3b (fmp)
      call fftfront
c
c     make the scalar summation over reciprocal lattice
c
      ntot = nfft1 * nfft2 * nfft3
      pterm = (pi/aewald)**2
      volterm = pi * volbox
      nff = nfft1 * nfft2
      nf1 = (nfft1+1) / 2
      nf2 = (nfft2+1) / 2
      nf3 = (nfft3+1) / 2
      do i = 1, ntot-1
         k3 = i/nff + 1
         j = i - (k3-1)*nff
         k2 = j/nfft1 + 1
         k1 = j - (k2-1)*nfft1 + 1
         m1 = k1 - 1
         m2 = k2 - 1
         m3 = k3 - 1
         if (k1 .gt. nf1)  m1 = m1 - nfft1
         if (k2 .gt. nf2)  m2 = m2 - nfft2
         if (k3 .gt. nf3)  m3 = m3 - nfft3
         r1 = dble(m1)
         r2 = dble(m2)
         r3 = dble(m3)
         h1 = recip(1,1)*r1 + recip(1,2)*r2 + recip(1,3)*r3
         h2 = recip(2,1)*r1 + recip(2,2)*r2 + recip(2,3)*r3
         h3 = recip(3,1)*r1 + recip(3,2)*r2 + recip(3,3)*r3
         hsq = h1*h1 + h2*h2 + h3*h3
         term = -pterm * hsq
         expterm = 0.0d0
         if (term .gt. -50.0d0) then
            denom = volterm*hsq*bsmod1(k1)*bsmod2(k2)*bsmod3(k3)
            expterm = exp(term) / denom
            if (.not. use_bounds) then
               expterm = expterm * (1.0d0-cos(pi*xbox*sqrt(hsq)))
            else if (octahedron) then
               if (mod(m1+m2+m3,2) .ne. 0)  expterm = 0.0d0
            end if
         end if
         qfac(k1,k2,k3) = expterm
      end do
c
c     account for the zeroth grid point for a finite system
c
      qfac(1,1,1) = 0.0d0
      if (.not. use_bounds) then
         expterm = 0.5d0 * pi / xbox
         qfac(1,1,1) = expterm
      end if
c
c     complete the transformation of the charge grid
c
      do k = 1, nfft3
         do j = 1, nfft2
            do i = 1, nfft1
               term = qfac(i,j,k)
               qgrid(1,i,j,k) = term * qgrid(1,i,j,k)
               qgrid(2,i,j,k) = term * qgrid(2,i,j,k)
            end do
         end do
      end do
c
c     perform 3-D FFT backward transform and get potential
c
      call fftback
      call fphi_mpole_3b (fphi)
c
c     sum over multipoles and increment total multipole energy
c
      e = 0.0d0
      do l1 = 1, npole3b
         i = pnum(l1)
         do k = 1, 10
            e = e + fmp(k,i)*fphi(k,i)
         end do
      end do
      e = 0.5d0 * electric * e
      em = em + e
c
c     convert Cartesian induced dipoles to fractional coordinates
c
      if (use_polar) then
         do i = 1, 3
            a(1,i) = dble(nfft1) * recip(i,1)
            a(2,i) = dble(nfft2) * recip(i,2)
            a(3,i) = dble(nfft3) * recip(i,3)
         end do
         do l1 = 1, npole3b
            i = pnum(l1)
            do k = 1, 3
               fuind(k,i) = a(k,1)*uind(1,i) + a(k,2)*uind(2,i)
     &                         + a(k,3)*uind(3,i)
            end do
         end do
c
c     sum over induced dipoles and increment total induced energy
c
         e = 0.0d0
         do l1 = 1, npole3b
            i = pnum(l1)
            do k = 1, 3
               e = e + fuind(k,i)*fphi(k+1,i)
            end do
         end do
         e = 0.5d0 * electric * e
         ep = ep + e
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (fuind)
      deallocate (cmp)
      deallocate (fmp)
      deallocate (fphi)
      return
      end

c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine erecip3_3b   mpole Ewald reciprocal analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "erecip3_3b" evaluates the reciprocal space portion the regular
c     Ewald summation energy due to atomic multipole interactions
c     and dipole polarizability, and prints information about the
c     energy over the reciprocal lattice vectors
c
c     literature reference:
c
c     W. Smith, "Point Multipoles in the Ewald Summation (Revisited)",
c     CCP5 Newsletter, 46, 18-30, 1998  (see http://www.ccp5.org/)
c
c
      subroutine erecip3_3b
      implicit none
      include 'sizes.i'
      include 'action.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'chgpot.i'
      include 'energi.i'
      include 'ewald.i'
      include 'ewreg.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'units.i'
      include 'bound.i'
      include 'potent.i'
      include 'molcul.i'
      include 'combo.i'
      integer i,j,k,l,ii,l1
      integer jmin,jmax
      integer kmin,kmax
      integer lmin,lmax
      real*8 e,ei,f,term,cut
      real*8 eterm,expterm
      real*8 xfr,yfr,zfr
      real*8 rj,rk,rl
      real*8 h1,h2,h3
      real*8 hsq,hleng
      real*8 t1,t2,t3,t4
      real*8 ck,dk,qk,uk
      real*8 q1,q2,q3
      real*8 ckr,skr
      real*8 cm(maxatm)
      real*8 dm(3,maxatm)
      real*8 qm(9,maxatm)
      real*8 um(3,maxatm)
      real*8 cjk(maxatm)
      real*8 sjk(maxatm)
      real*8 frecip
      logical header,huge
c
c
c     return if the Ewald coefficient is zero
c
      frecip = 0.5d0

      if (aewald .lt. 1.0d-6)  return
      header = .true.
      f = electric / dielec
      term = -0.25d0 / aewald**2
      eterm = 4.0d0 * pi * f / volbox
c
c     set the number of vectors based on box dimensions
c
      cut = 4.0d0 * pi * pi * frecip
      jmin = 0
      kmin = 0
      lmin = 1
      jmax = min(maxvec,int(frecip/recip(1,1)))
      kmax = min(maxvec,int(frecip/recip(2,2)))
      lmax = min(maxvec,int(frecip/recip(3,3)))
c
c     copy the multipole moments into local storage areas
c
      do l1 = 1, npole3b
         i=pnum(l1)
         cm(i) = rpole(1,i)
         dm(1,i) = rpole(2,i)
         dm(2,i) = rpole(3,i)
         dm(3,i) = rpole(4,i)
         qm(1,i) = rpole(5,i)
         qm(2,i) = rpole(6,i)
         qm(3,i) = rpole(7,i)
         qm(4,i) = rpole(8,i)
         qm(5,i) = rpole(9,i)
         qm(6,i) = rpole(10,i)
         qm(7,i) = rpole(11,i)
         qm(8,i) = rpole(12,i)
         qm(9,i) = rpole(13,i)
         um(1,i) = uind(1,i)
         um(2,i) = uind(2,i)
         um(3,i) = uind(3,i)
      end do
c
c     calculate and store the exponential factors
c
      do l1 = 1, npole3b
         i=pnum(l1)
         ii = ipole(i)
         zfr = (z(ii)/gamma_term) / zbox
         yfr = ((y(ii)-zfr*zbox*beta_term)/gamma_sin) / ybox
         xfr = (x(ii)-yfr*ybox*gamma_cos-zfr*zbox*beta_cos) / xbox
         xfr = 2.0d0 * pi * xfr
         yfr = 2.0d0 * pi * yfr
         zfr = 2.0d0 * pi * zfr
         ejc(i,0) = 1.0d0
         ejs(i,0) = 0.0d0
         ekc(i,0) = 1.0d0
         eks(i,0) = 0.0d0
         elc(i,0) = 1.0d0
         els(i,0) = 0.0d0
         ejc(i,1) = cos(xfr)
         ejs(i,1) = sin(xfr)
         ekc(i,1) = cos(yfr)
         eks(i,1) = sin(yfr)
         elc(i,1) = cos(zfr)
         els(i,1) = sin(zfr)
         ekc(i,-1) = ekc(i,1)
         eks(i,-1) = -eks(i,1)
         elc(i,-1) = elc(i,1)
         els(i,-1) = -els(i,1)
         do j = 2, jmax
            ejc(i,j) = ejc(i,j-1)*ejc(i,1) - ejs(i,j-1)*ejs(i,1)
            ejs(i,j) = ejs(i,j-1)*ejc(i,1) + ejc(i,j-1)*ejs(i,1)
         end do
         do j = 2, kmax
            ekc(i,j) = ekc(i,j-1)*ekc(i,1) - eks(i,j-1)*eks(i,1)
            eks(i,j) = eks(i,j-1)*ekc(i,1) + ekc(i,j-1)*eks(i,1)
            ekc(i,-j) = ekc(i,j)
            eks(i,-j) = -eks(i,j)
         end do
         do j = 2, lmax
            elc(i,j) = elc(i,j-1)*elc(i,1) - els(i,j-1)*els(i,1)
            els(i,j) = els(i,j-1)*elc(i,1) + elc(i,j-1)*els(i,1)
            elc(i,-j) = elc(i,j)
            els(i,-j) = -els(i,j)
         end do
      end do
c
c     loop over all k vectors from the reciprocal lattice
c
      do j = jmin, jmax
         rj = 2.0d0 * pi * dble(j)
         do k = kmin, kmax
            rk = 2.0d0 * pi * dble(k)
            do l1 = 1, npole3b
               i=pnum(l1)
               cjk(i) = ejc(i,j)*ekc(i,k) - ejs(i,j)*eks(i,k)
               sjk(i) = ejs(i,j)*ekc(i,k) + ejc(i,j)*eks(i,k)
            end do
            do l = lmin, lmax
               rl = 2.0d0 * pi * dble(l)
               h1 = recip(1,1)*rj
               h2 = recip(2,1)*rj + recip(2,2)*rk
               h3 = recip(3,1)*rj + recip(3,2)*rk + recip(3,3)*rl
               hsq = h1*h1 + h2*h2 + h3*h3
               if (hsq .le. cut) then
                  t1 = 0.0d0
                  t2 = 0.0d0
                  t3 = 0.0d0
                  t4 = 0.0d0
                  do l1 = 1, npole3b
                     i=pnum(l1)
                     ckr = cjk(i)*elc(i,l) - sjk(i)*els(i,l)
                     skr = sjk(i)*elc(i,l) + cjk(i)*els(i,l)
c
c     alternative phase cosine and sine is slower but avoids storage
c
c                    ii = ipole(i)
c                    phi = x(i)*h1 + y(i)*h2 + z(i)*h3
c                    ckr = cos(phi)
c                    skr = sin(phi)
                     ck = cm(i)
                     dk = h1*dm(1,i) + h2*dm(2,i) + h3*dm(3,i)
                     q1 = h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i)
                     q2 = h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i)
                     q3 = h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i)
                     qk = h1*q1 + h2*q2 + h3*q3
                     uk = h1*um(1,i) + h2*um(2,i) + h3*um(3,i)
                     t1 = t1 + (ck-qk)*skr + dk*ckr
                     t2 = t2 + (ck-qk)*ckr - dk*skr
                     t3 = t3 + uk*ckr
                     t4 = t4 - uk*skr
                  end do
                  expterm = eterm * exp(term*hsq) / hsq
                  if (octahedron) then
                     if (mod(j+k+l,2) .ne. 0)  expterm = 0.0d0
                  end if
                  e = expterm * (t1*t1+t2*t2)
                  ei = expterm * (t1*t3+t2*t4)
                  nem = nem + 1
                  nep = nep + 1
                  em = em + e
                  ep = ep + ei
c
c     print a message if the energy of this interaction is large
c
                  huge = (max(abs(e),abs(ei)) .gt. 10.0d0)
                  if (debug .or. (verbose.and.huge)) then
                     if (header) then
                        header = .false.
                        write (iout,10)
   10                   format (/,' Reciprocal Space Multipole',
     &                             ' and Polarization Terms :',
     &                          //,' Type',15x,'K-Vector',9x,'Fraction',
     &                             6x,'Energies (MPole, Polar)',/)
                     end if
                     hleng = hsq / (4.0d0*pi*pi)
                     write (iout,20)  j,k,l,hleng,e,ei
   20                format (' M-Pole',7x,3i5,7x,f8.4,4x,2f12.4)
                  end if
               end if
            end do
            lmin = -lmax
         end do
         kmin = -kmax
      end do
      return
      end

