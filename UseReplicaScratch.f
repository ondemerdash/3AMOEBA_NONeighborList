c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
      if (use_replica) then
c
c     calculate interaction with other unit cells
c

c      do i = 1, npole
      do l1= 1, npole3b
         i=pnum(l1)
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         di(1) = rpole(2,i)
         di(2) = rpole(3,i)
         di(3) = rpole(4,i)
         qi(1) = rpole(5,i)
         qi(2) = rpole(6,i)
         qi(3) = rpole(7,i)
         qi(4) = rpole(8,i)
         qi(5) = rpole(9,i)
         qi(6) = rpole(10,i)
         qi(7) = rpole(11,i)
         qi(8) = rpole(12,i)
         qi(9) = rpole(13,i)
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
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
         end do
c         do k = i, npole
         do l3 = l1,npole3b
            k=pnum(l3)
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
            if (.not. proceed)  goto 20
            do jcell = 1, ncell
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call imager (xr,yr,zr,jcell)
            r2 = xr*xr + yr*yr + zr*zr
c            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dk(1) = rpole(2,k)
               dk(2) = rpole(3,k)
               dk(3) = rpole(4,k)
               qk(1) = rpole(5,k)
               qk(2) = rpole(6,k)
               qk(3) = rpole(7,k)
               qk(4) = rpole(8,k)
               qk(5) = rpole(9,k)
               qk(6) = rpole(10,k)
               qk(7) = rpole(11,k)
               qk(8) = rpole(12,k)
               qk(9) = rpole(13,k)
c
c     apply Thole polarization damping to scale factors
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               do j = 1, 3
                  ddsc3(j) = 0.0d0
                  ddsc5(j) = 0.0d0
                  ddsc7(j) = 0.0d0
               end do
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                     scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *expdamp
                     temp3 = -3.0d0 * damp * expdamp / r2
                     temp5 = -damp
                     temp7 = -0.2d0 - 0.6d0*damp
                     ddsc3(1) = temp3 * xr
                     ddsc3(2) = temp3 * yr
                     ddsc3(3) = temp3 * zr
                     ddsc5(1) = temp5 * ddsc3(1)
                     ddsc5(2) = temp5 * ddsc3(2)
                     ddsc5(3) = temp5 * ddsc3(3)
                     ddsc7(1) = temp7 * ddsc5(1)
                     ddsc7(2) = temp7 * ddsc5(2)
                     ddsc7(3) = temp7 * ddsc5(3)
                  end if
               end if
               scale3i = scale3
               scale5i = scale5
               scale7i = scale7
               dsc3 = scale3
               dsc5 = scale5
               dsc7 = scale7
               psc3 = scale3
               psc5 = scale5
               psc7 = scale7
               if (use_polymer) then
                  if (r2 .le. polycut2) then
                     scale3i = scale3i * uscale(kk)
                     scale5i = scale5i * uscale(kk)
                     scale7i = scale7i * uscale(kk)
                     dsc3 = dsc3 * dscale(kk)
                     dsc5 = dsc5 * dscale(kk)
                     dsc7 = dsc7 * dscale(kk)
                     psc3 = psc3 * pscale(kk)
                     psc5 = psc5 * pscale(kk)
                     psc7 = psc7 * pscale(kk)
                  end if
               end if
c
c     construct necessary auxiliary vectors
c

c               dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
c               dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
c               dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
               dixuk(1) = di(2)*uind(3,l3) - di(3)*uind(2,l3)
               dixuk(2) = di(3)*uind(1,l3) - di(1)*uind(3,l3)
               dixuk(3) = di(1)*uind(2,l3) - di(2)*uind(1,l3)
               dkxui(1) = dk(2)*uind(3,l1) - dk(3)*uind(2,l1)
               dkxui(2) = dk(3)*uind(1,l1) - dk(1)*uind(3,l1)
               dkxui(3) = dk(1)*uind(2,l1) - dk(2)*uind(1,l1)
               dixukp(1) = di(2)*uinp(3,l3) - di(3)*uinp(2,l3)
               dixukp(2) = di(3)*uinp(1,l3) - di(1)*uinp(3,l3)
               dixukp(3) = di(1)*uinp(2,l3) - di(2)*uinp(1,l3)
               dkxuip(1) = dk(2)*uinp(3,l1) - dk(3)*uinp(2,l1)
               dkxuip(2) = dk(3)*uinp(1,l1) - dk(1)*uinp(3,l1)
               dkxuip(3) = dk(1)*uinp(2,l1) - dk(2)*uinp(1,l1)
               dixr(1) = di(2)*zr - di(3)*yr
               dixr(2) = di(3)*xr - di(1)*zr
               dixr(3) = di(1)*yr - di(2)*xr
               dkxr(1) = dk(2)*zr - dk(3)*yr
               dkxr(2) = dk(3)*xr - dk(1)*zr
               dkxr(3) = dk(1)*yr - dk(2)*xr
               qir(1) = qi(1)*xr + qi(4)*yr + qi(7)*zr
               qir(2) = qi(2)*xr + qi(5)*yr + qi(8)*zr
               qir(3) = qi(3)*xr + qi(6)*yr + qi(9)*zr
               qkr(1) = qk(1)*xr + qk(4)*yr + qk(7)*zr
               qkr(2) = qk(2)*xr + qk(5)*yr + qk(8)*zr
               qkr(3) = qk(3)*xr + qk(6)*yr + qk(9)*zr
c               qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
c               qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
c               qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
c               qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
c               qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
c               qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
c               qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)
c     &                       - qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
c               qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)
c     &                       - qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
c               qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)
c     &                       - qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
               rxqir(1) = yr*qir(3) - zr*qir(2)
               rxqir(2) = zr*qir(1) - xr*qir(3)
               rxqir(3) = xr*qir(2) - yr*qir(1)
               rxqkr(1) = yr*qkr(3) - zr*qkr(2)
               rxqkr(2) = zr*qkr(1) - xr*qkr(3)
               rxqkr(3) = xr*qkr(2) - yr*qkr(1)
c               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
c               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
c               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
c               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
c               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
c               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
c               qkrxqir(1) = qkr(2)*qir(3) - qkr(3)*qir(2)
c               qkrxqir(2) = qkr(3)*qir(1) - qkr(1)*qir(3)
c               qkrxqir(3) = qkr(1)*qir(2) - qkr(2)*qir(1)
c               qidk(1) = qi(1)*dk(1) + qi(4)*dk(2) + qi(7)*dk(3)
c               qidk(2) = qi(2)*dk(1) + qi(5)*dk(2) + qi(8)*dk(3)
c               qidk(3) = qi(3)*dk(1) + qi(6)*dk(2) + qi(9)*dk(3)
c               qkdi(1) = qk(1)*di(1) + qk(4)*di(2) + qk(7)*di(3)
c               qkdi(2) = qk(2)*di(1) + qk(5)*di(2) + qk(8)*di(3)
c               qkdi(3) = qk(3)*di(1) + qk(6)*di(2) + qk(9)*di(3)
               qiuk(1) = qi(1)*uind(1,l3) + qi(4)*uind(2,l3)
     &                      + qi(7)*uind(3,l3)
               qiuk(2) = qi(2)*uind(1,l3) + qi(5)*uind(2,l3)
     &                      + qi(8)*uind(3,l3)
               qiuk(3) = qi(3)*uind(1,l3) + qi(6)*uind(2,l3)
     &                      + qi(9)*uind(3,l3)
               qkui(1) = qk(1)*uind(1,l1) + qk(4)*uind(2,l1)
     &                      + qk(7)*uind(3,l1)
               qkui(2) = qk(2)*uind(1,l1) + qk(5)*uind(2,l1)
     &                      + qk(8)*uind(3,l1)
               qkui(3) = qk(3)*uind(1,l1) + qk(6)*uind(2,l1)
     &                      + qk(9)*uind(3,l1)
               qiukp(1) = qi(1)*uinp(1,l3) + qi(4)*uinp(2,l3)
     &                       + qi(7)*uinp(3,l3)
               qiukp(2) = qi(2)*uinp(1,l3) + qi(5)*uinp(2,l3)
     &                       + qi(8)*uinp(3,l3)
               qiukp(3) = qi(3)*uinp(1,l3) + qi(6)*uinp(2,l3)
     &                       + qi(9)*uinp(3,l3)
               qkuip(1) = qk(1)*uinp(1,l1) + qk(4)*uinp(2,l1)
     &                       + qk(7)*uinp(3,l1)
               qkuip(2) = qk(2)*uinp(1,l1) + qk(5)*uinp(2,l1)
     &                       + qk(8)*uinp(3,l1)
               qkuip(3) = qk(3)*uinp(1,l1) + qk(6)*uinp(2,l1)
     &                       + qk(9)*uinp(3,l1)
c               dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
c               dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
c               dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
c               dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
c               dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
c               dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
               uixqkr(1) = uind(2,l1)*qkr(3) - uind(3,l1)*qkr(2)
               uixqkr(2) = uind(3,l1)*qkr(1) - uind(1,l1)*qkr(3)
               uixqkr(3) = uind(1,l1)*qkr(2) - uind(2,l1)*qkr(1)
               ukxqir(1) = uind(2,l3)*qir(3) - uind(3,l3)*qir(2)
               ukxqir(2) = uind(3,l3)*qir(1) - uind(1,l3)*qir(3)
               ukxqir(3) = uind(1,l3)*qir(2) - uind(2,l3)*qir(1)
               uixqkrp(1) = uinp(2,l1)*qkr(3) - uinp(3,l1)*qkr(2)
               uixqkrp(2) = uinp(3,l1)*qkr(1) - uinp(1,l1)*qkr(3)
               uixqkrp(3) = uinp(1,l1)*qkr(2) - uinp(2,l1)*qkr(1)
               ukxqirp(1) = uinp(2,l3)*qir(3) - uinp(3,l3)*qir(2)
               ukxqirp(2) = uinp(3,l3)*qir(1) - uinp(1,l3)*qir(3)
               ukxqirp(3) = uinp(1,l3)*qir(2) - uinp(2,l3)*qir(1)
c               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
c               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
c               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
c               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
c               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
c               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
               rxqiuk(1) = yr*qiuk(3) - zr*qiuk(2)
               rxqiuk(2) = zr*qiuk(1) - xr*qiuk(3)
               rxqiuk(3) = xr*qiuk(2) - yr*qiuk(1)
               rxqkui(1) = yr*qkui(3) - zr*qkui(2)
               rxqkui(2) = zr*qkui(1) - xr*qkui(3)
               rxqkui(3) = xr*qkui(2) - yr*qkui(1)
               rxqiukp(1) = yr*qiukp(3) - zr*qiukp(2)
               rxqiukp(2) = zr*qiukp(1) - xr*qiukp(3)
               rxqiukp(3) = xr*qiukp(2) - yr*qiukp(1)
               rxqkuip(1) = yr*qkuip(3) - zr*qkuip(2)
               rxqkuip(2) = zr*qkuip(1) - xr*qkuip(3)
               rxqkuip(3) = xr*qkuip(2) - yr*qkuip(1)
c
c     calculate scalar products for permanent components
c
               sc(2) = di(1)*dk(1) + di(2)*dk(2) + di(3)*dk(3)
               sc(3) = di(1)*xr + di(2)*yr + di(3)*zr
               sc(4) = dk(1)*xr + dk(2)*yr + dk(3)*zr
               sc(5) = qir(1)*xr + qir(2)*yr + qir(3)*zr
               sc(6) = qkr(1)*xr + qkr(2)*yr + qkr(3)*zr
               sc(7) = qir(1)*dk(1) + qir(2)*dk(2) + qir(3)*dk(3)
               sc(8) = qkr(1)*di(1) + qkr(2)*di(2) + qkr(3)*di(3)
               sc(9) = qir(1)*qkr(1) + qir(2)*qkr(2) + qir(3)*qkr(3)
               sc(10) = qi(1)*qk(1) + qi(2)*qk(2) + qi(3)*qk(3)
     &                     + qi(4)*qk(4) + qi(5)*qk(5) + qi(6)*qk(6)
     &                     + qi(7)*qk(7) + qi(8)*qk(8) + qi(9)*qk(9)
c
c     calculate scalar products for induced components
c
               sci(1) = uind(1,l1)*dk(1) + uind(2,l1)*dk(2)
     &                     + uind(3,l1)*dk(3) + di(1)*uind(1,l3)
     &                     + di(2)*uind(2,l3) + di(3)*uind(3,l3)
               sci(2) = uind(1,l1)*uind(1,l3) + uind(2,l1)*uind(2,l3)
     &                     + uind(3,l1)*uind(3,l3)
               sci(3) = uind(1,l1)*xr + uind(2,l1)*yr + uind(3,l1)*zr
               sci(4) = uind(1,l3)*xr + uind(2,l3)*yr + uind(3,l3)*zr
               sci(7) = qir(1)*uind(1,l3) + qir(2)*uind(2,l3)
     &                     + qir(3)*uind(3,l3)
               sci(8) = qkr(1)*uind(1,l1) + qkr(2)*uind(2,l1)
     &                     + qkr(3)*uind(3,l1)
               scip(1) = uinp(1,l1)*dk(1) + uinp(2,l1)*dk(2)
     &                      + uinp(3,l1)*dk(3) + di(1)*uinp(1,l3)
     &                      + di(2)*uinp(2,l3) + di(3)*uinp(3,l3)
               scip(2) = uind(1,l1)*uinp(1,l3)+uind(2,l1)*uinp(2,l3)
     &                     + uind(3,l1)*uinp(3,l3)+uinp(1,l1)*uind(1,l3)
     &                     + uinp(2,l1)*uind(2,l3)+uinp(3,l1)*uind(3,l3)
               scip(3) = uinp(1,l1)*xr + uinp(2,l1)*yr + uinp(3,l1)*zr
               scip(4) = uinp(1,l3)*xr + uinp(2,l3)*yr + uinp(3,l3)*zr
               scip(7) = qir(1)*uinp(1,l3) + qir(2)*uinp(2,l3)
     &                      + qir(3)*uinp(3,l3)
               scip(8) = qkr(1)*uinp(1,l1) + qkr(2)*uinp(2,l1)
     &                      + qkr(3)*uinp(3,l1)
c
c     calculate the gl functions for permanent components
c

c               gl(0) = ci*ck
c               gl(1) = ck*sc(3) - ci*sc(4)
c               gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
c               gl(3) = sc(3)*sc(6) - sc(4)*sc(5)
c               gl(4) = sc(5)*sc(6)
c               gl(5) = -4.0d0 * sc(9)
c               gl(6) = sc(2)
c               gl(7) = 2.0d0 * (sc(7)-sc(8))
c               gl(8) = 2.0d0 * sc(10)

c
c     calculate the gl functions for induced components
c
               gli(1) = ck*sci(3) - ci*sci(4)
               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
               gli(6) = sci(1)
               gli(7) = 2.0d0 * (sci(7)-sci(8))
               glip(1) = ck*scip(3) - ci*scip(4)
               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
               glip(6) = scip(1)
               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c     compute the energy contributions for this interaction
c
               ei = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
     &                   + rr5*(gli(2)+gli(7))*psc5
     &                   + rr7*gli(3)*psc7)
               ei = f * ei
c               if (use_group) then
c                 ei = ei * fgrp
c               end if
               if (ii .eq. kk) then
                  ei = 0.5d0 * ei
               end if
               eptemp = eptemp + ei
c
c     increment the total intermolecular energy
c

c               if (molcule(ii) .ne. molcule(kk)) then
c                  einter = einter + e + ei
c               end if

c
c     intermediate variables for the permanent components
c

c               gf(1) = rr3*gl(0) + rr5*(gl(1)+gl(6))
c     &                    + rr7*(gl(2)+gl(7)+gl(8))
c     &                    + rr9*(gl(3)+gl(5)) + rr11*gl(4)
c               gf(2) = -ck*rr3 + sc(4)*rr5 - sc(6)*rr7
c               gf(3) = ci*rr3 + sc(3)*rr5 + sc(5)*rr7
c               gf(4) = 2.0d0 * rr5
c               gf(5) = 2.0d0 * (-ck*rr5+sc(4)*rr7-sc(6)*rr9)
c               gf(6) = 2.0d0 * (-ci*rr5-sc(3)*rr7-sc(5)*rr9)
c               gf(7) = 4.0d0 * rr7

c
c     intermediate variables for the induced components
c
               gfi(1) = 0.5d0 * rr5 * ((gli(1)+gli(6))*psc3
     &                                   + (glip(1)+glip(6))*dsc3
     &                                   + scip(2)*scale3i)
     &                + 0.5d0 * rr7 * ((gli(7)+gli(2))*psc5
     &                               + (glip(7)+glip(2))*dsc5
     &                      - (sci(3)*scip(4)+scip(3)*sci(4))*scale5i)
     &                + 0.5d0 * rr9 * (gli(3)*psc7+glip(3)*dsc7)
               gfi(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
               gfi(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
               gfi(4) = 2.0d0 * rr5
               gfi(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
               gfi(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)
c
c     get the permanent force components
c

c               ftm2(1) = gf(1)*xr + gf(2)*di(1) + gf(3)*dk(1)
c     &                      + gf(4)*(qkdi(1)-qidk(1)) + gf(5)*qir(1)
c     &                      + gf(6)*qkr(1) + gf(7)*(qiqkr(1)+qkqir(1))
c               ftm2(2) = gf(1)*yr + gf(2)*di(2) + gf(3)*dk(2)
c     &                      + gf(4)*(qkdi(2)-qidk(2)) + gf(5)*qir(2)
c     &                      + gf(6)*qkr(2) + gf(7)*(qiqkr(2)+qkqir(2))
c               ftm2(3) = gf(1)*zr + gf(2)*di(3) + gf(3)*dk(3)
c     &                      + gf(4)*(qkdi(3)-qidk(3)) + gf(5)*qir(3)
c     &                      + gf(6)*qkr(3) + gf(7)*(qiqkr(3)+qkqir(3))

c
c     get the induced force components
c
c
c     get the induced force components
c
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &           (- rr3*ck*(uind(1,l1)*psc3+uinp(1,l1)*dsc3)
     &            + rr5*sc(4)*(uind(1,l1)*psc5+uinp(1,l1)*dsc5)
     &            - rr7*sc(6)*(uind(1,l1)*psc7+uinp(1,l1)*dsc7))
     &            + (rr3*ci*(uind(1,l3)*psc3+uinp(1,l3)*dsc3)
     &            + rr5*sc(3)*(uind(1,l3)*psc5+uinp(1,l3)*dsc5)
     &            + rr7*sc(5)*(uind(1,l3)*psc7+uinp(1,l3)*dsc7))*0.5d0
     &            + rr5*scale5i*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &            + sci(3)*uinp(1,l3)+scip(3)*uind(1,l3))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
     &            + 0.5d0*gfi(4)*((qkui(1)-qiuk(1))*psc5
     &            + (qkuip(1)-qiukp(1))*dsc5)
     &            + gfi(5)*qir(1) + gfi(6)*qkr(1)
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &           (- rr3*ck*(uind(2,l1)*psc3+uinp(2,l1)*dsc3)
     &            + rr5*sc(4)*(uind(2,l1)*psc5+uinp(2,l1)*dsc5)
     &            - rr7*sc(6)*(uind(2,l1)*psc7+uinp(2,l1)*dsc7))
     &            + (rr3*ci*(uind(2,l3)*psc3+uinp(2,l3)*dsc3)
     &            + rr5*sc(3)*(uind(2,l3)*psc5+uinp(2,l3)*dsc5)
     &            + rr7*sc(5)*(uind(2,l3)*psc7+uinp(2,l3)*dsc7))*0.5d0
     &            + rr5*scale5i*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &            + sci(3)*uinp(2,l3)+scip(3)*uind(2,l3))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(2)
     &            + 0.5d0*gfi(4)*((qkui(2)-qiuk(2))*psc5
     &            + (qkuip(2)-qiukp(2))*dsc5)
     &            + gfi(5)*qir(2) + gfi(6)*qkr(2)
               ftm2i(3) = gfi(1)*zr  + 0.5d0*
     &           (- rr3*ck*(uind(3,l1)*psc3+uinp(3,l1)*dsc3)
     &            + rr5*sc(4)*(uind(3,l1)*psc5+uinp(3,l1)*dsc5)
     &            - rr7*sc(6)*(uind(3,l1)*psc7+uinp(3,l1)*dsc7))
     &            + (rr3*ci*(uind(3,l3)*psc3+uinp(3,l3)*dsc3)
     &            + rr5*sc(3)*(uind(3,l3)*psc5+uinp(3,l3)*dsc5)
     &            + rr7*sc(5)*(uind(3,l3)*psc7+uinp(3,l3)*dsc7))*0.5d0
     &            + rr5*scale5i*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &            + sci(3)*uinp(3,l3)+scip(3)*uind(3,l3))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
     &            + 0.5d0*gfi(4)*((qkui(3)-qiuk(3))*psc5
     &            + (qkuip(3)-qiukp(3))*dsc5)
     &            + gfi(5)*qir(3) + gfi(6)*qkr(3)
c
c     account for partially excluded induced interactions
c
               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(kk)
     &                                  +(glip(1)+glip(6))*dscale(kk))
               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(kk)
     &                                  +(glip(2)+glip(7))*dscale(kk))
               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(kk)
     &                                  +glip(3)*dscale(kk))
               fridmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
     &                        + temp7*ddsc7(1)
               fridmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
     &                        + temp7*ddsc7(2)
               fridmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
     &                        + temp7*ddsc7(3)
c
c     find some scaling terms for induced-induced force
c
               temp3 = 0.5d0 * rr3 * uscale(kk) * scip(2)
               temp5 = -0.5d0 * rr5 * uscale(kk)
     &                    * (sci(3)*scip(4)+scip(3)*sci(4))
               findmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
               findmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
               findmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
c
c     modify induced force for partially excluded interactions
c
               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)
c
c     correction to convert mutual to direct polarization force
c

c               if (poltyp .eq. 'DIRECT') then
c                  gfd = 0.5d0 * (rr5*scip(2)*scale3i
c     &                  - rr7*(scip(3)*sci(4)+sci(3)*scip(4))*scale5i)
c                  temp5 = 0.5d0 * rr5 * scale5i
c                  fdir(1) = gfd*xr + temp5
c     &                         * (sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
c     &                           +sci(3)*uinp(1,l3)+scip(3)*uind(1,l3))
c                  fdir(2) = gfd*yr + temp5
c     &                         * (sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
c     &                           +sci(3)*uinp(2,l3)+scip(3)*uind(2,l3))
c                  fdir(3) = gfd*zr + temp5
c     &                         * (sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
c     &                           +sci(3)*uinp(3,l3)+scip(3)*uind(3,l3))
c                  ftm2i(1) = ftm2i(1) - fdir(1) + findmp(1)
c                  ftm2i(2) = ftm2i(2) - fdir(2) + findmp(2)
c                  ftm2i(3) = ftm2i(3) - fdir(3) + findmp(3)
c               end if

c
c     intermediate terms for induced torque on multipoles
c
               gti(2) = 0.5d0 * rr5 * (sci(4)*psc5+scip(4)*dsc5)
               gti(3) = 0.5d0 * rr5 * (sci(3)*psc5+scip(3)*dsc5)
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
c
c     get the permanent torque components
c

c               ttm2(1) = -rr3*dixdk(1) + gf(2)*dixr(1) - gf(5)*rxqir(1)
c     &           + gf(4)*(dixqkr(1)+dkxqir(1)+rxqidk(1)-2.0d0*qixqk(1))
c     &           - gf(7)*(rxqikr(1)+qkrxqir(1))
c               ttm2(2) = -rr3*dixdk(2) + gf(2)*dixr(2) - gf(5)*rxqir(2)
c     &           + gf(4)*(dixqkr(2)+dkxqir(2)+rxqidk(2)-2.0d0*qixqk(2))
c     &           - gf(7)*(rxqikr(2)+qkrxqir(2))
c               ttm2(3) = -rr3*dixdk(3) + gf(2)*dixr(3) - gf(5)*rxqir(3)
c     &           + gf(4)*(dixqkr(3)+dkxqir(3)+rxqidk(3)-2.0d0*qixqk(3))
c     &           - gf(7)*(rxqikr(3)+qkrxqir(3))
c               ttm3(1) = rr3*dixdk(1) + gf(3)*dkxr(1) - gf(6)*rxqkr(1)
c     &           - gf(4)*(dixqkr(1)+dkxqir(1)+rxqkdi(1)-2.0d0*qixqk(1))
c     &           - gf(7)*(rxqkir(1)-qkrxqir(1))
c               ttm3(2) = rr3*dixdk(2) + gf(3)*dkxr(2) - gf(6)*rxqkr(2)
c     &           - gf(4)*(dixqkr(2)+dkxqir(2)+rxqkdi(2)-2.0d0*qixqk(2))
c     &           - gf(7)*(rxqkir(2)-qkrxqir(2))
c               ttm3(3) = rr3*dixdk(3) + gf(3)*dkxr(3) - gf(6)*rxqkr(3)
c     &           - gf(4)*(dixqkr(3)+dkxqir(3)+rxqkdi(3)-2.0d0*qixqk(3))
c     &           - gf(7)*(rxqkir(3)-qkrxqir(3))

c
c     get the induced torque components
c
               ttm2i(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
     &           + gti(2)*dixr(1) + gti(4)*((ukxqir(1)+rxqiuk(1))*psc5
     &           +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
     &           + gti(2)*dixr(2) + gti(4)*((ukxqir(2)+rxqiuk(2))*psc5
     &           +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
     &           + gti(2)*dixr(3) + gti(4)*((ukxqir(3)+rxqiuk(3))*psc5
     &           +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gti(5)*rxqir(3)
               ttm3i(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
     &           + gti(3)*dkxr(1) - gti(4)*((uixqkr(1)+rxqkui(1))*psc5
     &           +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0 - gti(6)*rxqkr(1)
               ttm3i(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
     &           + gti(3)*dkxr(2) - gti(4)*((uixqkr(2)+rxqkui(2))*psc5
     &           +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0 - gti(6)*rxqkr(2)
               ttm3i(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
     &           + gti(3)*dkxr(3) - gti(4)*((uixqkr(3)+rxqkui(3))*psc5
     &           +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0 - gti(6)*rxqkr(3)
c
c     handle the case where scaling is used
c
               do j = 1, 3
c                  ftm2(j) = f * ftm2(j)
                  ftm2i(j) = f * ftm2i(j)
c                  ttm2(j) = f * ttm2(j)
                  ttm2i(j) = f * ttm2i(j)
c                  ttm3(j) = f * ttm3(j)
                  ttm3i(j) = f * ttm3i(j)
               end do
c               if (use_polymer) then
c                  if (r2 .le. polycut2) then
c                     do j = 1, 3
c                        ftm2(j) = ftm2(j) * mscale(kk)
c                        ttm2(j) = ttm2(j) * mscale(kk)
c                        ttm3(j) = ttm3(j) * mscale(kk)
c                     end do
c                  end if
c               end if
c               if (use_group) then
c                  do j = 1, 3
c                     ftm2(j) = ftm2(j) * fgrp
c                     ttm2(j) = ttm2(j) * fgrp
c                     ttm3(j) = ttm3(j) * fgrp
c                    ftm2i(j) = ftm2i(j) * fgrp
c                    ttm2i(j) = ttm2i(j) * fgrp
c                    ttm3i(j) = ttm3i(j) * fgrp
c                  end do
c               end if
               if (ii .eq. kk) then
                  do j = 1, 3
c                     ftm2(j) = 0.5d0 * ftm2(j)
                     ftm2i(j) = 0.5d0 * ftm2i(j)
c                     ttm2(j) = 0.5d0 * ttm2(j)
                     ttm2i(j) = 0.5d0 * ttm2i(j)
c                     ttm3(j) = 0.5d0 * ttm3(j)
                     ttm3i(j) = 0.5d0 * ttm3i(j)
                  end do
               end if
c
c     increment gradient due to force and torque on first site
c

c               dem(1,ii) = dem(1,ii) + ftm2(1)
c               dem(2,ii) = dem(2,ii) + ftm2(2)
c               dem(3,ii) = dem(3,ii) + ftm2(3)
c               dep(1,ii) = dep(1,ii) + ftm2i(1)
c               dep(2,ii) = dep(2,ii) + ftm2i(2)
c               dep(3,ii) = dep(3,ii) + ftm2i(3)
               deptemp(1,i) = deptemp(1,i) + ftm2i(1)
               deptemp(2,i) = deptemp(2,i) + ftm2i(2)
               deptemp(3,i) = deptemp(3,i) + ftm2i(3)
c               call torque (i,ttm2,ttm2i,frcxi,frcyi,frczi)
               call torque_3b (deptemp,i,
     &          ttm2,ttm2i,frcxi,frcyi,frczi)
c
c     increment gradient due to force and torque on second site
c

c               dem(1,kk) = dem(1,kk) - ftm2(1)
c               dem(2,kk) = dem(2,kk) - ftm2(2)
c               dem(3,kk) = dem(3,kk) - ftm2(3)
c               dep(1,kk) = dep(1,kk) - ftm2i(1)
c               dep(2,kk) = dep(2,kk) - ftm2i(2)
c               dep(3,kk) = dep(3,kk) - ftm2i(3)
               deptemp(1,k) = deptemp(1,k) - ftm2i(1)
               deptemp(2,k) = deptemp(2,k) - ftm2i(2)
               deptemp(3,k) = deptemp(3,k) - ftm2i(3)
c               call torque (k,ttm3,ttm3i,frcxk,frcyk,frczk)
               call torque_3b (deptemp,k,
     &            ttm3,ttm3i,frcxk,frcyk,frczk)

c
c     increment the internal virial tensor components
c
               iaz = iz
               iax = ix
               iay = iy
               kaz = kz
               kax = kx
               kay = ky
               if (iaz .eq. 0)  iaz = ii
               if (iax .eq. 0)  iax = ii
               if (iay .eq. 0)  iay = ii
               if (kaz .eq. 0)  kaz = kk
               if (kax .eq. 0)  kax = kk
               if (kay .eq. 0)  kay = kk
               xiz = x(iaz) - x(ii)
               yiz = y(iaz) - y(ii)
               ziz = z(iaz) - z(ii)
               xix = x(iax) - x(ii)
               yix = y(iax) - y(ii)
               zix = z(iax) - z(ii)
               xiy = x(iay) - x(ii)
               yiy = y(iay) - y(ii)
               ziy = z(iay) - z(ii)
               xkz = x(kaz) - x(kk)
               ykz = y(kaz) - y(kk)
               zkz = z(kaz) - z(kk)
               xkx = x(kax) - x(kk)
               ykx = y(kax) - y(kk)
               zkx = z(kax) - z(kk)
               xky = x(kay) - x(kk)
               yky = y(kay) - y(kk)
               zky = z(kay) - z(kk)
c               vxx = -xr*(ftm2(1)+ftm2i(1)) + xix*frcxi(1)
c     &                  + xiy*frcyi(1) + xiz*frczi(1) + xkx*frcxk(1)
c     &                  + xky*frcyk(1) + xkz*frczk(1)
c               vyx = -yr*(ftm2(1)+ftm2i(1)) + yix*frcxi(1)
c     &                  + yiy*frcyi(1) + yiz*frczi(1) + ykx*frcxk(1)
c     &                  + yky*frcyk(1) + ykz*frczk(1)
c               vzx = -zr*(ftm2(1)+ftm2i(1)) + zix*frcxi(1)
c     &                  + ziy*frcyi(1) + ziz*frczi(1) + zkx*frcxk(1)
c     &                  + zky*frcyk(1) + zkz*frczk(1)
c               vyy = -yr*(ftm2(2)+ftm2i(2)) + yix*frcxi(2)
c     &                  + yiy*frcyi(2) + yiz*frczi(2) + ykx*frcxk(2)
c     &                  + yky*frcyk(2) + ykz*frczk(2)
c               vzy = -zr*(ftm2(2)+ftm2i(2)) + zix*frcxi(2)
c     &                  + ziy*frcyi(2) + ziz*frczi(2) + zkx*frcxk(2)
c     &                  + zky*frcyk(2) + zkz*frczk(2)
c               vzz = -zr*(ftm2(3)+ftm2i(3)) + zix*frcxi(3)
c     &                  + ziy*frcyi(3) + ziz*frczi(3) + zkx*frcxk(3)
c     &                  + zky*frcyk(3) + zkz*frczk(3)
c               vir(1,1) = vir(1,1) + vxx
c               vir(2,1) = vir(2,1) + vyx
c               vir(3,1) = vir(3,1) + vzx
c               vir(1,2) = vir(1,2) + vyx
c               vir(2,2) = vir(2,2) + vyy
c               vir(3,2) = vir(3,2) + vzy
c               vir(1,3) = vir(1,3) + vzx
c               vir(2,3) = vir(2,3) + vzy
c               vir(3,3) = vir(3,3) + vzz
c            end if
               vxx = -xr*(ftm2i(1)) + xix*frcxi(1)
     &                  + xiy*frcyi(1) + xiz*frczi(1) + xkx*frcxk(1)
     &                  + xky*frcyk(1) + xkz*frczk(1)
               vyx = -yr*(ftm2i(1)) + yix*frcxi(1)
     &                  + yiy*frcyi(1) + yiz*frczi(1) + ykx*frcxk(1)
     &                  + yky*frcyk(1) + ykz*frczk(1)
               vzx = -zr*(ftm2i(1)) + zix*frcxi(1)
     &                  + ziy*frcyi(1) + ziz*frczi(1) + zkx*frcxk(1)
     &                  + zky*frcyk(1) + zkz*frczk(1)
               vyy = -yr*(ftm2i(2)) + yix*frcxi(2)
     &                  + yiy*frcyi(2) + yiz*frczi(2) + ykx*frcxk(2)
     &                  + yky*frcyk(2) + ykz*frczk(2)
               vzy = -zr*(ftm2i(2)) + zix*frcxi(2)
     &                  + ziy*frcyi(2) + ziz*frczi(2) + zkx*frcxk(2)
     &                  + zky*frcyk(2) + zkz*frczk(2)
               vzz = -zr*(ftm2i(3)) + zix*frcxi(3)
     &                  + ziy*frcyi(3) + ziz*frczi(3) + zkx*frcxk(3)
     &                  + zky*frcyk(3) + zkz*frczk(3)

               virtemp(1,1) = virtemp(1,1) + vxx
               virtemp(2,1) = virtemp(2,1) + vyx
               virtemp(3,1) = virtemp(3,1) + vzx
               virtemp(1,2) = virtemp(1,2) + vyx
               virtemp(2,2) = virtemp(2,2) + vyy
               virtemp(3,2) = virtemp(3,2) + vzy
               virtemp(1,3) = virtemp(1,3) + vzx
               virtemp(2,3) = virtemp(2,3) + vzy
               virtemp(3,3) = virtemp(3,3) + vzz

            end do
   20       continue
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
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
            uscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
            uscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
            uscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
            uscale(ip14(j,ii)) = 1.0d0
         end do
      end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (mscale)
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      return
      end

