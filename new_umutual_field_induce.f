      subroutine udirect1_umutual_ewrp_3b(aewald3b,ewaldcut3b,
     & jmax3b,kmax3b,lmax3b,field,M_recip,npole3b,pnum)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'ewald.i'
c      include 'ewreg.i'
      include 'ewreg2.i'
      include 'math.i'
      include 'mpole.i'
c      include 'polar.i'
      include 'polar2.i'
      include 'units.i'
c      include 'combo.i'
      integer i,j,k,l,ii
      integer jmin,jmax
      integer kmin,kmax
      integer lmin,lmax
      integer l1,l3,i1,k1,i2,j2
      real*8 expterm,cut
      real*8 term,fterm
      real*8 xfr,yfr,zfr
      real*8 rj,rk,rl
      real*8 h1,h2,h3,hsq
      real*8 qf,t1,t2
      real*8 ck,dk,qk
      real*8 q1,q2,q3
      integer npole3b,pnum(*)
      real*8 field(3,*)
      real*8 cm(npole3b)
      real*8 dm(3,npole3b)
      real*8 qm(9,npole3b)
      real*8 frecip
      integer a,a1,a2,a3
      integer b,b1,b2,b3
      real*8 T1xx,T1xy,T1xz,T1yx,T1yy,T1yz,T1zx,T1zy,T1zz
      real*8 T2xx,T2xy,T2xz,T2yx,T2yy,T2yz,T2zx,T2zy,T2zz
      real*8 T1_rec(3*npole3b,3*npole3b)
      real*8 T2_rec(3*npole3b,3*npole3b)
      real*8 M_recip(3*npole3b,3*npole3b)
      real timediff,timearray(2)
      real*8 ckr(npole3b)
      real*8 skr(npole3b)
      real*8 cjk(npole3b)
      real*8 sjk(npole3b)
      real*8 ejc(npole3b,0:maxvec)
      real*8 ejs(npole3b,0:maxvec)
      real*8 ekc(npole3b,-maxvec:maxvec)
      real*8 eks(npole3b,-maxvec:maxvec)
      real*8 elc(npole3b,-maxvec:maxvec)
      real*8 els(npole3b,-maxvec:maxvec)
      real*8 aewald,aewald3b,ewaldcut3b
      integer jmax3b,kmax3b,lmax3b

      frecip = 0.5d0
      do l1 = 1, 3*npole3b
         do l3 = 1, 3*npole3b
            M_recip(l1,l3) = 0.0d0
            T1_rec(l1,l3) = 0.0d0
            T2_rec(l1,l3) = 0.0d0
         end do
      end do


      aewald=aewald3b

      if (aewald .lt. 1.0d-6)  return
      term = -0.25d0 / aewald**2
      fterm = 8.0d0 * pi / volbox
c
c     set the number of vectors based on box dimensions
c
      frecip = 0.5d0
      cut = 4.0d0 * pi * pi * frecip
      jmin = 0
      kmin = 0
      lmin = 1
c      jmax = min(maxvec,int(frecip/recip(1,1)))
c      kmax = min(maxvec,int(frecip/recip(2,2)))
c      lmax = min(maxvec,int(frecip/recip(3,3)))
      jmax=min(maxvec,jmax3b)
      kmax=min(maxvec,kmax3b)
      lmax=min(maxvec,lmax3b)
c
c     copy the multipole moments into local storage areas
c
      do l1 = 1, npole3b
         i = pnum(l1)
         cm(l1) = rpole(1,i)
         dm(1,l1) = rpole(2,i)
         dm(2,l1) = rpole(3,i)
         dm(3,l1) = rpole(4,i)
         qm(1,l1) = rpole(5,i)
         qm(2,l1) = rpole(6,i)
         qm(3,l1) = rpole(7,i)
         qm(4,l1) = rpole(8,i)
         qm(5,l1) = rpole(9,i)
         qm(6,l1) = rpole(10,i)
         qm(7,l1) = rpole(11,i)
         qm(8,l1) = rpole(12,i)
         qm(9,l1) = rpole(13,i)
      end do

      do l1 = 1, npole3b
         i = pnum(l1)
         ii = ipole(i)
         zfr = (z(ii)/gamma_term) / zbox
         yfr = ((y(ii)-zfr*zbox*beta_term)/gamma_sin) / ybox
         xfr = (x(ii)-yfr*ybox*gamma_cos-zfr*zbox*beta_cos) / xbox
         xfr = 2.0d0 * pi * xfr
         yfr = 2.0d0 * pi * yfr
         zfr = 2.0d0 * pi * zfr
         ejc(l1,0) = 1.0d0
         ejs(l1,0) = 0.0d0
         ekc(l1,0) = 1.0d0
         eks(l1,0) = 0.0d0
         elc(l1,0) = 1.0d0
         els(l1,0) = 0.0d0
         ejc(l1,1) = cos(xfr)
         ejs(l1,1) = sin(xfr)
         ekc(l1,1) = cos(yfr)
         eks(l1,1) = sin(yfr)
         elc(l1,1) = cos(zfr)
         els(l1,1) = sin(zfr)
         ekc(l1,-1) = ekc(l1,1)
         eks(l1,-1) = -eks(l1,1)
         elc(l1,-1) = elc(l1,1)
         els(l1,-1) = -els(l1,1)

         do j = 2, jmax
            ejc(l1,j) = ejc(l1,j-1)*ejc(l1,1) - ejs(l1,j-1)*ejs(l1,1)
            ejs(l1,j) = ejs(l1,j-1)*ejc(l1,1) + ejc(l1,j-1)*ejs(l1,1)
         end do
         do j = 2, kmax
            ekc(l1,j) = ekc(l1,j-1)*ekc(l1,1) - eks(l1,j-1)*eks(l1,1)
            eks(l1,j) = eks(l1,j-1)*ekc(l1,1) + ekc(l1,j-1)*eks(l1,1)
            ekc(l1,-j) = ekc(l1,j)
            eks(l1,-j) = -eks(l1,j)
         end do
         do j = 2, lmax
            elc(l1,j) = elc(l1,j-1)*elc(l1,1) - els(l1,j-1)*els(l1,1)
            els(l1,j) = els(l1,j-1)*elc(l1,1) + elc(l1,j-1)*els(l1,1)
            elc(l1,-j) = elc(l1,j)
            els(l1,-j) = -els(l1,j)
         end do
      end do

      do j = jmin, jmax
         rj = 2.0d0 * pi * dble(j)
         do k = kmin, kmax
            rk = 2.0d0 * pi * dble(k)
            do l1 = 1, npole3b
               i = pnum(l1)
c               cjk(i) = ejc(i,j)*ekc(i,k) - ejs(i,j)*eks(i,k)
c               sjk(i) = ejs(i,j)*ekc(i,k) + ejc(i,j)*eks(i,k)
               cjk(l1) = ejc(l1,j)*ekc(l1,k) - ejs(l1,j)*eks(l1,k)
               sjk(l1) = ejs(l1,j)*ekc(l1,k) + ejc(l1,j)*eks(l1,k)
            end do
            do l = lmin, lmax
               rl = 2.0d0 * pi * dble(l)
               h1 = recip(1,1)*rj
               h2 = recip(2,1)*rj + recip(2,2)*rk
               h3 = recip(3,1)*rj + recip(3,2)*rk + recip(3,3)*rl
               hsq = h1*h1 + h2*h2 + h3*h3
c               if (hsq .le. cut) then
                  t1 = 0.0d0
                  t2 = 0.0d0
                  do l1 = 1, npole3b
                     i = pnum(l1)
c                     ckr(i) = cjk(i)*elc(i,l) - sjk(i)*els(i,l)
c                     skr(i) = sjk(i)*elc(i,l) + cjk(i)*els(i,l)
                     ckr(l1) = cjk(l1)*elc(l1,l) - sjk(l1)*els(l1,l)
                     skr(l1) = sjk(l1)*elc(l1,l) + cjk(l1)*els(l1,l)
                     ck = cm(l1)
                     dk = h1*dm(1,l1) + h2*dm(2,l1) + h3*dm(3,l1)
                     q1 = h1*qm(1,l1) + h2*qm(4,l1) + h3*qm(7,l1)
                     q2 = h1*qm(2,l1) + h2*qm(5,l1) + h3*qm(8,l1)
                     q3 = h1*qm(3,l1) + h2*qm(6,l1) + h3*qm(9,l1)

                     qk = h1*q1 + h2*q2 + h3*q3
c                     t1 = t1 + (ck-qk)*skr(i) + dk*ckr(i)
c                     t2 = t2 + (ck-qk)*ckr(i) - dk*skr(i)
                     t1 = t1 + (ck-qk)*skr(l1) + dk*ckr(l1)
                     t2 = t2 + (ck-qk)*ckr(l1) - dk*skr(l1)

                  end do
                  expterm = fterm * exp(term*hsq) / hsq
                  if (octahedron) then
                     if (mod(j+k+l,2) .ne. 0)  expterm = 0.0d0
                  end if
                  do l1 = 1, npole3b
                     a1 = 3*(l1-1)
                     qf = expterm * (skr(l1)*t2-ckr(l1)*t1)
                     field(1,l1) = field(1,l1) + h1*qf
                     field(2,l1) = field(2,l1) + h2*qf
                     field(3,l1) = field(3,l1) + h3*qf
                     do l3 = l1, npole3b
                        b1 = 3*(l3-1)
                        T1xx = ckr(l1)*ckr(l3)*h1*h1
                        T1xy = ckr(l1)*ckr(l3)*h1*h2
                        T1xz = ckr(l1)*ckr(l3)*h1*h3
                        T1yy = ckr(l1)*ckr(l3)*h2*h2
                        T1yz = ckr(l1)*ckr(l3)*h2*h3
                        T1zz = ckr(l1)*ckr(l3)*h3*h3

                        T1yx = T1xy
                        T1zx = T1xz
                        T1zy = T1yz

                        T2xx = -skr(l1)*skr(l3)*h1*h1
                        T2xy = -skr(l1)*skr(l3)*h1*h2
                        T2xz = -skr(l1)*skr(l3)*h1*h3
                        T2yy = -skr(l1)*skr(l3)*h2*h2
                        T2yz = -skr(l1)*skr(l3)*h2*h3
                        T2zz = -skr(l1)*skr(l3)*h3*h3

                        T2yx = T2xy
                        T2zx = T2xz
                        T2zy = T2yz

                        M_recip(a1+1,b1+1) =  M_recip(a1+1,b1+1) +
     &                   expterm*(T2xx - T1xx)
                        M_recip(a1+1,b1+2) =  M_recip(a1+1,b1+2) +
     &                   expterm*(T2xy - T1xy)
                        M_recip(a1+1,b1+3) =  M_recip(a1+1,b1+3) +
     &                   expterm*(T2xz - T1xz)
                        M_recip(a1+2,b1+1) =  M_recip(a1+2,b1+1) +
     &                   expterm*(T2yx - T1yx)
                        M_recip(a1+2,b1+2) =  M_recip(a1+2,b1+2) +
     &                   expterm*(T2yy - T1yy)
                        M_recip(a1+2,b1+3) =  M_recip(a1+2,b1+3) +
     &                   expterm*(T2yz - T1yz)
                        M_recip(a1+3,b1+1) =  M_recip(a1+3,b1+1) +
     &                   expterm*(T2zx - T1zx)
                        M_recip(a1+3,b1+2) =  M_recip(a1+3,b1+2) +
     &                   expterm*(T2zy - T1zy)
                        M_recip(a1+3,b1+3) =  M_recip(a1+3,b1+3) +
     &                   expterm*(T2zz - T1zz)

                     end do
                  end do
c               end if
            end do
            lmin = -lmax
         end do
         kmin = -kmax
      end do
      do l1 = 1 , 3*npole3b-1
         do l3 = l1+1, 3*npole3b
            M_recip(l3,l1) = M_recip(l1,l3)
         end do
      end do
      return
      end


      subroutine udirect2a_umutual_ewrl_3b(aewald3b,ewaldcut3b,
     &  field,fieldp,M_real,npole3b,pnum)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'ewald.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar2.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'units.i'
      integer i,j,k,m
c      integer ii,kk,l1,l2
      integer ii,kk
      real*8 xr,yr,zr
      real*8 r,r2,rr3,rr5
      real*8 erfc,bfac,exp2a
      real*8 drr3,drr5,drr7
      real*8 prr3,prr5,prr7
      real*8 ci,dix,diy,diz
      real*8 qixx,qiyy,qizz
      real*8 qixy,qixz,qiyz
      real*8 ck,dkx,dky,dkz
      real*8 qkxx,qkyy,qkzz
      real*8 qkxy,qkxz,qkyz
      real*8 dir,dkr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 ralpha
      real*8 alsq2,alsq2n
      real*8 pdi,pti,pgamma
      real*8 damp,expdamp
      real*8 scale3_dir,scale5_dir
      real*8 scale7_dir
      real*8 scale3_mut,scale5_mut
      real*8 dsc3,dsc5,dsc7
      real*8 psc3,psc5,psc7
      real*8 bn(0:3)
      real*8 fimd(3),fkmd(3)
      real*8 fimp(3),fkmp(3)
      real*8 fim(3),fkm(3)
      real*8 fid(3),fkd(3)
      real*8 fip(3),fkp(3)
      real*8, allocatable :: dscale_dir(:)
      real*8, allocatable :: dscale_mut(:)
      real*8, allocatable :: pscale(:)
      real*8 field(3,*)
      real*8 fieldp(3,*)
      integer npole3b,pnum(*)
      character*6 mode
      external erfc
      integer l1,l2,l3,i2,k2
      real*8 M_real(3*npole3b,3*npole3b)
      real*8 Txx,Txy,Txz,Tyx
      real*8 Tyy,Tyz,Tzx,Tzy,Tzz
      real*8 fid_xx,fid_xy,fid_xz
      real*8 fid_yx,fid_yy,fid_yz
      real*8 fid_zx,fid_zy,fid_zz
      real*8 fimd_xx,fimd_xy,fimd_xz
      real*8 fimd_yx,fimd_yy,fimd_yz
      real*8 fimd_zx,fimd_zy,fimd_zz
      real*8 Mk_x,Mk_y,Mk_z,fkx,fky,fkz
      real*8 Mk0_x,Mk0_y,Mk0_z
      real*8 aewald,aewald3b,ewaldcut3b,ewaldcut3b_2
c
c
c     check for multipoles and set cutoff coefficients
c
      aewald=aewald3b
      if (npole3b .eq. 0)  return
c      mode = 'EWALD'
c      call switch (mode)

c      mode = 'EWALD3B'
c      call switch (mode)
      ewaldcut3b_2=ewaldcut3b*ewaldcut3b
c
c     perform dynamic allocation of some local arrays
c
      allocate (dscale_dir(n))
      allocate (dscale_mut(n))
      allocate (pscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      do i = 1, n
         pscale(i) = 1.0d0
         dscale_dir(i) = 1.0d0
         dscale_mut(i) = 1.0d0
      end do
c
c     compute the real space portion of the Ewald summation
c
      do l1 = 1, npole3b-1
         i = pnum(l1)
         i2 = 3*(l1-1)
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
         do j = 1, np11(ii)
            dscale_mut(ip11(j,ii)) = u1scale
            dscale_dir(ip11(j,ii)) = d1scale
         end do
         do j = 1, np12(ii)
            dscale_mut(ip12(j,ii)) = u2scale
            dscale_dir(ip12(j,ii)) = d2scale
         end do
         do j = 1, np13(ii)
            dscale_mut(ip13(j,ii)) = u3scale
            dscale_dir(ip13(j,ii)) = d3scale
         end do
         do j = 1, np14(ii)
            dscale_mut(ip14(j,ii)) = u4scale
            dscale_dir(ip14(j,ii)) = d4scale
         end do

         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do l3 = l1+1, npole3b
            k = pnum(l3)
            k2 = 3*(l3-1)
            kk = ipole(k)
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr* yr + zr*zr
c            if (r2 .le. cut2) then
            if (r2 .le. ewaldcut3b_2) then

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
               do j = 1, 3
                  bfac = dble(j+j-1)
                  alsq2n = alsq2 * alsq2n
                  bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
               end do
c
c     compute the error function scaled and unscaled terms
c

               scale3_dir = 1.0d0
               scale5_dir = 1.0d0
               scale7_dir = 1.0d0
               scale3_mut = dscale_mut(kk)
               scale5_mut = dscale_mut(kk)
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3_dir = 1.0d0 - expdamp
                     scale5_dir = 1.0d0 - expdamp*(1.0d0-damp)
                     scale7_dir = 1.0d0 - expdamp
     &                           *(1.0d0-damp+0.6d0*damp**2)

                     scale3_mut=scale3_mut * (1.0d0-expdamp)
                     scale5_mut=scale5_mut*(1.0d0-(1.0d0-damp)*expdamp)

                  end if
               end if
               rr3 = (1.0d0-scale3_mut) / (r*r2)
               rr5 = 3.0d0 * (1.0d0-scale5_mut) / (r*r2*r2)

               dsc3 = scale3_dir * dscale_dir(kk)
               dsc5 = scale5_dir * dscale_dir(kk)
               dsc7 = scale7_dir * dscale_dir(kk)
               psc3 = scale3_dir * pscale(kk)
               psc5 = scale5_dir * pscale(kk)
               psc7 = scale7_dir * pscale(kk)
               drr3 = (1.0d0-dsc3) / (r*r2)
               drr5 = 3.0d0 * (1.0d0-dsc5) / (r*r2*r2)
               drr7 = 15.0d0 * (1.0d0-dsc7) / (r*r2*r2*r2)
               prr3 = (1.0d0-psc3) / (r*r2)
               prr5 = 3.0d0 * (1.0d0-psc5) / (r*r2*r2)
               prr7 = 15.0d0 * (1.0d0-psc7) / (r*r2*r2*r2)
               dir = dix*xr + diy*yr + diz*zr
               qix = qixx*xr + qixy*yr + qixz*zr
               qiy = qixy*xr + qiyy*yr + qiyz*zr
               qiz = qixz*xr + qiyz*yr + qizz*zr
               qir = qix*xr + qiy*yr + qiz*zr
               dkr = dkx*xr + dky*yr + dkz*zr
               qkx = qkxx*xr + qkxy*yr + qkxz*zr
               qky = qkxy*xr + qkyy*yr + qkyz*zr
               qkz = qkxz*xr + qkyz*yr + qkzz*zr
               qkr = qkx*xr + qky*yr + qkz*zr
               fim(1) = -xr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                     - bn(1)*dkx + 2.0d0*bn(2)*qkx
               fim(2) = -yr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                     - bn(1)*dky + 2.0d0*bn(2)*qky
               fim(3) = -zr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                     - bn(1)*dkz + 2.0d0*bn(2)*qkz
               fkm(1) = xr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                     - bn(1)*dix - 2.0d0*bn(2)*qix
               fkm(2) = yr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                     - bn(1)*diy - 2.0d0*bn(2)*qiy
               fkm(3) = zr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                     - bn(1)*diz - 2.0d0*bn(2)*qiz
               fid(1) = -xr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                     - drr3*dkx + 2.0d0*drr5*qkx
               fid(2) = -yr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                     - drr3*dky + 2.0d0*drr5*qky
               fid(3) = -zr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                     - drr3*dkz + 2.0d0*drr5*qkz
               fkd(1) = xr*(drr3*ci+drr5*dir+drr7*qir)
     &                     - drr3*dix - 2.0d0*drr5*qix
               fkd(2) = yr*(drr3*ci+drr5*dir+drr7*qir)
     &                     - drr3*diy - 2.0d0*drr5*qiy
               fkd(3) = zr*(drr3*ci+drr5*dir+drr7*qir)
     &                     - drr3*diz - 2.0d0*drr5*qiz
               fip(1) = -xr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                     - prr3*dkx + 2.0d0*prr5*qkx
               fip(2) = -yr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                     - prr3*dky + 2.0d0*prr5*qky
               fip(3) = -zr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                     - prr3*dkz + 2.0d0*prr5*qkz
               fkp(1) = xr*(prr3*ci+prr5*dir+prr7*qir)
     &                     - prr3*dix - 2.0d0*prr5*qix
               fkp(2) = yr*(prr3*ci+prr5*dir+prr7*qir)
     &                     - prr3*diy - 2.0d0*prr5*qiy
               fkp(3) = zr*(prr3*ci+prr5*dir+prr7*qir)
     &                     - prr3*diz - 2.0d0*prr5*qiz

               fid_xx = -rr3 + xr*xr*rr5
               fid_xy =  xr*yr*rr5
               fid_xz =  xr*zr*rr5
               fid_yx =  yr*xr*rr5
               fid_yy = -rr3 + yr*yr*rr5
               fid_yz =  yr*zr*rr5
               fid_zx =  zr*xr*rr5
               fid_zy =  zr*yr*rr5
               fid_zz = -rr3 + zr*zr*rr5

               fimd_xx = -bn(1) + xr*xr*bn(2)
               fimd_xy =  xr*yr*bn(2)
               fimd_xz =  xr*zr*bn(2)
               fimd_yx =  yr*xr*bn(2)
               fimd_yy = -bn(1) + yr*yr*bn(2)
               fimd_yz =  yr*zr*bn(2)
               fimd_zx =  zr*xr*bn(2)
               fimd_zy =  zr*yr*bn(2)
               fimd_zz = -bn(1) + zr*zr*bn(2)

               Txx = (fimd_xx - fid_xx)
               Txy = (fimd_xy - fid_xy)
               Txz = (fimd_xz - fid_xz)
               Tyx = (fimd_yx - fid_yx)
               Tyy = (fimd_yy - fid_yy)
               Tyz = (fimd_yz - fid_yz)
               Tzx = (fimd_zx - fid_zx)
               Tzy = (fimd_zy - fid_zy)
               Tzz = (fimd_zz - fid_zz)

               M_real(i2+1,k2+1) = Txx
               M_real(i2+1,k2+2) = Txy
               M_real(i2+1,k2+3) = Txz
               M_real(i2+2,k2+1) = Tyx
               M_real(i2+2,k2+2) = Tyy
               M_real(i2+2,k2+3) = Tyz
               M_real(i2+3,k2+1) = Tzx
               M_real(i2+3,k2+2) = Tzy
               M_real(i2+3,k2+3) = Tzz

               M_real(k2+1,i2+1) = Txx
               M_real(k2+1,i2+2) = Txy
               M_real(k2+1,i2+3) = Txz
               M_real(k2+2,i2+1) = Tyx
               M_real(k2+2,i2+2) = Tyy
               M_real(k2+2,i2+3) = Tyz
               M_real(k2+3,i2+1) = Tzx
               M_real(k2+3,i2+2) = Tzy
               M_real(k2+3,i2+3) = Tzz

c
c     increment the field at each site due to this interaction
c
               do j = 1, 3
                  field(j,l1) = field(j,l1) + fim(j) - fid(j)
                  field(j,l3) = field(j,l3) + fkm(j) - fkd(j)
                  fieldp(j,l1) = fieldp(j,l1) + fim(j) - fip(j)
                  fieldp(j,l3) = fieldp(j,l3) + fkm(j) - fkp(j)
               end do
            end if
         end do

         do j = 1, n12(ii)
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale_dir(ip11(j,ii)) = 1.0d0
            dscale_mut(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale_dir(ip12(j,ii)) = 1.0d0
            dscale_mut(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale_dir(ip13(j,ii)) = 1.0d0
            dscale_mut(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale_dir(ip14(j,ii)) = 1.0d0
            dscale_mut(ip14(j,ii)) = 1.0d0
         end do
      end do

c
c     periodic boundary for lareg cutoffs via replicates method
c
      if (use_replica) then
        do l1 = 1, npole3b
          i = pnum(l1)
          i2 = 3*(l1-1)
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
            do j = 1, n12(ii)
               pscale(i12(j,ii)) = p2scale
            end do
            do j = 1, n13(ii)
               pscale(i13(j,ii)) = p3scale
            end do
            do j = 1, n14(ii)
               pscale(i14(j,ii)) = p4scale
            end do
            do j = 1, n15(ii)
               pscale(i15(j,ii)) = p5scale
            end do
            do j = 1, np11(ii)
               dscale_dir(ip11(j,ii)) = d1scale
               dscale_mut(ip11(j,ii)) = u1scale
            end do
            do j = 1, np12(ii)
               dscale_dir(ip12(j,ii)) = d2scale
               dscale_mut(ip12(j,ii)) = u2scale
            end do
            do j = 1, np13(ii)
               dscale_dir(ip13(j,ii)) = d3scale
               dscale_mut(ip13(j,ii)) = u3scale
            end do
            do j = 1, np14(ii)
               dscale_dir(ip14(j,ii)) = d4scale
               dscale_mut(ip14(j,ii)) = u4scale
            end do
            do l3 = l1, npole3b
               k = pnum(l3)
               k2 = 3*(l3-1)
               kk = ipole(k)
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
               do m = 1,ncell
                  xr = x(kk) - x(ii)
                  yr = y(kk) - y(ii)
                  zr = z(kk) - z(ii)
                  call imager (xr,yr,zr,m)
                  r2 = xr*xr + yr* yr + zr*zr
c                  if (r2 .le. cut2) then
                  if (r2 .le. ewaldcut3b_2) then
                     r = sqrt(r2)
                     ralpha = aewald * r
                     bn(0) = erfc(ralpha) / r
                     alsq2 = 2.0d0 * aewald**2
                     alsq2n = 0.0d0
                     if (aewald .gt. 0.0d0)
     &                  alsq2n = 1.0d0 / (sqrtpi*aewald)
                     exp2a = exp(-ralpha**2)
                     do j = 1, 3
                        bfac = dble(j+j-1)
                        alsq2n = alsq2 * alsq2n
                        bn(j) = (bfac*bn(j-1)+alsq2n*exp2a) / r2
                     end do
c
c     compute the error function scaled and unscaled terms
c
                     scale3_dir = 1.0d0
                     scale5_dir = 1.0d0
                     scale7_dir = 1.0d0
                     scale3_mut = 1.0d0
                     scale5_mut = 1.0d0
                     damp = pdi * pdamp(k)
                     if (damp .ne. 0.0d0) then
                        pgamma = min(pti,thole(k))
                        damp = -pgamma * (r/damp)**3
                        if (damp .gt. -50.0d0) then
                           expdamp = exp(damp)
                           scale3_dir = 1.0d0 - expdamp
                           scale5_dir = 1.0d0 - expdamp*(1.0d0-damp)
                           scale7_dir = 1.0d0 - expdamp
     &                                 *(1.0d0-damp+0.6d0*damp**2)

                           scale3_mut=scale3_mut*(1.0d0-expdamp)
                      scale5_mut=scale5_mut*(1.0d0-(1.0d0-damp)*expdamp)
                        end if
                     end if
                     rr3 = (1.0d0-scale3_mut) / (r*r2)
                     rr5 = 3.0d0 * (1.0d0-scale5_mut) / (r*r2*r2)
                     dsc3 = scale3_dir
                     dsc5 = scale5_dir
                     dsc7 = scale7_dir
                     psc3 = scale3_dir
                     psc5 = scale5_dir
                     psc7 = scale7_dir
                     if (use_polymer) then
                        if (r2 .le. polycut2) then
                           dsc3 = scale3_dir * dscale_dir(kk)
                           dsc5 = scale5_dir * dscale_dir(kk)
                           dsc7 = scale7_dir * dscale_dir(kk)
                           psc3 = scale3_dir * pscale(kk)
                           psc5 = scale5_dir * pscale(kk)
                           psc7 = scale7_dir * pscale(kk)
                           scale3_mut = scale3_mut * dscale_mut(kk)
                           scale5_mut = scale5_mut *dscale_mut(kk)
                        end if
                     end if

               fid_xx = -rr3 + xr*xr*rr5
               fid_xy =  xr*yr*rr5
               fid_xz =  xr*zr*rr5
               fid_yx =  yr*xr*rr5
               fid_yy = -rr3 + yr*yr*rr5
               fid_yz =  yr*zr*rr5
               fid_zx =  zr*xr*rr5
               fid_zy =  zr*yr*rr5
               fid_zz = -rr3 + zr*zr*rr5

               fimd_xx = -bn(1) + xr*xr*bn(2)
               fimd_xy =  xr*yr*bn(2)
               fimd_xz =  xr*zr*bn(2)
               fimd_yx =  yr*xr*bn(2)
               fimd_yy = -bn(1) + yr*yr*bn(2)
               fimd_yz =  yr*zr*bn(2)
               fimd_zx =  zr*xr*bn(2)
               fimd_zy =  zr*yr*bn(2)
               fimd_zz = -bn(1) + zr*zr*bn(2)

               Txx = (fimd_xx - fid_xx)
               Txy = (fimd_xy - fid_xy)
               Txz = (fimd_xz - fid_xz)
               Tyx = (fimd_yx - fid_yx)
               Tyy = (fimd_yy - fid_yy)
               Tyz = (fimd_yz - fid_yz)
               Tzx = (fimd_zx - fid_zx)
               Tzy = (fimd_zy - fid_zy)
               Tzz = (fimd_zz - fid_zz)

               M_real(i2+1,k2+1) = M_real(i2+1,k2+1)+Txx
               M_real(i2+1,k2+2) = M_real(i2+1,k2+2)+Txy
               M_real(i2+1,k2+3) = M_real(i2+1,k2+3)+Txz
               M_real(i2+2,k2+1) = M_real(i2+2,k2+1)+Tyx
               M_real(i2+2,k2+2) = M_real(i2+2,k2+2)+Tyy
               M_real(i2+2,k2+3) = M_real(i2+2,k2+3)+Tyz
               M_real(i2+3,k2+1) = M_real(i2+3,k2+1)+Tzx
               M_real(i2+3,k2+2) = M_real(i2+3,k2+2)+Tzy
               M_real(i2+3,k2+3) = M_real(i2+3,k2+3)+Tzz

               if (ii.ne.kk) then
                 M_real(k2+1,i2+1) = M_real(k2+1,i2+1)+Txx
                 M_real(k2+1,i2+2) = M_real(k2+1,i2+2)+Txy
                 M_real(k2+1,i2+3) = M_real(k2+1,i2+3)+Txz
                 M_real(k2+2,i2+1) = M_real(k2+2,i2+1)+Tyx
                 M_real(k2+2,i2+2) = M_real(k2+2,i2+2)+Tyy
                 M_real(k2+2,i2+3) = M_real(k2+2,i2+3)+Tyz
                 M_real(k2+3,i2+1) = M_real(k2+3,i2+1)+Tzx
                 M_real(k2+3,i2+2) = M_real(k2+3,i2+2)+Tzy
                 M_real(k2+3,i2+3) = M_real(k2+3,i2+3)+Tzz
               end if
                     drr3 = (1.0d0-dsc3) / (r*r2)
                     drr5 = 3.0d0 * (1.0d0-dsc5) / (r*r2*r2)
                     drr7 = 15.0d0 * (1.0d0-dsc7) / (r*r2*r2*r2)
                     prr3 = (1.0d0-psc3) / (r*r2)
                     prr5 = 3.0d0 * (1.0d0-psc5) / (r*r2*r2)
                     prr7 = 15.0d0 * (1.0d0-psc7) / (r*r2*r2*r2)
                     dir = dix*xr + diy*yr + diz*zr
                     qix = qixx*xr + qixy*yr + qixz*zr
                     qiy = qixy*xr + qiyy*yr + qiyz*zr
                     qiz = qixz*xr + qiyz*yr + qizz*zr
                     qir = qix*xr + qiy*yr + qiz*zr
                     dkr = dkx*xr + dky*yr + dkz*zr
                     qkx = qkxx*xr + qkxy*yr + qkxz*zr
                     qky = qkxy*xr + qkyy*yr + qkyz*zr
                     qkz = qkxz*xr + qkyz*yr + qkzz*zr
                     qkr = qkx*xr + qky*yr + qkz*zr
                     fim(1) = -xr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                           - bn(1)*dkx + 2.0d0*bn(2)*qkx
                     fim(2) = -yr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                           - bn(1)*dky + 2.0d0*bn(2)*qky
                     fim(3) = -zr*(bn(1)*ck-bn(2)*dkr+bn(3)*qkr)
     &                           - bn(1)*dkz + 2.0d0*bn(2)*qkz
                     fkm(1) = xr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                           - bn(1)*dix - 2.0d0*bn(2)*qix
                     fkm(2) = yr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                           - bn(1)*diy - 2.0d0*bn(2)*qiy
                     fkm(3) = zr*(bn(1)*ci+bn(2)*dir+bn(3)*qir)
     &                           - bn(1)*diz - 2.0d0*bn(2)*qiz
                     fid(1) = -xr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                           - drr3*dkx + 2.0d0*drr5*qkx
                     fid(2) = -yr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                           - drr3*dky + 2.0d0*drr5*qky
                     fid(3) = -zr*(drr3*ck-drr5*dkr+drr7*qkr)
     &                           - drr3*dkz + 2.0d0*drr5*qkz
                     fkd(1) = xr*(drr3*ci+drr5*dir+drr7*qir)
     &                           - drr3*dix - 2.0d0*drr5*qix
                     fkd(2) = yr*(drr3*ci+drr5*dir+drr7*qir)
     &                           - drr3*diy - 2.0d0*drr5*qiy
                     fkd(3) = zr*(drr3*ci+drr5*dir+drr7*qir)
     &                           - drr3*diz - 2.0d0*drr5*qiz
                     fip(1) = -xr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                           - prr3*dkx + 2.0d0*prr5*qkx
                     fip(2) = -yr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                           - prr3*dky + 2.0d0*prr5*qky
                     fip(3) = -zr*(prr3*ck-prr5*dkr+prr7*qkr)
     &                           - prr3*dkz + 2.0d0*prr5*qkz
                     fkp(1) = xr*(prr3*ci+prr5*dir+prr7*qir)
     &                           - prr3*dix - 2.0d0*prr5*qix
                     fkp(2) = yr*(prr3*ci+prr5*dir+prr7*qir)
     &                           - prr3*diy - 2.0d0*prr5*qiy
                     fkp(3) = zr*(prr3*ci+prr5*dir+prr7*qir)
     &                           - prr3*diz - 2.0d0*prr5*qiz
c
c     increment the field at each site due to this interaction
c
                     do j = 1, 3
                        field(j,l1) = field(j,l1) + fim(j) - fid(j)
                        fieldp(j,l1) = fieldp(j,l1) + fim(j) - fip(j)

                        if (ii .ne. kk) then
                           field(j,l3) = field(j,l3) + fkm(j) - fkd(j)
                           fieldp(j,l3) = fieldp(j,l3) + fkm(j) - fkp(j)
                        end if
                     end do
                  end if
               end do
            end do
c
c     reset interaction scaling coefficients for connected atoms
c
            do j = 1, n12(ii)
               pscale(i12(j,ii)) = 1.0d0
            end do
            do j = 1, n13(ii)
               pscale(i13(j,ii)) = 1.0d0
            end do
            do j = 1, n14(ii)
               pscale(i14(j,ii)) = 1.0d0
            end do
            do j = 1, n15(ii)
               pscale(i15(j,ii)) = 1.0d0
            end do
            do j = 1, np11(ii)
               dscale_dir(ip11(j,ii)) = 1.0d0
               dscale_mut(ip11(j,ii)) = 1.0d0
            end do
            do j = 1, np12(ii)
               dscale_dir(ip12(j,ii)) = 1.0d0
               dscale_mut(ip12(j,ii)) = 1.0d0
            end do
            do j = 1, np13(ii)
               dscale_dir(ip13(j,ii)) = 1.0d0
               dscale_mut(ip13(j,ii)) = 1.0d0
            end do
            do j = 1, np14(ii)
               dscale_dir(ip14(j,ii)) = 1.0d0
               dscale_mut(ip14(j,ii)) = 1.0d0
            end do
        end do
      end if
c
c     perform deallocation of some local arrays
c
      deallocate (dscale_dir)
      deallocate (dscale_mut)
      deallocate (pscale)
      return
      end

      subroutine induce0c_3b_new(aewald3b,ewaldcut3b,
     & jmax3b,kmax3b,lmax3b,npole3b,pnum,uind,uinp)
      implicit none
      include 'sizes.i'
c      include 'atoms.i'
      include 'boxes.i'
c      include 'combo.i'
      include 'ewald.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar2.i'
      include 'polpot.i'
      include 'potent.i'
      include 'units.i'
c      include 'uprior.i'
c      include 'molcul.i'
c      include 'pme.i'
      integer i,j,k,ii,i1,k1
      integer l1,l3,i2,j2,k2
c      integer iter,maxiter
c      real*8 eps,term
c      real*8 epsd,epsp
c      real*8 epsold
c      real*8 udsum,upsum
      real*8 ucell(3)
      real*8 ucellp(3)
      real*8 field(3,npole3b)
      real*8 fieldp(3,npole3b)
      integer npole3b,pnum(*)
      real*8 M_recip(3*npole3b,3*npole3b)
      real*8 M_real(3*npole3b,3*npole3b)
      real*8 M_tot(3*npole3b,3*npole3b)
      real timediff,timearray(2)
      logical done, use_liam
      real*8 uind(3,*)
      real*8 uinp(3,*)
      real*8 term,aewald,ewaldcut3b,aewald3b
      integer jmax3b,kmax3b,lmax3b

      aewald=aewald3b
c
c     zero out the induced dipole and the field at each site
c
      do i = 1, npole3b
         do j = 1, 3
            uind(j,i) = 0.0d0
            uinp(j,i) = 0.0d0
         end do
      end do
      do l1 = 1, 3*npole3b
         do l3 = 1, 3*npole3b
            M_recip(l1,l3) = 0.0d0
            M_real(l1,l3) = 0.0d0
            M_tot(l1,l3) = 0.0d0
         end do
      end do
      if (.not. use_polar)  return

        do i = 1, npole3b
          do j = 1, 3
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
          end do
        end do

       call udirect1_umutual_ewrp_3b(aewald3b,ewaldcut3b,jmax3b,kmax3b,
     & lmax3b,field,M_recip,npole3b,pnum)

      do l1 = 1, npole3b
        i=pnum(l1)
        do j = 1, 3
c          fieldp(j,i) = field(j,i)
          fieldp(j,l1)=field(j,l1)
        end do
      end do

       call udirect2a_umutual_ewrl_3b(aewald3b,ewaldcut3b,field,fieldp,
     &  M_real,npole3b,pnum)

      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do l1 = 1, npole3b
        i=pnum(l1)
         do j = 1, 3
            field(j,l1)=field(j,l1) + term*rpole(j+1,i)
            fieldp(j,l1)=fieldp(j,l1) + term*rpole(j+1,i)
         end do
      end do

      do l1 = 1, 3*npole3b-1
        do l3 = l1+1, 3*npole3b
           M_tot(l1,l3) = - (M_real(l1,l3) + M_recip(l1,l3))
           M_tot(l3,l1) = - (M_real(l3,l1) + M_recip(l3,l1))
        end do
      end do

      term = (4.0d0/3.0d0) * aewald**3 / sqrtpi
      do l1 = 1, npole3b
        i = pnum(l1)
        do j = 1, 3
          i1 = 3*(l1-1)+j
          M_tot(i1,i1) = - (M_real(i1,i1) + M_recip(i1,i1)) -
     &                      term + 1.0d0/polarity(i)
        end do
      end do

      call invert(3*npole3b,M_tot)

      do l1 = 1, npole3b
        i = pnum(l1)
        do i1 = 1, 3
          i2 = 3*(l1-1) + i1
          do l3 = 1, npole3b
             k = pnum(l3)
             k2 = 3*(l3-1)
c             uind(i1,i) = uind(i1,i)+M_tot(i2,k2+1)*field(1,k)+
c     &                               M_tot(i2,k2+2)*field(2,k)+
c     &                               M_tot(i2,k2+3)*field(3,k)
c             uind(i1,i) = uind(i1,i)+M_tot(i2,k2+1)*field(1,l3)+
c     &                               M_tot(i2,k2+2)*field(2,l3)+
c     &                               M_tot(i2,k2+3)*field(3,l3)

             uind(i1,l1) = uind(i1,l1)+M_tot(i2,k2+1)*field(1,l3)+
     &                               M_tot(i2,k2+2)*field(2,l3)+
     &                               M_tot(i2,k2+3)*field(3,l3)

          end do
        end do
      end do

      do l1 = 1, npole3b
c         i = pnum(l1)
         do j = 1, 3
            uinp(j,l1)=uind(j,l1)
         end do
      end do

      return
      end
