c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine induce0c  --  Ewald induced dipoles via loop  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "induce0c" computes the induced dipole moments at polarizable
c     sites using particle mesh Ewald summation a double loop
c
c
      subroutine induce0c_NoPol_Returnfield(field,fieldp)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'combo.i'
      include 'ewald.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'potent.i'
      include 'units.i'
      include 'uprior.i'
      integer i,j,k,ii
      integer iter,maxiter
      real*8 eps,term
      real*8 epsd,epsp
      real*8 epsold
      real*8 udsum,upsum
      real*8 ucell(3)
      real*8 ucellp(3)
      real*8 field(3,*)
      real*8 fieldp(3,*)
c      real*8, allocatable :: udir(:,:)
c      real*8, allocatable :: udirp(:,:)
c      real*8, allocatable :: uold(:,:)
c      real*8, allocatable :: uoldp(:,:)
c      real*8, allocatable :: field(:,:)
c      real*8, allocatable :: fieldp(:,:)
      logical done

c
c     get the reciprical space part of the electrostatic field
c
      call udirect1_reg_ewald (field)
      do i = 1, npole
         do j = 1, 3
            fieldp(j,i) = field(j,i)
         end do
      end do

c
c     get the real space portion of the electrostatic field
c
      call udirect2a (field,fieldp)
c
c     get the self-energy portion of the electrostatic field
c
      term = (4.0d0/3.0d0) * aewaldPerm**3 / sqrtpi
      do i = 1, npole3b
         do j = 1, 3
            field(j,i) = field(j,i) + term*rpole(j+1,i)
            fieldp(j,i) = fieldp(j,i) + term*rpole(j+1,i)
         end do
      end do
c
c     compute the cell dipole boundary correction to field
c
      if (boundary .eq. 'VACUUM') then
         do i = 1, 3
            ucell(i) = 0.0d0
         end do
         do i = 1, npole
            ii = ipole(i)
            ucell(1) = ucell(1) + rpole(2,i) + rpole(1,i)*x(ii)
            ucell(2) = ucell(2) + rpole(3,i) + rpole(1,i)*y(ii)
            ucell(3) = ucell(3) + rpole(4,i) + rpole(1,i)*z(ii)
         end do
         term = (4.0d0/3.0d0) * pi/volbox
         do i = 1, npole
            do j = 1, 3
               field(j,i) = field(j,i) - term*ucell(j)
               fieldp(j,i) = fieldp(j,i) - term*ucell(j)
            end do
         end do
      end if
      return
      end
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine induce0c  --  Ewald induced dipoles via loop  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "induce0c" computes the induced dipole moments at polarizable
c     sites using particle mesh Ewald summation a double loop
c
c
      subroutine induce0c_NoPol_Returnfield2(field,fieldp)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'combo.i'
      include 'ewald.i'
      include 'inform.i'
      include 'iounit.i'
      include 'math.i'
      include 'mpole.i'
      include 'polar.i'
      include 'polpot.i'
      include 'potent.i'
      include 'units.i'
      include 'uprior.i'
      integer i,j,k,ii
      integer iter,maxiter
      real*8 eps,term
      real*8 epsd,epsp
      real*8 epsold
      real*8 udsum,upsum
      real*8 ucell(3)
      real*8 ucellp(3)
      real*8 field(3,*)
      real*8 fieldp(3,*)
c      real*8, allocatable :: udir(:,:)
c      real*8, allocatable :: udirp(:,:)
c      real*8, allocatable :: uold(:,:)
c      real*8, allocatable :: uoldp(:,:)
c      real*8, allocatable :: field(:,:)
c      real*8, allocatable :: fieldp(:,:)
      logical done

c
c     get the reciprical space part of the electrostatic field
c
      call udirect1_reg_ewald (field)
      do i = 1, npole
         do j = 1, 3
            fieldp(j,i) = field(j,i)
         end do
      end do

c
c     get the real space portion of the electrostatic field
c
      call udirect2a2 (field,fieldp)
c
c     get the self-energy portion of the electrostatic field
c
      term = (4.0d0/3.0d0) * aewaldPerm**3 / sqrtpi
      do i = 1, npole3b
         do j = 1, 3
            field(j,i) = field(j,i) + term*rpole(j+1,i)
            fieldp(j,i) = fieldp(j,i) + term*rpole(j+1,i)
         end do
      end do
c
c     compute the cell dipole boundary correction to field
c
      if (boundary .eq. 'VACUUM') then
         do i = 1, 3
            ucell(i) = 0.0d0
         end do
         do i = 1, npole3b
            ii = ipole(i)
            ucell(1) = ucell(1) + rpole(2,i) + rpole(1,i)*x(ii)
            ucell(2) = ucell(2) + rpole(3,i) + rpole(1,i)*y(ii)
            ucell(3) = ucell(3) + rpole(4,i) + rpole(1,i)*z(ii)
         end do
         term = (4.0d0/3.0d0) * pi/volbox
         do i = 1, npole3b
            do j = 1, 3
               field(j,i) = field(j,i) - term*ucell(j)
               fieldp(j,i) = fieldp(j,i) - term*ucell(j)
            end do
         end do
      end if
      return
      end
c
c     #################################################################
c     ##                                                             ##
c     ##subroutine udirect1_reg_ewald - Ewald recip direct induced fd #
c     #
c     ##                                                             ##
c     #################################################################
c
c
c     "udirect1_reg_ewald" computes the reciprocal space contribution of
c
c     permanent atomic multipole moments to the electrostatic field
c     for use in finding the direct induced dipole moments via a
c     regular Ewald summation
c
c
      subroutine udirect1_reg_ewald (field)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'ewald.i'
      include 'ewreg.i'
      include 'math.i'
      include 'mpole.i'
      include 'units.i'
      integer i,j,k,l,ii
      integer jmin,jmax
      integer kmin,kmax
      integer lmin,lmax
      real*8 expterm,cut
      real*8 term,fterm
      real*8 xfr,yfr,zfr
      real*8 rj,rk,rl
      real*8 h1,h2,h3,hsq
      real*8 qf,t1,t2
      real*8 ck,dk,qk
      real*8 q1,q2,q3
      real*8 ckr(maxatm)
      real*8 skr(maxatm)
      real*8 cjk(maxatm)
      real*8 sjk(maxatm)
      real*8 cm(maxatm)
      real*8 dm(3,maxatm)
      real*8 qm(9,maxatm)
      real*8 field(3,maxatm)
      real*8 frecip
c
c     return if the Ewald coefficient is zero
c
      if (aewaldPerm .lt. 1.0d-6)  return
      term = -0.25d0 / aewaldPerm**2
      fterm = 8.0d0 * pi / volbox
c
c     set the number of vectors based on box dimensions
c
      frecip = 0.5d0
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
      do i = 1, npole
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
      end do
c
c     calculate and store the exponential factors
c
      do i = 1, npole
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
            do i = 1, npole
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
                  do i = 1, npole
                     ckr(i) = cjk(i)*elc(i,l) - sjk(i)*els(i,l)
                     skr(i) = sjk(i)*elc(i,l) + cjk(i)*els(i,l)
                     ck = cm(i)
                     dk = h1*dm(1,i) + h2*dm(2,i) + h3*dm(3,i)
                     q1 = h1*qm(1,i) + h2*qm(4,i) + h3*qm(7,i)
                     q2 = h1*qm(2,i) + h2*qm(5,i) + h3*qm(8,i)
                     q3 = h1*qm(3,i) + h2*qm(6,i) + h3*qm(9,i)
                     qk = h1*q1 + h2*q2 + h3*q3
                     t1 = t1 + (ck-qk)*skr(i) + dk*ckr(i)
                     t2 = t2 + (ck-qk)*ckr(i) - dk*skr(i)
                  end do
                  expterm = fterm * exp(term*hsq) / hsq
                  if (octahedron) then
                     if (mod(j+k+l,2) .ne. 0)  expterm = 0.0d0
                  end if
                  do i = 1, npole
                     qf = expterm * (skr(i)*t2-ckr(i)*t1)
                     field(1,i) = field(1,i) + h1*qf
                     field(2,i) = field(2,i) + h2*qf
                     field(3,i) = field(3,i) + h3*qf
                  end do
               end if
            end do
            lmin = -lmax
         end do
         kmin = -kmax
      end do
      return
      end

