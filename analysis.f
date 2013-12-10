c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine analysis  --  energy components and analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "analysis" calls the series of routines needed to calculate
c     the potential energy and perform energy partitioning analysis
c     in terms of type of interaction or atom number
c
c
      subroutine analysis (energy)
      implicit none
      include 'sizes.i'
      include 'analyz.i'
      include 'atoms.i'
      include 'bound.i'
      include 'cutoff.i'
      include 'combo.i'
      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'iounit.i'
      include 'potent.i'
      include 'vdwpot.i'
      include 'deriv.i'
      include 'mpole.i'
      include 'virial.i'
      integer i,j,i1,i3,ii,k
      real*8 energy,step
      real*8 cutoff
      real*8 epliam(3),forcp1(3)
      real*8 ethg(3),fthg(3)
      real*8 epArr(3)
      real*8 virtest(3,3)
c
c
c     zero out each of the potential energy components
c
      step=0.0000001d0
      do i3=1,3
c        x(3)=x(3)+step
      eb = 0.0d0
      ea = 0.0d0
      eba = 0.0d0
      eub = 0.0d0
      eaa = 0.0d0
      eopb = 0.0d0
      eopd = 0.0d0
      eid = 0.0d0
      eit = 0.0d0
      et = 0.0d0
      ept = 0.0d0
      ebt = 0.0d0
      ett = 0.0d0
      ev = 0.0d0
      ec = 0.0d0
      ecd = 0.0d0
      ed = 0.0d0
      em = 0.0d0
      ep = 0.0d0
      er = 0.0d0
      es = 0.0d0
      elf = 0.0d0
      eg = 0.0d0
      ex = 0.0d0
c
c     zero out energy partitioning components for each atom
c
      do i = 1, n
        do j=1,3
            deb(j,i) = 0.0d0
            dea(j,i) = 0.0d0
            deba(j,i) = 0.0d0
            deub(j,i) = 0.0d0
            deaa(j,i) = 0.0d0
            deopb(j,i) = 0.0d0
            deopd(j,i) = 0.0d0
            deid(j,i) = 0.0d0
            deit(j,i) = 0.0d0
            det(j,i) = 0.0d0
            dept(j,i) = 0.0d0
            debt(j,i) = 0.0d0
            dett(j,i) = 0.0d0
            dev(j,i) = 0.0d0
            dec(j,i) = 0.0d0
            decd(j,i) = 0.0d0
            ded(j,i) = 0.0d0
            dem(j,i) = 0.0d0
            dep(j,i) = 0.0d0
            der(j,i) = 0.0d0
            des(j,i) = 0.0d0
            delf(j,i) = 0.0d0
            deg(j,i) = 0.0d0
            dex(j,i) = 0.0d0
         end do
         aeb(i) = 0.0d0
         aea(i) = 0.0d0
         aeba(i) = 0.0d0
         aeub(i) = 0.0d0
         aeaa(i) = 0.0d0
         aeopb(i) = 0.0d0
         aeopd(i) = 0.0d0
         aeid(i) = 0.0d0
         aeit(i) = 0.0d0
         aet(i) = 0.0d0
         aept(i) = 0.0d0
         aebt(i) = 0.0d0
         aett(i) = 0.0d0
         aev(i) = 0.0d0
         aec(i) = 0.0d0
         aecd(i) = 0.0d0
         aed(i) = 0.0d0
         aem(i) = 0.0d0
         aep(i) = 0.0d0
         aer(i) = 0.0d0
         aes(i) = 0.0d0
         aelf(i) = 0.0d0
         aeg(i) = 0.0d0
         aex(i) = 0.0d0
      end do

      do i=1,3
         do j=1,3
           vir(i,j)=0.0d0
         end do 
      end do 
c
c     zero out the total intermolecular energy
c
      einter = 0.0d0
c
c     maintain any periodic boundary conditions
c
      if (use_bounds .and. .not.use_group)  call bounds

        x(3)=x(3)+step

c
c     remove any previous use of the replicates method
c
      cutoff = 0.0d0
      call replica (cutoff)
c
c     update the pairwise interaction neighbor lists
c
      if (use_list)  call nblist
c
c     many implicit solvation models require Born radii
c
      if (use_born)  call born
c
c     alter bond and torsion constants for pisystem
c
      if (use_orbit)  call picalc
c
c     call the local geometry energy component routines
c
      if (use_bond)  call ebond1
      if (use_angle)  call eangle1
      if (use_strbnd)  call estrbnd1
      if (use_urey)  call eurey1
      if (use_angang)  call eangang1
      if (use_opbend)  call eopbend1
      if (use_opdist)  call eopdist1
      if (use_improp)  call eimprop1
      if (use_imptor)  call eimptor1
      if (use_tors)  call etors1
      if (use_pitors)  call epitors1
      if (use_strtor)  call estrtor1
      if (use_tortor)  call etortor1
c
c     call the van der Waals energy component routines
c
      if (use_vdw) then
         if (vdwtyp .eq. 'LENNARD-JONES')  call elj1
         if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck1
         if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb1
         if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal1
         if (vdwtyp .eq. 'GAUSSIAN')  call egauss1
      end if
c
c     call the electrostatic energy component routines
c
      if (use_charge)  call echarge1
      if (use_chgdpl)  call echgdpl1
      if (use_dipole)  call edipole1
c      if (use_mpole .or. use_polar)  call empole1c_3b
      if (use_mpole .or. use_polar)  call empole1c_1a_3b
      if (use_rxnfld)  call erxnfld1
c
c     call any miscellaneous energy component routines
c
      if (use_solv)  call esolv1
      if (use_metal)  call emetal1
      if (use_geom)  call egeom1
      if (use_extra)  call extra1
c
c     sum up to give the total potential energy
c
      esum = eb + ea + eba + eub + eaa + eopb + eopd + eid + eit
     &          + et + ept + ebt + ett + ev + ec + ecd + ed + em
     &          + ep + er + es + elf + eg + ex
      energy = esum
c
c     sum up to give the total potential energy per atom
c
      do i = 1, n
         aesum(i) = aeb(i) + aea(i) + aeba(i) + aeub(i) + aeaa(i)
     &                 + aeopb(i) + aeopd(i) + aeid(i) + aeit(i)
     &                 + aet(i) + aept(i) + aebt(i) + aett(i) + aev(i)
     &                 + aec(i) + aecd(i) + aed(i) + aem(i) + aep(i)
     &                 + aer(i) + aes(i) + aelf(i) + aeg(i) + aex(i)
      end do
      do i = 1, n
         do j = 1, 3
            desum(j,i) = deb(j,i) + dea(j,i) + deba(j,i)
     &                      + deub(j,i) + deaa(j,i) + deopb(j,i)
     &                      + deopd(j,i) + deid(j,i) + deit(j,i)
     &                      + det(j,i) + dept(j,i) + debt(j,i)
     &                      + dett(j,i) + dev(j,i) + dec(j,i)
     &                      + decd(j,i) + ded(j,i) + dem(j,i)
     &                      + dep(j,i) + der(j,i) + des(j,i)
     &                      + delf(j,i) + deg(j,i) + dex(j,i)
         end do
      end do
      ethg(i3)=esum
      fthg(i3)=desum(1,3)
      epArr(i3)=ep
      end do
      
      
         do j=1,3
            do k=1,3
              virtest(j,k)=0.0d0
            end do
         end do

      do i=1,n
          virtest(1,1)=virtest(1,1)+x(i)*desum(1,i)
          virtest(1,2)=virtest(1,2)+x(i)*desum(2,i)
          virtest(1,3)=virtest(1,3)+x(i)*desum(3,i)
          virtest(2,1)=virtest(2,1)+y(i)*desum(1,i)
          virtest(2,2)=virtest(2,2)+y(i)*desum(2,i)
          virtest(2,3)=virtest(2,3)+y(i)*desum(3,i)
          virtest(3,1)=virtest(3,1)+z(i)*desum(1,i)
          virtest(3,2)=virtest(3,2)+z(i)*desum(2,i)
          virtest(3,3)=virtest(3,3)+z(i)*desum(3,i)
      end do

c      do i = 1, n
c          print*,'dem in analysis x', i, dem(1,i)
c          print*,'dem in analysis y', i, dem(2,i)
c          print*,'dem in analysis z', i, dem(3,i)
c      end do
      do i = 1, n
          print*,'dep in analysis x', i, dep(1,i)
          print*,'dep in analysis y', i, dep(2,i)
          print*,'dep in analysis z', i, dep(3,i)
      end do
      write(6,*)(ethg(3)-ethg(1))/2.0d0/step,fthg(2)
      write(6,*)'Total PE',ec,ecd,ed,em,ep
      print*,'PolEng each iter.',epArr(1),epArr(2),epArr(3)
      print*, 'Virialxx=',vir(1,1)
      print*, 'Virialxy=',vir(1,2)
      print*, 'Virialxz=',vir(1,3)
      print*, 'Virialyx=',vir(2,1)
      print*, 'Virialyy=',vir(2,2)
      print*, 'Virialyz=',vir(2,3)
      print*, 'Virialzx=',vir(3,1)
      print*, 'Virialzy=',vir(3,2)
      print*, 'Virialzz=',vir(3,3)

c      print*, 'VirialTESTxx=',virtest(1,1)
c      print*, 'VirialTESTxy=',virtest(1,2)
c      print*, 'VirialTESTxz=',virtest(1,3)
c      print*, 'VirialTESTyx=',virtest(2,1)
c      print*, 'VirialTESTyy=',virtest(2,2)
c      print*, 'VirialTESTyz=',virtest(2,3)
c      print*, 'VirialTESTzx=',virtest(3,1)
c      print*, 'VirialTESTzy=',virtest(3,2)
c      print*, 'VirialTESTzz=',virtest(3,3)

c     check for an illegal value for the total energy
c
      if (isnan(esum)) then
         write (iout,10)
   10    format (/,' ANALYSIS  --  Illegal Value for the Total',
     &              ' Potential Energy')
         call fatal
      end if
      return
      end
