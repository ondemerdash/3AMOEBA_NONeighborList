c
c
c     ################################################################
c     ##  COPYRIGHT (C) 2011 by Teresa Head-Gordon & Jay W. Ponder  ##
c     ##                     All Rights Reserved                    ##
c     ################################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine nose  --  Nose-Hoover NPT molecular dynamics  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "nose" performs a single molecular dynamics time step via
c     a Nose-Hoover extended system isothermal-isobaric algorithm
c
c     literature reference:
c
c     G. J. Martyna, M. E. Tuckerman, D. J. Tobias and M. L. Klein,
c     "Explicit Reversible Integrators for Extended Systems Dynamics",
c     Molecular Physics, 87, 1117-1157 (1996)
c
c     original version written by Teresa Head-Gordon, November 2011
c
c
c      subroutine nose (istep,dt)
c      implicit none
c      include 'sizes.i'
c      include 'atmtyp.i'
c      include 'atoms.i'
c      include 'bath.i'
c      include 'boxes.i'
c      include 'freeze.i'
c      include 'mdstuf.i'
c      include 'moldyn.i'
c      include 'units.i'
c      include 'usage.i'
c      include 'virial.i'
c      integer i,j,istep
c      real*8 dt,dt_2
c      real*8 epot,etot
c      real*8 eksum,temp
c      real*8 pres,press
c      real*8 poly,factor
c      real*8 term,expterm
c      real*8 term2,eterm2
c      real*8 e2,e4,e6,e8
c      real*8 ekin(3,3)
c      real*8 stress(3,3)
c      real*8, allocatable :: derivs(:,:)
c      save press

c
c
c     set some time values for the dynamics integration
c

      dt_2 = 0.5d0 * dt
      if (istep .eq. 1)  press = atmsph
c
c     update thermostat and barostat values, scale atomic velocities
c
      call hoover (dt,press)
c
c     get half-step velocities via Verlet recursion
c
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               v(j,i) = v(j,i) + a(j,i)*dt_2
            end do
         end if
      end do
c
c     update atomic positions via coupling to barostat
c
      term = vbar * dt_2
      term2 = term * term
      expterm = exp(term)
      eterm2 = expterm * expterm
      e2 = 1.0d0 / 6.0d0
      e4 = e2 / 20.0d0
      e6 = e4 / 42.0d0
      e8 = e6 / 72.0d0
      poly = 1.0d0 + term2*(e2+term2*(e4+term2*(e6+term2*e8)))
      poly = expterm * poly * dt
      do i = 1, n
         if (use(i)) then
            x(i) = x(i)*eterm2 + v(1,i)*poly
            y(i) = y(i)*eterm2 + v(2,i)*poly
            z(i) = z(i)*eterm2 + v(3,i)*poly
         end if
      end do
c
c     constraints under NH-NPT require the ROLL algorithm
c
      if (use_rattle)  call fatal
c
c     update the periodic box size and total volume
c
      xbox = xbox * eterm2
      ybox = ybox * eterm2
      zbox = zbox * eterm2
      call lattice
c
c     perform dynamic allocation of some local arrays
c
      allocate (derivs(3,n))
c
c     get the potential energy and atomic forces
c

c      call gradient (epot,derivs)
      if (taskid.eq.master) then

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
      ep3b2=0.0d0
      ep3b3=0.0d0

c
c     zero out each of the first derivative components
c
      do i = 1, n
         do j = 1, 3
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
            dep3b2(j,i) = 0.0d0
            dep3b3(j,i) = 0.0d0
         end do
      end do
c
c     zero out the virial and the intermolecular energy
c
      do i = 1, 3
         do j = 1, 3
            virep3b2(i,j)=0.0d0
            virep3b3(i,j)=0.0d0
            vir(j,i) = 0.0d0
         end do
      end do
      einter = 0.0d0

      end if

c      if (use_bounds .and. .not.use_rigid)  call bounds

      if (taskid.eq.master) then
         if (use_list)  call nblist
         call mollist2body
         call mollist3body
      end if

         call mpi_bcast(nmollst,nmol,mpi_integer,master,
     &   mpi_comm_world,ierr)
          call mpi_bcast(mollst,100*nmol,mpi_integer,master,
     &   mpi_comm_world,ierr)
         call mpi_bcast(nmollst3mod,nmol,mpi_integer,master,
     &   mpi_comm_world,ierr)
         call mpi_bcast(mollst3mod,800*nmol,mpi_integer,master,
     &   mpi_comm_world,ierr)
c         print*,"After mollst bcast"
         call mpi_barrier(mpi_comm_world,ierr)
        count3=0
         do k=1,nmol
          if(nmollst3mod(k).ne.0) then
            count3=count3+1
            mol3new(count3)=k
          end if
         end do

      if (use_bounds .and. .not.use_rigid)  call bounds
      cutoff = 0.0d0
      call replica (cutoff)

      if (taskid.eq.master) then

         if (use_born)  call born
c
c     alter bond and torsion constants for pisystem
c
         if (use_orbit)  call picalc
c
c     call the local geometry energy and gradient routines
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
c     call the van der Waals energy and gradient routines
c
         if (use_vdw) then
           if (vdwtyp .eq. 'LENNARD-JONES')  call elj1
           if (vdwtyp .eq. 'BUCKINGHAM')  call ebuck1
           if (vdwtyp .eq. 'MM3-HBOND')  call emm3hb1
           if (vdwtyp .eq. 'BUFFERED-14-7')  call ehal1
           if (vdwtyp .eq. 'GAUSSIAN')  call egauss1
         end if
c
c     call the electrostatic energy and gradient routines
c
         if (use_charge)  call echarge1
         if (use_chgdpl)  call echgdpl1
         if (use_dipole)  call edipole1
c      if (use_mpole .or. use_polar)  call empole1
         if (use_rxnfld)  call erxnfld1
c
c     call any miscellaneous energy and gradient routines
c
         if (use_solv)  call esolv1
         if (use_metal)  call emetal1
         if (use_geom)  call egeom1
         if (use_extra)  call extra1
          call empole1c_3b_Perm
      end if

         call mpi_bcast(rpole,13*maxatm,mpi_real8,master,
     &   mpi_comm_world,ierr)
         call mpi_barrier(mpi_comm_world,ierr)

                ep3bt_tot2=0.0d0
                ep3bt_tot3=0.0d0
                do i = 1, npole
                   do j = 1, 3
                      dep3bt_tot2(j,i) = 0.0d0
                      dep3bt_tot3(j,i) = 0.0d0
                   end do
                end do

                do i=1,3
                   do j=1,3
                      virep3bt_tot2(i,j)=0.0d0
                      virep3bt_tot3(i,j)=0.0d0
                   end do
                end do

c   
          mode = 'MPOLE'
         call switch (mode)


          offset=int(nmol/numtasks)
          remainder=mod(nmol,numtasks)
             if(taskid.le.remainder-1) then
               start=taskid*offset+1
               do moli1 =start,start+offset-1
                  call Innerloop2(moli1,ep3bt,virep3bt,dep3bt)
                 do i=1,3
                   do j=1,3
                 virep3bt_tot2(i,j)=virep3bt_tot2(i,j)+virep3bt(i,j)
                   end do
                 end do

                 do i=1,npole
                    do j=1,3
                 dep3bt_tot2(j,i)=dep3bt_tot2(j,i)+dep3bt(j,i)
                    end do
                 end do

                 ep3bt_tot2=ep3bt_tot2+ep3bt
               end do

                moli1=numtasks*offset+taskid+1
                  call Innerloop2(moli1,ep3bt,virep3bt,dep3bt)
                 do i=1,3
                   do j=1,3
                   virep3bt_tot2(i,j)=virep3bt_tot2(i,j)+virep3bt(i,j)
                   end do
                 end do
                 do i=1,npole
                    do j=1,3
                    dep3bt_tot2(j,i)=dep3bt_tot2(j,i)+dep3bt(j,i)
                    end do
                 end do
                ep3bt_tot2=ep3bt_tot2+ep3bt

                  call mpi_reduce(ep3bt_tot2,ep3b2,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt_tot2,dep3b2,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt_tot2,virep3b2,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)

             else
               start=taskid*offset+1
               do moli1 =start,start+offset-1
                  call Innerloop2(moli1,ep3bt,virep3bt,dep3bt)

                 do i=1,3
                   do j=1,3
                    virep3bt_tot2(i,j)=virep3bt_tot2(i,j)+virep3bt(i,j)
                   end do
                 end do

                 do i=1,npole
                    do j=1,3
                    dep3bt_tot2(j,i)=dep3bt_tot2(j,i)+dep3bt(j,i)
                    end do
                 end do

                ep3bt_tot2=ep3bt_tot2+ep3bt
               end do
                  call mpi_reduce(ep3bt_tot2,ep3b2,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt_tot2,dep3b2,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt_tot2,virep3b2,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
             end if

             offset=int(count3/numtasks)
             remainder=mod(count3,numtasks)

             if(taskid.le.remainder-1) then
               start=taskid*offset+1
c               do moli1 =start,start+offset-1
               do k1 =start,start+offset-1
                 moli1=mol3new(k1)
                  call Innerloop3(moli1,ep3bt,virep3bt,dep3bt)
                 do i=1,3
                   do j=1,3
                 virep3bt_tot3(i,j)=virep3bt_tot3(i,j)+virep3bt(i,j)
                   end do
                 end do

                 do i=1,npole
                    do j=1,3
                 dep3bt_tot3(j,i)=dep3bt_tot3(j,i)+dep3bt(j,i)
                    end do
                 end do

                 ep3bt_tot3=ep3bt_tot3+ep3bt
               end do

c                moli1=numtasks*offset+taskid+1
                k1=numtasks*offset+taskid+1
                moli1=mol3new(k1)
                  call Innerloop3(moli1,ep3bt,virep3bt,dep3bt)
                 do i=1,3
                   do j=1,3
                   virep3bt_tot3(i,j)=virep3bt_tot3(i,j)+virep3bt(i,j)
                   end do
                 end do
                 do i=1,npole
                    do j=1,3
                    dep3bt_tot3(j,i)=dep3bt_tot3(j,i)+dep3bt(j,i)
                    end do
                 end do
                ep3bt_tot3=ep3bt_tot3+ep3bt

                  call mpi_reduce(ep3bt_tot3,ep3b3,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt_tot3,dep3b3,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt_tot3,virep3b3,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)

             else
               start=taskid*offset+1
c               do moli1 =start,start+offset-1
               do k1=start,start+offset-1
                  moli1=mol3new(k1)

                  call Innerloop3(moli1,ep3bt,virep3bt,dep3bt)

                 do i=1,3
                   do j=1,3
                    virep3bt_tot3(i,j)=virep3bt_tot3(i,j)+virep3bt(i,j)
                   end do
                 end do

                 do i=1,npole
                    do j=1,3
                    dep3bt_tot3(j,i)=dep3bt_tot3(j,i)+dep3bt(j,i)
                    end do
                 end do

                ep3bt_tot3=ep3bt_tot3+ep3bt
               end do
                  call mpi_reduce(ep3bt_tot3,ep3b3,1,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(dep3bt_tot3,dep3b3,npole*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
                  call mpi_reduce(virep3bt_tot3,virep3b3,3*3,mpi_real8,
     &          mpi_sum,master,mpi_comm_world,ierr)
             end if

           call mpi_barrier(mpi_comm_world,ierr)

           if(taskid.eq.master) then
             ep = ep+ep3b2+ep3b3

             do i = 1, npole
               do j = 1, 3
                 dep(j,i) =dep(j,i)+dep3b2(j,i)+dep3b3(j,i)
               end do
             end do

             do i=1,3
               do j=1,3
                 vir(i,j)=vir(i,j)+virep3b2(i,j)+virep3b3(i,j)
               end do
             end do

c
c     sum up to get the total energy and first derivatives
c
              esum = eb + ea + eba + eub + eaa + eopb + eopd + eid + eit
     &          + et + ept + ebt + ett + ev + ec + ecd + ed + em
     &          + ep + er + es + elf + eg + ex
c              energy = esum
              epot=esum
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
                   derivs(j,i) = desum(j,i)
                 end do
              end do
c
c     check for an illegal value for the total energy
c
              if (isnan(esum)) then
                  write (iout,10)
   10            format (/,' GRADIENT  --  Illegal Value for the Total',
     &              ' Potential Energy')
                 call fatal
              end if

           end if

           call mpi_bcast(derivs,3*maxatm,mpi_real8,master,
     &   mpi_comm_world,ierr)
           call mpi_bcast(vir,3*3,mpi_real8,master,
     &   mpi_comm_world,ierr)
           call mpi_bcast(epot,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
           call mpi_barrier(mpi_comm_world,ierr)



c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using the Verlet recursion
c
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               a(j,i) = -convert * derivs(j,i) / mass(i)
               v(j,i) = v(j,i) + a(j,i)*dt_2
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (derivs)
c
c     constraints under NH-NPT require the ROLL algorithm
c
      if (use_rattle)  call fatal
c
c     update thermostat and barostat values, scale atomic velocities
c
      call hoover (dt,press)
c
c     set isotropic pressure to the average of tensor diagonal
c
      factor = prescon / volbox
      do i = 1, 3
         do j = 1, 3
            stress(j,i) = factor * (-vir(j,i))
         end do
      end do
      press = (stress(1,1)+stress(2,2)+stress(3,3)) / 3.0d0
c
c     accumulate the kinetic energy and its outer product
c
      call kinetic (eksum,ekin)
c
c     calculate the stress tensor for anisotropic systems
c
      do i = 1, 3
         do j = 1, 3
            stress(j,i) = factor * (2.0d0*ekin(j,i)-vir(j,i))
         end do
      end do
      pres = (stress(1,1)+stress(2,2)+stress(3,3)) / 3.0d0
c
c     get the instantaneous temperature from the kinetic energy
c
      temp = 2.0d0 * eksum / (dble(nfree) * gasconst)
      etot = epot + eksum
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      if(taskid.eq.master) then
        call mdsave (istep,dt,epot)
      end if

      return
      end
c
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine hoover  --  Nose-Hoover thermostat/barostat  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "hoover" applies a combined thermostat and barostat via a
c     Nose-Hoover chain algorithm
c
c
      subroutine hoover (dt,press)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bath.i'
      include 'boxes.i'
      include 'mdstuf.i'
      include 'moldyn.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,k
      integer nc,ns
      real*8 dt,dtc,dts
      real*8 dt2,dt4,dt8
      real*8 ekt,eksum
      real*8 df,odnf,gn1kt
      real*8 press,dpress
      real*8 expterm,scale
      real*8 w(3),ekin(3,3)
c
c
c     find kinetic energy and set an initial scale factor
c
      call kinetic (eksum,ekin)
      ekt = gasconst * kelvin
      nc = 5
      ns = 3
      dtc = dt / dble(nc)
      w(1) = 1.0d0 / (2.0d0-2.0d0**(1.0d0/3.0d0))
      w(2) = 1.0d0 - 2.0d0*w(1)
      w(3) = w(1)
      df = dble(nfree)
      odnf = 1.0d0 + 3.0d0/df
      gn1kt = (1.0d0+df) * ekt
      dpress = (press-atmsph) / prescon
      scale = 1.0d0
c
c     use multiple time steps to apply thermostat and barostat
c
      do k = 1, nc
         do j = 1, ns
            dts = w(j) * dtc
            dt2 = 0.5d0 * dts
            dt4 = 0.25d0 * dts
            dt8 = 0.125d0 * dts
c
c     update thermostat and barostat velocities and forces
c
            gnh(4) = (qnh(3)*vnh(3)*vnh(3)-ekt) / qnh(4)
            vnh(4) = vnh(4) + gnh(4)*dt4
            gnh(3) = (qnh(2)*vnh(2)*vnh(2)-ekt) / qnh(3)
            expterm = exp(-vnh(4)*dt8)
            vnh(3) = expterm * (vnh(3)*expterm+gnh(3)*dt4)
            gnh(2) = (qnh(1)*vnh(1)*vnh(1)-ekt) / qnh(2)
            expterm = exp(-vnh(3)*dt8)
            vnh(2) = expterm * (vnh(2)*expterm+gnh(2)*dt4)
            gnh(1) = (2.0d0*eksum+qbar*vbar*vbar-gn1kt) / qnh(1)
            expterm = exp(-vnh(2)*dt8)
            vnh(1) = expterm * (vnh(1)*expterm+gnh(1)*dt4)
            gbar = (2.0d0*eksum*odnf+3.0d0*volbox*dpress) / qbar
            expterm = exp(-vnh(1)*dt8)
            vbar = expterm * (vbar*expterm+gbar*dt4)
c
c     find velocity scale factor and update kinetic energy
c
            expterm = exp(-(vnh(1)+vbar*odnf)*dt2)
            scale = scale * expterm
            eksum = eksum * expterm * expterm 
c
c     update barostat and thermostat velocities and forces
c
            gbar = (2.0d0*eksum*odnf+3.0d0*volbox*dpress) / qbar
            expterm = exp(-vnh(1)*dt8)
            vbar = expterm * (vbar*expterm+gbar*dt4)
            gnh(1) = (2.0d0*eksum+qbar*vbar*vbar-gn1kt) / qnh(1)
            expterm = exp(-vnh(2)*dt8)
            vnh(1) = expterm * (vnh(1)*expterm+gnh(1)*dt4)
            gnh(2) = (qnh(1)*vnh(1)*vnh(1)-ekt) / qnh(2)
            expterm = exp(-vnh(3)*dt8)
            vnh(2) = expterm * (vnh(2)*expterm+gnh(2)*dt4)
            gnh(3) = (qnh(2)*vnh(2)*vnh(2)-ekt) / qnh(3)
            expterm = exp(-vnh(4)*dt8)
            vnh(3) = expterm * (vnh(3)*expterm+gnh(3)*dt4)
            gnh(4) = (qnh(3)*vnh(3)*vnh(3)-ekt) / qnh(4)
            vnh(4) = vnh(4) + gnh(4)*dt4
         end do
      end do
c
c     use scale factor to update the atomic velocities
c
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               v(j,i) = scale * v(j,i)
            end do
         end if
      end do
      return
      end
