c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine beeman  --  Beeman molecular dynamics step  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "beeman" performs a single molecular dynamics time step
c     via the Beeman multistep recursion formula; uses original
c     coefficients or Bernie Brooks' "Better Beeman" values
c
c     literature references:
c
c     D. Beeman, "Some Multistep Methods for Use in Molecular
c     Dynamics Calculations", Journal of Computational Physics,
c     20, 130-139 (1976)
c
c     B. R. Brooks, "Algorithms for Molecular Dynamics at Constant
c     Temperature and Pressure", DCRT Report, NIH, April 1988
c
c
      subroutine beeman (istep,dt)
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'freeze.i'
      include 'mdstuf.i'
      include 'moldyn.i'
      include 'units.i'
      include 'usage.i'
      integer i,j,istep
      real*8 dt,dt_x,factor
      real*8 etot,eksum,epot
      real*8 temp,pres
      real*8 part1,part2
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      real*8, allocatable :: derivs(:,:)
c
c
c     set time values and coefficients for Beeman integration
c
      factor = dble(bmnmix)
      dt_x = dt / factor
      part1 = 0.5d0*factor + 1.0d0
      part2 = part1 - 2.0d0
c
c     make half-step temperature and pressure corrections
c
      call temper (dt)
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
      allocate (derivs(3,n))
c
c     store the current atom positions, then find half-step
c     velocities and full-step positions via Beeman recursion
c
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               v(j,i) = v(j,i) + (part1*a(j,i)-aalt(j,i))*dt_x
            end do
            xold(i) = x(i)
            yold(i) = y(i)
            zold(i) = z(i)
            x(i) = x(i) + v(1,i)*dt
            y(i) = y(i) + v(2,i)*dt
            z(i) = z(i) + v(3,i)*dt
         end if
      end do
c
c     get constraint-corrected positions and half-step velocities
c
      if (use_rattle)  call rattle (dt,xold,yold,zold)
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
c     find the full-step velocities using the Beeman recursion
c
      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               aalt(j,i) = a(j,i)
               a(j,i) = -convert * derivs(j,i) / mass(i)
               v(j,i) = v(j,i) + (part2*a(j,i)+aalt(j,i))*dt_x
            end do
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)
      deallocate (derivs)
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     make full-step temperature and pressure corrections
c
      call temper2 (dt,eksum,ekin,temp)
      call pressure (dt,epot,ekin,temp,pres,stress)
c
c     total energy is sum of kinetic and potential energies
c
      etot = eksum + epot
c
c     compute statistics and save trajectory for this step
c
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      if(taskid.eq.master) then
        call mdsave (istep,dt,epot)
      end if
      call mdrest (istep)

      return
      end
