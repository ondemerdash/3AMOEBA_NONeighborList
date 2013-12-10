c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  program dynamic  --  run molecular or stochastic dynamics  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "dynamic" computes a molecular or stochastic dynamics trajectory
c     in one of the standard statistical mechanical ensembles and using
c     any of several possible integration methods
c
c
      program dynamic_parallel_nolist
      implicit none
      include 'sizes.i'
      include 'atmtyp.i'

      include 'atoms.i'
      include 'freeze.i'
      include 'moldyn.i'
      include 'units.i'
      include 'molcul.i'
      include 'bath.i'
      include 'boxes.i'
      include 'mpif.h'
      include 'mpole.i'
      include 'bond.i'
      include 'bound.i'
      include 'inform.i'
      include 'iounit.i'
      include 'keys.i'
      include 'mdstuf.i'
      include 'potent.i'
      include 'solute.i'
      include 'stodyn.i'
      include 'usage.i'
      include 'virial.i'
      include 'cutoff.i'
      include 'deriv.i'
      include 'energi.i'
      include 'inter.i'
      include 'rigid.i'
      include 'vdwpot.i'
      integer i,istep,nstep
      integer mode,next,k1,k
      real*8 dt,dtdump
      logical exist,query
      character*20 keyword
      character*120 record
      character*120 string
      integer ierr,taskid,numtasks,master,moli1
      integer stat(MPI_STATUS_SIZE)
      integer j,offset,remainder,start,last
      real*8 cutoff
      real*8 dt_2
      real*8 etot,epot
      real*8 eksum
      real*8 wall,cpu
      real*8 temp,pres
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      real*8, allocatable :: derivs(:,:)
      real*8 press
      real*8 poly,factor
      real*8 term,expterm
      real*8 term2,eterm2
      real*8 e2,e4,e6,e8
      real*8 dt_x
      real*8 part1,part2
      real*8 ep3bt,virep3bt(3,3),ep3bt_tot2,virep3bt_tot2(3,3)
      real*8 ep3bt_tot3,virep3bt_tot3(3,3)
      real*8 ep3b2,ep3b3,virep3b2(3,3),virep3b3(3,3)
      real*8, allocatable :: dep3bt(:,:)
      real*8, allocatable :: dep3bt_tot2(:,:)
      real*8, allocatable :: dep3b2(:,:)
      real*8, allocatable :: dep3bt_tot3(:,:)
      real*8, allocatable :: dep3b3(:,:)
      character*6 switchmode
      save press

c
c
c     set up the structure and molecular mechanics calculation
c
      call mpi_init(ierr)
      call mpi_comm_rank(mpi_comm_world,taskid,ierr)
      call mpi_comm_size(mpi_comm_world,numtasks,ierr)

      master=0

      if(taskid.eq.master) then
      call settime
      end if


      call initial
      call getxyz
      call mechanic
      allocate (dep3bt(3,npole))
      allocate (dep3b2(3,npole))
      allocate (dep3b3(3,npole))
      allocate (dep3bt_tot2(3,npole))
      allocate (dep3bt_tot3(3,npole))

c
c     initialize the temperature, pressure and coupling baths
c
      kelvin = 0.0d0
      atmsph = 0.0d0
      isothermal = .false.
      isobaric = .false.
c
c     check for keywords containing any altered parameters
c
      integrate = 'BEEMAN'
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         string = record(next:120)
         if (keyword(1:11) .eq. 'INTEGRATOR ') then
            call getword (record,integrate,next)
            call upcase (integrate)
         end if
      end do
c
c     initialize the simulation length as number of time steps
c
      query = .true.
      call nextarg (string,exist)
      if (exist) then
         read (string,*,err=10,end=10)  nstep
         query = .false.
      end if
   10 continue
      if (query) then
         write (iout,20)
   20    format (/,' Enter the Number of Dynamics Steps to be',
     &              ' Taken :  ',$)
         read (input,30)  nstep
   30    format (i10)
      end if
c
c     get the length of the dynamics time step in picoseconds
c
      dt = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=40,end=40)  dt
   40 continue
      do while (dt .lt. 0.0d0)
         write (iout,50)
   50    format (/,' Enter the Time Step Length in Femtoseconds',
     &              ' [1.0] :  ',$)
         read (input,60,err=70)  dt
   60    format (f20.0)
         if (dt .le. 0.0d0)  dt = 1.0d0
   70    continue
      end do
      dt = 0.001d0 * dt
c
c     enforce bounds on thermostat and barostat coupling times
c
      tautemp = max(tautemp,dt)
      taupres = max(taupres,dt)
c
c     set the time between trajectory snapshot coordinate dumps
c
      dtdump = -1.0d0
      call nextarg (string,exist)
      if (exist)  read (string,*,err=80,end=80)  dtdump
   80 continue
      do while (dtdump .lt. 0.0d0)
         write (iout,90)
   90    format (/,' Enter Time between Dumps in Picoseconds',
     &              ' [0.1] :  ',$)
         read (input,100,err=110)  dtdump
  100    format (f20.0)
         if (dtdump .le. 0.0d0)  dtdump = 0.1d0
  110    continue
      end do
      iwrite = nint(dtdump/dt)
c
c     get choice of statistical ensemble for periodic system
c
      if (use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=120,end=120)  mode
  120    continue
         do while (mode.lt.1 .or. mode.gt.4)
            write (iout,130)
  130       format (/,' Available Statistical Mechanical Ensembles :',
     &              //,4x,'(1) Microcanonical (NVE)',
     &              /,4x,'(2) Canonical (NVT)',
     &              /,4x,'(3) Isoenthalpic-Isobaric (NPH)',
     &              /,4x,'(4) Isothermal-Isobaric (NPT)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,140,err=150)  mode
  140       format (i10)
            if (mode .le. 0)  mode = 1
  150       continue
         end do
         if (integrate.eq.'BUSSI' .or. integrate.eq.'NOSE-HOOVER'
     &                .or. integrate.eq.'GHMC') then
            if (mode .ne. 4) then
               mode = 4
               write (iout,160)
  160          format (/,' Switching to NPT Ensemble as Required',
     &                    ' by Chosen Integrator')
            end if
         end if
         if (mode.eq.2 .or. mode.eq.4) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=170,end=170)  kelvin
  170       continue
            do while (kelvin .lt. 0.0d0)
               write (iout,180)
  180          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,190,err=200)  kelvin
  190          format (f20.0)
               if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  200          continue
            end do
         end if
         if (mode.eq.3 .or. mode.eq.4) then
            isobaric = .true.
            atmsph = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=210,end=210)  atmsph
  210       continue
            do while (atmsph .lt. 0.0d0)
               write (iout,220)
  220          format (/,' Enter the Desired Pressure in Atm',
     &                    ' [1.0] :  ',$)
               read (input,230,err=240)  atmsph
  230          format (f20.0)
               if (atmsph .le. 0.0d0)  atmsph = 1.0d0
  240          continue
            end do
         end if
      end if
c
c     use constant energy or temperature for nonperiodic system
c
      if (.not. use_bounds) then
         mode = -1
         call nextarg (string,exist)
         if (exist)  read (string,*,err=250,end=250)  mode
  250    continue
         do while (mode.lt.1 .or. mode.gt.2)
            write (iout,260)
  260       format (/,' Available Simulation Control Modes :',
     &              //,4x,'(1) Constant Total Energy Value (E)',
     &              /,4x,'(2) Constant Temperature via Thermostat (T)',
     &              //,' Enter the Number of the Desired Choice',
     &                 ' [1] :  ',$)
            read (input,270,err=280)  mode
  270       format (i10)
            if (mode .le. 0)  mode = 1
  280       continue
         end do
         if (mode .eq. 2) then
            isothermal = .true.
            kelvin = -1.0d0
            call nextarg (string,exist)
            if (exist)  read (string,*,err=290,end=290)  kelvin
  290       continue
            do while (kelvin .lt. 0.0d0)
               write (iout,300)
  300          format (/,' Enter the Desired Temperature in Degrees',
     &                    ' K [298] :  ',$)
               read (input,310,err=320)  kelvin
  310          format (f20.0)
               if (kelvin .le. 0.0d0)  kelvin = 298.0d0
  320          continue
            end do
         end if
      end if
c
c     initialize any holonomic constraints and setup dynamics
       call mpi_barrier(mpi_comm_world,ierr)
c
      call shakeup
       call mpi_barrier(mpi_comm_world,ierr)

      call mdinit
       call mpi_barrier(mpi_comm_world,ierr)

c
c     print out a header line for the dynamics computation
c
      if (integrate .eq. 'VERLET') then
         write (iout,330)
  330    format (/,' Molecular Dynamics Trajectory via',
     &              ' Velocity Verlet Algorithm')
      else if (integrate .eq. 'STOCHASTIC') then
         write (iout,340)
  340    format (/,' Stochastic Dynamics Trajectory via',
     &              ' Velocity Verlet Algorithm')
      else if (integrate .eq. 'BUSSI') then
         write (iout,350)
  350    format (/,' Molecular Dynamics Trajectory via',
     &              ' Bussi-Parrinello NPT Algorithm')
      else if (integrate .eq. 'NOSE-HOOVER') then
         write (iout,360)
  360    format (/,' Molecular Dynamics Trajectory via',
     &              ' Nose-Hoover NPT Algorithm')
      else if (integrate .eq. 'GHMC') then
         write (iout,370)
  370    format (/,' Stochastic Dynamics Trajectory via',
     &              ' Generalized Hybrid Monte Carlo')
      else if (integrate .eq. 'RIGIDBODY') then
         write (iout,380)
  380    format (/,' Molecular Dynamics Trajectory via',
     &              ' Rigid Body Algorithm')
      else if (integrate .eq. 'RESPA') then
         write (iout,390)
  390    format (/,' Molecular Dynamics Trajectory via',
     &              ' r-RESPA MTS Algorithm')
      else
         write (iout,400)
  400    format (/,' Molecular Dynamics Trajectory via',
     &              ' Modified Beeman Algorithm')
      end if

c       call mpi_barrier(mpi_comm_world,ierr)
c
c     integrate equations of motion to take a time step
c
      do istep = 1, nstep
         if (integrate .eq. 'VERLET') then
c            call verlet (istep,dt)
c
c
c     set some time values for the dynamics integration
c
      dt_2 = 0.5d0 * dt
c
c     make half-step temperature and pressure corrections
c

c      if (taskid.eq.master) then
        call temper (dt)
c      end if

c         call mpi_bcast(v,3*n,mpi_real8,master,
c     &   mpi_comm_world,ierr)
c         call mpi_bcast(a,3*n,mpi_real8,master,
c     &   mpi_comm_world,ierr)
c         call mpi_barrier(mpi_comm_world,ierr)

c
c     perform dynamic allocation of some local arrays
c

      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
      allocate (derivs(3,n))

      do i = 1, n
         if (use(i)) then
            do j = 1, 3
               v(j,i) = v(j,i) + a(j,i)*dt_2
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
c         call mollist2body
c         call mollist3body
      end if    

c         call mpi_bcast(nmollst,nmol,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c          call mpi_bcast(mollst,100*nmol,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c         call mpi_bcast(nmollst3mod,nmol,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c         call mpi_bcast(mollst3mod,800*nmol,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c         print*,"After mollst bcast"
c         call mpi_barrier(mpi_comm_world,ierr)

c        nonzerocount3=0
c         do k=1,nmol
c          if(nmollst3mod(k).ne.0) then
c            nonzerocount3=nonzerocount3+1
c            mol3new(nonzerocount3)=k
c          end if
c         end do

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
       print*,"After rpole mpi_bcast"
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
          switchmode = 'MPOLE'
         call switch (switchmode)
         call mpi_barrier(mpi_comm_world,ierr)


          offset=int((nmol-1)/numtasks)
          remainder=mod((nmol-1),numtasks)

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

c             offset=int(nonzerocount3/numtasks)
c             remainder=mod(nonzerocount3,numtasks)
             offset=int((nmol-2)/numtasks)
             remainder=mod((nmol-2),numtasks)

             if(taskid.le.remainder-1) then
               start=taskid*offset+1
               do moli1 =start,start+offset-1
c               do k1 =start,start+offset-1
c                 moli1=mol3new(k1)
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

                moli1=numtasks*offset+taskid+1
c                k1=numtasks*offset+taskid+1
c                moli1=mol3new(k1)
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
c               start=taskid*offset+1
               do moli1 =start,start+offset-1
c               do k1=start,start+offset-1
c                  moli1=mol3new(k1)

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
                  write (iout,14)
   14            format (/,' GRADIENT  --  Illegal Value for the Total',
     &              ' Potential Energy')
                 call fatal
              end if

           end if

           call mpi_bcast(derivs,3*n,mpi_real8,master,
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
c      if(taskid.eq.master) then
        do i = 1, n
          if (use(i)) then
            do j = 1, 3
               a(j,i) = -convert * derivs(j,i) / mass(i)
               v(j,i) = v(j,i) + a(j,i)*dt_2
            end do
          end if
        end do

c
c     find the constraint-corrected full-step velocities
c     
c   WORK TO BE DONE BY MASTER ONLY, I THINK    
        
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

      deallocate (derivs)

c      end if
c
c     perform deallocation of some local arrays
c
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)


         else if (integrate .eq. 'STOCHASTIC') then
            call sdstep (istep,dt)
         else if (integrate .eq. 'BUSSI') then
            call bussi (istep,dt)
         else if (integrate .eq. 'NOSE-HOOVER') then
c            call nose (istep,dt)
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
c         call mollist2body
c         call mollist3body
      end if

c         call mpi_bcast(nmollst,nmol,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c          call mpi_bcast(mollst,100*nmol,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c         call mpi_bcast(nmollst3mod,nmol,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c         call mpi_bcast(mollst3mod,800*nmol,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c         print*,"After mollst bcast"
c         call mpi_barrier(mpi_comm_world,ierr)
c        nonzerocount3=0
c         do k=1,nmol
c          if(nmollst3mod(k).ne.0) then
c            nonzerocount3=nonzerocount3+1
c            mol3new(nonzerocount3)=k
c          end if
c         end do

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
         switchmode = 'MPOLE'
         call switch (switchmode)
         call mpi_barrier(mpi_comm_world,ierr)


          offset=int((nmol-1)/numtasks)
          remainder=mod((nmol-1),numtasks)
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

             offset=int((nmol-2)/numtasks)
             remainder=mod((nmol-2),numtasks)

             if(taskid.le.remainder-1) then
               start=taskid*offset+1
               do moli1 =start,start+offset-1
c               do k1 =start,start+offset-1
c                 moli1=mol3new(k1)
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

                moli1=numtasks*offset+taskid+1
c                k1=numtasks*offset+taskid+1
c                moli1=mol3new(k1)
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
               do moli1 =start,start+offset-1
c               do k1=start,start+offset-1
c                  moli1=mol3new(k1)

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
                  write (iout,11)
   11            format (/,' GRADIENT  --  Illegal Value for the Total',
     &              ' Potential Energy')
                 call fatal
              end if

           end if

           call mpi_bcast(derivs,3*n,mpi_real8,master,
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


         else if (integrate .eq. 'GHMC') then
            call ghmcstep (istep,dt)
         else if (integrate .eq. 'RIGIDBODY') then
            call rgdstep (istep,dt)
         else if (integrate .eq. 'RESPA') then
            call respa (istep,dt)
         else
c            call beeman (istep,dt)
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
           call mpi_barrier(mpi_comm_world,ierr)

c      print*,"After temper in beeman"
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
           call mpi_barrier(mpi_comm_world,ierr)
 
c      print*,"After rattle in beeman"
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
c         call mollist2body
c         call mollist3body
      end if

c         call mpi_bcast(nmollst,nmol,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c          call mpi_bcast(mollst,100*nmol,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c         call mpi_bcast(nmollst3mod,nmol,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c         call mpi_bcast(mollst3mod,800*nmol,mpi_integer,master,
c     &   mpi_comm_world,ierr)
c         print*,"After mollst bcast"
c         call mpi_barrier(mpi_comm_world,ierr)
c        nonzerocount3=0
c         do k=1,nmol
c          if(nmollst3mod(k).ne.0) then
c            nonzerocount3=nonzerocount3+1
c            mol3new(nonzerocount3)=k
c          end if
c         end do

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
c         print*,"After rpole mpi_bcast in beeman"
      switchmode = 'MPOLE'
      call switch (switchmode)
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


             offset=int((nmol-1)/numtasks)
             remainder=mod((nmol-1),numtasks)

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

             offset=int((nmol-2)/numtasks)
             remainder=mod((nmol-2),numtasks)

             if(taskid.le.remainder-1) then
               start=taskid*offset+1
               do moli1 =start,start+offset-1
c               do k1 =start,start+offset-1
c                 moli1=mol3new(k1)
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

                moli1=numtasks*offset+taskid+1
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
               do moli1 =start,start+offset-1

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
             print*,"Ep in master task in beeman=",ep
             print*,"Em in master task in beeman=",em
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
                  write (iout,13)
   13            format (/,' GRADIENT  --  Illegal Value for the Total',
     &              ' Potential Energy')
                 call fatal
              end if

           end if

           call mpi_bcast(derivs,3*n,mpi_real8,master,
     &   mpi_comm_world,ierr)
           call mpi_bcast(vir,3*3,mpi_real8,master,
     &   mpi_comm_world,ierr)
           call mpi_bcast(epot,1,mpi_real8,master,
     &   mpi_comm_world,ierr)
c          print*,"After derivs mpi_bcast in beeman"
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

           call mpi_barrier(mpi_comm_world,ierr)

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


         end if
      end do
c
c     perform any final tasks before program exit
c
      call final
      if(taskid.eq.master) then
      call gettime (wall,cpu)
      print*,"Wall Time=", wall, "CPU Time=",cpu
      end if

      deallocate(dep3bt)
      deallocate(dep3bt_tot2)
      deallocate(dep3bt_tot3)
      deallocate(dep3b2)
      deallocate(dep3b3)
      call mpi_finalize(ierr)

      end
