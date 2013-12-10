c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine empole3c_3bi  --  double loop multipole analysis
c
c     ##                                                           ##
c     ###############################################################
c
c
c     "empole3c" calculates the atomic multipole and dipole
c     polarizability interaction energy using a particle mesh
c     Ewald summation and double loop, and partitions the energy
c     among the atoms
c
c
      subroutine empole3c_3bi_PolelecOnly
      implicit none
      include 'sizes.i'
      include 'energi.i'
      include 'atoms.i'
      include 'molcul.i'
      include 'polar.i'
      include 'combo.i'
      include 'mpole.i'
      real*8 em1(nmol),em2(nmol,nmol),em3(nmol,nmol,nmol)
      real*8 ep1(nmol),ep2(nmol,nmol),ep3(nmol,nmol,nmol)
      real*8 emliam,epliam
      real*8, allocatable :: field(:,:)
      real*8, allocatable :: fieldp(:,:)
      integer i,ii,b1,b2,b3,j,l1,counter(nmol,nmol),counts
      emliam = 0
      epliam = 0
      delta1 = 0
      delta2 = 0
      delta3 = 0
      b1 = 1
      b2 = 2
      b3 = 3


      allocate (field(3,npole))
      allocate (fieldp(3,npole))

      do i = 1, npole
         do j = 1, 3
c            uind(j,i) = 0.0d0
c            uinp(j,i) = 0.0d0
            field(j,i) = 0.0d0
            fieldp(j,i) = 0.0d0
         end do
      end do
      do i = 1, npole
         do j = 1, 3
            u1b(j,i) = 0.0d0
            u3b(j,i) = 0.0d0
            do count1 = 1, 100
               u_2b(j,i,count1) = 0.0d0
            end do
         end do
      end do

      call empole3c_PermelecOnly(field,fieldp)

      body1 = .true.
      if (nmol .ge. b1) then
         moli = 1
         do moli1 = 1, nmol
            call combo1
            call empole3c_3b_PolelecOnly(field,fieldp)
            ep1(moli) = ep
            epliam = epliam + ep1(moli)
            moli = moli + 1
         end do
      end if

      count1 = 1 
      count2 = 1
      count3 = 1
      counts = 1

      body1 = .false.
      body2 = .true.
      if (nmol .ge. b2) then
         moli = 1
         do moli1 = 1, nmol-1
            call combo1
            do moli2 = moli1+1, nmol
               call combo2
               call findr2
               call empole3c_3b_PolelecOnly(field,fieldp)
               ep2(moli1,moli2)=ep-ep1(moli1)-ep1(moli2)
               epliam = epliam + ep2(moli1,moli2)
               counts = counts + 1
               counter(moli1,moli2) = count1
               moli = moli + 1
            end do
         end do
      end if


      moli = 1
      body2 = .false.
      body3 = .true.
      if (nmol .ge. b3) then
         do moli1 = 1, nmol-2
            call combo1
            do moli2 = moli1+1, nmol-1
               call combo2
               do moli3 = moli2+1, nmol
                  count1 = counter(moli1,moli2)
                  count2 = counter(moli1,moli3)
                  count3 = counter(moli2,moli3)
                  call combo3
                  call findr3
                  call empole3c_3b_PolelecOnly(field,fieldp)
                  ep3(moli1,moli2,moli3)=ep-ep2(moli1,moli2)
     &                                - ep2(moli1,moli3)
     &                                - ep2(moli2,moli3)
     &                                - ep1(moli1)
     &                                - ep1(moli2)
     &                                - ep1(moli3)
c                emliam = emliam + em3(moli1,moli2,moli3) 
                  epliam = epliam + ep3(moli1,moli2,moli3)
c                ep3analyze(moli) = ep3(moli1,moli2,moli3)
                  moli = moli + 1
               end do
            end do
         end do
      end if
      print*,epliam
      em = em
      ep = epliam
      deallocate (field)
      deallocate (fieldp)
      return
      end

