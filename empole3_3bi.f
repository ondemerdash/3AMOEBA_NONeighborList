c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine empole3a_3bi  --  double loop multipole analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "empole3a" calculates the atomic multipole and dipole
c     polarizability interaction energy using a double loop,
c     and partitions the energy among the atoms
c
c
      subroutine empole3a_3bi
      implicit none
      include 'sizes.i'
      include 'energi.i'
      include 'atoms.i'
      include 'molcul.i'
      include 'combo.i'
      real*8 em1(nmol),em2(nmol,nmol),em3(nmol,nmol,nmol)
      real*8 ep1(nmol),ep2(nmol,nmol),ep3(nmol,nmol,nmol)
      real*8 emliam,epliam,r_123
      integer i,ii,b1,b2,b3
      emliam = 0
      epliam = 0
      delta1 = 0
      delta2 = 0
      delta3 = 0
      b1 = 1
      b2 = 2
      b3 = 3

      body1 = .true.
      if (nmol .ge. b1) then
          moli = 1
          do moli1 = 1, nmol
             call combo1
             call empole3a_3b
             em1(moli) = em
             ep1(moli) = ep
             emliam = emliam + em1(moli)
             epliam = epliam + ep1(moli)
             moli = moli + 1
          end do
       end if

       body1 = .false.
       body2 = .true.
    
       r_123=0.0d0
       if (nmol .ge. b2) then
          moli = 1
          do moli1 = 1, nmol-1 
             call combo1
             do moli2 = moli1+1, nmol
                call combo2
                call findr2 (r_123)
                if(r_123 .lt. 2.5d0) then

                   call empole3a_3b
                   em2(moli1,moli2)=em-em1(moli1)-em1(moli2)
                   ep2(moli1,moli2)=ep-ep1(moli1)-ep1(moli2)
                   emliam = emliam + em2(moli1,moli2)
                   epliam = epliam + ep2(moli1,moli2)
c                ep2analyze(moli) = ep2(moli1,moli2)
                else 
                   em2(moli1,moli2)=0
                   ep2(moli1,moli2)=0                   
                end if  

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
                   call combo3
                   call findr3 (r_123)
                   if(r_123 .lt. 2.5d0) then

                   call empole3a_3b                   
                   em3(moli1,moli2,moli3)=em-em2(moli1,moli2)
     &                                   - em2(moli1,moli3)
     &                                   - em2(moli2,moli3)
     &                                   - em1(moli1)
     &                                   - em1(moli2)
     &                                   - em1(moli3)
                   ep3(moli1,moli2,moli3)=ep-ep2(moli1,moli2)
     &                                   - ep2(moli1,moli3)
     &                                   - ep2(moli2,moli3)
     &                                   - ep1(moli1)
     &                                   - ep1(moli2)
     &                                   - ep1(moli3)
                   emliam = emliam + em3(moli1,moli2,moli3) 
                   epliam = epliam + ep3(moli1,moli2,moli3) 
c                   ep3analyze(moli) = ep3(moli1,moli2,moli3)
                   end if
                   moli = moli + 1
                end do
             end do 
          end do
        end if
        print*, "em3  =", emliam
        print*, "ep3  =", epliam
        em = emliam
        ep = epliam

      return
      end


c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine empole3b_3bi  --  double loop multipole analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "empole3b" calculates the atomic multipole and dipole
c     polarizability interaction energy using a neighbor list,
c     and partitions the energy among the atoms
c
c
      subroutine empole3b_3bi
      implicit none
      include 'sizes.i'
      include 'energi.i'
      include 'atoms.i'
      include 'molcul.i'
      include 'combo.i'
      real*8 em1(nmol),em2(nmol,nmol),em3(nmol,nmol,nmol)
      real*8 ep1(nmol),ep2(nmol,nmol),ep3(nmol,nmol,nmol)
      real*8 emliam,epliam
      integer i,ii,b1,b2,b3
      emliam = 0
      epliam = 0
      delta1 = 0
      delta2 = 0
      delta3 = 0
      b1 = 1
      b2 = 2
      b3 = 3

      if (nmol .ge. b1) then
          moli = 1
          do moli1 = 1, nmol
             call combo1
             call empole3b_3b
             em1(moli) = em
             ep1(moli) = ep
             emliam = emliam + em1(moli)
             epliam = epliam + ep1(moli)
             moli = moli + 1
          end do
       end if

       if (nmol .ge. b2) then
          moli = 1
          do moli1 = 1, nmol 
             call combo1
             do moli2 = moli1+1, nmol
                call combo2
                call empole3b_3b
                em2(moli1,moli2)=em-em1(moli1)-em1(moli2)
                ep2(moli1,moli2)=ep-ep1(moli1)-ep1(moli2)
                emliam = emliam + em2(moli1,moli2)
                epliam = epliam + ep2(moli1,moli2)
c                ep2analyze(moli) = ep2(moli1,moli2)
                moli = moli + 1
             end do
          end do
       end if
       moli = 1
       if (nmol .ge. b3) then
          do moli1 = 1, nmol
             call combo1
             do moli2 = moli1+1, nmol
                call combo2
                do moli3 = moli2+1, nmol
                   call combo3
                   call empole3b_3b
                   em3(moli1,moli2,moli3)=em-em2(moli1,moli2)
     &                                   - em2(moli1,moli3)
     &                                   - em2(moli2,moli3)
     &                                   - em1(moli1)
     &                                   - em1(moli2)
     &                                   - em1(moli3)
                   ep3(moli1,moli2,moli3)=ep-ep2(moli1,moli2)
     &                                   - ep2(moli1,moli3)
     &                                   - ep2(moli2,moli3)
     &                                   - ep1(moli1)
     &                                   - ep1(moli2)
     &                                   - ep1(moli3)
                   emliam = emliam + em3(moli1,moli2,moli3) 
                   epliam = epliam + ep3(moli1,moli2,moli3) 
c                   ep3analyze(moli) = ep3(moli1,moli2,moli3)
                   moli = moli + 1
                end do
             end do 
          end do
        end if
        print*, "em3  =", emliam
        print*, "ep3  =", epliam
        em = emliam
        ep = epliam

      return
      end


c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine empole3c_3bi  --  double loop multipole analysis  ##
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
      subroutine empole3c_3bi
      implicit none
      include 'sizes.i'
      include 'energi.i'
      include 'atoms.i'
      include 'molcul.i'
      include 'combo.i'
      real*8 em1(nmol),em2(nmol,nmol),em3(nmol,nmol,nmol)
      real*8 ep1(nmol),ep2(nmol,nmol),ep3(nmol,nmol,nmol)
      real*8 emliam,epliam
      real*8 r_123
      integer i,ii,b1,b2,b3
      emliam = 0
      epliam = 0
      delta1 = 0
      delta2 = 0
      delta3 = 0
      b1 = 1
      b2 = 2
      b3 = 3

      body1 = .true.
      if (nmol .ge. b1) then
         moli = 1
         do moli1 = 1, nmol
            call combo1
            call empole3c_3b
            em1(moli) = em
            ep1(moli) = ep
            emliam = emliam + em1(moli)
            epliam = epliam + ep1(moli)
            moli = moli + 1
         end do
      end if
     
      r_123 = 0.0d0

      body1 = .false.
      body2 = .true.
      if (nmol .ge. b2) then
         moli = 1
         do moli1 = 1, nmol-1 
            call combo1
            do moli2 = moli1+1, nmol
               call combo2
               r_123 = 0.0d0
               call findr2 (r_123)
                if(r_123 .lt. 2.5d0) then
                   call empole3c_3b
                   em2(moli1,moli2)=em-em1(moli1)-em1(moli2)
                   ep2(moli1,moli2)=ep-ep1(moli1)-ep1(moli2)
                   emliam = emliam + em2(moli1,moli2)
                   epliam = epliam + ep2(moli1,moli2)
c               ep2analyze(moli) = ep2(moli1,moli2)
                else
c                   em = 0.0d0
c                   ep = 0.0d0
c                   em2(moli1,moli2)=em-em1(moli1)-em1(moli2)
c                   ep2(moli1,moli2)=ep-ep1(moli1)-ep1(moli2)
                   em2(moli1,moli2)=0.0d0
                   ep2(moli1,moli2)=0.0d0
                   emliam = emliam + em2(moli1,moli2)
                   epliam = epliam + ep2(moli1,moli2)
                end if
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
                  call combo3
                  r_123 = 0.0d0
                  call findr3 (r_123)
                  if(r_123 .lt. 2.5d0) then
                    call empole3c_3b
                    em3(moli1,moli2,moli3)=em-em2(moli1,moli2)
     &                                   - em2(moli1,moli3)
     &                                   - em2(moli2,moli3)
     &                                   - em1(moli1)
     &                                   - em1(moli2)
     &                                   - em1(moli3)
                     ep3(moli1,moli2,moli3)=ep-ep2(moli1,moli2)
     &                                   - ep2(moli1,moli3)
     &                                   - ep2(moli2,moli3)
     &                                   - ep1(moli1)
     &                                   - ep1(moli2)
     &                                   - ep1(moli3)
                     emliam = emliam + em3(moli1,moli2,moli3) 
                     epliam = epliam + ep3(moli1,moli2,moli3) 
c                   ep3analyze(moli) = ep3(moli1,moli2,moli3)
                  else
c                     em = 0.0d0
c                     ep = 0.0d0
c                     em3(moli1,moli2,moli3)=em-em2(moli1,moli2)
c     &                                   - em2(moli1,moli3)
c     &                                   - em2(moli2,moli3)
c     &                                   - em1(moli1)
c     &                                   - em1(moli2)
c     &                                   - em1(moli3)
c                     ep3(moli1,moli2,moli3)=ep-ep2(moli1,moli2)
c     &                                   - ep2(moli1,moli3)
c     &                                   - ep2(moli2,moli3)
c     &                                   - ep1(moli1)
c     &                                   - ep1(moli2)
c     &                                   - ep1(moli3)
                    em3(moli1,moli2,moli3)=0.0d0
                    ep3(moli1,moli2,moli3)=0.0d0
                     emliam = emliam + em3(moli1,moli2,moli3)
                     epliam = epliam + ep3(moli1,moli2,moli3)
                  end if 
                   moli = moli + 1
               end do
            end do 
         end do
      end if
      print*, "em3  =", emliam
      print*, "ep3  =", epliam
      em = emliam
      ep = epliam

      return
      end


c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine empole3d_3bi  --  double loop multipole analysis  ##
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
      subroutine empole3d_3bi
      implicit none
      include 'sizes.i'
      include 'energi.i'
      include 'atoms.i'
      include 'molcul.i'
      include 'combo.i'
      real*8 em1(nmol),em2(nmol,nmol),em3(nmol,nmol,nmol)
      real*8 ep1(nmol),ep2(nmol,nmol),ep3(nmol,nmol,nmol)
      real*8 emliam,epliam
      integer i,ii,b1,b2,b3
      emliam = 0
      epliam = 0
      delta1 = 0
      delta2 = 0
      delta3 = 0
      b1 = 1
      b2 = 2
      b3 = 3

      if (nmol .ge. b1) then
          moli = 1
          do moli1 = 1, nmol
             call combo1
             call empole3d_3b
             em1(moli) = em
             ep1(moli) = ep
             emliam = emliam + em1(moli)
             epliam = epliam + ep1(moli)
             moli = moli + 1
          end do
       end if

       if (nmol .ge. b2) then
          moli = 1
          do moli1 = 1, nmol 
             call combo1
             do moli2 = moli1+1, nmol
                call combo2
                call empole3d_3b
                em2(moli1,moli2)=em-em1(moli1)-em1(moli2)
                ep2(moli1,moli2)=ep-ep1(moli1)-ep1(moli2)
                emliam = emliam + em2(moli1,moli2)
                epliam = epliam + ep2(moli1,moli2)
c                ep2analyze(moli) = ep2(moli1,moli2)
                moli = moli + 1
             end do
          end do
       end if
       moli = 1
       if (nmol .ge. b3) then
          do moli1 = 1, nmol
             call combo1
             do moli2 = moli1+1, nmol
                call combo2
                do moli3 = moli2+1, nmol
                   call combo3
                   call empole3d_3b
                   em3(moli1,moli2,moli3)=em-em2(moli1,moli2)
     &                                   - em2(moli1,moli3)
     &                                   - em2(moli2,moli3)
     &                                   - em1(moli1)
     &                                   - em1(moli2)
     &                                   - em1(moli3)
                   ep3(moli1,moli2,moli3)=ep-ep2(moli1,moli2)
     &                                   - ep2(moli1,moli3)
     &                                   - ep2(moli2,moli3)
     &                                   - ep1(moli1)
     &                                   - ep1(moli2)
     &                                   - ep1(moli3)
                   emliam = emliam + em3(moli1,moli2,moli3) 
                   epliam = epliam + ep3(moli1,moli2,moli3) 
c                   ep3analyze(moli) = ep3(moli1,moli2,moli3)
                   moli = moli + 1
                end do
             end do 
          end do
        end if
        print*, "em3  =", emliam
        print*, "ep3  =", epliam
        em = emliam
        ep = epliam

      return
      end
