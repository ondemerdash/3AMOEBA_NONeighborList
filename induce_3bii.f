c     #############################################################
c     ##      COPYRIGHT (C)2013 Liam D. O'Suilleabhain           ##
c     ##                   All Rights Reserved                   ##
c     #############################################################
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine induce_3bii  --  double loop multipole analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "induce1_3bii" calculates the induced dipole from 3 body effects
c
c
      subroutine induce_3bii
      implicit none
      include 'boxes.i'
      include 'sizes.i'
      include 'ewald.i'
      include 'math.i'
      include 'atoms.i'
      include 'mpole.i'
      include 'polar.i'
      include 'molcul.i'
      include 'combo.i'
c      real*8 u3b(3,maxatm),u2b(3,maxatm),u1b(3,maxatm)
      integer i,ii,j,l1,l2

      do i = 1, npole
         do j = 1, 3
            u1b(j,i) = 0.0d0
            u2b(j,i) = 0.0d0
            u3b(j,i) = 0.0d0
         end do
      end do

      body1 = .true.

       moli = 1
       do moli1 = 1, nmol
          call combo1
          call induce_3b

          do l1 = 1, np1
             i = pnum(l1)
             do j = 1, 3
                u1b(j,i) = uind(j,i)
                u3b(j,i) = u1b(j,i)
             end do
          end do
          moli = moli + 1
       end do

       body1 = .false.
       body2 = .true.

       moli = 1
       do moli1 = 1, nmol-1 
          call combo1
          do moli2 = moli1+1, nmol
             call combo2
             call induce_3b

             do l1 = 1, np1
                i = pnum(l1)
                do j = 1, 3
                   u_2b(j,i,moli2) = uind(j,i) - u1b(j,i)
                   u2b(j,i) = u_2b(j,i,moli2) 
                end do
             end do

             do l1 = np1+1, np2
                i = pnum(l1)
                do j = 1, 3
                   u_2b(j,i,moli1) = uind(j,i) - u1b(j,i)
                   u2b(j,i) = u_2b(j,i,moli1) 
                end do
             end do
             do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
                   u3b(j,i) = u3b(j,i) + u2b(j,i)
                end do
             end do
c             print*,moli
c             do l1 = 1, npole3b
c                i = pnum(l1)
c                do j = 1, 3
c                   print*,i,j,uind(j,i)
c                end do
c             end do


             moli = moli + 1
          end do
       end do
       body2 = .false.
       body3 = .true.

       do moli1 = 1, nmol-2
          call combo1
          do moli2 = moli1+1, nmol-1
             call combo2
             do moli3 = moli2+1, nmol
                call combo3
                call induce_3b

                do l1 = 1, np1
                   i = pnum(l1)
                   do j = 1, 3
                         u3b(j,i) = u3b(j,i) 
     &                          + uind(j,i) - u_2b(j,i,moli2)
     &                          - u_2b(j,i,moli3) - u1b(j,i)
                   end do
                end do

                do l1 = np1+1, np2
                   i = pnum(l1)
                   do j = 1, 3
                      u3b(j,i) = u3b(j,i) 
     &                          + uind(j,i) - u_2b(j,i,moli1)
     &                          - u_2b(j,i,moli3) - u1b(j,i)
                   end do
                end do

                do l1 = np2+1, np3
                   i = pnum(l1)
                   do j = 1, 3
                      u3b(j,i) = u3b(j,i) 
     &                          + uind(j,i) - u_2b(j,i,moli1)
     &                          - u_2b(j,i,moli2) - u1b(j,i)
                   end do
                end do
                moli = moli + 1
             end do
          end do
       end do

        do i = 1, npole
           do j = 1, 3
              uind(j,i) = u3b(j,i)
           end do
        end do

      return
      end 
