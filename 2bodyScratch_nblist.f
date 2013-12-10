
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine mollist -- build atom multipole neighbor lists  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "mollist" performs an update or a complete rebuild of the
c     body/molecule neighbor lists for 3-body approx of polariz. energy
c
c
      subroutine mollist2body
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'molcul.i'
      include 'combo.i'
      include 'atmtyp.i'
      integer i,j,k,ii,l1,lenmol
      real*8 xi,yi,zi,xcm1,ycm1,zcm1
      real*8 xr,yr,zr,M1
      real*8 radius,r2
c
c
c     neighbor list cannot be used with the replicates method
c
c      radius = sqrt(mbuf2)
c      call replica (radius)
c      if (use_replica) then
c         write (iout,10)
c   10    format (/,' MLIST  --  Pairwise Neighbor List cannot',
c     &              ' be used with Replicas')
c         call fatal
c      end if
c
c     perform a complete list build instead of an update
c
c      print*, "Hello from Mollist2body"
         mbuf2 = (mpolecut+lbuffer)**2
         mbufx = (mpolecut+2.0d0*lbuffer)**2

      molbuf2=25.0d0
      mol_lbuffer=4.0d0
      mol_lbuf2=(0.5d0*mol_lbuffer)**2
      molbufx=81.0d0

      if (domollst2bod) then
         domollst2bod = .false.
c         if (octahedron) then
c            do i = 1, npole
c               call mbuild (i)
c            end do
c         else
c          print*, "Within Mollist2body if domollst2bod"
            call molfull2body
c         end if
         return
      end if
c
c     test each site for displacement exceeding half the buffer
c
      do i = 1, nmol
c         ii = ipole(i)

         xcm1 = 0.0d0
         ycm1 = 0.0d0
         zcm1 = 0.0d0
         M1 = 0.0d0
         lenmol=imol(2,i)-imol(1,i)+1

         do l1 = 1, lenmol
c             j = pnum(l1)
            k=l1-1
            j=imol(1,i)+k
            xcm1 = xcm1 + (mass(j)*x(j))
            ycm1 = ycm1 + (mass(j)*y(j))
            zcm1 = zcm1 + (mass(j)*z(j))
            M1 = M1 + mass(j)
         end do

         xcm1 = xcm1/M1
         ycm1 = ycm1/M1
         zcm1 = zcm1/M1

         xr = xcm1 - xmolold(i)
         yr = ycm1 - ymolold(i)
         zr = zcm1 - zmolold(i)
         call imagen (xr,yr,zr)
         r2 = xr*xr + yr*yr + zr*zr
         if (r2 .ge. mol_lbuf2) then
c
c     store coordinates to reflect update of this site
c
            xmolold(i) = xcm1
            ymolold(i) = ycm1
            zmolold(i) = zcm1
c
c     rebuild the higher numbered neighbors of this site
c
            j = 0
            do k = i+1, nmol
               xr = xcm1 - xmolold(k)
               yr = ycm1 - ymolold(k)
               zr = zcm1 - zmolold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. molbuf2) then
                  j = j + 1
                  mollst(j,i) = k
               end if
            end do
            nmollst(i) = j
c
c     adjust lists of lower numbered neighbors of this site
c
            do k = 1, i-1
               xr = xcm1 - xmolold(k)
               yr = ycm1 - ymolold(k)
               zr = zcm1 - zmolold(k)
               call imagen (xr,yr,zr)
               r2 = xr*xr + yr*yr + zr*zr
               if (r2 .le. molbuf2) then
                  do j = 1, nmollst(k)
                     if (mollst(j,k) .eq. i)  goto 20
                  end do
                  nmollst(k) = nmollst(k) + 1
                  mollst(nmollst(k),k) = i
   20             continue
               else if (r2 .le. molbufx) then
                  do j = 1, nmollst(k)
                     if (mollst(j,k) .eq. i) then
                        mollst(j,k) = mollst(nmollst(k),k)
                        nmollst(k) = nmollst(k) - 1
                        goto 30
                     end if
                  end do
   30             continue
               end if
            end do
         end if
      end do
      return
      end

