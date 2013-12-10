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
      subroutine mollist3body
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'neigh.i'
      include 'molcul.i'
c      include 'combo.i'
      include 'atmtyp.i'
      integer i,j,k,k1,ii,l1,lenmol,listtemp(800),tempcount
      real*8 xi,yi,zi,x1,y1,z1
      real*8 xr,yr,zr,xr1,yr1,zr1,xr2,yr2,zr2,xr3,yr3,zr3,M1
      real*8 radius,r1,r3,r2,r2outer,molbuf_cobar
      real*8 mol_lbuffer,mol_lbuf2,shellsum
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
c      print*, "Hello from Mollist3body"

c         mbuf2 = (mpolecut+lbuffer)**2
c         mbufx = (mpolecut+2.0d0*lbuffer)**2


      molbuf_cobar=8.0d0
c      molbuf_cobar2=16.0d0
      mol_lbuffer=4.0d0
      mol_lbuf2=(0.5d0*mol_lbuffer)**2
c      molbuf_cobarx=256.0d0

      if (domollst3bod) then
         domollst3bod = .false.
c         if (octahedron) then
c            do i = 1, npole
c               call mbuild (i)
c            end do
c         else
c          print*, "Within Mollist3body if domollst3bod"
c            call molfull3body
c            call molfull3bodymod2
c            call molfull3bodymod3
            call molfull3bodycobar
c            call molfull3bodymod4nominimage
c         end if
         return
      end if
c
c     test each site for displacement exceeding half the buffer
c
      do i = 1, nmol
c         ii = ipole(i)
         lenmol=imol(2,i)-imol(1,i)+1

         do l1 = 1, lenmol
c             j = pnum(l1)
            k=l1-1
            j=imol(1,i)+k
           if(name(j).eq.'O') then
            x1 = x(j)
            y1 = y(j)
            z1 = z(j)
           end if
         end do

         xr = x1 - xmolold3(i)
         yr = y1 - ymolold3(i)
         zr = z1 - zmolold3(i)
         call imagen (xr,yr,zr)
         r2outer = xr*xr + yr*yr + zr*zr
         if (r2outer .ge. mol_lbuf2) then
c
c     store coordinates to reflect update of this site
c
            xmolold3(i) = x1
            ymolold3(i) = y1
            zmolold3(i) = z1
c
c     rebuild the higher numbered neighbors of this site
c
            j = 0
            do k = i+1, nmol-1
               xr = x1 - xmolold3(k)
               yr = y1 - ymolold3(k)
               zr = z1 - zmolold3(k)
               call imagen (xr,yr,zr)
               xr1=xr
               yr1=yr
               zr1=zr
               do k1=k+1, nmol
                   xr = x1 - xmolold3(k1)
                   yr = y1 - ymolold3(k1)
                   zr = z1 - zmolold3(k1)
                   call imagen (xr,yr,zr)
                   xr2=xr
                   yr2=yr
                   zr2=zr
                   xr = xmolold3(k) - xmolold3(k1)
                   yr = ymolold3(k) - ymolold3(k1)
                   zr = zmolold3(k) - zmolold3(k1)
                   call imagen (xr,yr,zr)
                   xr3=xr
                   yr3=yr
                   zr3=zr

                  r1=sqrt(xr1*xr1 + yr1*yr1 + zr1*zr1)
                  r2=sqrt(xr2*xr2 + yr2*yr2 + zr2*zr2)
                  r3=sqrt(xr3*xr3 + yr3*yr3 + zr3*zr3)

                  if( (r1.lt.r3).and.(r2.lt.r3) ) then
                     shellsum=r1+r2
                  else if ( (r1.lt.r2).and.(r3.lt.r2)) then
                     shellsum=r1+r3
                  else if ( (r2.lt.r1).and.(r3.lt.r1)) then
                     shellsum=r2+r3
                  end if

                  if (shellsum .le. molbuf_cobar) then
                       j =j+1
                       mollst3mod(j,i)=k
                       j =j+1
                       mollst3mod(j,i)=k1
                  end if
               end do  
            end do
            nmollst3mod(i) = j
c
c     adjust lists of lower numbered neighbors of this site
c
            do k = 1, i-1
               do k1=k+1,nmol
                 if (k1.ne.i) then
                   xr = x1 - xmolold3(k)
                   yr = y1 - ymolold3(k)
                   zr = z1 - zmolold3(k)
                   call imagen (xr,yr,zr)
                   xr1=xr
                   yr1=yr
                   zr1=zr
c                   r2 = xr*xr + yr*yr + zr*zr
                   xr = x1 - xmolold3(k1)
                   yr = y1 - ymolold3(k1)
                   zr = z1 - zmolold3(k1)
                   call imagen (xr,yr,zr)
                   xr2=xr
                   yr2=yr
                   zr2=zr

                   xr=xmolold3(k)-xmolold3(k1)
                   yr=ymolold3(k)-ymolold3(k1)
                   zr=zmolold3(k)-zmolold3(k1)
                   call imagen (xr,yr,zr)
                   xr3=xr
                   yr3=yr
                   zr3=zr
                
                  r1=sqrt(xr1*xr1 + yr1*yr1 + zr1*zr1)
                  r2=sqrt(xr2*xr2 + yr2*yr2 + zr2*zr2)
                  r3=sqrt(xr3*xr3 + yr3*yr3 + zr3*zr3)

                  if( (r1.lt.r3).and.(r2.lt.r3) ) then
                     shellsum=r1+r2
                  else if ( (r1.lt.r2).and.(r3.lt.r2)) then
                     shellsum=r1+r3
                  else if ( (r2.lt.r1).and.(r3.lt.r1)) then
                     shellsum=r2+r3
                  end if


      LEFT OFF HERE
c                   if (r2 .le. molbuf2) then
                   if(shellsum .le. molbuf_cobar) then
                     do j = 1, nmollst3mod(k),2
                       if(((mollst3mod(j-1,k).eq.k1).and.
     &(mollst3mod(j,k).eq.i)).or.((mollst3mod(j-1,k).eq.i).and.
     &                  (mollst3mod(j,k).eq.k1))) goto 20
                     end do
                     nmollst3mod(k) = nmollst3mod(k) + 2
                     if(k1.lt.i) then
                       mollst3mod(nmollst3mod(k)-1,k)=k1
                       mollst3mod(nmollst3mod(k),k)=i
                     else
                       mollst3mod(nmollst3mod(k)-1,k)=i
                       mollst3mod(nmollst3mod(k),k)=k1
                     end if
   20                continue
c                   else if (r2 .le. molbufx) then
                   else
                     do j = 1, nmollst3mod(k),2
                       if(((mollst3mod(j-1,k).eq.k1).and.
     &(mollst3mod(j,k).eq.i)).or.((mollst3mod(j-1,k).eq.i).and.
     &(mollst3mod(j,k).eq.k1))) then
                          mollst3mod(j-1,k)=0
                          mollst3mod(j,k)=0
c                         mollst(j,k) = mollst(nmollst(k),k)
c                         nmollst(k) = nmollst(k) - 1
                         goto 30
                       end if
                     end do
   30                continue
                   end if
                 end if
               end do

               do j=1,800
                 listtemp(j)=0
               end do 
               tempcount=0  
               do j=1,nmollst3mod(k)
                  if(mollst3mod(j,k).ne.0) then
                     tempcount=tempcount+1
                    listtemp(tempcount)=mollst3mod(j,k)
                  end if
               end do  
               nmollst3mod(k)=0
               do j=1,tempcount
                 nmollst3mod(k)=nmollst3mod(k)+1
                 mollst3mod(nmollst3mod(k),k)=listtemp(j)
               end do 
            end do
         end if
      end do
      return
      end

