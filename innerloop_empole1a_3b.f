c
c     ###############################################################
c     ##                                                           ##
c     ##               Subroutine empole1c_3b                      ##
c     ##  Liam O-Suilleabhain, Omar Demerdash, Teresa Head-Gordon  ##
c     ##                 Spring 2013                               ##
c     ###############################################################
c
c
c     "empole1c_3b" calculates the atomic multipole and dipole
c     polarizability interaction energy using the 3-body approximation
c
c
      subroutine Innerloop1(moli1,ep3bt,virep3bt,dep3bt) 
      implicit none
      include 'sizes.i'
      include 'energi.i'
      include 'molcul.i'
      include 'deriv.i'
      include 'mpole.i'
      include 'virial.i'
      include 'mpif.h'
      include 'neigh.i'
      real*8  ep1moli1,ep1moli2,ep2moli12,ep1moli3
      real*8  dep1moli1(3,30),dep1moli2(3,30)
      real*8  dep2moli12(3,30),dep1moli3(3,30)
      integer i,ii,j,l1,i1,i2,i3,k
      integer k1,k2,k5
      real*8 eptemp,deptemp(3,npole)
      real*8 vir1moli1(3,3),vir2moli12(3,3),vir1moli3(3,3)
      real*8 vir1moli2(3,3),virtemp(3,3)
      integer pnum(30),npole3b,moli1,moli2,moli3
      real*8 ep3bt,dep3bt(3,*),virep3bt(3,3)
      logical do2
                ep3bt=0.0d0
c                ep3b=0.0d0
                do i = 1, npole
                   do j = 1, 3
                      dep3bt(j,i) = 0.0d0
c                      dep3b(j,i) = 0.0d0
                   end do
                end do

                do i=1,3
                   do j=1,3
c                      virep3b(i,j)=0.0d0
                      virep3bt(i,j)=0.0d0
                   end do
                end do

      do i = 1, 30
        do j = 1, 3
          dep1moli1(j,i)=0.0d0
          dep2moli12(j,i)=0.0d0
          dep1moli2(j,i)=0.0d0
        end do
      end do

      do i=1,3
         do j=1,3
           vir1moli1(i,j)=0.0d0
           vir2moli12(i,j)=0.0d0
           vir1moli2(i,j)=0.0d0
         end do
      end do 

        do k1=1,nmollst(moli1)
          moli2=mollst(k1,moli1)
          
          pnum(1)=imol(1,moli1)
          pnum(2)=imol(1,moli1)+1
          pnum(3)=imol(2,moli1)
          npole3b=6
          pnum(4)=imol(1,moli2)
          pnum(5)=imol(1,moli2)+1
          pnum(6)=imol(2,moli2)
        call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)

          ep3bt = ep3bt + eptemp
          ep2moli12=eptemp
          do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
              dep3bt(j,i) = dep3bt(j,i)+deptemp(j,i)
              dep2moli12(j,l1)=deptemp(j,i)
            end do
          end do
          do i=1,3
             do j=1,3
                virep3bt(i,j)=virep3bt(i,j)+virtemp(i,j)
                vir2moli12(i,j)=virtemp(i,j)
             end do
          end do

c          do k2=1,nmollst3mod2(moli2,moli1) 
c            done(moli2)=.true.
          do k2=1,nmollst(moli2)
             moli3=mollst(k2,moli2)
c            moli3=mollst3mod2(k2,moli2,moli1)

            npole3b=9
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli2)
            pnum(5)=imol(1,moli2)+1
            pnum(6)=imol(2,moli2)
            pnum(7)=imol(1,moli3)
            pnum(8)=imol(1,moli3)+1
            pnum(9)=imol(2,moli3)
            call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
            ep3bt = ep3bt + eptemp
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
                dep3bt(j,i) = dep3bt(j,i)+deptemp(j,i)
              end do
            end do
            do i=1,3
               do j=1,3
                 virep3bt(i,j)=virep3bt(i,j)+virtemp(i,j) 
               end do
            end do
            
            npole3b=6
            ep3bt = ep3bt - ep2moli12
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
                dep3bt(j,i) = dep3bt(j,i)-dep2moli12(j,l1)
              end do
            end do
            do i=1,3
               do j=1,3
                virep3bt(i,j)=virep3bt(i,j)-vir2moli12(i,j)
               end do
            end do

            npole3b=6
            pnum(1)=imol(1,moli2)
            pnum(2)=imol(1,moli2)+1
            pnum(3)=imol(2,moli2)
            pnum(4)=imol(1,moli3)
            pnum(5)=imol(1,moli3)+1
            pnum(6)=imol(2,moli3)
        call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
            ep3bt = ep3bt - eptemp
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
                dep3bt(j,i) = dep3bt(j,i)-deptemp(j,i)
              end do
            end do
            do i=1,3
               do j=1,3
                virep3bt(i,j)=virep3bt(i,j)-virtemp(i,j)
               end do
            end do

            do2=.false.
            do k5=1,nmollst(moli1)
               if(mollst(k5,moli1).eq.moli3) then
                 do2=.true.
                 goto 31
               end if
            end do

   31 continue
          if(do2.eq..true.) then
            npole3b=6
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli3)
            pnum(5)=imol(1,moli3)+1
            pnum(6)=imol(2,moli3)
            call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
            ep3bt = ep3bt - eptemp
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
                dep3bt(j,i) = dep3bt(j,i)-deptemp(j,i)
              end do
            end do
            do i=1,3
               do j=1,3
                virep3bt(i,j)=virep3bt(i,j)-virtemp(i,j)
               end do
            end do
          end if


          end do
        end do 
 
c      print*,"InnerloopMoli1= ep3bt=",moli1,ep3bt
      return
      end


      subroutine getlocalsum(moli1,localsum)
      integer localsum,moli1

      localsum=0
      localsum=localsum+1
      print*,"Moli1 localsum",moli1,localsum
      return
      end

c     #############################################################
c     ##                                                         ##
c     ##  subroutine Innerloop2                                  ##
c     ##                                                         ##
c     #############################################################
c
c     Innerloop2 calculates the 2-body polarization energy,gradient, and
c     virial for a single iteration using a neighbor list.


      subroutine Innerloop2(moli1,ep3bt,virep3bt,dep3bt)
      implicit none
      include 'sizes.i'
      include 'energi.i'
      include 'atoms.i'
      include 'atmtyp.i'
      include 'molcul.i'
      include 'deriv.i'
      include 'mpole.i'
      include 'virial.i'
      include 'mpif.h'
c      include 'neigh.i'
      integer i,ii,j,l1,i1,i2,i3,k
      integer k1,k2,k5,np1,np2
      real*8 eptemp,deptemp(3,npole)
      real*8 virtemp(3,3)
      integer pnum(30),npole3b,moli1,moli2,moli3
      real*8 ep3bt,dep3bt(3,*),virep3bt(3,3)
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2
      real*8 xcm3,ycm3,zcm3,xr,yr,zr,r_123
      real*8 M1,M2,M3,xr1,yr1,zr1,xr2,yr2,zr2,xr3,yr3,zr3
      logical do2
                ep3bt=0.0d0
                do i = 1, npole
                   do j = 1, 3
                      dep3bt(j,i) = 0.0d0
                   end do
                end do

                do i=1,3
                   do j=1,3
                      virep3bt(i,j)=0.0d0
                   end do
                end do
c        do k1=1,nmollst(moli1)
c          moli2=mollst(k1,moli1)
        do moli2=moli1+1,nmol
          np1=3
          np2=6
          pnum(1)=imol(1,moli1)
          pnum(2)=imol(1,moli1)+1
          pnum(3)=imol(2,moli1)
          npole3b=6
          pnum(4)=imol(1,moli2)
          pnum(5)=imol(1,moli2)+1
          pnum(6)=imol(2,moli2)
          xcm1 = 0.0d0
          ycm1 = 0.0d0
          zcm1 = 0.0d0
          xcm2 = 0.0d0
          ycm2 = 0.0d0
          zcm2 = 0.0d0
          M1 = 0.0d0
          M2 = 0.0d0
          do l1 = 1, np1
            i = pnum(l1)
            xcm1 = xcm1 + (mass(i)*x(i))
            ycm1 = ycm1 + (mass(i)*y(i))
            zcm1 = zcm1 + (mass(i)*z(i))
            M1 = M1 + mass(i)
          end do
          do l1 = np1+1, np2
            i = pnum(l1)
            xcm2 = xcm2 + (mass(i)*x(i))
            ycm2 = ycm2 + (mass(i)*y(i))
            zcm2 = zcm2 + (mass(i)*z(i))
            M2 = M2 + mass(i)
          end do
          xcm1 = xcm1/M1
          ycm1 = ycm1/M1
          zcm1 = zcm1/M1
          xcm2 = xcm2/M2
          ycm2 = ycm2/M2
          zcm2 = zcm2/M2

          xr = xcm1 - xcm2
          yr = ycm1 - ycm2
          zr = zcm1 - zcm2

          call image(xr,yr,zr)
          r_123=xr*xr + yr*yr + zr*zr
         if(r_123 .lt. 25.0d0) then
        call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
          ep3bt = ep3bt + eptemp
          do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
              dep3bt(j,i) = dep3bt(j,i)+deptemp(j,i)
            end do
          end do
          do i=1,3
             do j=1,3
                virep3bt(i,j)=virep3bt(i,j)+virtemp(i,j)
             end do
          end do
         end if
        end do   
      return
      end

c     #############################################################
c     ##                                                         ##
c     ##  subroutine Innerloop3                                  ##
c     ##                                                         ##
c     #############################################################
c
c     Innerloop3 calculates the 3-body polarization energy,gradient, and
c     virial for a single iteration using a neighbor list.

      subroutine Innerloop3(moli1,ep3bt,virep3bt,dep3bt)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'atmtyp.i' 
      include 'energi.i'
      include 'molcul.i'
      include 'deriv.i'
      include 'mpole.i'
      include 'virial.i'
      include 'mpif.h'
      include 'neigh.i'
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2
      real*8 xcm3,ycm3,zcm3,xr,yr,zr,r_123
      real*8 M1,M2,M3,xr1,yr1,zr1,xr2,yr2,zr2,xr3,yr3,zr3
      integer np1,np2,np3
      real*8 x1,y1,z1,x2,y2,z2
      real*8 x3,y3,z3,r1,r2,r3
      real*8 shellsum,shellsum_mod
      integer i,ii,j,l1,i1,i2,i3,k
      integer k1,k2,k5
      real*8 eptemp,deptemp(3,npole)
      real*8 virtemp(3,3)
      integer pnum(30),npole3b,moli1,moli2,moli3
      real*8 ep3bt,dep3bt(3,*),virep3bt(3,3)
      logical do2
                ep3bt=0.0d0
                do i = 1, npole
                   do j = 1, 3
                      dep3bt(j,i) = 0.0d0
                   end do
                end do

                do i=1,3
                   do j=1,3
                      virep3bt(i,j)=0.0d0
                   end do
                end do
c       triplecount=0
c        do k1=1,nmollst3mod(moli1),2
c           moli2=mollst3mod(k1,moli1)
c           moli3=mollst3mod(k1+1,moli1)
c           triplecount=triplecount+1
      do moli2 = moli1+1,nmol-1 
        do moli3=moli2+1,nmol
            npole3b=9
              np1=3
              np2=6
              np3=9
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli2)
            pnum(5)=imol(1,moli2)+1
            pnum(6)=imol(2,moli2)
            pnum(7)=imol(1,moli3)
            pnum(8)=imol(1,moli3)+1
            pnum(9)=imol(2,moli3)
            x1 = 0.0d0
            y1 = 0.0d0
            z1 = 0.0d0
            x2 = 0.0d0
            y2 = 0.0d0
            z2 = 0.0d0
            x3 = 0.0d0
            y3 = 0.0d0
            z3 = 0.0d0
c
c Find center of mass of each molecule
c
            do l1 = 1, np1
              i = pnum(l1)
              if(name(i).eq.'O') then
                x1 = x(i)
                y1 = y(i)
                z1 = z(i)
              end if
            end do
            do l1 = np1+1, np2
              i = pnum(l1)
              if(name(i).eq.'O') then
                x2 = x(i)
                y2 = y(i)
                z2 = z(i)
              end if
            end do
            do l1 = np2+1, np3
              i = pnum(l1)
              if(name(i).eq.'O') then
                x3 = x(i)
                y3 = y(i)
                z3 = z(i)
              end if
            end do

            xr = x1 - x2
            yr = y1 - y2
            zr = z1 - z2
            call image(xr,yr,zr)
            xr1=xr
            yr1=yr
            zr1=zr
c           print*,"xr1 yr1 zr1",xr1,yr1,zr1,moli1,moli2,moli3
            xr = x1 - x3
            yr = y1 - y3
            zr = z1 - z3
            call image(xr,yr,zr)
            xr2=xr
            yr2=yr
            zr2=zr
c           print*,"xr2 yr2 zr2",xr2,yr2,zr2,moli1,moli2,moli3
            xr = x2 - x3
            yr = y2 - y3
            zr = z2 - z3
            call image(xr,yr,zr)

            xr3=xr
            yr3=yr
            zr3=zr

            r1=sqrt(xr1*xr1 + yr1*yr1 + zr1*zr1)
            r2=sqrt(xr2*xr2 + yr2*yr2 + zr2*zr2)
            r3=sqrt(xr3*xr3 + yr3*yr3 + zr3*zr3)

            if( (r1.lt.r3).and.(r2.lt.r3) ) then
c              shell1=int(r1/3.1d0)+1
c              shell2=int(r2/3.1d0)+1
c              shell3=int(r3/3.1d0)+1
c        print*,"shell1,r,shell2,r",shell1,r1/3.1d0,shell2,r2/3.1d0
c              shellsum=shell1+shell2
              shellsum=r1+r2
              shellsum_mod=r1+r2+r3
c              shellsum_mod=shell1+shell2+shell3
            else if ( (r1.lt.r2).and.(r3.lt.r2)) then
c              shell1=int(r1/3.1d0)+1
c              shell2=int(r3/3.1d0)+1
c              shell3=int(r2/3.1d0)+1
c        print*,"shell1,r,shell2,r",shell1,r1/3.1d0,shell2,r3/3.1d0
c              shellsum=shell1+shell2           
              shellsum=r1+r3
              shellsum_mod=r1+r2+r3
c              shellsum_mod=shell1+shell2+shell3   
            else if ( (r2.lt.r1).and.(r3.lt.r1)) then
c              shell1=int(r2/3.1d0)+1
c              shell2=int(r3/3.1d0)+1
c              shell3=int(r1/3.1d0)+1
c        print*,"shell1,r,shell2,r",shell1,r2/3.1d0,shell2,r3/3.1d0
c              shellsum=shell1+shell2
              shellsum=r2+r3
              shellsum_mod=r1+r2+r3
c              shellsum_mod=shell1+shell2+shell3
            end if

          if (shellsum .le. 8.0d0) then
            call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
            ep3bt = ep3bt + eptemp
            do l1 = 1, npole3b
              i = pnum(l1)
              do j = 1, 3
                dep3bt(j,i) = dep3bt(j,i)+deptemp(j,i)
              end do
            end do
            do i=1,3
               do j=1,3
                 virep3bt(i,j)=virep3bt(i,j)+virtemp(i,j)
               end do
            end do


c            do2=.false.
c            do k5=1,nmollst(moli1)
c               if(mollst(k5,moli1).eq.moli2) then
c                 do2=.true.
c                 goto 33
c               end if
c            end do

c   33 continue

c          if(do2.eq..true.) then
            np1=3
            np2=6
            npole3b=6
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli2)
            pnum(5)=imol(1,moli2)+1
            pnum(6)=imol(2,moli2)
            xcm1 = 0.0d0
            ycm1 = 0.0d0
            zcm1 = 0.0d0
            xcm2 = 0.0d0
            ycm2 = 0.0d0
            zcm2 = 0.0d0
            M1 = 0.0d0
            M2 = 0.0d0
            do l1 = 1, np1
              i = pnum(l1)
              xcm1 = xcm1 + (mass(i)*x(i))
              ycm1 = ycm1 + (mass(i)*y(i))
              zcm1 = zcm1 + (mass(i)*z(i))
              M1 = M1 + mass(i)
            end do
            do l1 = np1+1, np2
              i = pnum(l1)
              xcm2 = xcm2 + (mass(i)*x(i))
              ycm2 = ycm2 + (mass(i)*y(i))
              zcm2 = zcm2 + (mass(i)*z(i))
              M2 = M2 + mass(i)
            end do
            xcm1 = xcm1/M1
            ycm1 = ycm1/M1
            zcm1 = zcm1/M1
            xcm2 = xcm2/M2
            ycm2 = ycm2/M2
            zcm2 = zcm2/M2

            xr = xcm1 - xcm2
            yr = ycm1 - ycm2
            zr = zcm1 - zcm2

            call image(xr,yr,zr)
            r_123=xr*xr + yr*yr + zr*zr

            if(r_123 .lt. 25.0d0) then
             call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
              ep3bt=ep3bt-eptemp
              do l1 = 1, npole3b
                i = pnum(l1)
                do j = 1, 3
                 dep3bt(j,i) = dep3bt(j,i)-deptemp(j,i)
                end do
              end do
              do i=1,3
               do j=1,3
                virep3bt(i,j)=virep3bt(i,j)-virtemp(i,j)
               end do
              end do
            end if

c            do2=.false.
c            do k5=1,nmollst(moli2)
c               if(mollst(k5,moli2).eq.moli3) then
c                 do2=.true.
c                 goto 34
c               end if
c            end do

c   34 continue
c          if(do2.eq..true.) then
            np1=3
            np2=6
            npole3b=6
            pnum(1)=imol(1,moli2)
            pnum(2)=imol(1,moli2)+1
            pnum(3)=imol(2,moli2)
            pnum(4)=imol(1,moli3)
            pnum(5)=imol(1,moli3)+1
            pnum(6)=imol(2,moli3)
            xcm1 = 0.0d0
            ycm1 = 0.0d0
            zcm1 = 0.0d0
            xcm2 = 0.0d0
            ycm2 = 0.0d0
            zcm2 = 0.0d0
            M1 = 0.0d0
            M2 = 0.0d0
            do l1 = 1, np1
              i = pnum(l1)
              xcm1 = xcm1 + (mass(i)*x(i))
              ycm1 = ycm1 + (mass(i)*y(i))
              zcm1 = zcm1 + (mass(i)*z(i))
              M1 = M1 + mass(i)
            end do
            do l1 = np1+1, np2
              i = pnum(l1)
              xcm2 = xcm2 + (mass(i)*x(i))
              ycm2 = ycm2 + (mass(i)*y(i))
              zcm2 = zcm2 + (mass(i)*z(i))
              M2 = M2 + mass(i)
            end do
            xcm1 = xcm1/M1
            ycm1 = ycm1/M1
            zcm1 = zcm1/M1
            xcm2 = xcm2/M2
            ycm2 = ycm2/M2
            zcm2 = zcm2/M2

            xr = xcm1 - xcm2
            yr = ycm1 - ycm2
            zr = zcm1 - zcm2

            call image(xr,yr,zr)
            r_123=xr*xr + yr*yr + zr*zr
          
            if(r_123 .lt. 25.0d0) then
             call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
              ep3bt = ep3bt - eptemp
             do l1 = 1, npole3b
               i = pnum(l1)
               do j = 1, 3
                dep3bt(j,i) = dep3bt(j,i)-deptemp(j,i)
               end do
             end do
             do i=1,3
               do j=1,3
                virep3bt(i,j)=virep3bt(i,j)-virtemp(i,j)
               end do
             end do
            end if

c            do2=.false.
c            do k5=1,nmollst(moli1)
c               if(mollst(k5,moli1).eq.moli3) then
c                 do2=.true.
c                 goto 31
c               end if
c            end do

c   31 continue
c          if(do2.eq..true.) then
            np1=3
            np2=6
            npole3b=6
            pnum(1)=imol(1,moli1)
            pnum(2)=imol(1,moli1)+1
            pnum(3)=imol(2,moli1)
            pnum(4)=imol(1,moli3)
            pnum(5)=imol(1,moli3)+1
            pnum(6)=imol(2,moli3)
            xcm1 = 0.0d0
            ycm1 = 0.0d0
            zcm1 = 0.0d0
            xcm2 = 0.0d0
            ycm2 = 0.0d0
            zcm2 = 0.0d0
            M1 = 0.0d0
            M2 = 0.0d0
            do l1 = 1, np1
              i = pnum(l1)
              xcm1 = xcm1 + (mass(i)*x(i))
              ycm1 = ycm1 + (mass(i)*y(i))
              zcm1 = zcm1 + (mass(i)*z(i))
              M1 = M1 + mass(i)
            end do
            do l1 = np1+1, np2
              i = pnum(l1)
              xcm2 = xcm2 + (mass(i)*x(i))
              ycm2 = ycm2 + (mass(i)*y(i))
              zcm2 = zcm2 + (mass(i)*z(i))
              M2 = M2 + mass(i)
            end do
            xcm1 = xcm1/M1
            ycm1 = ycm1/M1
            zcm1 = zcm1/M1
            xcm2 = xcm2/M2
            ycm2 = ycm2/M2
            zcm2 = zcm2/M2
            
            xr = xcm1 - xcm2
            yr = ycm1 - ycm2
            zr = zcm1 - zcm2
            
            call image(xr,yr,zr)
            r_123=xr*xr + yr*yr + zr*zr
           
            if(r_123 .lt. 25.0d0) then
             call empole1a_3b_Polar(npole3b,pnum,eptemp,deptemp,virtemp)
             ep3bt = ep3bt - eptemp
             do l1 = 1, npole3b
               i = pnum(l1)
               do j = 1, 3
                 dep3bt(j,i) = dep3bt(j,i)-deptemp(j,i)
               end do
             end do
             do i=1,3
               do j=1,3
                virep3bt(i,j)=virep3bt(i,j)-virtemp(i,j)
               end do
             end do
            end if

          end if
        end do 
      end do
      return
      end
