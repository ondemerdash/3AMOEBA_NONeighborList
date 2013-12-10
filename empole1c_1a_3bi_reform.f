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
      subroutine empole1c_1a_3b
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'atmtyp.i'
      include 'energi.i'
      include 'molcul.i'
      include 'deriv.i'
      include 'mpole.i'
      include 'virial.i'
      include 'openmp.i'
c      real*8, allocatable :: field(:,:)
c      real*8, allocatable :: fieldp(:,:)
c      real*8  ep3b,dep3b(3,npole)
      real*8  ep2moli13,ep2moli12,ep2moli23,ep3b,ep2b
      real*8  dep1moli1(3,30),dep1moli2(3,30),dep1moli3(3,30)
      real*8  dep2moli12(3,30)
      integer i,ii,j,l1,i1,i2,i3
c      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2
c      real*8 xcm3,ycm3,zcm3,xcm,ycm,zcm
c      real*8 M,M1,M2,M3,r_123
c      real*8 virep3b(3,3),eptemp,deptemp(3,30)

      real*8 eptemp,deptemp(3,npole)
c      real*8 eptemp,deptemp(3,30)
      real*8 vir1moli1(3,3),vir2moli12(3,3),vir1moli3(3,3)
      real*8 vir1moli2(3,3),virtemp(3,3)
      integer pnum(30),npole3b,moli1,moli2,moli3
      integer pnum2(npole),nmol2,np1,np2,np3 
      real*8 ep3bt,dep3bt(3,npole),virep3bt(3,3)
      integer omp_get_num_threads,nthread2,tid
      integer omp_get_thread_num,triplecount,paircount
      real*8 findr3_minimage2func,r_123,cutoff
      real*8 findr2_minimage2func,ewaldcut3b
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2
      real*8 xcm3,ycm3,zcm3,xr,yr,zr,off3b
      real*8 M1,M2,M3,xr1,yr1,zr1,xr2,yr2,zr2,xr3,yr3,zr3
      character*6 mode
      external findr3_minimage2func
      external findr2_minimage2func
      paircount=0 
      triplecount=0
c      ep3b = 0.0d0
      ep3bt=0.0d0 
      ep3b=0.0d0
      ep2b=0.0d0
c      allocate (field(3,npole))
c      allocate (fieldp(3,npole))
c      ncell=7
      do i = 1, npole
        do j = 1, 3
c          field(j,i) = 0.0d0
c          fieldp(j,i) = 0.0d0
c          dep3b(j,i) = 0.0d0
          dep3bt(j,i) = 0.0d0
c          dep1moli1(j,i)=0.0d0
c          dep2moli12(j,i)=0.0d0
c          dep1moli2(j,i)=0.0d0
        end do
      end do

c        do i=1,npole
c          do j=1,3
c            dep1moli1(j,i)=0.0d0
c          end do 
c        end do   

      do i=1,3
         do j=1,3
c           virep3b(i,j)=0.0d0
           virep3bt(i,j)=0.0d0
c           vir1moli1(i,j)=0.0d0
c           vir2moli12(i,j)=0.0d0
c           vir1moli2(i,j)=0.0d0
         end do
      end do 
c
c    Usual 2-body permanent electrostatics but without polarization
c

      call empole1c_3b_Perm
c      cutoff = 0.0d0
c      call replica (cutoff)
      mode = 'MPOLE'
      call switch (mode)

c      call just_chk_rotpole

c
c    Begin 2-body and 3-body approximation
c
ccccc !$OMP PARALLEL DO default(private) shared(nmol2,imol)
ccccc !$OMP& reduction(+:ep3bt,virep3bt,dep3bt) schedule(dynamic)

!$OMP PARALLEL default(private) shared(nmol,imol,ep3bt,
!$OMP& virep3bt,dep3bt,x,y,z,mass)
!$OMP DO reduction(+:ep3bt,virep3bt,dep3bt)
!$OMP& schedule(dynamic)
      do moli1 = 1, nmol
c
        do moli2 = moli1+1, nmol

          npole3b=6
          np1=3
          np2=6
          pnum(1)=imol(1,moli1)
          pnum(2)=imol(1,moli1)+1
          pnum(3)=imol(2,moli1)
          pnum(4)=imol(1,moli2)
          pnum(5)=imol(1,moli2)+1
          pnum(6)=imol(2,moli2)
c          print*,"pnum(6)",pnum(6)
c          r_123=findr2_minimage2func(pnum,np1,np2)
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
c            print*,"x(ii)",x(ii)
c            print*,"x(i)",x(i)
c            print*,"mass(i)",mass(i)
          end do
c          print*,"M1",M1
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
c          print*,"xcm1 xcm2",xcm1,xcm2
c          print*,"ycm1 ycm2",ycm1,ycm2
c          print*,"zcm1 zcm2",zcm1,zcm2
          
          xr = xcm1 - xcm2
          yr = ycm1 - ycm2
          zr = zcm1 - zcm2
c          print*,"xr yr zr before image",xr,yr,zr
          call image(xr,yr,zr)
          r_123=xr*xr + yr*yr + zr*zr
c         off3b=sqrt(xr*xr + yr*yr + zr*zr)
          off3b=r_123        
c         if(r_123 .lt. 25.0d0) then
c          print*,"moli1 moli2 r_123",moli1,moli2,r_123
c          print*,"xr yr zr",moli1,moli2,xr,yr,zr
          paircount=paircount+1
          pnum(1)=imol(1,moli1)
          pnum(2)=imol(1,moli1)+1
          pnum(3)=imol(2,moli1)

          npole3b=6
          pnum(4)=imol(1,moli2)
          pnum(5)=imol(1,moli2)+1
          pnum(6)=imol(2,moli2)

          call empole1a_3b_Polar(npole3b,pnum,eptemp,
     &    deptemp,virtemp)
          ep3bt = ep3bt + eptemp
          ep2b=ep2b+eptemp
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

c         else
c          npole3b=6
c          pnum(1)=imol(1,moli1)
c          pnum(2)=imol(1,moli1)+1
c          pnum(3)=imol(2,moli1)
c          pnum(4)=imol(1,moli2)
c          pnum(5)=imol(1,moli2)+1
c          pnum(6)=imol(2,moli2)
c          ep2moli12=0.0d0
c          do l1 = 1, npole3b
c            i = pnum(l1)
c            do j = 1, 3
c              dep2moli12(j,l1)=0.0d0
c            end do
c          end do
c          do i=1,3
c             do j=1,3
c                vir2moli12(i,j)=0.0d0
c             end do
c          end do
c         end if

          do moli3 = moli2+1, nmol
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
            xcm1 = 0.0d0
            ycm1 = 0.0d0
            zcm1 = 0.0d0
            xcm2 = 0.0d0
            ycm2 = 0.0d0
            zcm2 = 0.0d0
            xcm3 = 0.0d0
            ycm3 = 0.0d0
            zcm3 = 0.0d0
            M1 = 0.0d0
            M2 = 0.0d0
            M3 = 0.0d0
c
c Find center of mass of each molecule
c
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
            do l1 = np2+1, np3
              i = pnum(l1)
              xcm3 = xcm3 + (mass(i)*x(i))
              ycm3 = ycm3 + (mass(i)*y(i))
              zcm3 = zcm3 + (mass(i)*z(i))
              M3 = M3 + mass(i)
            end do
            xcm1 = xcm1/M1
            ycm1 = ycm1/M1
            zcm1 = zcm1/M1
            xcm2 = xcm2/M2
            ycm2 = ycm2/M2
            zcm2 = zcm2/M2
            xcm3 = xcm3/M3
            ycm3 = ycm3/M3
            zcm3 = zcm3/M3

            xr = xcm1 - xcm2
            yr = ycm1 - ycm2
            zr = zcm1 - zcm2
            call image(xr,yr,zr)
            xr1=xr
            yr1=yr
            zr1=zr
            xr = xcm1 - xcm3
            yr = ycm1 - ycm3
            zr = zcm1 - zcm3
            call image(xr,yr,zr)
            xr2=xr
            yr2=yr
            zr2=zr
            xr = xcm2 - xcm3
            yr = ycm2 - ycm3
            zr = zcm2 - zcm3
            call image(xr,yr,zr)
 
            xr3=xr
            yr3=yr
            zr3=zr

           r_123 = (xr1*xr1 + yr1*yr1 + zr1*zr1+
     &        xr2*xr2 + yr2*yr2 + zr2*zr2 +
     &        xr3*xr3 + yr3*yr3 + zr3*zr3)/3.0d0
c           off3b=(sqrt(xr1*xr1 + yr1*yr1 + zr1*zr1)+
c     &     sqrt(xr2*xr2 + yr2*yr2 + zr2*zr2)+
c     &     sqrt(xr3*xr3 + yr3*yr3 + zr3*zr3))/3.0d0
c           print*,"r_123",r_123,moli1,moli2,moli3
           off3b=r_123
c           if(r_123 .lt. 25.0d0) then
            triplecount=triplecount+1
c           print*,"moli1 moli2 moli3 r_123",moli1,moli2,moli3,r_123
c           print*,"moli1 moli2 moli3 xr1 yr1,zr1",moli1,moli2,moli3,xr1,
c     &      yr1,zr1
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
c            ewaldcut3b=r_123+10.0d0
            call empole1a_3b_Polar(npole3b,pnum,eptemp,
     &      deptemp,virtemp)
            ep3b=ep3b+eptemp
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


c            if(ep2moli12.ne.0.0d0) then
              pnum(1)=imol(1,moli1)
              pnum(2)=imol(1,moli1)+1
              pnum(3)=imol(2,moli1)
              pnum(4)=imol(1,moli2)
              pnum(5)=imol(1,moli2)+1
              pnum(6)=imol(2,moli2)

              npole3b=6
              ep3bt = ep3bt - ep2moli12
              ep3b = ep3b-ep2moli12
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
c            end if

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
c            off3b=sqrt(xr*xr + yr*yr + zr*zr)
            off3b=r_123            
c            if(r_123 .lt. 25.0d0) then
              npole3b=6
              pnum(1)=imol(1,moli1)
              pnum(2)=imol(1,moli1)+1
              pnum(3)=imol(2,moli1)
              pnum(4)=imol(1,moli3)
              pnum(5)=imol(1,moli3)+1
              pnum(6)=imol(2,moli3)
c              ewaldcut3b=r_123+10.0d0
              call empole1a_3b_Polar(npole3b,pnum,eptemp,
     &        deptemp,virtemp)
              ep3bt = ep3bt - eptemp
              ep3b=ep3b-eptemp
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
c            end if  
            np1=3
            np2=6
            npole3b=6
            pnum(1)=imol(1,moli2)
            pnum(2)=imol(1,moli2)+1
            pnum(3)=imol(2,moli2)
            pnum(4)=imol(1,moli3)
            pnum(5)=imol(1,moli3)+1
            pnum(6)=imol(2,moli3)
c            r_123=findr2_minimage2func(pnum,np1,np2)

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
c            off3b=sqrt(xr*xr + yr*yr + zr*zr)
            off3b=r_123          
c            if(r_123 .lt. 25.0d0) then
              npole3b=6
              pnum(1)=imol(1,moli2)
              pnum(2)=imol(1,moli2)+1
              pnum(3)=imol(2,moli2)
              pnum(4)=imol(1,moli3)
              pnum(5)=imol(1,moli3)+1
              pnum(6)=imol(2,moli3)

c              ewaldcut3b=r_123+10.0d0
              call empole1a_3b_Polar(npole3b,pnum,eptemp,
     &        deptemp,virtemp)
              ep3bt = ep3bt - eptemp
              ep3b=ep3b-eptemp
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

c            end if

c           end if

          end do
        end do
      end do
!$OMP END DO
!$OMP END PARALLEL


cccccc!$OMP END PARALLEL DO

ccc 100  continue
c      ep3b=ep3b+ep3bt
      ep = ep+ep3bt
      print*,"2body poleng",ep2b
      print*,"3body poleng",ep3b

      write(6,*)ep
      print*,"Triplecount=",triplecount
      print*,"Paircount=",paircount
      do i = 1, npole
         do j = 1, 3
            dep(j,i) =dep(j,i)+dep3bt(j,i)
         end do
      end do

      do i=1,3
         do j=1,3
           vir(i,j)=vir(i,j)+virep3bt(i,j)
         end do
      end do 
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine empole1a_3b_Polar                           ##
c     ##                                                         ##
c     #############################################################
c
c     enoike1a_3b_Polar calculates the polarization energy,gradient, and
c     virial using ordinary Coulomb electrostatics with an infinite
c     cutoff


      subroutine empole1a_3b_Polar(npole3b,pnum,eptemp,
     &   deptemp,virtemp)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'bound.i'
      include 'boxes.i'
      include 'cell.i'
      include 'chgpot.i'
      include 'couple.i'
      include 'cutoff.i'
c      include 'deriv.i'
c      include 'energi.i'
      include 'group.i'
      include 'inter.i'
      include 'molcul.i'
      include 'mplpot.i'
      include 'mpole.i'
c      include 'polar.i'
      include 'polar2.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'shunt.i'
      include 'usage.i'
c      include 'virial.i'
      integer i,j,k
      integer ii,kk,jcell
      integer ix,iy,iz
      integer kx,ky,kz
      integer iax,iay,iaz
      integer kax,kay,kaz
      real*8 e,ei,f,fgrp,gfd
      real*8 damp,expdamp
      real*8 pdi,pti,pgamma
      real*8 scale3,scale3i
      real*8 scale5,scale5i
      real*8 scale7,scale7i
      real*8 temp3,temp5,temp7
      real*8 psc3,psc5,psc7
      real*8 dsc3,dsc5,dsc7
      real*8 xr,yr,zr
      real*8 xix,yix,zix
      real*8 xiy,yiy,ziy
      real*8 xiz,yiz,ziz
      real*8 xkx,ykx,zkx
      real*8 xky,yky,zky
      real*8 xkz,ykz,zkz
      real*8 r,r2,rr1,rr3
      real*8 rr5,rr7,rr9,rr11
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 ci,di(3),qi(9)
      real*8 ck,dk(3),qk(9)
      real*8 frcxi(3),frcxk(3)
      real*8 frcyi(3),frcyk(3)
      real*8 frczi(3),frczk(3)
      real*8 fridmp(3),findmp(3)
      real*8 ftm2(3),ftm2i(3)
      real*8 ttm2(3),ttm3(3)
      real*8 ttm2i(3),ttm3i(3)
      real*8 dixdk(3),fdir(3)
      real*8 dixuk(3),dkxui(3)
      real*8 dixukp(3),dkxuip(3)
      real*8 uixqkr(3),ukxqir(3)
      real*8 uixqkrp(3),ukxqirp(3)
      real*8 qiuk(3),qkui(3)
      real*8 qiukp(3),qkuip(3)
      real*8 rxqiuk(3),rxqkui(3)
      real*8 rxqiukp(3),rxqkuip(3)
      real*8 qidk(3),qkdi(3)
      real*8 qir(3),qkr(3)
      real*8 qiqkr(3),qkqir(3)
      real*8 qixqk(3),rxqir(3)
      real*8 dixr(3),dkxr(3)
      real*8 dixqkr(3),dkxqir(3)
      real*8 rxqkr(3),qkrxqir(3)
      real*8 rxqikr(3),rxqkir(3)
      real*8 rxqidk(3),rxqkdi(3)
      real*8 ddsc3(3),ddsc5(3)
      real*8 ddsc7(3)
      real*8 gl(0:8),gli(7),glip(7)
      real*8 sc(10),sci(8),scip(8)
      real*8 gf(7),gfi(6),gti(6)
      real*8, allocatable :: mscale(:)
      real*8, allocatable :: pscale(:)
      real*8, allocatable :: dscale(:)
      real*8, allocatable :: uscale(:)
      real*8 eptemp,deptemp(3,npole),virtemp(3,3)
      real*8 uind(3,npole3b)
      real*8 uinp(3,npole3b)
      real*8 off3b,eptemp_userep
      integer npole3b,pnum(*),l1,l3
      logical proceed,usei,usek
      character*6 mode

      eptemp = 0.0d0
      eptemp_userep =0.0d0

c   Zero out temporary gradient of polarization energy
      do l1 = 1, npole3b
        i = pnum(l1)
        do j = 1, 3
          deptemp(j,i) = 0.0d0
        end do
      end do

c   Zero out temporary virial
      do i=1,3
         do j=1,3
           virtemp(i,j)=0.0d0
         end do
      end do


      call induce0a_3b_PolelecOnly(npole3b,pnum,uind,uinp)

c      do l1 = 1, npole3b
c        i = pnum(l1)
c         print*,"uindx",uind(1,l1)
c         print*,"uindy",uind(2,l1)
c         print*,"uindz",uind(3,l1)
c      end do 
          

c     perform dynamic allocation of some local arrays
c
      allocate (mscale(n))
      allocate (pscale(n))
      allocate (dscale(n))
      allocate (uscale(n))
c
c     set arrays needed to scale connected atom interactions
c
      if (npole .eq. 0)  return
      do i = 1, n
         mscale(i) = 1.0d0
         pscale(i) = 1.0d0
         dscale(i) = 1.0d0
         uscale(i) = 1.0d0
      end do
c
c     set conversion factor, cutoff and switching coefficients
c
      f = electric / dielec
c      mode = 'MPOLE'
c      call switch (mode)

c
c     set scale factors for permanent multipole and induced terms
c
c      if(use_standard) then
c      if(use_replica.eq..false.) then

      do l1 = 1, npole3b-1
         i=pnum(l1)
         ii = ipole(i)
         iz = zaxis(i)
         ix = xaxis(i)
         iy = yaxis(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         di(1) = rpole(2,i)
         di(2) = rpole(3,i)
         di(3) = rpole(4,i)
         qi(1) = rpole(5,i)
         qi(2) = rpole(6,i)
         qi(3) = rpole(7,i)
         qi(4) = rpole(8,i)
         qi(5) = rpole(9,i)
         qi(6) = rpole(10,i)
         qi(7) = rpole(11,i)
         qi(8) = rpole(12,i)
         qi(9) = rpole(13,i)
         usei = (use(ii) .or. use(iz) .or. use(ix) .or. use(iy))
c
c     set interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = m2scale
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = m3scale
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = m4scale
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
                if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = m5scale
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = d1scale
            uscale(ip11(j,ii)) = u1scale
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = d2scale
            uscale(ip12(j,ii)) = u2scale
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = d3scale
            uscale(ip13(j,ii)) = u3scale
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = d4scale
            uscale(ip14(j,ii)) = u4scale
         end do
         do l3 = l1+1, npole3b
            k=pnum(l3)
            kk = ipole(k)
            kz = zaxis(k)
            kx = xaxis(k)
            ky = yaxis(k)
            usek = (use(kk) .or. use(kz) .or. use(kx) .or. use(ky))
            proceed = .true.
            if (use_group)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (.not. use_intra)  proceed = .true.
            if (proceed)  proceed = (usei .or. usek)
            if (.not. proceed)  goto 10
            xr = x(kk) - x(ii)
            yr = y(kk) - y(ii)
            zr = z(kk) - z(ii)
            call image (xr,yr,zr)
            r2 = xr*xr + yr*yr + zr*zr
c            print*,"r2 in stdloop empole1a_3b_Polar",r2
c            print*,"xr yr zr in stdloop empole1a_3b_Polar",xr,yr,zr

c            if (r2 .le. off2) then
               r = sqrt(r2)
               ck = rpole(1,k)
               dk(1) = rpole(2,k)
               dk(2) = rpole(3,k)
               dk(3) = rpole(4,k)
               qk(1) = rpole(5,k)
               qk(2) = rpole(6,k)
               qk(3) = rpole(7,k)
               qk(4) = rpole(8,k)
               qk(5) = rpole(9,k)
               qk(6) = rpole(10,k)
               qk(7) = rpole(11,k)
               qk(8) = rpole(12,k)
               qk(9) = rpole(13,k)
c
c     apply Thole polarization damping to scale factors
c
               rr1 = 1.0d0 / r
               rr3 = rr1 / r2
               rr5 = 3.0d0 * rr3 / r2
               rr7 = 5.0d0 * rr5 / r2
               rr9 = 7.0d0 * rr7 / r2
               rr11 = 9.0d0 * rr9 / r2
               scale3 = 1.0d0
               scale5 = 1.0d0
               scale7 = 1.0d0
               do j = 1, 3
                  ddsc3(j) = 0.0d0
                  ddsc5(j) = 0.0d0
                  ddsc7(j) = 0.0d0
               end do
               damp = pdi * pdamp(k)
               if (damp .ne. 0.0d0) then
                  pgamma = min(pti,thole(k))
                  damp = -pgamma * (r/damp)**3
                  if (damp .gt. -50.0d0) then
                     expdamp = exp(damp)
                     scale3 = 1.0d0 - expdamp
                     scale5 = 1.0d0 - (1.0d0-damp)*expdamp
                     scale7 = 1.0d0 - (1.0d0-damp+0.6d0*damp**2)
     &                                       *expdamp
                     temp3 = -3.0d0 * damp * expdamp / r2
                     temp5 = -damp
                     temp7 = -0.2d0 - 0.6d0*damp
                     ddsc3(1) = temp3 * xr
                     ddsc3(2) = temp3 * yr
                     ddsc3(3) = temp3 * zr
                     ddsc5(1) = temp5 * ddsc3(1)
                     ddsc5(2) = temp5 * ddsc3(2)
                     ddsc5(3) = temp5 * ddsc3(3)
                     ddsc7(1) = temp7 * ddsc5(1)
                     ddsc7(2) = temp7 * ddsc5(2)
                     ddsc7(3) = temp7 * ddsc5(3)
                  end if
               end if
               scale3i = scale3 * uscale(kk)
               scale5i = scale5 * uscale(kk)
               scale7i = scale7 * uscale(kk)
               dsc3 = scale3 * dscale(kk)
               dsc5 = scale5 * dscale(kk)
               dsc7 = scale7 * dscale(kk)
               psc3 = scale3 * pscale(kk)
               psc5 = scale5 * pscale(kk)
               psc7 = scale7 * pscale(kk)
c
c     construct necessary auxiliary vectors
c

c               dixdk(1) = di(2)*dk(3) - di(3)*dk(2)
c               dixdk(2) = di(3)*dk(1) - di(1)*dk(3)
c               dixdk(3) = di(1)*dk(2) - di(2)*dk(1)
               dixuk(1) = di(2)*uind(3,l3) - di(3)*uind(2,l3)
               dixuk(2) = di(3)*uind(1,l3) - di(1)*uind(3,l3)
               dixuk(3) = di(1)*uind(2,l3) - di(2)*uind(1,l3)
               dkxui(1) = dk(2)*uind(3,l1) - dk(3)*uind(2,l1)
               dkxui(2) = dk(3)*uind(1,l1) - dk(1)*uind(3,l1)
               dkxui(3) = dk(1)*uind(2,l1) - dk(2)*uind(1,l1)
               dixukp(1) = di(2)*uinp(3,l3) - di(3)*uinp(2,l3)
               dixukp(2) = di(3)*uinp(1,l3) - di(1)*uinp(3,l3)
               dixukp(3) = di(1)*uinp(2,l3) - di(2)*uinp(1,l3)
               dkxuip(1) = dk(2)*uinp(3,l1) - dk(3)*uinp(2,l1)
               dkxuip(2) = dk(3)*uinp(1,l1) - dk(1)*uinp(3,l1)
               dkxuip(3) = dk(1)*uinp(2,l1) - dk(2)*uinp(1,l1)
               dixr(1) = di(2)*zr - di(3)*yr
               dixr(2) = di(3)*xr - di(1)*zr
               dixr(3) = di(1)*yr - di(2)*xr
               dkxr(1) = dk(2)*zr - dk(3)*yr
               dkxr(2) = dk(3)*xr - dk(1)*zr
               dkxr(3) = dk(1)*yr - dk(2)*xr
               qir(1) = qi(1)*xr + qi(4)*yr + qi(7)*zr
               qir(2) = qi(2)*xr + qi(5)*yr + qi(8)*zr
               qir(3) = qi(3)*xr + qi(6)*yr + qi(9)*zr
               qkr(1) = qk(1)*xr + qk(4)*yr + qk(7)*zr
               qkr(2) = qk(2)*xr + qk(5)*yr + qk(8)*zr
               qkr(3) = qk(3)*xr + qk(6)*yr + qk(9)*zr
c               qiqkr(1) = qi(1)*qkr(1) + qi(4)*qkr(2) + qi(7)*qkr(3)
c               qiqkr(2) = qi(2)*qkr(1) + qi(5)*qkr(2) + qi(8)*qkr(3)
c               qiqkr(3) = qi(3)*qkr(1) + qi(6)*qkr(2) + qi(9)*qkr(3)
c               qkqir(1) = qk(1)*qir(1) + qk(4)*qir(2) + qk(7)*qir(3)
c               qkqir(2) = qk(2)*qir(1) + qk(5)*qir(2) + qk(8)*qir(3)
c               qkqir(3) = qk(3)*qir(1) + qk(6)*qir(2) + qk(9)*qir(3)
c               qixqk(1) = qi(2)*qk(3) + qi(5)*qk(6) + qi(8)*qk(9)
c     &                       - qi(3)*qk(2) - qi(6)*qk(5) - qi(9)*qk(8)
c               qixqk(2) = qi(3)*qk(1) + qi(6)*qk(4) + qi(9)*qk(7)
c     &                       - qi(1)*qk(3) - qi(4)*qk(6) - qi(7)*qk(9)
c               qixqk(3) = qi(1)*qk(2) + qi(4)*qk(5) + qi(7)*qk(8)
c     &                       - qi(2)*qk(1) - qi(5)*qk(4) - qi(8)*qk(7)
               rxqir(1) = yr*qir(3) - zr*qir(2)
               rxqir(2) = zr*qir(1) - xr*qir(3)
               rxqir(3) = xr*qir(2) - yr*qir(1)
               rxqkr(1) = yr*qkr(3) - zr*qkr(2)
               rxqkr(2) = zr*qkr(1) - xr*qkr(3)
               rxqkr(3) = xr*qkr(2) - yr*qkr(1)
c               rxqikr(1) = yr*qiqkr(3) - zr*qiqkr(2)
c               rxqikr(2) = zr*qiqkr(1) - xr*qiqkr(3)
c               rxqikr(3) = xr*qiqkr(2) - yr*qiqkr(1)
c               rxqkir(1) = yr*qkqir(3) - zr*qkqir(2)
c               rxqkir(2) = zr*qkqir(1) - xr*qkqir(3)
c               rxqkir(3) = xr*qkqir(2) - yr*qkqir(1)
c               qkrxqir(1) = qkr(2)*qir(3) - qkr(3)*qir(2)
c               qkrxqir(2) = qkr(3)*qir(1) - qkr(1)*qir(3)
c               qkrxqir(3) = qkr(1)*qir(2) - qkr(2)*qir(1)
c               qidk(1) = qi(1)*dk(1) + qi(4)*dk(2) + qi(7)*dk(3)
c               qidk(2) = qi(2)*dk(1) + qi(5)*dk(2) + qi(8)*dk(3)
c               qidk(3) = qi(3)*dk(1) + qi(6)*dk(2) + qi(9)*dk(3)
c               qkdi(1) = qk(1)*di(1) + qk(4)*di(2) + qk(7)*di(3)
c               qkdi(2) = qk(2)*di(1) + qk(5)*di(2) + qk(8)*di(3)
c               qkdi(3) = qk(3)*di(1) + qk(6)*di(2) + qk(9)*di(3)
               qiuk(1) = qi(1)*uind(1,l3) + qi(4)*uind(2,l3)
     &                      + qi(7)*uind(3,l3)
               qiuk(2) = qi(2)*uind(1,l3) + qi(5)*uind(2,l3)
     &                      + qi(8)*uind(3,l3)
               qiuk(3) = qi(3)*uind(1,l3) + qi(6)*uind(2,l3)
     &                      + qi(9)*uind(3,l3)
               qkui(1) = qk(1)*uind(1,l1) + qk(4)*uind(2,l1)
     &                      + qk(7)*uind(3,l1)
               qkui(2) = qk(2)*uind(1,l1) + qk(5)*uind(2,l1)
     &                      + qk(8)*uind(3,l1)
               qkui(3) = qk(3)*uind(1,l1) + qk(6)*uind(2,l1)
     &                      + qk(9)*uind(3,l1)
               qiukp(1) = qi(1)*uinp(1,l3) + qi(4)*uinp(2,l3)
     &                       + qi(7)*uinp(3,l3)
               qiukp(2) = qi(2)*uinp(1,l3) + qi(5)*uinp(2,l3)
     &                       + qi(8)*uinp(3,l3)
               qiukp(3) = qi(3)*uinp(1,l3) + qi(6)*uinp(2,l3)
     &                       + qi(9)*uinp(3,l3)
               qkuip(1) = qk(1)*uinp(1,l1) + qk(4)*uinp(2,l1)
     &                       + qk(7)*uinp(3,l1)
               qkuip(2) = qk(2)*uinp(1,l1) + qk(5)*uinp(2,l1)
     &                       + qk(8)*uinp(3,l1)
               qkuip(3) = qk(3)*uinp(1,l1) + qk(6)*uinp(2,l1)
     &                       + qk(9)*uinp(3,l1)
c               dixqkr(1) = di(2)*qkr(3) - di(3)*qkr(2)
c               dixqkr(2) = di(3)*qkr(1) - di(1)*qkr(3)
c               dixqkr(3) = di(1)*qkr(2) - di(2)*qkr(1)
c               dkxqir(1) = dk(2)*qir(3) - dk(3)*qir(2)
c               dkxqir(2) = dk(3)*qir(1) - dk(1)*qir(3)
c               dkxqir(3) = dk(1)*qir(2) - dk(2)*qir(1)
               uixqkr(1) = uind(2,l1)*qkr(3) - uind(3,l1)*qkr(2)
               uixqkr(2) = uind(3,l1)*qkr(1) - uind(1,l1)*qkr(3)
               uixqkr(3) = uind(1,l1)*qkr(2) - uind(2,l1)*qkr(1)
               ukxqir(1) = uind(2,l3)*qir(3) - uind(3,l3)*qir(2)
               ukxqir(2) = uind(3,l3)*qir(1) - uind(1,l3)*qir(3)
               ukxqir(3) = uind(1,l3)*qir(2) - uind(2,l3)*qir(1)
               uixqkrp(1) = uinp(2,l1)*qkr(3) - uinp(3,l1)*qkr(2)
               uixqkrp(2) = uinp(3,l1)*qkr(1) - uinp(1,l1)*qkr(3)
               uixqkrp(3) = uinp(1,l1)*qkr(2) - uinp(2,l1)*qkr(1)
               ukxqirp(1) = uinp(2,l3)*qir(3) - uinp(3,l3)*qir(2)
               ukxqirp(2) = uinp(3,l3)*qir(1) - uinp(1,l3)*qir(3)
               ukxqirp(3) = uinp(1,l3)*qir(2) - uinp(2,l3)*qir(1)
c               rxqidk(1) = yr*qidk(3) - zr*qidk(2)
c               rxqidk(2) = zr*qidk(1) - xr*qidk(3)
c               rxqidk(3) = xr*qidk(2) - yr*qidk(1)
c               rxqkdi(1) = yr*qkdi(3) - zr*qkdi(2)
c               rxqkdi(2) = zr*qkdi(1) - xr*qkdi(3)
c               rxqkdi(3) = xr*qkdi(2) - yr*qkdi(1)
               rxqiuk(1) = yr*qiuk(3) - zr*qiuk(2)
               rxqiuk(2) = zr*qiuk(1) - xr*qiuk(3)
               rxqiuk(3) = xr*qiuk(2) - yr*qiuk(1)
               rxqkui(1) = yr*qkui(3) - zr*qkui(2)
               rxqkui(2) = zr*qkui(1) - xr*qkui(3)
               rxqkui(3) = xr*qkui(2) - yr*qkui(1)
               rxqiukp(1) = yr*qiukp(3) - zr*qiukp(2)
               rxqiukp(2) = zr*qiukp(1) - xr*qiukp(3)
               rxqiukp(3) = xr*qiukp(2) - yr*qiukp(1)
               rxqkuip(1) = yr*qkuip(3) - zr*qkuip(2)
               rxqkuip(2) = zr*qkuip(1) - xr*qkuip(3)
               rxqkuip(3) = xr*qkuip(2) - yr*qkuip(1)
c
c     calculate scalar products for permanent components
c
               sc(2) = di(1)*dk(1) + di(2)*dk(2) + di(3)*dk(3)
               sc(3) = di(1)*xr + di(2)*yr + di(3)*zr
               sc(4) = dk(1)*xr + dk(2)*yr + dk(3)*zr
               sc(5) = qir(1)*xr + qir(2)*yr + qir(3)*zr
               sc(6) = qkr(1)*xr + qkr(2)*yr + qkr(3)*zr
               sc(7) = qir(1)*dk(1) + qir(2)*dk(2) + qir(3)*dk(3)
               sc(8) = qkr(1)*di(1) + qkr(2)*di(2) + qkr(3)*di(3)
               sc(9) = qir(1)*qkr(1) + qir(2)*qkr(2) + qir(3)*qkr(3)
               sc(10) = qi(1)*qk(1) + qi(2)*qk(2) + qi(3)*qk(3)
     &                     + qi(4)*qk(4) + qi(5)*qk(5) + qi(6)*qk(6)
     &                     + qi(7)*qk(7) + qi(8)*qk(8) + qi(9)*qk(9)
c
c     calculate scalar products for induced components
c
               sci(1) = uind(1,l1)*dk(1) + uind(2,l1)*dk(2)
     &                     + uind(3,l1)*dk(3) + di(1)*uind(1,l3)
     &                     + di(2)*uind(2,l3) + di(3)*uind(3,l3)
               sci(2) = uind(1,l1)*uind(1,l3) + uind(2,l1)*uind(2,l3)
     &                     + uind(3,l1)*uind(3,l3)
               sci(3) = uind(1,l1)*xr + uind(2,l1)*yr + uind(3,l1)*zr
               sci(4) = uind(1,l3)*xr + uind(2,l3)*yr + uind(3,l3)*zr
               sci(7) = qir(1)*uind(1,l3) + qir(2)*uind(2,l3)
     &                     + qir(3)*uind(3,l3)
               sci(8) = qkr(1)*uind(1,l1) + qkr(2)*uind(2,l1)
     &                     + qkr(3)*uind(3,l1)
               scip(1) = uinp(1,l1)*dk(1) + uinp(2,l1)*dk(2)
     &                      + uinp(3,l1)*dk(3) + di(1)*uinp(1,l3)
     &                      + di(2)*uinp(2,l3) + di(3)*uinp(3,l3)
               scip(2) = uind(1,l1)*uinp(1,l3)+uind(2,l1)*uinp(2,l3)
     &                   + uind(3,l1)*uinp(3,l3)+uinp(1,l1)*uind(1,l3)
     &                   + uinp(2,l1)*uind(2,l3)+uinp(3,l1)*uind(3,l3)
               scip(3) = uinp(1,l1)*xr + uinp(2,l1)*yr + uinp(3,l1)*zr
               scip(4) = uinp(1,l3)*xr + uinp(2,l3)*yr + uinp(3,l3)*zr
               scip(7) = qir(1)*uinp(1,l3) + qir(2)*uinp(2,l3)
     &                      + qir(3)*uinp(3,l3)
               scip(8) = qkr(1)*uinp(1,l1) + qkr(2)*uinp(2,l1)
     &                      + qkr(3)*uinp(3,l1)
c
c     calculate the gl functions for permanent components
c

c               gl(0) = ci*ck
c               gl(1) = ck*sc(3) - ci*sc(4)
c               gl(2) = ci*sc(6) + ck*sc(5) - sc(3)*sc(4)
c               gl(3) = sc(3)*sc(6) - sc(4)*sc(5)
c               gl(4) = sc(5)*sc(6)
c               gl(5) = -4.0d0 * sc(9)
c               gl(6) = sc(2)
c               gl(7) = 2.0d0 * (sc(7)-sc(8))
c               gl(8) = 2.0d0 * sc(10)

c
c     calculate the gl functions for induced components
c
               gli(1) = ck*sci(3) - ci*sci(4)
               gli(2) = -sc(3)*sci(4) - sci(3)*sc(4)
               gli(3) = sci(3)*sc(6) - sci(4)*sc(5)
               gli(6) = sci(1)
               gli(7) = 2.0d0 * (sci(7)-sci(8))
               glip(1) = ck*scip(3) - ci*scip(4)
               glip(2) = -sc(3)*scip(4) - scip(3)*sc(4)
               glip(3) = scip(3)*sc(6) - scip(4)*sc(5)
               glip(6) = scip(1)
               glip(7) = 2.0d0 * (scip(7)-scip(8))
c
c     compute the energy contributions for this interaction
c
               ei = 0.5d0*(rr3*(gli(1)+gli(6))*psc3
     &                   + rr5*(gli(2)+gli(7))*psc5
     &                   + rr7*gli(3)*psc7)
               ei = f * ei
               eptemp = eptemp + ei
c             print*,"eptemp in empole3a=",eptemp
c
c     increment the total intermolecular energy
c

c               if (molcule(ii) .ne. molcule(kk)) then
c                  einter = einter + e + ei
c               end if

c
c     intermediate variables for the permanent components
c

c               gf(1) = rr3*gl(0) + rr5*(gl(1)+gl(6))
c     &                    + rr7*(gl(2)+gl(7)+gl(8))
c     &                    + rr9*(gl(3)+gl(5)) + rr11*gl(4)
c               gf(2) = -ck*rr3 + sc(4)*rr5 - sc(6)*rr7
c               gf(3) = ci*rr3 + sc(3)*rr5 + sc(5)*rr7
c               gf(4) = 2.0d0 * rr5
c               gf(5) = 2.0d0 * (-ck*rr5+sc(4)*rr7-sc(6)*rr9)
c               gf(6) = 2.0d0 * (-ci*rr5-sc(3)*rr7-sc(5)*rr9)
c               gf(7) = 4.0d0 * rr7

c
c     intermediate variables for the induced components
c
               gfi(1) = 0.5d0 * rr5 * ((gli(1)+gli(6))*psc3
     &                                   + (glip(1)+glip(6))*dsc3
     &                                   + scip(2)*scale3i)
     &                + 0.5d0 * rr7 * ((gli(7)+gli(2))*psc5
     &                               + (glip(7)+glip(2))*dsc5
     &                      - (sci(3)*scip(4)+scip(3)*sci(4))*scale5i)
     &                + 0.5d0 * rr9 * (gli(3)*psc7+glip(3)*dsc7)
               gfi(2) = -rr3*ck + rr5*sc(4) - rr7*sc(6)
               gfi(3) = rr3*ci + rr5*sc(3) + rr7*sc(5)
               gfi(4) = 2.0d0 * rr5
               gfi(5) = rr7 * (sci(4)*psc7+scip(4)*dsc7)
               gfi(6) = -rr7 * (sci(3)*psc7+scip(3)*dsc7)
c
c     get the permanent force components
c

c               ftm2(1) = gf(1)*xr + gf(2)*di(1) + gf(3)*dk(1)
c     &                      + gf(4)*(qkdi(1)-qidk(1)) + gf(5)*qir(1)
c     &                      + gf(6)*qkr(1) + gf(7)*(qiqkr(1)+qkqir(1))
c               ftm2(2) = gf(1)*yr + gf(2)*di(2) + gf(3)*dk(2)
c     &                      + gf(4)*(qkdi(2)-qidk(2)) + gf(5)*qir(2)
c     &                      + gf(6)*qkr(2) + gf(7)*(qiqkr(2)+qkqir(2))
c               ftm2(3) = gf(1)*zr + gf(2)*di(3) + gf(3)*dk(3)
c     &                      + gf(4)*(qkdi(3)-qidk(3)) + gf(5)*qir(3)
c     &                      + gf(6)*qkr(3) + gf(7)*(qiqkr(3)+qkqir(3))

c
c     get the induced force components
c
               ftm2i(1) = gfi(1)*xr + 0.5d0*
     &           (- rr3*ck*(uind(1,l1)*psc3+uinp(1,l1)*dsc3)
     &            + rr5*sc(4)*(uind(1,l1)*psc5+uinp(1,l1)*dsc5)
     &            - rr7*sc(6)*(uind(1,l1)*psc7+uinp(1,l1)*dsc7))
     &            + (rr3*ci*(uind(1,l3)*psc3+uinp(1,l3)*dsc3)
     &            + rr5*sc(3)*(uind(1,l3)*psc5+uinp(1,l3)*dsc5)
     &            + rr7*sc(5)*(uind(1,l3)*psc7+uinp(1,l3)*dsc7))*0.5d0
     &            + rr5*scale5i*(sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
     &            + sci(3)*uinp(1,l3)+scip(3)*uind(1,l3))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(1)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(1)
     &            + 0.5d0*gfi(4)*((qkui(1)-qiuk(1))*psc5
     &            + (qkuip(1)-qiukp(1))*dsc5)
     &            + gfi(5)*qir(1) + gfi(6)*qkr(1)
               ftm2i(2) = gfi(1)*yr + 0.5d0*
     &           (- rr3*ck*(uind(2,l1)*psc3+uinp(2,l1)*dsc3)
     &            + rr5*sc(4)*(uind(2,l1)*psc5+uinp(2,l1)*dsc5)
     &            - rr7*sc(6)*(uind(2,l1)*psc7+uinp(2,l1)*dsc7))
     &            + (rr3*ci*(uind(2,l3)*psc3+uinp(2,l3)*dsc3)
     &            + rr5*sc(3)*(uind(2,l3)*psc5+uinp(2,l3)*dsc5)
     &            + rr7*sc(5)*(uind(2,l3)*psc7+uinp(2,l3)*dsc7))*0.5d0
     &            + rr5*scale5i*(sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
     &            + sci(3)*uinp(2,l3)+scip(3)*uind(2,l3))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(2)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(2)
     &            + 0.5d0*gfi(4)*((qkui(2)-qiuk(2))*psc5
     &            + (qkuip(2)-qiukp(2))*dsc5)
     &            + gfi(5)*qir(2) + gfi(6)*qkr(2)
               ftm2i(3) = gfi(1)*zr  + 0.5d0*
     &           (- rr3*ck*(uind(3,l1)*psc3+uinp(3,l1)*dsc3)
     &            + rr5*sc(4)*(uind(3,l1)*psc5+uinp(3,l1)*dsc5)
     &            - rr7*sc(6)*(uind(3,l1)*psc7+uinp(3,l1)*dsc7))
     &            + (rr3*ci*(uind(3,l3)*psc3+uinp(3,l3)*dsc3)
     &            + rr5*sc(3)*(uind(3,l3)*psc5+uinp(3,l3)*dsc5)
     &            + rr7*sc(5)*(uind(3,l3)*psc7+uinp(3,l3)*dsc7))*0.5d0
     &            + rr5*scale5i*(sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
     &            + sci(3)*uinp(3,l3)+scip(3)*uind(3,l3))*0.5d0
     &            + 0.5d0*(sci(4)*psc5+scip(4)*dsc5)*rr5*di(3)
     &            + 0.5d0*(sci(3)*psc5+scip(3)*dsc5)*rr5*dk(3)
     &            + 0.5d0*gfi(4)*((qkui(3)-qiuk(3))*psc5
     &            + (qkuip(3)-qiukp(3))*dsc5)
     &            + gfi(5)*qir(3) + gfi(6)*qkr(3)
c
c     account for partially excluded induced interactions
c
               temp3 = 0.5d0 * rr3 * ((gli(1)+gli(6))*pscale(kk)
     &                                  +(glip(1)+glip(6))*dscale(kk))
               temp5 = 0.5d0 * rr5 * ((gli(2)+gli(7))*pscale(kk)
     &                                  +(glip(2)+glip(7))*dscale(kk))
               temp7 = 0.5d0 * rr7 * (gli(3)*pscale(kk)
     &                                  +glip(3)*dscale(kk))
               fridmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
     &                        + temp7*ddsc7(1)
               fridmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
     &                        + temp7*ddsc7(2)
               fridmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
     &                        + temp7*ddsc7(3)
c
c     find some scaling terms for induced-induced force
c
               temp3 = 0.5d0 * rr3 * uscale(kk) * scip(2)
               temp5 = -0.5d0 * rr5 * uscale(kk)
     &                    * (sci(3)*scip(4)+scip(3)*sci(4))
               findmp(1) = temp3*ddsc3(1) + temp5*ddsc5(1)
               findmp(2) = temp3*ddsc3(2) + temp5*ddsc5(2)
               findmp(3) = temp3*ddsc3(3) + temp5*ddsc5(3)
c
c     modify induced force for partially excluded interactions
c
               ftm2i(1) = ftm2i(1) - fridmp(1) - findmp(1)
               ftm2i(2) = ftm2i(2) - fridmp(2) - findmp(2)
               ftm2i(3) = ftm2i(3) - fridmp(3) - findmp(3)
c
c     correction to convert mutual to direct polarization force
c

c               if (poltyp .eq. 'DIRECT') then
c                  gfd = 0.5d0 * (rr5*scip(2)*scale3i
c     &                  - rr7*(scip(3)*sci(4)+sci(3)*scip(4))*scale5i)
c                  temp5 = 0.5d0 * rr5 * scale5i
c                  fdir(1) = gfd*xr + temp5
c     &                         * (sci(4)*uinp(1,l1)+scip(4)*uind(1,l1)
c     &                           +sci(3)*uinp(1,l3)+scip(3)*uind(1,l3))
c                  fdir(2) = gfd*yr + temp5
c     &                         * (sci(4)*uinp(2,l1)+scip(4)*uind(2,l1)
c     &                           +sci(3)*uinp(2,l3)+scip(3)*uind(2,l3))
c                  fdir(3) = gfd*zr + temp5
c     &                         * (sci(4)*uinp(3,l1)+scip(4)*uind(3,l1)
c     &                           +sci(3)*uinp(3,l3)+scip(3)*uind(3,l3))
c                  ftm2i(1) = ftm2i(1) - fdir(1) + findmp(1)
c                  ftm2i(2) = ftm2i(2) - fdir(2) + findmp(2)
c                  ftm2i(3) = ftm2i(3) - fdir(3) + findmp(3)
c               end if


c
c     intermediate terms for induced torque on multipoles
c
               gti(2) = 0.5d0 * rr5 * (sci(4)*psc5+scip(4)*dsc5)
               gti(3) = 0.5d0 * rr5 * (sci(3)*psc5+scip(3)*dsc5)
               gti(4) = gfi(4)
               gti(5) = gfi(5)
               gti(6) = gfi(6)
c
c     get the permanent torque components
c

c               ttm2(1) = -rr3*dixdk(1) + gf(2)*dixr(1) - gf(5)*rxqir(1)
c     &           + gf(4)*(dixqkr(1)+dkxqir(1)+rxqidk(1)-2.0d0*qixqk(1))
c     &           - gf(7)*(rxqikr(1)+qkrxqir(1))
c               ttm2(2) = -rr3*dixdk(2) + gf(2)*dixr(2) - gf(5)*rxqir(2)
c     &           + gf(4)*(dixqkr(2)+dkxqir(2)+rxqidk(2)-2.0d0*qixqk(2))
c     &           - gf(7)*(rxqikr(2)+qkrxqir(2))
c               ttm2(3) = -rr3*dixdk(3) + gf(2)*dixr(3) - gf(5)*rxqir(3)
c     &           + gf(4)*(dixqkr(3)+dkxqir(3)+rxqidk(3)-2.0d0*qixqk(3))
c     &           - gf(7)*(rxqikr(3)+qkrxqir(3))
c               ttm3(1) = rr3*dixdk(1) + gf(3)*dkxr(1) - gf(6)*rxqkr(1)
c     &           - gf(4)*(dixqkr(1)+dkxqir(1)+rxqkdi(1)-2.0d0*qixqk(1))
c     &           - gf(7)*(rxqkir(1)-qkrxqir(1))
c               ttm3(2) = rr3*dixdk(2) + gf(3)*dkxr(2) - gf(6)*rxqkr(2)
c     &           - gf(4)*(dixqkr(2)+dkxqir(2)+rxqkdi(2)-2.0d0*qixqk(2))
c     &           - gf(7)*(rxqkir(2)-qkrxqir(2))
c               ttm3(3) = rr3*dixdk(3) + gf(3)*dkxr(3) - gf(6)*rxqkr(3)
c     &           - gf(4)*(dixqkr(3)+dkxqir(3)+rxqkdi(3)-2.0d0*qixqk(3))
c     &           - gf(7)*(rxqkir(3)-qkrxqir(3))

c
c     get the induced torque components
c
               ttm2i(1) = -rr3*(dixuk(1)*psc3+dixukp(1)*dsc3)*0.5d0
     &           + gti(2)*dixr(1) + gti(4)*((ukxqir(1)+rxqiuk(1))*psc5
     &           +(ukxqirp(1)+rxqiukp(1))*dsc5)*0.5d0 - gti(5)*rxqir(1)
               ttm2i(2) = -rr3*(dixuk(2)*psc3+dixukp(2)*dsc3)*0.5d0
     &           + gti(2)*dixr(2) + gti(4)*((ukxqir(2)+rxqiuk(2))*psc5
     &           +(ukxqirp(2)+rxqiukp(2))*dsc5)*0.5d0 - gti(5)*rxqir(2)
               ttm2i(3) = -rr3*(dixuk(3)*psc3+dixukp(3)*dsc3)*0.5d0
     &           + gti(2)*dixr(3) + gti(4)*((ukxqir(3)+rxqiuk(3))*psc5
     &           +(ukxqirp(3)+rxqiukp(3))*dsc5)*0.5d0 - gti(5)*rxqir(3)
               ttm3i(1) = -rr3*(dkxui(1)*psc3+dkxuip(1)*dsc3)*0.5d0
     &           + gti(3)*dkxr(1) - gti(4)*((uixqkr(1)+rxqkui(1))*psc5
     &           +(uixqkrp(1)+rxqkuip(1))*dsc5)*0.5d0 - gti(6)*rxqkr(1)
               ttm3i(2) = -rr3*(dkxui(2)*psc3+dkxuip(2)*dsc3)*0.5d0
     &           + gti(3)*dkxr(2) - gti(4)*((uixqkr(2)+rxqkui(2))*psc5
     &           +(uixqkrp(2)+rxqkuip(2))*dsc5)*0.5d0 - gti(6)*rxqkr(2)
               ttm3i(3) = -rr3*(dkxui(3)*psc3+dkxuip(3)*dsc3)*0.5d0
     &           + gti(3)*dkxr(3) - gti(4)*((uixqkr(3)+rxqkui(3))*psc5
     &           +(uixqkrp(3)+rxqkuip(3))*dsc5)*0.5d0 - gti(6)*rxqkr(3)
c
c     handle the case where scaling is used
c
               do j = 1, 3
c                  ftm2(j) = f * ftm2(j) * mscale(kk)
                  ftm2i(j) = f * ftm2i(j)
c                  ttm2(j) = f * ttm2(j) * mscale(kk)
                  ttm2i(j) = f * ttm2i(j)
c                  ttm3(j) = f * ttm3(j) * mscale(kk)
                  ttm3i(j) = f * ttm3i(j)
               end do
c
c     increment gradient due to force and torque on first site
c

c               dem(1,ii) = dem(1,ii) + ftm2(1)
c               dem(2,ii) = dem(2,ii) + ftm2(2)
c               dem(3,ii) = dem(3,ii) + ftm2(3)
c               dep(1,ii) = dep(1,ii) + ftm2i(1)
c               dep(2,ii) = dep(2,ii) + ftm2i(2)
c               dep(3,ii) = dep(3,ii) + ftm2i(3)
c               call torque (i,ttm2,ttm2i,frcxi,frcyi,frczi)
               deptemp(1,i) = deptemp(1,i) + ftm2i(1)
               deptemp(2,i) = deptemp(2,i) + ftm2i(2)
               deptemp(3,i) = deptemp(3,i) + ftm2i(3)
               call torque_3b (deptemp,i,
     &          ttm2,ttm2i,frcxi,frcyi,frczi)

c
c     increment gradient due to force and torque on second site
c

c               dem(1,kk) = dem(1,kk) - ftm2(1)
c               dem(2,kk) = dem(2,kk) - ftm2(2)
c               dem(3,kk) = dem(3,kk) - ftm2(3)
c               dep(1,kk) = dep(1,kk) - ftm2i(1)
c               dep(2,kk) = dep(2,kk) - ftm2i(2)
c               dep(3,kk) = dep(3,kk) - ftm2i(3)
c               call torque (k,ttm3,ttm3i,frcxk,frcyk,frczk)
               deptemp(1,k) = deptemp(1,k) - ftm2i(1)
               deptemp(2,k) = deptemp(2,k) - ftm2i(2)
               deptemp(3,k) = deptemp(3,k) - ftm2i(3)
               call torque_3b (deptemp,k,
     &            ttm3,ttm3i,frcxk,frcyk,frczk)

c
c     increment the internal virial tensor components
c
               iaz = iz
               iax = iz
               iay = iy
               kaz = kz
               kax = kx
               kay = ky
               if (iaz .eq. 0)  iaz = ii
               if (iax .eq. 0)  iax = ii
               if (iay .eq. 0)  iay = ii
               if (kaz .eq. 0)  kaz = kk
               if (kax .eq. 0)  kax = kk
               if (kay .eq. 0)  kay = kk
               xiz = x(iaz) - x(ii)
               yiz = y(iaz) - y(ii)
               ziz = z(iaz) - z(ii)
               xix = x(iax) - x(ii)
               yix = y(iax) - y(ii)
               zix = z(iax) - z(ii)
               xiy = x(iay) - x(ii)
               yiy = y(iay) - y(ii)
               ziy = z(iay) - z(ii)
               xkz = x(kaz) - x(kk)
               ykz = y(kaz) - y(kk)
               zkz = z(kaz) - z(kk)
               xkx = x(kax) - x(kk)
               ykx = y(kax) - y(kk)
               zkx = z(kax) - z(kk)
               xky = x(kay) - x(kk)
               yky = y(kay) - y(kk)
               zky = z(kay) - z(kk)
c               vxx = -xr*(ftm2(1)+ftm2i(1)) + xix*frcxi(1)
c     &                  + xiy*frcyi(1) + xiz*frczi(1) + xkx*frcxk(1)
c     &                  + xky*frcyk(1) + xkz*frczk(1)
c               vyx = -yr*(ftm2(1)+ftm2i(1)) + yix*frcxi(1)
c     &                  + yiy*frcyi(1) + yiz*frczi(1) + ykx*frcxk(1)
c     &                  + yky*frcyk(1) + ykz*frczk(1)
c               vzx = -zr*(ftm2(1)+ftm2i(1)) + zix*frcxi(1)
c     &                  + ziy*frcyi(1) + ziz*frczi(1) + zkx*frcxk(1)
c     &                  + zky*frcyk(1) + zkz*frczk(1)
c               vyy = -yr*(ftm2(2)+ftm2i(2)) + yix*frcxi(2)
c     &                  + yiy*frcyi(2) + yiz*frczi(2) + ykx*frcxk(2)
c     &                  + yky*frcyk(2) + ykz*frczk(2)
c               vzy = -zr*(ftm2(2)+ftm2i(2)) + zix*frcxi(2)
c     &                  + ziy*frcyi(2) + ziz*frczi(2) + zkx*frcxk(2)
c     &                  + zky*frcyk(2) + zkz*frczk(2)
c               vzz = -zr*(ftm2(3)+ftm2i(3)) + zix*frcxi(3)
c     &                  + ziy*frcyi(3) + ziz*frczi(3) + zkx*frcxk(3)
c     &                  + zky*frcyk(3) + zkz*frczk(3)
c               vir(1,1) = vir(1,1) + vxx
c               vir(2,1) = vir(2,1) + vyx
c               vir(3,1) = vir(3,1) + vzx
c               vir(1,2) = vir(1,2) + vyx
c               vir(2,2) = vir(2,2) + vyy
c               vir(3,2) = vir(3,2) + vzy
c               vir(1,3) = vir(1,3) + vzx
c               vir(2,3) = vir(2,3) + vzy
c               vir(3,3) = vir(3,3) + vzz
c            end if
               vxx = -xr*(ftm2i(1)) + xix*frcxi(1)
     &                  + xiy*frcyi(1) + xiz*frczi(1) + xkx*frcxk(1)
     &                  + xky*frcyk(1) + xkz*frczk(1)
               vyx = -yr*(ftm2i(1)) + yix*frcxi(1)
     &                  + yiy*frcyi(1) + yiz*frczi(1) + ykx*frcxk(1)
     &                  + yky*frcyk(1) + ykz*frczk(1)
               vzx = -zr*(ftm2i(1)) + zix*frcxi(1)
     &                  + ziy*frcyi(1) + ziz*frczi(1) + zkx*frcxk(1)
     &                  + zky*frcyk(1) + zkz*frczk(1)
               vyy = -yr*(ftm2i(2)) + yix*frcxi(2)
     &                  + yiy*frcyi(2) + yiz*frczi(2) + ykx*frcxk(2)
     &                  + yky*frcyk(2) + ykz*frczk(2)
               vzy = -zr*(ftm2i(2)) + zix*frcxi(2)
     &                  + ziy*frcyi(2) + ziz*frczi(2) + zkx*frcxk(2)
     &                  + zky*frcyk(2) + zkz*frczk(2)
               vzz = -zr*(ftm2i(3)) + zix*frcxi(3)
     &                  + ziy*frcyi(3) + ziz*frczi(3) + zkx*frcxk(3)
     &                  + zky*frcyk(3) + zkz*frczk(3)

               virtemp(1,1) = virtemp(1,1) + vxx
               virtemp(2,1) = virtemp(2,1) + vyx
               virtemp(3,1) = virtemp(3,1) + vzx
               virtemp(1,2) = virtemp(1,2) + vyx
               virtemp(2,2) = virtemp(2,2) + vyy
               virtemp(3,2) = virtemp(3,2) + vzy
               virtemp(1,3) = virtemp(1,3) + vzx
               virtemp(2,3) = virtemp(2,3) + vzy
               virtemp(3,3) = virtemp(3,3) + vzz

   10       continue
         end do
c
c     reset interaction scaling coefficients for connected atoms
c
         do j = 1, n12(ii)
            mscale(i12(j,ii)) = 1.0d0
            pscale(i12(j,ii)) = 1.0d0
         end do
         do j = 1, n13(ii)
            mscale(i13(j,ii)) = 1.0d0
            pscale(i13(j,ii)) = 1.0d0
         end do
         do j = 1, n14(ii)
            mscale(i14(j,ii)) = 1.0d0
            pscale(i14(j,ii)) = 1.0d0
         end do
         do j = 1, n15(ii)
            mscale(i15(j,ii)) = 1.0d0
            pscale(i15(j,ii)) = 1.0d0
         end do
         do j = 1, np11(ii)
            dscale(ip11(j,ii)) = 1.0d0
            uscale(ip11(j,ii)) = 1.0d0
         end do
         do j = 1, np12(ii)
            dscale(ip12(j,ii)) = 1.0d0
            uscale(ip12(j,ii)) = 1.0d0
         end do
         do j = 1, np13(ii)
            dscale(ip13(j,ii)) = 1.0d0
            uscale(ip13(j,ii)) = 1.0d0
         end do
         do j = 1, np14(ii)
            dscale(ip14(j,ii)) = 1.0d0
            uscale(ip14(j,ii)) = 1.0d0
         end do
      end do
c
c     perform deallocation of some local arrays
c

c
c     for periodic boundary conditions with large cutoffs
c     neighbors must be found by the replicates method
c
c       use_replica=.true.
c   NOTE TO SELF:  REMOVED if (use_replica) BECAUSE NCELL WAS
c   ALWAYS EQUAL TO ZERO AFTER PERIODIC BOUNDARY CONDITIONS WERE
c   IMPLEMENTED CORRECTLY.  ALSO, AS AN ADDITIONAL CHECK, THE ENERGY
c   IN THIS SECTION WAS EQUAL TO ZERO.
c
c     perform deallocation of some local arrays
c

      deallocate (mscale)
      deallocate (pscale)
      deallocate (dscale)
      deallocate (uscale)
      return
      end


      subroutine induce0a_3b_PolelecOnly(npole3b,pnum,uind,uinp)
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'boxes.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
      include 'mpole.i'
      include 'polar2.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'potent.i'
      include 'shunt.i'
      include 'units.i'
      include 'uprior.i'
      integer i,j,k,m
      integer ii,kk
c      real*8 xr,yr,zr
c      real*8 fgrp,r,r2
c      real*8 rr3,rr5,rr7
c      real*8 ci,dix,diy,diz
c      real*8 duix,duiy,duiz
c      real*8 puix,puiy,puiz
c      real*8 qixx,qixy,qixz
c      real*8 qiyy,qiyz,qizz
c      real*8 ck,dkx,dky,dkz
c      real*8 dukx,duky,dukz
c      real*8 pukx,puky,pukz
c      real*8 qkxx,qkxy,qkxz
c      real*8 qkyy,qkyz,qkzz
c      real*8 dir,duir,puir
c      real*8 dkr,dukr,pukr
c      real*8 qix,qiy,qiz,qir
c      real*8 qkx,qky,qkz,qkr
c      real*8 udsum,upsum
c      real*8 damp,expdamp
c      real*8 scale3,scale5
c      real*8 scale7
c      real*8 pdi,pti,pgamma
c      real*8 fid(3),fkd(3)
c      real*8 fip(3),fkp(3)
c      real*8, allocatable :: dscale(:)
c      real*8, allocatable :: pscale(:)
c      real*8, allocatable :: field(:,:)
c      real*8, allocatable :: fieldp(:,:)
      real*8 field(3,npole3b),off3b
      real*8 fieldp(3,npole3b)
c      real*8, allocatable :: udir(:,:)
c      real*8, allocatable :: udirp(:,:)
c      real*8, allocatable :: uold(:,:)
c      real*8, allocatable :: uoldp(:,:)
      logical proceed,done,use_liam
      character*6 mode
      integer l1,l2,l3,k1,k2,i1,i2
c      real*8 fdir(3,npole)
c      real*8 liam(3,maxatm)
      real*8 M_tot(3*npole3b,3*npole3b)
      real*8 uind(3,npole3b)
      real*8 uinp(3,npole3b)
      integer npole3b,pnum(*)


      if (.not. use_polar)  return
        do l1=1,npole3b
           do j=1,3
              field(j,l1)=0.0d0
              fieldp(j,l1)=0.0d0
           end do
        end do

         do l1 = 1, 3*npole3b
            do l3 = 1, 3*npole3b
               M_tot(l1,l3) = 0.0d0
            end do
         end do

       call field_noewald_umutual_rl_3b (field,fieldp,M_tot,
     &   npole3b,pnum)

         do l1 = 1, npole3b
c            i = pnum(l1)
            do j = 1, 3
               uind(j,l1) = 0.0d0
               uinp(j,l1) = 0.0d0
            end do
         end do

c         call umutual_rl (M_tot)

         do l1 = 1, npole3b
            i = pnum(l1)
            do j = 1, 3
               i1 = 3*(l1-1)+j
c               M_tot(i1,i1) = 1.0d0/polarity(i)
               M_tot(i1,i1) = M_tot(i1,i1)+1.0d0/polarity(i)
            end do
         end do
c         do l1 = 1, 3*npole3b
c            do l3 = 1, 3*npole3b
c               print*,"M_tot",M_tot(l1,l3)
c            end do
c         end do
c         do l1=1,npole3b
c           do j=1,3
c              print*,"Direct field l1 j",l1,j,field(j,l1)
c           end do
c        end do

         call invert(3*npole3b,M_tot)

         do l1 = 1, npole3b
            i = pnum(l1)
            do i1 = 1, 3
               i2 = 3*(l1-1) + i1
               do l3 = 1, npole3b
                  k = pnum(l3)
                  k2 = 3*(l3-1)
c                  uind(i1,i) = uind(i1,i) + M_tot(i2,k2+1)*fdir(1,k) +
c     &                                      M_tot(i2,k2+2)*fdir(2,k) +
c     &                                      M_tot(i2,k2+3)*fdir(3,k)
                  uind(i1,l1)=uind(i1,l1)+ M_tot(i2,k2+1)*field(1,l3) +
     &                                      M_tot(i2,k2+2)*field(2,l3) +
     &                                      M_tot(i2,k2+3)*field(3,l3)

               end do
c               print*,i,i1,uind(i1,i)
            end do
         end do
      do l1 = 1, npole3b
c         i = pnum(l1)
         do j = 1, 3
            uinp(j,l1)=uind(j,l1)
         end do
      end do

c      if(npole3b .eq. 9) then
c         do l1 = 1, npole3b
c            i = pnum(l1)
c            do j = 1, 3
c      print*,"moli1 moli2 moli3 uindx=",uind(1,i)
c      print*,"moli1 moli2 moli3 uindy=",uind(2,i)
c      print*,"moli1 moli2 moli3 uindz=",uind(3,i)
c            end do
c         end do
c      end if

      return
      end
c
c  field_noewald_umutual_rl_3b (field,fieldp,M_tot,npole3b,pnum)
c
      subroutine field_noewald_umutual_rl_3b(field,fieldp,M_tot,
     &   npole3b,pnum)
      implicit none
      include 'boxes.i'
      include 'sizes.i'
c      include 'ewald.i'
      include 'math.i'
      include 'atoms.i'
      include 'mpole.i'
c      include 'combo.i'
      include 'bound.i'
      include 'cell.i'
      include 'couple.i'
      include 'group.i'
      include 'inform.i'
      include 'iounit.i'
c      include 'polar.i'
      include 'polar2.i'
      include 'polgrp.i'
      include 'polpot.i'
      include 'potent.i'
      include 'shunt.i'
      include 'units.i'
      include 'uprior.i'
      include 'molcul.i'
c      real*8 term
c      real*8 ucell(3)
      real*8 field(3,*),off3b
      real*8 fieldp(3,*)
      real*8, allocatable :: dscale_dir(:)
      real*8, allocatable :: dscale_mut(:)
      real*8, allocatable :: pscale(:)
      real*8 xr,yr,zr
      real*8 fgrp,r,r2
      real*8 rr3_dir,rr5_dir,rr7_dir
      real*8 rr3,rr5
      real*8 ci,dix,diy,diz
      real*8 duix,duiy,duiz
      real*8 puix,puiy,puiz
      real*8 qixx,qixy,qixz
      real*8 qiyy,qiyz,qizz
      real*8 ck,dkx,dky,dkz
      real*8 dukx,duky,dukz
      real*8 pukx,puky,pukz
      real*8 qkxx,qkxy,qkxz
      real*8 qkyy,qkyz,qkzz
      real*8 dir,duir,puir
      real*8 dkr,dukr,pukr
      real*8 qix,qiy,qiz,qir
      real*8 qkx,qky,qkz,qkr
      real*8 damp,expdamp
      real*8 scale3_dir,scale5_dir
      real*8 scale7_dir
      real*8 scale3_mut,scale5_mut
      real*8 pdi,pti,pgamma
      real*8 fid(3),fkd(3),M_tot(3*npole3b,3*npole3b)
      real*8 fip(3),fkp(3)
      logical proceed
      character*6 mode
      integer i,ii,j,l1,k,m,kk
      integer l2,l3,k1,k2,i1,i2
      integer pnum(*),npole3b
      real*8 Txx,Txy,Txz,Tyx
      real*8 Tyy,Tyz,Tzx,Tzy,Tzz
      

      allocate (dscale_dir(n))
      allocate (dscale_mut(n))
      allocate (pscale(n))
c
c     set the switching function coefficients
c

c      mode = 'MPOLE'
c      call switch (mode)

c      mode = 'MPOLE3'
c      call switch (mode,off3b)
c       call replica (off3b)
c       use_replica = .false.

c
c     compute the direct induced dipole moment at each atom
c
c      if(use_replica.eq..false.) then
      do l1 = 1, npole3b-1
         i = pnum(l1)
         i2 = 3*(l1-1)
         ii = ipole(i)
         pdi = pdamp(i)
         pti = thole(i)
         ci = rpole(1,i)
         dix = rpole(2,i)
         diy = rpole(3,i)
         diz = rpole(4,i)
         qixx = rpole(5,i)
         qixy = rpole(6,i)
         qixz = rpole(7,i)
         qiyy = rpole(9,i)
         qiyz = rpole(10,i)
         qizz = rpole(13,i)
         do j = i+1, npole
            dscale_dir(ipole(j)) = 1.0d0
            pscale(ipole(j)) = 1.0d0
            dscale_mut(ipole(j)) = 1.0d0
         end do
         do j = 1, n12(ii)
            pscale(i12(j,ii)) = p2scale
         end do
         do j = 1, n13(ii)
            pscale(i13(j,ii)) = p3scale
         end do
         do j = 1, n14(ii)
            pscale(i14(j,ii)) = p4scale
            do k = 1, np11(ii)
               if (i14(j,ii) .eq. ip11(k,ii))
     &            pscale(i14(j,ii)) = p4scale * p41scale
            end do
         end do
         do j = 1, n15(ii)
            pscale(i15(j,ii)) = p5scale
         end do
         do j = 1, np11(ii)
            dscale_dir(ip11(j,ii)) = d1scale
            dscale_mut(ip11(j,ii)) = u1scale

         end do
         do j = 1, np12(ii)
            dscale_dir(ip12(j,ii)) = d2scale
            dscale_mut(ip12(j,ii)) = u2scale

         end do
         do j = 1, np13(ii)
            dscale_dir(ip13(j,ii)) = d3scale
            dscale_mut(ip13(j,ii)) = u3scale

         end do
         do j = 1, np14(ii)
            dscale_dir(ip14(j,ii)) = d4scale
            dscale_mut(ip14(j,ii)) = u4scale
         end do
         do l3 = l1+1, npole3b
            k = pnum(l3)
            k2 = 3*(l3-1)
            kk = ipole(k)
            proceed = .true.
            if (use_intra)  call groups (proceed,fgrp,ii,kk,0,0,0,0)
            if (proceed) then
               xr = x(kk) - x(ii)
               yr = y(kk) - y(ii)
               zr = z(kk) - z(ii)
               call image (xr,yr,zr)
               r2 = xr*xr + yr* yr + zr*zr
c               if (r2 .le. off2) then
                  r = sqrt(r2)
                  ck = rpole(1,k)
                  dkx = rpole(2,k)
                  dky = rpole(3,k)
                  dkz = rpole(4,k)
                  qkxx = rpole(5,k)
                  qkxy = rpole(6,k)
                  qkxz = rpole(7,k)
                  qkyy = rpole(9,k)
                  qkyz = rpole(10,k)
                  qkzz = rpole(13,k)
                  scale3_dir = 1.0d0
                  scale5_dir = 1.0d0
                  scale7_dir = 1.0d0

                  scale3_mut = dscale_mut(kk)
                  scale5_mut = dscale_mut(kk)

                  damp = pdi * pdamp(k)
                  if (damp .ne. 0.0d0) then
                     pgamma = min(pti,thole(k))
                     damp = -pgamma * (r/damp)**3
                     if (damp .gt. -50.0d0) then
                        expdamp = exp(damp)
                        scale3_dir = 1.0d0 - expdamp
                        scale5_dir = 1.0d0 - expdamp*(1.0d0-damp)
                        scale7_dir = 1.0d0 - expdamp
     &                              *(1.0d0-damp+0.6d0*damp**2)
                        scale3_mut = scale3_mut * (1.0d0-expdamp)
                        scale5_mut = scale5_mut * (1.0d0-expdamp
     &                                        *(1.0d0-damp))

                     end if
                  end if
                  rr3_dir = scale3_dir / (r*r2)
                  rr5_dir = 3.0d0 * scale5_dir / (r*r2*r2)
                  rr7_dir = 15.0d0 * scale7_dir / (r*r2*r2*r2)

                  rr3 = scale3_mut / (r*r2)
                  rr5 = 3.0d0 * scale5_mut / (r*r2*r2)


                  Txx = -(-rr3 + xr*xr*rr5)
                  Txy = -(xr*yr*rr5)
                  Txz = -(xr*zr*rr5)
                  Tyx = Txy
                  Tyy = -(-rr3 + yr*yr*rr5)
                  Tyz = -(yr*zr*rr5)
                  Tzx = Txz
                  Tzy = Tyz
                  Tzz = -(-rr3 + zr*zr*rr5)

                  M_tot(i2+1,k2+1) = Txx
                  M_tot(i2+1,k2+2) = Txy
                  M_tot(i2+1,k2+3) = Txz
                  M_tot(i2+2,k2+1) = Tyx
                  M_tot(i2+2,k2+2) = Tyy
                  M_tot(i2+2,k2+3) = Tyz
                  M_tot(i2+3,k2+1) = Tzx
                  M_tot(i2+3,k2+2) = Tzy
                  M_tot(i2+3,k2+3) = Tzz

                  M_tot(k2+1,i2+1) = Txx
                  M_tot(k2+1,i2+2) = Txy
                  M_tot(k2+1,i2+3) = Txz
                  M_tot(k2+2,i2+1) = Tyx
                  M_tot(k2+2,i2+2) = Tyy
                  M_tot(k2+2,i2+3) = Tyz
                  M_tot(k2+3,i2+1) = Tzx
                  M_tot(k2+3,i2+2) = Tzy
                  M_tot(k2+3,i2+3) = Tzz

                  dir = dix*xr + diy*yr + diz*zr
                  qix = qixx*xr + qixy*yr + qixz*zr
                  qiy = qixy*xr + qiyy*yr + qiyz*zr
                  qiz = qixz*xr + qiyz*yr + qizz*zr
                  qir = qix*xr + qiy*yr + qiz*zr
                  dkr = dkx*xr + dky*yr + dkz*zr
                  qkx = qkxx*xr + qkxy*yr + qkxz*zr
                  qky = qkxy*xr + qkyy*yr + qkyz*zr
                  qkz = qkxz*xr + qkyz*yr + qkzz*zr
                  qkr = qkx*xr + qky*yr + qkz*zr
                  fid(1) = -xr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
     &                        - rr3_dir*dkx + 2.0d0*rr5_dir*qkx
                  fid(2) = -yr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
     &                        - rr3_dir*dky + 2.0d0*rr5_dir*qky
                  fid(3) = -zr*(rr3_dir*ck-rr5_dir*dkr+rr7_dir*qkr)
     &                        - rr3_dir*dkz + 2.0d0*rr5_dir*qkz
                  fkd(1) = xr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
     &                        - rr3_dir*dix - 2.0d0*rr5_dir*qix
                  fkd(2) = yr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
     &                        - rr3_dir*diy - 2.0d0*rr5_dir*qiy
                  fkd(3) = zr*(rr3_dir*ci+rr5_dir*dir+rr7_dir*qir)
     &                        - rr3_dir*diz - 2.0d0*rr5_dir*qiz
                  do j = 1, 3
c                     field(j,i) = field(j,i) + fid(j)*dscale(kk)
c                     field(j,k) = field(j,k) + fkd(j)*dscale(kk)
c                     fieldp(j,i) = fieldp(j,i) + fid(j)*pscale(kk)
c                     fieldp(j,k) = fieldp(j,k) + fkd(j)*pscale(kk)

                     field(j,l1) = field(j,l1) + fid(j)*dscale_dir(kk)
                     field(j,l3) = field(j,l3) + fkd(j)*dscale_dir(kk)
                     fieldp(j,l1) = fieldp(j,l1) + fid(j)*pscale(kk)
                     fieldp(j,l3) = fieldp(j,l3) + fkd(j)*pscale(kk)

                  end do
c               end if
            end if
         end do
      end do
c      else
c
c     periodic boundary for large cutoffs via replicates method
c
c       use_replica=.true.
c   NOTE TO SELF:  REMOVED if (use_replica) BECAUSE NCELL WAS
c   ALWAYS EQUAL TO ZERO AFTER PERIODIC BOUNDARY CONDITIONS WERE
c   IMPLEMENTED CORRECTLY. 

      deallocate (dscale_dir)
      deallocate (dscale_mut)
      deallocate (pscale)
      return
      end


      subroutine just_chk_rotpole
      call chkpole
c
c     rotate the multipole components into the global frame
c
      call rotpole
      return
      end

