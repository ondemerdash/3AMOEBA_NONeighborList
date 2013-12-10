c
c
c     ##########################################################
c     ##  COPYRIGHT (C)  2012  by  Liam Denis O'Suilleabhain  ##
c     ##              All Rights Reserved                     ##
c     ##########################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  combo.i  --  get the atom sites for m_pole interactions   ##
c     ##                                                            ##
c     ################################################################
c
c
c     moli      index for molecule number in empole
c     moli1     index for molecule number in empole
c     moli2     index for molecule number in empole
c     moli3 
c     latom     lower index of atom
c     hatom     higher index of atom
c     latom1    lower index of atom
c     hatom1    higher index of atom
c     latom2    lower index of atom
c     hatom2    higher index of atom
c     latom3    lower index of atom
c     hatom3    higher index of atom
c     lpole     lower index polarized atom
c     hpole     higher index polarized atom 
c     lpole1    lower index polarized atom
c     hpole1    higher index polarized atom 
c     lpole2    lower index polarized atom
c     hpole2    higher index polarized atom 
c     lpole3    lower index polarized atom
c     hpole3    higher index polarized atom 
c     num_pole  number of poles to evaluate
c     num_atom  number of atoms to evaluate
c     pnum      index of the polarizable site
c     anum      index of the atom site
c     ep2analyze
c     ep3analyze
c     r2b    
c     r3b_1
c     r3b_2
c     r3b_3
c     delta1
c
c    fd_2b has dims (3,npole3b,numbbodies)

      integer moli,moli1,moli2,moli3,
     &        count1,count2,count3,
     &        na1,na2,na3,np1,np2,np3,
     &        npole3b,natom3b,pnum,anum
      real*8 etest,delta1,delta2,r3b,delta3,
     &        total_field,u3b,u_2b,u2b,u1b
      logical body1,body2,body3
      common /combo/ moli,moli1,moli2,moli3,
     &               count1,count2,count3,
     &               np1,np2,np3,npole3b,natom3b,
     &               anum(10000),pnum(10000),
     &               etest,delta1,delta2,r3b,
     &               delta3,total_field(3,maxatm),
     &               u3b(3,maxatm),u_2b(3,1000,6000),
     &               u2b(3,maxatm),u1b(3,maxatm),
     &               body1,body2,body3
