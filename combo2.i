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

      integer moli,moli1,moli2,moli3,moli4,
     &        na1,na2,na3,na4,np1,np2,np3,np4,
     &        npole3b,natom3b,pnum,anum
      real*8  ep2analyze,ep3analyze,r2b,r3b,delta1,
     &        delta2,delta3,fd_2b,fd_1b,fd
      logical use_liam,body1,body2,body3
      common  /combo/ moli,moli1,moli2,moli3,moli4,
     &               na1,na2,na3,na4,np1,np2,np3,np4,
     &               npole3b,natom3b,pnum(11000000),
     &               anum(1700000),ep2analyze(650000),
     &               ep3analyze(11000000),r2b(maxatm),
     &               r3b(maxatm*maxatm),delta1,delta2,delta3,
     &               fd_2b(3,maxatm,220),fd_1b(3,maxatm),
     &               fd(3,maxatm),use_liam,body1,body2,body3

