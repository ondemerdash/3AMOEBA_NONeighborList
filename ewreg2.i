
c     ##############################################################
c     ##                                                          ##
c     ##  ewreg.i  --  exponential factors for regular Ewald sum  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     maxvec    maximum number of k-vectors per reciprocal axis
c
c     ejc       exponental factors for cosine along the j-axis
c     ejs       exponental factors for sine along the j-axis
c     ekc       exponental factors for cosine along the k-axis
c     eks       exponental factors for sine along the k-axis
c     elc       exponental factors for cosine along the l-axis
c     els       exponental factors for sine along the l-axis
c
c
      integer maxvec
      parameter (maxvec=15)
