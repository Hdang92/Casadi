C Simple example of use of MA57 package
      INTEGER LRHS,LFACT,LKEEP,LIFACT,LWORK
      PARAMETER (LKEEP=100,LRHS=10,LWORK=12,LFACT=1000,LIFACT=500)
      INTEGER IRN(10),JCN(10),IWORK(40),KEEP(LKEEP),IFACT(LIFACT)
      INTEGER ICNTL(40),INFO(40),N,NE,I,JOB
      REAL             A(30),WORK(LWORK),FACT(LFACT),RHS(LRHS,1),
     +                 CNTL(5),RINFO(20)
      EXTERNAL MA57I,MA57A,MA57B,MA57C


C Set default values for control parameters.
      CALL MA57I(CNTL,ICNTL)

C Ask for full printing from MA57 package
      ICNTL(5) = 4

C Read matrix and right-hand side
      READ (5,*) N,NE
      READ (5,*) (IRN(I),JCN(I),A(I),I=1,NE)
      READ (5,*) (RHS(I,1),I=1,N)

C Analyse sparsity pattern
      CALL MA57A(N,NE,IRN,JCN,LKEEP,KEEP,IWORK,ICNTL,INFO,RINFO)

C Factorize matrix
      CALL MA57B(N,NE,A,FACT,LFACT,IFACT,LIFACT,LKEEP,KEEP,
     +            IWORK,ICNTL,CNTL,INFO,RINFO)

C Solve the equations
      JOB = 1
      CALL MA57C(JOB,N,FACT,LFACT,IFACT,LIFACT,1,RHS,LRHS,
     +            WORK,LWORK,IWORK,ICNTL,INFO)
      END
