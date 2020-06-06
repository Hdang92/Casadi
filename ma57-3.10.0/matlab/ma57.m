%
% hsl_ma57 package is composed of two procedures:
% 	hsl_ma57fact
% 	hsl_ma57solv 
%  i.e.
%
% [struct[,info,rinfo {optional}]] = hsl_ma57fact(A[,control {optional}])
% [x[,info,rinfo {optional}]]=hsl_ma57solve(A,b,struct [,nitre {optional}])
% 
% The user must call the two  functions  separately  following  the rules  
% for the optional argin and argout described below. 
% The use of two separate mex functions allows the user  to
% solve several systems with the same  matrix but  different  right-hand 
% sides that are not available simultaneously and so it is not 
% possible to solve for a block of right-hand  sides.
% In  this case, the call to the factorization  is executed only once.  
% However, the  user  
%                    MUST  NOT MODIFY
% the  values of  struct after  the call to hsl_ma57fact 
% and should pass these variables unchanged to hsl_ma57solve.
%
% The  value  of  nitre  gives the maximum number of  steps that  can be  
% performed by the Iterative Refinement.  
% If  the  user  wants to  avoid the  use  of  
% Iterative  Refinement  then he or she can set nitre = 0  before the 
% call to hsl_ma57solve. 
% The new value of nitre will overwrite the default value of  10.
% The parameter nitre is optional and can be omitted.
% 
% hsl_ma57fact  performs an LDL'  factorization  of  a  permutation of a
% sparse A, using the  HSL_2011  routines  MA57AD, MA57BD, and MA57ED.
% These   routines   implement   a   multifrontal  algorithm   for   the
% factorization of A. The matrix L is a unit lower triangular matrix and
% D is block diagonal  with blocks of order 1 or 2.  hsl_ma57  uses only 
% the diagonal and the lower triangle of A. 
% The upper triangle is assumed to be  the transpose of  the lower. 
% Before  the factorization  starts, scaling is performed using the MC64  
% package of HSL_2011.
% In order to minimize the  fill-in  during  the  factorization process,  
% the user can choose to reorder the matrix A either by  Approximate 
% Minimum  Degree (AMD) (the default value) or, if the user wants to use 
% a permutation p, he/she can call the interfaces using A(p,p) and
% setting the control.order parameter to a value  > 1 ( see optional ARGIN
% below). 
%  
% For  stability, all the  pivots are  tested numerically.   However, in
% order to avoid excessive disruption of the AMD or the user's reorderings
% threshold pivoting is used. Therefore,  denoting  the  entries  to  be
% processed  at step  k of the  factorization by f_{ij}, a 1x1 numerical 
% pivot will be used if
%  
% |f_{kk} | \ge u \max_{j \ne k} |f_{kj} |,
%  
% or a 2x2 numerical pivot will be used if
%  
% (|f_{k,k}|   |f_{k,k+1}|  ) (\max_{j\ne k,k+1} |f_{k,j}|  )        (1)
% (                         ) (                             ) <= u^-1( )
% (|f_{k+1,k}| |f_{k+1,k+1}|) (\max_{j\ne k,k+1} |f_{k+1,j}|)        (1)
%  
% where u  is  a threshold  parameter given  in the control.thres 
% parameter (the default value is 0.01), see the ARGIN optional below.
% During stage k of a multifrontal implementation
% of Gaussian elimination, some entries may  not be eligible to be chosen 
% as a pivot.  In  this  case, the  algorithm  will  proceed  using   the 
% eligible variables only and  will choose  the rejected pivot at a later 
% stage when  the candidate  entry will  become eligible.   
% Additional details about the theory of the multifrontal method and the  
% numerical pivoting strategy can be found in [1] and [2].
%  
% [struct[,info,rinfo {optional}]] = hsl_ma57fact(A[,control {optional}])
%   ARGIN
%  
%       A  input matrix (only the lower triangular part will be used).
%
%   OPTIONAL ARGIN
%  
%   control.order  
%           Allows the user to choose a reordering strategy for the pivot
%           sequence:
%              control.order <= 0    ==> AMD ordering is used
%              control.order  > 1    ==> User's reordering p is used and 
%                                        the user must pass the 
%                                        permuted matrix A(p,p) in the 
%                                        input
%           Its default value is 0.
%  
%   control.thres  
%           The value of the threshold u (see ME57 Specification sheet).
%           Its default value is 0.01.
%  
%   control.stpv  
%           If stpv is set to a value greater than zero, static pivoting  
%           is invoked.   This  means that if, at  any   stage  of   the 
%           multifrontal elimination, the threshold criterion defined by  
%           thresh prevents us choosing a  pivot, then  we  first weaken 
%           the threshold by successively trying values one tenth of the 
%           previous until a  threshold of sqrt(thresh)*stvp is reached.   
%           If the pivot is still not chosen, but involves non  eligible  
%           variables, then  the  diagonal  entries are replaced by  
%           control.stpv times the largest modulus of an entry in the  
%           original matrix and are used as 1 by 1 pivots in the 
%           factorization. The default value is 0. 
%           We advise the user interested in this option to use  first a 
%           value close to sqrt(eps) for stpv. Iterative refinement must
%           be used but convergence is not guaranteed.
%
% If any of the fields of control is omitted the default value is 
% automatically set. Incorrect values or empty values produce an error 
% message.
%
%   OPTIONAL ARGOUT
%  
%   info and rinfo 
%                 hold the statistics of the analysis and factorization
%                 execution and other internal information: see the
%                 detailed description below.
%
% hsl_ma57solve uses the  HSL_2011  routines  MA57CD and MA57DD in order 
% to compute the solution(s) corresponding to the right-hand side(s) 
% held in b. 
%
% [x[,info,rinfo {optional}]]=hsl_ma57solve(A,b,struct [,nitre {optional}])
%
%       A  input matrix (only the lower triangular part will be used).
%
%       b  input right-hand side(s)
%
%  struct  input from the previous call to hsl_ma57fact
%
%   OPTIONAL ARGIN
%  
%   nitre  the user may choose the maximum number of iterative refinememnt
%          steps (the default value is 10). If the user sets nitre to 0 
%          iterative refinement is not used.
%          If b is a matrix with more than one column then nitre is 
%          automatically set to 0.
%  
%   ARGOUT
%  
%       x   Computed  solution(s).   If  b  has  several  columns, the 
%           corresponding  columns of  x will hold   the  corresponding
%           solutions.
%           If b is a vector, iterative refinement will be used to give
%           the best possible solution.
%  
%  
%   OPTIONAL ARGOUT
%  
%   info and rinfo 
%                 Hold the statistics of the analysis and factorization
%                 from the previous call to hsl_ma57fact plus the statistics
%                 of the solution phase (see description below).
% 
% 
% If A is not symmetric hsl_ma57fact and hsl_ma57solve will use only the 
% lower triangular part of A.
% 
% ERRORS and WARNINGS
% 
%   If A is not square or sparse an ERROR message is printed. If a 
%   fatal error occurs during the factorization an error message is
%   printed. Fatal errors include the wrong number of varargin and/or 
%   varargout parameters. If it is impossible to allocate  memory 
%   a fatal error message will be printed.
%  
%   If the numerical rank < length(A) a WARNING message is printed.
% 
% 
%
%  Description of the info and rinfo values:
%
%  info(1)  will be 0 after a successful execution; 
%           a positive value indicates that the matrix A has 
%           numerical pseudo rank < length(A). The solution can be 
%           incorrect, and info(25) holds the numerical pseudo rank 
%           value.
%  info(2)  to info(4)  are only used internally.
%  info(5)  Forecast number of real in the factors.
%  info(6)  Forecast number of integers in the factors.
%  info(7)  Forecast maximum frontal matrix size.
%  info(8)  Forecast number of nodes in the assembly tree.
%  info(9)  Forecast of the length of fact array(real)  
%               without numerical pivoting.
%  info(10) Forecast of the length of ifact array(integer) 
%           without numerical pivoting.
%  info(11) Length of fact required for a successful
%               completion of the numerical phase allowing
%           data compression (without numerical pivoting).
%  info(12) Length of ifact required for a successful 
%           completion of the numerical phase allowing
%           data compression (without numerical pivoting).
%  info(13) Number of data compresses.
%  info(14) Number of entries in factors.
%  info(15) Storage for real data in factors.
%  info(16) Storage for integer data in factors.
%  info(17) Minimum length of fact required for a successful
%           completion of the numerical phase.
%  info(18) Minimum length of ifact required for a successful 
%           completion of the numerical phase.
%  info(19) Minimum length of fact required for a successful
%           completion of the numerical phase without 
%           data compression.
%  info(20) Minimum length of ifact required for a successful
%           completion of the numerical phyase without data compression.
%  info(21) Order of the largest frontal matrix.
%  info(22) Number of 2x2  pivots.
%  info(23) Total number of fully-summed variables passed to
%           the father node because of numerical pivoting.
%  info(24) Number of negative eigenvalues.
%  info(25) Rank of the factorization of the matrix.
%  info(26) Only used internally
%  info(27) Pivot step where static pivoting commences and
%           it is set to zero if no modification performed.
%  info(28) Number of data compresses on real data structures.
%  info(29) Number of data compresses on integer data structures.
%  info(30) Number of steps performed by Iterative Refinement.
%  info(31) Number of block pivots in factors (see info(8)).
%  info(32) Number of zeros in the triangle of the factors.
%  info(33) Number of zeros in the rectangle of the factors.
%  info(34) Number of zero columns in rectangle of the factors.
%  info(35) Number of pivots chosen in the static pivoting phase.
%  info(36) to info(40) are only used internally.
% 
% rinfo(1)  Forecast number of floating_point additions
%           for the assembly processes. 
% rinfo(2)  Forecast number of floating_point operations
%           to perform the elimination operations 
%           counting a multiply-add pair as two operations 
% rinfo(3)  Number of floating_point additions 
%           for the assembly processes.
% rinfo(4)  Number of floating_point operations
%           to perform the elimination operations
%           counting multiply-add pair as two operations.
% rinfo(5)  Only used internally
% rinfo(6)  Omega1 (backward error) (see MA57 spec).
% rinfo(7)  Omega2 (backward error) (see MA57 spec).
% rinfo(8) to rinfo(10) are only used internally
% rinfo(11) Cond number for Omega1
% rinfo(12) Cond number for Omega2
% rinfo(13) Estimated forward error
% rinfo(14) to rinfo(15) are only used internally.
% rinfo(16) Minimum value of the scaling factor (MC64).
% rinfo(17) Maximum value of the scaling factor (MC64).
% rinfo(18) Maximum modulus of matrix entry.
% rinfo(19) and rinfo(20) are only used internally.
% rinfo(21) Analysis phase cputime (sec).
% rinfo(22) Factorization phase cputime (sec).
% rinfo(23) Solution phase cputime (sec) without iterative refinement.
% rinfo(24) Solution phase cputime (sec) with iterative refinement.
% 
% 
% References
%  
% [1] I.S. Duff ``MA57: a new  code for the solution of sparse symmetric
% indefinite   systems'',    Technical   Report   RAL-TR-2002-024,   ACM
% Trans. Math. Software 30 (2004), 118-144.
%  
% [2] I.S. Duff and S. Pralet, ``Strategies for scaling and pivoting for
% sparse symmetric indefinite problems'', SIAM. J. Matrix Anal. & Appl. 27, 
% (2005), 313-340.
%
