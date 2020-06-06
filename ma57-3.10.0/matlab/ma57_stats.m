function hsl_ma57_stats(order,info,rinfo,file)


if (nargin < 4),
    fid = 1;
else
    fid = fopen(file,'a+');
end

if (order >= 1)
    fprintf(fid,'User''s reordering in use\n\n');
end
if (order <= 0)
    fprintf(fid,'AMD reordering in use\n\n');
end



if (info(1) > 0)
    fprintf(fid,'Numerical rank                             = %i \n',info(25));
end


fprintf(fid,' \n \nhsl_ma57 Statistics of the Analysis phase \n \n');

fprintf(fid,'Analysis Phase cputime                              : %9.3e sec \n', rinfo(21));
fprintf(fid,'Main error (if zero ma57 run successfully)          =  %i \n',info(1));
fprintf(fid,'Forecast number of real in the factors              =  %i \n',info(5));
fprintf(fid,'Forecast number of integers in the factors          =  %i \n',info(6));
fprintf(fid,'Forecast maximum frontal matrix size                =  %i \n',info(7));
fprintf(fid,'Forecast number of nodes in the assembly tree       =  %i \n',info(8));
fprintf(fid,'Forecast of the length of fact array(real) \n');
fprintf(fid,'without numerical pivot                             =  %i \n',info(9));
fprintf(fid,'Forecast of the length of ifact array(integer) \n');
fprintf(fid,'without numerical pivot                             =  %i \n',info(10));
fprintf(fid,'Length of fact required for a successful \n');
fprintf(fid,'completion of the numerical phase allowing\n');
fprintf(fid,'data compression (without numerical pivoting)       =  %i \n',info(11));
fprintf(fid,'Length of ifact required for a successful \n');
fprintf(fid,'completion of the numerical phase allowing\n');
fprintf(fid,'data compression (without numerical pivoting)       =  %i \n',info(12));
fprintf(fid,'Number of data compresses                           =  %i \n',info(13));
fprintf(fid,'Forecast number of floating point additions\n');
fprintf(fid,'for the assembly processes                          =  %9.3e \n',rinfo(1));
fprintf(fid,'Forecast number of floating point operations \n');
fprintf(fid,'to perform the elimination operations \n');
fprintf(fid,'counting multiply-add pair as two operations        =  %9.3e \n',rinfo(2));
fprintf(fid,'\n \n');

fprintf(fid,' \n \n hsl_ma57 Statistics of the numerical factorization phase \n \n');

fprintf(fid,'Factorization cputime                               : %9.3e sec \n',rinfo(22));

fprintf(fid,'Number of entries in factors                        = %i \n',info(14));
fprintf(fid,'Storage for real data in factors                    = %i \n',info(15));
fprintf(fid,'Storage for integer data in factors                 = %i \n',info(16));
fprintf(fid,'Minimum length of fact required for a successful \n');
fprintf(fid,'completion of the numerical phase                   = %i \n',info(17));
fprintf(fid,'Minimum length of ifact required for a successful \n');
fprintf(fid,'completion of the numerical phase                 = %i \n',info(18));
fprintf(fid,'Minimum length of fact required for a successful \n');
fprintf(fid,'completion of the numerical phase without \n');
fprintf(fid,'data compression                                    = %i \n',info(19));
fprintf(fid,'Minimum length of ifact required for a successful \n');
fprintf(fid,'completion of the numerical phyase without \n');
fprintf(fid,'data compression                                    = %i \n',info(20));
fprintf(fid,'Order of the largest frontal matrix                 = %i \n',info(21));
fprintf(fid,'Number of 2x2 numerical pivots                      = %i \n',info(22));
fprintf(fid,'Total number of fully-summed variables passed to \n');
fprintf(fid,'the father node because of numerical pivot          = %i \n',info(23));
fprintf(fid,'Number of negative eigenvalues                      = %i \n',info(24));
fprintf(fid,'Rank of factorization of the matrix                 = %i \n',info(25));
fprintf(fid,'Pivot step where static pivot commences and \n');
fprintf(fid,'it is set to zero if no modification performed      = %i \n',info(27));
fprintf(fid,'Number of data compresses on real \n');
fprintf(fid,'data structures                                     = %i \n',info(28));
fprintf(fid,'Number of data compresses on integer \n');
fprintf(fid,'data structures                                     = %i \n',info(29));
fprintf(fid,'Number of block pivots in factors                   = %i \n',info(31));
fprintf(fid,'Number of zeros in the triangle of the factors      = %i \n',info(32));
fprintf(fid,'Number of zeros in the rectangle of the factors     = %i \n',info(33));
fprintf(fid,'Number of zero column in rectangle of the factors   = %i \n',info(34));
fprintf(fid,'Number of pivots chosen in static pivoting phase    = %i \n',info(35));
fprintf(fid,'Number of floating point additions \n');
fprintf(fid,'for the assembly processes                          = %9.3e \n',rinfo(3));
fprintf(fid,'Number of floating point operations \n');
fprintf(fid,'to perform the elimination operations \n');
fprintf(fid,'counting multiply-add pair as two operations        = %9.3e \n',rinfo(4));
fprintf(fid,'Minimum value of the scaling factor (MC64)          = %9.3e \n',rinfo(16));
fprintf(fid,'Maximum value of the scaling factor (MC64)          = %9.3e \n',rinfo(17));
fprintf(fid,'Maximum modulus of matrix entry                     = %9.3e \n\n\n',rinfo(18));

fprintf(fid,'\n\nhsl_ma57 statistics of the solution phase \n \n');

fprintf(fid,'Solution phase cputime without Iter. Ref.  : %9.3e sec \n', rinfo(23));
fprintf(fid,'Solution phase cputime with    Iter. Ref.  : %9.3e sec \n', rinfo(24));

fprintf(fid,'Number of steps performed by Iterative Refinement   = %i \n\n\n',info(30));
fprintf(fid,'Omega1 (backward error)                             = %9.3e \n',rinfo(6));
fprintf(fid,'Omega2 (backward error)                             = %9.3e \n',rinfo(7));
fprintf(fid,'Cond number for Omega1                              = %9.3e \n',rinfo(11));
fprintf(fid,'Cond number for Omega2                              = %9.3e \n',rinfo(12));
fprintf(fid,'Estimated forward error                             = %9.3e \n',rinfo(13));

if (nargin == 4),
    fclose(fid);
end;

