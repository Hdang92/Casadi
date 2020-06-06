clear all;
close;


%
% This is a Test Example that show the possible use of the
% matlab interfaces when User's is NOT available.
% At the end the user can find a simple example of how he can
% bypass the two matlab function ma57_F and ma57 and use
% the two mex functions ma57ab and ma57cd separately.
% We point out the in this case the user must not modify the
% parameters fact, ifact, and infoAB in output from ma57ab. 
% Any modification will cause an error that may be not detected or
% trapped by the function ma57cd and causing a segmentation error.
% The only values that can be safely changed between the call of ma57ab
% and the call of ma57cd is icntl(9). The value of icntl(9) gives the
% maximum number of steps that can be performed by the Iterative
% Refinement. If the user want to avoid the use of the Iterative Refinement
% in all is calls then he can set icntl(9) = 0 before the call of
% ma57cd.
%

control.order = 0; % AMD control.ordering
control.thres = 0.01;
control.stpv = 0;
nitre= 10;
load ma57_test1.mat; 
A = M;
xe = ones(length(A),1);
b  = A*xe;

% Error test 1
% matrix is NOT sparse
A1 = full(A);
% try
%     [y] = hsl_ma57(A1,b,control);
% catch 
%     errstr = lasterror();
%     str=strtrim(errstr.message);
%     str1=strtrim('The input matrix A must be square and sparse');
%     if (size(strfind(str,str1),2)==0)
%         error('Order - Failure at error test 1.1')
%     end
% end
try
    [struct,infoA,rinfoA] = ma57_factor(A1,control);

    [y,info,rinfo] = ma57_solve(A,b,struct);
catch
    errstr = lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('The input matrix A must be square and sparse');
    if (size(strfind(str,str1),2)==0)
        error('Order - Failure at error test 1.1')
    end
end
try
    [struct,infoA,rinfoA] = ma57_factor(A,control);

    [y,info,rinfo] = ma57_solve(A1,b,struct,nitre);
catch
    errstr = lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('Error in argument A. Expected sparse matrix.');
    if (size(strfind(str,str1),2)==0)
        error('Order - Failure at error test 1.2')
    end
end
clear control;
control.thres = 0.01;
control.stpv = 0;
control.nitre =10;
control.order = [];
try
    [struct,infoA,rinfoA] = ma57_factor(A,control);

    [y,info,rinfo] = ma57_solve(A1,b,struct,nitre);
catch
    errstr = lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('Error in argument control.order. Expected 1x1 matrix.');
    if (size(strfind(str,str1),2)==0)
        error('Order - Failure at error test 1.3')
    end
end
clear control;
control.thres = [];
control.stpv = 0;
control.order = 0;
try
    [struct,infoA,rinfoA] = ma57_factor(A,control);

    [y,info,rinfo] = ma57_solve(A1,b,struct,nitre);
catch
    errstr = lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('Error in argument control.thres. Expected 1x1 matrix.');
    if (size(strfind(str,str1),2)==0)
        error('Order - Failure at error test 1.4')
    end
end
clear control;
control.thres = 0.01;
control.stpv = [];
control.order = 0;
try
    [struct,infoA,rinfoA] = ma57_factor(A,control);

    [y,info,rinfo] = ma57_solve(A1,b,struct,nitre);
catch
    errstr = lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('Error in argument control.stpv. Expected 1x1 matrix.');
    if (size(strfind(str,str1),2)==0)
        error('Order - Failure at error test 1.5')
    end
end
% Error test 2
% number input output parameters

control.order = 0; % AMD control.ordering
control.thres = 0.01;
control.stpv = 0;
nitre= 10;
try
    [y,info,rinfo,z] = ma57_factor(A,b,control);
catch
    errstr = lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('Wrong # of arguments');
    if (size(strfind(str,str1),2)==0)
        error('Order - Failure at error test 2.1')
    end
end
try
    [struct,infoA,rinfoA] = ma57_factor(A,control);

    [y,info,rinfo,z] = ma57_solve(A,b,struct);
catch
    errstr = lasterror();
    str=strtrim(errstr.message);
    str1=strtrim(' Wrong # of output arguments ');
    if (size(strfind(str,str1),2)==0)
        error('Order - Failure at error test 2.2')
    end
end
try
    [struct,infoA,rinfoA] = ma57_factor(A,control);

    [y,info,rinfo] = ma57_solve(A,b,struct,nitre,b);
catch
    errstr = lasterror();
    str=strtrim(errstr.message);
    str1=strtrim(' Wrong # of input arguments ');
    if (size(strfind(str,str1),2)==0)
        error('Order - Failure at error test 2.3')
    end
end
try
    [struct] = ma57_factor(A,control);

%    [y,info,rinfo] = ma57_solve(A1,b,fact,ifact,cntl,icntl,infoA);
catch
    errstr = lasterror();
    str=strtrim(errstr.message);
    str1=strtrim(' Wrong # of arguments ');
    if (size(strfind(str,str1),2)==0)
        error('Order - Failure at error test 2.4')
    end
end
try
    [struct,infoA] = ma57_factor(A,control,b);
catch
    errstr = lasterror();
    str=strtrim(errstr.message);
    str1=strtrim(' Wrong # of arguments ');
    if (size(strfind(str,str1),2)==0)
        error('Order - Failure at error test 2.5')
    end
end
% Error test 3
% Incompatible rhs
b1 = b(1:10);
try
    [struct,infoA,rinfoA] = ma57_factor(A,control);

    [y,info,rinfo] = ma57_solve(A,b1,struct,nitre);
catch
    errstr = lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('The RHS has dimensions incompatible with the input matrix');
    if (size(strfind(str,str1),2)==0)
        error('Order - Failure at error test 1')
    end
end
% Error test 4
% Incompatible Struct is not a structure
b1 = b(1:10);stru = 1;
try
     [struct,info,rinfo] = ma57_factor(A,control);

    [y,info,rinfo] = ma57_solve(A,b,stru,nitre);
catch
    errstr = lasterror();
    str=strtrim(errstr.message);
    str1=strtrim('The third entry is not a structure');
    if (size(strfind(str,str1),2)==0)
        error('Order - Failure at error test 5.1')
    end
end
% WARNIBG messages
warning('off', 'all');
clear control;
control.order = 0; % AMD control.ordering
control.thres = 0.01;
control.stpv = 0;
nitre= 0;
    [struct] = ma57_factor(A,control);
    [y,info,rinfo] = ma57_solve(A,b,struct,nitre);

    msgstr = lastwarn();
    str=strtrim(msgstr);
    str1=strtrim('Iterative refinement is not executed');
    if (size(strfind(str,str1),2)==0)
        error('Order - Failure at error test 5.4')
    end

[n,m]=size(A);
A(:,n) = zeros(n,1);A(n,:) = zeros(1,n);

     [struct,info,rinfo] = ma57_factor(A,control);

%    [y,info,rinfo] = ma57_solve(A,b,struct,nitre);

    msgstr = lastwarn();
    str=strtrim(msgstr);
    str1=strtrim('The matrix is not full rank: look at info(25)');
    if (size(strfind(str,str1),2)==0)
        error('Order - Failure at error test 5.5')
    end

    nitre = 10;
  %  [struct,info,rinfo] = ma57_factor(A,control);

    [y,info,rinfo] = ma57_solve(A,b,struct,nitre);

    msgstr = lastwarn();
    str=strtrim(msgstr);
    str1=strtrim('The matrix is not full rank: look at info(25)');
    if (size(strfind(str,str1),2)==0)
        error('Order - Failure at error test 5.6')
    end

 A=M;
%      warning('on', 'all');
clear control;
%control.order = 0; % AMD control.ordering
control.thres = 0.01;
control.stpv = 0;
nitre= 10;
     [struct,info,rinfo] = ma57_factor(A,control);

%    [y,info,rinfo] = ma57_solve(A,b,struct,nitre);

    msgstr = lastwarn();
    str=strtrim(msgstr);
    str1=strtrim('control.order is not assigned and it is set to default value 0');
    if (size(strfind(str,str1),2)==0)
        error('Order - Failure at error test 5.1')
    end
    
clear control;
control.order = 0; % AMD control.ordering
%control.thres = 0.01;
control.stpv = 0;
nitre= 10;
     [struct,info,rinfo] = ma57_factor(A,control);

%    [y,info,rinfo] = ma57_solve(A,b,struct,nitre);

    msgstr = lastwarn();
    str=strtrim(msgstr);
    str1=strtrim('control.thres is not assigned and it is set to default value 0.01');
    if (size(strfind(str,str1),2)==0)
        error('Order - Failure at error test 5.2')
    end
clear control;
control.order = 0; % AMD control.ordering
control.thres = 0.01;
%control.stpv = 0;
nitre= 10;
     [struct,info,rinfo] = ma57_factor(A,control);

%    [y,info,rinfo] = ma57_solve(A,b,struct,nitre);

    msgstr = lastwarn();
    str=strtrim(msgstr);
    str1=strtrim('control.stpv is not assigned and it is set to default value 0.0');
    if (size(strfind(str,str1),2)==0)
        error('Order - Failure at error test 5.3')
    end
warning('on', 'all');    
fprintf('\nAll tests succeeded\n');
