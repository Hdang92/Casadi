%
% This is a Test Example that show the possible use of the
% matlab interfaces.
% 
% We point out the user MUST NOT MODIFY the values in struct
% in output from ma57_factor. 
% Any modification will cause an error that may be not detected or
% trapped by the function ma57_solve and causing a segmentation error.
% See last example.
% If the user want to avoid the use of the Iterative Refinement
% in all is calls then he can set nitre = 0 before the call of
% ma57_solve.
%
clear all;
close;

fprintf('\n\n Example 1\n');
fprintf('\nTEST THE SOLVER WITH AMD ORDERING AND ERROR ANALYSIS (M)\n');


control.order = 0; % AMD control.ordering
control.thres = 0.01;
control.stpv = 0;

%load M.mat;  
load ma57_test1.mat;
A = M;
xe = ones(length(A),1);
b  = A*xe;

fprintf('\n\n Example 2\n');
fprintf('\nTEST THE SOLVER WITH AMD ORDERING AND ERROR ANALYSIS (Wathen)\n');
C = gallery('wathen',2,2);
A = sparse(C);
xe = ones(length(A),1);
b  = A*xe;

control.order = 0; % AMD control.ordering
control.thres = 0.01;
control.stpv = 0;

[str,info,rinfo] = ma57_factor(A,control);

[y,info,rinfo] = ma57_solve(A,b,str);

fprintf('\n Relative Norm of the errors = %7.1e\n',norm(xe-y)/norm(xe));
fprintf('\n info = %i\n',info(1));

%
% Print Statistics on line
%
%ma57_stats(control.order,info,rinfo);

%load DUALC1_C_I.mat;
load ma57_test2.mat;
xe = ones(length(A),1);


fprintf('\n\n Example 3\n');
fprintf('\nTEST THE SOLVER WITH AMD ORDERING AND MultipleRHS\n');

control.order = 0; % AMD control.ordering
control.thres = 0.01;
control.stpv = 0;

Xe = [xe ,(xe+1)];

B1 = A*Xe;


[str,info,rinfo] = ma57_factor(A,control);

for k = 1:2
    y = ma57_solve(A,B1(:,k),str);
    X(:,k) = y;
end;
fprintf('\n Relative Norm of the errors = %7.1e\n',norm(Xe-X)/norm(Xe));


fprintf('\n\n Example 4\n');

fprintf('\nTEST THE SOLVER WITH AMD ORDERING AND ERROR ANALYSIS\n');

control.order = 0; % AMD control.ordering
control.thres = 0.01;
control.stpv = 0;
nitre = 5;
xe = ones(length(A),1);
b  = A*xe;


[str,info,rinfo] = ma57_factor(A,control);

[y,info,rinfo] = ma57_solve(A,b,str,nitre);

fprintf('\n Relative Norm of the errors = %7.1e\n',norm(xe-y)/norm(xe));
fprintf('\n info = %i\n',info(1));


fprintf('\n\n Example 5\n');

fprintf('\nTEST THE SOLVER WITH AMD ORDERING AND STATIC PIVOT\n');

%load cont201;
load ma57_test3.mat;
A = Problem.A;
xe = ones(length(A),1);
b  = A*xe;

control.order = 0; % AMD control.ordering
control.thres = 0.01;
control.stpv = 1e-08;


[str,info,rinfo] = ma57_factor(A,control);

[y,info,rinfo] = ma57_solve(A,b,str);
    
fprintf('\n Relative Norm of the errors = %7.1e\n',norm(xe-y)/norm(xe));



fprintf('\n\n Example 6\n');

fprintf('\nTEST THE SOLVER WITH USER ORDERING AND ERROR ANALYSIS\n');
fprintf('\nExample of use of ma57_stats\n'); 
load ma57_test2.mat;

control.order = 5; % user control.ordering
control.thres = 0.01;
control.stpv = 0;
p =symrcm(A);
A = A(p,p);
xe = ones(length(A),1);
b  = A*xe;


[str,info,rinfo] = ma57_factor(A,control);

[y,info,rinfo] = ma57_solve(A,b,str);

fprintf('\n Relative Norm of the errors = %7.1e\n',norm(xe-y)/norm(xe));
fprintf('\n info = %i\n',info(1));



%
% Print Statistics on file 'ma57_test.log'
%
ma57_stats(control.order,info,rinfo,'ma57_test.log');

fprintf('\n\n Example 7\n');
fprintf('\nTEST THE SOLVER WHEN WE SIMULATE A FAILURE IN ma57_factor\n');


control.order = 0; % AMD control.ordering
control.thres = 0.01;
control.stpv = 0;

%load M.mat;  
load ma57_test1.mat;
A = M;
xe = ones(length(A),1);
b  = A*xe;

% simulation of an error in the factorization phase
% the returning solution is NaN array

[str,info,rinfo] = ma57_factor(A,control); str.info(1) =-100;

try
   [y,info,rinfo] = ma57_solve(A,b,str);
catch
   errstr = lasterror();
   str=strtrim(errstr.message);
   str1=strtrim('Previous factorization failed');
   if (size(strfind(str,str1),2)==0)
      error('Unexpected error message')
   end
end
fprintf('Test OK\n')
