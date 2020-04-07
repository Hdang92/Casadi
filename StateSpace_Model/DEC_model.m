function [ ode, tc ] = DEC_model(tt,mp,H,J)

%% import casadi
    addpath('/home/hdang/Matlab/casadi-linux-matlabR2014b-v3.5.1')
    import casadi.*

%% initialize variables
    % Fix other parameters
        g0=9.81;
        alpha = 0.05;
        kp = mp*g0*H/180*pi;

    % Declare model variables
        x1=SX.sym('x1');
        x2=SX.sym('x2');
        x3=SX.sym('x3');
        x4=SX.sym('x4');
        x5=SX.sym('x5');
        x6=SX.sym('x6');
        u_bar=SX.sym('u_bar');
        x=[x1;x2;x3;x4;x5;x6;u_bar];
        u=SX.sym('u'); 

    % Process noise
        b1=SX.sym('b1');
        b7=SX.sym('b7');
        b=[b1;b7];

    % Parameters to estimate
        wfg=SX.sym('wfg');
        wff=SX.sym('wff');
        lambda=SX.sym('lambda');
        sg=SX.sym('sg');
        kd=SX.sym('kd');
        tau=SX.sym('tau');
        kf=SX.sym('kf');
        fs=SX.sym('fs');
        params_d=[wfg;wff;lambda;sg;kd;tau;kf;fs];

    % Internal variables
        uu = -(wfg+sg)*(x1+b1) -wff*(x6) +sg*u +(x3);
        tc = Function('tc',{x4,x5,kd,tau},{kp*(x5) + kd*((x4)-(x5))/tau});

%% Model equations
    xdot = [(x2);
            kd/J/tau*(x4) + 1/J*(kp-kd/tau)*(x5) + mp*g0*H/180*pi/J*(x1+b1); % unity conversion here!
            kf/fs*(kp-kd/tau)*(x5) + kf*kd/tau/fs*(x4) - 1/fs*(x3);
            uu/tau - (x4)/tau;
            ((x4)-(x5))/tau;
            1/2*sqrt( ((u-(u_bar+b7))/tt-lambda)^2 + alpha*lambda^2 ) - 1/2*sqrt( ((u-(u_bar+b7))/tt+lambda)^2 + alpha*lambda^2 ) + (u-(u_bar+b7))/tt;
            (u-(u_bar+b7))/tt];

%% construct dae structure
    ode = struct('x',x,'p',[u;params_d;b],'ode',xdot);
        