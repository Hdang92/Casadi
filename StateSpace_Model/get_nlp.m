function [ nlp, m, lbw, ubw, lbg, ubg, w0, R ] = get_nlp(ode, u1, y1, out, n_amp, N, tt)

% Inputs: - ode:
%         - u1: 
%         - y1: 
%         - out: 
%         - n_amp: number of amplitudes
%         - N: number of discretization 
%         - tt:
% Outputs:- nlp:
%         - m:
%         - lbw:
%         - ubw:
%         - lbg:
%         - ubg:
%         - w0:
%         - R:
%% import casadi
    addpath('/home/haidang/Matlab/casadi-linux-matlabR2014b-v3.4.5')
    import casadi.*

%% Formulation of NLP
    % Implementation of system dynamics in Casadi
        opts = struct('tf', tt,'number_of_finite_elements',15, 'jit', true);
        F = integrator('F', 'rk', ode, opts);
               
    nu = 1; % number of inputs
    nx = 7; % number of model variables
    nb = 2; % number of state noise terms

% New NLP variable for the parameters
    Wfg=MX.sym('Wfg');
    Wff=MX.sym('Wff');
    Lambda=MX.sym('Lambda');
    SG=MX.sym('SG');
    Kd=MX.sym('Kd');
    Tau=MX.sym('Tau');
    Kf=MX.sym('Kf');
    Fs=MX.sym('Fs');
    params_D=[Wfg; Wff; Lambda; SG; Kd; Tau; Kf; Fs];

% Start with an empty NLP
    w={};     % collocation of NLP across all time points
    w0 = [];  % initial conditions for NLP
    lbw = []; % lower bound w
    ubw = []; % upper bound w
    g={};     % equality constraints
    lbg = []; % lower bound g
    ubg = []; % upper bound g

    x1 = {};
    x_0 = {};
    x = {};
    B = {};
    U = {};
% initialization of sample time count
    n = 0;
    
% "Lift" initial conditions
    %J=0;
    
% Formulate the NLP
for i=1:n_amp
    % drop state memory when new amplitude starts
        Xk = MX.sym(['X' num2str(n)], nx); % state memory
        Bk = MX.sym(['B' num2str(n)] ,nb); % noise memory
        

        w = [w, {Xk}, {Bk}];                        % symbolic nlp vector
        x1 = [x1,{Xk(1,1)}];
        
        x_0 = [x_0,{Xk(2,1)},{Xk(3,1)},{Xk(4,1)},{Xk(5,1)},{Xk(6,1)},{Xk(7,1)}];
        
        B = [B,{Bk}];
        %%%%%%%%%% RP - changed initial X values from zero to inf
        lbw = [lbw; -inf*ones(nx,1); zeros(nb,1)]; % lower bounds of initial state and noise values
        ubw = [ubw; inf*ones(nx,1); zeros(nb,1)];  % upper bounds of initial state and noise values
        
        w0 = [w0; zeros(nx,1); zeros(nb,1)];       % numeric initial state and noise values

    for k=0:N-1
        % New NLP variable for the control
            Uk = MX.sym(['U_' num2str(n)]);
            w = [w, {Uk}];
            U = [U,{Uk}];
            lbw = [lbw; -inf];
            ubw = [ubw;  inf];
            w0 = [w0;  u1(i,k+1)];

        % Integrate till the end of the interval; for continuity condition
            Fk = F('x0', Xk, 'p', [Uk; params_D; Bk]);
            Xk_end = Fk.xf;
            
            %%%%%%%%%%%% RP - changed collocation 'counting' to improve readability
            n = n+1; % increasing time sample count
            %%%%%%%%%%%%
            
            Xk = MX.sym(['X_' num2str(n)], nx); % New NLP variable for state at end of interval
            Bk = MX.sym(['B_' num2str(n)], nb); % New NLP variable for the process noise
            w = [w, {Xk}, {Bk}];
            x1 = [x1, {Xk(1,1)}];
            B = [B, {Bk}];
            %x = [x,{Xk(2:7,1)}];

            lbw = [lbw; -inf*ones(nx,1) ; -inf*ones(nb,1)];
            ubw = [ubw;  inf*ones(nx,1); inf*ones(nb,1)];
            
            sim_states(:,1) = out(:,i,k+2);     % initial values for this sample time point from simulation
            w0 = [w0; sim_states; zeros(nb,1)]; %%%%%%%%%%%%% RP - changed initial state values from zero to previously simulated values

        % Add equality constraint 
            g = [g, {Xk_end-Xk}];
            lbg = [lbg; zeros(nx,1)];
            ubg = [ubg; zeros(nx,1)];
    end
    % increase sample time count for next cycles initial conditions
        n = n+1; % increasing time sample count

end


%% Add parameters constraints
    lbw = [lbw; 0; 0; 0; 0; 0; 0.025; 0; 0.025];
    ubw = [ubw; 1; 1; inf; inf; inf; 60.5; inf; 60.5];

% Decision variables and their initialization
    params=[vertcat(w{:});params_D];

% Weights; more trust into component --> higher weight
% weight is proportional to the inverse of the variance
    wy = 5;
    wu = 5;
    wb = 0.5;

% input and output measurements
    % measurement values as row vector
    y1 = reshape(y1', [(N+1)*n_amp,1]);
    y1 = y1';
    u1 = u1(:, 1:end-1);
    u1 = reshape(u1', [N*n_amp,1]);
    u1 = u1';
    m=[u1,y1];                                % Concatenate measurements values

    % symbolic measuremnt row vector
    U1=MX.sym('U1',nu,N*n_amp);
    Y1=MX.sym('Y1',1,(N+1)*n_amp);
    
    M=[U1,Y1];                                % Concatenate measurements in symbolic

% Objective
    ey = abs(y1'-vertcat(x1{:}));
    eu = abs(u1'- vertcat(U{:}));
    eb = abs(vertcat(B{:}));
    ex = abs(vertcat(x_0{:}));
    R = [sqrt(wu).*eu;sqrt(wy).*ey;sqrt(wb).*eb;ex];
% Objective function
    J = 0.5*wy*dot(ey,ey) + 0.5*wu*dot(eu,eu) + 0.5*wb*dot(eb,eb)+ 0.1*dot(ex,ex);
% create nlp structure
    nlp = struct('f', J, 'x', params, 'g', vertcat(g{:}), 'p', M);
    
    
    
    