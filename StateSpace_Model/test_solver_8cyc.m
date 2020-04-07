%% initialize
addpath(genpath('/home/hdang/Matlab/LA_toolbox'))
addpath('/home/hdang/Matlab/casadi-linux-matlabR2014b-v3.5.1')

close all
clc

import casadi.*

%% get input and output data
    load('cycles_eyes_closed.mat')
    n=0;
    for k=1:4
        for m=2:3
            n=n+1;
%             u_t = mean(data(k).FS(m,:),1);
%             u1(n,:) = u_t - mean(u_t);
            y_t = mean(data(k).BS(m,:),1);
            y1(n,:) = y_t - mean(y_t);
        end
    end
    clear u_t n k m data
    sr = 100; % sampling rate
    
    [x,xi]=pseudogen3(5,[0 0 1 -1 -1],[2 0 1 0 2],25);
    u1(1,:) = ( xi/(max(xi)-min(xi)) ) * 1;
    u1(2,:) = ( xi/(max(xi)-min(xi)) ) * 1;
    u1(3,:) = ( xi/(max(xi)-min(xi)) ) * 2;
    u1(4,:) = ( xi/(max(xi)-min(xi)) ) * 2;
    u1(5,:) = ( xi/(max(xi)-min(xi)) ) * 4;
    u1(6,:) = ( xi/(max(xi)-min(xi)) ) * 4;
    u1(7,:) = ( xi/(max(xi)-min(xi)) ) * 8;
    u1(8,:) = ( xi/(max(xi)-min(xi)) ) * 8;
    
%     figure(10)
%     t=(1:length(u1))/100;
%     plot(t,u1'); hold on

%% reduction of data points
    NSamp = 25; % downsampling factor

    u_t=u1; y_t=y1; clear u1 y1
    for k=1:size(u_t,1);
        u1(k,:) = u_t(k,1:NSamp:end);
        y1(k,:) = y_t(k,1:NSamp:end);
    end
    clear u_t k
    
    % define data structure and time intervals
        tt = NSamp/sr; ff=1/tt;
        N = size(u1,2)-1;
        T = 60.5-tt;
        n_cyc = size(u1,1);
        
%     t2 = (0:length(u1)-1)/ff;
%     figure(10); plot(t2,u1','x'); title('input')


%% define body parameters
        mp = 69;
        H = 0.93;
        J = (72.8)*pi/180;
        
%% get model
    [ dae, tc ] = DEC_model(tt,mp,H,J);
    
%% Simulate model with fixed parameters
    % Implementation of system dynamics in Casadi
        opts = struct('tf', tt,'number_of_finite_elements',15, 'jit', true);
        F = integrator('F', 'rk', dae, opts);
        
        params_Start = [0.55; 0.25; 0.3; 1; 0.3*mp*9.81*H/180*pi; 0.08; 0.1; 20];
                        % Wfg; Wff; Lambda; SG; Kd; Tau; Kf; Fs;
    % evaluate F for every time step
        out = zeros(7,size(u1,1),size(u1,2));
        for j=1:size(u1,1) % number of cycles
            X = zeros(length(dae.x),1); % set all states to zero (start of new cycle)
            for i=1:N % samples within a cycle
                x_out = F('x0',X,'p',[u1(j,i); params_Start; zeros(2,1)]); 
                X = x_out.xf; % memorize state for next iteration
                out(1:7,j,i+1) = full(X);
            end
        end
%         y1(:,:) = out(1,:,:);
        clear opts F X x_out i j
        
%         figure(11); plot(y1'); title('simulated output')
    
%% create nlp
    disp('Creating nlp')
    [ nlp, m, lbw, ubw, lbg, ubg, w0 ] = get_nlp(dae, u1, y1, out, n_cyc, N, tt);
    params_init=[w0; params_Start];

    
%% Create an NLP solver
    disp('Creating nlp solver')
    options.ipopt.max_iter=1000;                                             
    options.ipopt.linear_solver='ma57';
    %options.hess_lag = hess_lag;
tic
    solver = nlpsol('solver', 'ipopt', nlp, options);
toc

%% Solve the NLP
    disp('solving nlp')
tic
    sol = solver('x0', params_init, 'lbx', lbw, 'ubx', ubw,...
                'lbg', lbg, 'ubg', ubg, 'p',m);
toc
        sol.x = full(sol.x);

%% Plot the results
% Parameters estimated
    np = 8;%length(params_d); 
    params_opt = sol.x(end-np+1:end,:)
    % Wfg; Wff; Lambda; SG; Kd; Tau; Kf; Fs;

    nu = 1;
    nx = 7;
    nb = 2;
    s = nx+nu+nb;    

% Plot the Optimized COM
x1_opt = [];
for i= 1:n_cyc
    x1i_opt = sol.x( 1+(i-1)*( N*s+(s-1) ) :s: 1+(i-1)*( N*s+(s-1) ) + N*s);
    x1_opt = [x1_opt; x1i_opt];
end

x4_opt = [];
for i= 1:n_cyc
    x4i_opt = sol.x( 1+(i-1)*( N*s+(s-1) ) +3:s: 1+(i-1)*( N*s+(s-1) ) + N*s +3);
    x4_opt = [x4_opt; x4i_opt];
end

x5_opt = [];
for i= 1:n_cyc
    x5i_opt = sol.x( 1+(i-1)*( N*s+(s-1) ) +4 :s: 1+(i-1)*( N*s+(s-1) ) + N*s +4);
    x5_opt = [x5_opt; x5i_opt];
end

b1_opt = [];
for i= 1:n_cyc
    b1i_opt = sol.x( 1+(i-1)*( N*s+(s-1) ) + nx:s: 1+(i-1)*( N*s+(s-1) ) + N*s + nx);
    b1_opt = [b1_opt; b1i_opt];
end

b7_opt = [];
for i= 1:n_cyc
    b7i_opt = sol.x( 1+(i-1)*( N*s+(s-1) ) + nx +1:s: 1+(i-1)*( N*s+(s-1) ) + N*s + nx+1);
    b7_opt = [b7_opt; b7i_opt];
end

u_opt = [];
for i= 1:n_cyc
    ui_opt =  sol.x( 1+(i-1)*( N*s+(s-1) ) +nx+nb :s: 1+(i-1)*( N*s+(s-1) ) +nx+nb  + (N-1)*s);
    u_opt = [u_opt; ui_opt];
end

    t_sim = [0:1/ff:n_cyc*T+(n_cyc-1)*1/ff];
    y1 = reshape(y1', [(N+1)*n_cyc,1]);
    y1 = y1';
    u1 = u1(:, 1:end-1);
    u1 = reshape(u1', [N*n_cyc,1]);
    u1 = u1';

%%

figure(2)
hold on
plot(t_sim, y1)
plot(t_sim, x1_opt)
% plot(t_sim, X_int(1,:))
legend('Measured','Optimized', 'Location', 'Northeast') %,'Simulated'
title('COM')
xlabel('Time (s)')
ylabel('COM sway ({\circ})')

t_sim_u = [];
for i=0:n_cyc-1

t_sim_ui = [i*(T+1/ff) :1/ff: i*(T+1/ff)+T-1/ff];
t_sim_u = [t_sim_u, t_sim_ui];

end

figure(3)
hold on
plot(t_sim_u,u1)
plot(t_sim_u, u_opt')
legend('Applied','Optimized')
title('Tilt stimulus')
xlabel('Time (s)')
ylabel('Tilt angle ({\circ})')

% % calculate and plot ankle torque from function tc
% Tc=tc(x4_opt,x5_opt,params_opt(5),params_opt(6));
% Tc=full(Tc);
% 
% figure(4)
% hold on
% plot(t_sim, Tc)
% plot(t_sim, x1_opt)
% legend('Torque','COM')
% 

figure(5)
hold on
plot(t_sim, b1_opt)
plot(t_sim, b7_opt)
legend('b1', 'b7')
title('Process noise')
xlabel('Time (s)')
%axis([0 60.5 -0.02 0.02])

