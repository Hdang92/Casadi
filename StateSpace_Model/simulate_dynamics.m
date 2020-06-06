%% initialize
addpath(genpath('/home/hdang/Matlab/LA_toolbox'))
addpath('/home/hdang/Matlab/casadi-linux-matlabR2014b-v3.4.5')

close all
clc


import casadi.*

%% get input and output data
    load('cycles_eyes_closed.mat')
    n=0;
    for k=1:4
        for m=2
            n=n+1;
            u_t = mean(data(k).FS(m,:),1);
            u1(n,:) = u_t - mean(u_t);
            y_t = mean(data(k).BS(m,:),1);
            y1(n,:) = y_t - mean(y_t);
        end
    end
    clear k m n u_t y_t
    sr = 100; % sampling rate

%% reduction of data points
    NSamp = 5; % downsampling factor

    u_t=u1; y_t=y1; clear u1 y1
    for k=1:size(u_t,1);
        u1(k,:) = u_t(k,1:NSamp:end);
        y1(k,:) = y_t(k,1:NSamp:end);
    end
    clear u_t y_t k
    
    % define data structure and time intervals
        tt = NSamp/sr; ff=1/tt;
        N = size(u1,2)-1;
        T = 60.5-tt;
        n_cyc = size(u1,1);
        
        
%% define body parameters
        mp = 69;
        g0 = 9.81;
        H = 0.93;
        J = (72.8)*pi/180;
        
%% get model
    [ dae, tc ] = DEC_model(tt,mp,H,J);

%% Simulate model with fixed parameters
    % Implementation of system dynamics in Casadi
        opts = struct('tf', tt,'number_of_finite_elements',15, 'jit', true);
        F = integrator('F', 'rk', dae, opts);
        
        params_Start = [0.55; 0.25; 0.3; 1; 0.3*mp*g0*H/180*pi; 0.08; 0.1; 20];
        
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
        
%%

% % check the plot
y1_exp = reshape(y1', 1,size(y1,1)*size(y1,2));
t_sim = (1:length(y1_exp))*tt;
out_t(:,:) = out(1,:,:);
y1_sim(:,:) = reshape(out_t',1,size(out,2)*size(out,3));
% % check the plot

figure(1)
hold on
plot(t_sim,y1_exp)
plot(t_sim,y1_sim)
legend('Measured','Simulated')
