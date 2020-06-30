%% initialize
addpath(genpath('/home/haidang/Matlab/LA_toolbox'))
addpath('/home/haidang/Matlab/casadi-linux-matlabR2014b-v3.4.5')


clear all
  
import casadi.*

%% get input and output data
sr = 100; % sampling rate
N_cyc = 1;
N_amp = 2;
amp = 1:N_amp;%[1,2,4,8]
%[x,xi]=pseudogen3(2,[1 1],[0 2],25);
[x,xi]=pseudogen3(3,[1 0 -1],[1 0 2],25);
%[x,xi]=pseudogen3(5,[0 0 1 -1 -1],[2 0 1 0 2],25);


%% reduction of data points

NSamp = 25; % downsampling factor
% define data structure and time intervals
tt = NSamp/sr; 
ff=1/tt;
%T = 60.5-tt;
T = 6.5-tt;

% define noise for realistic data
fco = [0.02, 0.2]; % [Hz] Cut off Frequenzen fuer den bandpass filter

%% define body parameters
mp = 69;
H = 0.93;
J = (72.8)*pi/180;

%% get model
[ ode, tc ] = DEC_model(tt,mp,H,J);

%% Simulate model with fixed parameters
% Implementation of system dynamics in Casadi
opts = struct('tf', tt,'number_of_finite_elements',15, 'jit', false);
F = integrator('F', 'rk', ode, opts);
params_Start = [0.55; 0.25; 0.3; 1; 0.3*mp*9.81*H/180*pi; 0.08; 0.1; 20];
                % Wfg; Wff; Lambda; SG; Kd; Tau; Kf; Fs;
                

%% initialize result cell array
Q = cell(N_amp,N_cyc);


for j_cyc=1:N_cyc
    for i_amp=1:N_amp
        u = zeros(i_amp,size(xi,2)*j_cyc);

        % simulate N_cyc cycles
        for k=1:i_amp
            u(k,:) = repmat( amp(1,k)*( xi/(max(xi)-min(xi)) ),1,j_cyc);
        end

        u = u(:,1:NSamp:end,1);
        Q{i_amp,j_cyc,4} = u;
        N = size(u,2)-1;
        
        % evaluate F for every time step

        out = zeros(7,size(u,1),size(u,2));
        for j=1:i_amp % number of different amplitudes
            X = zeros(length(ode.x),1); % set all states to zero (start of new amplitude)
            for i=2:N % samples within an amplitude
                x_out = F('x0',X,'p',[u(j,i); params_Start; zeros(2,1)]); 
                X = x_out.xf; % memorize state for next iteration
                out(1:7,j,i+1) = full(X);
            end
        end

        % get first state trajectory (which is body sway)
        y(:,:) = out(1,:,:);
        
        Q{i_amp,j_cyc,2} = y;

        %plot(y'); title('simulated output')
        
        % add sway variability to simulated sway response
        noise = wgn(size(y,1),size(y,2),0); % white noise Signal
        fnyq = NSamp/2;
        [b,a] = butter(2,fco/fnyq,'bandpass'); % bandpass filter
        noise = filter(b,a,noise);
        
        y = y + noise;
        
        hold on
        %plot(y');
        Q{i_amp,j_cyc,3} = y;
        %% create nlp
        disp('Creating nlp')
        [ nlp, m, lbw, ubw, lbg, ubg, w0 ] = get_nlp(ode, u, y, out, i_amp, N, tt);

        %% Create an NLP solver
        disp('Creating nlp solver')
        options.ipopt.max_iter=5;                                             
        %options.ipopt.linear_solver='ma57';
        %options.hess_lag = hess_lag;

        tic
        solver = nlpsol('solver', 'ipopt', nlp, options);
        %Hess = solver.get_function('nlp_hess_l');
        toc

        %% Solve the NLP
        disp('solving nlp')

        % params_init=[w0; params_Start]; % start with ideal parameters
        params_init=[w0;0.55; 1.25; 0.8; 2; 0.3*mp*9.81*H/180*pi; 0.2; 0.5; 44]; % vary starting parameters from ideal ones
                        % Wfg; Wff; Lambda; SG; Kd; Tau; Kf; Fs;
        tic
        sol = solver('x0', params_init, 'lbx', lbw, 'ubx', ubw,...
                    'lbg', lbg, 'ubg', ubg, 'p',m);
        toc
        sol.x = full(sol.x);
        Q{i_amp,j_cyc,1}=sol.x(end-7:end,:);

        clear y w0
    end
end
%% Get Jacobian
    nu = 1;
    nx = 7;
    nb = 2;
    s = nx+nu+nb;          
y_in = reshape(Q{2,1,3}',[(N+1)*N_amp,1]);
u_in = reshape(Q{2,1,4}',[(N+1)*N_amp,1]);

for i = N_amp:-1:1 
    u_in(i*(N+1)) = [];
end

y_opt = [];
for i= 1:N_amp*N_cyc
    y_opt = [y_opt; sol.x( 1+(i-1)*( N*s+(s-1) ) :s: 1+(i-1)*( N*s+(s-1) ) + N*s)];
end
u_opt = [];
for i= 1:N_amp*N_cyc
    u_opt = [u_opt; sol.x( 1+(i-1)*( N*s+(s-1) ) +nx+nb :s: 1+(i-1)*( N*s+(s-1) ) +nx+nb  + (N-1)*s)];
end
b1_opt = [];
for i= 1:N_amp*N_cyc
    b1_opt = [b1_opt; sol.x( 1+(i-1)*( N*s+(s-1) ) + nx:s: 1+(i-1)*( N*s+(s-1) ) + N*s + nx)];
end
b7_opt = [];
for i= 1:N_amp*N_cyc
    b7_opt = [b7_opt;sol.x( 1+(i-1)*( N*s+(s-1) ) + nx +1:s: 1+(i-1)*( N*s+(s-1) ) + N*s + nx+1)];
end
res_u = abs(u_opt - u_in);
res_y = abs(y_opt - y_in);
res_b1 = abs(b1_opt);
res_b7 = abs(b7_opt);
wy = 10;
wu = 10;
wb = 0.1;
r = [sqrt(wu).*res_u; sqrt(wy).*res_y; sqrt(wb).*res_b1;sqrt(wb).*res_b7];
jac_f_fun = Function('jac_f_fun',{nlp.x, nlp.p},{jacobian(nlp.f,nlp.x)}); 
jac_g_fun = Function('jac_g_fun',{nlp.x, nlp.p},{jacobian(nlp.g,nlp.x)}); 
jac_r_fun = Function('jac_r_opt',{nlp.x},{jacobian(r,nlp.x)});
jac_f_opt = jac_f_fun(sol.x, m);
jac_g_opt = jac_g_fun(sol.x, m);
jac_r_opt = jac_r_fun(sol.x);
Df_full = full(jac_r_opt);
Dg_full = full(jac_g_opt);
Df_fullmat = Df_full'*Df_full;
det(Df_fullmat)
ng = size(nlp.g,1);
nx = size(nlp.x,1);
if rank(Dg_full) ~= ng %|| rank([Df_full; Dg_full]) ~= nx
    disp('rank check failed')
else
    disp('rank check succesfull')
end
 

%% Plot the results
% Solution with k cycle and l amplitude
k=1;
l=1;
% for j=2:k
%     
% figure(j)
% plot(Q{l,j,2})
% hold on
% plot(Q{l,j,3})
% legend('Simulated body sway','Body sway with noises', 'Location', 'Northeast') %,'Simulated'
% title('Body sway trajectory')
% xlabel('Time (s)')
% ylabel('Body sway ({\circ})')
% end
% 
% 
% P=Q(:,:,1);
% P=cell2mat(P);
% P=reshape(P,8,N_amp*N_cyc);
% % Calculate error
% Error = zeros(N_amp*N_cyc,1);
% for i=1:N_amp*N_cyc
% Error(i) = norm(P(:,i)-params_Start,2);
% end
% plot(Error(2:4))
% Parameters estimated
np = 8;%length(params_d); 
params_opt = sol.x(end-np+1:end,:);
% Wfg; Wff; Lambda; SG; Kd; Tau; Kf; Fs;

% nu = 1;
% nx = 7;
% nb = 2;
% s = nx+nu+nb;    

