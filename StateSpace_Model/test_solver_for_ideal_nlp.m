%% initialize
addpath(genpath('/home/haidang/Matlab/LA_toolbox'))
addpath('/home/haidang/Matlab/casadi-linux-matlabR2014b-v3.4.5')


clear all
  
import casadi.*

%% get input and output data
sr = 100; % sampling rate
N_cyc = 3;
N_amp = 3;

% Turn noise test on=1 or off else
noise_test = 1;
if noise_test==1
    noise_amp = 4; % Compute with different noises
    Q = cell(noise_amp,6);
     
elseif noise_test==0
    noise_amp = 1;
    Q = cell(N_amp,N_cyc,6); 
end

amp = 1:N_amp;%[1,2,4,8]
%[x,xi]=pseudogen3(2,[1 1],[0 2],25);
[x,xi]=pseudogen3(3,[1 0 -1],[1 0 2],25);
%[x,xi]=pseudogen3(5,[0 0 1 -1 -1],[2 0 1 0 2],25);
ndata = N_amp*size(x,2);
np = 8;
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
             % [K_{grav};K_{th};lambda; K_{s}; K_{d}; tau;K_{f};tau_f
                

%% initialize result cell array


if noise_test == 1
    for i_noise=1:noise_amp
        u = zeros(N_amp,size(xi,2)*N_cyc);

        % simulate N_cyc cycles
        for k=1:N_amp
            u(k,:) = repmat( amp(1,k)*( xi/(max(xi)-min(xi)) ),1,N_cyc);
        end

        u = u(:,1:NSamp:end,1);
        Q{i_noise,1} = u; % Save stimulus
        N = size(u,2)-1;
        
        % evaluate F for every time step

        out = zeros(7,size(u,1),size(u,2));
        for j=1:N_amp % number of different amplitudes
            X = zeros(length(ode.x),1); % set all states to zero (start of new amplitude)
            for i=1:N % samples within an amplitude
                x_out = F('x0',X,'p',[u(j,i); params_Start; zeros(2,1)]); 
                X = x_out.xf; % memorize state for next iteration
                out(1:7,j,i+1) = full(X);
            end
        end

        % get first state trajectory (which is body sway)
        y(:,:) = out(1,:,:);
        Q{i_noise,2} = y;   % Save body sway trajectory
        
        % add sway variability to simulated sway response
        noise = wgn(size(y,1),size(y,2),0); % white noise Signal
        fnyq = NSamp/2;
        [b,a] = butter(2,fco/fnyq,'bandpass'); % bandpass filter
        noise = filter(b,a,noise);
        
        y = y + 10*i_noise*noise;

        Q{i_noise,3} = y;   % Save body sway trajectory after including noise
        %% create nlp
        disp('Creating nlp')
        [ nlp, m, lbw, ubw, lbg, ubg, w0, R ] = get_nlp(ode, u, y, out, N_amp, N, tt);

        %% Create an NLP solver
        disp('Creating nlp solver')
        options.ipopt.max_iter=100;                                             
        %options.ipopt.linear_solver='ma57';

        tic
        solver = nlpsol('solver', 'ipopt', nlp, options);
        toc

        %% Solve the NLP
        disp('solving nlp')

        % params_init=[w0; params_Start]; % start with ideal parameters
        params_init=[w0;2; 0.2; 0.2; 4; 0.2*mp*9.81*H/180*pi; 0.7; 0.2; 33]; % vary starting parameters from ideal ones
                  % [w0;K_{grav};K_{th};lambda; K_{s}; K_{d}; tau;K_{f};tau_f
        tic
        sol = solver('x0', params_init, 'lbx', lbw, 'ubx', ubw,...
                    'lbg', lbg, 'ubg', ubg, 'p',m);
        toc
        sol.x = full(sol.x);
        Q{i_noise,4}=sol.x(end-7:end,:);    % Save solution of the parameter
        params_opt = sol.x(end-np+1:end,:);
        ng = size(nlp.g,1);
        nx = size(nlp.x,1);
        % Define derivative of the equality constraints and the residuals
        jac_g_fun = Function('jac_g_fun',{nlp.x, nlp.p},{jacobian(nlp.g,nlp.x)}); 
        jac_r_fun = Function('jac_r_fun',{nlp.x},{jacobian(R,nlp.x)});

        % Evaluate the derivatives at the solution 
        jac_g_opt = jac_g_fun(sol.x, m);
        jac_r_opt = jac_r_fun(sol.x);

        % Convert to full matrices
        Dr_full = full(jac_r_opt);
        Dg_full = full(jac_g_opt);

        %% Build the M matrix to extract the covariance matrix via LU decomposition
        Dr_fullmat = Dr_full'*Dr_full;
        M = sparse([jac_r_opt'*jac_r_opt jac_g_opt'; jac_g_opt zeros(ng)]); 
        M_decomposed = decomposition(M,'LU');
        XZ = M_decomposed \ [eye(nx); zeros(ng, nx)];
        X  = XZ(1:nx, 1:nx);
        C_p = X(end-np+1:end, end-np+1:end);

        % Rank check for more information see documentary
        r_g = rank(Dg_full);
        r_r = rank([Dr_full; Dg_full]);
        sprintf('Rank check with N_amp= %d and noise_amp= %d',N_amp,noise_amp)
        if r_g ~= ng || r_r ~= nx
            sprintf('rank check failed : Rank g = %d and ng = %d; Rank r = %d and nx = %d',r_g,ng,r_r,nx)
        else
            disp('rank check succesfull')
        end

        %% Calculate confidence interval
        %conf_b = [params_opt-(1.96.*(sqrt(diag(C_p)./(ndata)))) params_opt+(1.96.*(sqrt(diag(C_p)./(ndata))))];
        conf_b = [params_opt-(2.*np/(ndata-np)).*1.96.*(sqrt(diag(C_p))) params_opt+(2.*np/(ndata-np)).*1.96.*(sqrt(diag(C_p)))]; 
        conf_err = 1.96.*(sqrt(diag(C_p)./(ndata)));
        Q{i_noise,5} = conf_b;  % Save confident bounds
        Q{i_noise,6} = conf_err;

        clear y w0 conf_b conf_err

    end

% If noise test is turned off
elseif noise_test ==0    
    for i_amp=1:N_amp
        for j_cyc=1:N_cyc
            u = zeros(i_amp,size(xi,2)*j_cyc);

            % simulate N_cyc cycles
            for k=1:i_amp
                u(k,:) = repmat( amp(1,k).*( xi./(max(xi)-min(xi)) ),1,j_cyc);
            end

            u = u(:,1:NSamp:end,1);
            Q{i_amp,j_cyc,1} = u;   % Save stimulus
            N = size(u,2)-1;

            % evaluate F for every time step

            out = zeros(7,size(u,1),size(u,2));
            for j=1:i_amp % number of different amplitudes
                X = zeros(length(ode.x),1); % set all states to zero (start of new amplitude)
                for i=1:N % samples within an amplitude
                    x_out = F('x0',X,'p',[u(j,i); params_Start; zeros(2,1)]); 
                    X = x_out.xf; % memorize state for next iteration
                    out(1:7,j,i+1) = full(X);
                end
            end

            % get first state trajectory (which is body sway)
            y(:,:) = out(1,:,:);
            Q{i_amp,j_cyc,2} = y;   % Save body sway trajectory

            % add sway variability to simulated sway response
            noise = wgn(size(y,1),size(y,2),0); % white noise Signal
            fnyq = NSamp/2;
            [b,a] = butter(2,fco/fnyq,'bandpass'); % bandpass filter
            noise = filter(b,a,noise);

            y = y + noise_amp*noise;

            Q{i_amp,j_cyc,3} = y;   % Save body sway trajectory after including noise
            %% create nlp
            disp('Creating nlp')
            [ nlp, m, lbw, ubw, lbg, ubg, w0, R ] = get_nlp(ode, u, y, out, i_amp, N, tt);

            %% Create an NLP solver
            disp('Creating nlp solver')
            options.ipopt.max_iter=100;                                             
            %options.ipopt.linear_solver='ma57';

            tic
            solver = nlpsol('solver', 'ipopt', nlp, options);
            %Hess = solver.get_function('nlp_hess_l');
            toc

            %% Solve the NLP
            disp('solving nlp')

            % params_init=[w0; params_Start]; % start with ideal parameters
            params_init=[w0;2; 0.2; 0.2; 4; 0.2*mp*9.81*H/180*pi; 0.7; 0.2; 33]; % vary starting parameters from ideal ones
                            % [K_{grav};K_{th};lambda; K_{s}; K_{d}; tau;K_{f};tau_f
            tic
            sol = solver('x0', params_init, 'lbx', lbw, 'ubx', ubw,...
                        'lbg', lbg, 'ubg', ubg, 'p',m);
            toc
            sol.x = full(sol.x);
            Q{i_amp,j_cyc,4}=sol.x(end-7:end,:);    % Save solution of the parameter
            params_opt = sol.x(end-np+1:end,:);
            ng = size(nlp.g,1);
            nx = size(nlp.x,1);
            % Define derivative of the equality constraints and the residuals
            jac_g_fun = Function('jac_g_fun',{nlp.x, nlp.p},{jacobian(nlp.g,nlp.x)}); 
            jac_r_fun = Function('jac_r_fun',{nlp.x},{jacobian(R,nlp.x)});

            % Evaluate the derivatives at the solution 
            jac_g_opt = jac_g_fun(sol.x, m);
            jac_r_opt = jac_r_fun(sol.x);

            % Convert to full matrices
            Dr_full = full(jac_r_opt);
            Dg_full = full(jac_g_opt);

            %% Build the M matrix to extract the covariance matrix via LU decomposition
            Dr_fullmat = Dr_full'*Dr_full;
            M = sparse([jac_r_opt'*jac_r_opt jac_g_opt'; jac_g_opt zeros(ng)]); 
            M_decomposed = decomposition(M,'LU');
            XZ = M_decomposed \ [eye(nx); zeros(ng, nx)];
            X  = XZ(1:nx, 1:nx);
            C_p = X(end-np+1:end, end-np+1:end);

            % Rank check for more information see documentary
            r_g = rank(Dg_full);
            r_r = rank([Dr_full; Dg_full]);
            sprintf('Rank check with i_amp= %d and j_cyc= %d',i_amp,j_cyc)
            if r_g ~= ng || r_r ~= nx
                sprintf('rank check failed : Rank g = %d and ng = %d; Rank r = %d and nx = %d',r_g,ng,r_r,nx)
            else
                disp('rank check succesfull')
            end

            %% Calculate confidence interval
            %conf_b = [params_opt-(2.*np/(ndata-np)).*1.96.*(sqrt(diag(C_p))) params_opt+(2.*np/(ndata-np)).*1.96.*(sqrt(diag(C_p)))]; 
            conf_b = [params_opt-(1.96.*(sqrt(diag(C_p)./(ndata)))) params_opt+(1.96.*(sqrt(diag(C_p)./(ndata))))];
            conf_err = 1.96.*(sqrt(diag(C_p)./(ndata)));
            Q{i_amp,j_cyc,5} = conf_b; % Save confidence bounds
            Q{i_amp,j_cyc,6} = conf_err;
            clear y w0 conf_b conf_err
        end
    end
end




if noise_test == 1
    %% Plot results testing noise amplitude

    bound = reshape(cell2mat(Q(:,5)),np,noise_amp*2);
    bound_error = reshape(cell2mat(Q(:,6)),np,noise_amp);
    est_param = reshape(cell2mat(Q(:,4)),np,noise_amp);
    p_names = ["K_{grav}"; "K_{th}"; "\lambda"; "K_{s}"; "K_{d}"; "\tau"; "K_{f}"; "\tau_f"];
    figure(1)
    for i=1:8
        subplot(4,2,i);    
        plot_real = plot(1:noise_amp,params_Start(i).*ones(noise_amp,1),'s','linewidth',2);
        plot_real.LineStyle = 'none';
        plot_real.MarkerSize = 8;
        plot_real.MarkerEdgeColor = 'blue';
        plot_real.MarkerFaceColor = [0.2 .6 .6]; 
        hold on
        plot_CI1 = errorbar(1:noise_amp,est_param(i,:)',bound_error(i,:),'linewidth',1.5);
        plot_CI1.LineStyle = 'none';
        plot_CI1.Marker = 'o';
        plot_CI1.MarkerSize = 5;
        plot_CI1.MarkerEdgeColor = 'red';
        plot_CI1.MarkerFaceColor = [1 .6 .6];    
        set(gca, 'XLim',[0 noise_amp+1],'XTick', 1:noise_amp)    

        ylim= ([min(bound(i,:))-0.1*abs(min(bound(i,:))) max(bound(i,:))+0.2*abs(max(bound(i,:)))]);

        title(p_names(i) )
        xlabel('Amplitude number')
        grid
    end
elseif noise_test==0
    %% Plot results testing duration and amplitude 
    cyc_ind = 2;
    bound = reshape(cell2mat(Q(:,cyc_ind,5)),np,N_amp*2);
    bound_error = reshape(cell2mat(Q(:,cyc_ind,6)),np,N_amp);
	est_param = reshape(cell2mat(Q(:,cyc_ind,4)),np,N_amp);
    p_names = ["K_{grav}"; "K_{th}"; "\lambda"; "K_{s}"; "K_{d}"; "\tau"; "K_{f}"; "\tau_f"];
    figure(2)
    for i=1:8
        subplot(4,2,i);    
        plot_real = plot(1:N_amp,params_Start(i).*ones(N_amp,1),'s','linewidth',2);
        plot_real.LineStyle = 'none';
        plot_real.MarkerSize = 8;
        plot_real.MarkerEdgeColor = 'blue';
        plot_real.MarkerFaceColor = [0.2 .6 .6]; 
        hold on
        plot_CI1 = errorbar(1:N_amp,est_param(i,:)',bound_error(i,:),'linewidth',1.5);
        plot_CI1.LineStyle = 'none';
        plot_CI1.Marker = 'o';
        plot_CI1.MarkerSize = 5;
        plot_CI1.MarkerEdgeColor = 'red';
        plot_CI1.MarkerFaceColor = [1 .6 .6];    
        set(gca, 'XLim',[0 N_amp+1],'XTick', 1:N_amp)    
        ylim= ([min(bound(i,:))-0.1*abs(min(bound(i,:))) max(bound(i,:))+0.2*abs(max(bound(i,:)))]);

        title(p_names(i) )
        xlabel('Amplitude number')
        grid
    end
end
%%
% close all
% figure(2)
% Cycle index
% cyc_ind = 2;
% Parameter index
% param_ind = 1;
% 
% real_param = params_Start(param_ind).*ones(N_amp,1);
% bound = reshape(cell2mat(Q(:,cyc_ind,5)),np,N_amp*2);
% est_param = reshape(cell2mat(Q(:,cyc_ind,4)),np,N_amp);
% plot_CI = errorbar(1:N_amp,est_param(param_ind,:)',bound(param_ind,1:N_amp)',bound(param_ind,N_amp+1:end)');
% 
% plot_CI.LineStyle = 'none';
% xticks([ 1 2 ])
% xticklabels({'1','2'})
% axis([0 N_amp+1 min(bound(param_ind,:))-1 max(bound(param_ind,:))+1])
% plot_CI.Marker = 'o';
% plot_CI.MarkerSize = 5;
% plot_CI.MarkerEdgeColor = 'red';
% plot_CI.MarkerFaceColor = [1 .6 .6];
% title(['Confidence bounds for parameter ' ,num2str(param_ind), ' with ', num2str(cyc_ind),' cycles' ])
% xlabel('Amplitude number')
% grid on

%%
% Solution with k cycle and l amplitude
% k=1;
% l=1;
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



% nu = 1;
% nx = 7;
% nb = 2;
% s = nx+nu+nb;    
% 
% for i = N_amp:-1:1 
%     u_in(i*(N+1)) = [];
% end
% %% Get Jacobian       
% y_in = reshape(Q{2,1,3}',[(N+1)*N_amp,1]);
% u_in = reshape(Q{2,1,4}',[(N+1)*N_amp,1]);
% res_u = abs(u_opt - u_in);
% res_y = abs(y_opt - y_in);
% res_b1 = abs(b1_opt);
% res_b7 = abs(b7_opt);

% y_opt = [];
% for i= 1:N_amp*N_cyc
%     y_opt = [y_opt; sol.x( 1+(i-1)*( N*s+(s-1) ) :s: 1+(i-1)*( N*s+(s-1) ) + N*s)];
% end
% u_opt = [];
% for i= 1:N_amp*N_cyc
%     u_opt = [u_opt; sol.x( 1+(i-1)*( N*s+(s-1) ) +nx+nb :s: 1+(i-1)*( N*s+(s-1) ) +nx+nb  + (N-1)*s)];
% end
% b1_opt = [];
% for i= 1:N_amp*N_cyc
%     b1_opt = [b1_opt; sol.x( 1+(i-1)*( N*s+(s-1) ) + nx:s: 1+(i-1)*( N*s+(s-1) ) + N*s + nx)];
% end
% b7_opt = [];
% for i= 1:N_amp*N_cyc
%     b7_opt = [b7_opt;sol.x( 1+(i-1)*( N*s+(s-1) ) + nx +1:s: 1+(i-1)*( N*s+(s-1) ) + N*s + nx+1)];
% end

