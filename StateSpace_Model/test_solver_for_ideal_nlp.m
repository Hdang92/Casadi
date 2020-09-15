%% initialize
addpath(genpath('/home/haidang/Matlab/LA_toolbox'))
addpath('/home/haidang/Matlab/casadi-linux-matlabR2014b-v3.4.5')


clear all
  
import casadi.*

%% get input and output data
sr = 100; % sampling rate
n_cyc = 2;
n_stimuli = 2;

% Turn noise test on=1 or off else
noise_test = 1;
if noise_test==1
    noise_amp = 4; % Compute with different noises
    Q = cell(noise_amp,6);
     
elseif noise_test==0
    noise_amp = 1;
    Q = cell(n_stimuli,n_cyc,6); 
end

amp = 1:n_stimuli;%[1,2,4,8]
%[x,xi]=pseudogen3(2,[1 1],[0 2],25);
[x,xi]=pseudogen3(3,[1 0 -1],[1 0 2],25);
%[x,xi]=pseudogen3(5,[0 0 1 -1 -1],[2 0 1 0 2],25);
ndata = n_stimuli*size(x,2);
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
    for n_noise=1:noise_amp
        u = zeros(n_stimuli,size(xi,2)*n_cyc);

        % simulate N_cyc cycles
        for k=1:n_stimuli
            u(k,:) = repmat( amp(1,k)*( xi/(max(xi)-min(xi)) ),1,n_cyc);
        end

        u = u(:,1:NSamp:end,1);
        Q{n_noise,1} = u; % Save stimulus
        N = size(u,2)-1;
        
        % evaluate F for every time step

        out = zeros(7,size(u,1),size(u,2));
        for j=1:n_stimuli % number of different amplitudes
            X = zeros(length(ode.x),1); % set all states to zero (start of new amplitude)
            for i=1:N % samples within an amplitude
                x_out = F('x0',X,'p',[u(j,i); params_Start; zeros(2,1)]); 
                X = x_out.xf; % memorize state for next iteration
                out(1:7,j,i+1) = full(X);
            end
        end

        % get first state trajectory (which is body sway)
        y(:,:) = out(1,:,:);
        Q{n_noise,2} = y;   % Save body sway trajectory
        
        % add sway variability to simulated sway response
        noise = wgn(1,size(y,2),0); % white noise Signal
        fnyq = NSamp/2;
        [b,a] = butter(2,fco/fnyq,'bandpass'); % bandpass filter
        noise = filter(b,a,noise);
        %noise = repmat(noise,n_stimuli,1);
        y = y + n_noise*noise;
        Q{n_noise,7} = noise(1,:);
        Q{n_noise,3} = y;   % Save body sway trajectory after including noise
        %% create nlp
        disp('Creating nlp')
        [ nlp, m, lbw, ubw, lbg, ubg, w0, R ] = get_nlp(ode, u, y, out, n_stimuli, N, tt);

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
        Q{n_noise,4}=sol.x(end-7:end,:);    % Save solution of the parameter
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
        sprintf('Rank check with N_amp= %d and noise_amp= %d',n_stimuli,noise_amp)
        if r_g ~= ng || r_r ~= nx
            sprintf('rank check failed : Rank g = %d and ng = %d; Rank r = %d and nx = %d',r_g,ng,r_r,nx)
        else
            disp('rank check succesfull')
        end

        %% Calculate confidence interval
        %conf_b = [params_opt-(1.96.*(sqrt(diag(C_p)./(ndata)))) params_opt+(1.96.*(sqrt(diag(C_p)./(ndata))))];
        conf_b = [params_opt-(2.*np/(ndata-np)).*1.96.*(sqrt(diag(C_p))) params_opt+(2.*np/(ndata-np)).*1.96.*(sqrt(diag(C_p)))]; 
        conf_err = 1.96.*(sqrt(diag(C_p)./(ndata)));
        Q{n_noise,5} = conf_b;  % Save confident bounds
        Q{n_noise,6} = conf_err;

        clear y w0 conf_b conf_err

    end

% If noise test is turned off
elseif noise_test ==0    
    for i_amp=1:n_stimuli
        for j_cyc=1:n_cyc
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
            %noise = repmat(noise,i_amp,1);
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
            Q{i_amp,j_cyc,5} = conf_b; % Save confident bounds
            Q{i_amp,j_cyc,6} = conf_err;
            clear y w0 conf_b conf_err
        end
    end
end


%% plot stimuli strength


    
%% Plot results testing noise amplitude
if noise_test == 1

    bound = reshape(cell2mat(Q(:,5)),np,noise_amp*2);
    bound_error = reshape(cell2mat(Q(:,6)),np,noise_amp);
    est_param = reshape(cell2mat(Q(:,4)),np,noise_amp);
         

    p_names = ["K_{grav}"; "K_{th}"; "\lambda"; "K_{s}"; "K_{d}"; "\tau"; "K_{f}"; "\tau_f"];
    figure(1) 
    set(gcf, 'Position',  [200, 200, 600, 800])  
    for i=1:np
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

        %ylim= ([min(bound(i,:))-0.1*abs(min(bound(i,:))) max(bound(i,:))+0.2*abs(max(bound(i,:)))]);

        title(p_names(i) )
        xlabel('number of noise (n_{noise})')
        grid
        matlab2tikz('5.tex')
    end
    %%
    noi = cell2mat(Q(1,7));
    figure(2)
    for i=1:n_noise
        
        plot(i.*noi(1,:));
        hold on
        xlabel('time (s)')
        ylabel('tilt angle (°)')
    end
        legend('pp 1° (n_{noise}=1)','pp 2° (n_{noise}=1)','pp 3° (n_{noise}=3)')%,'pp 4° (n_{noise}=4)'
%% Plot results testing duration and amplitude

elseif noise_test==0
    % Choose which cycle should be plotted
    cyc_ind = 4;
    
    bound_error = reshape(cell2mat(Q(:,cyc_ind,6)),np,n_stimuli);
    bound_error_old = bound_error;	
    est_param = reshape(cell2mat(Q(:,cyc_ind,4)),np,n_stimuli);
    % Scale plots by cutting errorbars which are 5 times larger then second
    % element
    index_cut = zeros(np,1);
    for j= 1:np
        if bound_error(j,1) > 5*bound_error(j,2)
            index_cut(j) = 1;
            %bound_error(j,1) = 2*bound_error(j,2);
        end
    end
    
    p_names = ["K_{grav}"; "K_{th}"; "\lambda"; "K_{s}"; "K_{d}"; "\tau"; "K_{f}"; "\tau_f"];
    figure(1)        
    set(gcf, 'Position',  [200, 200, 600, 800])    
    for i=1:np
        subplot(4,2,i);    
        plot_real = plot(1:n_stimuli,params_Start(i).*ones(n_stimuli,1),'s','linewidth',2);
        plot_real.LineStyle = 'none';
        plot_real.MarkerSize = 8;
        plot_real.MarkerEdgeColor = 'blue';
        plot_real.MarkerFaceColor = [0.2 .6 .6]; 
        hold on
        plot_CI1 = errorbar(1:n_stimuli,est_param(i,:)',bound_error(i,:),'linewidth',1.5);
        plot_CI1.LineStyle = 'none';
        plot_CI1.Marker = 'o';
        plot_CI1.MarkerSize = 5;
        plot_CI1.MarkerEdgeColor = 'red';
        plot_CI1.MarkerFaceColor = [1 .6 .6]; 
    
        if index_cut(i) ==1
            set(gca, 'XLim',[0 n_stimuli+1],'XTick', 1:n_stimuli, 'YLim',[est_param(i,2)-1.3*bound_error(i,2) bound_error(i,2)+1.3*est_param(i,2)])

        else
            set(gca, 'XLim',[0 n_stimuli+1],'XTick', 1:n_stimuli, 'YLim',[min(est_param(i,:)-bound_error(i,:))-0.1*abs(min(est_param(i,:)-bound_error(i,:))) max(est_param(i,:)+bound_error(i,:))+0.2*abs(min(est_param(i,:)+bound_error(i,:)))])
        end
        title(p_names(i) )
        xlabel('number of stimuli (n_{stimuli})')
        grid
        %matlab2tikz('4.tex')
    end
    %%
    stim = Q{n_stimuli,n_cyc,1};
    figure(2)
    for i = 1:n_stimuli
        
        plot(stim(i,:));
        hold on
        xlabel('time (s)')
        ylabel('tilt angle (°)')
        legend('pp 1°','pp 2°','pp 3°','pp 4°')
    end
end



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
