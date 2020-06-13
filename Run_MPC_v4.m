%% Load model parameters
load('SEIQRDP_ZA_v5.mat');

%%
dt = 1; % time step
Tadd = 500;
time1 = datetime(t_start):dt:datetime(datestr(floor(now)+datenum(Tadd)));
N = numel(time1);
t = (0:N-1).*dt;

%Beta2 = [0.35 0.34 0.0025];
Beta2 = [0.25 0.364 0.00265];

[S,E,I,Q,R,D,P] = ...
     SEIQRDP_v3(alpha,Beta2,gamma,delta,Lambda,Kappa,Npop,E0,I0,Q0,R0,D0,t,lambdaFun);
 
 plot(t,Q)

%% Define model
% Evaluate time-based parameters to get constant values - can expand to
% time-variant model later
kMax = 1000;
t_beta = 0:(2*kMax);

t_eval = 1;
beta = betaFun(Beta2,t_beta);
lambda = lambdaFun(Lambda,t_eval);
kappa = kappaFun(Kappa,t_eval);
A = [-alpha   0      0             0         0 0 0;
    0    -gamma   0             0         0 0 0;
    0     gamma -delta          0         0 0 0;
    0      0     delta  (-kappa - lambda) 0 0 0;
    0      0      0          lambda       0 0 0;
    0      0      0          kappa        0 0 0;
    alpha  0      0             0         0 0 0];

C = eye(7);

%% Run simulation with no control
dT = 1;

% States are:
% S - Susceptible population
% E - Exposed
% I - Infectious
% Q - Quarantined (confirmed infectious)
% R - Recovered
% D - Deceased
% P - Insusceptible
X0 = [Npop-Q0-E0-R0-D0-I0, E0, I0, Q0, R0, D0, 0];

kMax = 1000;
t = 0:kMax;

% Run sim
Y = zeros(kMax+1,7);
Y(1,:) = X0;
X = X0';

for k = 2:(kMax+1)
    
    B = [-beta(k)/Npop beta(k)/Npop 0 0 0 0 0]';
    
    U = Y(k-1,1)*Y(k-1,3);
    
    % Runge-Kutta method
    k_1 = A*X + B*U;
    k_2 = A*(X+0.5*dT*k_1) + B*U;
    k_3 = A*(X+0.5*dT*k_2) + B*U;
    k_4 = A*(X+dT*k_3) + B*U;
    % States
    X = X + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dT;
    
    Y(k,:) = C*X;
end

figure(1)
plot(t,Y(:,2:7))
legend_entries = {'S','E','I','Q','R','D','P'};
legend(legend_entries(2:7));

%% Run control simulation where input is lockdown level (integer value)

% Beta multipliers
% Level 5 - 1
% Level 4 - 1.25
% Level 3 - 1.75
% Level 2 - 2.25
% Level 1 - 2.5
%beta_levels = [2.5 2.0 1.6 1.25 1];
beta_levels = [1/0.6 1/0.6*0.95 1.5 1.25 1];

%Initialize control parameters
dT = 1;
kMax = 500;
tsim = 0:kMax;

% Control quarantined cases Y(4)
y_hi = 1e7;
y_lo = 0;

% Initialze simulation in/outputs
Yc = zeros(kMax+1,7);
Yc(1,:) = X0;
Xc = X0';
Umpc = zeros(kMax+1,1);
Umpc(1) = 1;

% Control config variables
Np = 300;
Nc = 30; % When using blocking the control horizon is not used
Nb = [10 10 10]; % Blocking
S = 1e-6;
R = 1e3;
T = 1e2;

%options = optimoptions('ga','Display', 'off', 'FitnessLimit', 0, 'MaxGenerations', 100, 'MaxTime', 5, 'FunctionTolerance', 1);

% Run control
% Setup waitbar that shows simulation progress
h = waitbar(0,'Simulating...');

for k = 2:(kMax+1)
    tic;
    
    B = [-beta(k)/Npop beta(k)/Npop 0 0 0 0 0]';
    % Uncontrolled portion of input
    Un = Yc(k-1,1)*Yc(k-1,3);
    
    % Run controller on first iteration and every 10th iteration after that
    if (mod(k,10) - 2) == 0
        
        % On the first run, start without an InitialPopulationMatrix, on each
        % successive run initialize with the final population from the previous
        if k == 2
            [Uc,fval,exitflag,output,final_pop] = ga(@(Ucon) FitnessFunc_v2(Ucon,Np,Nc,Nb,A,B,C,y_lo,y_hi,S,R,T,beta_levels,k,beta,Npop,dT,Xc),3,[],[],[],[],[1 1 1],[5 5 5],[],[1 2 3]);
        else
            [Uc,fval,exitflag,output,final_pop] = ga(@(Ucon) FitnessFunc_v2(Ucon,Np,Nc,Nb,A,B,C,y_lo,y_hi,S,R,T,beta_levels,k,beta,Npop,dT,Xc),3,[],[],[],[],[1 1 1],[5 5 5],[],[1 2 3]);
        end
        %options = optimoptions('ga','InitialPopulationMatrix', final_pop, 'Display', 'off', 'FitnessLimit', 0, 'MaxGenerations', 50, 'MaxTime', 0.5, 'FunctionTolerance', 1);
        %options = optimoptions('ga','InitialPopulationMatrix', final_pop, 'Display', 'off', 'FitnessLimit', 0, 'MaxGenerations', 100, 'MaxTime', 5, 'FunctionTolerance', 1);
        
        Uinp = Uc(1);
        
    end
    
    Umpc(k) = Uinp;
    
    U = Un*beta_levels(Uinp);
    % Runge-Kutta method
    k_1 = A*Xc + B*U;
    k_2 = A*(Xc+0.5*dT*k_1) + B*U;
    k_3 = A*(Xc+0.5*dT*k_2) + B*U;
    k_4 = A*(Xc+dT*k_3) + B*U;
    % States
    Xc = Xc + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dT;
    
    Yc(k,:) = C*Xc;
    tt = toc;
    % Display time per step or comment out
    fprintf('Run %d took %d seconds \n',k,tt)
    waitbar(k/kMax)
end

close(h)

% Plot results of control simulation

figure(2)
subplot(2,1,1)
plot(tsim,Yc(:,2:7))
legend('E','I','Q','R','D','P')
subplot(2,1,2)
plot(tsim,Umpc)
legend('U')

%% Run simulation with own control

kMax = 1000;
t = 0:kMax;

% Run sim
Yo = zeros(kMax+1,7);
Yo(1,:) = X0;
Xo = X0';

%U_pre_def = ones(kMax+1,1);
%U_pre_def = beta_levels(Umpc);
U_pre_def = 2.5.*ones(kMax+1,1);
U_pre_def(30:300) = 1;
U_pre_def(301:end) = 1.9;

for k = 2:(kMax+1)
    
    B = [-beta(k)/Npop beta(k)/Npop 0 0 0 0 0]';
    
    Uo = Yo(k-1,1)*Yo(k-1,3)*U_pre_def(k);
    
    % Runge-Kutta method
    k_1 = A*Xo + B*Uo;
    k_2 = A*(Xo+0.5*dT*k_1) + B*Uo;
    k_3 = A*(Xo+0.5*dT*k_2) + B*Uo;
    k_4 = A*(Xo+dT*k_3) + B*Uo;
    % States
    Xo = Xo + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dT;
    
    Yo(k,:) = C*Xo;
end

figure(3)
subplot(2,1,1)
plot(t,Yo(:,2:7))
legend_entries = {'S','E','I','Q','R','D','P'};
legend(legend_entries(2:7));
subplot(2,1,2)
plot(t,U_pre_def)
legend('U')

%% Run control simulation where input level is a real value

%Initialize control parameters
dT = 1;
kMax = 500;
tsim = 0:kMax;

% Control quarantined cases Y(4)
y_hi = 1e7;
y_lo = 0;

% Initialze simulation in/outputs
Yc_o = zeros(kMax+1,7);
Yc_o(1,:) = X0;
Xc = X0';
Umpc_o = zeros(kMax+1,1);
Umpc_o(1) = 1;

% Control config variables
Np = 300;
Nc = 30; % When using blocking the control horizon is not used
Nb = [10 10 10]; % Blocking
S = 1e-6;
R = 1e3;
T = 1e1;

% Run control
% Setup waitbar that shows simulation progress
h = waitbar(0,'Simulating...');

for k = 2:(kMax+1)
    tic;
    
    B = [-beta(k)/Npop beta(k)/Npop 0 0 0 0 0]';
    % Uncontrolled portion of input
    Un = Yc_o(k-1,1)*Yc_o(k-1,3)*2.5;
    
    % Run controller on first iteration and every 10th iteration after that
    if (mod(k,10) - 2) == 0
        
        % On the first run, start without an InitialPopulationMatrix, on each
        % successive run initialize with the final population from the previous
        if k == 2
            [Uc,fval,exitflag,output,final_pop] = fmincon(@(Ucon) FitnessFunc_v3(Ucon,Np,Nc,Nb,A,B,C,y_lo,y_hi,S,R,T,k,beta,Npop,dT,Xc),[1 1 1],[],[],[],[],[0.4 0.4 0.4],[1 1 1]);
        else
            [Uc,fval,exitflag,output,final_pop] = fmincon(@(Ucon) FitnessFunc_v3(Ucon,Np,Nc,Nb,A,B,C,y_lo,y_hi,S,R,T,k,beta,Npop,dT,Xc),[1 1 1],[],[],[],[],[0.4 0.4 0.4],[1 1 1]);
        end
        
        Uinp = Uc(1);
        
    end
    
    Umpc_o(k) = Uinp;
    
    U = Un*Uinp;
    % Runge-Kutta method
    k_1 = A*Xc + B*U;
    k_2 = A*(Xc+0.5*dT*k_1) + B*U;
    k_3 = A*(Xc+0.5*dT*k_2) + B*U;
    k_4 = A*(Xc+dT*k_3) + B*U;
    % States
    Xc = Xc + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dT;
    
    Yc_o(k,:) = C*Xc;
    tt = toc;
    % Display time per step or comment out
    fprintf('Run %d took %d seconds \n',k,tt)
    waitbar(k/kMax)
end

close(h)

% Plot results of control simulation

figure(4)
subplot(2,1,1)
plot(tsim,Yc_o(:,2:7))
legend('E','I','Q','R','D','P')
subplot(2,1,2)
plot(tsim,(1-Umpc_o)*4/0.6+1)
legend('U')