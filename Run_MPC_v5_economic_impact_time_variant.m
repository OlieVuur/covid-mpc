%% Load model parameters
load('SEIQRDP_ZA_v5.mat');

%% Define model
% Evaluate time-based parameters to get constant values - can expand to
% time-variant model later
kMax = 1000;
t_beta = 0:(2*kMax);

beta = betaFun(Beta2,t_beta);
lambda = lambdaFun(Lambda,t_beta);
kappa = kappaFun(Kappa,t_beta);

C = eye(7);

%% Run control simulation

% Beta multipliers
% Level 5 - 1
% Level 4 - 1.25
% Level 3 - 1.75
% Level 2 - 2.25
% Level 1 - 2.5
%beta_levels = [2.5 2.0 1.6 1.25 1];
beta_levels = [1/0.6 1/0.6*0.95 1.5 1.25 1];

economic_loss_by_level = [0 1 2 3 4];
cumulative_economic_loss = 0;

%Initialize control parameters
dT = 1;
X0 = [Npop-Q0-E0-R0-D0-I0, E0, I0, Q0, R0, D0, 0];
kMax = 1000;
tsim = 0:kMax;

% Control quarantined cases Y(4)
y_hi = 2.1e5;
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
R = 1e6;
T = 1e6;
E = 1e5;

%options = optimoptions('ga','Display', 'off', 'FitnessLimit', 0, 'MaxGenerations', 100, 'MaxTime', 5, 'FunctionTolerance', 1);

% Run control
% Setup waitbar that shows simulation progress
h = waitbar(0,'Simulating...');

for k = 2:(kMax+1)
    tic;
    
    A = [-alpha   0      0             0         0 0 0;
           0    -gamma   0             0         0 0 0;
           0     gamma -delta          0         0 0 0;
           0     0     delta  (-kappa(k) - lambda(k)) 0 0 0;
           0      0      0          lambda(k)       0 0 0;
           0      0      0          kappa(k)        0 0 0;
           alpha  0      0             0            0 0 0];
    
    B = [-beta(k)/Npop beta(k)/Npop 0 0 0 0 0]';
    % Uncontrolled portion of input
    Un = Yc(k-1,1)*Yc(k-1,3);
    
    % Run controller on first iteration and every 10th iteration after that
    if (mod(k,10) - 2) == 0
        
        % On the first run, start without an InitialPopulationMatrix, on each
        % successive run initialize with the final population from the previous
        if k == 2
            [Uc,fval,exitflag,output,final_pop] = ga(@(Ucon) FitnessFunc_v4_time_variant(Ucon,Np,Nc,Nb,alpha,gamma,delta,lambda,kappa,B,C,E,economic_loss_by_level,cumulative_economic_loss,y_lo,y_hi,S,R,T,beta_levels,k,beta,Npop,dT,Xc),3,[],[],[],[],[1 1 1],[5 5 5],[],[1 2 3]);
        else
            [Uc,fval,exitflag,output,final_pop] = ga(@(Ucon) FitnessFunc_v4_time_variant(Ucon,Np,Nc,Nb,alpha,gamma,delta,lambda,kappa,B,C,E,economic_loss_by_level,cumulative_economic_loss,y_lo,y_hi,S,R,T,beta_levels,k,beta,Npop,dT,Xc),3,[],[],[],[],[1 1 1],[5 5 5],[],[1 2 3]);
        end
        
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
    
    cumulative_economic_loss = cumulative_economic_loss + E*economic_loss_by_level(Uinp);
    
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
hold on;
plot(tsim,1e7.*ones(size(tsim)),'k--')
legend('E','I','Q','R','D','P')
subplot(2,1,2)
plot(tsim,Umpc)
legend('U')

%% Export fig
time1 = datetime(t_start):datetime(datestr(datenum(t_start)+kMax));
figure(3)
subplot(2,1,1)
plot(time1,Yc(:,4),time1,y_hi.*ones(length(Yc),1),'k--','LineWidth',1.2)
xlim([min(time1),max(time1)])
xtickformat('dd-MMM-yy')
subplot(2,1,2)
plot(time1,Umpc,'LineWidth',1.2)
xlim([min(time1),max(time1)])
xtickformat('dd-MMM-yy')

%%
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,'HMPC_economic_impact','-dpdf','-r0')