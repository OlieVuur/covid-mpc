%% Load model parameters
load('SEIQRDP_ZA_v4.mat');

%% Define model
% Evaluate time-based parameters to get constant values - can expand to
% time-variant model later
t = 1;
beta = lambdaFun(Beta2,t);
lambda = lambdaFun(Lambda,t);
kappa = kappaFun(Kappa,t);
A = [-alpha   0      0             0         0 0 0;
       0    -gamma   0             0         0 0 0;
       0     gamma -delta          0         0 0 0;
       0      0     delta  (-kappa - lambda) 0 0 0;
       0      0      0          lambda       0 0 0;
       0      0      0          kappa        0 0 0;
       alpha  0      0             0         0 0 0];
   
B = [-beta/Npop beta/Npop 0 0 0 0 0]';
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
legend('E','I','Q','R','D','P')

%% Run control simulation 

%Initialize control parameters
dT = 1;
kMax = 1000;
tsim = 0:kMax;

% Control quarantined cases Y(4)
y_hi = 1e6;
y_lo = 0;

% For now the controlled variable is the "social distancing effectiveness"
% which is a value between 0 and 100. 0 is no social distancing, 100 is
% full social distancing (no contact between people)
% In future the idea should be to connect a lockdown level to the control
% as input
u_hi = 100;
u_lo = 0;

% Initialze simulation in/outputs
Yc = zeros(kMax+1,7);
Yc(1,:) = X0;
Xc = X0';
Umpc = zeros(kMax+1,1);
Umpc(1) = 0;

% Control config variables
Np = 200;
Nc = 30; % When using blocking the control horizon is not used
Nb = [10 10 10]; % Blocking
S = 1e-6;
R = 1e3;

options = optimoptions('ga','Display', 'off', 'FitnessLimit', 0, 'MaxGenerations', 100, 'MaxTime', 5, 'FunctionTolerance', 1);

% Run control
% Setup waitbar that shows simulation progress
h = waitbar(0,'Simulating...');

for k = 2:(kMax+1)
    tic;
    % Uncontrolled portion of input
    Un = Yc(k-1,1)*Yc(k-1,3);
    
    % Run controller on first iteration and every 10th iteration after that
    if (mod(k,10) - 2) == 0
        
        % On the first run, start without an InitialPopulationMatrix, on each
        % successive run initialize with the final population from the previous
        if k == 2
            [Uc,fval,exitflag,output,final_pop] = ga(@(Ucon) FitnessFunc(Ucon,Np,Nc,Nb,A,B,C,y_lo,y_hi,S,R,dT,Xc),3,[],[],[],[],[0 0 0],[100 100 100],[],[1 2 3]);
        else
            [Uc,fval,exitflag,output,final_pop] = ga(@(Ucon) FitnessFunc(Ucon,Np,Nc,Nb,A,B,C,y_lo,y_hi,S,R,dT,Xc),3,[],[],[],[],[0 0 0],[100 100 100],[],[1 2 3],options);
        end
        %options = optimoptions('ga','InitialPopulationMatrix', final_pop, 'Display', 'off', 'FitnessLimit', 0, 'MaxGenerations', 50, 'MaxTime', 0.5, 'FunctionTolerance', 1);
        %options = optimoptions('ga','InitialPopulationMatrix', final_pop, 'Display', 'off', 'FitnessLimit', 0, 'MaxGenerations', 100, 'MaxTime', 5, 'FunctionTolerance', 1);
        Uc = 5.*Uc;
        
        Uinp = Uc(1);
        
    end
    
    Umpc(k) = Uinp;
    
    U = Un*(100 - Uinp)./100;
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

%% Plot results of control simulation

figure(1)
subplot(2,1,1)
plot(tsim,Yc(:,2:7))
legend('E','I','Q','R','D','P')
subplot(2,1,2)
plot(tsim,Umpc)
legend('U')