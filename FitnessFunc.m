function obj = FitnessFunc(Ucon,Np,Nc,Nb,A,B,C,yLo,yHi,S,R,dT,X0)

    % Initialize control states and outputs
    Y_eval = zeros(size(C,1),Np-1);
    X_eval = zeros(size(A,1),Np-1);
    X_eval(:,1) = X0;
    Y_eval(:,1) = X0;
    
    Ucon = 1.*Ucon;
    % Repeat input according to blocking vector
    NMPC_U = zeros(size(B,2),Np-1);
    NMPC_U(:,1:Nb(1)) = repeat(Ucon(:,1),Nb(1),1);
    for k = 2:length(Nb)
        ind = sum(Nb(1:(k-1)));
        NMPC_U(:,(ind+1):(ind+Nb(k))) = repeat(Ucon(:,k),Nb(k),1);
    end
    ind = sum(Nb);
    NMPC_U(:,(ind+1):end) = NMPC_U(:,ind)*ones(1,Np-ind-1);

    % Run simulation over prediction horizon
    i = 2;

    while (i < Np)
        
        Ustep = Y_eval(1,i-1)*Y_eval(3,i-1)*(100 - NMPC_U(i-1))./100;
        
        % Runge Kutta 4th order for simulation
        k_1 = A*X_eval(:,i-1) + B*Ustep;
        k_2 = A*(X_eval(:,i-1)+0.5*dT*k_1) + B*Ustep;
        k_3 = A*(X_eval(:,i-1)+0.5*dT*k_2) + B*Ustep;
        k_4 = A*(X_eval(:,i-1)+dT*k_3) + B*Ustep;
        % States
        X_eval(:,i) = X_eval(:,i-1) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dT;
        Y_eval(:,i) = C*X_eval(:,i);
        
        i = i + 1;
    end

    % Objective function evaluation
    
    % Penalize high/low limit violations
    high_lim_viol = max(Y_eval(4,:) - yHi,0);
    low_lim_viol = max(yLo - Y_eval(4,:),0);
    lim_viol = high_lim_viol + low_lim_viol;
    
    objLim = 0;
    objR = 0;
    for j = 1:(Np-1)
        objLim = objLim + lim_viol(j)'*S*lim_viol(j); % Limit violations weighted with slack variables
        objR = objR + (NMPC_U(j))'*R*(NMPC_U(j)); % Penalize control move away from 100 weighted
    end
    
    obj = objLim + objR;