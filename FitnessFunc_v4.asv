function obj = FitnessFunc_v4(Ucon,Np,Nc,Nb,A,B,C,yLo,yHi,S,R,T,E,Betas,economic_loss,cumulative_eco_loss,k,beta,Npop,dT,X0)

    % Initialize control states and outputs
    Y_eval = zeros(size(C,1),Np-1);
    X_eval = zeros(size(A,1),Np-1);
    X_eval(:,1) = X0;
    Y_eval(:,1) = X0;
    
    
    % Repeat input according to blocking vector
    NMPC_U = zeros(size(B,2),Np-1);
    NMPC_U(:,1:Nb(1)) = repeat(Ucon(:,1),Nb(1),1);
    for j = 2:length(Nb)
        ind = sum(Nb(1:(j-1)));
        NMPC_U(:,(ind+1):(ind+Nb(j))) = repeat(Ucon(:,j),Nb(j),1);
    end
    ind = sum(Nb);
    NMPC_U(:,(ind+1):end) = NMPC_U(:,ind)*ones(1,Np-ind-1);

    % Run simulation over prediction horizon
    i = 2;

    while (i < Np)
        
        B = [-beta(k+i-1)/Npop beta(k+i-1)/Npop 0 0 0 0 0]';
        Ustep = Y_eval(1,i-1)*Y_eval(3,i-1)*Betas(NMPC_U(i));
        
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
    
    % Penalize control move changes
    objDeltaU = sum(T.*diff(Ucon).^2);
    
    % Penalize control value
    objR = sum(R.*(Ucon.^2));
    
    % Cumulative economic loss
    objEco = cumulative_eco_loss + sum(E.*economic_loss(NMPC_U(1:sum(Nb))));
    
    objLim = 0;
    for j = 1:(Np-1)
        objLim = objLim + lim_viol(j)'*S*lim_viol(j); % Limit violations weighted with slack variables
    end
    
    if k == 102
        aa = 1;
    end
    
    obj = objLim + objR + objDeltaU + objEco;