% DiscretizeEarningsDynamicsModel

%%
% For some reason the Params seems to be 'accumulating' kappa_j across
% runs. I have implemented this 'clear' to make sure there are no inputs
% except what is supposed to be there.
temp_agejshifter=Params.agejshifter;
temp_J=Params.J;
clear Params
Params.agejshifter=temp_agejshifter;
Params.J=temp_J;



%% Parameters

% Discount rate
Params.beta = 0.96;
% Preferences
Params.sigma = 2; % Coeff of relative risk aversion (curvature of consumption)
Params.eta = 1.5; % Curvature of leisure (This will end up being 1/Frisch elasty)
Params.psi = 10; % Weight on leisure

% Prices
Params.w=1; % Wage
Params.r=0.05; % Interest rate (0.05 is 5%)

% Demographics
Params.agej=1:1:Params.J; % Is a vector of all the agej: 1,2,3,...,J
Params.Jr=66-Params.agejshifter; % Age 65 is last working age, age 66 is retired

% Pensions
Params.pension=0.3;

% Unemployment benefits
Params.unempbenefit=0.1;

% Conditional survival probabilities: sj is the probability of surviving to be age j+1, given alive at age j
% Most countries have calculations of these (as they are used by the government departments that oversee pensions)
% In fact I will here get data on the conditional death probabilities, and then survival is just 1-death.
% Here I just use them for the US, taken from "National Vital Statistics Report, volume 58, number 10, March 2010."
% I took them from first column (qx) of Table 1 (Total Population)
% Conditional death probabilities
Params.dj=[0.006879, 0.000463, 0.000307, 0.000220, 0.000184, 0.000172, 0.000160, 0.000149, 0.000133, 0.000114, 0.000100, 0.000105, 0.000143, 0.000221, 0.000329, 0.000449, 0.000563, 0.000667, 0.000753, 0.000823,...
    0.000894, 0.000962, 0.001005, 0.001016, 0.001003, 0.000983, 0.000967, 0.000960, 0.000970, 0.000994, 0.001027, 0.001065, 0.001115, 0.001154, 0.001209, 0.001271, 0.001351, 0.001460, 0.001603, 0.001769, 0.001943, 0.002120, 0.002311, 0.002520, 0.002747, 0.002989, 0.003242, 0.003512, 0.003803, 0.004118, 0.004464, 0.004837, 0.005217, 0.005591, 0.005963, 0.006346, 0.006768, 0.007261, 0.007866, 0.008596, 0.009473, 0.010450, 0.011456, 0.012407, 0.013320, 0.014299, 0.015323,...
    0.016558, 0.018029, 0.019723, 0.021607, 0.023723, 0.026143, 0.028892, 0.031988, 0.035476, 0.039238, 0.043382, 0.047941, 0.052953, 0.058457, 0.064494,...
    0.071107, 0.078342, 0.086244, 0.094861, 0.104242, 0.114432, 0.125479, 0.137427, 0.150317, 0.164187, 0.179066, 0.194979, 0.211941, 0.229957, 0.249020, 0.269112, 0.290198, 0.312231, 1.000000]; 
% dj covers Ages 0 to 100 (since it starts from 0, we actually use from Params.agejshifter+1+1 on
Params.sj=1-Params.dj(Params.agejshifter+1+1:101); % Conditional survival probabilities
Params.sj(end)=0; % In the present model the last period (j=J) value of sj is actually irrelevant

% Warm glow of bequest
Params.warmglow1=0.3; % (relative) importance of bequests
Params.warmglow2=3; % bliss point of bequests (essentially, the target amount)
Params.warmglow3=Params.sigma; % By using the same curvature as the utility of consumption it makes it much easier to guess appropraite parameter values for the warm glow

%% The life-cycle models of GKOS2021: Set up all six discretizations

% The only part of the AR(1) process that depends on age are the mixture
% probabilities of the gaussian-mixture innovations.

% kappa_j, deterministic age-dependent (log) earnings
% z, the persistent AR(1) with gaussian-mixture innovations
% eta, the gaussian-mixture innovations to z

% alpha, the fixed effect
% kappa_beta, heterogenous income profiles (only used in model 6)

%% 
if useModel==1 || useModel==11 || useModel==12 % The three are about using different methods to discretize z
    fellagallipolipanoptions=struct();
    kirkbyoptions=struct();
    farmertodaoptions=struct();
    alphaoptions=struct();
    
    %% Create the discretized earnings for ages 25 to 65 (as this is what GKOS2021 used)
    % Then at the end fill all the retirment age values in with zeros
    
    % Age-dependent labor productivity units
    Params.kappa_j=0.740+0.337*(Params.agej/10)+0.070*(Params.agej/10).^2;
    % Note: agej=age-24, which is what GKOS2021 use. [kappa_j is what they denote g(t)]
    Params.kappa_j(Params.Jr:end)=0; % Now fill in the rest of the retirement ages with zero productivity
    
    %% Life-cycle AR(1) process with gaussian innovations
    Params.rho=1.005*ones(1,Params.Jr-1);
    Params.sigma_eta=0.134*ones(1,Params.Jr-1);
    Params.sigma_z0=0.343;
    
    if useModel==1 % Extended Tauchen (extended Rouwenhorst cannot handle rho>1)
        fellagallipolipanoptions.nSigmas=nSigmaz;
        fellagallipolipanoptions.initialj0sigma_z=Params.sigma_z0;
        % FellaGallipoliPan is extended Rouwenhorst
        [z_grid_J, pi_z_J,jequaloneDistz,otheroutputs_z] = discretizeLifeCycleAR1_FellaGallipoliPanTauchen(Params.rho,Params.sigma_eta,n_z,Params.Jr-1,fellagallipolipanoptions);
    elseif useModel==11 % Extended Farmer-Toda targetting 2 moments
        kirkbyoptions.nSigmas=nSigmaz;
        kirkbyoptions.initialj0sigma_z=Params.sigma_z0;
        kirkbyoptions.nMoments=2;
        % Kirkby is extended Farmer-Toda
        [z_grid_J, pi_z_J,jequaloneDistz,otheroutputs_z] = discretizeLifeCycleAR1_Kirkby(Params.rho,Params.sigma_eta,n_z,Params.Jr-1,kirkbyoptions); % z_grid_J is n_z-by-J, so z_grid_J(:,j) is the grid for age j
        % pi_z_J is n_z-by-n_z-by-J, so pi_z_J(:,:,j) is the transition matrix for age j
        
        % Take a look at how many moments we matched for transitions from each grid point
        figure(20)
        h = heatmap(otheroutputs_z.nMoments_grid);
        saveas(h,['./SavedOutput/EvaluateDiscretization/discretization_momentsheatmap_model',num2str(useModel),'_nz',num2str(n_z),'_nSigmaz',num2str(nSigmaz),'.pdf'])
        % Can see from the heat map that it gets 2 moments in the vast majority of grid points
        [sum(sum(otheroutputs_z.nMoments_grid==2)),numel(otheroutputs_z.nMoments_grid)] % Over 90% of grid points match all 2 moments
    elseif useModel==12 % Extended Farmer-Toda targetting 4 moments
        kirkbyoptions.nSigmas=nSigmaz;
        kirkbyoptions.initialj0sigma_z=Params.sigma_z0;
        kirkbyoptions.nMoments=4; % This is anyway the default
        % Kirkby is extended Farmer-Toda
        [z_grid_J, pi_z_J,jequaloneDistz,otheroutputs_z] = discretizeLifeCycleAR1_Kirkby(Params.rho,Params.sigma_eta,n_z,Params.Jr-1,kirkbyoptions); % z_grid_J is n_z-by-J, so z_grid_J(:,j) is the grid for age j
        
        % Take a look at how many moments we matched for transitions from each grid point
        figure(20)
        h = heatmap(otheroutputs_z.nMoments_grid);
        saveas(h,['./SavedOutput/EvaluateDiscretization/discretization_momentsheatmap_model',num2str(useModel),'_nz',num2str(n_z),'_nSigmaz',num2str(nSigmaz),'.pdf'])
        % Can see from the heat map that it gets 4 moments in the vast majority of grid points
        [sum(sum(otheroutputs_z.nMoments_grid==4)),numel(otheroutputs_z.nMoments_grid)] % Over 90% of grid points match all 4 moments
    end
    
    %% Even though there is no non-employment shock upsilon, I will set it up as a single grid point because this makes it easier to use the same codes as for the other models
    upsilon_grid=0; % Since we have (1-upsilon)*ln(y), upsilon=0 is essentially turning it off
    % pi_upsilon=1;
    upsilon_grid_J=zeros(1,Params.Jr-1);
    pi_zupsilon_J=pi_z_J;
    
    jequaloneDistzupsilon=jequaloneDistz;
    
    %% Transitory gaussian iid shocks
    farmertodaoptions.nSigmas=nSigmaz;
    Params.sigma_epsilon=0.696;
    [epsilon_grid,pi_epsilon] = discretizeAR1_FarmerToda(0,0,Params.sigma_epsilon,n_epsilon,farmertodaoptions);

    % Note that epsilon is iid, so
    pi_epsilon=pi_epsilon(1,:)';
    
    jequaloneDistzupsilonepsilon=zeros(n_z*1*n_epsilon,1);
    for e_c=1:n_epsilon
        jequaloneDistzupsilonepsilon((1:1:1*n_z)+(n_z*n_upsilon)*(e_c-1))=jequaloneDistzupsilon*pi_epsilon(e_c);
    end
    
    %% Fill in all the retirement age values
    % Fill in the rest of the retirement ages with zero productivity
    z_grid_J=[z_grid_J,zeros(n_z,Params.J-Params.Jr+1)];
    upsilon_grid_J=[upsilon_grid_J,zeros(1,Params.J-Params.Jr+1)];
    epsilon_grid_J=[epsilon_grid.*ones(1,Params.Jr-1),zeros(n_epsilon,Params.J-Params.Jr+1)];
    pi_epsilon_J=pi_epsilon*ones(1,Params.J);
    
    zupsilon_grid_J=[z_grid_J; upsilon_grid_J];
    
    % Fill in the retirement ages with uniform transition probabilities (these are anyway irrelevant)
    temp=pi_zupsilon_J;
    pi_zupsilon_J=ones(n_z,n_z,Params.J)/n_z;
    pi_zupsilon_J(:,:,1:Params.Jr-1)=temp;
    % Fill in the retirement ages for z with identity matrix (these are anyway irrelevant)
    pi_z_J2=pi_z_J;
    pi_z_J=repmat(eye(n_z,n_z),1,1,Params.J);
    pi_z_J(:,:,1:Params.Jr-1)=pi_z_J2;
    
    %% Fixed effect
    alphaoptions.nSigmas=nSigma_alpha;
    Params.sigma_alpha=1.182;
    [alpha_grid,pi_alpha] = discretizeAR1_FarmerToda(0,0,Params.sigma_alpha,n_alpha,alphaoptions);
    
    statdist_alpha=pi_alpha(1,:)';
    
    Params.alpha=alpha_grid; 
    
    % No kappabeta
    Params.kappabeta=0; % 0 makes it effectively disappear from model
        
    Params.statdist_alpha=statdist_alpha;
    PTypeDistParamNames={'statdist_alpha'};
    
elseif useModel==2 || useModel==21 || useModel==22 % The three are about using different methods to discretize z
    kirkbyoptions=struct();
    fellagallipolipanoptions=struct();
    farmertodaoptions=struct();
    alphaoptions=struct();
    
    %% Create the discretized earnings for ages 25 to 65 (as this is what GKOS2021 used)
    % Then at the end fill all the retirment age values in with zeros
    
    % Age-dependent labor productivity units
    Params.kappa_j=2.569+0.766*(Params.agej/10)-0.152*(Params.agej/10).^2;
    % Note: agej=age-24, which is what GKOS2021 use. [kappa_j is what they denote g(t)]
    Params.kappa_j(Params.Jr:end)=0; % Now fill in the rest of the retirement ages with zero productivity
    
    %% Life-cycle AR(1) process with gaussian innovations
    Params.rho=0.967*ones(1,Params.Jr-1);
    Params.sigma_eta=0.197*ones(1,Params.Jr-1);
    Params.sigma_z0=0.563;
    if useModel==2 % Extended Rouwenhorst
        fellagallipolipanoptions.nSigmas=nSigmaz;
        fellagallipolipanoptions.initialj0sigma_z=Params.sigma_z0;
        % FellaGallipoliPan is extended Rouwenhorst
        [z_grid_J, pi_z_J,jequaloneDistz,otheroutputs_z] = discretizeLifeCycleAR1_FellaGallipoliPan(Params.rho,Params.sigma_eta,n_z,Params.Jr-1,fellagallipolipanoptions);
    elseif useModel==21 % Extended Farmer-Toda targetting 2 moments
        kirkbyoptions.nSigmas=nSigmaz;
        kirkbyoptions.initialj0sigma_z=Params.sigma_z0;
        kirkbyoptions.nMoments=2;
        % Kirkby is extended Farmer-Toda
        [z_grid_J, pi_z_J,jequaloneDistz,otheroutputs_z] = discretizeLifeCycleAR1_Kirkby(Params.rho,Params.sigma_eta,n_z,Params.Jr-1,kirkbyoptions); % z_grid_J is n_z-by-J, so z_grid_J(:,j) is the grid for age j
        % pi_z_J is n_z-by-n_z-by-J, so pi_z_J(:,:,j) is the transition matrix for age j
        
        % Take a look at how many moments we matched for transitions from each grid point
        figure(20)
        h = heatmap(otheroutputs_z.nMoments_grid);
        saveas(h,['./SavedOutput/EvaluateDiscretization/discretization_momentsheatmap_model',num2str(useModel),'_nz',num2str(n_z),'_nSigmaz',num2str(nSigmaz),'.pdf'])
        % Can see from the heat map that it gets 2 moments in the vast majority of grid points
        [sum(sum(otheroutputs_z.nMoments_grid==2)),numel(otheroutputs_z.nMoments_grid)] % Over 90% of grid points match all 2 moments
    elseif useModel==22 % Extended Farmer-Toda targetting 4 moments
        kirkbyoptions.nSigmas=nSigmaz;
        kirkbyoptions.initialj0sigma_z=Params.sigma_z0;
        kirkbyoptions.nMoments=4; % This is anyway the default
        % Kirkby is extended Farmer-Toda
        [z_grid_J, pi_z_J,jequaloneDistz,otheroutputs_z] = discretizeLifeCycleAR1_Kirkby(Params.rho,Params.sigma_eta,n_z,Params.Jr-1,kirkbyoptions); % z_grid_J is n_z-by-J, so z_grid_J(:,j) is the grid for age j
        % pi_z_J is n_z-by-n_z-by-J, so pi_z_J(:,:,j) is the transition matrix for age j
        
        % Take a look at how many moments we matched for transitions from each grid point
        figure(20)
        h = heatmap(otheroutputs_z.nMoments_grid);
        saveas(h,['./SavedOutput/EvaluateDiscretization/discretization_momentsheatmap_model',num2str(useModel),'_nz',num2str(n_z),'_nSigmaz',num2str(nSigmaz),'.pdf'])
        % Can see from the heat map that it gets 4 moments in the vast majority of grid points
        [sum(sum(otheroutputs_z.nMoments_grid==4)),numel(otheroutputs_z.nMoments_grid)] % Over 90% of grid points match all 4 moments
    end
    
    %% Non-Employment Shocks
    % Note that the non-employment shocks depend on z
    Params.lambda=0.03;
%     Params.xi=@(agej,z) -3.036-0.917*agej-5.397*z-4.442*agej*z;
    Params.xi=@(agej,z) -3.036-0.917*(agej/10)-5.397*z-4.442*(agej/10)*z;
    Params.prob_upsilon=@(xi) exp(xi)/(1+exp(xi));
    % Need two states, 0 and min(1,exp(lambda))
    % Note that the probabilities, xi, are based on agej and z in the same
    % period (which is next period in the transition matrix)
    % Create grid
    upsilon_grid_J=[zeros(1,Params.Jr-1); min(1,exp(Params.lambda))*ones(1,Params.Jr-1)];
    % pi_upsilon_J cannot be defined independent of pi_z_J, so create the
    % joint transition matrix for (z, upsilon)
    pi_zupsilon_J=zeros(n_z*2,n_z*2,Params.Jr-1);
    for jj=1:Params.Jr-2
        for z_c=1:n_z
            xi=Params.xi(jj+1,z_grid_J(z_c,jj+1));
            prob_upsilon=Params.prob_upsilon(xi);
            % Note all that matters for (next period) upsilon is next period z and next period age
            pi_zupsilon_J(1:n_z,z_c,jj)=pi_z_J(:,z_c,jj)*(1-prob_upsilon); % Corresponds to upsilon=0
            pi_zupsilon_J(n_z+1:2*n_z,z_c,jj)=pi_z_J(:,z_c,jj)*(1-prob_upsilon);  % Corresponds to upsilon=0
            pi_zupsilon_J(1:n_z,n_z+z_c,jj)=pi_z_J(:,z_c,jj)*prob_upsilon; % Corresponds to upsilon=min(1,exp(lambda))
            pi_zupsilon_J(n_z+1:2*n_z,n_z+z_c,jj)=pi_z_J(:,z_c,jj)*prob_upsilon; % Corresponds to upsilon=min(1,exp(lambda))
        end
    end
    pi_zupsilon_J(:,:,Params.Jr-1)=ones(n_z*2,n_z*2)/(n_z*2); % Note that the agej=Jr-1 transition is irrelevant in any case
    
    %
    jequaloneDistzupsilon=zeros(n_z*2,1);
    for z_c=1:n_z
        xi=Params.xi(1,z_grid_J(z_c,1));
        prob_upsilon=Params.prob_upsilon(xi);
        % Note all that matters for (next period) upsilon is next period z and next period age
        jequaloneDistzupsilon(z_c)=jequaloneDistz(z_c)*(1-prob_upsilon); % Corresponds to upsilon=0
        jequaloneDistzupsilon(n_z+z_c)=jequaloneDistz(z_c)*prob_upsilon; % Corresponds to upsilon=min(1,exp(lambda))
    end
    
    %% Transitory gaussian iid shocks
    farmertodaoptions.nSigmas=nSigmaz;
    Params.sigma_epsilon=0.163;
    [epsilon_grid,pi_epsilon] = discretizeAR1_FarmerToda(0,0,Params.sigma_epsilon,n_epsilon,farmertodaoptions);

    % Note that epsilon is iid, so
    pi_epsilon=pi_epsilon(1,:)';
    
    jequaloneDistzupsilonepsilon=zeros(n_z*2*n_epsilon,1);
    for e_c=1:n_epsilon
        jequaloneDistzupsilonepsilon((1:1:2*n_z)+(n_z*n_upsilon)*(e_c-1))=jequaloneDistzupsilon*pi_epsilon(e_c);
    end
    
    %% Fill in all the retirement age values
    % Fill in the rest of the retirement ages with zero productivity
    z_grid_J=[z_grid_J,zeros(n_z,Params.J-Params.Jr+1)];
    upsilon_grid_J=[upsilon_grid_J,zeros(2,Params.J-Params.Jr+1)];
    epsilon_grid_J=[epsilon_grid.*ones(1,Params.Jr-1),zeros(n_epsilon,Params.J-Params.Jr+1)];
    pi_epsilon_J=pi_epsilon*ones(1,Params.J);

    zupsilon_grid_J=[z_grid_J; upsilon_grid_J];
    
    % Fill in the retirement ages with uniform transition probabilities (these are anyway irrelevant)
    temp=pi_zupsilon_J;
    pi_zupsilon_J=ones(n_z*2,n_z*2,Params.J)/(n_z*2);
    pi_zupsilon_J(:,:,1:Params.Jr-1)=temp;
    % Fill in the retirement ages for z with identity matrix (these are anyway irrelevant)
    pi_z_J2=pi_z_J;
    pi_z_J=repmat(eye(n_z,n_z),1,1,Params.J);
    pi_z_J(:,:,1:Params.Jr-1)=pi_z_J2;
    
    %% Fixed effect
    alphaoptions.nSigmas=nSigma_alpha;
    Params.sigma_alpha=0.655;
    [alpha_grid,pi_alpha] = discretizeAR1_FarmerToda(0,0,Params.sigma_alpha,n_alpha,alphaoptions);

    statdist_alpha=pi_alpha(1,:)';
    
    Params.alpha=alpha_grid; 
    
    % No kappabeta
    Params.kappabeta=0; % 0 makes it effectively disappear from model
        
    Params.statdist_alpha=statdist_alpha;
    PTypeDistParamNames={'statdist_alpha'};
    
elseif useModel==3
    kirkbyoptions=struct();
    farmertodaoptions=struct();
    alphaoptions=struct();
    
    % Create the discretized earnings for ages 25 to 65 (as this is what GKOS2021 used)
    % Then at the end fill all the retirment age values in with zeros
    
    %% Age-dependent labor productivity units
    Params.kappa_j=2.547-0.144*(Params.agej/10)-0.059*(Params.agej/10).^2;
    % Note: agej=age-24, which is what GKOS2021 use. [kappa_j is what they denote g(t)]
    Params.kappa_j(Params.Jr:end)=0; % Now fill in the rest of the retirement ages with zero productivity

    %% Persistent AR(1) process with gaussian-mixture innovations
    Params.rho=1.010*ones(1,Params.Jr-1);
    Params.mixprobs_eta=[0.05;1-0.05].*ones(1,Params.Jr-1);
    % Gaussian-mixture innovations
    Params.mu_eta1=-1*ones(1,Params.Jr-1); 
    % Note mew_eta2 is determined to make mean(eta)=0 (from mew_et1 and mixture-probabilities)
    kirkbyoptions.setmixturemutoenforcezeromean=1; % If missing mean for one of the gaussian mixture innovations this is used to set it so that mean equals zero (assumes the missing one is the 'last' one)
    Params.sigma_eta1=1.421*ones(1,Params.Jr-1);
    Params.sigma_eta2=0.010*ones(1,Params.Jr-1);
    % Initial agej=0 distribution of the life-cycle AR(1) process with gaussian-mixture innovations
    Params.sigma_z0=0.213; % N(0,sigma_z0)
    kirkbyoptions.initialj0sigma_z=Params.sigma_z0;
    
    % Now we can set up grids and transtition probabilities for all ages
    Params.sigma_eta=[Params.sigma_eta1;Params.sigma_eta2];
    Params.mu_eta=Params.mu_eta1; % Using kirkbyoptions.setmixturemutoenforcezeromean=1
    
    kirkbyoptions.nSigmas=nSigmaz;
    [z_grid_J, pi_z_J,jequaloneDistz,otheroutputs_z] = discretizeLifeCycleAR1wGM_Kirkby(Params.rho,Params.mixprobs_eta,Params.mu_eta,Params.sigma_eta,n_z,Params.Jr-1,kirkbyoptions); % z_grid_J is n_z-by-J, so z_grid_J(:,j) is the grid for age j
    % pi_z_J is n_z-by-n_z-by-J, so pi_z_J(:,:,j) is the transition matrix for age j    
    
    % Take a look at how many moments we matched for transitions from each grid point
    figure(20)
    h = heatmap(otheroutputs_z.nMoments_grid);
    saveas(h,['./SavedOutput/EvaluateDiscretization/discretization_momentsheatmap_model',num2str(useModel),'_nz',num2str(n_z),'_nSigmaz',num2str(nSigmaz),'.pdf'])
    % Can see from the heat map that it gets 4 moments in the vast majority of grid points
    [sum(sum(otheroutputs_z.nMoments_grid==4)),numel(otheroutputs_z.nMoments_grid)] % Over 90% of grid points match all 4 moments
        
    %% Even though there is no non-employment shock upsilon, I will set it up as a single grid point because this makes it easier to use the same codes as for the other models
    upsilon_grid=0; % Since we have (1-upsilon)*ln(y), upsilon=0 is essentially turning it off
    % pi_upsilon=1;
    upsilon_grid_J=zeros(1,Params.Jr-1);
    pi_zupsilon_J=pi_z_J;
    
    jequaloneDistzupsilon=jequaloneDistz;
    
    %% Transitory shock
    Params.p_epsilon=[0.118;1-0.118];
    Params.mu_epsilon1=-0.826;
    Params.mu_epsilon2=-(Params.p_epsilon(1)*Params.mu_epsilon1)/Params.p_epsilon(2); % Rearrange: p.*mu=0
    Params.mu_epsilon=[Params.mu_epsilon1;Params.mu_epsilon2];
    Params.sigma_epsilon1=1.549;
    Params.sigma_epsilon2=0.020;
    Params.sigma_epsilon=[Params.sigma_epsilon1;Params.sigma_epsilon2];
    
    farmertodaoptions.nSigmas=nSigmaz;
    farmertodaoptions.method='GMQ'; % Did not succeed in hitting all 4 moments when using evenly spaced grid; GMQ grid can hit all 4
    [epsilon_grid,pi_epsilon, otheroutputs_epsilon] = discretizeAR1wGM_FarmerToda(0,0,Params.p_epsilon,Params.mu_epsilon,Params.sigma_epsilon,n_epsilon,farmertodaoptions);
    % otheroutputs.nMoments_grid shows that only succeed in matching first three moments
    
    % Note that epsilon is iid, so
    pi_epsilon=pi_epsilon(1,:)';
    
    jequaloneDistzupsilonepsilon=zeros(n_z*1*n_epsilon,1);
    for e_c=1:n_epsilon
        jequaloneDistzupsilonepsilon((1:1:1*n_z)+(n_z*n_upsilon)*(e_c-1))=jequaloneDistzupsilon*pi_epsilon(e_c);
    end
    
    %% Fill in all the retirement age values
    % Fill in the rest of the retirement ages with zero productivity
    z_grid_J=[z_grid_J,zeros(n_z,Params.J-Params.Jr+1)];
    upsilon_grid_J=[upsilon_grid_J,zeros(1,Params.J-Params.Jr+1)];
    epsilon_grid_J=[epsilon_grid.*ones(1,Params.Jr-1),zeros(n_epsilon,Params.J-Params.Jr+1)];
    pi_epsilon_J=pi_epsilon*ones(1,Params.J);

    zupsilon_grid_J=[z_grid_J; upsilon_grid_J];
    
    % Fill in the retirement ages with uniform transition probabilities (these are anyway irrelevant)
    temp=pi_zupsilon_J;
    pi_zupsilon_J=ones(n_z,n_z,Params.J)/n_z;
    pi_zupsilon_J(:,:,1:Params.Jr-1)=temp;
    % Fill in the retirement ages for z with identity matrix (these are anyway irrelevant)
    pi_z_J2=pi_z_J;
    pi_z_J=repmat(eye(n_z,n_z),1,1,Params.J);
    pi_z_J(:,:,1:Params.Jr-1)=pi_z_J2;
    
    %% Fixed effect alpha
    alphaoptions.nSigmas=nSigma_alpha;
    Params.sigma_alpha=0.273;
    [alpha_grid,pi_alpha] = discretizeAR1_FarmerToda(0,0,Params.sigma_alpha,n_alpha,alphaoptions);

    statdist_alpha=pi_alpha(1,:)';
    
    Params.alpha=alpha_grid; 
    
    % No kappabeta
    Params.kappabeta=0; % 0 makes it effectively disappear from model
    
    Params.statdist_alpha=statdist_alpha;
    PTypeDistParamNames={'statdist_alpha'};
    
elseif useModel==4
    kirkbyoptions=struct();
    farmertodaoptions=struct();
    alphaoptions=struct();
    
    %% Age-dependent labor productivity units
    Params.kappa_j=2.176+0.169*(Params.agej/10)-0.1*(Params.agej/10).^2;
    % Note: agej=age-24, which is what GKOS2021 use. [kappa_j is what they denote g(t)]
    Params.kappa_j(Params.Jr:end)=0; % Now fill in the rest of the retirement ages with zero productivity

    %% Life-cycle AR(1) process with gaussian-mixture innovations
    kirkbyoptions.customGKOSmodel4=1;
    
    Params.rho=0.992*ones(1,Params.Jr-1);
    % Gaussian-mixture innovations
    Params.mu_eta1=-1*ones(1,Params.Jr-1); % Note mew_eta2 is below determined to make mean=0 (from mew_et1 and mixture-probabilities)
    Params.sigma_eta1=1.07*ones(1,Params.Jr-1);
    Params.sigma_eta2=0.032*ones(1,Params.Jr-1);
    xi_fn=@(agej,zlag) (-0.474+1.961*(agej/10)-3.183*zlag-0.187*(agej/10)*zlag);
    Params.mixprobs_eta=@(agej,zlag) [(exp(xi_fn(agej,zlag)))/(1+exp(xi_fn(agej,zlag))); 1-(exp(xi_fn(agej,zlag)))/(1+exp(xi_fn(agej,zlag)))]; % GKOS2021 denote this p_z
    % Note: mixprobs_eta is logistic function of xi: e(xi)/(1+e(xi))
    
    % Initial agej=1 distribution of the life-cycle AR(1) process with gaussian-mixture innovations
    Params.sigma_z0=0.446; % N(0,sigma_z0)
    % Need to set up the age j=1 grid
    kirkbyoptions.initialj0sigma_z=Params.sigma_z0;
    
    % Now we can set up grids and transtition probabilities for all ages
    Params.sigma_eta=[Params.sigma_eta1;Params.sigma_eta2];
    Params.mu_eta=Params.mu_eta1;
    kirkbyoptions.setmixturemutoenforcezeromean=1; % If missing mean for one of the gaussian mixture innovations this is used to set it so that mean equals zero (assumes the missing one is the 'last' one)

    kirkbyoptions.nSigmas=nSigmaz;
    [z_grid_J, pi_z_J,jequaloneDistz,otheroutputs_z] = discretizeLifeCycleAR1wGM_KirkbyFull(Params.rho,Params.mixprobs_eta,Params.mu_eta,Params.sigma_eta,n_z,Params.Jr-1,kirkbyoptions); % z_grid_J is n_z-by-J, so z_grid_J(:,j) is the grid for age j
    % pi_z_J is n_z-by-n_z-by-J, so pi_z_J(:,:,j) is the transition matrix for age j
    
    % For debugging purposes:
    rho=Params.rho;
    mixprobs_i=Params.mixprobs_eta; mu_i=Params.mu_eta; 
    sigma_i=Params.sigma_eta;
    znum=n_z;J=Params.Jr-1;
    
    % Take a look at how many moments we matched for transitions from each grid point
    figure(20)
    h = heatmap(otheroutputs_z.nMoments_grid);
    saveas(h,['./SavedOutput/EvaluateDiscretization/discretization_momentsheatmap_model',num2str(useModel),'_nz',num2str(n_z),'_nSigmaz',num2str(nSigmaz),'.pdf'])
    % Check how often it hits four grid points
    [sum(sum(otheroutputs_z.nMoments_grid==4)),numel(otheroutputs_z.nMoments_grid)] % Roughly 80% of grid points match all 4 moments
    
    %% Non-employment shock (not used for Model 4)
    n_upsilon=1;
    upsilon_grid_J=zeros(1,Params.Jr-1);
    pi_upsilon=1;
    
    pi_zupsilon_J=pi_z_J;
    jequaloneDistzupsilon=jequaloneDistz;

    disp('Mass of jequaloneDistz')
    sum(jequaloneDistz(:))
    
    %% transitory shock epsilon, gaussian-mixture
    Params.p_epsilon=[0.088; 1-0.088];
    Params.mu_epsilon1=0.311;
    Params.mu_epsilon2=-(Params.p_epsilon(1)*Params.mu_epsilon1)/Params.p_epsilon(2); % Rearrange: p.*mu=0
    Params.mu_epsilon=[Params.mu_epsilon1;Params.mu_epsilon2];
    Params.sigma_epsilon1=0.795;
    Params.sigma_epsilon2=0.020;
    Params.sigma_epsilon=[Params.sigma_epsilon1;Params.sigma_epsilon2];
    
    farmertodaoptions.nSigmas=nSigmaz;
    farmertodaoptions.method='GMQ'; % Did not succeed in hitting all 4 moments when using evenly spaced grid; GMQ grid can hit all 4
    [epsilon_grid,pi_epsilon,otheroutputs_epsilon] = discretizeAR1wGM_FarmerToda(0,0,Params.p_epsilon,Params.mu_epsilon,Params.sigma_epsilon,n_epsilon,farmertodaoptions);
    % otheroutputs.nMoments_grid shows that only succeed in matching first three moments
    
    % Note that epsilon is iid, so
    pi_epsilon=pi_epsilon(1,:)';
    
    jequaloneDistzupsilonepsilon=zeros(n_z*n_epsilon,1);
    for e_c=1:n_epsilon
        jequaloneDistzupsilonepsilon((1:1:n_z)+(n_z)*(e_c-1))=jequaloneDistzupsilon*pi_epsilon(e_c);
    end
        
    %% Fill in all the retirement age values
    % Fill in the rest of the retirement ages with zero productivity
    z_grid_J=[z_grid_J,zeros(n_z,Params.J-Params.Jr+1)];
    upsilon_grid_J=[upsilon_grid_J,zeros(1,Params.J-Params.Jr+1)];
    epsilon_grid_J=[epsilon_grid.*ones(1,Params.Jr-1),zeros(n_epsilon,Params.J-Params.Jr+1)];
    pi_epsilon_J=pi_epsilon*ones(1,Params.J);

    zupsilon_grid_J=[z_grid_J; upsilon_grid_J];
    
    % Fill in the retirement ages with uniform transition probabilities (these are anyway irrelevant)
    temp=pi_zupsilon_J;
    pi_zupsilon_J=ones(n_z,n_z,Params.J)/(n_z);
    pi_zupsilon_J(:,:,1:Params.Jr-1)=temp;
    % Fill in the retirement ages for z with identity matrix (these are anyway irrelevant)
    pi_z_J2=pi_z_J;
    pi_z_J=repmat(eye(n_z,n_z),1,1,Params.J);
    pi_z_J(:,:,1:Params.Jr-1)=pi_z_J2;
    
    %% Fixed effect alpha
    alphaoptions.nSigmas=nSigma_alpha;
    Params.sigma_alpha=0.473;
    [alpha_grid,pi_alpha] = discretizeAR1_FarmerToda(0,0,Params.sigma_alpha,n_alpha,alphaoptions);

    statdist_alpha=pi_alpha(1,:)';
    
    Params.alpha=alpha_grid; 
    
    % No kappabeta
    Params.kappabeta=0; % 0 makes it effectively disappear from model
    
    Params.statdist_alpha=statdist_alpha;
    PTypeDistParamNames={'statdist_alpha'};
    
elseif useModel==5
    farmertodaoptions=struct();
    kirkbyoptions=struct();
    alphaoptions=struct();
    
	% Create the discretized earnings for ages 25 to 65 (as this is what GKOS2021 used)
    % Then at the end fill all the retirment age values in with zeros
    
    %% Age-dependent labor productivity units
    Params.kappa_j=2.746+0.624*(Params.agej/10)-0.167*(Params.agej/10).^2;
    % Note: agej=age-24, which is what GKOS2021 use. [kappa_j is what they denote g(t)]
    Params.kappa_j(Params.Jr:end)=0; % Now fill in the rest of the retirement ages with zero productivity

    %% Persistent AR(1) process with gaussian-mixture innovations
    
    Params.rho=0.991*ones(1,Params.Jr-1);
    Params.mixprobs_eta=[0.176;1-0.176].*ones(1,Params.Jr-1);
    % Gaussian-mixture innovations
    Params.mu_eta1=-0.524*ones(1,Params.Jr-1); 
    % Note mew_eta2 is determined to make mean(eta)=0 (from mew_et1 and mixture-probabilities)
    kirkbyoptions.setmixturemutoenforcezeromean=1; % If missing mean for one of the gaussian mixture innovations this is used to set it so that mean equals zero (assumes the missing one is the 'last' one)
    Params.sigma_eta1=0.113*ones(1,Params.Jr-1);
    Params.sigma_eta2=0.046*ones(1,Params.Jr-1);
    % Initial agej=0 distribution of the life-cycle AR(1) process with gaussian-mixture innovations
    Params.sigma_z0=0.450; % N(0,sigma_z0)
    kirkbyoptions.initialj0sigma_z=Params.sigma_z0;
    
    % Now we can set up grids and transtition probabilities for all ages
    Params.sigma_eta=[Params.sigma_eta1;Params.sigma_eta2];
    Params.mu_eta=Params.mu_eta1; % Using kirkbyoptions.setmixturemutoenforcezeromean=1
    
    kirkbyoptions.nSigmas=nSigmaz;
    [z_grid_J, pi_z_J,jequaloneDistz,otheroutputs_z] = discretizeLifeCycleAR1wGM_Kirkby(Params.rho,Params.mixprobs_eta,Params.mu_eta,Params.sigma_eta,n_z,Params.Jr-1,kirkbyoptions); % z_grid_J is n_z-by-J, so z_grid_J(:,j) is the grid for age j
    % pi_z_J is n_z-by-n_z-by-J, so pi_z_J(:,:,j) is the transition matrix for age j    
    
    % Take a look at how many moments we matched for transitions from each grid point
    figure(20)
    h = heatmap(otheroutputs_z.nMoments_grid);
    saveas(h,['./SavedOutput/EvaluateDiscretization/discretization_momentsheatmap_model',num2str(useModel),'_nz',num2str(n_z),'_nSigmaz',num2str(nSigmaz),'.pdf'])
    % Can see from the heat map that it gets 4 moments in the vast majority of grid points
    [sum(sum(otheroutputs_z.nMoments_grid==4)),numel(otheroutputs_z.nMoments_grid)] % Over 90% of grid points match all 4 moments
        
    %% Non-Employment Shocks
    % Note that the non-employment shocks depend on z
    Params.lambda=0.016;
%     Params.xi=@(agej,z) -2.495-1.037*agej-5.051*z-1.087*agej*z;
    Params.xi=@(agej,z) -2.495-1.037*(agej/10)-5.051*z-1.087*(agej/10)*z;
    Params.prob_upsilon=@(xi) exp(xi)/(1+exp(xi));
    % Need two states, 0 and min(1,exp(lambda))
    % Note that the probabilities, xi, are based on agej and z in the same
    % period (which is next period in the transition matrix)
    % Create grid
    upsilon_grid_J=[zeros(1,Params.Jr-1); min(1,exp(Params.lambda))*ones(1,Params.Jr-1)];
    % pi_upsilon_J cannot be defined independent of pi_z_J, so create the
    % joint transition matrix for (z, upsilon)
    pi_zupsilon_J=zeros(n_z*2,n_z*2,Params.Jr-1);
    for jj=1:Params.Jr-2
        for z_c=1:n_z
            xi=Params.xi(jj+1,z_grid_J(z_c,jj+1));
            prob_upsilon=Params.prob_upsilon(xi);
            % Note all that matters for (next period) upsilon is next period z and next period age
            pi_zupsilon_J(1:n_z,z_c,jj)=pi_z_J(:,z_c,jj)*(1-prob_upsilon); % Corresponds to upsilon=0
            pi_zupsilon_J(n_z+1:2*n_z,z_c,jj)=pi_z_J(:,z_c,jj)*(1-prob_upsilon);  % Corresponds to upsilon=0
            pi_zupsilon_J(1:n_z,n_z+z_c,jj)=pi_z_J(:,z_c,jj)*prob_upsilon; % Corresponds to upsilon=min(1,exp(lambda))
            pi_zupsilon_J(n_z+1:2*n_z,n_z+z_c,jj)=pi_z_J(:,z_c,jj)*prob_upsilon; % Corresponds to upsilon=min(1,exp(lambda))
        end
    end
    pi_zupsilon_J(:,:,Params.Jr-1)=ones(n_z*2,n_z*2)/(n_z*2); % Note that the agej=Jr-1 transition is irrelevant in any case
    
    %
    jequaloneDistzupsilon=zeros(n_z*2,1);
    for z_c=1:n_z
        xi=Params.xi(1,z_grid_J(z_c,1));
        prob_upsilon=Params.prob_upsilon(xi);
        % Note all that matters for (next period) upsilon is next period z and next period age
        jequaloneDistzupsilon(z_c)=jequaloneDistz(z_c)*(1-prob_upsilon); % Corresponds to upsilon=0
        jequaloneDistzupsilon(n_z+z_c)=jequaloneDistz(z_c)*prob_upsilon; % Corresponds to upsilon=min(1,exp(lambda))
    end
    
    %% Transitory shock
    Params.p_epsilon=[0.044;1-0.044];
    Params.mu_epsilon1=0.134;
    Params.mu_epsilon2=-(Params.p_epsilon(1)*Params.mu_epsilon1)/Params.p_epsilon(2); % Rearrange: p.*mu=0
    Params.mu_epsilon=[Params.mu_epsilon1;Params.mu_epsilon2];
    Params.sigma_epsilon1=0.762;
    Params.sigma_epsilon2=0.055;
    Params.sigma_epsilon=[Params.sigma_epsilon1;Params.sigma_epsilon2];
    
    farmertodaoptions.nSigmas=nSigmaz;
    farmertodaoptions.method='GMQ'; % Did not succeed in hitting all 4 moments when using evenly spaced grid; GMQ grid can hit all 4
    [epsilon_grid,pi_epsilon, otheroutputs_epsilon] = discretizeAR1wGM_FarmerToda(0,0,Params.p_epsilon,Params.mu_epsilon,Params.sigma_epsilon,n_epsilon,farmertodaoptions);
    % otheroutputs.nMoments_grid shows that only succeed in matching first three moments
    
    % Note that epsilon is iid, so
    pi_epsilon=pi_epsilon(1,:)';
    
    jequaloneDistzupsilonepsilon=zeros(n_z*2*n_epsilon,1);
    for e_c=1:n_epsilon
        jequaloneDistzupsilonepsilon((1:1:2*n_z)+(n_z*n_upsilon)*(e_c-1))=jequaloneDistzupsilon*pi_epsilon(e_c);
    end
    
    %% Fill in all the retirement age values
    % Fill in the rest of the retirement ages with zero productivity
    z_grid_J=[z_grid_J,zeros(n_z,Params.J-Params.Jr+1)];
    upsilon_grid_J=[upsilon_grid_J,zeros(2,Params.J-Params.Jr+1)];
    epsilon_grid_J=[epsilon_grid.*ones(1,Params.Jr-1),zeros(n_epsilon,Params.J-Params.Jr+1)];
    pi_epsilon_J=pi_epsilon*ones(1,Params.J);

    zupsilon_grid_J=[z_grid_J; upsilon_grid_J];
    
    % Fill in the retirement ages with uniform transition probabilities (these are anyway irrelevant)
    temp=pi_zupsilon_J;
    pi_zupsilon_J=ones(n_z*2,n_z*2,Params.J)/(n_z*2);
    pi_zupsilon_J(:,:,1:Params.Jr-1)=temp;
    % Fill in the retirement ages for z with identity matrix (these are anyway irrelevant)
    pi_z_J2=pi_z_J;
    pi_z_J=repmat(eye(n_z,n_z),1,1,Params.J);
    pi_z_J(:,:,1:Params.Jr-1)=pi_z_J2;
    
    %% Fixed effect alpha
    alphaoptions.nSigmas=nSigma_alpha;
    Params.sigma_alpha=0.472;
    [alpha_grid,pi_alpha] = discretizeAR1_FarmerToda(0,0,Params.sigma_alpha,n_alpha,alphaoptions);

    statdist_alpha=pi_alpha(1,:)';
    
    Params.alpha=alpha_grid; 
    
    % No kappabeta
    Params.kappabeta=0; % 0 makes it effectively disappear from model
    
    Params.statdist_alpha=statdist_alpha;
    PTypeDistParamNames={'statdist_alpha'};
    
elseif useModel==6
    kirkbyoptions=struct();
    farmertodaoptions=struct();
    alphaoptions=struct();
    
	% Create the discretized earnings for ages 25 to 65 (as this is what GKOS2021 used)
    % Then at the end fill all the retirment age values in with zeros
    
    %% Age-dependent labor productivity units
    Params.kappa_j=2.581+0.812*(Params.agej/10)-0.185*(Params.agej/10).^2;
    % Note: agej=age-24, which is what GKOS2021 use. [kappa_j is what they denote g(t)]
    Params.kappa_j(Params.Jr:end)=0; % Now fill in the rest of the retirement ages with zero productivity
    
    %% Persistent AR(1) process with gaussian-mixture innovations
    Params.rho=0.959*ones(1,Params.Jr-1);
    Params.mixprobs_eta=[0.407;1-0.407].*ones(1,Params.Jr-1);
    % Gaussian-mixture innovations
    Params.mu_eta1=-0.085*ones(1,Params.Jr-1); 
    % Note mew_eta2 is determined to make mean(eta)=0 (from mew_et1 and mixture-probabilities)
    kirkbyoptions.setmixturemutoenforcezeromean=1; % If missing mean for one of the gaussian mixture innovations this is used to set it so that mean equals zero (assumes the missing one is the 'last' one)
    Params.sigma_eta1=0.364*ones(1,Params.Jr-1);
    Params.sigma_eta2=0.069*ones(1,Params.Jr-1);
    % Initial agej=0 distribution of the life-cycle AR(1) process with gaussian-mixture innovations
    Params.sigma_z0=0.714; % N(0,sigma_z0)
    kirkbyoptions.initialj0sigma_z=Params.sigma_z0;
    
    % Now we can set up grids and transtition probabilities for all ages
    Params.sigma_eta=[Params.sigma_eta1;Params.sigma_eta2];
    Params.mu_eta=Params.mu_eta1; % Using kirkbyoptions.setmixturemutoenforcezeromean=1
    
    kirkbyoptions.nSigmas=nSigmaz;
    [z_grid_J, pi_z_J,jequaloneDistz,otheroutputs_z] = discretizeLifeCycleAR1wGM_Kirkby(Params.rho,Params.mixprobs_eta,Params.mu_eta,Params.sigma_eta,n_z,Params.Jr-1,kirkbyoptions); % z_grid_J is n_z-by-J, so z_grid_J(:,j) is the grid for age j
    % pi_z_J is n_z-by-n_z-by-J, so pi_z_J(:,:,j) is the transition matrix for age j    
    
    % Take a look at how many moments we matched for transitions from each grid point
    figure(20)
    h = heatmap(otheroutputs_z.nMoments_grid);
    saveas(h,['./SavedOutput/EvaluateDiscretization/discretization_momentsheatmap_model',num2str(useModel),'_nz',num2str(n_z),'_nSigmaz',num2str(nSigmaz),'.pdf'])
    % Can see from the heat map that it gets 4 moments in the vast majority of grid points
    [sum(sum(otheroutputs_z.nMoments_grid==4)),numel(otheroutputs_z.nMoments_grid)] % Over 90% of grid points match all 4 moments
    
    %% Non-Employment Shocks
    % Note that the non-employment shocks depend on z
    Params.lambda=0.0001;
%     Params.xi=@(agej,z) -3.353-0.859*agej-5.034*z-2.895*agej*z;
    Params.xi=@(agej,z) -3.353-0.859*(agej/10)-5.034*z-2.895*(agej/10)*z;
    Params.prob_upsilon=@(xi) exp(xi)/(1+exp(xi));
    % Need two states, 0 and min(1,exp(lambda))
    % Note that the probabilities, xi, are based on agej and z in the same
    % period (which is next period in the transition matrix)
    % Create grid
    upsilon_grid_J=[zeros(1,Params.Jr-1); min(1,exp(Params.lambda))*ones(1,Params.Jr-1)];
    % pi_upsilon_J cannot be defined independent of pi_z_J, so create the
    % joint transition matrix for (z, upsilon)
    pi_zupsilon_J=zeros(n_z*2,n_z*2,Params.Jr-1);
    for jj=1:Params.Jr-2
        for z_c=1:n_z
            xi=Params.xi(jj+1,z_grid_J(z_c,jj+1));
            prob_upsilon=Params.prob_upsilon(xi);
            % Note all that matters for (next period) upsilon is next period z and next period age
            pi_zupsilon_J(1:n_z,z_c,jj)=pi_z_J(:,z_c,jj)*(1-prob_upsilon); % Corresponds to upsilon=0
            pi_zupsilon_J(n_z+1:2*n_z,z_c,jj)=pi_z_J(:,z_c,jj)*(1-prob_upsilon);  % Corresponds to upsilon=0
            pi_zupsilon_J(1:n_z,n_z+z_c,jj)=pi_z_J(:,z_c,jj)*prob_upsilon; % Corresponds to upsilon=min(1,exp(lambda))
            pi_zupsilon_J(n_z+1:2*n_z,n_z+z_c,jj)=pi_z_J(:,z_c,jj)*prob_upsilon; % Corresponds to upsilon=min(1,exp(lambda))
        end
    end
    pi_zupsilon_J(:,:,Params.Jr-1)=ones(n_z*2,n_z*2)/(n_z*2); % Note that the agej=Jr-1 transition is irrelevant in any case
    
    %
    jequaloneDistzupsilon=zeros(n_z*2,1);
    for z_c=1:n_z
        xi=Params.xi(1,z_grid_J(z_c,1));
        prob_upsilon=Params.prob_upsilon(xi);
        % Note all that matters for (next period) upsilon is next period z and next period age
        jequaloneDistzupsilon(z_c)=jequaloneDistz(z_c)*(1-prob_upsilon); % Corresponds to upsilon=0
        jequaloneDistzupsilon(n_z+z_c)=jequaloneDistz(z_c)*prob_upsilon; % Corresponds to upsilon=min(1,exp(lambda))
    end
    
    %% Transitory shock
    
    Params.p_epsilon=[0.130;1-0.130];
    Params.mu_epsilon1=0.271;
    Params.mu_epsilon2=-(Params.p_epsilon(1)*Params.mu_epsilon1)/Params.p_epsilon(2); % Rearrange: p.*mu=0
    Params.mu_epsilon=[Params.mu_epsilon1;Params.mu_epsilon2];
    Params.sigma_epsilon1=0.795;
    Params.sigma_epsilon2=0.020;
    Params.sigma_epsilon=[Params.sigma_epsilon1;Params.sigma_epsilon2];
    
    farmertodaoptions.nSigmas=nSigmaz;
    farmertodaoptions.method='GMQ'; % Did not succeed in hitting all 4 moments when using evenly spaced grid; GMQ grid can hit all 4
    [epsilon_grid,pi_epsilon, otheroutputs_epsilon] = discretizeAR1wGM_FarmerToda(0,0,Params.p_epsilon,Params.mu_epsilon,Params.sigma_epsilon,n_epsilon,farmertodaoptions);
    % otheroutputs.nMoments_grid shows that only succeed in matching first three moments
    
    % Note that epsilon is iid, so
    pi_epsilon=pi_epsilon(1,:)';
    
    jequaloneDistzupsilonepsilon=zeros(n_z*2*n_epsilon,1);
    for e_c=1:n_epsilon
        jequaloneDistzupsilonepsilon((1:1:2*n_z)+(n_z*n_upsilon)*(e_c-1))=jequaloneDistzupsilon*pi_epsilon(e_c);
    end
    
    %% Fill in all the retirement age values
    % Fill in the rest of the retirement ages with zero productivity
    z_grid_J=[z_grid_J,zeros(n_z,Params.J-Params.Jr+1)];
    upsilon_grid_J=[upsilon_grid_J,zeros(2,Params.J-Params.Jr+1)];
    epsilon_grid_J=[epsilon_grid.*ones(1,Params.Jr-1),zeros(n_epsilon,Params.J-Params.Jr+1)];
    pi_epsilon_J=pi_epsilon*ones(1,Params.J);

    zupsilon_grid_J=[z_grid_J; upsilon_grid_J];
    
    % Fill in the retirement ages with uniform transition probabilities (these are anyway irrelevant)
    temp=pi_zupsilon_J;
    pi_zupsilon_J=ones(n_z*2,n_z*2,Params.J)/(n_z*2);
    pi_zupsilon_J(:,:,1:Params.Jr-1)=temp;
    % Fill in the retirement ages for z with identity matrix (these are anyway irrelevant)
    pi_z_J2=pi_z_J;
    pi_z_J=repmat(eye(n_z,n_z),1,1,Params.J);
    pi_z_J(:,:,1:Params.Jr-1)=pi_z_J2;

    
    %% Fixed effect and Heterogenous Income Profiles
    alphaoptions.nSigmas=nSigma_alpha;
    % alpha and kappabeta are correlated, and modelled as permanent types
    % Parameters for the fixed effect (Model 6 of GKOS2021)
    Params.sigma_alpha=0.3;
    % Parameters for the heterogenous income profile (Model 6 of GKOS2021)
    Params.sigma_kappabeta=0.196; % GKOS2021 call this beta, I call it kappa_beta (as it modifies kappa, and because beta is the discount factor)
    Params.corr_alpha_kappabeta=0.798;
    % Both alpha and kappabeta are treated as 'Permanent Types'

    % Create discrete quadrature grids for alpha and kappa_beta, which are joint normally distributed
    Params.covar_alpha_kappabeta=Params.corr_alpha_kappabeta*(Params.sigma_alpha*Params.sigma_kappabeta);
    VarCovMatrix_alpha_kappabeta=[Params.sigma_alpha^2, Params.covar_alpha_kappabeta; Params.covar_alpha_kappabeta, Params.sigma_kappabeta]; % SigmaSq
    [alpha_kappabeta_grid,pi_alpha_kappabeta] = discretizeVAR1_FarmerToda(zeros(2,1),zeros(2,2),VarCovMatrix_alpha_kappabeta,n_alpha,alphaoptions);
    % Note: n_alpha=n_kappabeta is required by discretizeVAR1_FarmerToda()
    % We have no actual use for the transition matrix, as (alpha,kappabeta) is iid
    statdist_alpha_kappabeta=pi_alpha_kappabeta(1,:)';
    
    % alpha_kappabeta_grid is jointly-determined
    Params.alpha=alpha_kappabeta_grid(:,1); % first row is alpha
    Params.kappabeta=alpha_kappabeta_grid(:,2).*[linspace(1,Params.Jr-1,Params.Jr-1)/10,zeros(1,Params.J-Params.Jr+1)]; % second row is kappabeta
    % Note: the linspace(1,Params.Jr-1,Params.Jr-1)/10 term is 't=(age-24)/10' in GKOS2021 notation
    
    Params.statdist_alpha_kappabeta=statdist_alpha_kappabeta;
    PTypeDistParamNames={'statdist_alpha_kappabeta'};
    
elseif useModel==13
    % Table D.V: Additional Specifications
    % Reports a Model 13, which is like the baseline Model 6, but with
    % addition of age-dependence (and lag dependence) of the innovations of
    % the persistent AR(1) with gaussian-mixture shocks.
    % But it was not estimated using the same weights on the target moments
    % (see footnote of same table), and nor are the deterministic life-cycle
    % profile nor heterogenous income profile parameters reported.
end




%% The following is unrelated to the discretization. It is simply to create a Table containing model parameter values.

% Table of life-cycle model parameters (except the earnings process)
FID = fopen('./SavedOutput/LatexInputs/Table_LifeCycleModelParameters.tex', 'w');
fprintf(FID, 'Parameters of the Life-Cycle Model \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcc} \n \\hline \\hline \n');
fprintf(FID, ' Total periods & $J$ & %8.2f \\\\ \n', Params.J);
fprintf(FID, ' Retirement period & $Jr$ & %8.2f \\\\ \n', Params.Jr);
fprintf(FID, ' Age in years at j=0 & $Jr$ & %8.2f \\\\ \n', Params.agejshifter);
fprintf(FID, ' \\multicolumn{3}{l}{Preferences} \\\\ \n');
fprintf(FID, ' Discount factor & $\\beta$ & %8.2f \\\\ \n', Params.beta);
fprintf(FID, ' CRRA & $\\sigma$ & %8.2f \\\\ \n', Params.sigma);
fprintf(FID, ' Curv. of leisure & $\\eta$ & %8.2f \\\\ \n', Params.eta);
fprintf(FID, ' Weight on leisure & $\\psi$ & %8.2f \\\\ \n', Params.psi);
fprintf(FID, ' \\multicolumn{3}{l}{Prices and policies} \\\\ \n');
fprintf(FID, ' Wage & $w$ & %8.2f \\\\ \n', Params.w);
fprintf(FID, ' Interest rate & $r$ & %8.2f \\\\ \n', Params.r);
fprintf(FID, ' Pension & $pension$ & %8.2f \\\\ \n', Params.pension);
fprintf(FID, ' Unemp. benefits & $unempbenefit$ & %8.2f \\\\ \n', Params.unempbenefit);
fprintf(FID, ' \\multicolumn{3}{l}{Warm glow of bequests} \\\\ \n');
fprintf(FID, ' relative importance of bequests & $warmglow1$ & %8.1f \\\\ \n', Params.warmglow1);
fprintf(FID, ' bliss point of bequests & $warmglow2$ & %8.1f \\\\ \n', Params.warmglow2);
fprintf(FID, ' curvature for bequests & $warmglow3$ & %8.1f \\\\ \n', Params.warmglow3);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
% fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% fprintf(FID, 'Note: notice. \n');
% fprintf(FID, '}} \\end{minipage}');
fclose(FID);


jstr='';
istr='';
cstr='';
for jj=1:10
    jstr=[jstr,'& %8.6f '];
    istr=[istr,'& %i '];
    cstr=[cstr,'c'];
end
jstrend='';
istrend='';
for jj=71:Params.J
    jstrend=[jstrend,'& %8.6f '];
    istrend=[istrend,'& %i '];
end
for jj=Params.J+1:80
    jstrend=[jstrend,'& '];
    istrend=[istrend,'& '];
end

temp_agej=1:1:Params.J;
temp_ageyear=Params.agejshifter+(1:1:Params.J);
temp_dj=Params.dj(Params.agejshifter+1+1:101);
temp_sj=Params.sj;
temp_cumsj=cumprod(Params.sj);

% Table of age conditional death and survival probabilities
FID = fopen('./SavedOutput/LatexInputs/Table_LifeCycleModelSurvivalProbs.tex', 'w');
fprintf(FID, 'Age-Conditional Survival/Death Probabilities for the Life-Cycle Model \\\\ \n');
fprintf(FID, ['\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lc',cstr,'} \n \\hline \\hline \n']);
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(1:10));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(1:10));
fprintf(FID, [' Death & $d_j$ ',jstr,' \\\\ \n'], temp_dj(1:10));
fprintf(FID, [' Survival & $s_j=1-d_j$ ',jstr,' \\\\ \n'], temp_sj(1:10));
fprintf(FID, [' Cumulative survival & $\\Pi_{j=1}^J s_j$ ',jstr,' \\\\ \\hline \n'], temp_cumsj(1:10));
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(11:20));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(11:20));
fprintf(FID, [' Death & $d_j$ ',jstr,' \\\\ \n'], temp_dj(11:20));
fprintf(FID, [' Survival & $s_j=1-d_j$ ',jstr,' \\\\ \n'], temp_sj(11:20));
fprintf(FID, [' Cumulative survival & $\\Pi_{j=1}^J s_j$ ',jstr,' \\\\ \\hline \n'], temp_cumsj(11:20));
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(21:30));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(21:30));
fprintf(FID, [' Death & $d_j$ ',jstr,' \\\\ \n'], temp_dj(21:30));
fprintf(FID, [' Survival & $s_j=1-d_j$ ',jstr,' \\\\ \n'], temp_sj(21:30));
fprintf(FID, [' Cumulative survival & $\\Pi_{j=1}^J s_j$ ',jstr,' \\\\ \\hline \n'], temp_cumsj(21:30));
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(31:40));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(31:40));
fprintf(FID, [' Death & $d_j$ ',jstr,' \\\\ \n'], temp_dj(31:40));
fprintf(FID, [' Survival & $s_j=1-d_j$ ',jstr,' \\\\ \n'], temp_sj(31:40));
fprintf(FID, [' Cumulative survival & $\\Pi_{j=1}^J s_j$ ',jstr,' \\\\ \\hline \n'], temp_cumsj(31:40));
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(41:50));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(41:50));
fprintf(FID, [' Death & $d_j$ ',jstr,' \\\\ \n'], temp_dj(41:50));
fprintf(FID, [' Survival & $s_j=1-d_j$ ',jstr,' \\\\ \n'], temp_sj(41:50));
fprintf(FID, [' Cumulative survival & $\\Pi_{j=1}^J s_j$ ',jstr,' \\\\ \\hline \n'], temp_cumsj(41:50));
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(51:60));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(51:60));
fprintf(FID, [' Death & $d_j$ ',jstr,' \\\\ \n'], temp_dj(51:60));
fprintf(FID, [' Survival & $s_j=1-d_j$ ',jstr,' \\\\ \n'], temp_sj(51:60));
fprintf(FID, [' Cumulative survival & $\\Pi_{j=1}^J s_j$ ',jstr,' \\\\ \\hline \n'], temp_cumsj(51:60));
fprintf(FID, ['Model age & $j$ ',istr,' \\\\ \n'], temp_agej(61:70));
fprintf(FID, ['Age in Years & ',istr,' \\\\ \n'], temp_ageyear(61:70));
fprintf(FID, [' Death & $d_j$ ',jstr,' \\\\ \n'], temp_dj(61:70));
fprintf(FID, [' Survival & $s_j=1-d_j$ ',jstr,' \\\\ \n'], temp_sj(61:70));
fprintf(FID, [' Cumulative survival & $\\Pi_{j=1}^J s_j$ ',jstr,' \\\\ \\hline \n'], temp_cumsj(61:70));
fprintf(FID, ['Model age & $j$ ',istrend,' \\\\ \n'], temp_agej(71:end));
fprintf(FID, ['Age in Years & ',istrend,' \\\\ \n'], temp_ageyear(71:end));
fprintf(FID, [' Death & $d_j$ ',jstrend,' \\\\ \n'], temp_dj(71:end));
fprintf(FID, [' Survival & $s_j=1-d_j$ ',jstrend,' \\\\ \n'], temp_sj(71:end));
fprintf(FID, [' Cumulative survival & $\\Pi_{j=1}^J s_j$ ',jstrend,' \\\\ \\hline \n'], temp_cumsj(71:end));
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Source: $s_j$ is the probability of surviving to be age $j+1$, given alive at age $j$. Age-conditional death probabilites for the US are taken from "National Vital Statistics Report, volume 58, number 10, March 2010.", first column (qx) of Table 1 (Total Population). Conditional survival probabilities are calculated as then $s_j=1-d_j$. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);


