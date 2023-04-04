%% Discretize Earnings Dynamics 
% We create a standard exogenous labor life-cycle model based on the
% Earnings Dynamics of Guvenen, Karahan, Ozkan & Song (2021).

useModel=2 % Can be: 1,2,3,4,5,6
% There is also useModel 21 and 22, these are just same as 2, except using extended Farmer-Toda to target 2 and 4 moments respectively
nSigmaz=2; % Use: 2,3,4 (number of standard deviations for the max and min points of grid used for z)
% Note: nSigmaz is really just intended for use with EvaluateDiscretization.
nSigma_alpha=1; % This was originally 2, but seems to give weird things for model 6, so reduced to 1
LoadDiscertization=1 % If 0, then does the discretization, if 1 then just loads the discretization (must be run with 0 before can run with 1)
% Typically I always just load the discretization. (A code called 'JustDiscretizeAndEvaluate' creates all the discretizations)

Setmedianwage=1 % This must equal zero the first time this is run % Set the wage parameter so that the median earnings for employed in the model equals

SolveV=0 % 0 means skip value fn and get straight to stationary dist
SolveStationaryDist=0 % 0 means skip stationary dist and get straight to life-cycle profiles
StatDistIterate=0; % Can iterate (=1), or simulate (=0) [normally I always iterate, but exogenous shocks in this model are large so ran out of memory and so simulating instead]
CalculateStatistics=[0,0,0,0,0,1]; % 0 means skip, 1 means calculate (divided it into six parts to make it easier to just run specific bits)
 % 1st, 2nd, 3rd are life-cycle profiles and similar, 4th is simul panel data, 5th consumption insurance using panel data, 6th is welfare calculation
simoptions.numbersims=10^5; % Note: does numbersims per PType (I use 10^5, but this can be set lower for speed when just testing)

% A line I needed for running on the Server
addpath(genpath('./MatlabToolkits/'))
% gpuDevice(1) % reset gpu to clear out memory

% Because there are a lot of permanent types we will use option that while
% things are calculated on gpu they get stored on cpu. (Otherwise run out of gpu memory)
% Results in minor speed reduction, but managable.
vfoptions.ptypestorecpu=1;
simoptions.ptypestorecpu=1;

if CalculateStatistics(4)==0 && Setmedianwage==0
    error('Cannot calculate the median median wage as not calculating the panel data that is needed for this')
end

%% Begin setting up to use VFI Toolkit to solve
% Lets model agents from age 25 to age 100, so 76 periods

Params.agejshifter=24; % Age 25 minus one. Makes keeping track of actual age easy in terms of model age
Params.J=100-Params.agejshifter; % Number of period in life-cycle

% Grid sizes to use
n_d=0; % Labor is exogenous
n_a=2501; %2501; % Endogenous asset holdings
% (I let asset grid be 0 to 5*10^6; with 2501 points the top two points are $6000 apart.)

% Exogenous states
n_z=51; % Persistent earnings shock
if useModel==2 || useModel==5 || useModel==6 || useModel==21 || useModel==22
    n_upsilon=2; % non-employment shock
else
    n_upsilon=1; % unused (no non-employment shock)
end
n_epsilon=15; % Transitory earnings shock
% n_zupsilonepsilon=[n_z, n_upsilon, n_epsilon];
n_zupsilon=[n_z, n_upsilon];

% Permanent types
n_alpha=5; % Fixed effect
if useModel==6
    n_kappabeta=n_alpha;  % Heterogenous income profile (GKOS2021 call this beta; set below depending on model)
    % Note: discretization routine requires that n_kappabeta=n_alpha
else
    n_kappabeta=1; % unused
end

% There are n_alpha*n_kappabeta permanent types
N_i=n_alpha*n_kappabeta;
% Note that n_kappabeta=1 when there is not permanent types of kappabeta (no Heterogeneous Income Profiles)

N_j=Params.J; % Number of periods in finite horizon

save(['./SavedOutput/Main/BasicOutput',num2str(useModel),'.mat'], 'n_d','n_a','n_z','n_upsilon','n_epsilon','n_zupsilon','n_alpha','n_kappabeta','N_i','N_j')

%% Discretized earnings dynamics 
% Creates Params
% Creates exogenous states
% Creates permanent types
if LoadDiscertization==0
    DiscretizeEarningsDynamicsModel
    EvaluateDiscretization
    save(['./SavedOutput/BasicDiscretization/DiscretizedEarnings',num2str(useModel),'nSigmaz',num2str(nSigmaz),'nSigma_alpha',num2str(nSigma_alpha),'.mat']', 'zupsilon_grid_J','epsilon_grid', 'pi_zupsilon_J', 'pi_epsilon', 'jequaloneDistzupsilonepsilon', 'jequaloneDistzupsilon', 'PTypeDistParamNames', 'Params','z_grid_J','pi_z_J','epsilon_grid_J','pi_epsilon_J','Epi_upsilon_J','-v7.3');
elseif LoadDiscertization==1
    load(['./SavedOutput/BasicDiscretization/DiscretizedEarnings',num2str(useModel),'nSigmaz',num2str(nSigmaz),'nSigma_alpha',num2str(nSigma_alpha),'.mat'])
end

%%
% To use exogenous shocks that depend on age you have to add them to vfoptions and simoptions
vfoptions.z_grid_J=zupsilon_grid_J; % Note: naming of vfoptions.z_grid_J has to be exactly as is.
vfoptions.pi_z_J=pi_zupsilon_J; % Note: naming of vfoptions.z_grid_J has to be exactly as is.
vfoptions.n_e=n_epsilon;
vfoptions.e_grid_J=epsilon_grid_J; % You could just use vfoptions.e_grid=epsilon_grid, I don't do this purely as gives cleaner results this way when plotting life-cycle profiles (epsilon does not depend on age, but it is irrelevant to retirement so replace it with zeros during retirement to make things look nicer)
vfoptions.pi_e_J=pi_epsilon_J;

simoptions.z_grid_J=vfoptions.z_grid_J; % Note: naming of vfoptions.z_grid_J has to be exactly as is.
simoptions.pi_z_J=vfoptions.pi_z_J; % Note: naming of vfoptions.z_grid_J has to be exactly as is.
simoptions.n_e=n_epsilon;
simoptions.e_grid_J=epsilon_grid_J; % You can use vfoptions.e_grid_J and vfoptions.pi_e_J, but not needed here
simoptions.pi_e_J=pi_epsilon_J;


% You then just pass a 'placeholder' for z_grid and pi_z, and the commands
% will ignore these and will only use what is in vfoptions/simoptions
zupsilon_grid=zupsilon_grid_J(:,1); % Not actually used
pi_zupsilon=pi_zupsilon_J(:,:,1); % Not actually used

%% Grids
% The ^3 means that there are more points near 0 than near 1. We know from theory that the value function will be more 'curved' near zero assets,
% and putting more points near curvature (where the derivative changes the most) increases accuracy of results.
a_grid=(5*10^6)*(linspace(0,1,n_a).^3)'; % The ^3 means most points are near zero, which is where the derivative of the value fn changes most.
% Earnings for Model 1 are up to 10000 per period, for model 6 are up to 250000 per period.
% (Presumably the original estimates were just $-per-year)

if n_a==201 % This is just a thing I used to test and debug so everything was faster
    a_grid=(5*10^2)*(linspace(0,1,n_a).^3)';
end

d_grid=[];

%% Now, create the return function 
DiscountFactorParamNames={'beta','sj'};

ReturnFn=@(aprime,a,z,upsilon,epsilon,alpha,kappabeta,w,sigma,agej,Jr,pension,incomefloor,r,kappa_j,warmglow1,warmglow2,warmglow3,beta,sj,eta1,eta2,CEV) EarningsDynamics_ReturnFn(aprime,a,z,upsilon,epsilon,alpha,kappabeta,w,sigma,agej,Jr,pension,incomefloor,r,kappa_j,warmglow1,warmglow2,warmglow3,beta,sj,eta1,eta2,CEV);

%% CEV is only used for consumption equivalent variation calculations, so turn it 'off' be setting to zero
Params.CEV=0;

%% Set up taxes and income floor
% IncomeTax=eta1+eta2*log(Income)*Income, where $IncomeTax$ is the amount paid by a household with $Income$.
% This functional form is found to have a good emprical fit to the US income tax system by GunerKaygusuvVentura2014.
Params.eta1=0.099;
Params.eta2=0.035;

Params.incomefloor=0.1;
% Note that pension>incomefloor (pension=0.4 when not setting median wage)

%% Set the wage to be equal to the median for US in 2010 (GKOS2021 use 2010 as their reference year, to which they convert all nominal amounts using the PCE index)
% $26,363.55 (which is the median earnings for US in 2010, the year for which the dollars in GKOS2021 estimates are for: https://www.ssa.gov/cgi-bin/netcomp.cgi?year=2010)
targetmedianwage=26363.55;
if Setmedianwage==1
    load(['./SavedOutput/Main/MedianEarnings',num2str(useModel),'_wageequal1.mat'], 'modelmedianearnings_conditionalonemployed_wageequal1')
    Params.w=targetmedianwage/modelmedianearnings_conditionalonemployed_wageequal1;
end
% And also set the income floor relative to this
targetincomefloor=6400; % Hubbard, Skinner & Zeldes (1995) estimate the total value of various government pro-
% grams (food stamps, AFDC, housing subsidies, etc.) for a female-headed family with
% two children and no outside earnings and assets. They obtain a value of about $7,000 in
% 1984 dollars—or about $12,800 in 2010 dollars. Using the OECD equivalence scale for
% one adult plus two children, this comes to $6,400 per adult person. (this
% paragraph copied from pg 7 of Guvenen, Karahan & Ozkan (2018WP) "Consumption and Savings Under Non-Gaussian Income Risk"
targetpension=15400; 
% $15,434 was the median income of aged 65+ with no earnings and one retirement benefit in US in 2010
% https://www.ssa.gov/policy/docs/chartbooks/income_aged/2010/iac10.html
% --> Income Sources --> Median income, by receipt of earnings and retirement benefits
if Setmedianwage==1
    Params.incomefloor=targetincomefloor;
    Params.pension=targetpension;
    
    % Because there is a target final assets (for warm glow of bequests), I
    % also scale this with w else it becomes irrelevant
    Params.warmglow2=Params.warmglow2*Params.w; % Note: original w=1, so this is scaling up by the increase in w
end


%% Now solve the value function iteration problem, just to check that things are working before we go to General Equilbrium
vfoptions.verbose=1;
vfoptions.lowmemory=1;
vfoptions.paroverz=1;

if SolveV==1
    disp('Solving the Value Fn and Policy Fn')
    tic;
    % Note: z_grid and pi_z, this will be ignored due to presence of vfoptions.z_grid_J and vfoptions.pi_z_J
    [V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_zupsilon,N_j,N_i, d_grid, a_grid, zupsilon_grid, pi_zupsilon, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
    time1=toc
    
    save(['./SavedOutput/Main/VSave',num2str(useModel),'.mat'],'V','-v7.3')
    save(['./SavedOutput/Main/PolicySave',num2str(useModel),'.mat'],'Policy','-v7.3')
    save(['./SavedOutput/Main/RestSave',num2str(useModel),'.mat'],'n_zupsilon','n_epsilon','zupsilon_grid','epsilon_grid','pi_zupsilon','pi_epsilon','time1')
    
else
    load(['./SavedOutput/Main/VSave',num2str(useModel),'.mat'])
    load(['./SavedOutput/Main/PolicySave',num2str(useModel),'.mat'])
    load(['./SavedOutput/Main/RestSave',num2str(useModel),'.mat'])
end

%% Now, we want to graph Life-Cycle Profiles

%% Initial distribution of agents at birth (j=1)

% Dist over z: jequaloneDistz'
% Dist over upsilon: 
% Dist over epsilon: statdist_epsilon
% Dist over permanent type: statdist_alpha_kappabeta

% Before we plot the life-cycle profiles we have to define how agents are at age j=1. We will give them all zero assets.
PTypeDist=Params.(PTypeDistParamNames{1}); % PTypeDist equals statdist_alpha_kappabeta or statdist_alpha depending on model
jequaloneDist=zeros(n_a,prod(n_zupsilon)*n_epsilon,'gpuArray');
jequaloneDist(1,:)=jequaloneDistzupsilonepsilon';  % All agents start with zero assets
jequaloneDist=reshape(jequaloneDist,[n_a,n_zupsilon,n_epsilon]);

% With permanent types there are two ways to set up jequaloneDist; first as
% just a distribution that does not depend on i and is assumed to be same
% for all the permanent types. Second, as a structure with fields for each i.


%% We now compute the 'stationary distribution' of households
% Start with a mass of one at initial age, use the conditional survival
% probabilities sj to calculate the mass of those who survive to next
% period, repeat. Once done for all ages, normalize to one
Params.mewj=ones(1,Params.J); % Marginal distribution of households over age
for jj=2:length(Params.mewj)
    Params.mewj(jj)=Params.sj(jj-1)*Params.mewj(jj-1);
end
Params.mewj=Params.mewj./sum(Params.mewj); % Normalize to one
AgeWeightsParamNames={'mewj'}; % So VFI Toolkit knows which parameter is the mass of agents of each age
simoptions.verbose=1;

if StatDistIterate==1
    %% iterate on agent dist
    if SolveStationaryDist==1
        disp('Solving the Agent Dist')
        simoptions.parallel=6;
        tic;
        StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_zupsilon,N_j,N_i,pi_zupsilon,Params,simoptions);
        time2=toc
        save(['./SavedOutput/Main/StationaryDist',num2str(useModel),'.mat'], 'StationaryDist','jequaloneDist','time2','-v7.3')
    else
        load(['./SavedOutput/Main/StationaryDist',num2str(useModel),'.mat'])
    end
else
    %% simulate agent dist
    if SolveStationaryDist==1
        disp('Solving the Agent Dist')
        simoptions.iterate=0;
        tic;
        StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_zupsilon,N_j,N_i,pi_zupsilon,Params,simoptions);
        time2=toc
        save(['./SavedOutput/Main/StationaryDist_sim',num2str(useModel),'.mat'], 'StationaryDist','jequaloneDist','time2','-v7.3')
    else
        load(['./SavedOutput/Main/StationaryDist_sim',num2str(useModel),'.mat'])
    end
end


%% Test that Policy is not trying to leave the top of the grid
Names_i=fieldnames(V);
test_PolicyLeavingTopOfGrid=zeros(N_i,5);
for ii=1:N_i
    temp=Policy.(Names_i{ii});
    tempdist=StationaryDist.(Names_i{ii});
    
    test_PolicyLeavingTopOfGrid(ii,1)=sum(sum(sum(sum(sum(sum(temp==n_a))))));
    test_PolicyLeavingTopOfGrid(ii,2)=sum(sum(sum(sum(sum(tempdist(end-49:end,:,:,:,:))))));
    test_PolicyLeavingTopOfGrid(ii,3)=sum(sum(sum(sum(sum(tempdist(end-99:end,:,:,:,:))))));
    test_PolicyLeavingTopOfGrid(ii,4)=sum(sum(sum(sum(sum(tempdist(end-199:end,:,:,:,:))))));
    test_PolicyLeavingTopOfGrid(ii,5)=sum(sum(sum(sum(sum(tempdist(end-299:end,:,:,:,:))))));
end
save(['./SavedOutput/Main/TestPolicy',num2str(useModel),'.mat'],'test_PolicyLeavingTopOfGrid')

%% FnsToEvaluate are how we say what we want to simulate and graph
% Note: everything involving earnings needs to include a *(agej<Jr) term
% Like with return function, we have to include (aprime,a,z) as first inputs, then just any relevant parameters.
FnsToEvaluate.earnings=@(aprime,a,z,upsilon,epsilon,w,kappa_j,alpha,kappabeta,agej,Jr) w*(1-upsilon)*exp(kappa_j+alpha+kappabeta+z+epsilon)*(agej<Jr); % labor earnings
FnsToEvaluate.income=@(aprime,a,z,upsilon,epsilon,w,kappa_j,alpha,kappabeta,r,agej,Jr,incomefloor) max(w*(1-upsilon)*exp(kappa_j+alpha+kappabeta+z+epsilon)*(agej<Jr)+r*a,incomefloor); % labor earnings+ r*a
FnsToEvaluate.assets=@(aprime,a,z,upsilon,epsilon) a; % a is the current asset holdings
FnsToEvaluate.consumption=@(aprime,a,z,upsilon,epsilon,alpha,kappabeta,w,agej,Jr,pension,incomefloor,r,kappa_j,eta1,eta2) EarningsDynamics_ConsumptionFn(aprime,a,z,upsilon,epsilon,alpha,kappabeta,w,agej,Jr,pension,incomefloor,r,kappa_j,eta1,eta2); % consumption
FnsToEvaluate.z=@(aprime,a,z,upsilon,epsilon) z; % AR(1) with gaussian-mixture innovations
FnsToEvaluate.upsilon=@(aprime,a,z,upsilon,epsilon) upsilon; % non-employment shock
FnsToEvaluate.epsilon=@(aprime,a,z,upsilon,epsilon) epsilon; % transitiory shock
FnsToEvaluate.deterministicearnings=@(aprime,a,z,upsilon,epsilon,w,kappa_j,alpha,kappabeta,agej,Jr) w*exp(kappa_j+alpha+kappabeta)*(agej<Jr); % the part of earnings which is just deterministic (in terms of age and permanent type)
FnsToEvaluate.potentialearnings=@(aprime,a,z,upsilon,epsilon,w,kappa_j,alpha,kappabeta,agej,Jr) w*exp(kappa_j+alpha+kappabeta+z+epsilon)*(agej<Jr); % what earnings would be without non-employment shocks
FnsToEvaluate.disposableincome=@(aprime,a,z,upsilon,epsilon,alpha,kappabeta,w,agej,Jr,pension,incomefloor,r,kappa_j,eta1,eta2) EarningsDynamics_DisposableIncomeFn(aprime,a,z,upsilon,epsilon,alpha,kappabeta,w,agej,Jr,pension,incomefloor,r,kappa_j,eta1,eta2);

% We also use log consumption as that way the variance has a nice interpretation as a percentage
FnsToEvaluate.logconsumption=@(aprime,a,z,upsilon,epsilon,alpha,kappabeta,w,agej,Jr,pension,incomefloor,r,kappa_j,eta1,eta2) log(EarningsDynamics_ConsumptionFn(aprime,a,z,upsilon,epsilon,alpha,kappabeta,w,agej,Jr,pension,incomefloor,r,kappa_j,eta1,eta2)); % log(consumption)
FnsToEvaluate.agej=@(aprime,a,z,upsilon,epsilon,agej) agej;

% For debugging purposes
FnsToEvaluate.nextperiodassets=@(aprime,a,z,upsilon,epsilon) aprime; % aprime is next period asset holdings
FnsToEvaluate.agej=@(aprime,a,z,upsilon,epsilon,agej) agej; % check that age dependent parameters are handled correctly
FnsToEvaluate.alpha=@(aprime,a,z,upsilon,epsilon,alpha) alpha; % check that ptype dependent parameters are handled correctly

% notice that we have called these earnings, assets and consumption

%% Calculate the life-cycle profiles: by ptype
if CalculateStatistics(1)==1
    disp('Calculate the life-cycle profiles')
    simoptions.groupptypesforstats=0;
    tic;
    AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params,n_d,n_a,n_zupsilon,N_j,N_i,d_grid, a_grid, zupsilon_grid, simoptions);
    time3=toc
    
    save(['./SavedOutput/Main/AgeConditionalStats',num2str(useModel),'.mat'], 'AgeConditionalStats','-v7.3')
else
    load(['./SavedOutput/Main/AgeConditionalStats',num2str(useModel),'.mat'])
end

%% Calculate the life-cycle profiles for the working-age population by using the simoptions.agegroupings
% This is directly comparable to GKOS2021 as their sample is just the working age
if CalculateStatistics(2)==1
    disp('Calculate the life-cycle profiles for age-groupings')

    simoptions.agegroupings=[1,Params.Jr]; % Will look at stats for age-group 1 to Params.Jr-1, which is everyone of working age (Params.Jr is first year of retirement)
    simoptions.groupptypesforstats=1;
    tic;
    AgeConditionalStats_Grouped_AgeGroupings=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params,n_d,n_a,n_zupsilon,N_j,N_i,d_grid, a_grid, zupsilon_grid, simoptions);
    time4=toc
    
    save(['./SavedOutput/Main/AgeConditionalStats_Grouped_AgeGroupings',num2str(useModel),'.mat'], 'AgeConditionalStats_Grouped_AgeGroupings','-v7.3')
    
    simoptions=rmfield(simoptions,'agegroupings');
else
    load(['./SavedOutput/Main/AgeConditionalStats_Grouped_AgeGroupings',num2str(useModel),'.mat'])
end

if CalculateStatistics(3)==1
    disp('Calculate AllStats on agent dist')

    tic;
    AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params,n_d,n_a,n_zupsilon,N_j,N_i,d_grid, a_grid, zupsilon_grid, simoptions);
    time5=toc
    save(['./SavedOutput/Main/VariousStats',num2str(useModel),'.mat'], 'AllStats','-v7.3')
end


%% Create some individual household simulations
if CalculateStatistics(4)==1
    disp('Create some individual household simulations')
    tic;
    SimPanelValues=SimPanelValues_FHorz_Case1_PType(jequaloneDist,PTypeDistParamNames,Policy,FnsToEvaluate,Params,n_d,n_a,n_zupsilon,N_j,N_i,d_grid,a_grid,zupsilon_grid,pi_zupsilon, simoptions);
    time6=toc
    
    simoptions.npoints=5; % Report transition probabilities for quintiles
    RankTransitionProbabilities_Grouped=EvalPanelData_RankTransProbs(SimPanelValues,simoptions);
        
    save(['./SavedOutput/Main/SimPanelValues',num2str(useModel),'.mat'], 'SimPanelValues','RankTransitionProbabilities_Grouped','-v7.3')
    if Setmedianwage==0
        % Params.w=1;
        keep=(SimPanelValues.earnings>0);
        modelearnings_conditionalonemployed=SimPanelValues.earnings(logical(keep));
        modelmedianearnings_conditionalonemployed_wageequal1=median(modelearnings_conditionalonemployed);
        save(['./SavedOutput/Main/MedianEarnings',num2str(useModel),'_wageequal1.mat'], 'modelmedianearnings_conditionalonemployed_wageequal1')
    end
else
    load(['./SavedOutput/Main/SimPanelValues',num2str(useModel),'.mat'])
end


%% Consumption insurance calculations (based on the simulated panel data)
if CalculateStatistics(5)==1
    disp('Some calculations based on the simulated panel data')
    % Calculate consumption insurance (two coefficients following BPP2008)
    ConsumptionInsurance_BPP2008coeffs=struct();
    % Need to use t-2 to t+1: so think of 3 as first period and end-1 as last period for time t
    % Hence 4:end is t+1; 2:end-2 is t-1, and 1:end-3 is t-2
    deltalogct=SimPanelValues.logconsumption(3:end-1,:)-SimPanelValues.logconsumption(2:end-2,:); % log c_t - log c_{t-1}
    deltalogyt=log(SimPanelValues.disposableincome(3:end-1,:))-log(SimPanelValues.disposableincome(2:end-2,:));
    diff_ytp1_ytm2=log(SimPanelValues.disposableincome(4:end,:))-log(SimPanelValues.disposableincome(1:end-3,:)); % log y_{t+1} - log y_{t-2}
    delta_ytp1=log(SimPanelValues.disposableincome(4:end,:))-log(SimPanelValues.disposableincome(3:end-1,:));
    
    % Need to drop the zero income observations from this sample before calculating the BPP coefficients
    keep=logical(isfinite(deltalogyt).*isfinite(diff_ytp1_ytm2).*isfinite(delta_ytp1)); % Note: c is always >0 so deltalogct is fine
    deltalogct=deltalogct(keep);
    deltalogyt=deltalogyt(keep);
    diff_ytp1_ytm2=diff_ytp1_ytm2(keep);
    delta_ytp1=delta_ytp1(keep);
    
    covmatrix1=cov(deltalogct,diff_ytp1_ytm2);
    covmatrix2=cov(deltalogyt,diff_ytp1_ytm2);
    covmatrix3=cov(deltalogct,delta_ytp1);
    covmatrix4=cov(deltalogyt,delta_ytp1);
    BPP_coeff1=1-covmatrix1(1,2)/covmatrix2(1,2);     
    % 1- Cov(delta log c_t, log y_{t+1} - log y_{t-2})/Cov(delta log y_t, log y_{t+1} - log y_{t-2})
    BPP_coeff2=1-covmatrix3(1,2)/covmatrix4(1,2);  
    % 1- Cov(delta log c_t, delta log y_{t+1})/Cov(delta log y_t, delta log y_{t+1})
    % Note that y here is disposable income, c is consumption
    
    ConsumptionInsurance_BPP2008coeffs.persistent=BPP_coeff1;
    ConsumptionInsurance_BPP2008coeffs.transitory=BPP_coeff2;
    
    % Want to do a plot of age-conditional variance(log(consumption)) and variance(log(earnings)), need to use the restricted sample dropping zero income observations
    logC=SimPanelValues.logconsumption;
    logY=log(SimPanelValues.disposableincome);
    age=SimPanelValues.agej;
    VarLogDispIncome_LCP=zeros(Params.J,1);
    VarLogCons_LCP=zeros(Params.J,1);
    for jj=1:Params.J
        logC_j=logC(logical(age==jj));
        logY_j=logY(logical(age==jj));
        VarLogDispIncome_LCP(jj)=var(logY_j);
        VarLogCons_LCP(jj)=var(logC_j);
    end
    
    save(['./SavedOutput/Main/ConsumptionInsurance',num2str(useModel),'.mat'], 'ConsumptionInsurance_BPP2008coeffs','VarLogDispIncome_LCP','VarLogCons_LCP')
    
    % Calculate lifetime earnings inequality, following GKSW2022
    % First, restrict the sample
    earningsmin=1885; % 2013 dollars
    earningsmin=earningsmin/52250;  % According to US Census: 2013 U.S. median house-hold income was $52,250
    KeepIndicator=ones(1,simoptions.numbersims);
    for ii=1:simoptions.numbersims
        currentearningssim=SimPanelValues.earnings(:,ii);
        count=0; % GKSW2022 include households with earnings over a minimum amount (earningsmin)
        for jj=(25-Params.agejshifter):(55-Params.agejshifter) % Ages 25 to 55
            if currentearningssim(jj)>earningsmin
                count=count+1;
            end
        end
        if count<15 % GKSW2022 drop those without at least 15 observations meeting the earnings minimum
            KeepIndicator(ii)=0; % Drop
        end
    end
    varnames=fieldnames(SimPanelValues);
    for vv=1:length(varnames)
        temp=SimPanelValues.(varnames{vv});
        temp=temp(:,logical(KeepIndicator)); % Drop those with KeepIndicator=0
        SimPanelValues.(varnames{vv})=temp;
    end
    % Now, compute the lifetime earnings
    LifetimeEarningsSample=sum(SimPanelValues.earnings,1)/31;
    % Now, some inequality measures
    LorenzCurve_LifetimeEarnings=LorenzCurve_FromSampleObs(LifetimeEarningsSample);
    % GKSW2022, Figure 8, provide the std dev of log, and the interquartile range
    stddev_logLifetimeEarnings=std(log(LifetimeEarningsSample));
    LifetimeEarningsPercentiles=prctile(LifetimeEarningsSample,[10,25,50,75,90]);
    LifetimeEarnings_P75P25ratio=LifetimeEarningsPercentiles(4)/LifetimeEarningsPercentiles(2);    
    LifetimeEarnings_P90P50ratio=LifetimeEarningsPercentiles(5)/LifetimeEarningsPercentiles(3);
    LifetimeEarnings_P50P10ratio=LifetimeEarningsPercentiles(3)/LifetimeEarningsPercentiles(1);

    save(['./SavedOutput/Main/LifeTimeEarnings',num2str(useModel),'.mat'], 'LorenzCurve_LifetimeEarnings','LifetimeEarningsSample','stddev_logLifetimeEarnings','LifetimeEarnings_P75P25ratio','LifetimeEarnings_P90P50ratio','LifetimeEarnings_P50P10ratio')
end


%%
% Save various other things about the model
Names_i=fieldnames(V);
save(['./SavedOutput/Main/GeneralOutput',num2str(useModel),'.mat'], 'Params','vfoptions','simoptions','a_grid','jequaloneDist','PTypeDist','AgeWeightsParamNames','Names_i')


%% Solve deterministic version (without shocks) to use for welfare calculations.

% To compare to DeNardi, Fella & Paz-Pardo (2018) results on welfare, we need the value function at age 1.
% We also need to solve a version without shocks and keep the value function at age 1

if useModel==3
    save ./SavedOutput/DebugCEV.mat
    save('./SavedOutput/DebugCEV2.mat','Policy','StationaryDist','V','-v7.3')
end

% We turn off all the shocks, and replace the deterministic earnings
% profiles with the age-conditional mean (for each permanent type)
if CalculateStatistics(6)==1
    if useModel==6
        PTypemassvec=Params.statdist_alpha_kappabeta;
    else
        PTypemassvec=Params.statdist_alpha;
    end
    
      
    %% Before doing the deterministic versions, solve the risky one lots of time for different possible CEV
    
    % Loop over possible values of CEV and figure out which is appropriate
    CEV_vec=0:0.01:3; % Consider 0% to 300% (to the nearest percentage points)
    W_vec=zeros(length(CEV_vec),1);
    W_ii_vec=zeros(length(CEV_vec),N_i);
    for cc=1:length(CEV_vec)
        Params.CEV=CEV_vec(cc);
        W=0;

        V_cc=ValueFnFromPolicy_Case1_FHorz_PType(Policy,n_d,n_a,n_zupsilon,N_j,Names_i,d_grid,a_grid,zupsilon_grid, pi_zupsilon, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
        for ii=1:length(Names_i)
            V_ii=V_cc.(Names_i{ii});
            StationaryDist_ii=StationaryDist.(Names_i{ii});
            
            V_CEV=sum(sum(sum(sum(V_ii(:,:,:,:,1).*StationaryDist_ii(:,:,:,:,1)))));
            
            W=W+PTypemassvec(ii)*V_CEV;
            W_ii_vec(cc,ii)=V_CEV;
        end
        W_vec(cc)=W;
        
        if cc==1
            save('./SavedOutput/DebugCEV5.mat','V_cc','-v7.3')
        end
    end
    Params.CEV=0;
    
    % To check that these are doing the right thing, plot the age-conditional earnings which should be unchanged
    FnsToEvaluate_CEVcheck.earnings=@(aprime,a,z,upsilon,epsilon,w,kappa_j,alpha,kappabeta,agej,Jr) w*(1-upsilon)*exp(kappa_j+alpha+kappabeta+z+epsilon)*(agej<Jr); % labor earnings
    FnsToEvaluate_CEVcheck.upsilon=@(aprime,a,z,upsilon,epsilon,w,kappa_j,alpha,kappabeta,agej,Jr) upsilon; % non-employment
    simoptions.groupptypesforstats=0;
    AgeConditionalStats_CEVcheck=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate_CEVcheck, Params,n_d,n_a,n_zupsilon,N_j,N_i,d_grid, a_grid, zupsilon_grid, simoptions);
    
    %% Now solve the deterministic versions of the model
    
    % This problem is simple enough we can just use the default options
    vfoptions=struct();
    simoptions=struct();
    
    % Turn off z, epsilon
    n_z=1;
    z_grid=0;
    pi_z=1;
    
    n_epsilon=1;
    epsilon_grid=0;
    pi_epsilon=1;
    
    %% For those models with non-employment shocks I first do a version in which just z and epsilon are turned off (will be CEV3)
    if n_upsilon==2
        modelhasnonemployment=1;
        upsilon_grid=[0;1]; % =zupsilon_grid_J(end-1:end,1);
    else
        modelhasnonemployment=0;        
    end
    
    if modelhasnonemployment==1
        % Note that since epsilon is no longer being put in as iid using vfoptions, we need to include it in the standard exogenous shocks processes
        n_zupsilonepsilon=[n_z,n_upsilon,n_epsilon];
        zupsilonepsilon_grid=[z_grid; upsilon_grid; epsilon_grid];
        pi_zupsilonepsilon_J=zeros(2,2,Params.J); % because z and epsilon are just single points we get that the transition matrix is just that for upsilon
        for jj=1:Params.Jr-1
            pi_zupsilonepsilon_J(:,:,jj)=Epi_upsilon_J(:,:,jj);
        end
        for jj=Params.Jr:Params.J
            pi_zupsilonepsilon_J(1,1,jj)=1; % always the zero value (which is the first point on grid); note that in retirement upsilon is anyway not used
            pi_zupsilonepsilon_J(2,1,jj)=1; % always the zero value (which is the first point on grid); note that in retirement upsilon is anyway not used
        end
        vfoptions.pi_z_J=pi_zupsilonepsilon_J;
        simoptions.pi_z_J=pi_zupsilonepsilon_J;
        vfoptions.z_grid_J=zupsilonepsilon_grid.*ones(1,Params.J);
        simoptions.z_grid_J=zupsilonepsilon_grid.*ones(1,Params.J);    
        % Need a placeholder
        pi_zupsilonepsilon=pi_zupsilonepsilon_J(:,:,1);
        
        % Note: eliminate the permament types (alpha and kappa_beta) as these will be captured in the different deterministic earnings profiles.
        % Turn off the permanent types
        Params.alpha=0;
        Params.kappabeta=0;
        % Instead, these are captured by the deterministic earnings profiles
        Params=rmfield(Params,'kappa_j');
        for ii=1:N_i
            % Earnings for each age/fraction employed for each age
            % So that earnings conditional on employed are unchanged
            Params.kappa_j.(Names_i{ii})=log(AgeConditionalStats.earnings.(Names_i{ii}).Mean./(Params.w*(1-AgeConditionalStats.upsilon.(Names_i{ii}).Mean))); % log() as in model it is exp(kappa_j+...)
        end
        [V_deterministic3, Policy_deterministic3]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_zupsilonepsilon,N_j,N_i, d_grid, a_grid, zupsilonepsilon_grid, pi_zupsilonepsilon, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
        
        simoptions.parallel=2;
        jequaloneDist_deterministic3=reshape(sum(sum(jequaloneDist,4),2),[n_a,1,n_upsilon,1]);
        % Try iterating since now no longer have z and epsilon
        StationaryDist_deterministic3=StationaryDist_Case1_FHorz_PType(jequaloneDist_deterministic3,AgeWeightsParamNames,PTypeDistParamNames,Policy_deterministic3,n_d,n_a,n_zupsilonepsilon,N_j,N_i,pi_zupsilonepsilon,Params,simoptions);
        
        simoptions.groupptypesforstats=0;
        AgeConditionalStats_CEVcheck3=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_deterministic3, Policy_deterministic3, FnsToEvaluate_CEVcheck, Params,n_d,n_a,n_zupsilonepsilon,N_j,N_i,d_grid, a_grid, zupsilonepsilon_grid, simoptions);

        % Turn the options back off again
        vfoptions=struct();
        simoptions=struct();
    end
    
    
    %% Now, turn of z, epsilon and upsilon (will be CEV1)
    % Turn off upsilon
    n_upsilon=1;
    upsilon_grid=0;
    pi_upsilon=1;
    
    % Note that since epsilon is no longer being put in as iid using vfoptions, we need to include it in the standard exogenous shocks processes
    n_zupsilonepsilon=[n_z,n_upsilon,n_epsilon];
    zupsilonepsilon_grid=[z_grid; upsilon_grid; epsilon_grid];
    pi_zupsilonepsilon=1; % because z and epsilon are just single points we get that the transition matrix is just that for upsilon
    
    % Note: eliminate alpha and kappa_beta as these will be captured in the different deterministic earnings profiles.
    % Actually keeping the permanent types for now, they are just done via kappa_j
    Params.alpha=0;
    Params.kappabeta=0;
    % Instead, these are captured by the deterministic earnings profiles
    Params=rmfield(Params,'kappa_j');
    for ii=1:N_i
        Params.kappa_j.(Names_i{ii})=log(AgeConditionalStats.earnings.(Names_i{ii}).Mean/Params.w); % log() as in model it is exp(kappa_j+...)
    end
    [V_deterministic, Policy_deterministic]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_zupsilonepsilon,N_j,N_i, d_grid, a_grid, zupsilonepsilon_grid, pi_zupsilonepsilon, ReturnFn, Params, DiscountFactorParamNames, vfoptions);

    simoptions.parallel=2;
    jequaloneDist_deterministic=reshape(sum(sum(sum(jequaloneDist,4),3),2),[n_a,1]);
    % Try iterating since now no longer have z and epsilon
    StationaryDist_deterministic=StationaryDist_Case1_FHorz_PType(jequaloneDist_deterministic,AgeWeightsParamNames,PTypeDistParamNames,Policy_deterministic,n_d,n_a,n_zupsilonepsilon,N_j,N_i,pi_zupsilonepsilon,Params,simoptions);

    simoptions.groupptypesforstats=0;
    AgeConditionalStats_CEVcheck1=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_deterministic, Policy_deterministic, FnsToEvaluate_CEVcheck, Params,n_d,n_a,n_zupsilonepsilon,N_j,N_i,d_grid, a_grid, zupsilonepsilon_grid, simoptions);

    %% Now do a version eliminating z, epsilon, upsilon and permanent types (will be CEV2)
    % Do an alternative version where we also eliminate permanent types (as from the perspective of the 'behind the veil' this is also an uninsurable risk).
    Params=rmfield(Params,'kappa_j');
    Params.kappa_j=log(AgeConditionalStats.earnings.Mean/Params.w);  % log() as in model it is exp(kappa_j+...)
    % Note: no need to remove statdist_alpha (or statdist_alpha_kappabeta) as it is anyway not used by following line
    % Note: following line is not PType
    [V_deterministic2, Policy_deterministic2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_zupsilonepsilon,N_j, d_grid, a_grid, zupsilonepsilon_grid, pi_zupsilonepsilon, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
    
    simoptions.parallel=2;
    jequaloneDist_deterministic2=reshape(sum(sum(sum(jequaloneDist,4),3),2),[n_a,1]);
    % Try iterating since now no longer have z and epsilon
    StationaryDist_deterministic2=StationaryDist_FHorz_Case1(jequaloneDist_deterministic2,AgeWeightsParamNames,Policy_deterministic2,n_d,n_a,n_zupsilonepsilon,N_j,pi_zupsilonepsilon,Params,simoptions);

    simoptions.groupptypesforstats=0;
    AgeConditionalStats_CEVcheck2=LifeCycleProfiles_FHorz_Case1(StationaryDist_deterministic2, Policy_deterministic2, FnsToEvaluate_CEVcheck, [], Params,n_d,n_a,n_zupsilonepsilon,N_j,d_grid, a_grid, zupsilonepsilon_grid, simoptions);

    %% Do the welfare calculations
    W=0;
    Wd=0;
    Wd2=0;
    Wd3=0;
    
    % For debugging purposes, keep each ptype one
    W_ii=zeros(1,N_i);
    Wd_ii=zeros(1,N_i);
    Wd2_ii=zeros(1,N_i);
    Wd3_ii=zeros(1,N_i);
    
    % Note: first, get the age 1 values and average across idiosyncratic shocks
    for ii=1:length(Names_i)
        V_ii=V.(Names_i{ii});
        V_deterministic_ii=V_deterministic.(Names_i{ii});
        
        % We just want age 1
        V_ii=V_ii(:,:,:,:,1);
        V_deterministic_ii=V_deterministic_ii(:,:,:,:,1);
        V_deterministic2_ii=V_deterministic2(:,:,:,:,1); % There are no permanent types in deterministic2 so just the same for all

        % Get the relevant stationary dist for age 1 (only applies to V_ii)
        StationaryDist_ii=StationaryDist.(Names_i{ii});
        % For the rest we have turned the shocks off
        distage1=sum(sum(sum(StationaryDist_ii(:,:,:,:,1),4),3),2); % age 1, sum over, z, upsilon and epsilon

        V_ii=sum(sum(sum(sum(V_ii.*StationaryDist_ii(:,:,:,:,1)))));

        V_deterministic_ii=sum(V_deterministic_ii.*distage1); % sums over the asset dimensions
        V_deterministic2_ii=sum(V_deterministic2_ii.*distage1); % sums over the asset dimensions
        
        W=W+PTypemassvec(ii)*V_ii;
        Wd=Wd+PTypemassvec(ii)*V_deterministic_ii;
        Wd2=Wd2+PTypemassvec(ii)*V_deterministic2_ii;
        
        W_ii(ii)=V_ii;
        Wd_ii(ii)=V_deterministic_ii;
        Wd2_ii(ii)=V_deterministic2_ii;
        
        if modelhasnonemployment==1
            V_deterministic3_ii=V_deterministic3.(Names_i{ii});
            V_deterministic3_ii=V_deterministic3_ii(:,:,:,:,1);
            upsilondistage1=sum(sum(StationaryDist_ii(:,:,:,:,1),4),2); % age 1, sum over assets, z and epsilon
            V_deterministic3_ii=sum(sum(V_deterministic3_ii.*upsilondistage1)); % sum over the asset and upsilon dimensions

            Wd3=Wd3+PTypemassvec(ii)*V_deterministic3_ii;
            
            Wd3_ii(ii)=V_deterministic3_ii;
        end
    end
    
    
    CEV3=0; % Placeholder for models without upsilon
    
    [~,CEV1_index]=min(abs(W_vec-Wd));
    CEV1=CEV_vec(CEV1_index);
    [~,CEV2_index]=min(abs(W_vec-Wd2));
    CEV2=CEV_vec(CEV2_index);
    if modelhasnonemployment==1
        [~,CEV3_index]=min(abs(W_vec-Wd3));
        CEV3=CEV_vec(CEV3_index);
    end
    
    save(['./SavedOutput/Main/WelfareCEV',num2str(useModel),'.mat'],'CEV1','CEV2','CEV3','W','Wd','Wd2','Wd3','W_vec')
    % load (['./SavedOutput/Main/WelfareCEV',num2str(useModel),'.mat'])
    
    if useModel==2
        save('./SavedOutput/DebugWelfareCEV.mat','CEV1','CEV2','CEV3','W','Wd','Wd2','Wd3','W_vec')
        save('./SavedOutput/DebugCEV3.mat','V_deterministic','V_deterministic2','V_deterministic3','W_ii','Wd_ii','Wd2_ii','Wd3_ii','-v7.3')
        save('./SavedOutput/DebugCEV4.mat','AgeConditionalStats_CEVcheck', 'AgeConditionalStats_CEVcheck1', 'AgeConditionalStats_CEVcheck2', 'AgeConditionalStats_CEVcheck3')
        save('./SavedOutput/DebugCEV6.mat','StationaryDist_deterministic', 'StationaryDist_deterministic2', 'StationaryDist_deterministic3')
    elseif useModel==3
        save('./SavedOutput/DebugWelfareCEV.mat','CEV1','CEV2','W','Wd','Wd2','W_vec')
        save('./SavedOutput/DebugCEV3.mat','V_deterministic','V_deterministic2','W_ii','Wd_ii','Wd2_ii','-v7.3')
        save('./SavedOutput/DebugCEV4.mat','AgeConditionalStats_CEVcheck', 'AgeConditionalStats_CEVcheck1', 'AgeConditionalStats_CEVcheck2')
        save('./SavedOutput/DebugCEV6.mat','StationaryDist_deterministic', 'StationaryDist_deterministic2')
    end
    
    % Check that the age-conditional earnings are not really changing across these
    disp('Check age-conditional earnings')
    check11=AgeConditionalStats_CEVcheck.earnings.ptype001.Mean-AgeConditionalStats_CEVcheck1.earnings.ptype001.Mean;
    check12=AgeConditionalStats_CEVcheck.earnings.ptype002.Mean-AgeConditionalStats_CEVcheck1.earnings.ptype002.Mean;
    check13=AgeConditionalStats_CEVcheck.earnings.ptype003.Mean-AgeConditionalStats_CEVcheck1.earnings.ptype003.Mean;
    check14=AgeConditionalStats_CEVcheck.earnings.ptype004.Mean-AgeConditionalStats_CEVcheck1.earnings.ptype004.Mean;
    check15=AgeConditionalStats_CEVcheck.earnings.ptype005.Mean-AgeConditionalStats_CEVcheck1.earnings.ptype005.Mean;
    [max(abs(check11(:))),max(abs(check12(:))),max(abs(check13(:))),max(abs(check14(:))),max(abs(check15(:)))]
    
    check2=AgeConditionalStats_CEVcheck.earnings.Mean-AgeConditionalStats_CEVcheck2.earnings.Mean;
    max(abs(check2(:)))
    
    if modelhasnonemployment==1
        check31=AgeConditionalStats_CEVcheck.earnings.ptype001.Mean-AgeConditionalStats_CEVcheck3.earnings.ptype001.Mean;
        check32=AgeConditionalStats_CEVcheck.earnings.ptype002.Mean-AgeConditionalStats_CEVcheck3.earnings.ptype002.Mean;
        check33=AgeConditionalStats_CEVcheck.earnings.ptype003.Mean-AgeConditionalStats_CEVcheck3.earnings.ptype003.Mean;
        check34=AgeConditionalStats_CEVcheck.earnings.ptype004.Mean-AgeConditionalStats_CEVcheck3.earnings.ptype004.Mean;
        check35=AgeConditionalStats_CEVcheck.earnings.ptype005.Mean-AgeConditionalStats_CEVcheck3.earnings.ptype005.Mean;
        [max(abs(check31(:))),max(abs(check32(:))),max(abs(check33(:))),max(abs(check34(:))),max(abs(check35(:)))]
    end

end

