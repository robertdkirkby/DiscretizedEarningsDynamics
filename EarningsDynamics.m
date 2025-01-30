%% Discretize Earnings Dynamics 
% We create a standard exogenous labor life-cycle model based on the
% Earnings Dynamics of Guvenen, Karahan, Ozkan & Song (2021).

% A line I needed for running on the Server
addpath(genpath('./MatlabToolkits/'))
% gpuDevice(1) % reset gpu to clear out memory
gpuDeviceCount % For some reason server failed to find gpu without this

useModel=6 % Can be: 1,2,3,4,5,6
% There is also useModel 21 and 22, these are just same as 2, except using extended Farmer-Toda to target 2 and 4 moments respectively
nSigmaz=2; % Use: 2,3,4 (number of standard deviations for the max and min points of grid used for z)
% Note: nSigmaz is really just intended for use with EvaluateDiscretization.
nSigma_alpha=1; % This was originally 2, but seems to give weird things for model 6, so reduced to 1

% Typically I always just load the discretization. (A code called 'JustDiscretizeAndEvaluate' creates all the discretizations)
LoadDiscertization=1 % If 0, then does the discretization, if 1 then just loads the discretization (must be run with 0 before can run with 1)
% First, some calibration
preCalib=0 % Pre-calibrate some parameters (note: before this preCalib, the asset grid is ridiculous, but doesn't matter since earnings is exogenous as the precalib only needs earnings)
doCalib=0 % Calibrate some parameters
% Next, solve the model
SolveVandStationaryDist=0 % if 0, loads saved versions from previous run
% Now, go through all sorts of model outputs
CalculateStatistics=[1,1,0,0,0,0] % 0 means skip, 1 means calculate (divided it into six parts to make it easier to just run specific bits)
% 1st, 2nd, 3rd are life-cycle profiles and similar, 4th is simul panel data, 5th consumption insurance using panel data, 6th is welfare calculation

% Because there are a lot of permanent types we will use option that while
% things are calculated on gpu they get stored on cpu. (Otherwise run out of gpu memory)
% Results in minor speed reduction, but managable.
vfoptions.ptypestorecpu=1;
simoptions.ptypestorecpu=1;
% Settings for the panel data
simoptions.numbersims=10^5; % Note: does numbersims per PType (I use 10^5, but this can be set lower for speed when just testing)
% Problem is monotone, so to solve faster use ddivide-and-conquer
vfoptions.divideandconquer=1;
vfoptions.level1n=25;


%% Begin setting up to use VFI Toolkit to solve
% Lets model agents from age 25 to age 100, so 76 periods

Params.agejshifter=24; % Age 25 minus one. Makes keeping track of actual age easy in terms of model age
Params.J=100-Params.agejshifter; % Number of period in life-cycle

% Grid sizes to use
n_d=0; % Labor is exogenous
n_a=2501; % Endogenous asset holdings
% (I let asset grid be 0 to 5*10^6; with 2501 points the top two points are $6000 apart.)

% Exogenous states
n_z=21; % 51 % Persistent earnings shock
if useModel==2 || useModel==5 || useModel==6 || useModel==21 || useModel==22
    n_upsilon=2; % non-employment shock
else
    n_upsilon=1; % unused (no non-employment shock)
end
n_epsilon=9; % 15 % Transitory earnings shock
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

%% Discretized earnings dynamics 
% Creates Params (many are just pre-calibration initial values)
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
vfoptions.n_e=n_epsilon;
vfoptions.e_grid=epsilon_grid_J; % You could just use vfoptions.e_grid=epsilon_grid, I don't do this purely as gives cleaner results this way when plotting life-cycle profiles (epsilon does not depend on age, but it is irrelevant to retirement so replace it with zeros during retirement to make things look nicer)
vfoptions.pi_e=pi_epsilon_J;

simoptions.n_e=n_epsilon;
simoptions.e_grid=epsilon_grid_J; % You can use vfoptions.e_grid_J and vfoptions.pi_e_J, but not needed here
simoptions.pi_e=pi_epsilon_J;

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

ReturnFn=@(aprime,a,z,upsilon,epsilon,alpha,kappabeta,w,sigma,agej,Jr,pension,incomefloor,r,kappa_j,warmglow1,warmglow2,warmglow3,beta,sj,eta1,eta2,Jbeq,CEV) EarningsDynamics_ReturnFn(aprime,a,z,upsilon,epsilon,alpha,kappabeta,w,sigma,agej,Jr,pension,incomefloor,r,kappa_j,warmglow1,warmglow2,warmglow3,beta,sj,eta1,eta2,Jbeq,CEV);

%% CEV is only used for consumption equivalent variation calculations, so turn it 'off' be setting to zero
Params.CEV=0;

% Following parameters were created in DiscretizeEarningsDynamicsModel.m
% (along with the parameters that are obviously about the shocks themselves)
% Demographics
% Params.agej=1:1:Params.J; % Is a vector of all the agej: 1,2,3,...,J
% Params.Jr=66-Params.agejshifter; % Age 65 is last working age, age 66 is retired

%% Parameters (except those for earnings dynamics where are declared in 'DiscretizeEarningsDynamicsModel.m')
% Discount rate
Params.beta = 0.96;
% Preferences
Params.sigma = 2; % Coeff of relative risk aversion (curvature of consumption)

% Prices
Params.w=1; % Wage (just an initial guess)
Params.r=0.05; % Interest rate (0.05 is 5%)

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

%% Bequests parameters
Params.Jbeq=80-Params.agejshifter; % consider people dying ages 80+ to leave bequests

% Warm glow of bequest (using functional form of De Nardi, 2004)
Params.warmglow1=2; % (relative) importance of bequests
Params.warmglow2=0.6; % extent to which bequests are a 'luxury' good
Params.warmglow3=Params.sigma; % By using the same curvature as the utility of consumption it makes it much easier to guess appropraite parameter values for the warm glow
% Note: originally I had warmglow1=0.3 and warmglow2=3, these were updated
% to something closer to what the calibration delivers so as to speed the codes.

%% Set up taxes, income floor and pensions
% IncomeTax=eta1+eta2*log(Income)*Income, where $IncomeTax$ is the amount paid by a household with $Income$.
% This functional form is found to have a good emprical fit to the US income tax system by GunerKaygusuvVentura2014.
Params.eta1=0.099;
Params.eta2=0.035;

% Pensions
Params.pension=0.3; % Just an initial guess
% Income floor (income in non-employment state during working age)
Params.incomefloor=0.1; % Just an initial guess

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


%% Calibrate the life-cycle model based on calibration targets (part 1, precalibration)
% Because this is an exogenous labor model, we can just solve once, and based on 
% this we can accurately hit some of our calibration targets directly.
% We do this for: median wage, pension, incomefloor
FnsToEvaluate.earnings=@(aprime,a,z,upsilon,epsilon,w,kappa_j,alpha,kappabeta,agej,Jr) w*(1-upsilon)*exp(kappa_j+alpha+kappabeta+z+epsilon)*(agej<Jr); % labor earnings
if preCalib==1
    % First, solve the model to get median wage
    % Then set w, incomefloor, and pension, based on this
    disp('Value Fn')
    [~, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_zupsilon,N_j,N_i, d_grid, a_grid, zupsilon_grid_J, pi_zupsilon_J, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
    disp('Agent Dist')
    StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_zupsilon,N_j,N_i,pi_zupsilon_J,Params,simoptions);

    disp('All Stats')
    % We want the median earnings, conditional on being employed
    simoptions.conditionalrestrictions.employed=@(aprime,a,z,upsilon,epsilon,agej,Jr) (1-upsilon)*(agej<Jr); % 'employed' (non-employment upsilon is equal to zero; and working age)
    AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params,n_d,n_a,n_zupsilon,N_j,N_i,d_grid, a_grid, zupsilon_grid_J, simoptions);
    simoptions=rmfield(simoptions,'conditionalrestrictions');

    % Set the wage to be equal to the median for US in 2010 (GKOS2021 use 2010 as their reference year, to which they convert all nominal amounts using the PCE index)
    % $26,363.55 (which is the median earnings for US in 2010, the year for which the dollars in GKOS2021 estimates are for: https://www.ssa.gov/cgi-bin/netcomp.cgi?year=2010)
    targetmedianwage=26363.55;

    Params.w=targetmedianwage/AllStats.employed.earnings.Median; % Note, this median is with the conditional restriction, so is conditional on being employed
    
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

    Params.incomefloor=targetincomefloor;
    Params.pension=targetpension;
    
    % Note: asset grids were just silly prior to now, but since this is an
    % exogenous earnings model that made no difference to the earnings
    % which is the only thing we used until now.

    % clean up large objects
    clear Policy StationaryDist

    % I also want to keep
    PreCalib_MeanEarnings=Params.w*AllStats.earnings.Mean;

    save(['./SavedOutput/Main/PreCalib',num2str(useModel),'.mat'],'Params')
else
    load(['./SavedOutput/Main/PreCalib',num2str(useModel),'.mat'],'Params')    
end


%% Calibrate the life-cycle model based on calibration targets (part 2)
% Set up FnsToEvaluate
FnsToEvaluate.bequest=@(aprime,a,z,upsilon,epsilon,agej,Jbeq,sj) (1-sj)*aprime*(agej>=Jbeq); % value of bequests
FnsToEvaluate.bequestamount=@(aprime,a,z,upsilon,epsilon,agej,Jbeq) aprime*(agej>=Jbeq); % value of bequests (conditional on age and dying)
FnsToEvaluate.assets=@(aprime,a,z,upsilon,epsilon) a; % a is the current asset holdings
if doCalib==1
    % calibrate beta to target aggregate asset-income ratio
    % calibrate two warm-glow parameters for three targets
    % put weight of 2 on beta target, weights of 1 on the three warm-glow targets 
    caliboptions.verbose=1;

    % Parameters to calibrate
    CalibParamNames={'beta','warmglow1','warmglow2'};
    % all three of these must be positive valued
    caliboptions.constrainpositive={'beta','warmglow1','warmglow2'};

    % Setup for calibration targets
    calibtarget_assettoearningsratio=(3/2)*3.9; % =assets to earnings ratio = capital to income ratio of 3.9, together with earnings is 2/3rd income
    % 3.9 capital ratio: from https://cepr.org/voxeu/columns/us-capital-glut-and-other-myths
    
    % Now set up these targets for the CalibrateLifeCycleModel_PType() command
    TargetMoments.AllStats.assets.Mean=calibtarget_assettoearningsratio*PreCalib_MeanEarnings; % target for beta

    % Only two warm-glow parameters (the third was just set equal to have same curvature as utility)

    % Nishiyama (2002) uses “The calibration in this paper thus uses the annual flow of bequests, 1.00% of net worth, estimated by Gale and Scholz (1994). (See Table III.)”
    % We know total earnings (as it is exogenous), and we know that total
    % asset is targeted to be X times earnings. So we just need total
    % bequests to be X/100 times earnings.
    % Note: we compute the total (weighted) bequests, then divide it by the
    % probabilities of death at each age, and this is what targets X/100 times earnings
    TargetMoments.AllStats.bequest.Mean=(1/100)*calibtarget_assettoearningsratio*PreCalib_MeanEarnings;

    % “The mean estate value is $94,469, but the median is half as much, $50,000.” (Hurd & Smith, 1999)
    % All the above are 1994 dollars, so need to multiply by 1.485 to get 2010 dollars [https://fred.stlouisfed.org/series/CPIAUCSL#0]

    % Rather than actually target this, I instead target that the Mean-to-Median ratio averages 94,469/50,000 for assets during the
    % bequest ages (this is not the same things as it ignores that each of these ages should be weighted by the probability of dying as that age; but is close enough).
    simoptions.conditionalrestrictions.bequestages=@(aprime,a,z,upsilon,epsilon,agej,Jbeq) (agej>=Jbeq); % ages at which bequests are left
    TargetMoments.AllStats.bequestages.bequestamount.RatioMeanToMedian=94469/50000; % Note, can just leave them in 1994 dollars, as anyway just using ratio
    % Because of how median works, this needs to use bequestamount, not bequest
    
    ParametrizePTypeFn=[]; % not needed

    if useModel==6
        % To speed things up (this model timed out as server has 24hr job limit, rest were fast enough)
        Params.beta=1;
        Params.warmglow1=8500;
        Params.warmglow2=0.001;
    end
    
    % Perform the calibration
    caliboptions.metric='sum_logratiosquared';
    [CalibParams,calibsummary]=CalibrateLifeCycleModel_PType(CalibParamNames,TargetMoments,n_d,n_a,n_zupsilon,N_j,N_i,d_grid, a_grid, zupsilon_grid_J, pi_zupsilon_J, ReturnFn, Params, DiscountFactorParamNames, jequaloneDist,AgeWeightsParamNames, PTypeDistParamNames, ParametrizePTypeFn, FnsToEvaluate, caliboptions, vfoptions,simoptions);
    % Store the calibration in Params
    for pp=1:length(CalibParamNames)
        Params.(CalibParamNames{pp})=CalibParams.(CalibParamNames{pp});
    end
    % Calculate moments that were targetted in the calibration
    [~, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_zupsilon,N_j,N_i, d_grid, a_grid, zupsilon_grid_J, pi_zupsilon_J, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
    StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_zupsilon,N_j,N_i,pi_zupsilon_J,Params,simoptions);
    AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params,n_d,n_a,n_zupsilon,N_j,N_i,d_grid, a_grid, zupsilon_grid_J, simoptions);

    % Calibrated moments
    calibresults=nan(2,3);
    calibresults(1,1)=AllStats.assets.Mean;
    calibresults(2,1)=TargetMoments.AllStats.assets.Mean;
    calibresults(1,2)=AllStats.bequest.Mean;
    calibresults(2,2)=TargetMoments.AllStats.bequest.Mean;
    calibresults(1,3)=AllStats.bequestages.bequestamount.RatioMeanToMedian;
    calibresults(2,3)=TargetMoments.AllStats.bequestages.bequestamount.RatioMeanToMedian;
    
    % Delete some things we no longer need
    simoptions=rmfield(simoptions,'conditionalrestrictions');
    
    save(['./SavedOutput/Main/Calib',num2str(useModel),'.mat'],'Params','CalibParams','calibsummary','calibresults','FnsToEvaluate')
else
    load(['./SavedOutput/Main/Calib',num2str(useModel),'.mat'],'Params','CalibParams','calibsummary','calibresults','FnsToEvaluate')
end

%% Model is now properly calibrated
fprintf('Finished Calibration! \n')

%% Now solve the value function iteration problem and stationary dist
% Same some minor things
save(['./SavedOutput/Main/BasicOutput',num2str(useModel),'.mat'], 'n_d','n_a','n_z','n_upsilon','n_epsilon','n_zupsilon','n_alpha','n_kappabeta','N_i','N_j')
if SolveVandStationaryDist==1
    disp('Solving the Value Fn and Policy Fn')
    tic;
    [V, Policy]=ValueFnIter_Case1_FHorz_PType(n_d,n_a,n_zupsilon,N_j,N_i, d_grid, a_grid, zupsilon_grid_J, pi_zupsilon_J, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
    time1=toc
    disp('Solving the Agent Dist')
    tic;
    StationaryDist=StationaryDist_Case1_FHorz_PType(jequaloneDist,AgeWeightsParamNames,PTypeDistParamNames,Policy,n_d,n_a,n_zupsilon,N_j,N_i,pi_zupsilon_J,Params,simoptions);
    time2=toc

    Names_i=fieldnames(Policy);
    % Don't actually need V for anything, only the first period of it, so just save that and clear the rest to free up space
    V_age1=struct();
    for ii=1:N_i
        V_ii=V.(Names_i{ii});
        V_age1.(Names_i{ii})=V_ii(:,:,:,:,1);
    end
    save(['./SavedOutput/Main/Vage1Save',num2str(useModel),'.mat'],'V_age1','-v7.3')
    % save(['./SavedOutput/Main/VSave',num2str(useModel),'.mat'],'V','-v7.3')
    clear V V_ii
    
    save(['./SavedOutput/Main/PolicySave',num2str(useModel),'.mat'],'Policy','-v7.3')
    save(['./SavedOutput/Main/RestSave',num2str(useModel),'.mat'],'n_zupsilon','n_epsilon','zupsilon_grid_J','epsilon_grid','pi_zupsilon_J','pi_epsilon','Names_i')
    save(['./SavedOutput/Main/StationaryDist',num2str(useModel),'.mat'], 'StationaryDist','jequaloneDist','-v7.3')

    %% Test that Policy is not trying to leave the top of the grid
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
    clear temp tempdist % free up memory
else
    load(['./SavedOutput/Main/PolicySave',num2str(useModel),'.mat'])
    load(['./SavedOutput/Main/RestSave',num2str(useModel),'.mat'])
    load(['./SavedOutput/Main/StationaryDist',num2str(useModel),'.mat'])
end
% NOTE: I no longer keep a saved copy of StationaryDist because it uses too
% much harddrive space and is anyway easy to create. This does mean you
% CANNOT set SolveVandStationaryDist=0, but no big deal.


%% Now, we want to graph Life-Cycle Profiles


%% FnsToEvaluate are how we say what we want to simulate and graph
% Note: everything involving earnings needs to include a *(agej<Jr) term
% Like with return function, we have to include (aprime,a,z) as first inputs, then just any relevant parameters.
FnsToEvaluate.income=@(aprime,a,z,upsilon,epsilon,w,kappa_j,alpha,kappabeta,r,agej,Jr,incomefloor) max(w*(1-upsilon)*exp(kappa_j+alpha+kappabeta+z+epsilon)*(agej<Jr)+r*a,incomefloor); % labor earnings+ r*a
FnsToEvaluate.consumption=@(aprime,a,z,upsilon,epsilon,alpha,kappabeta,w,agej,Jr,pension,incomefloor,r,kappa_j,eta1,eta2) EarningsDynamics_ConsumptionFn(aprime,a,z,upsilon,epsilon,alpha,kappabeta,w,agej,Jr,pension,incomefloor,r,kappa_j,eta1,eta2); % consumption
FnsToEvaluate.z=@(aprime,a,z,upsilon,epsilon) z; % AR(1) with gaussian-mixture innovations
FnsToEvaluate.upsilon=@(aprime,a,z,upsilon,epsilon) upsilon; % non-employment shock
FnsToEvaluate.epsilon=@(aprime,a,z,upsilon,epsilon) epsilon; % transitiory shock
FnsToEvaluate.deterministicearnings=@(aprime,a,z,upsilon,epsilon,w,kappa_j,alpha,kappabeta,agej,Jr) w*exp(kappa_j+alpha+kappabeta)*(agej<Jr); % the part of earnings which is just deterministic (in terms of age and permanent type)
FnsToEvaluate.potentialearnings=@(aprime,a,z,upsilon,epsilon,w,kappa_j,alpha,kappabeta,agej,Jr) w*exp(kappa_j+alpha+kappabeta+z+epsilon)*(agej<Jr); % what earnings would be without non-employment shocks
FnsToEvaluate.disposableincome=@(aprime,a,z,upsilon,epsilon,alpha,kappabeta,w,agej,Jr,pension,incomefloor,r,kappa_j,eta1,eta2) EarningsDynamics_DisposableIncomeFn(aprime,a,z,upsilon,epsilon,alpha,kappabeta,w,agej,Jr,pension,incomefloor,r,kappa_j,eta1,eta2);

% The following four were all set above, just include them here for convience
% FnsToEvaluate.earnings=@(aprime,a,z,upsilon,epsilon,w,kappa_j,alpha,kappabeta,agej,Jr) w*(1-upsilon)*exp(kappa_j+alpha+kappabeta+z+epsilon)*(agej<Jr); % labor earnings
% FnsToEvaluate.assets=@(aprime,a,z,upsilon,epsilon) a; % a is the current asset holdings
% FnsToEvaluate.bequest=@(aprime,a,z,upsilon,epsilon,agej,Jbeq,sj) (1-sj)*aprime*(agej>=Jbeq); % value of bequests
% FnsToEvaluate.bequestamount=@(aprime,a,z,upsilon,epsilon,agej,Jbeq) aprime*(agej>=Jbeq); % value of bequests (conditional on age and dying)

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

    simoptions.lowmemory=1; % this thing keeps running out of memory, so use lowmemory
    simoptions.whichstats=ones(1,7); % need to use lower memory versions with such large age-groupings
    tic;
    AgeConditionalStats=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params,n_d,n_a,n_zupsilon,N_j,N_i,d_grid, a_grid, zupsilon_grid_J, simoptions);
    time3=toc
    simoptions=rmfield(simoptions,'whichstats');
    simoptions=rmfield(simoptions,'lowmemory');

    save(['./SavedOutput/Main/AgeConditionalStats',num2str(useModel),'.mat'], 'AgeConditionalStats','-v7.3')
else
    load(['./SavedOutput/Main/AgeConditionalStats',num2str(useModel),'.mat'])
end


%% Calculate the life-cycle profiles for the working-age population by using the simoptions.agegroupings
% This is directly comparable to GKOS2021 as their sample is just the working age
if CalculateStatistics(2)==1
    disp('Calculate the life-cycle profiles for age-groupings')

    simoptions.agegroupings=[1,Params.Jr]; % Will look at stats for age-group 1 to Params.Jr-1, which is everyone of working age (Params.Jr is first year of retirement)
    simoptions.whichstats=ones(1,7); % need to use lower memory versions with such large age-groupings
    simoptions.lowmemory=1; % this thing keeps running out of memory, so use lowmemory
    tic;
    AgeConditionalStats_Grouped_AgeGroupings=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params,n_d,n_a,n_zupsilon,N_j,N_i,d_grid, a_grid, zupsilon_grid_J, simoptions);
    time4=toc
    
    save(['./SavedOutput/Main/AgeConditionalStats_Grouped_AgeGroupings',num2str(useModel),'.mat'], 'AgeConditionalStats_Grouped_AgeGroupings','-v7.3')
    
    simoptions=rmfield(simoptions,'agegroupings');
    simoptions=rmfield(simoptions,'whichstats');
    simoptions=rmfield(simoptions,'lowmemory');
else
    load(['./SavedOutput/Main/AgeConditionalStats_Grouped_AgeGroupings',num2str(useModel),'.mat'])
end


%% Some stats on the entire agent dist
if CalculateStatistics(3)==1
    disp('Calculate AllStats on agent dist')

    tic;
    AllStats=EvalFnOnAgentDist_AllStats_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate, Params,n_d,n_a,n_zupsilon,N_j,N_i,d_grid, a_grid, zupsilon_grid_J, simoptions);
    time5=toc
    save(['./SavedOutput/Main/VariousStats',num2str(useModel),'.mat'], 'AllStats','-v7.3')

    % only use: earnings, income, assets, consumption
end


%% Create some individual household simulations
if CalculateStatistics(4)==1
    disp('Create some individual household simulations')
    simoptions.lowmemory=1;
    
    tic;
    SimPanelValues=SimPanelValues_FHorz_Case1_PType(jequaloneDist,PTypeDistParamNames,Policy,FnsToEvaluate,Params,n_d,n_a,n_zupsilon,N_j,N_i,d_grid,a_grid,zupsilon_grid_J,pi_zupsilon_J, simoptions);
    time6=toc
    
    tic;
    simoptions.npoints=5; % Report transition probabilities for quintiles
    RankTransitionProbabilities_Grouped=EvalPanelData_RankTransProbs(SimPanelValues,simoptions);
    time7=toc

    simoptions=rmfield(simoptions,'lowmemory');
    simoptions=rmfield(simoptions,'npoints');

    save(['./SavedOutput/Main/SimPanelValues',num2str(useModel),'.mat'], 'SimPanelValues','RankTransitionProbabilities_Grouped','-v7.3')
else
    load(['./SavedOutput/Main/SimPanelValues',num2str(useModel),'.mat'])
end


%% Consumption insurance calculations (based on the simulated panel data)
if CalculateStatistics(5)==1
    tic;
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

    time8=toc
end


%%
% Save various other things about the model
Names_i=fieldnames(Policy);
save(['./SavedOutput/Main/GeneralOutput',num2str(useModel),'.mat'], 'Params','vfoptions','simoptions','a_grid','jequaloneDist','PTypeDist','AgeWeightsParamNames','Names_i')


%% Solve deterministic version (without shocks) to use for welfare calculations.

% To compare to DeNardi, Fella & Paz-Pardo (2018) results on welfare, we need the value function at age 1.
% We also need to solve a version without shocks and keep the value function at age 1

% We turn off all the shocks, and replace the deterministic earnings
% profiles with the age-conditional mean (for each permanent type)
if CalculateStatistics(6)==1
    tic;
    disp('Welfare costs of shocks')
    % load(['./SavedOutput/Main/VSave',num2str(useModel),'.mat'])
    % We don't need V, just the period 1 part
    load(['./SavedOutput/Main/Vage1Save',num2str(useModel),'.mat'],'V_age1')


    if useModel==6
        PTypemassvec=Params.statdist_alpha_kappabeta;
    else
        PTypemassvec=Params.statdist_alpha;
    end

    time9a=toc
    
      
    %% Before doing the deterministic versions, solve the risky one lots of time for different possible CEV
    tic;
    % Loop over possible values of CEV and figure out which is appropriate
    CEV_vec=0:0.01:3; % Consider 0% to 300% (to the nearest percentage points)
    W_vec=zeros(length(CEV_vec),1);
    W_ii_vec=zeros(length(CEV_vec),N_i);
    for cc=1:length(CEV_vec)
        Params.CEV=CEV_vec(cc);
        W=0;

        V_cc=ValueFnFromPolicy_Case1_FHorz_PType(Policy,n_d,n_a,n_zupsilon,N_j,Names_i,d_grid,a_grid,zupsilon_grid_J, pi_zupsilon_J, ReturnFn, Params, DiscountFactorParamNames, vfoptions);
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
    AgeConditionalStats_CEVcheck=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist, Policy, FnsToEvaluate_CEVcheck, Params,n_d,n_a,n_zupsilon,N_j,N_i,d_grid, a_grid, zupsilon_grid_J, simoptions);
    
    time9b=toc


    %% Now solve the deterministic versions of the model
    tic;
    % This problem is simple enough we can just use the default options
    vfoptions=struct();
    simoptions=struct();
    
    % Turn off z, epsilon
    n_z=1;
    z_grid=0;
    pi_zupsilon=1;
    
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
        pi_zupsilonepsilon=pi_zupsilonepsilon_J;
        zupsilonepsilon_grid=zupsilonepsilon_grid.*ones(1,Params.J);
        
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
        
        % simoptions.parallel=2;
        jequaloneDist_deterministic3=reshape(sum(sum(jequaloneDist,4),2),[n_a,1,n_upsilon,1]);
        % Try iterating since now no longer have z and epsilon
        StationaryDist_deterministic3=StationaryDist_Case1_FHorz_PType(jequaloneDist_deterministic3,AgeWeightsParamNames,PTypeDistParamNames,Policy_deterministic3,n_d,n_a,n_zupsilonepsilon,N_j,N_i,pi_zupsilonepsilon,Params,simoptions);
        
        AgeConditionalStats_CEVcheck3=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_deterministic3, Policy_deterministic3, FnsToEvaluate_CEVcheck, Params,n_d,n_a,n_zupsilonepsilon,N_j,N_i,d_grid, a_grid, zupsilonepsilon_grid, simoptions);

        % Turn the options back off again
        vfoptions=struct();
        simoptions=struct();
    end
    
    time9c=toc
    
    %% Now, turn of z, epsilon and upsilon (will be CEV1)
    tic;
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

    jequaloneDist_deterministic=reshape(sum(sum(sum(jequaloneDist,4),3),2),[n_a,1]);
    % Try iterating since now no longer have z and epsilon
    StationaryDist_deterministic=StationaryDist_Case1_FHorz_PType(jequaloneDist_deterministic,AgeWeightsParamNames,PTypeDistParamNames,Policy_deterministic,n_d,n_a,n_zupsilonepsilon,N_j,N_i,pi_zupsilonepsilon,Params,simoptions);

    AgeConditionalStats_CEVcheck1=LifeCycleProfiles_FHorz_Case1_PType(StationaryDist_deterministic, Policy_deterministic, FnsToEvaluate_CEVcheck, Params,n_d,n_a,n_zupsilonepsilon,N_j,N_i,d_grid, a_grid, zupsilonepsilon_grid, simoptions);

    time9d=toc

    %% Now do a version eliminating z, epsilon, upsilon and permanent types (will be CEV2)
    tic;
    % Do an alternative version where we also eliminate permanent types (as from the perspective of the 'behind the veil' this is also an uninsurable risk).
    Params=rmfield(Params,'kappa_j');
    Params.kappa_j=log(AgeConditionalStats.earnings.Mean/Params.w);  % log() as in model it is exp(kappa_j+...)
    % Note: no need to remove statdist_alpha (or statdist_alpha_kappabeta) as it is anyway not used by following line
    % Note: following line is not PType
    [V_deterministic2, Policy_deterministic2]=ValueFnIter_Case1_FHorz(n_d,n_a,n_zupsilonepsilon,N_j, d_grid, a_grid, zupsilonepsilon_grid, pi_zupsilonepsilon, ReturnFn, Params, DiscountFactorParamNames, [], vfoptions);
    
    % simoptions.parallel=2;
    jequaloneDist_deterministic2=reshape(sum(sum(sum(jequaloneDist,4),3),2),[n_a,1]);
    % Try iterating since now no longer have z and epsilon
    StationaryDist_deterministic2=StationaryDist_FHorz_Case1(jequaloneDist_deterministic2,AgeWeightsParamNames,Policy_deterministic2,n_d,n_a,n_zupsilonepsilon,N_j,pi_zupsilonepsilon,Params,simoptions);

    AgeConditionalStats_CEVcheck2=LifeCycleProfiles_FHorz_Case1(StationaryDist_deterministic2, Policy_deterministic2, FnsToEvaluate_CEVcheck, Params,[],n_d,n_a,n_zupsilonepsilon,N_j,d_grid, a_grid, zupsilonepsilon_grid, simoptions);

    time9e=toc

    %% Do the welfare calculations
    tic;

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
        V_ii=V_age1.(Names_i{ii});
        % V_ii=V.(Names_i{ii});
        V_deterministic_ii=V_deterministic.(Names_i{ii});
        
        % We just want age 1
        % V_ii=V_ii(:,:,:,:,1);
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
    
    % if useModel==2
    %     save('./SavedOutput/DebugWelfareCEV.mat','CEV1','CEV2','CEV3','W','Wd','Wd2','Wd3','W_vec')
    %     save('./SavedOutput/DebugCEV3.mat','V_deterministic','V_deterministic2','V_deterministic3','W_ii','Wd_ii','Wd2_ii','Wd3_ii','-v7.3')
    %     save('./SavedOutput/DebugCEV4.mat','AgeConditionalStats_CEVcheck', 'AgeConditionalStats_CEVcheck1', 'AgeConditionalStats_CEVcheck2', 'AgeConditionalStats_CEVcheck3')
    %     save('./SavedOutput/DebugCEV6.mat','StationaryDist_deterministic', 'StationaryDist_deterministic2', 'StationaryDist_deterministic3')
    % elseif useModel==3
    %     save('./SavedOutput/DebugWelfareCEV.mat','CEV1','CEV2','W','Wd','Wd2','W_vec')
    %     save('./SavedOutput/DebugCEV3.mat','V_deterministic','V_deterministic2','W_ii','Wd_ii','Wd2_ii','-v7.3')
    %     save('./SavedOutput/DebugCEV4.mat','AgeConditionalStats_CEVcheck', 'AgeConditionalStats_CEVcheck1', 'AgeConditionalStats_CEVcheck2')
    %     save('./SavedOutput/DebugCEV6.mat','StationaryDist_deterministic', 'StationaryDist_deterministic2')
    % end
    
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

    time9f=toc

    time9=time9a+time9b+time9c+time9d+time9e+time9f
end


%% All the runtimes so I can see which parts take more/less time
fprintf('All the runtimes! \n')
if SolveVandStationaryDist==1
    fprintf('value fn:   %3.2f \n', time1)
    fprintf('agent dist: %3.2f \n', time2)
end
if CalculateStatistics(1)==1
    fprintf('Life-Cycle Profiles: %3.2f \n', time3)
end
if CalculateStatistics(2)==1
    fprintf('Life-Cycle Profiles, group ages: %3.2f \n', time4)
end
if CalculateStatistics(3)==1
    fprintf('AllStats: %3.2f \n', time5)
end
if CalculateStatistics(4)==1
    fprintf('sim panel: %3.2f \n', time6)
    fprintf('rank transitions: %3.2f \n', time7)
end
if CalculateStatistics(5)==1
    fprintf('consumption insurance: %3.2f \n', time8)
end
if CalculateStatistics(6)==1
    fprintf('welfare eval: %3.2f \n', time9)
    fprintf('         (breakdown: %3.2f, %3.2f, %3.2f, %3.2f, %3.2f, %3.2f \n', time9a, time9b, time9c, time9d, time9e, time9f)
end


%% Create Tables of the calibrated parametesr
load(['./SavedOutput/Main/Calib',num2str(useModel),'.mat'],'Params')

% Table of life-cycle model parameters (except the earnings process)
FID = fopen(['./SavedOutput/LatexInputs/Table_LifeCycleModelParameters_Model',num2str(useModel),'.tex'], 'w');
fprintf(FID, 'Parameters of the Life-Cycle Model \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcc} \n \\hline \\hline \n');
fprintf(FID, ' Total periods & $J$ & %8.2f \\\\ \n', Params.J);
fprintf(FID, ' Retirement period & $Jr$ & %8.2f \\\\ \n', Params.Jr);
fprintf(FID, ' Leave bequests after age & $Jbeq$ & %8.2f \\\\ \n', Params.Jbeq);
fprintf(FID, ' Age in years at j=1 & $ $ & %8.2f \\\\ \n', Params.agejshifter+1);
fprintf(FID, ' \\multicolumn{3}{l}{Preferences} \\\\ \n');
fprintf(FID, ' Discount factor & $\\beta$ & %8.2f \\\\ \n', Params.beta);
fprintf(FID, ' CRRA & $\\sigma$ & %8.2f \\\\ \n', Params.sigma);
fprintf(FID, ' \\multicolumn{3}{l}{Prices and policies} \\\\ \n');
fprintf(FID, ' Wage & $w$ & %8.2f \\\\ \n', Params.w);
fprintf(FID, ' Interest rate & $r$ & %8.2f \\\\ \n', Params.r);
fprintf(FID, ' Pension & $b$ & %8.2f \\\\ \n', Params.pension);
fprintf(FID, ' Income floor & $\\theta$ & %8.2f \\\\ \n', Params.incomefloor);
fprintf(FID, ' \\multicolumn{3}{l}{Warm glow of bequests} \\\\ \n');
fprintf(FID, ' relative importance of bequests & $\\phi_1$ & %8.1f \\\\ \n', Params.warmglow1);
fprintf(FID, ' bequests as a luxury good & $\\phi_2$ & %8.1f \\\\ \n', Params.warmglow2);
fprintf(FID, ' curvature for bequests & $\\phi_3$ & %8.1f \\\\ \n', Params.warmglow3);
fprintf(FID, '\\hline \n \\end{tabular*} \n');
% fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
% fprintf(FID, 'Note: notice. \n');
% fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%% Create a Table of the conditional survival probabilities (these are same across all models, so just overwrite the file each time it is run)
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


