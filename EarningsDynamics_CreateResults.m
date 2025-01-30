%% Create all the results for paper
% This code is to be run after 'EarningsDynamics.m' (which needs to be run
% 6 times, for useModel=1,2,3,4,5,6)
% It takes all the output generated and saved by that script and uses it to
% create the results reported in the paper.

CreateFiguresAndTables=1 % 0 means skip, 1 creates all the figures and tables from the results
figure_c=0;

useModel_vec=[1,2,3,4,5,6] %[1,2,3,4,5,6]

nSigmaz_vec=[2,2,2,0.3,2,2]; % number of standard deviations for the max and min points of grid used for z
nSigma_alpha=1; % number of standard deviations for the max and min points of grid used for alpha (and kappa_beta)

%% Bunch of results for each model
for useModel=useModel_vec
    % There is also useModel 21 and 22, these are just same as 2, except using extended Farmer-Toda to target 2 and 4 moments respectively
    
    useModel % Just to make it easy to see the progress
    
    nSigmaz=nSigmaz_vec(useModel);
    
    %% Dicretized earnings dynamics
    load(['./SavedOutput/Main/BasicOutput',num2str(useModel),'.mat']) %, 'n_d','n_a','n_z','n_upsilon','n_epsilon','n_zupsilon','n_alpha','n_kappabeta','N_i','N_j')
    load(['./SavedOutput/BasicDiscretization/DiscretizedEarnings',num2str(useModel),'nSigmaz',num2str(nSigmaz),'nSigma_alpha',num2str(nSigma_alpha),'.mat'])
    
    % Load everything about the model (that was created and saved by EarningsDynamics.m)
%     load(['./SavedOutput/Main/PolicySave',num2str(useModel),'.mat'],'Policy') % Not needed at this stage
%     load(['./SavedOutput/Main/VSave',num2str(useModel),'.mat'],'V') % Not needed for anything
%     load(['./SavedOutput/Main/StationaryDist',num2str(useModel),'.mat']) % Not needed at this stage

    load(['./SavedOutput/Main/RestSave',num2str(useModel),'.mat'])

    load(['./SavedOutput/Main/GeneralOutput',num2str(useModel),'.mat'], 'Params','vfoptions','simoptions','a_grid','jequaloneDist','PTypeDist','AgeWeightsParamNames','Names_i')
    
    load(['./SavedOutput/Main/AgeConditionalStats',num2str(useModel),'.mat'])
    load(['./SavedOutput/Main/AgeConditionalStats_Grouped_AgeGroupings',num2str(useModel),'.mat'])
    load(['./SavedOutput/Main/VariousStats',num2str(useModel),'.mat'])
    load(['./SavedOutput/Main/SimPanelValues',num2str(useModel),'.mat'])
    
    %% Plot figures
    Params.ageyears=Params.agej+Params.agejshifter;   
    
    % Names_i is loaded from GeneralOutput
    
    Names_i2=cell(2*N_i,1); % Used for some graphs
    for ii=1:N_i
        Names_i2{2*ii-1}=Names_i{ii};
        Names_i2{2*ii}='';
    end
    
    if CreateFiguresAndTables==1
        %% Some life-cycle profiles
        % Life-cycle profiles of the mean and variance of earnings by permanent type, and grouped
        figure_c=figure_c+1;
        figure(figure_c)
        subplot(3,1,1)
        hold on
        for ii=1:N_i
            plot(Params.ageyears,AgeConditionalStats.earnings.(Names_i{ii}).Mean)
        end
        hold off
        legend(Names_i)
        title('Life-Cycle Profiles of mean earnings')
        subplot(3,1,2)
        hold on
        for ii=1:N_i
            plot(Params.ageyears,AgeConditionalStats.earnings.(Names_i{ii}).Mean)
            plot(Params.ageyears,AgeConditionalStats.earnings.(Names_i{ii}).StdDeviation,'.')
        end
        hold off
        legend(Names_i2)
        title('Life-Cycle Profiles of mean earnings and standard deviation')
        subplot(3,1,3)
        hold on
        plot(Params.ageyears,AgeConditionalStats.earnings.Mean)
        plot(Params.ageyears,AgeConditionalStats.earnings.StdDeviation,'.')
        hold off
        title('Life-Cycle Profiles of mean earnings and standard deviation for whole population')
        saveas(figure_c,['./SavedOutput/Graphs/EarningsProfiles_Model',num2str(useModel),'.png'])
        
        % Life-cycle profiles of the mean and variance of assets by permanent type
        figure_c=figure_c+1;
        figure(figure_c)
        subplot(3,1,1)
        hold on
        for ii=1:N_i
            plot(Params.ageyears,AgeConditionalStats.assets.(Names_i{ii}).Mean)
        end
        hold off
        legend(Names_i)
        title('Life-Cycle Profiles of mean assets')
        subplot(3,1,2)
        hold on
        for ii=1:N_i
            plot(Params.ageyears,AgeConditionalStats.assets.(Names_i{ii}).Mean)
            plot(Params.ageyears,AgeConditionalStats.assets.(Names_i{ii}).StdDeviation,'.')
        end
        hold off
        legend(Names_i2)
        title('Life-Cycle Profiles of mean assets and standard deviation')
        subplot(3,1,3)
        hold on
        plot(Params.ageyears,AgeConditionalStats.assets.Mean)
        plot(Params.ageyears,AgeConditionalStats.assets.StdDev,'.')
        hold off
        title('Life-Cycle Profiles of mean assets and standard deviation for whole population')
        saveas(figure_c,['./SavedOutput/Graphs/AssetsProfiles_Model',num2str(useModel),'.png'])
        
        % Life-cycle profiles of the mean and variance of consumption by permanent type
        figure_c=figure_c+1;
        figure(figure_c)
        subplot(3,1,1)
        hold on
        for ii=1:N_i
            plot(Params.ageyears,AgeConditionalStats.consumption.(Names_i{ii}).Mean)
        end
        hold off
        legend(Names_i)
        title('Life-Cycle Profiles of mean consumption')
        subplot(3,1,2)
        hold on
        for ii=1:N_i
            plot(Params.ageyears,AgeConditionalStats.consumption.(Names_i{ii}).Mean)
            plot(Params.ageyears,AgeConditionalStats.consumption.(Names_i{ii}).StdDeviation,'.')
        end
        hold off
        legend(Names_i2)
        title('Life-Cycle Profiles of mean consumption and standard deviation')
        subplot(3,1,3)
        hold on
        plot(Params.ageyears,AgeConditionalStats.consumption.Mean)
        plot(Params.ageyears,AgeConditionalStats.consumption.StdDev,'.')
        hold off
        title('Life-Cycle Profiles of mean consumption and standard deviation for whole population')
        saveas(figure_c,['./SavedOutput/Graphs/ConsumptionProfiles_Model',num2str(useModel),'.png'])
        
        % Life-cycle profiles of the means of earnings, assets and consumption by permanent type, and grouped
        figure_c=figure_c+1;
        figure(figure_c)
        subplot(3,1,1)
        plot(Params.ageyears,AgeConditionalStats.consumption.Mean)
        hold on
        for ii=1:N_i
            plot(Params.ageyears,AgeConditionalStats.consumption.(Names_i{ii}).Mean)
        end
        hold off
        legend({'Grouped',Names_i{:}})
        title('Life-Cycle Profiles of mean consumption')
        subplot(3,1,2)
        plot(Params.ageyears,AgeConditionalStats.earnings.Mean)
        hold on
        for ii=1:N_i
            plot(Params.ageyears,AgeConditionalStats.earnings.(Names_i{ii}).Mean)
        end
        hold off
        legend({'Grouped',Names_i{:}})
        title('Life-Cycle Profiles of mean earnings')
        subplot(3,1,3)
        plot(Params.ageyears,AgeConditionalStats.assets.Mean)
        hold on
        for ii=1:N_i
            plot(Params.ageyears,AgeConditionalStats.assets.(Names_i{ii}).Mean)
        end
        hold off
        legend({'Grouped',Names_i{:}})
        title('Life-Cycle Profiles of mean assets')
        saveas(figure_c,['./SavedOutput/Graphs/MeanLifeCycleProfiles_Model',num2str(useModel),'.png'])
        
        
        %% Some individual simulations
        
        figure_c=figure_c+1;
        figure(figure_c)
        indexestoplot=round(simoptions.numbersims*rand(1,10));% 10 households chosen at random
        subplot(4,2,1); plot(1:1:N_j, SimPanelValues.earnings(:,indexestoplot)) % Just plot the first 10
        title('Earnings')
        subplot(4,2,2); plot(1:1:N_j, SimPanelValues.assets(:,indexestoplot))
        title('Assets')
        subplot(4,2,3); plot(1:1:N_j, SimPanelValues.consumption(:,indexestoplot))
        title('Consumption')
        subplot(4,2,4); plot(1:1:N_j, SimPanelValues.z(:,indexestoplot))
        title('z')
        subplot(4,2,5); plot(1:1:N_j, SimPanelValues.upsilon(:,indexestoplot))
        title('upsilon (non-employment shock)')
        subplot(4,2,6); plot(1:1:N_j, SimPanelValues.epsilon(:,indexestoplot))
        title('epsilon')
        subplot(4,2,7); plot(1:1:N_j, SimPanelValues.deterministicearnings(:,indexestoplot))
        title('Deterministic Earnings')
        subplot(4,2,8); plot(1:1:N_j, SimPanelValues.potentialearnings(:,indexestoplot))
        title('Potential Earnings')
        saveas(figure_c,['./SavedOutput/Graphs/SimIndividuals_Model',num2str(useModel),'.png'])
        
        clear SimPanelValues
    end
    
    %% Miscellaneous things
    
    % To compare to US inequality statistics we want to keep
    InequalityStats_TableRaw.(['model',num2str(useModel)])=AllStats;
    
    % To compare to US working age inequality statistics we want to keep
    WorkingAgeInequalityStats_TableRaw.(['model',num2str(useModel)])=AgeConditionalStats_Grouped_AgeGroupings;
    
    % To compare to the Geweke & Keane (2000) results on the persistence of poverty we want to keep
    RankTransitionProbablities_TableRaw.(['model',num2str(useModel)])=RankTransitionProbabilities_Grouped;
    
    ParamsForModel.(['model',num2str(useModel)])=Params;
    
    
    if useModel==6
        % Create a plot of the heterogenous income profiles
        figure_c=figure_c+1;
        figure(figure_c)
        plot(Params.kappa_j'+Params.kappabeta')
        title('Model 6: Heterogenous income profiles, exp(kappa+kappabeta)')
        saveas(figure_c,['./SavedOutput/Graphs/EarningsDynamics_Model',num2str(useModel),'_HeteroIncomeProfiles.png'])
    end
    
end
save ./SavedOutput/Main/InequalityStats.mat InequalityStats_TableRaw WorkingAgeInequalityStats_TableRaw
save ./SavedOutput//Main/RankTransitionProbablities.mat RankTransitionProbablities_TableRaw
save ./SavedOutput/Main/ParamsForModel.mat ParamsForModel



%% Plot of agent distribution
for useModel=useModel_vec
    if CreateFiguresAndTables==1
        load(['./SavedOutput/Main/StationaryDist',num2str(useModel),'.mat'])
        load(['./SavedOutput/Main/BasicOutput',num2str(useModel),'.mat'],'N_i')
        load(['./SavedOutput/Main/GeneralOutput',num2str(useModel),'.mat'],'Names_i','a_grid')
        
        % Graph the pdf and cdf of assets
        % Note, looking at the following graph is also a way to check for people running in to top of grid
        figure_c=figure_c+1;
        figure(figure_c)
        subplot(2,1,1);
        semilogx(a_grid,sum(sum(sum(sum(StationaryDist.(Names_i{1}),2),3),4),5) ) % have to do this one before hold on so as to set the x-axis as log-spacing
        hold on
        for ii=2:N_i
            plot(a_grid,sum(sum(sum(sum(StationaryDist.(Names_i{ii}),2),3),4),5) ) % sums are eliminating all the dimensions other than assets
        end
        hold off
        title('pdf of assets')
        subplot(2,1,2);
        semilogx(a_grid,cumsum(sum(sum(sum(sum(StationaryDist.(Names_i{1}),2),3),4),5)) )
        hold on
        for ii=2:N_i
            plot(a_grid,cumsum(sum(sum(sum(sum(StationaryDist.(Names_i{ii}),2),3),4),5)) )
        end
        hold off
        title('cdf of assets')
        saveas(figure_c,['./SavedOutput/Graphs/Assetpdfcdf_Model',num2str(useModel),'.png'])
        
        clear StationaryDist
    end
end


%% A check that top of grid is not being reached
% Idea is to make sure that there are zero agents near the top of the asset grid
% Should print out lots of zeros :)
for useModel=useModel_vec
    load(['./SavedOutput/Main/TestPolicy',num2str(useModel),'.mat'],'test_PolicyLeavingTopOfGrid')

    load(['./SavedOutput/Main/GeneralOutput',num2str(useModel),'.mat'],'Names_i')
    N_i=length(Names_i);
    
    % Check for leaving the top of asset grid
    fprintf('For model %i \n',useModel)
    fprintf(['Following are to check that max value of grid on assets is large enough for model number ',num2str(useModel),': \n'])
    for ii=1:N_i
        fprintf('Number of grid points choosing max assets for PType %i is %i \n',ii, test_PolicyLeavingTopOfGrid(ii,1) )
        fprintf('  Mass of agents on top 50 asset points for PType %i is %8.12f \n', ii, test_PolicyLeavingTopOfGrid(ii,2) )
        fprintf('  Mass of agents on top 100 asset points for PType %i is %8.12f \n', ii, test_PolicyLeavingTopOfGrid(ii,3) )
        fprintf('  Mass of agents on top 200 asset points for PType %i is %8.12f \n', ii, test_PolicyLeavingTopOfGrid(ii,4) )
        fprintf('  Mass of agents on top 300 asset points for PType %i is %8.12f \n', ii, test_PolicyLeavingTopOfGrid(ii,5) )
    end
    fprintf('\n')

    
end




%% Create table about inequality statistics
load ./SavedOutput/Main/InequalityStats.mat InequalityStats_TableRaw WorkingAgeInequalityStats_TableRaw

EarningsInequality=zeros(9,6);
IncomeInequality=zeros(9,6);
AssetsInequality=zeros(9,6);
ConsumptionInequality=zeros(9,6);
for useModel=useModel_vec
    LorenzCurve_Earnings=InequalityStats_TableRaw.(['model',num2str(useModel)]).earnings.LorenzCurve;
    EarningsInequality(1,useModel)=Gini_from_LorenzCurve(LorenzCurve_Earnings); % Gini Coefficient
    EarningsInequality(2:6,useModel)=[LorenzCurve_Earnings(20); LorenzCurve_Earnings(40)-LorenzCurve_Earnings(20); LorenzCurve_Earnings(60)-LorenzCurve_Earnings(40); LorenzCurve_Earnings(80)-LorenzCurve_Earnings(60); LorenzCurve_Earnings(100)-LorenzCurve_Earnings(80)]; % Quintiles
    EarningsInequality(7,useModel)=LorenzCurve_Earnings(95)-LorenzCurve_Earnings(90); % 90-95
    EarningsInequality(8,useModel)=LorenzCurve_Earnings(99)-LorenzCurve_Earnings(95); % 95-99
    EarningsInequality(9,useModel)=LorenzCurve_Earnings(100)-LorenzCurve_Earnings(99); % 99-100

    LorenzCurve_Income=InequalityStats_TableRaw.(['model',num2str(useModel)]).income.LorenzCurve;
    IncomeInequality(1,useModel)=Gini_from_LorenzCurve(LorenzCurve_Income); % Gini Coefficient
    IncomeInequality(2:6,useModel)=[LorenzCurve_Income(20); LorenzCurve_Income(40)-LorenzCurve_Income(20); LorenzCurve_Income(60)-LorenzCurve_Income(40); LorenzCurve_Income(80)-LorenzCurve_Income(60); LorenzCurve_Income(100)-LorenzCurve_Income(80)]; % Quintiles
    IncomeInequality(7,useModel)=LorenzCurve_Income(95)-LorenzCurve_Income(90); % 90-95
    IncomeInequality(8,useModel)=LorenzCurve_Income(99)-LorenzCurve_Income(95); % 95-99
    IncomeInequality(9,useModel)=LorenzCurve_Income(100)-LorenzCurve_Income(99); % 99-100

    LorenzCurve_Assets=InequalityStats_TableRaw.(['model',num2str(useModel)]).assets.LorenzCurve;
    AssetsInequality(1,useModel)=Gini_from_LorenzCurve(LorenzCurve_Assets); % Gini Coefficient
    AssetsInequality(2:6,useModel)=[LorenzCurve_Assets(20); LorenzCurve_Assets(40)-LorenzCurve_Assets(20); LorenzCurve_Assets(60)-LorenzCurve_Assets(40); LorenzCurve_Assets(80)-LorenzCurve_Assets(60); LorenzCurve_Assets(100)-LorenzCurve_Assets(80)]; % Quintiles
    AssetsInequality(7,useModel)=LorenzCurve_Assets(95)-LorenzCurve_Assets(90); % 90-95
    AssetsInequality(8,useModel)=LorenzCurve_Assets(99)-LorenzCurve_Assets(95); % 95-99
    AssetsInequality(9,useModel)=LorenzCurve_Assets(100)-LorenzCurve_Assets(99); % 99-100

    LorenzCurve_Consumption=InequalityStats_TableRaw.(['model',num2str(useModel)]).consumption.LorenzCurve;
    ConsumptionInequality(1,useModel)=Gini_from_LorenzCurve(LorenzCurve_Consumption); % Gini Coefficient
    ConsumptionInequality(2:6,useModel)=[LorenzCurve_Consumption(20); LorenzCurve_Consumption(40)-LorenzCurve_Consumption(20); LorenzCurve_Consumption(60)-LorenzCurve_Consumption(40); LorenzCurve_Consumption(80)-LorenzCurve_Consumption(60); LorenzCurve_Consumption(100)-LorenzCurve_Consumption(80)]; % Quintiles
    ConsumptionInequality(7,useModel)=LorenzCurve_Consumption(95)-LorenzCurve_Consumption(90); % 90-95
    ConsumptionInequality(8,useModel)=LorenzCurve_Consumption(99)-LorenzCurve_Consumption(95); % 95-99
    ConsumptionInequality(9,useModel)=LorenzCurve_Consumption(100)-LorenzCurve_Consumption(99); % 99-100
end

% Change the units (5% as 5 rather than 0.05)
EarningsInequality(2:end,:)=100*EarningsInequality(2:end,:); % (1,:) is the Gini coeff
IncomeInequality(2:end,:)=100*IncomeInequality(2:end,:); % (1,:) is the Gini coeff
AssetsInequality(2:end,:)=100*AssetsInequality(2:end,:); % (1,:) is the Gini coeff
ConsumptionInequality(2:end,:)=100*ConsumptionInequality(2:end,:); % (1,:) is the Gini coeff

%Table: Inequality Version 1: By concept
FID = fopen('./SavedOutput/LatexInputs/EarningsDynamics_InequalityV1.tex', 'w');
fprintf(FID, 'Inequality: Share of X held by quintile/top percent \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccccccccc} \n \\hline \\hline \n');
fprintf(FID, '& & & & & & & \\multicolumn{3}{c}{Top Groups} \\\\ \\cline{8-10} \n');
fprintf(FID, '& & \\multicolumn{5}{c}{Quintile} & \\multicolumn{3}{c}{Percentile} \\\\ \\cline{3-7} \\cline{8-10} \n');
fprintf(FID, '              & Gini  & First  & Second & Third & Fourth & Fifth & 90th-95th & 95th-99th & 99th-100th \\\\ \n \\hline \n');
fprintf(FID, '\\multicolumn{10}{l}{\textbf{Earnings:}} \\\\ \n');
fprintf(FID, 'US Data       & 0.67  & -0.1  & 3.0   & 10.4  & 20.2  & 66.5  & 12.4  & 18.4  & 18.8  \\\\ \n');
fprintf(FID, 'Model 1       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,1)');
fprintf(FID, 'Model 2       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,2)');
fprintf(FID, 'Model 3       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,3)');
fprintf(FID, 'Model 4       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,4)');
fprintf(FID, 'Model 5       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,5)');
fprintf(FID, 'Model 6       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,6)');
fprintf(FID, '\\multicolumn{10}{l}{Income:} \\\\ \n');
fprintf(FID, 'US Data       & 0.58  & 3.0   & 6.5   & 10.9  & 18.1  & 61.4  & 10.8  & 16.5  & 19.7  \\\\ \n');
fprintf(FID, 'Model 1       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,1)');
fprintf(FID, 'Model 2       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,2)');
fprintf(FID, 'Model 3       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,3)');
fprintf(FID, 'Model 4       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,4)');
fprintf(FID, 'Model 5       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,5)');
fprintf(FID, 'Model 6       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,6)');
fprintf(FID, '\\multicolumn{10}{l}{Wealth:} \\\\ \n');
fprintf(FID, 'US Data       & 0.85  & -0.7  & 0.6   & 3.2   & 9.8   & 87.0  & 12.1  & 27.4  & 35.5  \\\\ \n');
fprintf(FID, 'Model 1       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,1)');
fprintf(FID, 'Model 2       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,2)');
fprintf(FID, 'Model 3       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,3)');
fprintf(FID, 'Model 4       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,4)');
fprintf(FID, 'Model 5       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,5)');
fprintf(FID, 'Model 6       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,6)');
fprintf(FID, '\\multicolumn{10}{l}{Consumption:} \\\\ \n');
fprintf(FID, 'US Data       & 0.32  & --  & --   & --   & --   & --  & --  & --  & --  \\\\ \n');
fprintf(FID, 'Model 1       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,1)');
fprintf(FID, 'Model 2       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,2)');
fprintf(FID, 'Model 3       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,3)');
fprintf(FID, 'Model 4       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,4)');
fprintf(FID, 'Model 5       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,5)');
fprintf(FID, 'Model 6       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,6)');
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: Income equals earnings plus interest income on assets. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%Table: Inequality Version 2: By Model
FID = fopen('./SavedOutput/LatexInputs/EarningsDynamics_InequalityV2.tex', 'w');
fprintf(FID, 'Inequality:  Share of X held by quintile/top percent \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccccccccc} \n \\hline \\hline \n');
fprintf(FID, '& & & & & & & \\multicolumn{3}{c}{Top Groups} \\\\ \\cline{8-10} \n');
fprintf(FID, '& & \\multicolumn{5}{c}{Quintile} & \\multicolumn{3}{c}{(Percentile)} \\\\ \\cline{3-7} \\cline{8-10} \n');
fprintf(FID, '              & Gini  & First  & Second & Third & Fourth & Fifth & 90th-95th & 95th-99th & 99th-100th \\\\ \n \\hline \n');
fprintf(FID, '\\multicolumn{10}{l}{Model 1:} \\\\ \n');
fprintf(FID, 'Earnings      & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,1)');
fprintf(FID, 'Income        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,1)');
fprintf(FID, 'Wealth        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,1)');
fprintf(FID, 'Consumption   & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,1)');
fprintf(FID, '\\multicolumn{10}{l}{Model 2:} \\\\ \n');
fprintf(FID, 'Earnings      & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,2)');
fprintf(FID, 'Income        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,2)');
fprintf(FID, 'Wealth        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,2)');
fprintf(FID, 'Consumption   & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,2)');
fprintf(FID, '\\multicolumn{10}{l}{Model 3:} \\\\ \n');
fprintf(FID, 'Earnings      & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,3)');
fprintf(FID, 'Income        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,3)');
fprintf(FID, 'Wealth        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,3)');
fprintf(FID, 'Consumption   & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,3)');
fprintf(FID, '\\multicolumn{10}{l}{Model 4:} \\\\ \n');
fprintf(FID, 'Earnings      & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,4)');
fprintf(FID, 'Income        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,4)');
fprintf(FID, 'Wealth        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,4)');
fprintf(FID, 'Consumption   & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,4)');
fprintf(FID, '\\multicolumn{10}{l}{Model 5:} \\\\ \n');
fprintf(FID, 'Earnings      & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,5)');
fprintf(FID, 'Income        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,5)');
fprintf(FID, 'Wealth        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,5)');
fprintf(FID, 'Consumption   & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,5)');
fprintf(FID, '\\multicolumn{10}{l}{Model 6:} \\\\ \n');
fprintf(FID, 'Earnings      & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,6)');
fprintf(FID, 'Income        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,6)');
fprintf(FID, 'Wealth        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,6)');
fprintf(FID, 'Consumption   & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,6)');
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: Income equals earnings plus interest income on assets. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%% Inequality again, but this time just working age, not whole population.

EarningsInequality=zeros(9,6);
IncomeInequality=zeros(9,6);
AssetsInequality=zeros(9,6);
ConsumptionInequality=zeros(9,6);
for useModel=useModel_vec
    LorenzCurve_Earnings=WorkingAgeInequalityStats_TableRaw.(['model',num2str(useModel)]).earnings.LorenzCurve(:,1); % The 1 is working age
    EarningsInequality(1,useModel)=Gini_from_LorenzCurve(LorenzCurve_Earnings); % Gini Coefficient
    EarningsInequality(2:6,useModel)=[LorenzCurve_Earnings(20); LorenzCurve_Earnings(40)-LorenzCurve_Earnings(20); LorenzCurve_Earnings(60)-LorenzCurve_Earnings(40); LorenzCurve_Earnings(80)-LorenzCurve_Earnings(60); LorenzCurve_Earnings(100)-LorenzCurve_Earnings(80)]; % Quintiles
    EarningsInequality(7,useModel)=LorenzCurve_Earnings(95)-LorenzCurve_Earnings(90); % 90-95
    EarningsInequality(8,useModel)=LorenzCurve_Earnings(99)-LorenzCurve_Earnings(95); % 95-99
    EarningsInequality(9,useModel)=LorenzCurve_Earnings(100)-LorenzCurve_Earnings(99); % 99-100

    LorenzCurve_Income=WorkingAgeInequalityStats_TableRaw.(['model',num2str(useModel)]).income.LorenzCurve(:,1);
    IncomeInequality(1,useModel)=Gini_from_LorenzCurve(LorenzCurve_Income); % Gini Coefficient
    IncomeInequality(2:6,useModel)=[LorenzCurve_Income(20); LorenzCurve_Income(40)-LorenzCurve_Income(20); LorenzCurve_Income(60)-LorenzCurve_Income(40); LorenzCurve_Income(80)-LorenzCurve_Income(60); LorenzCurve_Income(100)-LorenzCurve_Income(80)]; % Quintiles
    IncomeInequality(7,useModel)=LorenzCurve_Income(95)-LorenzCurve_Income(90); % 90-95
    IncomeInequality(8,useModel)=LorenzCurve_Income(99)-LorenzCurve_Income(95); % 95-99
    IncomeInequality(9,useModel)=LorenzCurve_Income(100)-LorenzCurve_Income(99); % 99-100

    LorenzCurve_Assets=WorkingAgeInequalityStats_TableRaw.(['model',num2str(useModel)]).assets.LorenzCurve(:,1);
    AssetsInequality(1,useModel)=Gini_from_LorenzCurve(LorenzCurve_Assets); % Gini Coefficient
    AssetsInequality(2:6,useModel)=[LorenzCurve_Assets(20); LorenzCurve_Assets(40)-LorenzCurve_Assets(20); LorenzCurve_Assets(60)-LorenzCurve_Assets(40); LorenzCurve_Assets(80)-LorenzCurve_Assets(60); LorenzCurve_Assets(100)-LorenzCurve_Assets(80)]; % Quintiles
    AssetsInequality(7,useModel)=LorenzCurve_Assets(95)-LorenzCurve_Assets(90); % 90-95
    AssetsInequality(8,useModel)=LorenzCurve_Assets(99)-LorenzCurve_Assets(95); % 95-99
    AssetsInequality(9,useModel)=LorenzCurve_Assets(100)-LorenzCurve_Assets(99); % 99-100

    LorenzCurve_Consumption=WorkingAgeInequalityStats_TableRaw.(['model',num2str(useModel)]).consumption.LorenzCurve(:,1);
    ConsumptionInequality(1,useModel)=Gini_from_LorenzCurve(LorenzCurve_Consumption); % Gini Coefficient
    ConsumptionInequality(2:6,useModel)=[LorenzCurve_Consumption(20); LorenzCurve_Consumption(40)-LorenzCurve_Consumption(20); LorenzCurve_Consumption(60)-LorenzCurve_Consumption(40); LorenzCurve_Consumption(80)-LorenzCurve_Consumption(60); LorenzCurve_Consumption(100)-LorenzCurve_Consumption(80)]; % Quintiles
    ConsumptionInequality(7,useModel)=LorenzCurve_Consumption(95)-LorenzCurve_Consumption(90); % 90-95
    ConsumptionInequality(8,useModel)=LorenzCurve_Consumption(99)-LorenzCurve_Consumption(95); % 95-99
    ConsumptionInequality(9,useModel)=LorenzCurve_Consumption(100)-LorenzCurve_Consumption(99); % 99-100
end

% Change the units (5% as 5 rather than 0.05)
EarningsInequality(2:end,:)=100*EarningsInequality(2:end,:); % (1,:) is the Gini coeff
IncomeInequality(2:end,:)=100*IncomeInequality(2:end,:); % (1,:) is the Gini coeff
AssetsInequality(2:end,:)=100*AssetsInequality(2:end,:); % (1,:) is the Gini coeff
ConsumptionInequality(2:end,:)=100*ConsumptionInequality(2:end,:); % (1,:) is the Gini coeff

%Table: Working Age Inequality Version 1: By concept
FID = fopen('./SavedOutput/LatexInputs/EarningsDynamics_WorkingAgeInequalityV1.tex', 'w');
fprintf(FID, 'Working Age Inequality:  Share of X held by quintile/top percent \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccccccccc} \n \\hline \\hline \n');
fprintf(FID, '& & & & & & & \\multicolumn{3}{c}{Top Groups} \\\\ \\cline{8-10} \n');
fprintf(FID, '& & \\multicolumn{5}{c}{Quintile} & \\multicolumn{3}{c}{(Percentile)} \\\\ \\cline{3-7} \\cline{8-10} \n');
fprintf(FID, '              & Gini  & First  & Second & Third & Fourth & Fifth & 90th-95th & 95th-99th & 99th-100th \\\\ \n \\hline \n');
fprintf(FID, '\\multicolumn{10}{l}{Earnings:} \\\\ \n');
fprintf(FID, 'Model 1       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,1)');
fprintf(FID, 'Model 2       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,2)');
fprintf(FID, 'Model 3       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,3)');
fprintf(FID, 'Model 4       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,4)');
fprintf(FID, 'Model 5       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,5)');
fprintf(FID, 'Model 6       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,6)');
fprintf(FID, '\\multicolumn{10}{l}{Income:} \\\\ \n');
fprintf(FID, 'Model 1       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,1)');
fprintf(FID, 'Model 2       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,2)');
fprintf(FID, 'Model 3       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,3)');
fprintf(FID, 'Model 4       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,4)');
fprintf(FID, 'Model 5       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,5)');
fprintf(FID, 'Model 6       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,6)');
fprintf(FID, '\\multicolumn{10}{l}{Wealth:} \\\\ \n');
fprintf(FID, 'Model 1       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,1)');
fprintf(FID, 'Model 2       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,2)');
fprintf(FID, 'Model 3       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,3)');
fprintf(FID, 'Model 4       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,4)');
fprintf(FID, 'Model 5       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,5)');
fprintf(FID, 'Model 6       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,6)');
fprintf(FID, '\\multicolumn{10}{l}{Consumption:} \\\\ \n');
fprintf(FID, 'Model 1       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,1)');
fprintf(FID, 'Model 2       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,2)');
fprintf(FID, 'Model 3       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,3)');
fprintf(FID, 'Model 4       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,4)');
fprintf(FID, 'Model 5       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,5)');
fprintf(FID, 'Model 6       & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,6)');
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: Income equals earnings plus interest income on assets. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

%Table: Working Age Inequality Version 2: By Model
FID = fopen('./SavedOutput/LatexInputs/EarningsDynamics_WorkingAgeInequalityV2.tex', 'w');
fprintf(FID, 'Working Age Inequality: Share of X held by quintile/top percent \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lccccccccc} \n \\hline \\hline \n');
fprintf(FID, '& & & & & & & \\multicolumn{3}{c}{Top Groups} \\\\ \\cline{8-10} \n');
fprintf(FID, '& & \\multicolumn{5}{c}{Quintile} & \\multicolumn{3}{c}{(Percentile)} \\\\ \\cline{3-7} \\cline{8-10} \n');
fprintf(FID, '              & Gini  & First  & Second & Third & Fourth & Fifth & 90th-95th & 95th-99th & 99th-100th \\\\ \n \\hline \n');
fprintf(FID, '\\multicolumn{10}{l}{Model 1:} \\\\ \n');
fprintf(FID, 'Earnings      & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,1)');
fprintf(FID, 'Income        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,1)');
fprintf(FID, 'Wealth        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,1)');
fprintf(FID, 'Consumption   & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,1)');
fprintf(FID, '\\multicolumn{10}{l}{Model 2:} \\\\ \n');
fprintf(FID, 'Earnings      & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,2)');
fprintf(FID, 'Income        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,2)');
fprintf(FID, 'Wealth        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,2)');
fprintf(FID, 'Consumption   & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,2)');
fprintf(FID, '\\multicolumn{10}{l}{Model 3:} \\\\ \n');
fprintf(FID, 'Earnings      & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,3)');
fprintf(FID, 'Income        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,3)');
fprintf(FID, 'Wealth        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,3)');
fprintf(FID, 'Consumption   & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,3)');
fprintf(FID, '\\multicolumn{10}{l}{Model 4:} \\\\ \n');
fprintf(FID, 'Earnings      & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,4)');
fprintf(FID, 'Income        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,4)');
fprintf(FID, 'Wealth        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,4)');
fprintf(FID, 'Consumption   & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,4)');
fprintf(FID, '\\multicolumn{10}{l}{Model 5:} \\\\ \n');
fprintf(FID, 'Earnings      & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,5)');
fprintf(FID, 'Income        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,5)');
fprintf(FID, 'Wealth        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,5)');
fprintf(FID, 'Consumption   & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,5)');
fprintf(FID, '\\multicolumn{10}{l}{Model 6:} \\\\ \n');
fprintf(FID, 'Earnings      & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', EarningsInequality(:,6)');
fprintf(FID, 'Income        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', IncomeInequality(:,6)');
fprintf(FID, 'Wealth        & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', AssetsInequality(:,6)');
fprintf(FID, 'Consumption   & %8.2f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', ConsumptionInequality(:,6)');
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: Income equals earnings plus interest income on assets. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);


%% Look at lifetime earnings inequality
% To be able to compare to Guvenen, Kaplan, Song & Weidner (2022)

% GKSW2022 do not report lorenz curves, but seem like they would be nice to draw
figure_c=figure_c+1;
figure(figure_c)
legendstr={}; cc=0;
for useModel=useModel_vec
    load(['./SavedOutput/Main/LifeTimeEarnings',num2str(useModel),'.mat'], 'LorenzCurve_LifetimeEarnings')
    if useModel==1
        plot(LorenzCurve_LifetimeEarnings)
        hold on
    else
        plot(LorenzCurve_LifetimeEarnings)
    end
    cc=cc+1; % Just a counter in case not doing all the models
    legendstr{cc}=['Model ',num2str(useModel)];
end
hold off
legend(legendstr{:},'Location','northwest')
saveas(figure_c,'./SavedOutput/Graphs/EarningsDynamics_LifetimeEarningsInequality.png')

% GKSW2022 do report interquartile ratio and the std dev of log, so create table of these
Table_GKSW2022=zeros(6,4); % six models by four measures
for useModel=useModel_vec
    load(['./SavedOutput/Main/LifeTimeEarnings',num2str(useModel),'.mat'],'stddev_logLifetimeEarnings','LifetimeEarnings_P75P25ratio','LifetimeEarnings_P90P50ratio','LifetimeEarnings_P50P10ratio')

    Table_GKSW2022(useModel,1)=stddev_logLifetimeEarnings;
    Table_GKSW2022(useModel,2)=LifetimeEarnings_P75P25ratio;
    Table_GKSW2022(useModel,3)=LifetimeEarnings_P90P50ratio;
    Table_GKSW2022(useModel,4)=LifetimeEarnings_P50P10ratio;
end
%Table: Lifetime earnings inequality (compare to GKSW2022, Figure 8 and Figure 9)
FID = fopen('./SavedOutput/LatexInputs/EarningsDynamics_LifetimeEarningsInequality.tex', 'w');
fprintf(FID, 'Lifetime Earnings Inequality \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcccc} \n \\hline \\hline \n');
fprintf(FID, '  & Std dev of log & interquartile ratio (P75/P25)  & P90/P50 ratio & P50/P10 ratio \\\\ \\hline \n');
fprintf(FID, 'US Data & 0.80  & 2.75  & 2.50  & 2.90  \\\\ \n');
fprintf(FID, 'Model 1 & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table_GKSW2022(1,1), Table_GKSW2022(1,2), Table_GKSW2022(1,4), Table_GKSW2022(1,4));
fprintf(FID, 'Model 2 & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table_GKSW2022(2,1), Table_GKSW2022(2,2), Table_GKSW2022(2,4), Table_GKSW2022(2,4));
fprintf(FID, 'Model 3 & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table_GKSW2022(3,1), Table_GKSW2022(3,2), Table_GKSW2022(3,4), Table_GKSW2022(3,4));
fprintf(FID, 'Model 4 & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table_GKSW2022(4,1), Table_GKSW2022(4,2), Table_GKSW2022(4,4), Table_GKSW2022(4,4));
fprintf(FID, 'Model 5 & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table_GKSW2022(5,1), Table_GKSW2022(5,2), Table_GKSW2022(5,4), Table_GKSW2022(5,4));
fprintf(FID, 'Model 6 & %8.2f & %8.2f & %8.2f & %8.2f \\\\ \n', Table_GKSW2022(6,1), Table_GKSW2022(6,2), Table_GKSW2022(6,4), Table_GKSW2022(6,4));
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: Reports measures of lifetime earnings inequality following Guvenen, Kaplan, Song \\& Weidner (2022). The interquartile ratio is the ratio of the 75th percentile to the 25th percentile. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);


%% Create table about persisitence of poverty
load ./SavedOutput/Main/RankTransitionProbablities.mat RankTransitionProbablities_TableRaw
load ./SavedOutput/Main/ParamsForModel.mat ParamsForModel
Params.Jr=ParamsForModel.model1.Jr; % Jr is the same for all models

Table_PersistenceOfPovery=zeros(14,6);

% By definition the probability of bottom quintile is 0.2
Table_PersistenceOfPovery(1,:)=0.2*ones(1,6);
Table_PersistenceOfPovery(2,:)=0.8*ones(1,6);
% Now we do the two period version
for useModel=useModel_vec
    TransProbs_Quintiles=RankTransitionProbablities_TableRaw.(['model',num2str(useModel)]).earnings;
    % Use just the working ages
    TransProbs_Quintiles=TransProbs_Quintiles(:,:,1:Params.Jr-2); % The transitions essentially end at this point
    % Average over ages
    TransProbs_Quintiles=mean(TransProbs_Quintiles,3);
    TransProbs_Quintiles=TransProbs_Quintiles./sum(TransProbs_Quintiles,2);
    % We are just interested in bottom quintile vs all others
    temp=mean(TransProbs_Quintiles(2:5,:),1);
    temp=[temp(1,1),sum(temp(1,2:5))];
    TransProbs=[TransProbs_Quintiles(1,1),sum(TransProbs_Quintiles(1,2:5)); temp];
    % TransProbs is now 2x2, representing bottom quintile and all other quintiles

    Table_PersistenceOfPovery(3,useModel)=Table_PersistenceOfPovery(1,useModel)*TransProbs(1,1); % --
    Table_PersistenceOfPovery(4,useModel)=Table_PersistenceOfPovery(1,useModel)*TransProbs(1,2); % -+
    Table_PersistenceOfPovery(5,useModel)=Table_PersistenceOfPovery(2,useModel)*TransProbs(2,1); % +-
    Table_PersistenceOfPovery(6,useModel)=Table_PersistenceOfPovery(2,useModel)*TransProbs(2,2); % ++

    Table_PersistenceOfPovery(7,useModel)=Table_PersistenceOfPovery(3,useModel)*TransProbs(1,1); % ---
    Table_PersistenceOfPovery(8,useModel)=Table_PersistenceOfPovery(3,useModel)*TransProbs(1,2); % --+
    Table_PersistenceOfPovery(9,useModel)=Table_PersistenceOfPovery(4,useModel)*TransProbs(2,1); % -+-
    Table_PersistenceOfPovery(10,useModel)=Table_PersistenceOfPovery(4,useModel)*TransProbs(2,2); % -++
    Table_PersistenceOfPovery(11,useModel)=Table_PersistenceOfPovery(5,useModel)*TransProbs(1,1); % +--
    Table_PersistenceOfPovery(12,useModel)=Table_PersistenceOfPovery(5,useModel)*TransProbs(1,2); % +-+
    Table_PersistenceOfPovery(13,useModel)=Table_PersistenceOfPovery(6,useModel)*TransProbs(2,1); % ++-
    Table_PersistenceOfPovery(14,useModel)=Table_PersistenceOfPovery(6,useModel)*TransProbs(2,2); % +++
end

%Table: Persistence of Poverty
FID = fopen('./SavedOutput/LatexInputs/EarningsDynamics_PersistenceOfPoverty.tex', 'w');
fprintf(FID, 'Persistence of Poverty: Probabilities of Earnings Quintile Sequences \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcccccc} \n \\hline \\hline \n');
fprintf(FID, ' & \\multicolumn{6}{l}{Model} \\\\ \n');
fprintf(FID, ' Sequence & 1 & 2 & 3 & 4 & 5 & 6 \\\\ \\hline \n');
fprintf(FID, '$-$   & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table_PersistenceOfPovery(1,:) );
fprintf(FID, '$+$   & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table_PersistenceOfPovery(2,:) );
fprintf(FID, '$--$   & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table_PersistenceOfPovery(3,:) );
fprintf(FID, '$-+$   & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table_PersistenceOfPovery(4,:) );
fprintf(FID, '$+-$   & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table_PersistenceOfPovery(5,:) );
fprintf(FID, '$++$   & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table_PersistenceOfPovery(6,:) );
fprintf(FID, '$---$   & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table_PersistenceOfPovery(7,:) );
fprintf(FID, '$--+$   & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table_PersistenceOfPovery(8,:) );
fprintf(FID, '$-+-$   & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table_PersistenceOfPovery(9,:) );
fprintf(FID, '$-++$   & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table_PersistenceOfPovery(10,:) );
fprintf(FID, '$+--$   & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table_PersistenceOfPovery(11,:) );
fprintf(FID, '$+-+$   & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table_PersistenceOfPovery(12,:) );
fprintf(FID, '$++-$   & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table_PersistenceOfPovery(13,:) );
fprintf(FID, '$+++$   & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f & %8.3f \\\\ \n', Table_PersistenceOfPovery(14,:) );
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: In the column headed Sequence, - indicates the person is in the bottom quintile of the earnings distribution, and + indicates they are not. \n');
fprintf(FID,'A --- indicates the person is in the bottom quintile of the earnings distribution for three consectutive years. \n');
fprintf(FID, 'By definition the first two rows are 0.2 and 0.8. In Geweke \\& Keane (2000) the first two rows relate to subgroups of the population and so need not be 0.2 and 0.8; their numbers are also based on a sample with a different distribution of ages and so are not directly comparable. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);

% Maybe the issue is about the transition probability of upsilon


%% Create table about welfare impacts
% To compare to DeNardi, Fella & Paz-Pardo (2018) results on welfare

Table_WelfareCEV=zeros(2,6);

for useModel=useModel_vec
    load(['./SavedOutput/Main/BasicOutput',num2str(useModel),'.mat'],'n_upsilon')
    if n_upsilon==1
        load(['./SavedOutput/Main/WelfareCEV',num2str(useModel),'.mat'],'CEV1')
        CEV2=CEV1;
    else
        load(['./SavedOutput/Main/WelfareCEV',num2str(useModel),'.mat'],'CEV1','CEV2')
    end
    Table_WelfareCEV(1,useModel)=CEV1; % Turn of upsilon, z, and epsilon
    Table_WelfareCEV(2,useModel)=CEV2; % Additionally turn off alpha (can kappa_beta)
end

%Table: Welfare cost of idiosyncratic shocks
FID = fopen('./SavedOutput/LatexInputs/EarningsDynamics_Welfare.tex', 'w');
fprintf(FID, 'Consumption Equivalance measure of welfare cost of shocks \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcccccc} \n \\hline \\hline \n');
fprintf(FID, ' \\multicolumn{7}{c}{Eliminating idiosyncratic shocks} \\\\ \n');
fprintf(FID, 'Model: & 1 & 2 & 3 & 4 & 5 & 6 \\\\ \\hline \n');
fprintf(FID, 'CEV (\\%) :   & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\  \\hline \n', 100*Table_WelfareCEV(1,:) );
fprintf(FID, ' \\multicolumn{7}{c}{Further eliminating permanent types} \\\\ \n');
fprintf(FID, 'Model: & 1 & 2 & 3 & 4 & 5 & 6 \\\\ \\hline \n');
fprintf(FID, 'CEV (\\%) :   & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f & %8.1f \\\\ \n', 100*Table_WelfareCEV(2,:) );
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: The deterministic model is different for each model. For all models: Set z=0 and epsilon=0. Set upsilon=0 and replace $\\kappa_j$ with the age-conditional mean earnings (so turning off all shocks leaves age-conditional mean earnings unchanged, conditional on permanent type). When further eliminating permanent type all agents get the age-conditional mean earnings (unconditional).  \n');
fprintf(FID, 'A double-check showed that age-conditional mean earnings are unchanged (change was less than 10^(-9)). For Model 4 we only tried as high as 300\\%, actual value is higher. \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);


%% Look at consumption insurance against income shocks

% Life-cycle profiles for the variance of log consumption and log disposible income
load ./SavedOutput/Main/ParamsForModel.mat ParamsForModel
figure_c=figure_c+1;
figure(figure_c)
for useModel=useModel_vec
    load(['./SavedOutput/Main/ConsumptionInsurance',num2str(useModel),'.mat'],'VarLogDispIncome_LCP','VarLogCons_LCP')
    Params=ParamsForModel.(['model',num2str(useModel)]);
    
    subplot(3,2,useModel); plot(Params.ageyears,VarLogCons_LCP,'.')
    hold on
    subplot(3,2,useModel); plot(Params.ageyears,VarLogDispIncome_LCP,'.')
    hold off
    title(['Model', num2str(useModel)])
end
legend('Var[ln(c)]','Var[ln(disp. income)]')
saveas(figure_c,'./SavedOutput/Graphs/EarningsDynamics_Var_lnc_lndispincome.png')

% Redo the previous graph, but this time just showing the working age.
% (The full graph looks a bit silly because of some silly variance in
% the earnings process at the end of working life in some of the earnings processes)
figure_c=figure_c+1;
figure(figure_c)
for useModel=useModel_vec
    load(['./SavedOutput/Main/ConsumptionInsurance',num2str(useModel),'.mat'],'VarLogDispIncome_LCP','VarLogCons_LCP')
    Params=ParamsForModel.(['model',num2str(useModel)]);
    
    subplot(3,2,useModel); plot(Params.ageyears(1:Params.Jr-1),VarLogCons_LCP(1:Params.Jr-1),'.')
    hold on
    subplot(3,2,useModel); plot(Params.ageyears(1:Params.Jr-1),VarLogDispIncome_LCP(1:Params.Jr-1),'.')
    hold off
    title(['Model', num2str(useModel)])
end
legend('Var[ln(c)]','Var[ln(disp. income)]')
saveas(figure_c,'./SavedOutput/Graphs/EarningsDynamics_Var_lnc_lndispincome_WorkingAge.png')


% Life-cycle profiles for the variance of consumption and disposible income (to show
% that decreasing variance of consumption in model 2 is an artifact of the log)
% Note: When just comparing graphs it could also be an artifact of the
% sample restriction (that earnings>0, needed to be able to take logs), but
% more detailed investigation (commented out code below in this section)
% shows it is the log.
figure_c=figure_c+1;
figure(figure_c)
for useModel=useModel_vec
    load(['./SavedOutput/Main/AgeConditionalStats',num2str(useModel),'.mat'])
    subplot(3,2,useModel); plot(Params.ageyears,AgeConditionalStats.consumption.Variance,'.')
    hold on
    subplot(3,2,useModel); plot(Params.ageyears,AgeConditionalStats.disposableincome.Variance,'.')
    hold off
    title(['Model', num2str(useModel)])
end
legend('Var[c]','Var[disp. income]')
saveas(figure_c,'./SavedOutput/Graphs/EarningsDynamics_Var_c_dispincome.png')


% Blundell, Pistaferri & Preston (2008) coefficients
Table_BPPcoeffs=zeros(6,2); % six models by two BPP coeffs
for useModel=useModel_vec
    load(['./SavedOutput/Main/ConsumptionInsurance',num2str(useModel),'.mat'], 'ConsumptionInsurance_BPP2008coeffs')

    Table_BPPcoeffs(useModel,1)=ConsumptionInsurance_BPP2008coeffs.persistent;
    Table_BPPcoeffs(useModel,2)=ConsumptionInsurance_BPP2008coeffs.transitory;
end
%Table: Consumption insurance against income shocks
FID = fopen('./SavedOutput/LatexInputs/EarningsDynamics_BPPcoeffs.tex', 'w');
fprintf(FID, 'Consumption Insurance against Income Shocks \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lcc} \n \\hline \\hline \n');
fprintf(FID, ' & \\multicolumn{2}{c}{BPP2008 Consumption Insurance Coeff.s} \\\\ \n');
fprintf(FID, ' & Persistent ($\\phi^p$) & Transitory ($\\phi^{tr}$) \\\\ \\hline \n');
fprintf(FID, 'US Data & 0.36  & 0.95 \\\\ \n');
fprintf(FID, 'Model 1 & %8.2f & %8.2f \\\\ \n', Table_BPPcoeffs(1,1), Table_BPPcoeffs(1,2));
fprintf(FID, 'Model 2 & %8.2f & %8.2f \\\\ \n', Table_BPPcoeffs(2,1), Table_BPPcoeffs(2,2));
fprintf(FID, 'Model 3 & %8.2f & %8.2f \\\\ \n', Table_BPPcoeffs(3,1), Table_BPPcoeffs(3,2));
fprintf(FID, 'Model 4 & %8.2f & %8.2f \\\\ \n', Table_BPPcoeffs(4,1), Table_BPPcoeffs(4,2));
fprintf(FID, 'Model 5 & %8.2f & %8.2f \\\\ \n', Table_BPPcoeffs(5,1), Table_BPPcoeffs(5,2));
fprintf(FID, 'Model 6 & %8.2f & %8.2f \\\\ \n', Table_BPPcoeffs(6,1), Table_BPPcoeffs(6,2));
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, 'Note: Reports consumption insurance coefficients following Blundell, Pistaferri \\& Preston (2008) measuring insurance against persisent and transitory shocks, respectively. Range from 0, representing no insurance, to 1 representing perfect insurance (these are actually 1- the BPP verions, to ease interpretation). \n');
fprintf(FID, '}} \\end{minipage}');
fclose(FID);


% % To better understand what is going on with consumption in a specific model
% useModel=2;
% load ./SavedOutput/Main/ParamsForModel.mat ParamsForModel
% Params=ParamsForModel.(['model',num2str(useModel)]);
% load(['./SavedOutput/Main/SimPanelValues',num2str(useModel),'.mat'])
% % Look at the age-condtional cdf for consumption
% figure_c=figure_c+1;
% figure(figure_c)
% c_jj=SimPanelValues.consumption(1,:);
% subplot(2,1,1); cdfplot(c_jj)
% xlim([min(min(SimPanelValues.consumption)),max(max(SimPanelValues.consumption))])
% hold on
% for jj=[11,21,31,41,51]
%     c_jj=SimPanelValues.consumption(jj,:);
%     subplot(2,1,1); cdfplot(c_jj)
% end
% hold off
% legend('agej=1','agej=11','agej=21','agej=31','agej=41','agej=51')
% 
% e_jj=SimPanelValues.consumption(1,:);
% subplot(2,1,2); cdfplot(e_jj)
% xlim([min(min(SimPanelValues.earnings)),max(max(SimPanelValues.earnings))])
% hold on
% for jj=[11,21,31,41,51]
%     e_jj=SimPanelValues.earnings(jj,:);
%     subplot(2,1,2); cdfplot(e_jj)
% end
% hold off
% legend('agej=1','agej=11','agej=21','agej=31','agej=41','agej=51')
% 
% % What if I look at variance of consumption rather than variance of log consumption?
% figure_c=figure_c+1;
% figure(figure_c)
% load(['./SavedOutput/Main/AgeConditionalStats',num2str(useModel),'.mat'])
% plot(Params.ageyears,AgeConditionalStats.consumption.Variance,'.')
% hold on
% plot(Params.ageyears,AgeConditionalStats.earnings.Variance,'.')
% hold off
% legend('var(c)', 'var(earnings)')
% % Without restricting the sample
% Csample=SimPanelValues.consumption;
% Ysample=SimPanelValues.earnings;
% age=SimPanelValues.agej;
% for jj=1:Params.J
%     Csample_j=Csample(logical(age==jj));
%     Ysample_j=Ysample(logical(age==jj));
%     VarEarnings_LCP(jj)=var(Ysample_j);
%     VarCons_LCP(jj)=var(Csample_j);
%     VarLogCons_LCP3(jj)=var(log(Csample_j));
% end
% % With the sample restriction (earnings>0, which is required to make it possible to log())
% keep=logical(SimPanelValues.earnings>0); % for model 2, keeps roughly 46% of observations
% Csample=SimPanelValues.consumption(keep);
% Ysample=SimPanelValues.earnings(keep);
% age=SimPanelValues.agej(keep);
% for jj=1:Params.J
%     Csample_j=Csample(logical(age==jj));
%     Ysample_j=Ysample(logical(age==jj));
%     VarEarnings_LCP2(jj)=var(Ysample_j);
%     VarCons_LCP2(jj)=var(Csample_j);
%     VarLogEarnings_LCP2(jj)=var(log(Ysample_j));
%     VarLogCons_LCP2(jj)=var(log(Csample_j));
% end
% % Looking at: VarEarnings_LCP, VarEarnings_LCP2, VarCons_LCP, and VarCons_LCP2
% % Seems likely that the odd graph is a product of log() somehow.
% % Note that: VarCons_LCP2 increase with age (restricted sample, without
% % log), but that VarLogCons_LCP2 decreases with age (restricted sample,
% % with log). VarLogCons_LCP3 (unrestricted sample, with log) displays the
% % decreasing in age, but VarCons_LCP (unrestricted sample, without log) is
% % increasing in age (until retirement).
% 
% % Note, to calculate the age-conditional variance of consumption and
% % earnings I have to use a restricted sample, try replotting the ecdf using
% % the restricted sample.
% keep=logical(SimPanelValues.earnings>0); % for model 2, keeps roughly 46% of observations
% logC=SimPanelValues.logconsumption(keep);
% logY=SimPanelValues.logearnings(keep);
% age=SimPanelValues.agej(keep);
% VarLogEarnings_LCP=zeros(Params.J,1);
% VarLogCons_LCP=zeros(Params.J,1);
% for jj=1:Params.J
%     logC_j=logC(logical(age==jj));
%     logY_j=logY(logical(age==jj));
%     VarLogEarnings_LCP(jj)=var(logY_j);
%     VarLogCons_LCP(jj)=var(logC_j);
% end



%% Compare the shock processes on z and epsilon

% Plot some of the life-cycle profiles for the variance of the shocks
load ./SavedOutput/Main/ParamsForModel.mat ParamsForModel
figure_c=figure_c+1;
figure(figure_c)
legendstr={}; cc=0;
for useModel=useModel_vec
    load(['./SavedOutput/Main/AgeConditionalStats',num2str(useModel),'.mat'])
    Params=ParamsForModel.(['model',num2str(useModel)]);
    
    hold on
    subplot(4,1,1); plot(Params.ageyears,AgeConditionalStats.z.StdDev)
    hold off
    hold on
    subplot(4,1,2); plot(Params.ageyears,AgeConditionalStats.epsilon.StdDev)
    hold off
    hold on
    subplot(4,1,3); plot(Params.ageyears,AgeConditionalStats.upsilon.Mean)
    hold off
    hold on
    subplot(4,1,4); plot(Params.ageyears,Params.kappa_j)
    hold off
    
    cc=cc+1; % Just a counter in case not doing all the models
    legendstr{cc}=['Model ',num2str(useModel)];
end
subplot(4,1,1); title('Age-Conditional Standard Deviation of z')
subplot(4,1,2); title('Age-Conditional Standard Deviation of epsilon')
subplot(4,1,3); title('Age-Conditional Mean of upsilon')
subplot(4,1,4); title('Age-Conditional deterministic earnings, kappa')
legend(legendstr)
saveas(figure_c,'./SavedOutput/Graphs/EarningsDynamics_CompareShocks.png')


%% Look at how prob(upsilon=1) varies with z (and age)
% GKOS2021 use that: prob(upsilon=1)=c0+c1*z+c2*t+c3+z+t (where t=(age-24)/10)
% Model 2: -3.036,  -0.917,   -5.397,   -4.442
% Model 5: -2.495,  -1.037,   -5.051,   -1.087
% Model 6: -3.353,  -0.859,   -5.034,   -2.895
% So in all cases the probability decreases as z and t increase

% Note that the probability is independent of lag of z, and of lag of upsilon
% So will just plot the probability that upsilon is 1, conditional on z,
% for various ages

agevec=[1,11,21,31,41];
legendstr={};
for useModel=[2,5,6]
    load(['./SavedOutput/BasicDiscretization/DiscretizedEarnings',num2str(useModel),'nSigmaz',num2str(nSigmaz),'nSigma_alpha',num2str(nSigma_alpha),'.mat'])
    figure_c=figure_c+1;
    figure(figure_c)
    legendstr{1}=num2str(agevec(1)+Params.agejshifter);
    plot(z_grid_J(:,1),pi_zupsilon_J(1,n_z+1:end,agevec(1))); 
    % Arbitrarily use first row, because it is anyway independent of z and upsilon
    % n_z+1:end corresponds to upsilon=1
    hold on
    for jj=2:length(agevec)
        legendstr{jj}=num2str(agevec(jj)+Params.agejshifter);
        plot(z_grid_J(:,jj),pi_zupsilon_J(1,n_z+1:end,agevec(jj)));
    end
    hold off
    legend(legendstr)
    xlabel('z')
    ylim([0,1])
    title('Prob(upsilon=1) in terms of z (for various ages)')
    saveas(figure_c,['./SavedOutput/Graphs/EarningsDynamics_ProbUnemp_useModel',num2str(useModel),'.png'])
end

% Following codes print out a bunch of zeros, confirming that all rows are identical as expected
% for jj=1:size(pi_zupsilon_J,3)
%     for ii=1:size(pi_zupsilon_J,1)
%         max(max(pi_zupsilon_J(1,:,jj)-pi_zupsilon_J(ii,:,jj)))
%     end
% end



