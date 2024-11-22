% Compare discretiziation of the AR1wGM from model 3 using different grid sizes
% (Model three has the smallest/simplest AR1wGM)

% This script simply loops over the DiscretizeEarningsDynamicsModel.m and
% EvaluateDiscretization.m files to create and evaluate all the discretizations.

% Lets model agents from age 25 to age 100, so 76 periods

Params.agejshifter=24; % Age 25 minus one. Makes keeping track of actual age easy in terms of model age
Params.J=100-Params.agejshifter; % Number of period in life-cycle

% Note: grid on assets is not relavant to the discretization of the earnings process

n_epsilon=15; % Transitory earnings shock

% Permanent types
n_alpha=5; % Fixed effect

N_j=Params.J; % Number of periods in finite horizon

% Exogenous states

%% 
nzvec=[101,75,51,41,31,21,11];
nSigmazvec=[2,3];
useModelvec=[3,5,6];

% Table about how much of z grid is being used in various specifications
Table_truerho=zeros(41,length(nzvec),length(nSigmazvec),3); % last dimension is the model
Table_discrho=zeros(41,length(nzvec),length(nSigmazvec),3); % last dimension is the model
Table_true_uncondlmean=zeros(41,length(nzvec),length(nSigmazvec),3); % last dimension is the model
Table_true_uncondlvar=zeros(41,length(nzvec),length(nSigmazvec),3); % last dimension is the model
Table_disc_uncondlmean=zeros(41,length(nzvec),length(nSigmazvec),3); % last dimension is the model
Table_disc_uncondlvar=zeros(41,length(nzvec),length(nSigmazvec),3); % last dimension is the model
Table_disc_fraction4=zeros(length(nzvec),length(nSigmazvec),3); % last dimension is the model
% note: 41 is number of periods for working age

%%
nSigma_alpha=1; % Not relevant to what currently doing
for i1=1:length(nzvec)
    n_z=nzvec(i1);
    for i2=1:length(nSigmazvec)
        nSigmaz=nSigmazvec(i2);
        % Use: 2,3 (number of standard deviations for the max and min points of grid used for the discretized shocks: z, epsilon
        % Note: nSigmaz is really just intended for use with EvaluateDiscretization.
        % Use a different number of standard deviations for the max and min points of grids used for the discretization of the fixed-effect (and HIP): alpha (and kappa_beta)
        for i3=1:length(useModelvec)
            useModel=useModelvec(i3);
            % Models 3,5,6 are the ones that use AR1wGM, which is what we want to evaluate the accuracy of the discretization for.

            fprintf('Currently discretizing and evaluating useModel=%i and nSigmaz=%i \n',useModel,nSigmaz)

            if useModel==2 || useModel==5 || useModel==6 || useModel==21 || useModel==22
                n_upsilon=2; % non-employment shock
            else
                n_upsilon=1; % unused (no non-employment shock)
            end
            n_zupsilon=[n_z, n_upsilon];

            if useModel==6
                n_kappabeta=n_alpha;  % Heterogenous income profile (GKOS2021 call this beta; set below depending on model)
                % Note: discretization routine requires that n_kappabeta=n_alpha
            else
                n_kappabeta=1; % unused
            end

            % There are n_alpha*n_kappabeta permanent types
            N_i=n_alpha*n_kappabeta;
            % Note that n_kappabeta=1 when there is not permanent types of kappabeta (no Heterogeneous Income Profiles)

            %% Discretized earnings dynamics
            % Creates Params
            % Creates exogenous states
            % Creates permanent types
            DiscretizeEarningsDynamicsModel
            EvaluateDiscretization2
            
            %% I want to create a Table that reports accuracy on the
            % autocorrelation
            Table_true_rho(:,i1,i2,i3)=Params.rho(1:41)';
            Table_disc_rho(:,i1,i2,i3)=corr_disc_z';
            % unconditional moments
            Table_true_uncondlmean(:,i1,i2,i3)=zeros(41,1);
            Table_true_uncondlvar(:,i1,i2,i3)=otheroutputs_z.sigma_z';
            Table_disc_uncondlmean(:,i1,i2,i3)=mean_disc_z(1:41)';
            Table_disc_uncondlvar(:,i1,i2,i3)=variance_disc_z(1:41)';
            % fraction of transitions for which hit 4 moments
            Table_disc_fraction4(i1,i2,i3)=sum(sum(otheroutputs_z.nMoments_grid==4))/numel(otheroutputs_z.nMoments_grid);
        end
    end
end

save ./SavedOutput/CompareGridSize/CompareAccuracyGridSizes.mat

%% Now use these things to compute accuracy
% autocorrelation
autocorrelation_accuracy=(Table_disc_rho(2:41,:,:,:)-Table_true_rho(2:41,:,:,:))./Table_true_rho(2:41,:,:,:); % autocorrelation in 1st period is not observable
autocorrelation_accuracy=squeeze(mean(autocorrelation_accuracy,1));
% (n_z, nSigmaz, useModel) are the three dimensions
autocorrelation_accuracy(:,:,3)

% unconditional mean
unconditionalmean_accuracy=(Table_disc_uncondlmean-Table_true_uncondlmean); % true is often zero, so cannot measure as fraction
unconditionalmean_accuracy=squeeze(mean(unconditionalmean_accuracy,1));
unconditionalmean_accuracy(:,:,3)
% unconditional variance
unconditionalvar_accuracy=(Table_disc_uncondlvar-Table_true_uncondlvar)./Table_true_uncondlvar;
unconditionalvar_accuracy=squeeze(mean(unconditionalvar_accuracy,1));
unconditionalvar_accuracy(:,:,3)

% fraction of transitions for which hit 4 moments
Table_disc_fraction4(:,:,3)

%% Create a table
% Table 6
FID = fopen('./SavedOutput/LatexInputs/CompareAccuracyofGrids.tex', 'w');
fprintf(FID, 'Evaluating Quadrature for Different Size Grids \\\\ \n');
fprintf(FID, '\\begin{tabular*}{1.00\\textwidth}{@{\\extracolsep{\\fill}}lllll} \n \\hline \\hline \n');
fprintf(FID, ' \\# grid points ($n\\_{}z$) & Autocorr. & Uncondl Mean & Uncondl Variance & Fraction hit 4 moments \\\\ \\hline \n');
fprintf(FID, ' \\multicolumn{5}{c}{Model 3 with nSigmaz=2} \\\\ \\hline \n');
for i1=1:length(nzvec)
    fprintf(FID, ' %i & %1.3f & %1.3f & %1.3f & %1.3f \\\\ \n',nzvec(i1), autocorrelation_accuracy(i1,1,3), unconditionalmean_accuracy(i1,1,3), unconditionalvar_accuracy(i1,1,3), Table_disc_fraction4(i1,1,3));
end
% fprintf(FID, '\\multicolumn{5}{c}{Model 3 with nSigmaz=3} \\\\ \\hline \n');
% for i1=1:length(nzvec)
%     fprintf(FID, ' %i & %1.3f & %1.3f & %1.3f & %1.3f \\\\ \n',nzvec(i1), autocorrelation_accuracy(i1,2,3), unconditionalmean_accuracy(i1,2,3), unconditionalvar_accuracy(i1,2,3), Table_disc_fraction4(i1,2,3));
% end
% fprintf(FID, '\\multicolumn{5}{c}{Model 5 with nSigmaz=2} \\\\ \\hline \n');
% for i1=1:length(nzvec)
%     fprintf(FID, ' %i & %1.3f & %1.3f & %1.3f & %1.3f \\\\ \n',nzvec(i1), autocorrelation_accuracy(i1,1,5), unconditionalmean_accuracy(i1,1,5), unconditionalvar_accuracy(i1,1,5), Table_disc_fraction4(i1,1,5));
% end
% fprintf(FID, '\\multicolumn{5}{c}{Model 5 with nSigmaz=3} \\\\ \\hline \n');
% for i1=1:length(nzvec)
%     fprintf(FID, ' %i & %1.3f & %1.3f & %1.3f & %1.3f \\\\ \n',nzvec(i1), autocorrelation_accuracy(i1,2,5), unconditionalmean_accuracy(i1,2,5), unconditionalvar_accuracy(i1,2,5), Table_disc_fraction4(i1,2,5));
% end
% fprintf(FID, '\\multicolumn{5}{c}{Model 6 with nSigmaz=2} \\\\ \\hline \n');
% for i1=1:length(nzvec)
%     fprintf(FID, ' %i & %1.3f & %1.3f & %1.3f & %1.3f \\\\ \n',nzvec(i1), autocorrelation_accuracy(i1,1,6), unconditionalmean_accuracy(i1,1,6), unconditionalvar_accuracy(i1,1,6), Table_disc_fraction4(i1,1,6));
% end
% fprintf(FID, '\\multicolumn{5}{c}{Model 6 with nSigmaz=3} \\\\ \\hline \n');
% for i1=1:length(nzvec)
%     fprintf(FID, ' %i & %1.3f & %1.3f & %1.3f & %1.3f \\\\ \n',nzvec(i1), autocorrelation_accuracy(i1,2,6), unconditionalmean_accuracy(i1,2,6), unconditionalvar_accuracy(i1,2,6), Table_disc_fraction4(i1,2,6));
% end
fprintf(FID, '\\hline \n \\end{tabular*} \n');
fprintf(FID, '\\begin{minipage}[t]{1.00\\textwidth}{\\baselineskip=.5\\baselineskip \\vspace{.3cm} \\footnotesize{ \n');
fprintf(FID, ['Note: Accuracy of the autocorrelation coefficient, the unconditional mean and the unconditional variance are evaluated as the absolute deviation from the true value as a fraction of the true value, averaged over ages. ' ...
    'Fraction hit four moments is the number of grid points for which the transition probabilities were able to be chosen to hit all four target conditional moments.' ...
    'Model 3, 5 and 6 are the three models with gaussian-mixture innovations The results reported relate only to the discretization' ...
    'of the life-cycle AR(1) process. \n']);
fprintf(FID, '}} \\end{minipage}');
fclose(FID);



