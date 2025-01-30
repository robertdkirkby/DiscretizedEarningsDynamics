% This script simply loops over the DiscretizeEarningsDynamicsModel.m and
% EvaluateDiscretization.m files to create and evaluate all the discretizations.

% Lets model agents from age 25 to age 100, so 76 periods

Params.agejshifter=24; % Age 25 minus one. Makes keeping track of actual age easy in terms of model age
Params.J=100-Params.agejshifter; % Number of period in life-cycle

% Note: grid on assets is not relavant to the discretization of the earnings process

% Permanent types
n_alpha=5; % Fixed effect

N_j=Params.J; % Number of periods in finite horizon

% Exogenous states
% n_z=51; % Persistent earnings shock
% n_epsilon=21; % Transitory earnings shock

% Table about how much of z grid is being used in various specifications
TableAboutzgrid=struct();

% The one used in the paper is
% n_z=21; nSigmaz=2; nSigma_alpha=1; n_epsilon=9;
% nzvec=21;
% nSigmazvec=2;
% nSigma_alphavec=1;
n_epsilon=9; % Transitory earnings shock


nzvec=[101,75,51,31,21];
% Note: 21 is done last so that the save files just contain the 21 point discretization (it is the one I use for the Life-Cycle models)
% Because it would be weird to use this same ordering of n_z in the tables
% I also create a version for use in creating tables
nzvectable=sort(nzvec);
nSigmazvec=[2,3,4];
nSigma_alphavec=[1,2];
for n_z=nzvec
    % First, do Models 2 and 6, and consider different nSigmaz
    for useModel=[1,11,12,2,21,22,3,4,5,6] % Can be: 1,2,3,4,5,6
        % There is also useModel 11 and 12, these are just same as 1, except using extended Farmer-Toda to target 2 and 4 moments respectively
        % There is also useModel 21 and 22, these are just same as 2, except using extended Farmer-Toda to target 2 and 4 moments respectively
        for nSigmaz=nSigmazvec
            % Use: 2,3,4 (number of standard deviations for the max and min points of grid used for the discretized shocks: z, epsilon
            % Note: nSigmaz is really just intended for use with EvaluateDiscretization.
            for nSigma_alpha=nSigma_alphavec
                % Use a different number of standard deviations for the max and min points of grids used for the discretization of the fixed-effect (and HIP): alpha (and kappa_beta)
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

                if useModel==4
                    nSigmaz=0.3;
                end

                %% Dicretized earnings dynamics
                % Creates Params
                % Creates exogenous states
                % Creates permanent types
                DiscretizeEarningsDynamicsModel
                EvaluateDiscretization
                save(['./SavedOutput/BasicDiscretization/DiscretizedEarnings',num2str(useModel),'nSigmaz',num2str(nSigmaz),'nSigma_alpha',num2str(nSigma_alpha),'.mat']', 'zupsilon_grid_J','epsilon_grid', 'pi_zupsilon_J', 'pi_epsilon', 'jequaloneDistzupsilonepsilon', 'jequaloneDistzupsilon', 'PTypeDistParamNames', 'Params','z_grid_J','pi_z_J','epsilon_grid_J','pi_epsilon_J','Epi_upsilon_J','-v7.3');

                % TableAboutzgridtemp is created by 'EvaluateDiscretization'
                TableAboutzgrid.(['useModel',num2str(useModel)]).(['nSigmaz',num2str(nSigmaz*10)]).(['nz',num2str(n_z)])=TableAboutzgridtemp; % Note that these are all independent of nSigma_alpha anyway
                TableAboutpiupsilon.(['useModel',num2str(useModel)]).(['nSigmaz',num2str(nSigmaz*10)]).(['nz',num2str(n_z)])=Epi_upsilon_J;
                TableAboutUsingz_z_grid_J.(['useModel',num2str(useModel)]).(['nSigmaz',num2str(nSigmaz*10)]).(['nz',num2str(n_z)])=z_grid_J;
                TableAboutUsingz_statdist_z.(['useModel',num2str(useModel)]).(['nSigmaz',num2str(nSigmaz*10)]).(['nz',num2str(n_z)])=statdist_z;
            end
        end
    end
end



%% Create a 'Table' that summarizes how much of the z_grid ends up actually being used by the discretization
% Note: I don't end up using this. Instead EvaluateDiscretization is
% creating a seperate one for each of the useModel/nSimgaz/n_z combos and these are used.
jjvec=[1,11,21,31,41];
FID = fopen(['./SavedOutput/EvaluateDiscretization/Table_HowMuchzgridUsed.tex'], 'w');
fprintf(FID, 'How much of the z grid is being used? \\\\ \n');
useModel=2;
for n_z=nzvectable
    for nSigmaz=nSigmazvec
        TableAboutzgridtemp=TableAboutzgrid.(['useModel',num2str(useModel)]).(['nSigmaz',num2str(nSigmaz*10)]).(['nz',num2str(n_z)]);
        fprintf(FID, 'Model 2: Rouwenhorst, nSigma=%i, $n_z$=%i \\\\ \n',nSigmaz,n_z);
        for jjc=1:5
            jj=jjvec(jjc);
            fprintf(FID, '  Number of points with mass greater than 10 to -9: %i of %i at age j=%i \\\\ \n',TableAboutzgridtemp(jjc,1) ,n_z,jj);
        end
        for jjc=1:5
            jj=jjvec(jjc);
        fprintf(FID,'  Bottom 10 points sum to mass %8.6f, top 10 points sum to mass %8.6f, at age j=%i \\\\ \n',TableAboutzgridtemp(jjc,2),TableAboutzgridtemp(jjc,3),jj);
        end
        for jjc=1:5
            jj=jjvec(jjc);
        fprintf(FID,'  Bottom 5 points sum to mass %8.6f, top 5 points sum to mass %8.6f, at age j=%i \\\\ \n',TableAboutzgridtemp(jjc,4),TableAboutzgridtemp(jjc,5),jj);
        end
        for jjc=1:5
            jj=jjvec(jjc);
            fprintf(FID,'  Min grid point %8.6f, max grid point %8.6f, at age j=%i  \\\\ \n',TableAboutzgridtemp(jjc,6),TableAboutzgridtemp(jjc,7),jj);
        end
    end
end
useModel=21;
for n_z=nzvectable
    for nSigmaz=nSigmazvec
        TableAboutzgridtemp=TableAboutzgrid.(['useModel',num2str(useModel)]).(['nSigmaz',num2str(nSigmaz*10)]).(['nz',num2str(n_z)]);
        fprintf(FID, 'Model 2: Farmer-Toda targetting 2 moments, nSigma=%i, $n_z$=%i \\\\ \n',nSigmaz,n_z);
        for jjc=1:5
            jj=jjvec(jjc);
            fprintf(FID, '  Number of points with mass greater than 10 to -9: %i of %i at age j=%i \\\\ \n',TableAboutzgridtemp(jjc,1) ,n_z,jj);
        end
        for jjc=1:5
            jj=jjvec(jjc);
        fprintf(FID,'  Bottom 10 points sum to mass %8.6f, top 10 points sum to mass %8.6f, at age j=%i \\\\ \n',TableAboutzgridtemp(jjc,2),TableAboutzgridtemp(jjc,3),jj);
        end
        for jjc=1:5
            jj=jjvec(jjc);
        fprintf(FID,'  Bottom 5 points sum to mass %8.6f, top 5 points sum to mass %8.6f, at age j=%i \\\\ \n',TableAboutzgridtemp(jjc,4),TableAboutzgridtemp(jjc,5),jj);
        end
        for jjc=1:5
            jj=jjvec(jjc);
            fprintf(FID,'  Min grid point %8.6f, max grid point %8.6f, at age j=%i  \\\\ \n',TableAboutzgridtemp(jjc,6),TableAboutzgridtemp(jjc,7),jj);
        end
    end
end
useModel=22;
for n_z=nzvectable
    for nSigmaz=nSigmazvec
        TableAboutzgridtemp=TableAboutzgrid.(['useModel',num2str(useModel)]).(['nSigmaz',num2str(nSigmaz*10)]).(['nz',num2str(n_z)]);
        fprintf(FID, 'Model 2: Farmer-Toda targetting 4 moments, nSigma=%i, $n_z$=%i \\\\ \n',nSigmaz,n_z);
        for jjc=1:5
            jj=jjvec(jjc);
            fprintf(FID, '  Number of points with mass greater than 10 to -9: %i of %i at age j=%i \\\\ \n',TableAboutzgridtemp(jjc,1) ,n_z,jj);
        end
        for jjc=1:5
            jj=jjvec(jjc);
        fprintf(FID,'  Bottom 10 points sum to mass %8.6f, top 10 points sum to mass %8.6f, at age j=%i \\\\ \n',TableAboutzgridtemp(jjc,2),TableAboutzgridtemp(jjc,3),jj);
        end
        for jjc=1:5
            jj=jjvec(jjc);
        fprintf(FID,'  Bottom 5 points sum to mass %8.6f, top 5 points sum to mass %8.6f, at age j=%i \\\\ \n',TableAboutzgridtemp(jjc,4),TableAboutzgridtemp(jjc,5),jj);
        end
        for jjc=1:5
            jj=jjvec(jjc);
            fprintf(FID,'  Min grid point %8.6f, max grid point %8.6f, at age j=%i  \\\\ \n',TableAboutzgridtemp(jjc,6),TableAboutzgridtemp(jjc,7),jj);
        end
    end
end
useModel=6;
for n_z=nzvectable
    for nSigmaz=nSigmazvec
        TableAboutzgridtemp=TableAboutzgrid.(['useModel',num2str(useModel)]).(['nSigmaz',num2str(nSigmaz*10)]).(['nz',num2str(n_z)]);
        fprintf(FID, 'Model 6: Farmer-Toda targetting 4 moments, nSigma=%i, $n_z$=%i \\\\ \n',nSigmaz,n_z);
        for jjc=1:5
            jj=jjvec(jjc);
            fprintf(FID, '  Number of points with mass greater than 10 to -9: %i of %i at age j=%i \\\\ \n',TableAboutzgridtemp(jjc,1) ,n_z,jj);
        end
        for jjc=1:5
            jj=jjvec(jjc);
        fprintf(FID,'  Bottom 10 points sum to mass %8.6f, top 10 points sum to mass %8.6f, at age j=%i \\\\ \n',TableAboutzgridtemp(jjc,2),TableAboutzgridtemp(jjc,3),jj);
        end
        for jjc=1:5
            jj=jjvec(jjc);
        fprintf(FID,'  Bottom 5 points sum to mass %8.6f, top 5 points sum to mass %8.6f, at age j=%i \\\\ \n',TableAboutzgridtemp(jjc,4),TableAboutzgridtemp(jjc,5),jj);
        end
        for jjc=1:5
            jj=jjvec(jjc);
            fprintf(FID,'  Min grid point %8.6f, max grid point %8.6f, at age j=%i  \\\\ \n',TableAboutzgridtemp(jjc,6),TableAboutzgridtemp(jjc,7),jj);
        end
    end
end
fclose(FID);

%% Since Model 2 is normally distributed we know exactly what the mass should be for each additional +-sigma.
% Create a table that looks at this.
FID = fopen(['./SavedOutput/EvaluateDiscretization/Table_Evaluate_Normalz_stddev.tex'], 'w');
fprintf(FID, 'For a normal distribution +-1$\\sigma$ covers 0.68 of the probability distribution \\\\ \n');
fprintf(FID, 'For a normal distribution +-2$\\sigma$ covers 0.954 of the probability distribution \\\\ \n');
fprintf(FID, 'For a normal distribution +-3$\\sigma$ covers 0.997 of the probability distribution \\\\ \n');
fprintf(FID, 'For a normal distribution +-4$\\sigma$ covers over 0.9999 of the probability distribution \\\\ \n');
fprintf(FID, 'For Model 2, at age j=21, we get the following \n');
for useModel=[2,21,22]
    for n_z=nzvectable
        for nSigmaz=nSigmazvec
            z_grid_J=TableAboutUsingz_z_grid_J.(['useModel',num2str(useModel)]).(['nSigmaz',num2str(nSigmaz*10)]).(['nz',num2str(n_z)]);
            statdist_z=TableAboutUsingz_statdist_z.(['useModel',num2str(useModel)]).(['nSigmaz',num2str(nSigmaz*10)]).(['nz',num2str(n_z)]);

            z_grid_temp=z_grid_J(:,21);
            statdistz_temp=statdist_z(:,21);
            if useModel==2
                fprintf(FID, 'Model 2: Rouwenhorst, nSigma=%i, $n_z$=%i \\\\ \n',nSigmaz,n_z);
            elseif useModel==21
                fprintf(FID, 'Model 2: Farmer-Toda targetting 2 moments, nSigma=%i, $n_z$=%i \\\\ \n',nSigmaz,n_z);
            elseif useModel==22
                fprintf(FID, 'Model 2: Farmer-Toda targetting 4 moments, nSigma=%i, $n_z$=%i \\\\ \n',nSigmaz,n_z);
            end
            if nSigmaz==3
                fprintf(FID,'   +-1/2/3 standard deviation has mass %8.4f / %8.4f / %8.4f \\\\ \n', sum(statdistz_temp(abs(z_grid_temp)<=z_grid_temp(end)/3)), sum(statdistz_temp(abs(z_grid_temp)<=2*z_grid_temp(end)/3)), sum(statdistz_temp(abs(z_grid_temp)<=3*z_grid_temp(end)/3)));
            elseif nSigmaz==4
                fprintf(FID,'   +-1/2/3/4 standard deviation has mass %8.4f \\\\ \n', sum(statdistz_temp(abs(z_grid_temp)<=z_grid_temp(end)/4)), sum(statdistz_temp(abs(z_grid_temp)<=2*z_grid_temp(end)/4)), sum(statdistz_temp(abs(z_grid_temp)<=3*z_grid_temp(end)/4)), sum(statdistz_temp(abs(z_grid_temp)<=4*z_grid_temp(end)/4)));
            end
        end
    end
end
fprintf(FID,'For Model 2, at age j=41, we get the following \n');
for useModel=[2,21,22]
    for n_z=nzvectable
        for nSigmaz=nSigmazvec
            z_grid_temp=z_grid_J(:,41);
            statdistz_temp=statdist_z(:,41);
            if useModel==2
                fprintf(FID, 'Model 2: Rouwenhorst, nSigma=%i, $n_z$=%i \\\\ \n',nSigmaz,n_z);
            elseif useModel==21
                fprintf(FID, 'Model 2: Farmer-Toda targetting 2 moments, nSigma=%i, $n_z$=%i \\\\ \n',nSigmaz,n_z);
            elseif useModel==22
                fprintf(FID, 'Model 2: Farmer-Toda targetting 4 moments, nSigma=%i, $n_z$=%i \\\\ \n',nSigmaz,n_z);
            end
            if nSigmaz==3
                fprintf(FID,'   +-1/2/3 standard deviation has mass %8.4f / %8.4f / %8.4f \\\\ \n', sum(statdistz_temp(abs(z_grid_temp)<=z_grid_temp(end)/3)), sum(statdistz_temp(abs(z_grid_temp)<=2*z_grid_temp(end)/3)), sum(statdistz_temp(abs(z_grid_temp)<=3*z_grid_temp(end)/3)));
            elseif nSigmaz==4
                fprintf(FID,'   +-1/2/3/4 standard deviation has mass %8.4f \\\\ \n', sum(statdistz_temp(abs(z_grid_temp)<=z_grid_temp(end)/4)), sum(statdistz_temp(abs(z_grid_temp)<=2*z_grid_temp(end)/4)), sum(statdistz_temp(abs(z_grid_temp)<=3*z_grid_temp(end)/4)), sum(statdistz_temp(abs(z_grid_temp)<=4*z_grid_temp(end)/4)));
            end
        end
    end
end
fclose(FID);


%% Create a table of pi_upsilon for models 2 and 6
FID = fopen(['./SavedOutput/EvaluateDiscretization/Table_Epi_upsilon.tex'], 'w');
useModel=6;
for n_z=nzvectable
    for nSigmaz=nSigmazvec
        Epi_upsilon_J=TableAboutpiupsilon.(['useModel',num2str(useModel)]).(['nSigmaz',num2str(nSigmaz*10)]).(['nz',num2str(n_z)]);
        fprintf(FID, 'Model 6: Farmer-Toda targetting 4 moments, nSigma=%i, $n_z$=%i \\\\ \n',nSigmaz,n_z);
        for jjc=1:5
            jj=jjvec(jjc);
            fprintf(FID, '  Epi upsilon=[%8.4f, %8.4f; %8.4f, %8.4f] at age j=%i \\\\ \n',Epi_upsilon_J(:,:,jj),jj);
        end
    end
end
useModel=2;
for n_z=nzvectable
    for nSigmaz=nSigmazvec
        Epi_upsilon_J=TableAboutpiupsilon.(['useModel',num2str(useModel)]).(['nSigmaz',num2str(nSigmaz*10)]).(['nz',num2str(n_z)]);
        fprintf(FID, 'Model 2: Rouwenhorst, nSigma=%i, $n_z$=%i \\\\ \n',nSigmaz,n_z);
        for jjc=1:5
            jj=jjvec(jjc);
            fprintf(FID, '  Epi upsilon=[%8.4f, %8.4f; %8.4f, %8.4f] at age j=%i \\\\ \n',Epi_upsilon_J(:,:,jj),jj);
        end
    end
end
fclose(FID);








