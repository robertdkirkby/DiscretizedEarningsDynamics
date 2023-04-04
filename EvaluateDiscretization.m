% EvaluateDiscretization

% z_grid_J, pi_z_J: check that simulating them gives right autocorrelation and the right gaussian-mixture innovations
% jequaloneDistz: check it is normal dist
% heatmaps already show the number of moments so no need to check those

% Note: Plotting upsilon does not seem likely to be useful. Elsewhere I plot the life-cycle profile of mean of upsilon using model.

N=10^5; % Number of simulations to use for z, to calculate autocorrelations and plot innovations (pdf & cdf)

% To use exogenous shocks that depend on age you have to add them to vfoptions and simoptions
vfoptions.z_grid_J=zupsilon_grid_J; % Note: naming of vfoptions.z_grid_J has to be exactly as is.
vfoptions.pi_z_J=pi_zupsilon_J; % Note: naming of vfoptions.z_grid_J has to be exactly as is.
vfoptions.n_e=n_epsilon;
vfoptions.e_grid=epsilon_grid; % You can use vfoptions.e_grid_J and vfoptions.pi_e_J, but not needed here
vfoptions.pi_e=pi_epsilon;

simoptions.z_grid_J=vfoptions.z_grid_J; % Note: naming of vfoptions.z_grid_J has to be exactly as is.
simoptions.pi_z_J=vfoptions.pi_z_J; % Note: naming of vfoptions.z_grid_J has to be exactly as is.
simoptions.n_e=n_epsilon;
simoptions.e_grid=epsilon_grid; % You can use vfoptions.e_grid_J and vfoptions.pi_e_J, but not needed here
simoptions.pi_e=pi_epsilon;

%%
if useModel==1
    %% kappa_j
    fig=figure(1);
    plot((1:1:Params.J)+Params.agejshifter,Params.kappa_j)
    title('Deterministic labor efficiency units')
    xlabel('Age (in Years)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model1_kappaj.pdf')
    %% epsilon
    fig=figure(2);
    x = linspace(1.5*min(epsilon_grid),1.5*max(epsilon_grid),201);
    y = normpdf(x,0,Params.sigma_epsilon);
    subplot(1,2,1); plot(x,y)
    hold on
    scatter(epsilon_grid,pi_epsilon)
    hold off
    title('i.i.d shock: normal, pdf')
    legend('Normal pdf','Discretization')
    y = normcdf(x,0,Params.sigma_epsilon);
    subplot(1,2,2); plot(x,y)
    hold on
    scatter(epsilon_grid,cumsum(pi_epsilon))
    hold off
    title('i.i.d. shock: normal, cdf')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model1_epsilon.pdf')
    
    %% z
    markovsimoptions.nsims=10^5;
    z_panel=MarkovChainSimulation_FHorz(z_grid_J,pi_z_J,jequaloneDistz,markovsimoptions);
    z_panel=z_panel(:,1:Params.Jr-1); % Drop the retirement age as is not relevant to the actual discretization of z
    % Calculate the correlation for each age: compare to Params.rho_z
    fig=figure(3);
    discretized_rho_z=nan(1,Params.Jr-1);
    for jj=2:Params.Jr-1
        discretized_rho_z(jj)=corr(z_panel(:,jj),z_panel(:,jj-1));
    end
    plot((1:1:(Params.Jr-1))+Params.agejshifter,Params.rho,(1:1:(Params.Jr-1))+Params.agejshifter,discretized_rho_z)
    title('Autocorrelation coefficient of z')
    legend('Original','Discretized')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model1_zcorr.pdf')
    % Calculate the innovations and plot them
    fig=figure(4);
    % z_innov=zeros(N,Params.Jr-1-1); % Extra minus one as first difference
    z_innov=z_panel(:,2:end)-z_panel(:,1:end-1);
    aa=floor((Params.Jr-1)/10); % Create 10 subplots per row
    bb=ceil((Params.Jr-1)/aa);  
    fig=figure(4);
    for djj=1:Params.Jr-2 % djj=1 will give transition from period 1 to 2
        x = linspace(1.5*min(z_innov(:,djj)),1.5*max(z_innov(:,djj)),201)';
        y = normpdf(x,0,Params.sigma_eta(djj));
        subplot(bb,aa,djj); plot(x,y)
        hold on
        histogram(z_innov(:,djj),n_z,'Normalization','probability')
        hold off        
    end
    legend('Normal pdf','Discretization')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model1_zinnov_pdf.pdf')
    fig=figure(5);
    for djj=1:Params.Jr-2 % djj=1 will give transition from period 1 to 2
        x = linspace(1.5*min(z_innov(:,djj)),1.5*max(z_innov(:,djj)),201)';
        y = normcdf(x,0,Params.sigma_eta(djj));
        subplot(bb,aa,djj); plot(x,y)
        hold on
        histogram(z_innov(:,djj),n_z,'Normalization','cdf')
        hold off        
    end
    legend('Normal cdf','Discretization')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model1_zinnov_cdf.pdf')
    
    % Does all of the z_grid get used? Plot the age conditional pdfs and cdfs
    z_statdist=zeros(n_z,Params.Jr);
    z_statdist(:,1)=jequaloneDistz;
    for jj=1:Params.Jr-1
        z_statdist(:,jj+1)=pi_z_J(:,:,jj)'*z_statdist(:,1);
    end
    aa=floor((Params.Jr)/10); % Create 10 subplots per row
    bb=ceil((Params.Jr)/aa);  
    fig=figure(6);
    for jj=1:Params.Jr-1
        subplot(bb,aa,jj); plot(z_grid_J(:,jj),z_statdist(:,jj))
        xlim([min(z_grid_J(:,jj)),max(z_grid_J(:,jj))])
    end
    legend('Discretized z pdf')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model1_z_pdf.pdf')
    fig=figure(7);
    for jj=1:Params.Jr-1
        subplot(bb,aa,jj); plot(z_grid_J(:,jj),cumsum(z_statdist(:,jj)))
        xlim([min(z_grid_J(:,jj)),max(z_grid_J(:,jj))])
    end
    legend('Discretized z cdf')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model1_z_cdf.pdf')
    
    % Calculate the standard deviation of the discretized z process and
    % compare to the theortical ones (which were used to create the grid for z)
    [mean_disc_z,variance_disc_z,corr_disc_z,statdist_z]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J,jequaloneDistz);
    fig=figure(8);
    subplot(2,1,1);  plot((1:1:(Params.Jr-1))+Params.agejshifter,sqrt(variance_disc_z(1:Params.Jr-1)), (1:1:(Params.Jr-1))+Params.agejshifter,otheroutputs_z.sigma_z)
    legend('Discretized','Theory')
    title('Standard deviation of z (from age to next age)')
    subplot(2,1,2); plot((1:1:Params.Jr-1)+Params.agejshifter,mean_disc_z(1:Params.Jr-1))
    ylim([-0.1,0.1])
    legend('Discretized')
    title('Mean of z (conditional on age)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model1_z_stddev.pdf')
    % Do another check of the autocorrelation coefficient
    fig=figure(9);
    plot((1:1:(Params.Jr-1))+Params.agejshifter,Params.rho,(1:1:(Params.Jr-1))+Params.agejshifter,discretized_rho_z,(2:1:(Params.Jr-1))+Params.agejshifter,corr_disc_z(2:Params.Jr-1))
    ylim([0.9,1])
    title('Autocorrelation coefficient of z')
    legend('Original','Discretized','Discretized double-check')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model1_zcorr2.pdf')
    
    % Plot the transition matrix at certain ages as 2D-pdfs
    fig=figure(10);
    subplot(2,2,1); surf(pi_z_J(:,:,1))
    title('Transition matrix pdf for j=1 (age 25)')
    subplot(2,2,2); surf(pi_z_J(:,:,11))
    title('Transition matrix pdf for j=11 (age 35)')
    subplot(2,2,3); surf(pi_z_J(:,:,21))
    title('Transition matrix pdf for j=21 (age 45)')
    subplot(2,2,4); surf(pi_z_J(:,:,31))
    title('Transition matrix pdf for j=31 (age 55)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model1_surf_pi_z_J.pdf')
    % Plot the transition matrix at certain ages as cdfs for each z
    fig=figure(11);
    subplot(2,2,1); plot(cumsum(pi_z_J(:,:,1),2)')
    title('Transition matrix cdf for j=1 (age 25)')
    subplot(2,2,2); plot(cumsum(pi_z_J(:,:,11),2)')
    title('Transition matrix cdf for j=11 (age 35)')
    subplot(2,2,3); plot(cumsum(pi_z_J(:,:,21),2)')
    title('Transition matrix cdf for j=21 (age 45)')
    subplot(2,2,4); plot(cumsum(pi_z_J(:,:,31),2)')
    title('Transition matrix cdf for j=31 (age 55)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model1_cdf_pi_z_J.pdf')
    
    % Plotting upsilon does not seem likely to be useful. Elsewhere I plot the life-cycle profile of mean of upsilon using model.
    
    %% alpha
    fig=figure(8);
    x = linspace(1.5*min(Params.alpha),1.5*max(Params.alpha),201); % Params.alpha is alpha_grid
    y = normpdf(x,0,Params.sigma_alpha);
    subplot(1,2,1); plot(x,y)
    hold on
    scatter(Params.alpha,Params.statdist_alpha)
    hold off
    title('Fixed effect: i.i.d. normal, pdf')
    legend('Normal','Discretization')
    y = normcdf(x,0,Params.sigma_alpha);
    subplot(1,2,2); plot(x,y)
    hold on
    scatter(Params.alpha,cumsum(Params.statdist_alpha))
    hold off
    title('Fixed effect: i.i.d. normal, cdf')
%     legend('Normal','Discretization')
    saveas(fig,['./SavedOutput/EvaluateDiscretization/Model1_alpha_nSigma',num2str(nSigma_alpha),'.pdf'])

elseif useModel==2
    %% kappa_j
    fig=figure(1);
    plot((1:1:Params.J)+Params.agejshifter,Params.kappa_j)
    title('Deterministic labor efficiency units')
    xlabel('Age (in Years)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model2_kappaj.pdf')
    %% epsilon
    fig=figure(2);
    x = linspace(1.5*min(epsilon_grid),1.5*max(epsilon_grid),201);
    y = normpdf(x,0,Params.sigma_epsilon);
    subplot(1,2,1); plot(x,y)
    hold on
    scatter(epsilon_grid,pi_epsilon)
    hold off
    title('i.i.d shock: normal, pdf')
    legend('Normal pdf','Discretization')
    y = normcdf(x,0,Params.sigma_epsilon);
    subplot(1,2,2); plot(x,y)
    hold on
    scatter(epsilon_grid,cumsum(pi_epsilon))
    hold off
    title('i.i.d. shock: normal, cdf')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model2_epsilon.pdf')
    
    %% z
    markovsimoptions.nsims=10^5;
    z_panel=MarkovChainSimulation_FHorz(z_grid_J,pi_z_J,jequaloneDistz,markovsimoptions);
    z_panel=z_panel(:,1:Params.Jr-1); % Drop the retirement age as is not relevant to the actual discretization of z
    % Calculate the correlation for each age: compare to Params.rho_z
    fig=figure(3);
    discretized_rho_z=nan(1,Params.Jr-1);
    for jj=2:Params.Jr-1
        discretized_rho_z(jj)=corr(z_panel(:,jj),z_panel(:,jj-1));
    end
    plot((1:1:(Params.Jr-1))+Params.agejshifter,Params.rho,(1:1:(Params.Jr-1))+Params.agejshifter,discretized_rho_z)
    title('Autocorrelation coefficient of z')
    legend('Original','Discretized')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model2_zcorr.pdf')
    % Calculate the innovations and plot them
    fig=figure(4);
    % z_innov=zeros(N,Params.Jr-1-1); % Extra minus one as first difference
    z_innov=z_panel(:,2:end)-z_panel(:,1:end-1);
    aa=floor((Params.Jr-1)/10); % Create 10 subplots per row
    bb=ceil((Params.Jr-1)/aa);  
    fig=figure(4);
    for djj=1:Params.Jr-2 % djj=1 will give transition from period 1 to 2
        x = linspace(1.5*min(z_innov(:,djj)),1.5*max(z_innov(:,djj)),201)';
        y = normpdf(x,0,Params.sigma_eta(djj));
        subplot(bb,aa,djj); plot(x,y)
        hold on
        histogram(z_innov(:,djj),n_z,'Normalization','probability')
        hold off        
    end
    legend('Normal pdf','Discretization')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model2_zinnov_pdf.pdf')
    fig=figure(5);
    for djj=1:Params.Jr-2 % djj=1 will give transition from period 1 to 2
        x = linspace(1.5*min(z_innov(:,djj)),1.5*max(z_innov(:,djj)),201)';
        y = normcdf(x,0,Params.sigma_eta(djj));
        subplot(bb,aa,djj); plot(x,y)
        hold on
        histogram(z_innov(:,djj),n_z,'Normalization','cdf')
        hold off        
    end
    legend('Normal cdf','Discretization')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model2_zinnov_cdf.pdf')
    
    % Does all of the z_grid get used? Plot the age conditional pdfs and cdfs
    z_statdist=zeros(n_z,Params.Jr);
    z_statdist(:,1)=jequaloneDistz;
    for jj=1:Params.Jr-1
        z_statdist(:,jj+1)=pi_z_J(:,:,jj)'*z_statdist(:,1);
    end
    aa=floor((Params.Jr)/10); % Create 10 subplots per row
    bb=ceil((Params.Jr)/aa);  
    fig=figure(6);
    for jj=1:Params.Jr-1
        subplot(bb,aa,jj); plot(z_grid_J(:,jj),z_statdist(:,jj))
        xlim([min(z_grid_J(:,jj)),max(z_grid_J(:,jj))])
    end
    legend('Discretized z pdf')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model2_z_pdf.pdf')
    fig=figure(7);
    for jj=1:Params.Jr-1
        subplot(bb,aa,jj); plot(z_grid_J(:,jj),cumsum(z_statdist(:,jj)))
        xlim([min(z_grid_J(:,jj)),max(z_grid_J(:,jj))])
    end
    legend('Discretized z cdf')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model2_z_cdf.pdf')
    
    % Calculate the standard deviation of the discretized z process and
    % compare to the theortical ones (which were used to create the grid for z)
    [mean_disc_z,variance_disc_z,corr_disc_z,statdist_z]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J,jequaloneDistz);
    fig=figure(8);
    subplot(2,1,1);  plot((1:1:(Params.Jr-1))+Params.agejshifter,sqrt(variance_disc_z(1:Params.Jr-1)), (1:1:(Params.Jr-1))+Params.agejshifter,otheroutputs_z.sigma_z)
    legend('Discretized','Theory')
    title('Standard deviation of z (from age to next age)')
    subplot(2,1,2); plot((1:1:Params.Jr-1)+Params.agejshifter,mean_disc_z(1:Params.Jr-1))
    ylim([-0.1,0.1])
    legend('Discretized')
    title('Mean of z (conditional on age)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model2_z_stddev.pdf')
    % Do another check of the autocorrelation coefficient
    fig=figure(9);
    plot((1:1:(Params.Jr-1))+Params.agejshifter,Params.rho,(1:1:(Params.Jr-1))+Params.agejshifter,discretized_rho_z,(2:1:(Params.Jr-1))+Params.agejshifter,corr_disc_z(2:Params.Jr-1))
    ylim([0.9,1])
    title('Autocorrelation coefficient of z')
    legend('Original','Discretized','Discretized double-check')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model2_zcorr2.pdf')
    
    % Plot the transition matrix at certain ages as 2D-pdfs
    fig=figure(10);
    subplot(2,2,1); surf(pi_z_J(:,:,1))
    title('Transition matrix pdf for j=1 (age 25)')
    subplot(2,2,2); surf(pi_z_J(:,:,11))
    title('Transition matrix pdf for j=11 (age 35)')
    subplot(2,2,3); surf(pi_z_J(:,:,21))
    title('Transition matrix pdf for j=21 (age 45)')
    subplot(2,2,4); surf(pi_z_J(:,:,31))
    title('Transition matrix pdf for j=31 (age 55)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model2_surf_pi_z_J.pdf')
    % Plot the transition matrix at certain ages as cdfs for each z
    fig=figure(11);
    subplot(2,2,1); plot(cumsum(pi_z_J(:,:,1),2)')
    title('Transition matrix cdf for j=1 (age 25)')
    subplot(2,2,2); plot(cumsum(pi_z_J(:,:,11),2)')
    title('Transition matrix cdf for j=11 (age 35)')
    subplot(2,2,3); plot(cumsum(pi_z_J(:,:,21),2)')
    title('Transition matrix cdf for j=21 (age 45)')
    subplot(2,2,4); plot(cumsum(pi_z_J(:,:,31),2)')
    title('Transition matrix cdf for j=31 (age 55)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model2_cdf_pi_z_J.pdf')
    
    % Plotting upsilon does not seem likely to be useful. Elsewhere I plot the life-cycle profile of mean of upsilon using model.
    
    %% alpha
    fig=figure(8);
    x = linspace(1.5*min(Params.alpha),1.5*max(Params.alpha),201); % Params.alpha is alpha_grid
    y = normpdf(x,0,Params.sigma_alpha);
    subplot(1,2,1); plot(x,y)
    hold on
    scatter(Params.alpha,Params.statdist_alpha)
    hold off
    title('Fixed effect: i.i.d. normal, pdf')
    legend('Normal','Discretization')
    y = normcdf(x,0,Params.sigma_alpha);
    subplot(1,2,2); plot(x,y)
    hold on
    scatter(Params.alpha,cumsum(Params.statdist_alpha))
    hold off
    title('Fixed effect: i.i.d. normal, cdf')
%     legend('Normal','Discretization')
    saveas(fig,['./SavedOutput/EvaluateDiscretization/Model2_alpha_nSigma',num2str(nSigma_alpha),'.pdf'])

    
elseif useModel==3
    %% kappa_j
    fig=figure(1);
    plot((1:1:Params.J)+Params.agejshifter,Params.kappa_j)
    title('Deterministic labor efficiency units')
    xlabel('Age (in Years)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model3_kappaj.pdf')
    %% epsilon
    gm = gmdistribution(Params.mu_epsilon,shiftdim(Params.sigma_epsilon,-2),Params.p_epsilon); % mu and sigma need to be row vectors
    fig=figure(2);
    x = linspace(1.5*min(epsilon_grid),1.5*max(epsilon_grid),201)';
    ypdf = pdf(gm,x);
    subplot(1,2,2); plot(x,ypdf)
    hold on
    scatter(epsilon_grid,pi_epsilon)
    hold off
    title('i.i.d. shock: gaussian-mixture, pdf')
    legend('Gaussian-mixture','Discretization')
    ycdf = cdf(gm,x);
    subplot(1,2,1); plot(x,ycdf)
    hold on
    scatter(epsilon_grid,cumsum(pi_epsilon))
    hold off
    title('i.i.d. shock: gaussian-mixture, cdf')
%     legend('Gaussian-mixture','Discretization')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model3_epsilon.pdf')
    
    %% z
    markovsimoptions.nsims=10^5;
    z_panel=MarkovChainSimulation_FHorz(z_grid_J,pi_z_J,jequaloneDistz,markovsimoptions);
    z_panel=z_panel(:,1:Params.Jr-1); % Drop the retirement age as is not relevant to the actual discretization of z
    % Calculate the correlation for each age: compare to Params.rho_z
    fig=figure(3);
    discretized_rho_z=nan(1,Params.Jr-1);
    for jj=2:Params.Jr-1
        discretized_rho_z(jj)=corr(z_panel(:,jj),z_panel(:,jj-1));
    end
    plot((1:1:(Params.Jr-1))+Params.agejshifter,Params.rho,(1:1:(Params.Jr-1))+Params.agejshifter,discretized_rho_z)
    title('Autocorrelation coefficient of z')
    legend('Original','Discretized')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model3_zcorr.pdf')
    % Calculate the innovations and plot them
    fig=figure(4);
    % z_innov=zeros(N,Params.Jr-1-1); % Extra minus one as first difference
    z_innov=z_panel(:,2:end)-z_panel(:,1:end-1);
    aa=floor((Params.Jr-1)/10); % Create 10 subplots per row
    bb=ceil((Params.Jr-1)/aa);  
    fig=figure(4);
    for djj=1:Params.Jr-2 % djj=1 will give transition from period 1 to 2
        x = linspace(1.5*min(z_innov(:,djj)),1.5*max(z_innov(:,djj)),201)';
        % Params.mu_eta is just the first, the second is determined to make E[eta]=0.
        mu_eta=[Params.mu_eta(djj); -(Params.mu_eta(djj).*Params.mixprobs_eta(1,djj))/Params.mixprobs_eta(2,djj)]; % Simple rearrangement of mu_i(:,jj).*mixprob_i(:,jj)=0, which is the requirement that mean of gaussian-mixture innovations=0
        gm = gmdistribution(mu_eta,shiftdim(Params.sigma_eta(:,djj),-2),Params.mixprobs_eta(:,djj));
        y = pdf(gm,x);
        subplot(bb,aa,djj); plot(x,y)
        hold on
        histogram(z_innov(:,djj),n_z,'Normalization','probability')
        hold off        
    end
    legend('Gaussian-mixture pdf','Discretization')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model3_zinnov_pdf.pdf')
    fig=figure(5);
    for djj=1:Params.Jr-2 % djj=1 will give transition from period 1 to 2
        x = linspace(1.5*min(z_innov(:,djj)),1.5*max(z_innov(:,djj)),201)';
        % Params.mu_eta is just the first, the second is determined to make E[eta]=0.
        mu_eta=[Params.mu_eta(djj); -(Params.mu_eta(djj).*Params.mixprobs_eta(1,djj))/Params.mixprobs_eta(2,djj)]; % Simple rearrangement of mu_i(:,jj).*mixprob_i(:,jj)=0, which is the requirement that mean of gaussian-mixture innovations=0
        gm = gmdistribution(mu_eta,shiftdim(Params.sigma_eta(:,djj),-2),Params.mixprobs_eta(:,djj));
        y = cdf(gm,x);
        subplot(bb,aa,djj); plot(x,y)
        hold on
        histogram(z_innov(:,djj),n_z,'Normalization','cdf')
        hold off        
    end
    legend('Gaussian-mixture cdf','Discretization')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model3_zinnov_cdf.pdf')
    
    % Does all of the z_grid get used? Plot the age conditional pdfs and cdfs
    z_statdist=zeros(n_z,Params.Jr);
    z_statdist(:,1)=jequaloneDistz;
    for jj=1:Params.Jr-1
        z_statdist(:,jj+1)=pi_z_J(:,:,jj)'*z_statdist(:,1);
    end
    aa=floor((Params.Jr)/10); % Create 10 subplots per row
    bb=ceil((Params.Jr)/aa);  
    fig=figure(6);
    for jj=1:Params.Jr-1
        subplot(bb,aa,jj); plot(z_grid_J(:,jj),z_statdist(:,jj))
        xlim([min(z_grid_J(:,jj)),max(z_grid_J(:,jj))])
    end
    legend('Discretized z pdf')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model3_z_pdf.pdf')
    fig=figure(7);
    for jj=1:Params.Jr-1
        subplot(bb,aa,jj); plot(z_grid_J(:,jj),cumsum(z_statdist(:,jj)))
        xlim([min(z_grid_J(:,jj)),max(z_grid_J(:,jj))])
    end
    legend('Discretized z cdf')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model3_z_cdf.pdf')
    
    % Calculate the standard deviation of the discretized z process and
    % compare to the theortical ones (which were used to create the grid for z)
    [mean_disc_z,variance_disc_z,corr_disc_z,statdist_z]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J,jequaloneDistz);
    fig=figure(8);
    subplot(2,1,1);  plot((1:1:(Params.Jr-1))+Params.agejshifter,sqrt(variance_disc_z(1:Params.Jr-1)), (1:1:(Params.Jr-1))+Params.agejshifter,otheroutputs_z.sigma_z)
    legend('Discretized','Theory')
    title('Standard deviation of z (from age to next age)')
    subplot(2,1,2); plot((1:1:Params.Jr-1)+Params.agejshifter,mean_disc_z(1:Params.Jr-1))
    ylim([-0.1,0.1])
    legend('Discretized')
    title('Mean of z (conditional on age)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model3_z_stddev.pdf')
    % Do another check of the autocorrelation coefficient
    fig=figure(9);
    plot((1:1:(Params.Jr-1))+Params.agejshifter,Params.rho,(1:1:(Params.Jr-1))+Params.agejshifter,discretized_rho_z,(2:1:(Params.Jr-1))+Params.agejshifter,corr_disc_z(2:Params.Jr-1))
    ylim([0.9,1])
    title('Autocorrelation coefficient of z')
    legend('Original','Discretized','Discretized double-check')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model3_zcorr2.pdf')
    
    % Plot the transition matrix at certain ages as 2D-pdfs
    fig=figure(10);
    subplot(2,2,1); surf(pi_z_J(:,:,1))
    title('Transition matrix pdf for j=1 (age 25)')
    subplot(2,2,2); surf(pi_z_J(:,:,11))
    title('Transition matrix pdf for j=11 (age 35)')
    subplot(2,2,3); surf(pi_z_J(:,:,21))
    title('Transition matrix pdf for j=21 (age 45)')
    subplot(2,2,4); surf(pi_z_J(:,:,31))
    title('Transition matrix pdf for j=31 (age 55)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model3_surf_pi_z_J.pdf')
    % Plot the transition matrix at certain ages as cdfs for each z
    fig=figure(11);
    subplot(2,2,1); plot(cumsum(pi_z_J(:,:,1),2)')
    title('Transition matrix cdf for j=1 (age 25)')
    subplot(2,2,2); plot(cumsum(pi_z_J(:,:,11),2)')
    title('Transition matrix cdf for j=11 (age 35)')
    subplot(2,2,3); plot(cumsum(pi_z_J(:,:,21),2)')
    title('Transition matrix cdf for j=21 (age 45)')
    subplot(2,2,4); plot(cumsum(pi_z_J(:,:,31),2)')
    title('Transition matrix cdf for j=31 (age 55)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model3_cdf_pi_z_J.pdf')
    
    % Plotting upsilon does not seem likely to be useful. 
    % Elsewhere I plot the life-cycle profile of mean of upsilon using model.
    
    %% alpha
    fig=figure(8);
    x = linspace(1.5*min(Params.alpha),1.5*max(Params.alpha),201); % Params.alpha is alpha_grid
    y = normpdf(x,0,Params.sigma_alpha);
    subplot(1,2,1); plot(x,y)
    hold on
    scatter(Params.alpha,Params.statdist_alpha)
    hold off
    title('Fixed effect: i.i.d. normal')
    legend('Normal pdf','Discretization')
    y = normcdf(x,0,Params.sigma_alpha);
    subplot(1,2,2); plot(x,y)
    hold on
    scatter(Params.alpha,cumsum(Params.statdist_alpha))
    hold off
    title('Fixed effect: i.i.d. normal')
%     legend('Normal cdf','Discretization')
    saveas(fig,['./SavedOutput/EvaluateDiscretization/Model3_alpha_nSigma',num2str(nSigma_alpha),'.pdf'])
    
elseif useModel==4
    %% kappa_j
    fig=figure(1);
    plot((1:1:Params.J)+Params.agejshifter,Params.kappa_j)
    title('Deterministic labor efficiency units')
    xlabel('Age (in Years)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model4_kappaj.pdf')
    %% epsilon
    gm = gmdistribution(Params.mu_epsilon,shiftdim(Params.sigma_epsilon,-2),Params.p_epsilon); % mu and sigma need to be row vectors
    fig=figure(2);
    x = linspace(1.5*min(epsilon_grid),1.5*max(epsilon_grid),201)';
    ypdf = pdf(gm,x);
    subplot(1,2,2); plot(x,ypdf)
    hold on
    scatter(epsilon_grid,pi_epsilon)
    hold off
    title('i.i.d. shock: gaussian-mixture, pdf')
    legend('Gaussian-mixture','Discretization')
    ycdf = cdf(gm,x);
    subplot(1,2,1); plot(x,ycdf)
    hold on
    scatter(epsilon_grid,cumsum(pi_epsilon))
    hold off
    title('i.i.d. shock: gaussian-mixture, cdf')
%     legend('Gaussian-mixture','Discretization')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model4_epsilon.pdf')
    %% z
    % NOTE: FOLLOWING NOT YET CORRECT AS DEPENDS ON ZLAG
    markovsimoptions.nsims=10^5;
    z_panel=MarkovChainSimulation_FHorz(z_grid_J,pi_z_J,jequaloneDistz,markovsimoptions);
    z_panel=z_panel(:,1:Params.Jr-1); % Drop the retirement age as is not relevant to the actual discretization of z
    % Calculate the correlation for each age: compare to Params.rho_z
    fig=figure(3);
    discretized_rho_z=nan(1,Params.Jr-1);
    for jj=2:Params.Jr-1
        discretized_rho_z(jj)=corr(z_panel(:,jj),z_panel(:,jj-1));
    end
    plot((1:1:(Params.Jr-1))+Params.agejshifter,Params.rho,(1:1:(Params.Jr-1))+Params.agejshifter,discretized_rho_z)
    title('Autocorrelation coefficient of z')
    legend('Original','Discretized')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model4_zcorr.pdf')
    % NOTE: FOLLOWING NOT YET CORRECT AS DEPENDS ON ZLAG
%     % Calculate the innovations and plot them
%     fig=figure(4)
%     % z_innov=zeros(N,Params.Jr-1-1); % Extra minus one as first difference
%     z_innov=z_panel(:,2:end)-z_panel(:,1:end-1);
%     aa=floor((Params.Jr-1)/10); % Create 10 subplots per row
%     bb=ceil((Params.Jr-1)/aa);    
%     fig=figure(4);
%     for djj=1:Params.Jr-2 % djj=1 will give transition from period 1 to 2
%         x = linspace(1.5*min(z_innov(:,djj)),1.5*max(z_innov(:,djj)),201)';
%         % Params.mu_eta is just the first, the second is determined to make E[eta]=0.
%         mu_eta=[Params.mu_eta(djj); -(Params.mu_eta(djj).*Params.mixprobs_eta(1,djj))/Params.mixprobs_eta(2,djj)]; % Simple rearrangement of mu_i(:,jj).*mixprob_i(:,jj)=0, which is the requirement that mean of gaussian-mixture innovations=0
%         gm = gmdistribution(mu_eta,shiftdim(Params.sigma_eta(:,djj),-2),Params.mixprobs_eta(:,djj));
%         y = pdf(gm,x);
%         subplot(bb,aa,djj); plot(x,y)
%         hold on
%         histogram(z_innov(:,djj),n_z,'Normalization','probability')
%         hold off        
%     end
%     legend('Gaussian-mixture pdf','Discretization')
%     saveas(fig,'./SavedOutput/EvaluateDiscretization/Model4_zinnov_pdf.pdf')
%     fig=figure(5);
%     for djj=1:Params.Jr-2 % djj=1 will give transition from period 1 to 2
%         x = linspace(1.5*min(z_innov(:,djj)),1.5*max(z_innov(:,djj)),201)';
%         % Params.mu_eta is just the first, the second is determined to make E[eta]=0.
%         mu_eta=[Params.mu_eta(djj); -(Params.mu_eta(djj).*Params.mixprobs_eta(1,djj))/Params.mixprobs_eta(2,djj)]; % Simple rearrangement of mu_i(:,jj).*mixprob_i(:,jj)=0, which is the requirement that mean of gaussian-mixture innovations=0
%         gm = gmdistribution(mu_eta,shiftdim(Params.sigma_eta(:,djj),-2),Params.mixprobs_eta(:,djj));
%         y = cdf(gm,x);
%         subplot(bb,aa,djj); plot(x,y)
%         hold on
%         histogram(z_innov(:,djj),n_z,'Normalization','cdf')
%         hold off        
%     end
%     legend('Gaussian-mixture cdf','Discretization')
%     saveas(fig,'./SavedOutput/EvaluateDiscretization/Model4_zinnov_cdf.pdf')
    
    % Does all of the z_grid get used? Plot the age conditional pdfs and cdfs
    z_statdist=zeros(n_z,Params.Jr);
    z_statdist(:,1)=jequaloneDistz;
    for jj=1:Params.Jr-1
        z_statdist(:,jj+1)=pi_z_J(:,:,jj)'*z_statdist(:,1);
    end
    aa=floor((Params.Jr)/10); % Create 10 subplots per row
    bb=ceil((Params.Jr)/aa);  
    fig=figure(6);
    for jj=1:Params.Jr-1
        subplot(bb,aa,jj); plot(z_grid_J(:,jj),z_statdist(:,jj))
        xlim([min(z_grid_J(:,jj)),max(z_grid_J(:,jj))])
    end
    legend('Discretized z pdf')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model4_z_pdf.pdf')
    fig=figure(7);
    for jj=1:Params.Jr-1
        subplot(bb,aa,jj); plot(z_grid_J(:,jj),cumsum(z_statdist(:,jj)))
        xlim([min(z_grid_J(:,jj)),max(z_grid_J(:,jj))])
    end
    legend('Discretized z cdf')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model4_z_cdf.pdf')
    
    
    % Calculate the standard deviation of the discretized z process and
    % compare to the theortical ones (which were used to create the grid for z)
    [mean_disc_z,variance_disc_z,corr_disc_z,statdist_z]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J,jequaloneDistz);
    fig=figure(8);
    subplot(2,1,1);  plot((1:1:(Params.Jr-1))+Params.agejshifter,sqrt(variance_disc_z(1:Params.Jr-1)), (1:1:(Params.Jr-1))+Params.agejshifter,otheroutputs_z.sigma_z)
    legend('Discretized','Theory')
    title('Standard deviation of z (from age to next age)')
    subplot(2,1,2); plot((1:1:Params.Jr-1)+Params.agejshifter,mean_disc_z(1:Params.Jr-1))
    ylim([-0.1,0.1])
    legend('Discretized')
    title('Mean of z (conditional on age)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model4_z_stddev.pdf')
    % Do another check of the autocorrelation coefficient
    fig=figure(9);
    plot((1:1:(Params.Jr-1))+Params.agejshifter,Params.rho,(1:1:(Params.Jr-1))+Params.agejshifter,discretized_rho_z,(2:1:(Params.Jr-1))+Params.agejshifter,corr_disc_z(2:Params.Jr-1))
    ylim([0.9,1])
    title('Autocorrelation coefficient of z')
    legend('Original','Discretized','Discretized double-check')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model4_zcorr2.pdf')
    
    % Plot the transition matrix at certain ages as 2D-pdfs
    fig=figure(10);
    subplot(2,2,1); surf(pi_z_J(:,:,1))
    title('Transition matrix pdf for j=1 (age 25)')
    subplot(2,2,2); surf(pi_z_J(:,:,11))
    title('Transition matrix pdf for j=11 (age 35)')
    subplot(2,2,3); surf(pi_z_J(:,:,21))
    title('Transition matrix pdf for j=21 (age 45)')
    subplot(2,2,4); surf(pi_z_J(:,:,31))
    title('Transition matrix pdf for j=31 (age 55)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model4_surf_pi_z_J.pdf')
    % Plot the transition matrix at certain ages as cdfs for each z
    fig=figure(11);
    subplot(2,2,1); plot(cumsum(pi_z_J(:,:,1),2)')
    title('Transition matrix cdf for j=1 (age 25)')
    subplot(2,2,2); plot(cumsum(pi_z_J(:,:,11),2)')
    title('Transition matrix cdf for j=11 (age 35)')
    subplot(2,2,3); plot(cumsum(pi_z_J(:,:,21),2)')
    title('Transition matrix cdf for j=21 (age 45)')
    subplot(2,2,4); plot(cumsum(pi_z_J(:,:,31),2)')
    title('Transition matrix cdf for j=31 (age 55)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model4_cdf_pi_z_J.pdf')

    % Plotting upsilon does not seem likely to be useful. Elsewhere I plot the life-cycle profile of mean of upsilon using model.
        
    %% alpha
    fig=figure(8);
    x = linspace(1.5*min(Params.alpha),1.5*max(Params.alpha),201); % Params.alpha is alpha_grid
    y = normpdf(x,0,Params.sigma_alpha);
    subplot(1,2,1); plot(x,y)
    hold on
    scatter(Params.alpha,Params.statdist_alpha)
    hold off
    legend('Normal','Discretization')
    title('Fixed effect: i.i.d. normal, pdf')
    y = normcdf(x,0,Params.sigma_alpha);
    subplot(1,2,2); plot(x,y)
    hold on
    scatter(Params.alpha,cumsum(Params.statdist_alpha))
    hold off
    title('Fixed effect: i.i.d. normal, cdf')
%     legend('Normal','Discretization')
    saveas(fig,['./SavedOutput/EvaluateDiscretization/Model4_alpha_nSigma',num2str(nSigma_alpha),'.pdf'])


elseif useModel==5
    %% kappa_j
    fig=figure(1);
    plot((1:1:Params.J)+Params.agejshifter,Params.kappa_j)
    title('Deterministic labor efficiency units')
    xlabel('Age (in Years)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model5_kappaj.pdf')
    %% epsilon
    gm = gmdistribution(Params.mu_epsilon,shiftdim(Params.sigma_epsilon,-2),Params.p_epsilon); % mu and sigma need to be row vectors
    fig=figure(2);
    x = linspace(1.5*min(epsilon_grid),1.5*max(epsilon_grid),201)';
    ypdf = pdf(gm,x);
    subplot(1,2,2); plot(x,ypdf)
    hold on
    scatter(epsilon_grid,pi_epsilon)
    hold off
    title('i.i.d. shock: gaussian-mixture, pdf')
    legend('Gaussian-mixture','Discretization')
    ycdf = cdf(gm,x);
    subplot(1,2,1); plot(x,ycdf)
    hold on
    scatter(epsilon_grid,cumsum(pi_epsilon))
    hold off
    title('i.i.d. shock: gaussian-mixture, cdf')
%     legend('Gaussian-mixture','Discretization')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model5_epsilon.pdf')
    
    %% z
    markovsimoptions.nsims=10^5;
    z_panel=MarkovChainSimulation_FHorz(z_grid_J,pi_z_J,jequaloneDistz,markovsimoptions);
    z_panel=z_panel(:,1:Params.Jr-1); % Drop the retirement age as is not relevant to the actual discretization of z
    % Calculate the correlation for each age: compare to Params.rho_z
    fig=figure(3);
    discretized_rho_z=nan(1,Params.Jr-1);
    for jj=2:Params.Jr-1
        discretized_rho_z(jj)=corr(z_panel(:,jj),z_panel(:,jj-1));
    end
    plot((1:1:(Params.Jr-1))+Params.agejshifter,Params.rho,(1:1:(Params.Jr-1))+Params.agejshifter,discretized_rho_z)
    title('Autocorrelation coefficient of z')
    legend('Original','Discretized')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model5_zcorr.pdf')
    % Calculate the innovations and plot them
    fig=figure(4);
    % z_innov=zeros(N,Params.Jr-1-1); % Extra minus one as first difference
    z_innov=z_panel(:,2:end)-z_panel(:,1:end-1);
    aa=floor((Params.Jr-1)/10); % Create 10 subplots per row
    bb=ceil((Params.Jr-1)/aa);  
    fig=figure(4);
    for djj=1:Params.Jr-2 % djj=1 will give transition from period 1 to 2
        x = linspace(1.5*min(z_innov(:,djj)),1.5*max(z_innov(:,djj)),201)';
        % Params.mu_eta is just the first, the second is determined to make E[eta]=0.
        mu_eta=[Params.mu_eta(djj); -(Params.mu_eta(djj).*Params.mixprobs_eta(1,djj))/Params.mixprobs_eta(2,djj)]; % Simple rearrangement of mu_i(:,jj).*mixprob_i(:,jj)=0, which is the requirement that mean of gaussian-mixture innovations=0
        gm = gmdistribution(mu_eta,shiftdim(Params.sigma_eta(:,djj),-2),Params.mixprobs_eta(:,djj));
        y = pdf(gm,x);
        subplot(bb,aa,djj); plot(x,y)
        hold on
        histogram(z_innov(:,djj),n_z,'Normalization','probability')
        hold off        
    end
    legend('Gaussian-mixture pdf','Discretization')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model5_zinnov_pdf.pdf')
    fig=figure(5);
    for djj=1:Params.Jr-2 % djj=1 will give transition from period 1 to 2
        x = linspace(1.5*min(z_innov(:,djj)),1.5*max(z_innov(:,djj)),201)';
        % Params.mu_eta is just the first, the second is determined to make E[eta]=0.
        mu_eta=[Params.mu_eta(djj); -(Params.mu_eta(djj).*Params.mixprobs_eta(1,djj))/Params.mixprobs_eta(2,djj)]; % Simple rearrangement of mu_i(:,jj).*mixprob_i(:,jj)=0, which is the requirement that mean of gaussian-mixture innovations=0
        gm = gmdistribution(mu_eta,shiftdim(Params.sigma_eta(:,djj),-2),Params.mixprobs_eta(:,djj));
        y = cdf(gm,x);
        subplot(bb,aa,djj); plot(x,y)
        hold on
        histogram(z_innov(:,djj),n_z,'Normalization','cdf')
        hold off        
    end
    legend('Gaussian-mixture cdf','Discretization')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model5_zinnov_cdf.pdf')
    
    % Does all of the z_grid get used? Plot the age conditional pdfs and cdfs
    z_statdist=zeros(n_z,Params.Jr);
    z_statdist(:,1)=jequaloneDistz;
    for jj=1:Params.Jr-1
        z_statdist(:,jj+1)=pi_z_J(:,:,jj)'*z_statdist(:,1);
    end
    aa=floor((Params.Jr)/10); % Create 10 subplots per row
    bb=ceil((Params.Jr)/aa);  
    fig=figure(6);
    for jj=1:Params.Jr-1
        subplot(bb,aa,jj); plot(z_grid_J(:,jj),z_statdist(:,jj))
        xlim([min(z_grid_J(:,jj)),max(z_grid_J(:,jj))])
    end
    legend('Discretized z pdf')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model5_z_pdf.pdf')
    fig=figure(7);
    for jj=1:Params.Jr-1
        subplot(bb,aa,jj); plot(z_grid_J(:,jj),cumsum(z_statdist(:,jj)))
        xlim([min(z_grid_J(:,jj)),max(z_grid_J(:,jj))])
    end
    legend('Discretized z cdf')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model5_z_cdf.pdf')
    
    % Calculate the standard deviation of the discretized z process and
    % compare to the theortical ones (which were used to create the grid for z)
    [mean_disc_z,variance_disc_z,corr_disc_z,statdist_z]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J,jequaloneDistz);
    fig=figure(8);
    subplot(2,1,1);  plot((1:1:(Params.Jr-1))+Params.agejshifter,sqrt(variance_disc_z(1:Params.Jr-1)), (1:1:(Params.Jr-1))+Params.agejshifter,otheroutputs_z.sigma_z)
    legend('Discretized','Theory')
    title('Standard deviation of z (from age to next age)')
    subplot(2,1,2); plot((1:1:Params.Jr-1)+Params.agejshifter,mean_disc_z(1:Params.Jr-1))
    ylim([-0.1,0.1])
    legend('Discretized')
    title('Mean of z (conditional on age)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model5_z_stddev.pdf')
    % Do another check of the autocorrelation coefficient
    fig=figure(9);
    plot((1:1:(Params.Jr-1))+Params.agejshifter,Params.rho,(1:1:(Params.Jr-1))+Params.agejshifter,discretized_rho_z,(2:1:(Params.Jr-1))+Params.agejshifter,corr_disc_z(2:Params.Jr-1))
    ylim([0.9,1])
    title('Autocorrelation coefficient of z')
    legend('Original','Discretized','Discretized double-check')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model5_zcorr2.pdf')
    
    % Plot the transition matrix at certain ages as 2D-pdfs
    fig=figure(10);
    subplot(2,2,1); surf(pi_z_J(:,:,1))
    title('Transition matrix pdf for j=1 (age 25)')
    subplot(2,2,2); surf(pi_z_J(:,:,11))
    title('Transition matrix pdf for j=11 (age 35)')
    subplot(2,2,3); surf(pi_z_J(:,:,21))
    title('Transition matrix pdf for j=21 (age 45)')
    subplot(2,2,4); surf(pi_z_J(:,:,31))
    title('Transition matrix pdf for j=31 (age 55)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model5_surf_pi_z_J.pdf')
    % Plot the transition matrix at certain ages as cdfs for each z
    fig=figure(11);
    subplot(2,2,1); plot(cumsum(pi_z_J(:,:,1),2)')
    title('Transition matrix cdf for j=1 (age 25)')
    subplot(2,2,2); plot(cumsum(pi_z_J(:,:,11),2)')
    title('Transition matrix cdf for j=11 (age 35)')
    subplot(2,2,3); plot(cumsum(pi_z_J(:,:,21),2)')
    title('Transition matrix cdf for j=21 (age 45)')
    subplot(2,2,4); plot(cumsum(pi_z_J(:,:,31),2)')
    title('Transition matrix cdf for j=31 (age 55)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model5_cdf_pi_z_J.pdf')
    
    % Plotting upsilon does not seem likely to be useful. 
    % Elsewhere I plot the life-cycle profile of mean of upsilon using model.
    
    %% alpha
    fig=figure(8);
    x = linspace(1.5*min(Params.alpha),1.5*max(Params.alpha),201); % Params.alpha is alpha_grid
    y = normpdf(x,0,Params.sigma_alpha);
    subplot(1,2,1); plot(x,y)
    hold on
    scatter(Params.alpha,Params.statdist_alpha)
    hold off
    title('Fixed effect: i.i.d. normal')
    legend('Normal pdf','Discretization')
    y = normcdf(x,0,Params.sigma_alpha);
    subplot(1,2,2); plot(x,y)
    hold on
    scatter(Params.alpha,cumsum(Params.statdist_alpha))
    hold off
    title('Fixed effect: i.i.d. normal')
%     legend('Normal cdf','Discretization')
    saveas(fig,['./SavedOutput/EvaluateDiscretization/Model5_alpha_nSigma',num2str(nSigma_alpha),'.pdf'])

elseif useModel==6
    
    %% kappa_j
    % kappa_beta means that there are heterogeneous income profiles
    fig=figure(1);
    plot((1:1:Params.J)+Params.agejshifter,Params.kappa_j+Params.kappabeta(1,:))
    hold on
    for ii=2:n_kappabeta
        plot((1:1:Params.J)+Params.agejshifter,Params.kappa_j+Params.kappabeta(ii,:));
    end
    hold off
    title('Deterministic labor efficiency units (hetero income profiles)')
    xlabel('Age (in Years)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model6_kappaj.pdf')
    % NOTE: alpha and kappa_beta are jointly determined, so this may be
    % plotting all the same kappa_beta (with different values of alpha)
    
    %% epsilon
    gm = gmdistribution(Params.mu_epsilon,shiftdim(Params.sigma_epsilon,-2),Params.p_epsilon); % mu and sigma need to be row vectors
    fig=figure(2);
    x = linspace(1.5*min(epsilon_grid),1.5*max(epsilon_grid),201)';
    ypdf = pdf(gm,x);
    subplot(1,2,2); plot(x,ypdf)
    hold on
    scatter(epsilon_grid,pi_epsilon)
    hold off
    title('i.i.d. shock: gaussian-mixture, pdf')
    legend('Gaussian-mixture','Discretization')
    ycdf = cdf(gm,x);
    subplot(1,2,1); plot(x,ycdf)
    hold on
    scatter(epsilon_grid,cumsum(pi_epsilon))
    hold off
    title('i.i.d. shock: gaussian-mixture, cdf')
%     legend('Gaussian-mixture','Discretization')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model6_epsilon.pdf')
    
    %% z
    markovsimoptions.nsims=10^5;
    z_panel=MarkovChainSimulation_FHorz(z_grid_J,pi_z_J,jequaloneDistz,markovsimoptions);
    z_panel=z_panel(:,1:Params.Jr-1); % Drop the retirement age as is not relevant to the actual discretization of z
    % Calculate the correlation for each age: compare to Params.rho_z
    fig=figure(3);
    discretized_rho_z=nan(1,Params.Jr-1);
    for jj=2:Params.Jr-1
        discretized_rho_z(jj)=corr(z_panel(:,jj),z_panel(:,jj-1));
    end
    plot((1:1:(Params.Jr-1))+Params.agejshifter,Params.rho,(1:1:(Params.Jr-1))+Params.agejshifter,discretized_rho_z)
    title('Autocorrelation coefficient of z')
    legend('Original','Discretized')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model6_zcorr.pdf')
    % Calculate the innovations and plot them
    % z_innov=zeros(N,Params.Jr-1-1); % Extra minus one as first difference
    z_innov=z_panel(:,2:end)-z_panel(:,1:end-1);
    aa=floor((Params.Jr-1)/10); % Create 10 subplots per row
    bb=ceil((Params.Jr-1)/aa);  
    fig=figure(4);
    for djj=1:Params.Jr-2 % djj=1 will give transition from period 1 to 2
        x = linspace(1.5*min(z_innov(:,djj)),1.5*max(z_innov(:,djj)),201)';
        % Params.mu_eta is just the first, the second is determined to make E[eta]=0.
        mu_eta=[Params.mu_eta(djj); -(Params.mu_eta(djj).*Params.mixprobs_eta(1,djj))/Params.mixprobs_eta(2,djj)]; % Simple rearrangement of mu_i(:,jj).*mixprob_i(:,jj)=0, which is the requirement that mean of gaussian-mixture innovations=0
        gm = gmdistribution(mu_eta,shiftdim(Params.sigma_eta(:,djj),-2),Params.mixprobs_eta(:,djj));
        y = pdf(gm,x);
        subplot(bb,aa,djj); plot(x,y)
        hold on
        histogram(z_innov(:,djj),n_z,'Normalization','probability')
        hold off        
    end
    legend('Gaussian-mixture pdf','Discretization')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model6_zinnov_pdf.pdf')
    fig=figure(5);
    for djj=1:Params.Jr-2 % djj=1 will give transition from period 1 to 2
        x = linspace(1.5*min(z_innov(:,djj)),1.5*max(z_innov(:,djj)),201)';
        % Params.mu_eta is just the first, the second is determined to make E[eta]=0.
        mu_eta=[Params.mu_eta(djj); -(Params.mu_eta(djj).*Params.mixprobs_eta(1,djj))/Params.mixprobs_eta(2,djj)]; % Simple rearrangement of mu_i(:,jj).*mixprob_i(:,jj)=0, which is the requirement that mean of gaussian-mixture innovations=0
        gm = gmdistribution(mu_eta,shiftdim(Params.sigma_eta(:,djj),-2),Params.mixprobs_eta(:,djj));
        y = cdf(gm,x);
        subplot(bb,aa,djj); plot(x,y)
        hold on
        histogram(z_innov(:,djj),n_z,'Normalization','cdf')
        hold off        
    end
    legend('Gaussian-mixture cdf','Discretization')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model6_zinnov_cdf.pdf')
    
    % Does all of the z_grid get used? Plot the age conditional pdfs and cdfs
    z_statdist=zeros(n_z,Params.Jr);
    z_statdist(:,1)=jequaloneDistz;
    for jj=1:Params.Jr-1
        z_statdist(:,jj+1)=pi_z_J(:,:,jj)'*z_statdist(:,1);
    end
    aa=floor((Params.Jr)/10); % Create 10 subplots per row
    bb=ceil((Params.Jr)/aa);  
    fig=figure(6);
    for jj=1:Params.Jr-1
        subplot(bb,aa,jj); plot(z_grid_J(:,jj),z_statdist(:,jj))
        xlim([min(z_grid_J(:,jj)),max(z_grid_J(:,jj))])
    end
    legend('Discretized z pdf')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model6_z_pdf.pdf')
    fig=figure(7);
    for jj=1:Params.Jr-1
        subplot(bb,aa,jj); plot(z_grid_J(:,jj),cumsum(z_statdist(:,jj)))
        xlim([min(z_grid_J(:,jj)),max(z_grid_J(:,jj))])
    end
    legend('Discretized z cdf')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model6_z_cdf.pdf')

    % Calculate the standard deviation of the discretized z process and
    % compare to the theortical ones (which were used to create the grid for z)
    [mean_disc_z,variance_disc_z,corr_disc_z,statdist_z]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J,jequaloneDistz);
    fig=figure(8);
    subplot(2,1,1);  plot((1:1:(Params.Jr-1))+Params.agejshifter,sqrt(variance_disc_z(1:Params.Jr-1)), (1:1:(Params.Jr-1))+Params.agejshifter,otheroutputs_z.sigma_z)
    legend('Discretized','Theory')
    title('Standard deviation of z (from age to next age)')
    subplot(2,1,2); plot((1:1:Params.Jr-1)+Params.agejshifter,mean_disc_z(1:Params.Jr-1))
    ylim([-0.1,0.1])
    legend('Discretized')
    title('Mean of z (conditional on age)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model6_z_stddev.pdf')
    % Do another check of the autocorrelation coefficient
    fig=figure(9);
    plot((1:1:(Params.Jr-1))+Params.agejshifter,Params.rho,(1:1:(Params.Jr-1))+Params.agejshifter,discretized_rho_z,(2:1:(Params.Jr-1))+Params.agejshifter,corr_disc_z(2:Params.Jr-1))
    ylim([0.9,1])
    title('Autocorrelation coefficient of z')
    legend('Original','Discretized','Discretized double-check')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model6_zcorr2.pdf')
    
    % Plot the transition matrix at certain ages as 2D-pdfs
    fig=figure(10);
    subplot(2,2,1); surf(pi_z_J(:,:,1))
    title('Transition matrix pdf for j=1 (age 25)')
    subplot(2,2,2); surf(pi_z_J(:,:,11))
    title('Transition matrix pdf for j=11 (age 35)')
    subplot(2,2,3); surf(pi_z_J(:,:,21))
    title('Transition matrix pdf for j=21 (age 45)')
    subplot(2,2,4); surf(pi_z_J(:,:,31))
    title('Transition matrix pdf for j=31 (age 55)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model6_surf_pi_z_J.pdf')
    % Plot the transition matrix at certain ages as cdfs for each z
    fig=figure(11);
    subplot(2,2,1); plot(cumsum(pi_z_J(:,:,1),2)')
    title('Transition matrix cdf for j=1 (age 25)')
    subplot(2,2,2); plot(cumsum(pi_z_J(:,:,11),2)')
    title('Transition matrix cdf for j=11 (age 35)')
    subplot(2,2,3); plot(cumsum(pi_z_J(:,:,21),2)')
    title('Transition matrix cdf for j=21 (age 45)')
    subplot(2,2,4); plot(cumsum(pi_z_J(:,:,31),2)')
    title('Transition matrix cdf for j=31 (age 55)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model6_cdf_pi_z_J.pdf')
    
    % Plotting upsilon does not seem likely to be useful. 
    % Elsewhere I plot the life-cycle profile of mean of upsilon using model.
    % I also plot percentiles of age-contional earnings
        
    %% alpha and kappa_beta
    % They are joint normally distributed
    
    % Start by finding the marginal grid and marginal distribution
    marginal_alpha=unique(alpha_kappabeta_grid(1,:));
    marginaldist_alpha=zeros(1,length(marginal_alpha));
    for ii=1:length(marginal_alpha)
        for jj=1:length(alpha_kappabeta_grid(1,:))
            if marginal_alpha(ii)==alpha_kappabeta_grid(1,jj)
                marginaldist_alpha(ii)=marginaldist_alpha(ii)+statdist_alpha_kappabeta(jj);
            end
        end
    end
    marginal_kappabeta=unique(alpha_kappabeta_grid(2,:));
    marginaldist_kappabeta=zeros(1,length(marginal_kappabeta));
    for ii=1:length(marginal_kappabeta)
        for jj=1:length(alpha_kappabeta_grid(2,:))
            if marginal_kappabeta(ii)==alpha_kappabeta_grid(2,jj)
                marginaldist_kappabeta(ii)=marginaldist_kappabeta(ii)+statdist_alpha_kappabeta(jj);
            end
        end
    end

    % Following commented out lines did not work as the resulting matrix
    % was very sparse and so surf() and mesh() of it just don't produce
    % anything useful. I use scatter3() of the irregular grids instead
% %     % Because the actual grids (alpha_kappabeta_grid) are irregular we
% %     % cannot just use surf() or mesh() directly, we first need to create a
% %     % regular grid from these
% %     statdist_regular=nan(length(marginal_alpha),length(marginal_kappabeta));
% %     for ii=1:length(marginal_alpha)
% %         for kk=1:size(alpha_kappabeta_grid,2)
% %             if marginal_alpha(ii)==alpha_kappabeta_grid(1,kk)
% %                 for jj=1:length(marginal_kappabeta)
% %                     if marginal_kappabeta(jj)==alpha_kappabeta_grid(2,kk)
% %                         statdist_regular(ii,jj)=statdist_alpha_kappabeta(kk);
% %                     end
% %                 end
% %             end
% %         end
% %     end
    
    
    % Plot surface of joint pdfs
    fig=figure(13);
    subplot(1,2,1); scatter3(alpha_kappabeta_grid(:,1),alpha_kappabeta_grid(:,2), statdist_alpha_kappabeta,[],statdist_alpha_kappabeta,'filled','MarkerEdgeColor','k')
    xlabel('alpha_i')
    ylabel('kappa_beta_i')
    title('Fixed Effects: Discretized pdf')
    X1=linspace(-4*Params.sigma_alpha,4*Params.sigma_alpha,201)';
    X2=linspace(-4*Params.sigma_kappabeta,4*Params.sigma_kappabeta,201)';
    X=[kron(ones(201,1),X1),kron(X2,ones(201,1))];
    abMu=[0, 0];
    abSigma=[Params.sigma_alpha^2,  Params.covar_alpha_kappabeta; Params.covar_alpha_kappabeta, Params.sigma_kappabeta^2];
    y = mvnpdf(X,abMu,abSigma);
    y = reshape(y,[201,201]);
    subplot(1,2,2); surf(X1,X2,y)
    title('Fixed Effects: Joint-normal pdf')
    saveas(fig,['./SavedOutput/EvaluateDiscretization/Model6_alphabeta_joint_nSigma',num2str(nSigma_alpha),'.pdf'])

    % Plot the marginal pdfs
    fig=figure(12);
    subplot(2,2,1); plot(marginal_alpha,marginaldist_alpha)
    title('Marginal dist of discretized alpha')
    subplot(2,2,2); plot(marginal_kappabeta,marginaldist_kappabeta)
    title('Marginal dist of discretized kappa_beta')
    subplot(2,2,3); plot(X1,sum(y,2))
    title('Marginal dist of alpha')
    subplot(2,2,4); plot(X2,sum(y,1))
    title('Marginal dist of kappa_beta')
    saveas(fig,['./SavedOutput/EvaluateDiscretization/Model6_alphabeta_marginal_nSigma',num2str(nSigma_alpha),'.pdf'])
    
    % Plot the conditional pdfs
    % It is not obvious to me that plotting these is likely to be a good idea, 
    % I feel like they should look very different for different numbers of grid points 
    % and is not obvious what an 'accurate' approximation should look like.

    
end




%% That is the end of the main part of the code. Rest is just about useModel==21 and 22 (which are just useModel==2, but different method to discretize z)

%% Just do z for useModel 21 and 22, as these are otherwise just the same as useModel=2
if useModel==21
    %% z
    markovsimoptions.nsims=10^5;
    z_panel=MarkovChainSimulation_FHorz(z_grid_J,pi_z_J,jequaloneDistz,markovsimoptions);
    z_panel=z_panel(:,1:Params.Jr-1); % Drop the retirement age as is not relevant to the actual discretization of z
    % Calculate the correlation for each age: compare to Params.rho_z
    fig=figure(3);
    discretized_rho_z=nan(1,Params.Jr-1);
    for jj=2:Params.Jr-1
        discretized_rho_z(jj)=corr(z_panel(:,jj),z_panel(:,jj-1));
    end
    plot((1:1:(Params.Jr-1))+Params.agejshifter,Params.rho,(1:1:(Params.Jr-1))+Params.agejshifter,discretized_rho_z)
    title('Autocorrelation coefficient of z')
    legend('Original','Discretized')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model21_zcorr.pdf')
    % Calculate the innovations and plot them
    fig=figure(4);
    % z_innov=zeros(N,Params.Jr-1-1); % Extra minus one as first difference
    z_innov=z_panel(:,2:end)-z_panel(:,1:end-1);
    aa=floor((Params.Jr-1)/10); % Create 10 subplots per row
    bb=ceil((Params.Jr-1)/aa);  
    fig=figure(4);
    for djj=1:Params.Jr-2 % djj=1 will give transition from period 1 to 2
        x = linspace(1.5*min(z_innov(:,djj)),1.5*max(z_innov(:,djj)),201)';
        y = normpdf(x,0,Params.sigma_eta(djj));
        subplot(bb,aa,djj); plot(x,y)
        hold on
        histogram(z_innov(:,djj),n_z,'Normalization','probability')
        hold off        
    end
    legend('Normal pdf','Discretization')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model21_zinnov_pdf.pdf')
    fig=figure(5);
    for djj=1:Params.Jr-2 % djj=1 will give transition from period 1 to 2
        x = linspace(1.5*min(z_innov(:,djj)),1.5*max(z_innov(:,djj)),201)';
        y = normcdf(x,0,Params.sigma_eta(djj));
        subplot(bb,aa,djj); plot(x,y)
        hold on
        histogram(z_innov(:,djj),n_z,'Normalization','cdf')
        hold off        
    end
    legend('Normal cdf','Discretization')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model21_zinnov_cdf.pdf')
    
    % Does all of the z_grid get used? Plot the age conditional pdfs and cdfs
    z_statdist=zeros(n_z,Params.Jr);
    z_statdist(:,1)=jequaloneDistz;
    for jj=1:Params.Jr-1
        z_statdist(:,jj+1)=pi_z_J(:,:,jj)'*z_statdist(:,1);
    end
    aa=floor((Params.Jr)/10); % Create 10 subplots per row
    bb=ceil((Params.Jr)/aa);  
    fig=figure(6);
    for jj=1:Params.Jr-1
        subplot(bb,aa,jj); plot(z_grid_J(:,jj),z_statdist(:,jj))
        xlim([min(z_grid_J(:,jj)),max(z_grid_J(:,jj))])
    end
    legend('Discretized z pdf')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model21_z_pdf.pdf')
    fig=figure(7);
    for jj=1:Params.Jr-1
        subplot(bb,aa,jj); plot(z_grid_J(:,jj),cumsum(z_statdist(:,jj)))
        xlim([min(z_grid_J(:,jj)),max(z_grid_J(:,jj))])
    end
    legend('Discretized z cdf')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model21_z_cdf.pdf')
    
    % Calculate the standard deviation of the discretized z process and
    % compare to the theortical ones (which were used to create the grid for z)
    [mean_disc_z,variance_disc_z,corr_disc_z,statdist_z]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J,jequaloneDistz);
    fig=figure(8);
    subplot(2,1,1);  plot((1:1:(Params.Jr-1))+Params.agejshifter,sqrt(variance_disc_z(1:Params.Jr-1)), (1:1:(Params.Jr-1))+Params.agejshifter,otheroutputs_z.sigma_z)
    legend('Discretized','Theory')
    title('Standard deviation of z (from age to next age)')
    subplot(2,1,2); plot((1:1:Params.Jr-1)+Params.agejshifter,mean_disc_z(1:Params.Jr-1))
    ylim([-0.1,0.1])
    legend('Discretized')
    title('Mean of z (conditional on age)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model21_z_stddev.pdf')
    % Do another check of the autocorrelation coefficient
    fig=figure(9);
    plot((1:1:(Params.Jr-1))+Params.agejshifter,Params.rho,(1:1:(Params.Jr-1))+Params.agejshifter,discretized_rho_z,(2:1:(Params.Jr-1))+Params.agejshifter,corr_disc_z(2:Params.Jr-1))
    ylim([0.9,1])
    title('Autocorrelation coefficient of z')
    legend('Original','Discretized','Discretized double-check')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model21_zcorr2.pdf')
    
    % Plot the transition matrix at certain ages as 2D-pdfs
    fig=figure(10);
    subplot(2,2,1); surf(pi_z_J(:,:,1))
    title('Transition matrix pdf for j=1 (age 25)')
    subplot(2,2,2); surf(pi_z_J(:,:,11))
    title('Transition matrix pdf for j=11 (age 35)')
    subplot(2,2,3); surf(pi_z_J(:,:,21))
    title('Transition matrix pdf for j=21 (age 45)')
    subplot(2,2,4); surf(pi_z_J(:,:,31))
    title('Transition matrix pdf for j=31 (age 55)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model21_surf_pi_z_J.pdf')
    % Plot the transition matrix at certain ages as cdfs for each z
    fig=figure(11);
    subplot(2,2,1); plot(cumsum(pi_z_J(:,:,1),2)')
    title('Transition matrix cdf for j=1 (age 25)')
    subplot(2,2,2); plot(cumsum(pi_z_J(:,:,11),2)')
    title('Transition matrix cdf for j=11 (age 35)')
    subplot(2,2,3); plot(cumsum(pi_z_J(:,:,21),2)')
    title('Transition matrix cdf for j=21 (age 45)')
    subplot(2,2,4); plot(cumsum(pi_z_J(:,:,31),2)')
    title('Transition matrix cdf for j=31 (age 55)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model21_cdf_pi_z_J.pdf')
    
elseif useModel==22
    %% z
    markovsimoptions.nsims=10^5;
    z_panel=MarkovChainSimulation_FHorz(z_grid_J,pi_z_J,jequaloneDistz,markovsimoptions);
    z_panel=z_panel(:,1:Params.Jr-1); % Drop the retirement age as is not relevant to the actual discretization of z
    % Calculate the correlation for each age: compare to Params.rho_z
    fig=figure(3);
    discretized_rho_z=nan(1,Params.Jr-1);
    for jj=2:Params.Jr-1
        discretized_rho_z(jj)=corr(z_panel(:,jj),z_panel(:,jj-1));
    end
    plot((1:1:(Params.Jr-1))+Params.agejshifter,Params.rho,(1:1:(Params.Jr-1))+Params.agejshifter,discretized_rho_z)
    title('Autocorrelation coefficient of z')
    legend('Original','Discretized')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model22_zcorr.pdf')
    % Calculate the innovations and plot them
    % z_innov=zeros(N,Params.Jr-1-1); % Extra minus one as first difference
    z_innov=z_panel(:,2:end)-z_panel(:,1:end-1);
    aa=floor((Params.Jr-1)/10); % Create 10 subplots per row
    bb=ceil((Params.Jr-1)/aa);  
    fig=figure(4);
    for djj=1:Params.Jr-2 % djj=1 will give transition from period 1 to 2
        x = linspace(1.5*min(z_innov(:,djj)),1.5*max(z_innov(:,djj)),201)';
        y = normpdf(x,0,Params.sigma_eta(djj));
        subplot(bb,aa,djj); plot(x,y)
        hold on
        histogram(z_innov(:,djj),n_z,'Normalization','probability')
        hold off        
    end
    legend('Normal pdf','Discretization')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model22_zinnov_pdf.pdf')
    fig=figure(5);
    for djj=1:Params.Jr-2 % djj=1 will give transition from period 1 to 2
        x = linspace(1.5*min(z_innov(:,djj)),1.5*max(z_innov(:,djj)),201)';
        y = normcdf(x,0,Params.sigma_eta(djj));
        subplot(bb,aa,djj); plot(x,y)
        hold on
        histogram(z_innov(:,djj),n_z,'Normalization','cdf')
        hold off        
    end
    legend('Normal cdf','Discretization')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model22_zinnov_cdf.pdf')
    
    % Does all of the z_grid get used? Plot the age conditional pdfs and cdfs
    z_statdist=zeros(n_z,Params.Jr);
    z_statdist(:,1)=jequaloneDistz;
    for jj=1:Params.Jr-1
        z_statdist(:,jj+1)=pi_z_J(:,:,jj)'*z_statdist(:,1);
    end
    aa=floor((Params.Jr)/10); % Create 10 subplots per row
    bb=ceil((Params.Jr)/aa);  
    fig=figure(6);
    for jj=1:Params.Jr-1
        subplot(bb,aa,jj); plot(z_grid_J(:,jj),z_statdist(:,jj))
        xlim([min(z_grid_J(:,jj)),max(z_grid_J(:,jj))])
    end
    legend('Discretized z pdf')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model22_z_pdf.pdf')
    fig=figure(7);
    for jj=1:Params.Jr-1
        subplot(bb,aa,jj); plot(z_grid_J(:,jj),cumsum(z_statdist(:,jj)))
        xlim([min(z_grid_J(:,jj)),max(z_grid_J(:,jj))])
    end
    legend('Discretized z cdf')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model22_z_cdf.pdf')
    
    % Calculate the standard deviation of the discretized z process and
    % compare to the theortical ones (which were used to create the grid for z)
    [mean_disc_z,variance_disc_z,corr_disc_z,statdist_z]=MarkovChainMoments_FHorz(z_grid_J,pi_z_J,jequaloneDistz);
    fig=figure(8);
    subplot(2,1,1);  plot((1:1:(Params.Jr-1))+Params.agejshifter,sqrt(variance_disc_z(1:Params.Jr-1)), (1:1:(Params.Jr-1))+Params.agejshifter,otheroutputs_z.sigma_z)
    legend('Discretized','Theory')
    title('Standard deviation of z (from age to next age)')
    subplot(2,1,2); plot((1:1:Params.Jr-1)+Params.agejshifter,mean_disc_z(1:Params.Jr-1))
    ylim([-0.1,0.1])
    legend('Discretized')
    title('Mean of z (conditional on age)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model22_z_stddev.pdf')
    % Do another check of the autocorrelation coefficient
    fig=figure(9);
    plot((1:1:(Params.Jr-1))+Params.agejshifter,Params.rho,(1:1:(Params.Jr-1))+Params.agejshifter,discretized_rho_z,(2:1:(Params.Jr-1))+Params.agejshifter,corr_disc_z(2:Params.Jr-1))
    ylim([0.9,1])
    title('Autocorrelation coefficient of z')
    legend('Original','Discretized','Discretized double-check')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model22_zcorr2.pdf')
    
    % Plot the transition matrix at certain ages as 2D-pdfs
    fig=figure(10);
    subplot(2,2,1); surf(pi_z_J(:,:,1))
    title('Transition matrix pdf for j=1 (age 25)')
    subplot(2,2,2); surf(pi_z_J(:,:,11))
    title('Transition matrix pdf for j=11 (age 35)')
    subplot(2,2,3); surf(pi_z_J(:,:,21))
    title('Transition matrix pdf for j=21 (age 45)')
    subplot(2,2,4); surf(pi_z_J(:,:,31))
    title('Transition matrix pdf for j=31 (age 55)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model22_surf_pi_z_J.pdf')
    % Plot the transition matrix at certain ages as cdfs for each z
    fig=figure(11);
    subplot(2,2,1); plot(cumsum(pi_z_J(:,:,1),2)')
    title('Transition matrix cdf for j=1 (age 25)')
    subplot(2,2,2); plot(cumsum(pi_z_J(:,:,11),2)')
    title('Transition matrix cdf for j=11 (age 35)')
    subplot(2,2,3); plot(cumsum(pi_z_J(:,:,21),2)')
    title('Transition matrix cdf for j=21 (age 45)')
    subplot(2,2,4); plot(cumsum(pi_z_J(:,:,31),2)')
    title('Transition matrix cdf for j=31 (age 55)')
    saveas(fig,'./SavedOutput/EvaluateDiscretization/Model22_cdf_pi_z_J.pdf')
end



%% How much of the z grid is actually being used if nSigmaz=4?? 
% Following is to take a look.

% First, how many points have more than 10^(-9) of mass?
FID = fopen(['./SavedOutput/EvaluateDiscretization/Model',num2str(useModel),'_nz',num2str(n_z),'_nSigmaz',num2str(nSigmaz),'_statdist_z.tex'], 'w');
fprintf(FID,'How much of z grid is used? (1/3) \\\\ \n');
for jj=1:Params.Jr-1 % Only interested in discretization for working ages
    fprintf(FID,'  %i of %i points have more than 10\\^{}(-9) mass at age j=%i \\\\ \n',sum((statdist_z(:,jj)>10^(-9))),n_z,jj);
end
fprintf(FID,'How much of z grid is used? (2/3) \\\\ \n');
ipoints=10; % Look at bottom and top ipoints
for jj=1:Params.Jr-1 % Only interested in discretization for working ages
    fprintf(FID,'  Bottom %i points sum to mass %8.6f, top %i points sum to mass %8.6f, at age j=%i \\\\ \n',ipoints,sum(statdist_z(1:ipoints,jj)),ipoints,sum(statdist_z(end-ipoints+1:end,jj)),jj);
end
fprintf(FID,'How much of z grid is used? (3/3) \\\\ \n');
ipoints=5; % Look at bottom and top ipoints
for jj=1:Params.Jr-1 % Only interested in discretization for working ages
    fprintf(FID,'  Bottom %i points sum to mass %8.6f, top %i points sum to mass %8.6f, at age j=%i  \\\\ \n',ipoints,sum(statdist_z(1:ipoints,jj)),ipoints,sum(statdist_z(end-ipoints+1:end,jj)),jj);
end
fclose(FID);
% Based on the contents of this I decided to do a much shorter version
FID = fopen(['./SavedOutput/EvaluateDiscretization/Model',num2str(useModel),'_nz',num2str(n_z),'_nSigmaz',num2str(nSigmaz),'_statdist_z_shortversion.tex'], 'w');
fprintf(FID,'How much of z grid is used? (1/4) \\\\ \n');
for jj=1:Params.Jr-1 % Only interested in discretization for working ages
    fprintf(FID,'  %i of %i points have more than 10\\^{}(-9) mass at age j=%i \\\\ \n',sum((statdist_z(:,jj)>10^(-9))),n_z,jj);
end
fprintf(FID,'How much of z grid is used? (2/4) \\\\ \n');
ipoints=10; % Look at bottom and top ipoints
for jj=1:10:Params.Jr-1 % Only interested in discretization for working ages
    fprintf(FID,'  Bottom %i points sum to mass %8.6f, top %i points sum to mass %8.6f, at age j=%i \\\\ \n',ipoints,sum(statdist_z(1:ipoints,jj)),ipoints,sum(statdist_z(end-ipoints+1:end,jj)),jj);
end
fprintf(FID,'How much of z grid is used? (3/4) \\\\ \n');
ipoints=5; % Look at bottom and top ipoints
for jj=1:10:Params.Jr-1 % Only interested in discretization for working ages
    fprintf(FID,'  Bottom %i points sum to mass %8.6f, top %i points sum to mass %8.6f, at age j=%i  \\\\ \n',ipoints,sum(statdist_z(1:ipoints,jj)),ipoints,sum(statdist_z(end-ipoints+1:end,jj)),jj);
end
fprintf(FID,'Range of z grid? (4/4) \\\\ \n');
for jj=1:10:Params.Jr-1 % Only interested in discretization for working ages
    fprintf(FID,'  Min grid point %8.6f, max grid point %8.6f, at age j=%i  \\\\ \n',z_grid_J(1,jj),z_grid_J(end,jj),jj);
end
fclose(FID);

% Store these so that I can put together a big table about them that looks across many different useModel/n_z/nSigmaz
TableAboutzgridtemp=nan(5,7);  % fourages-sevennumbersofinterest (for models 2,21,22,6)
TableAboutzgridtemp(:,1)=[sum((statdist_z(:,1)>10^(-9))); sum((statdist_z(:,11)>10^(-9))); sum((statdist_z(:,21)>10^(-9))); sum((statdist_z(:,31)>10^(-9))); sum((statdist_z(:,41)>10^(-9)))];
ipoints=10;
TableAboutzgridtemp(:,2)=[sum(statdist_z(1:ipoints,1));sum(statdist_z(1:ipoints,11));sum(statdist_z(1:ipoints,21));sum(statdist_z(1:ipoints,31));sum(statdist_z(1:ipoints,41))];
TableAboutzgridtemp(:,3)=[sum(statdist_z(end-ipoints+1:end,1)); sum(statdist_z(end-ipoints+1:end,11)); sum(statdist_z(end-ipoints+1:end,21)); sum(statdist_z(end-ipoints+1:end,31)); sum(statdist_z(end-ipoints+1:end,41))];
ipoints=5;
TableAboutzgridtemp(:,4)=[sum(statdist_z(1:ipoints,1));sum(statdist_z(1:ipoints,11));sum(statdist_z(1:ipoints,21));sum(statdist_z(1:ipoints,31));sum(statdist_z(1:ipoints,41))];
TableAboutzgridtemp(:,5)=[sum(statdist_z(end-ipoints+1:end,1)); sum(statdist_z(end-ipoints+1:end,11)); sum(statdist_z(end-ipoints+1:end,21)); sum(statdist_z(end-ipoints+1:end,31)); sum(statdist_z(end-ipoints+1:end,41))];
TableAboutzgridtemp(:,6)=[z_grid_J(1,1); z_grid_J(1,11); z_grid_J(1,21); z_grid_J(1,31); z_grid_J(1,41)];
TableAboutzgridtemp(:,7)=[z_grid_J(end,1); z_grid_J(end,11); z_grid_J(end,21); z_grid_J(end,31); z_grid_J(end,41)];
fig=figure(13);
subplot(2,2,1); plot(cumsum(statdist_z(:,1)))
xlim([1,n_z])
title('Stationary distn of z cdf for j=1 (age 25)')
subplot(2,2,2); plot(cumsum(statdist_z(:,11)))
xlim([1,n_z])
title('Stationary distn of z cdf for j=11 (age 35)')
subplot(2,2,3); plot(cumsum(statdist_z(:,21)))
xlim([1,n_z])
title('Stationary distn of z cdf for j=21 (age 45)')
subplot(2,2,4); plot(cumsum(statdist_z(:,31)))
xlim([1,n_z])
title('Stationary distn of z cdf for j=31 (age 55)')
saveas(fig,['./SavedOutput/EvaluateDiscretization/Model',num2str(useModel),'_nz',num2str(n_z),'_nSigmaz',num2str(nSigmaz),'_cdf_statdist_z.pdf'])

% These transition matrix graphs were done earlier, but now just using different filenames to save them
% Plot the transition matrix at certain ages as 2D-pdfs
fig=figure(14);
subplot(2,2,1); surf(pi_z_J(:,:,1))
title('Transition matrix pdf for j=1 (age 25)')
subplot(2,2,2); surf(pi_z_J(:,:,11))
title('Transition matrix pdf for j=11 (age 35)')
subplot(2,2,3); surf(pi_z_J(:,:,21))
title('Transition matrix pdf for j=21 (age 45)')
subplot(2,2,4); surf(pi_z_J(:,:,31))
title('Transition matrix pdf for j=31 (age 55)')
saveas(fig,['./SavedOutput/EvaluateDiscretization/Model',num2str(useModel),'_nz',num2str(n_z),'_nSigmaz',num2str(nSigmaz),'_surf_pi_z_J.pdf'])
% Plot the transition matrix at certain ages as cdfs for each z
fig=figure(15);
subplot(2,2,1); plot(cumsum(pi_z_J(:,:,1),2)')
title('Transition matrix cdf for j=1 (age 25)')
subplot(2,2,2); plot(cumsum(pi_z_J(:,:,11),2)')
title('Transition matrix cdf for j=11 (age 35)')
subplot(2,2,3); plot(cumsum(pi_z_J(:,:,21),2)')
title('Transition matrix cdf for j=21 (age 45)')
subplot(2,2,4); plot(cumsum(pi_z_J(:,:,31),2)')
title('Transition matrix cdf for j=31 (age 55)')
saveas(fig,['./SavedOutput/EvaluateDiscretization/Model',num2str(useModel),'_nz',num2str(n_z),'_nSigmaz',num2str(nSigmaz),'_cdf_pi_z_J.pdf'])

if useModel==2 || useModel==21 || useModel==22
    % Calculate the innovations and plot them
    % z_innov=zeros(N,Params.Jr-1-1); % Extra minus one as first difference
    z_innov=z_panel(:,2:end)-z_panel(:,1:end-1);
    aa=floor((Params.Jr-1)/10); % Create 10 subplots per row
    bb=ceil((Params.Jr-1)/aa);
    fig=figure(16);
    djjvec=[1,11,21,31];
    for djjc=1:4 % djj=1 will give transition from period 1 to 2
        djj=djjvec(djjc);
        x = linspace(1.5*min(z_innov(:,djj)),1.5*max(z_innov(:,djj)),201)';
        y = normpdf(x,0,Params.sigma_eta(djj));
        
        subplot(2,2,djjc); plot(x,y)
        hold on
        histogram(z_innov(:,djj),n_z,'Normalization','probability')
        title(['pdf of innovations from period ',num2str(djj),' to period ',num2str(djj+1)])
        hold off
    end
    legend('Normal pdf','Discretization')
    saveas(fig,['./SavedOutput/EvaluateDiscretization/Model',num2str(useModel),'_nz',num2str(n_z),'_nSigmaz',num2str(nSigmaz),'_zinnov_pdf.pdf'])
    fig=figure(17);
    djjvec=[1,11,21,31];
    for djjc=1:4 % djj=1 will give transition from period 1 to 2
        djj=djjvec(djjc);
        x = linspace(1.5*min(z_innov(:,djj)),1.5*max(z_innov(:,djj)),201)';
        y = normcdf(x,0,Params.sigma_eta(djj));
        subplot(2,2,djjc); plot(x,y)
        hold on
        histogram(z_innov(:,djj),n_z,'Normalization','cdf')
        title(['cdf of innovations from period ',num2str(djj),' to period ',num2str(djj+1)])
        hold off
    end
    legend('Normal cdf','Discretization')
    saveas(fig,['./SavedOutput/EvaluateDiscretization/Model',num2str(useModel),'_nz',num2str(n_z),'_nSigmaz',num2str(nSigmaz),'_zinnov_cdf.pdf'])
elseif useModel==6
    % Calculate the innovations and plot them
    % z_innov=zeros(N,Params.Jr-1-1); % Extra minus one as first difference
    z_innov=z_panel(:,2:end)-z_panel(:,1:end-1);
    fig=figure(16);
    djjvec=[1,11,21,31];
    for djjc=1:4 % djj=1 will give transition from period 1 to 2
        djj=djjvec(djjc);
        x = linspace(1.5*min(z_innov(:,djj)),1.5*max(z_innov(:,djj)),201)';
        % Params.mu_eta is just the first, the second is determined to make E[eta]=0.
        mu_eta=[Params.mu_eta(djj); -(Params.mu_eta(djj).*Params.mixprobs_eta(1,djj))/Params.mixprobs_eta(2,djj)]; % Simple rearrangement of mu_i(:,jj).*mixprob_i(:,jj)=0, which is the requirement that mean of gaussian-mixture innovations=0
        gm = gmdistribution(mu_eta,shiftdim(Params.sigma_eta(:,djj),-2),Params.mixprobs_eta(:,djj));
        y = pdf(gm,x);
        subplot(2,2,djjc); plot(x,y)
        hold on
        histogram(z_innov(:,djj),n_z,'Normalization','probability')
        title(['pdf of innovations from period ',num2str(djj),' to period ',num2str(djj+1)])
        hold off
    end
    legend('Gaussian-mixture pdf','Discretization')
    saveas(fig,['./SavedOutput/EvaluateDiscretization/Model',num2str(useModel),'_nz',num2str(n_z),'_nSigmaz',num2str(nSigmaz),'_zinnov_pdf.pdf'])
    fig=figure(17);
    djjvec=[1,11,21,31];
    for djjc=1:4 % djj=1 will give transition from period 1 to 2
        djj=djjvec(djjc);
        x = linspace(1.5*min(z_innov(:,djj)),1.5*max(z_innov(:,djj)),201)';
        % Params.mu_eta is just the first, the second is determined to make E[eta]=0.
        mu_eta=[Params.mu_eta(djj); -(Params.mu_eta(djj).*Params.mixprobs_eta(1,djj))/Params.mixprobs_eta(2,djj)]; % Simple rearrangement of mu_i(:,jj).*mixprob_i(:,jj)=0, which is the requirement that mean of gaussian-mixture innovations=0
        gm = gmdistribution(mu_eta,shiftdim(Params.sigma_eta(:,djj),-2),Params.mixprobs_eta(:,djj));
        y = cdf(gm,x);
        subplot(2,2,djjc); plot(x,y)
        hold on
        histogram(z_innov(:,djj),n_z,'Normalization','cdf')
        title(['cdf of innovations from period ',num2str(djj),' to period ',num2str(djj+1)])
        hold off        
    end
    legend('Gaussian-mixture cdf','Discretization')
    saveas(fig,['./SavedOutput/EvaluateDiscretization/Model',num2str(useModel),'_nz',num2str(n_z),'_nSigmaz',num2str(nSigmaz),'_zinnov_cdf.pdf'])
end



%% Evaluate Epi_upsilon, the age-conditional transition probabilies for pi_upsilon that are the expectation over z (taken using stationary dist)
% If it is a model that uses upsilon
if useModel==2 || useModel==5 || useModel==6 || useModel==21 || useModel==22
   Epi_upsilon_J=zeros(2,2,Params.Jr-1);
   for jj=1:Params.Jr-1
       temp=sum(pi_zupsilon_J(1:n_z,1:n_z,jj),2); % upsilon=0 to upsilon'=0. Sum over z. 
       Epi_upsilon_J(1,1,jj)=sum(temp.*statdist_z(:,jj));
       temp=sum(pi_zupsilon_J(1:n_z,(n_z+1):(2*n_z),jj),2); % upsilon=0 to upsilon'=1. Sum over z. 
       Epi_upsilon_J(1,2,jj)=sum(temp.*statdist_z(:,jj));
       temp=sum(pi_zupsilon_J((n_z+1):(2*n_z),1:n_z,jj),2); % upsilon=1 to upsilon'=0. Sum over z. 
       Epi_upsilon_J(2,1,jj)=sum(temp.*statdist_z(:,jj));
       temp=sum(pi_zupsilon_J((n_z+1):(2*n_z),(n_z+1):(2*n_z),jj),2); % upsilon=1 to upsilon'=1. Sum over z. 
       Epi_upsilon_J(2,2,jj)=sum(temp.*statdist_z(:,jj));
   end
else
    Epi_upsilon_J=nan; % This is just so that I can save Epi_upsilon without having to worry if it exists
end


