clc
clear
cd '-----'

raw_data = readtable('DATI_FISCAL_CK.txt', 'Delimiter','	'); % load the dataset
data=table2array(raw_data(:,{'Var4','Var5','Var1'})); % select the variables: Tax, Gov, GDP

A_res=[1,0,-2.08;0,1,0;nan,nan,1]; % A matrix restrictions
B_res=[nan,nan,0;0,nan,0;0,0,nan]; % B matrix restrictions
det_comp='ct';                     % deterministic components including constant and trend)
lags=4;                            % number of lags of the VAR
HorizonIRF=20;                     % time horizon of the IRF
repetitions=1000;                  % bootstrap repetitions

% estimate structural SVAR using the SVAR_AB function (MLE)
[A, B, SE_A, SE_B, Sigma, IRF, PHI_Boot]=SVAR_AB(data, lags, det_comp, A_res, B_res, HorizonIRF, repetitions); 

%% bootstrap confidence intervals (90%)
PHI_Inf_Boot = prctile(PHI_Boot,5,4);  % lower 5th percentile across bootstrap repetitions
PHI_Sup_Boot = prctile(PHI_Boot,95,4); % upper 95th percentile across bootstrap repetitions

for h = 0 : HorizonIRF
    for i=1:3
        for j=1:3
            IRF_inf(h+1,3*(i-1)+j)=PHI_Inf_Boot(i,j,h+1); % fill lower confidence bound for (i,j) IRF
            IRF_sup(h+1,3*(i-1)+j)=PHI_Sup_Boot(i,j,h+1); % fill upper confidence bound for (i,j) IRF
        end
    end
end

%% IRFs plots
for i = 1:3 
    if i==1, var='Tax';
    elseif i==2, var='G';
    elseif i==3, var='GDP';
    end

% (1) response to the identified tax shock
figure;
col = 1+(i-1)*3;         % column index for tax shock IRF
x = 0:HorizonIRF;        % IRF horizon vector

fill([x fliplr(x)], [IRF_sup(:,col)' fliplr(IRF_inf(:,col)')], ...
     [0.7 0.9 0.9], 'EdgeColor','none','FaceAlpha',0.35);            % shaded 90% CI bond
hold on; 

median_line = median(PHI_Boot(i,1,1:HorizonIRF+1,:),4); 
plot(x, squeeze(median_line), 'Color',[0 0.5 0.5], 'LineWidth',0.8); % plot bootstrap median line
plot(x, IRF(:,col), 'k', 'LineWidth',0.8);                           % plot main IRF
plot(xlim,[0 0],'k--','LineWidth',0.5);                              % plot zero line

ax = gca; 
ax.YAxis.Exponent = 0;   % disable scientific notation
ax.XLabel = [];          % remove x-axis label       
grid on;                 % add grid for clearer reading

legend({'Bootstrap 90% CI','Bootstrap median','IRF'}, 'Orientation','horizontal', 'Location','southoutside', 'Box','off', 'FontSize',10);
title(['Response of ' var ' to the identified tax shock'], 'Interpreter','none','FontSize',10);
print(['graphs/IRF_' var '_to_Tax'],'-dpng','-r300')


% (2) response to the identified G shock
figure;
col = 2+(i-1)*3;         % column index for G shock IRF
x = 0:HorizonIRF;        % IRF horizon vector

fill([x fliplr(x)], [IRF_sup(:,col)' fliplr(IRF_inf(:,col)')], ...
     [0.7 0.9 0.9], 'EdgeColor','none','FaceAlpha',0.35);            % shaded 90% CI band
hold on;

median_line = median(PHI_Boot(i,2,1:HorizonIRF+1,:),4); 
plot(x, squeeze(median_line), 'Color',[0 0.5 0.5], 'LineWidth',0.8); % bootstrap median line
plot(x, IRF(:,col), 'k', 'LineWidth',0.8);                           % main IRF
plot(xlim,[0 0],'k--','LineWidth',0.5);                              % zero line

ax = gca;
ax.YAxis.Exponent = 0;   % disable scientific notation
ax.XLabel = [];          % remove x-axis label
grid on;                 % grid for readability

legend({'Bootstrap 90% CI','Bootstrap median','IRF'}, 'Orientation','horizontal', 'Location','southoutside', 'Box','off', 'FontSize',10);
title(['Response of ' var ' to the identified Gov shock'], 'Interpreter','none','FontSize',10);
print(['graphs/IRF_' var '_to_G'],'-dpng','-r300')


% (3) response to the identified GDP shock
figure;
col = 3+(i-1)*3;         % column index for GDP shock IRF
x = 0:HorizonIRF;        % IRF horizon vector

fill([x fliplr(x)], [IRF_sup(:,col)' fliplr(IRF_inf(:,col)')], ...
     [0.7 0.9 0.9], 'EdgeColor','none','FaceAlpha',0.35);            % shaded 90% CI band
hold on;

median_line = median(PHI_Boot(i,3,1:HorizonIRF+1,:),4); 
plot(x, squeeze(median_line), 'Color',[0 0.5 0.5], 'LineWidth',0.8); % bootstrap median line
plot(x, IRF(:,col), 'k', 'LineWidth',0.8);                           % main IRF
plot(xlim,[0 0],'k--','LineWidth',0.5);                              % zero line

ax = gca;
ax.YAxis.Exponent = 0;   % disable scientific notation
ax.XLabel = [];          % remove x-axis label
grid on;                 

legend({'Bootstrap 90% CI','Bootstrap median','IRF'}, 'Orientation','horizontal', 'Location','southoutside', 'Box','off', 'FontSize',10);
title(['Response of ' var ' to the identified GDP shock'], 'Interpreter','none','FontSize',10);
print(['graphs/IRF_' var '_to_GDP'],'-dpng','-r300')

end

%% multipliers
data_non_logged = exp(data);                               % convert logged data back to levels
ratio = [data_non_logged(:,3)./data_non_logged(:,1), ...
         data_non_logged(:,3)./data_non_logged(:,2)];      % compute GDP/Tax and GDP/G for each quarter
ratio_mean = mean(ratio);                                  % average GDP-to-tax and GDP-to-G ratios

% Tax multiplier
multiplier_tax = (IRF(:,7) / IRF(1,1)) * ratio_mean(1);    % GDP IRF to tax shock scaled by GDP/Tax ratio

% Gov multiplier
multiplier_G = (IRF(:,8) / IRF(1,5)) * ratio_mean(2);      % GDP IRF to G shock scaled by GDP/G ratio

%% multipliers bootstrap

% compute multipliers at each horizon for every bootstrap repetition
multiplier_tax_boot=[]; 
multiplier_G_boot=[];
for boot=1:repetitions
    for h=0:HorizonIRF
        multiplier_tax_boot(h+1,boot)=PHI_Boot(3,1,h+1,boot)/PHI_Boot(1,1,1,boot)*ratio_mean(1); % tax multiplier in bootstrap draw
        multiplier_G_boot(h+1,boot)=PHI_Boot(3,2,h+1,boot)/PHI_Boot(2,2,1,boot)*ratio_mean(2);   % G multiplier in bootstrap draw
    end
end 

% computing the 5th and 95th bootstrap percentiles
multiplier_tax_inf=prctile(multiplier_tax_boot,5,2);  % lower CI for tax multiplier
multiplier_tax_sup=prctile(multiplier_tax_boot,95,2); % upper CI for tax multiplier
multiplier_G_inf=prctile(multiplier_G_boot,5,2);      % lower CI for Gov multiplier
multiplier_G_sup=prctile(multiplier_G_boot,95,2);     % upper CI for Gov multiplier

%% multipliers plots

% (1) plot the tax multiplier
figure;
x = 0:HorizonIRF;

fill([x fliplr(x)], [multiplier_tax_sup' fliplr(multiplier_tax_inf')], ...
     [0.7 0.9 0.9], 'EdgeColor','none','FaceAlpha',0.35);                     % shaded 90% CI
hold on;

median_tax = median(multiplier_tax_boot,2);                                   % bootstrap median
plot(x, median_tax, 'Color',[0 0.5 0.5], 'LineWidth',0.8);                    % median line
plot(x, multiplier_tax, 'k', 'LineWidth',0.8);                                % main multiplier line
plot(xlim,[0 0],'k--','LineWidth',0.5);                                       % zero line

ax = gca;
ax.YAxis.Exponent = 0;    % disable scientific notation
ax.XLabel = [];           % remove x-axis label
grid on;

legend({'Bootstrap 90% CI','Bootstrap median','Multiplier'}, 'Orientation','horizontal','Location','southoutside','Box','off','FontSize',10);
title('Tax multiplier ', 'Interpreter','none','FontSize',10);
print('graphs/tax_multiplier','-dpng','-r300')


% (2) plot the Gov multiplier
figure;
x = 0:HorizonIRF;

fill([x fliplr(x)], [multiplier_G_sup' fliplr(multiplier_G_inf')], ...
     [0.7 0.9 0.9], 'EdgeColor','none','FaceAlpha',0.35);                     % shaded 90% CI
hold on;

median_G = median(multiplier_G_boot,2);                                       % bootstrap median
plot(x, median_G, 'Color',[0 0.5 0.5], 'LineWidth',0.8);                      % median line
plot(x, multiplier_G, 'k', 'LineWidth',0.8);                                  % main multiplier line
plot(xlim,[0 0],'k--','LineWidth',0.5);                                       % zero line

ax = gca;
ax.YAxis.Exponent = 0;    % disable scientific notation
ax.XLabel = [];           % remove x-axis label
grid on;

legend({'Bootstrap 90% CI','Bootstrap median','Multiplier'}, 'Orientation','horizontal','Location','southoutside','Box','off','FontSize',10);
title('Fiscal spending multiplier', 'Interpreter','none','FontSize',10);
print('graphs/G_multiplier','-dpng','-r300')

