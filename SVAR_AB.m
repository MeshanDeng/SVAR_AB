function [A, B, SE_A,SE_B,Sigma, IRF, PHI_Boot ] = SVAR_AB(data, lags, det_comp, A_res,B_res, HorizonIRF, repetitions)

% OUTPUTS
% A = estimated A matrix
% B = estimated B matrix
% SE_A = estimated standard errors of A
% SE_B = estimated standard errors of B
% Sigma = covariance matrix of reduced-form residuals
% IRF = impulse response functions
% PHI_Boot = bootstrap IRFs

% INPUTS
% data = data matrix
% lags = lags of the VAR
% det_comp = deterministic components
% A_res = A matrix restrictions
% B_res = B matrix restrictions
% HorizonIRF = IRF horizon
% repetitions = bootstrap repetitions


%% set up VAR model
    if ~(strcmp(det_comp, 'ct') || strcmp(det_comp, 't') || strcmp(det_comp, 'c') || strcmp(det_comp, 'none'))
        disp('Unknown deterministic component')  % check deterministic component
    end 

    T = size(data,1)-lags;  % effective sample size
    M = size(data,2);       % number of variables

    VAR_PI = {}; 
    for i = 1 : lags
        VAR_PI{i} = nan(M); % initialize coefficient matrices                   
    end 

    if strcmp(det_comp, 'c')        % for constant
        VAR_c = nan(M,1);
        VAR_t = zeros(M,1);
    elseif strcmp(det_comp, 'ct')   % for constant and trend
        VAR_c = nan(M,1);
        VAR_t = nan(M,1);
    elseif strcmp(det_comp, 't')    % for trend
        VAR_t = nan(M,1);
        VAR_c = zeros(M,1);
    elseif strcmp(det_comp, 'none') % for nothing
        VAR_c = zeros(M,1);
        VAR_t = zeros(M,1);
    end 
    
    VAR = varm('Constant',VAR_c,'AR',VAR_PI,'Trend',VAR_t); % build VAR model container

%% estimate reduced-form VAR
    [EstVAR,EstSE,logLikVAR,Residuals] = estimate(VAR,data,'Display',"full"); %  ML estimation
    Sigma = EstVAR.Covariance; % estimated covariance matrix

%% estimate structural SVAR parameters (MLE)
    StructuralParam = sum(sum(isnan(A_res)))+sum(sum(isnan(B_res))); % number of free structural parameters
    InitialValues=randn(StructuralParam,1)/10; % initial values for the likelihood maximization

    params_A=find(isnan(A_res));  % indices of free entries in A
    params_B=find(isnan(B_res));  % indices of free entries in B

    options = optimset('MaxFunEvals',200000,'TolFun',1e-500,'MaxIter',200000,'TolX',1e-100);
    fun = @(theta) Likelihood_SVAR_AB(theta, A_res, B_res, Sigma, T, M, params_A, params_B); 
    
    [StructuralParam_Estimation_MATRIX, n_logLik,exitflag,output,grad, Hessian_MATRIX] = fminunc(fun, InitialValues, options); % maximize log-likelihood (via minimizing its negative)
    SE_Hessian_MATRIX = diag(inv(Hessian_MATRIX)).^0.5; % standard errors from inverse Hessian

%% construct A and B, and test overidentification
    A = A_res; 
    SE_A = nan(size(Sigma,1)); 
    HSelection_A = zeros(M*M,size(params_A,1)); % selection matrix mapping theta elements to A

    for c_par = 1 : size(params_A,1) 
        A(params_A(c_par,1)) = StructuralParam_Estimation_MATRIX(c_par);
        SE_A(params_A(c_par,1)) = SE_Hessian_MATRIX(c_par);
        HSelection_A(params_A(c_par,1),c_par) = 1;
    end

    B = B_res; 
    SE_B = nan(size(Sigma,1));
    HSelection_B = zeros(M*M,size(params_B,1)); % selection matrix mapping theta elements to B

    for c_par = size(params_A,1)+1 : size(params_B,1)+size(params_A,1) 
        B(params_B(c_par-size(params_A,1),1))=StructuralParam_Estimation_MATRIX(c_par);
        SE_B(params_B(c_par-size(params_A,1),1)) = SE_Hessian_MATRIX(c_par);
        HSelection_B(params_B(c_par-size(params_A,1),1),c_par-size(params_A,1)) = 1;
    end

    for i = 1:M
        if B(i,i)<0, B(:,i)=-B(:,i); end % ensure diagonal of B is positive
        if A(i,i)<0, A(:,i)=-A(:,i); end % ensure diagonal of A is positive
    end

    logLik_SVAR_AB = -1*n_logLik;                      % convert negative log-likelihood to log-likelihood
    LR_test_overid = -2*(logLik_SVAR_AB - logLikVAR);  % LR test comparing structural vs reduced-form likelihood
    df = M*(M+1)/2 -size(params_A,1)-size(params_B,1); % degrees of freedom
    
    if df>0
        PVal = 1-chi2cdf(LR_test_overid,df); % p-value
        disp(['The p-value of the overidentification restrictions is ' num2str(PVal)])
    end 

%% compute IRFs
    J=[eye(M) zeros(M,M*(lags-1))]; % selects the first M rows of the companion form
    CompanionMatrix = [];           % initialize companion matrix container


    for p = 1:lags
        CompanionMatrix = [CompanionMatrix EstVAR.AR{p}]; % stack VAR coefficient matrices horizontally
    end 

    CompanionMatrix=[CompanionMatrix;eye(M*(lags-1)) zeros(M*(lags-1),M)]; % add lower block to complete companion form

    for h = 0 : HorizonIRF
        PHI(:,:,h+1)=J*CompanionMatrix^h*J'*pinv(A)*B;    % structural IRF: at horizon h to each 
    end 

    for h = 0 : HorizonIRF
        for i=1:M
            for j=1:M
                IRF(h+1,M*(i-1)+j)=PHI(i,j,h+1);          % flatten IRFs in row form 
            end
        end
    end 

%% Bootstrap IRFs

% (1) generate bootstrap data by resampling residuals
    Residuals_B=[];                                            % container for resampled residuals
    data_B=zeros(T+lags,M*repetitions);                        % container for bootstrap data matrix
    data_B(1:lags,:)=kron(ones(1,repetitions),data(1:lags,:)); % initialize bootstrap series with original values

    for boot = 1 : repetitions
        TBoot=datasample(1:T,T);                               % resample from 1 to T         
        Residuals_B=[Residuals_B Residuals(TBoot,:)];          % collect resampled residuals across repetitions
    end 

    for t=1+lags:T+lags
        data_h=zeros(1,size(data_B,2));                        % predicted value before adding bootstrap residuals
        for p=1:lags
            data_h=data_h+data_B(t-p,:)*kron(eye(repetitions),EstVAR.AR{p}'); % VAR prediction using previous values
        end 
        data_B(t,:)=data_h+kron(ones(1,repetitions),EstVAR.Constant')...      % add constant term
            +kron(ones(1,repetitions),EstVAR.Trend')*t...                     % add trend component
            +Residuals_B(t-lags,:);                                           % add bootstrap residuals
    end

% (2) split stacked bootstrap matrix into separate datasets
    data_B_all={}; 
    for boot = 1 : repetitions
        data_B_all{boot}=data_B(:,1+(boot-1)*M:M+(boot-1)*M);  
    end 

% (3) re-estimate in each bootstrap dataset and compute IRFs
    for boot = 1 : repetitions 
        disp(boot)
        data_B=data_B_all{boot};  % select bootstrap dataset

        [EstVAR_B] = estimate(VAR,data_B);  
        Sigma_B=EstVAR_B.Covariance;

        options = optimset('MaxFunEvals',200000,'TolFun',1e-500,'MaxIter',200000,'TolX',1e-100,'Display', 'off');   
        fun = @(theta) Likelihood_SVAR_AB(theta, A_res, B_res, Sigma_B, T, M, params_A, params_B);
        [StructuralParam_Estimation_MATRIX_B] = fminunc(fun, InitialValues, options);

        A_B = A_res;
        B_B = B_res;

        for c_par = 1 : size(params_A,1)
            A_B(params_A(c_par,1)) = StructuralParam_Estimation_MATRIX_B(c_par);
        end

        for c_par = size(params_A,1)+1 : size(params_B,1)+size(params_A,1)
            B_B(params_B(c_par-size(params_A,1),1))=StructuralParam_Estimation_MATRIX_B(c_par);
        end

        for i = 1:M
            if B_B(i,i)<0, B_B(:,i)=-B_B(:,i); end 
            if A_B(i,i)<0, A_B(:,i)=-A_B(:,i); end 
        end

        J=[eye(M) zeros(M,M*(lags-1))];
        CompanionMatrix_Boot = [];

        for p = 1:lags
            CompanionMatrix_Boot = [CompanionMatrix_Boot EstVAR_B.AR{p}];
        end

        CompanionMatrix_Boot=[CompanionMatrix_Boot;eye(M*(lags-1)) zeros(M*(lags-1),M)];

        for h = 0 : HorizonIRF
            PHI_Boot(:,:,h+1,boot)=J*CompanionMatrix_Boot^h*J'*pinv(A_B)*B_B;
        end
    end 
end 