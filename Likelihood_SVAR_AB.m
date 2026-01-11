function [n_logLik]=Likelihood_SVAR_AB(theta, A_res, B_res, Sigma, T, M, params_A, params_B)

% OTPUTS
% n_logLik = negative log likelihood

% INPUTS
% theta =vector containing all free parameters (A first, then B)
% A_res = A matrix restrictions
% B_res = B matrix restrictions
% T = period length
% M = number of variables
% params_A = indices of free entries in A
% params_B = indices of free entries in B

%% fill A with free parameters
    A = A_res;
    for c_par = 1 : size(params_A,1)
        A(params_A(c_par,1))=theta(c_par);
    end 

%% fill B with free parameters 
    B = B_res;
    for c_par = size(params_A,1)+1 : size(params_B,1)+size(params_A,1)
        B(params_B(c_par-size(params_A,1),1))=theta(c_par);
    end 

%% negative log-likelihood
    logLik = -0.5*T*M*(log(2*pi))...
             - 0.5*T*log((det(pinv(A)*(B*B')*pinv(A)')))...
             -0.5*T*trace((A'*pinv(B)'*pinv(B)*A)*Sigma); 
    n_logLik = -logLik;
end
