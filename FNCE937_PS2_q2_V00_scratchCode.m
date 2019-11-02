% For PS02, ex 02. Beginning

% matrix coefficients contains results:
%
% Housekeeping
clear
clc
%%
vRhos = [0.6 0.7 08];%[0.5 0.6 0.7 0.8 0.9];
% vRhos = [0.5 0.6 0.7 0.8 0.9];
vLambdas = [0.01 0.025 0.05]; %[ 0.01 0.02 0.025 0.03 0.05];
% vLambdas = [ 0.01 0.02 0.025 0.03 0.05];
numRhos = length(vRhos);
numLambdas = length(vLambdas);
numCoefs = 9; % rho, lambda, beta0-beta2, gamma0-gamma3
coefficients = zeros(numRhos*numLambdas,numCoefs);
ii = 0;
j = 0;
for rrho= vRhos
    ii = ii+1;
    coefficients(((ii-1)*numLambdas +1):(ii*numLambdas),1) = rrho;
    j = 0;
    for llambda = vLambdas
        j = j+1;
        coefficients((ii-1)*numLambdas+j,2) = llambda;
    end
end


%% Get the coefficients from code of 1.d)
ii = 0;
j = 0;
for rrho= vRhos
    ii = ii+1;
    j = 0;
    for llambda = vLambdas
        j = j+1;
        
        vectorCoefs = 1:7; % beta0-beta2, gamma0-gamma3
        coefficients((ii-1)*numLambdas+j,3:numCoefs) = vectorCoefs;
    end
end
coefficients