% For PS02, ex 02.
% Wharton @ Penn.
%   Shasha and Rodrigo M.

% matrix coefficients contains results:
%
% Housekeeping
loadedData = 1;
if loadedData ==1
    Vinit = Vrenewed;
else
    clear
end
clc
%someParameters
M = 0.99;
ddelta = 0.1;
aalphaK = 0.3;
aalphaL = 0.6;
W = 2; % wage
ttaoC=0.15; % corporate tax rate
Nk = 20;
kMin = 0.00001;
kSteadyState = 1;
% a_grid and a's transition matrix
m = 3; % parameter for tauchen method
Na = 5;
ssigma = 0.05;
Nb = 15;
tolerance=0.00001;
kGridLength = [15]; % number of points in grid for capital
numsteps = 1000;
%checking Tauchen matrix
% r = 1/M-1; % interest rate for notation convenience
% RRHO = 0.7;
% tempK = kSteadyState^((aalphaK+aalphaL-1)/(1-aalphaL));
% aMean = ((r+ddelta)/aalphaK/tempK)^(1-aalphaL)*(W/aalphaL)^aalphaL; 
% [grid_a_log,m_a_prob] = tauchen(Na,log(aMean),RRHO,ssigma,m);
%buildOtherCoefficients to use in the function calling the simulation
vOtherCoefs= [M ddelta aalphaK, aalphaL, W, ttaoC, Nk, kMin, kSteadyState];
vOtherCoefs= [vOtherCoefs m Na ssigma Nb tolerance]; %size 14
if loadedData ==1
    Vinit = Vrenewed;
else
    Vinit = ones(kGridLength(1),Nb,Na,Na); % initial guess
end
%First make sure that this is working for rho = 0.7 , lambda = 0.025
%[vectorCoefs, Vnew] = funVectorCoefs02(0.7,0.025, Vinit,vOtherCoefs,...
%    kGridLength,1, numsteps);
%vectorCoefs
%% check it should be quick!! :)
%[vectorCoefs, Vnew] = funVectorCoefs02(0.7,0.025, Vnew,vOtherCoefs,...
%    kGridLength,0, numsteps);
%%
%vRhos = 0.7*ones(1,3);
%vRhos = [0.6 0.7 08];%[0.5 0.6 0.7 0.8 0.9];
vRhos = [0.5 0.6 0.7 0.8 0.9];
%vLambdas = 0.025*ones(1,3);
%vLambdas = [0.01 0.025 0.05]; %[ 0.01 0.02 0.025 0.03 0.05];
vLambdas = [ 0.01 0.02 0.025 0.03 0.05];
numRhos = length(vRhos);
numLambdas = length(vLambdas);
numCoefs = 9; % rho, lambda, beta0-beta2, gamma0-gamma3
coefficients = zeros(numRhos*numLambdas,numCoefs);
%Generate the skeleton of the matrix coefficients:
 % [rho1, lambda1, [beta0-beta2, gamma0-gamma3];
    %[ rho1, lambda2, [beta0-..-gamma3]; ...
    %       ...]
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


%% Get the coefficients from code of 1.d) for each rho/lambda
ii = 0;
j = 0;
for rrho= vRhos
    ii = ii+1;
    j = 0;
    for llambda = vLambdas
        j = j+1;
        % HERE WE RUN THE CODE FROM 1.d with the given llambda and rho..
        %PUNCHLINE OF THIS EXERCISE, GET THE VECTOR OF BETAs..
        % call using a function:
        [vectorCoefs, Vnew] = funVectorCoefs02(rrho,llambda, Vinit,vOtherCoefs,...
            kGridLength,0, numsteps); %vectorCoefs = 1:7;
        Vinit = Vnew;
        coefficients((ii-1)*numLambdas+j,3:numCoefs) = vectorCoefs;
    end
end
%coefficients

%once we have the coefficients, we just see which one is the closest one
%       to the one given
betaObjective = [0.004 0.02 -0.15 -0.4 0.05]; %b1, b2, g1,g2,g3
%get rid of beta0 and gamma0
coefToCompare = coefficients(:,[4,5,7,8,9]);
errors = zeros(numRhos*numLambdas,1);
% find the best coefficients w.r.to the norm2 distance vs the given coefs.
errors(1) = norm(coefToCompare(1,:) - betaObjective,2);
minerror = errors(1);
minerrorii = 1;
for ii=2:(numRhos*numLambdas)
    errors(ii) = norm(coefToCompare(ii,:) - betaObjective,2);
    if errors(ii) < errors(ii-1)
        minerror = errors(ii);
        minerrorii = ii;
    end
end
display(['The best coefficients found are = [rho,lambda,betas,gammas] =']);
display([num2str(coefficients(minerrorii,:))])

saveData = 1;
if saveData ==1
    filename = '00_data_ex03_01';
    save(filename)
end
