%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPENN
% ECON 706, Metrics
% March 24, 2019.
% Playing with markov simulations
% By Rodrigo A Morales M :)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% https://www.mathworks.com/help/econ/dtmc.graphplot.html
% https://www.mathworks.com/help/econ/dtmc.simplot.html
clear;
clc;
cd /home/moralesmendozar/Dropbox/03_UPenn/classes/2019_fall/01_937-/PS02/Markov

N = 10;
m = rand(N,N);
%Generate the random matrix
%Generate the random matrix
P = (m ./( ones(N,1)*sum(m)))'; %Regular random matrix
P = zeros(N,N);
P(1,1) = 1;
P(N,N) = 1;
for i =2:(N-1)
    P(i,i-1) = rand();
    P(i, i+1) = 1-P(i,i-1);
end

mc = dtmc(P);
% mc = dtmc(P,'StateNames',['s01','s02','s03'])
figure(1)
graphplot(mc,'ColorNodes',true,'ColorEdges',true)

% Perform a simulation of 20 steps
numsteps = 20;
X = simulate(mc,numsteps);

figure(2)
simplot(mc,X)

% Perform three 10 step simulations for each state:
Nsim = 7;
%t0 = round(rand(Nsim,1)*(N)+0.5,0);
x0 = 100*ones(1,mc.NumStates); %states the number of simulations to begin in each state
x0 = [0 0 0 0 100 0 0 0 0 0];
%or
X0 = zeros(1,N);
X0(5) = 100;
numStep = 100;
X2 = simulate(mc,numStep, 'X0',x0);
%X2 = simulate(mc,numStep, 'X0',100*ones(1,N));

figure(3)
simplot(mc,X2)
%https://www.mathworks.com/help/econ/dtmc.graphplot.html