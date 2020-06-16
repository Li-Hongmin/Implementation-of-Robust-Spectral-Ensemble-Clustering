
%% config
addpath('./data')
addpath('./utils')
load iris_uni.mat
k = 3; % cluster number
N = 100; % base clustering number
%%
% generate base clusterings by kmeans
baseCls = RPS(X, k, N);

%% run RSEC
lambda1 = 0.1;
lambda2 = 0.01;

label = RSEC(baseCls, k, lambda1, lambda2);

%% evaluate

acc = Cal_ACC(Y, label)
nmi = Cal_NMI(Y, label)