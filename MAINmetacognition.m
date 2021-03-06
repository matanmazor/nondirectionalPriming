close all;
clc;
clear;
warning('off')

% Stelzer, J., Chen, Y., & Turner, R. (2013). 
% Statistical inference and multiple testing correction in 
% classification-based multi-voxel pattern analysis (MVPA): 
% random permutations and cluster size control. Neuroimage, 65, 69-82.
params.N_perm = 25; %number of label shuffling per participant
params.N_null = 10000; %number of samples in bootstrapped null distribution
params.rng = 1;

params.N_splits = 500; %for sign consistency analysis
params.control_for = ''; 
params.predict = 'Stimulus'; 
params.x = 'Confidence';
params.accuracy = 'Accuracy';

params.signConsistency = true;
params.directional = true;

params.statistic = @(conf,accuracy) getAUC(conf,accuracy);
filename = 'data_Mazor_2020_Discrimination';

analyzeMetacognition(params,filename);



