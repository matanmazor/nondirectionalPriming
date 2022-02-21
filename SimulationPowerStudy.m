close all;
clc;
clear;
warning('off')

%{
The point of this script is to identify cases where a nondirectional
analysis outperforms a directional analysis, and cases where the opposite
is true. 
In these simulations, we assume that the dependent variable (here we call it rt) is sampled from
one of two distributions: A or B (here termed incong and cong). A simplifying assumption here is that for
all subjects, the variance of these two distributions equals 1. The mean of
the A (incong) distribution is always 0. The mean of the B (cong) distribution, muB, differs
between subjects. Different simulations make different assumptions about
the distributions from which muB s sampled for different subjects.
Specifically, what is the mean of this distribution, and what is its
variance.
%}



% ANALYSIS PARAMS
params.N_perm = 25; %number of label shuffling per participant
params.N_null = 10000; %number of samples in bootstrapped null distribution
params.rng = 1;

params.N_splits = 100; %for sign consistency analysis
params.control_for = ''; 
params.predict = 'cong'; 
params.x = 'rt';
params.filter_column = '';
% params.inclusion_value = 1; %for filtering;

params.SVM = false;
params.signConsistency = true;
params.directional = true;

params.save=false;
params.plot=false;

params.statistic = @(x) mean(x);

% GENERATE DATA 
N_sub = 20;
N_trial = 100;

%manipulate_sd:
% population_sd = [0,exp(-3:0.5:-0.5)];
% population_mean = 0;

%manipulate mean:
population_sd = 0;
population_mean = 0;

N_iter = 1000;
consistency_p = nan(N_iter,length(population_sd));
directional_p = nan(N_iter,length(population_sd));

% for i_sd = 1:length(population_sd)
for i_mean = 1:length(population_mean)
    for i_iter = 1:N_iter

        i_iter
        rng(i_iter*3); %the analysis code always sets the seed to 1, so without this line all iterations produce identical data.
        subNum = [];
        cong = {};
        rt = [];
        Exp = {};

        for i_s = 1:N_sub

%             subj_effect = normrnd(population_mean,population_sd(i_sd));
            subj_effect = normrnd(population_mean(i_mean),population_sd);
            cong_s = binornd(1,0.5,N_trial,1);
            rt_s = normrnd(0,1,N_trial,1)+cong_s*subj_effect;

            subNum = [subNum; i_s*ones(N_trial,1)];
            for i_t=1:N_trial
               if cong_s(i_t)==1
                   cong{end+1}='cong';
               else
                   cong{end+1}='incong';
               end
               Exp{end+1}='sim';
            end

            rt = [rt; rt_s];
        end

        cong = cong';
        Exp = Exp';
        T = table(subNum,cong,rt,Exp);

        % ANALYZE
        resultsTable = analyzePriming(params,T);
        consistency_p(i_iter, i_mean)=resultsTable.consistency_p;
        directional_p(i_iter, i_mean)=resultsTable.directional_p;
    end
end

fig=figure;
hold on;
% plot(population_sd, mean(directional_p<0.05))
% plot(population_sd, mean(consistency_p<0.05))
plot(population_mean, mean(directional_p<0.05))
plot(population_mean, mean(consistency_p<0.05))
xlabel('Group-level mu');
xlim([0,population_mean(end)])
ylabel('Statistical power');
title('Statistical power for 100 trials x 20 participants')
legend('directional','nondirectional', 'Location','southeast');

s=hgexport('readstyle','presentation');
s.Format = 'png';
s.Width = 7;
s.Height = 3;

hgexport(fig,fullfile('analyzed','simulation_power_study','power_by_mu.png'),s);

