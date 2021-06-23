% code to reproduce the simulation results in Umakantha, Morina, Cowley et al., (Neuron, 2021)

clear; clc; close all;

%% add paths
addpath(genpath('helper_functions'))
addpath(genpath('main_functions'))

%% generate figure 3f
fprintf('Generating Figure 3f...\n');
fig_3f
drawnow
clear

% colorbar indicates the loading similarity

%% generate figure 3g
fprintf('Generating Figure 3g...\n');
fig_3g
drawnow
clear

%% generate figure 4
fprintf('Generating Figure 4...\n');
fig_4
clear