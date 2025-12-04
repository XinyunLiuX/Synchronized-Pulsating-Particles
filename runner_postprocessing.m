clc; clear all; close all;
folder = 'RES';
myfiles = dir(strcat(folder,'/res*.mat'));
parfor k = 1:length(myfiles)
    fileName = myfiles(k).name;
    pp_visualize_dynamics(folder, fileName)
end