clear
clc
close all
addpath('../source_func')
%%
cmap = customize_colormap(2);
load('Mar_v_test.mat')
imagesc(marmousi_velocity)
colormap(cmap)
colorbar