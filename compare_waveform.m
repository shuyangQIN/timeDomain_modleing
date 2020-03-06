clear
clc
close all
addpath('data')
load('sensor_data_fdm.mat')
load('sensor_data_kwave.mat')
index = randsample(size(sensor_data_fdm,1),1,1);
plot(sensor_data_kwave(index,:))
hold on
plot(sensor_data_fdm(index,:))
hold off
legend('kwave','fdm')
