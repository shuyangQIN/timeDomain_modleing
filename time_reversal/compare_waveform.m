clear
clc
close all
addpath('data')
load('sensor_data_fdm11.mat')
load('sensor_data_kwave1.mat')
index = randsample(size(sensor_data_fdm11,1),1,1);
plot(sensor_data_kwave1(index,:))
hold on
plot(sensor_data_fdm11(index,:))
hold off
legend('kwave','fdm')
