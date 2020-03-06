dt = 1/1000;
%h = @(s,r) exp(-(s-r).^2/2);
h = @(s,r) (1-abs(s-r)).*(abs(s-r)<=1);
t = (-4:dt:4)';
f = h(t,0.7);
g = h(t,-0.7);
figure
plot(t,f,'r-',t,g,'b-')


[wass_dist,~,~,dist] = wasserstein_distance_1d(f,g,t);
% figure
% plot(dist)
wass_dist