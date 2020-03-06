function [quad_wass_dist,shift,scale,dist] = wasserstein_distance_1d(f,g,t)
%% Compute the quadratic wasserstein distance between 2 discrete vectors
if size(f,1) ~= size(g,1) || size(f,1) ~= size(t,1)
    error('Dimension of 2 input vectors not matching!\n')
end
n = size(f,1);
%% preprocessing -- removing negative values by shifting
shift = [-1.1*min(f.*(f<0));-1.1*min(g.*(g<0))];
f = f + shift(1);
g = g + shift(2);
%% trapziodal interpolation to compute rescaling
weight = ones(n,1);
weight(1) = 1/2; weight(n) = 1/2;
f_scale = sum(f.*weight);
g_scale = sum(g.*weight);
scale = [f_scale;g_scale];

f_rescale = f.*weight/f_scale;
%g_rescale = g.*weight/g_scale;
%% find the integral value at grid points
f_int = zeros(n,1);
g_int = zeros(n,1);
for j = 2:n
    weight = ones(j,1);
    weight(1) = 1/2; weight(j) = 1/2;
    f_int(j) = sum(f(1:j).*weight)/f_scale;
    g_int(j) = sum(g(1:j).*weight)/g_scale;
end
%% find the corresponding points that maps to roughly the same cumulative point
temp = 1;
dist = zeros(n,1);
for j1 = 1:n
    for j2 = temp:n
        if f_int(j1) <= g_int(j2)
            temp = j2;
            dist(j1) = (t(j2) - t(j1))^2;
            break;
        end
    end
end
%% compute quadratic wasserstein distance
quad_wass_dist = sum(dist.*f_rescale);
end