function w = blackmanharris_derivative1(N)

a0 = 0.35875;
a1 = 0.48829;
a2 = 0.14128;
a3 = 0.01168;

n = (0:N-1)';
w = a0 - a1*cos(2*pi*n/(N-1)) + a2*cos(4*pi*n/(N-1)) - a3*cos(6*pi*n/(N-1));

end