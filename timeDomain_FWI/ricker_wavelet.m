function w = ricker_wavelet(t,t0,sigma)
w = (2/(sqrt(3*sigma)*pi^0.25))*(1-((t-t0)/sigma).^2).*exp(-(t-t0).^2/(2*sigma^2));
w = w/10;
%plot(w)
end