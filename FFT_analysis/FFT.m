function [f, P1] = FFT(vector, Fs)
% FFT is used to analysis the signal and get the frequency spectrum from
% the whole window;
L = length(vector);
Y = fft(vector);
P2 = abs(Y/L);
P1 = P2(1: round(L/2)+1);
P1(2 : end-1) = 2*P1(2 : end-1);
f = Fs*(0 : round(L/2))/L;

end