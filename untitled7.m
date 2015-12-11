clear all
close all
%%% Example of convolution theorem for finite support signals
nx = 40;
x = linspace(0,2*pi,nx)';
k = linspace(-nx/2,nx/2 - 1,nx)';
u = sin(x);
upad = [u;zeros(length(u) - 1,1)];
uhat = fftshift(fft(u))/nx;
uhatpad = fftshift(fft(upad));

Gf = heaviside(pi/0.1 - abs(k));
%%%% Convolution in frequency domain == fft of multiplication in physical
%%%% domain. Relavent for u^2 term

uhat2 = fftshift(fft(u.^3)/nx);
uhat2_freq = myconv(uhat,myconv(uhat,uhat,nx),nx); %myconv =  1./nx*fftshift(cconv(uhat,what,nx));
norm(uhat2_freq - uhat2)

%%% Now opposite. fft of convolution in physical domain = mult in freq
%%% domain. 
uc = fftshift(fft(conv(u,u,'same'))/nx);
uc2 = uhat.*uhat;
%uc = conv(u,u,'same');
%uc2 = ifft(uhatpad.*uhatpad);
plot(abs(uc))
hold on
plot(abs(uc2))
