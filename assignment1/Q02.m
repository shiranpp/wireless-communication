function [sigma p]=Q02
k=1.380000000000000e-23;
f=3.1623;
t=300;
n0=k*t*f;%we have power density of noise here
sigma=sqrt(n0);%sigma is the Square root of power density

fs=60000;%sampling frequency
wp=10000/fs*2*pi;% pass band in angle
ws=12000/fs*2*pi;% stop band in angle
wt = ws - wp;%trans band of window
N = ceil(12*pi/wt);% level of window
wc = (ws+wp)/2;
h = 100*fir1(N,wc/pi,blackman(N+1));% i use the blackman window in creat a filter. it can't reach 50db at 12khz but it is best i can find
N1=fs/10;%window size of fft with resolution of 10hz
hdft = fft(h,N1);% fouier transform
hdft = hdft(1:N1/2+1);% show the half part of filter 
freq = fs/2/length(hdft):fs/2/length(hdft):fs/2% the X axis of filter spectrum

figure('name','Q021');
title('Filter window')
plot(freq,10*log10(abs(hdft))); 
xlabel('Frequency/hz') ;
ylabel('/dB') ;

v = randn(100000,1) .* sigma;% produce the noise
y=filter(h,1,v);% filter the noise
N1=fs/10;%window size of fft with resolution of 10hz
ydft = fft(y,N1);% fouier transform
ydft = ydft(1:N1/2+1);
psdy = (1/(N1)) * abs(ydft).^2;% do the power density
psdy(2:end-1) = 2*psdy(2:end-1);
freq = fs/2/length(psdy):fs/2/length(psdy):fs/2% the X axis of filter spectrum

figure('name','Q022');
plot(freq,10*log10(psdy))
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)');

p=sum(abs(y).^2);

end
