function Q03 % i do the PSD of Doppler function and produce the Doppler filter with it. then i filter a white nose with it
%produce the psd and Corresponding correlation
fd=10;
fs=1000;
fx=-fs/2:1:fs/2;
S=zeros(1,length(fx));
n = randn(1,100000);%white Noise with sigama 1 mu 0
for i=-fs/2:1:fs/2
    if abs(i)<fd
        S(i+fs/2+1)=1/sqrt(1-(i/fd)^2);%Produce the PSD of Doppler function
    else
        S(i+fs/2+1)=0;
    end
end
S1=circshift(S,ceil(size(S)/2));
S2=sqrt(S);%frequency respose of Doppler filter
h=ifft(sqrt(S1))%Doppler filter
h=circshift(h,floor(size(h)/2));
T=length(h)/fs*1000;%time table in ms
fx1=T/length(h):T/length(h):T;
figure('name','Q031');
plot(fx1,h)
title('Impulse response of Doppler filter;')
xlabel('Time(ms)')
ylabel('signal magnitude');
fx2=-fs/2:1:fs/2;
figure('name','Q032');
plot(fx2,S2)
title('Frequency response of Doppler filter;')
xlabel('Frequency')
ylabel('signal magnitude');
y=filter(h,1,n); %filter the noise with doppler filter
m=100;%block number
N=length(y);
P=N/m;
S3=zeros(1,P)
for j=0:1:m-1
    S3=S3+ (1/(P*fs))*(abs(fft(y(1+j*P:(j+1)*P))).^2);%fft for each blcok
end
S3=S3/m;
S4=circshift(S3,floor(size(S3)/2));% PSD of the filtered noise
fx3=-fs/2:fs/(length(S4)-1):fs/2;
figure('name','Q033');
plot(fx3,S4)
title('PSD of Gain;')
xlabel('Frequency')
ylabel('Power/Frequency');
C=ifft(S3);% Corresponding correlation of filted signale

T=P/fs*1000;%time table in ms
fx4=T/length(C):T/length(C):T;
figure('name','Q033');
plot(fx4,C)
title('Corresponding correlation function')
xlabel('Time(ms)')
ylabel('correlation magnitude');

end
