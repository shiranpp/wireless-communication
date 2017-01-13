function shi_04_v1
clear all;
T=200;%total Times , increase to get more accurate result
fs=2000;
s = randi(4,1,fs*T);%generate a 2bits/symbol signal
t = [1+1i -1+1i -1-1i 1-1i];
a = t(s);%transfor it in to4QAM random sequence
SNR = 0:1:60;% Eb/N0 in db


%QAM
Eb = sum(abs(a).^2)/T/fs;% compute the power of signal
N0 = Eb./10.^(SNR/10);
sigma = sqrt( N0 /2);% compute the sigma of noise by SNR
no = randn(length(sigma),length(a)) + randn(length(sigma),length(a)) * 1i;
no = no.*repmat(sigma,length(a),1)';%complex noise
a1 = repmat(a,length(sigma),1)+no;% add the noise to signal
%make the reciver signal decision
s1 = 1*(real(a1)>0&imag(a1)>=0) + 2*(real(a1)<=0&imag(a1)>0) + 3*(real(a1)<0&imag(a1)<=0) + 4*(real(a1)>=0&imag(a1)<0);
ser = s1 -  repmat(s,length(sigma),1);
ser = (0~=ser);
ser = (sum(ser'))./fs/T;%Symbol error rate
figure('name','assignment4')
semilogy(SNR,ser,'yellow') 
hold on


%DIF QAM
c = [1 + 1i a/(1+1i)];
for l=2:length(c)
c(l) = c(l-1)*c(l);%generate the DIF signal of orginal signal
end
Eb = sum(abs(c).^2)/T/fs;% compute the power of signal
N0 = Eb./10.^(SNR/10);
sigma = sqrt( N0 /2);% compute the sigma of noise by SNR
no = randn(length(sigma),length(c)) + randn(length(sigma),length(c)) * 1i;
no = no.*repmat(sigma,length(c),1)';%complex noise
c1 = repmat(c,length(sigma),1)+no;% add the noise to signal
a1 = c1(:,2:length(c1))./c1(:,1:length(c1)-1).*(1+1i);% restore the orginal signal form DIF
%make the reciver signal decision
s1 = 1*(real(a1)>0&imag(a1)>=0) + 2*(real(a1)<=0&imag(a1)>0) + 3*(real(a1)<0&imag(a1)<=0) + 4*(real(a1)>=0&imag(a1)<0);
ser = s1 -  repmat(s,length(sigma),1);
ser = (0~=ser);
ser = (sum(ser'))./fs/T;%Symbol error rate
semilogy(SNR,ser,'magenta') 


%DIF QAM with rayleigh
color = {'cyan'  'red'  'green'};
fd = [20 60 200];%let fd = 20,60,200
for l=1:3
    g = genG( fd(l) , length(c));%generate the gain with doppler frequency fd
    c1 = repmat(c,length(sigma),1).*repmat(g,length(sigma),1)+no;% add noise and gain to signal
    a1 = c1(:,2:length(c1))./c1(:,1:length(c1)-1).*(1+1i);% restore the orginal signal form DIF
    %make the reciver signal decision
    s1 = 1*(real(a1)>0&imag(a1)>=0) + 2*(real(a1)<=0&imag(a1)>0) + 3*(real(a1)<0&imag(a1)<=0) + 4*(real(a1)>=0&imag(a1)<0);
    ser = s1 -  repmat(s,length(sigma),1);
    ser = (0~=ser);
    ser = (sum(ser'))./fs/T;%Symbol error rate
    semilogy(SNR,ser,color{l}) 
end


%theoretical ser in 4QAM with rayleigh fading
r=1.5*log2(4)/(4-1);
ser=3./(4.*r.*(10.^(SNR/10))); % for 4 QAM rayleigh fading ser 
semilogy(SNR,ser,'blue') ;
title('Symbol error rate of Rayleigh fading')
xlabel('SNR(dB)')
ylabel('Ser');
legend('QAM','DIF QAM','DIF QAM 20hz fd','DIF QAM 60hz fd','DIF QAM 200hz fd','QAM rayleigh theoretical ser')
end


function g = genG( fd , n1) %generate the gain with doppler frequency and sampling frequency and time T
fs=2000;
fx = - fs / 2 : 1 : fs / 2;
n = randn( 1 , n1 ) + randn( 1 , n1 ) * 1i;%white Noise with sigama 1 mu 0
S=zeros(1,length(fx));
for l=-fs/2:1:fs/2
    if abs(l)<fd
        S(l+fs/2+1)=6000/4/pi/fd/sqrt(1-(l/fd)^2);%Produce the PSD of Doppler function
    else
        S(l+fs/2+1)=0;
    end
end
S1 = circshift(S,ceil(size(S)/2));
h = ifft(sqrt(S1));%Doppler filter
h = circshift(h,floor(size(h)/2));
g = filter ( h , 1 , n );
end