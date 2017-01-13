 function shi_06_v1

%
%Part1
%
%Average branch SNR
aSNR = 10;
%Number of value in data
N = 10^6;

%X22M probability density function numerical 
figure('name','assignment5');
xlabel('SNR(dB)');
ylabel('Pobility');
title('X22M probability density With Average branch SNR 10dB');
color={'r','g','b'};
SNRc = zeros (1,N);
hold on
for x=1:3
    %rayleigh gain
    g = 1./sqrt(2) .* (randn(1,N) + 1i * randn(1,N));
    %combing n channel
    SNRc = SNRc + (abs(g) .^ 2) .* (aSNR);
    %get the psd
    [Pc, SNR] = histR(SNRc,1000);
    Pc = Pc / N / diff(SNR(1 : 2));
    plot(10 * log10(SNR),Pc,color{x});
end

% X22M probability density function Theoretical
color={'r*','g*','b*'};
SNRdb = -20:1:40;
SNR = 10.^(SNRdb./10);
for L=1:3 
    % Compute Theoretical probability density
    Pc = SNR.^(L-1).*exp(-SNR./aSNR)./((aSNR).^L)./factorialR(L-1);
    plot(10 * log10(SNR),Pc,color{L});
end
xlim([0 40]);
legend('numerical M=1','numerical M=2','numerical M=1','Theoretical M=1','Theoretical M=2','Theoretical M=3')

hold off


%
%Part2
%



%Number of symbols
N = 10^6;
%Generate a 2bits/symbol signal
s = randi(4,1,N);
t = [1+1i -1+1i -1-1i 1-1i];
%Transform it in to4QAM random sequence
a = t(s);
% Eb/N0 in db
Eb0db = -10:2:40;
% Eb/N0
Eb0 = 10.^(Eb0db/10);
M=4;
%SNR
SNR = Eb0 * log2(M);
%The power of signal
Es = sum(abs(a).^2)/N;
N0 = Es./SNR;
% compute the sigma of noise by SNR
sigma = sqrt(N0);
v = (1./sqrt(2) .* (randn(1,N) + 1i * randn(1,N))); 
%build figuie label
figure('name','assignment5');

Constant=2;
r = sqrt(1.5* log2(M)/(M-1));
D = sqrt(Eb0./(1+Eb0));
%SER for 1 channel
ser = Constant*((1-D)./2);
semilogy(Eb0db,ser,'k:') ;
%SER for 2 channel
ser = Constant*((1-D)./2).^2.*(1+(2*((1+D)./2)));
hold on
semilogy(Eb0db,ser,'m:') ;
%SER for 3 channel
ser = Constant*((1-D)./2).^3.*(1+(3*((1+D)./2)+(6*((1+D)./2).^2)));
semilogy(Eb0db,ser,'c:') ;

%simulation SER
ser = zeros(1,length(Eb0));
g= zeros (1,N);

%Equivalent gain of MCR with L channel
n1 = sqrt(0.5).*(randn(1,N) + 1i * randn(1,N));
g1 =  sqrt((abs(n1) .^ 2) .* (1));
for l=1:length(Eb0);
    %produce noise
    v1 = sigma (l).* v;
    %add noise and gain
    a1 = a .* g1 + v1;
    %decision decoding
    s1 = 1*(real(a1)>0&imag(a1)>=0) + 2*(real(a1)<=0&imag(a1)>0) + 3*(real(a1)<0&imag(a1)<=0) + 4*(real(a1)>=0&imag(a1)<0);
    %SER 
    ser(l) = (sum(s1~=s))/N;
end

semilogy(Eb0db,ser,'r') ;

%Equivalent gain of MCR with L channel
n1 = sqrt(0.5).*(randn(1,N) + 1i * randn(1,N));
n2 = sqrt(0.5).*(randn(1,N) + 1i * randn(1,N));
g1 =  sqrt((abs(n1) .^ 2) + (abs(n2) .^ 2));
for l=1:length(Eb0);
    %produce noise
    v1 = sigma (l).* v;
    %add noise and gain
    a1 = a .* g1 + v1;
    %decision decoding
    s1 = 1*(real(a1)>0&imag(a1)>=0) + 2*(real(a1)<=0&imag(a1)>0) + 3*(real(a1)<0&imag(a1)<=0) + 4*(real(a1)>=0&imag(a1)<0);
    %SER 
    ser(l) = (sum(s1~=s))/N;
end

semilogy(Eb0db,ser,'r') ;

%Equivalent gain of MCR with L channel
n1 = sqrt(0.5).*(randn(1,N) + 1i * randn(1,N));
n2 = sqrt(0.5).*(randn(1,N) + 1i * randn(1,N));
n3 = sqrt(0.5).*(randn(1,N) + 1i * randn(1,N));
g1 =  sqrt((abs(n1) .^ 2) + (abs(n2) .^ 2) + (abs(n3) .^ 2));
for l=1:length(Eb0);
    %produce noise
    v1 = sigma (l).* v;
    %add noise and gain
    a1 = a .* g1 + v1;
    %decision decoding
    s1 = 1*(real(a1)>0&imag(a1)>=0) + 2*(real(a1)<=0&imag(a1)>0) + 3*(real(a1)<0&imag(a1)<=0) + 4*(real(a1)>=0&imag(a1)<0);
    %SER 
    ser(l) = (sum(s1~=s))/N;
end

semilogy(Eb0db,ser,'r') ;

xlim([-10 40])
ylim([0.00001 1])
xlabel('SNR(dB)');
ylabel('Error rate');
title('SER of MRC');
legend('The L=1','The L=2','The L=3','Sim L=1','Sim L=2','Sim L=3');
hold off



end

%factorial
function a=factorialR(n)
b = 1 : 1 : n;
a = prod (b);
end

%histgrom function%
function [no,xo] = histR(a,b) %replace command hist to get histogram. a is data. b is the total block number.
maxa=max(a);%get max value
mina=min(a);%get min value
l=(maxa-mina)/b;%get the length of each block
xo= zeros(1,b);
for i= 1:1:b
    xo(i)=mina+l*(i-1);%produce the map of block
end
no= zeros(1,b);
for j= 1:1:length(a)
    t=round((a(j)-mina)/l)+1;%run the data add count to each block
    if t>b
        t=b;
    end
     if t==0
        t=1;
    end
    no(t)=no(t)+1;
end
end
