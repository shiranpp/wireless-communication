 function shi_06_v2
 
 
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

%Theoretical SER

Constant=1.75;
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
% for one branch
    for l=1:length(Eb0);
        %rayleigh gain
        n = 1./sqrt(2) .* (randn(1,N) + 1i * randn(1,N));
        %ai for each branch
        g =  sqrt((abs(n) .^ 2) .* (SNR(l)));
        %produce noise
        v1 = sigma (l).* v;
        %add noise and gain
        a1 = a .* g.*abs(n) + v1.*g;
        %decision decoding
        s1 = 1*(real(a1)>0&imag(a1)>=0) + 2*(real(a1)<=0&imag(a1)>0) + 3*(real(a1)<0&imag(a1)<=0) + 4*(real(a1)>=0&imag(a1)<0);
        %SER 
        ser(l) = (sum(s1~=s))/N;
    end
semilogy(Eb0db,ser,'r') ;
% for two branch
    for l=1:length(Eb0);
        %rayleigh gain
        n1 = 1./sqrt(2) .* (randn(1,N) + 1i * randn(1,N));
        n2 = 1./sqrt(2) .* (randn(1,N) + 1i * randn(1,N));
        %ai for each branch
        g1 =  sqrt((abs(n1) .^ 2) .* (SNR(l)));
        g2 =  sqrt((abs(n2) .^ 2) .* (SNR(l)));
        %produce noise
        v1 = sigma (l).* (1./sqrt(2) .* (randn(1,N) + 1i * randn(1,N)));
        v2 = sigma (l).* (1./sqrt(2) .* (randn(1,N) + 1i * randn(1,N)));
        %add noise and gain
        a1 = a .* g1.*abs(n1) + v1.*g1+a .* g2.*abs(n2) + v2.*g2;
        %decision decoding
        s1 = 1*(real(a1)>0&imag(a1)>=0) + 2*(real(a1)<=0&imag(a1)>0) + 3*(real(a1)<0&imag(a1)<=0) + 4*(real(a1)>=0&imag(a1)<0);
        %SER 
        ser(l) = (sum(s1~=s))/N;
    end
semilogy(Eb0db,ser,'g') ;
% for 3 branch
    for l=1:length(Eb0);
        %rayleigh gain
        n1 = 1./sqrt(2) .* (randn(1,N) + 1i * randn(1,N));
        n2 = 1./sqrt(2) .* (randn(1,N) + 1i * randn(1,N));
        n3 = 1./sqrt(2) .* (randn(1,N) + 1i * randn(1,N));
        %ai for each branch
        g1 =  sqrt((abs(n1) .^ 2) .* (SNR(l)));
        g2 =  sqrt((abs(n2) .^ 2) .* (SNR(l)));
        g3 =  sqrt((abs(n3) .^ 2) .* (SNR(l)));
        %produce noise
        v1 = sigma (l).* (1./sqrt(2) .* (randn(1,N) + 1i * randn(1,N)));
        v2 = sigma (l).* (1./sqrt(2) .* (randn(1,N) + 1i * randn(1,N)));
        v3 = sigma (l).* (1./sqrt(2) .* (randn(1,N) + 1i * randn(1,N)));
        %add noise and gain
        a1 = a .* g1.*abs(n1) + v1.*g1+a .* g2.*abs(n2) + v2.*g2+a .* g3.*abs(n3) + v3.*g3;
        %decision decoding
        s1 = 1*(real(a1)>0&imag(a1)>=0) + 2*(real(a1)<=0&imag(a1)>0) + 3*(real(a1)<0&imag(a1)<=0) + 4*(real(a1)>=0&imag(a1)<0);
        %SER 
        ser(l) = (sum(s1~=s))/N;
    end
semilogy(Eb0db,ser,'b') ;

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
