function shi_06_v2
aSNR = 10;%Average branch SNR
Pos = 1;%Power of signal
N = 10^6;%Number of value in compute

%X22M probability density function numerical 
SNRc = zeros (1,N);
hold on
for x=1:3
g = 1./sqrt(2) .* (randn(1,N) + 1i * randn(1,N)); 
SNRc = SNRc + (abs(g) .^ 2) .* (aSNR);
[Pc, SNR] = histR(SNRc,1000);
Pc = Pc / N / diff(SNR(1 : 2));
plot(10 * log10(SNR),Pc);
end

%X22M probability density function Theoretical
SNRdb = -20:1:40;
SNR = 10.^(SNRdb./10);
for M=1:3
Pc = SNR.^(M-1).*exp(-SNR./aSNR)./((aSNR).^M)./factorialR(M-1);
plot(10 * log10(SNR),Pc);
end
hold off
end

%factorial
function a=factorialR(n)
b = 1 : 1 : n;
a = prod (b);
end

%histgrom function%
function [no,xo] = histR(a,b) %replace command hist to get histogram. a is data. b is the total block number.
maxa=0;
mina=0;
maxa=max(a);%get max value
mina=min(a);%get min value
l=0;
l=(maxa-mina)/b;%get the length of each block
xo= zeros(1,b);
for i= 1:1:b
    xo(i)=mina+l*(i-1);%produce the map of block
end
no= zeros(1,b);
t=0;
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
