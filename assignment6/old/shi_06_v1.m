function Pc=shi_06_v1

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
xlabel('SNR(dB)');
ylabel('Error rate');
title('SER of MRC');

%Theoretical SER
%since we don't know what is constant in equation 
%i assume it is 2 because it have same ser with rayleigth fading in -
%-signal channel
r = sqrt(1.5* log2(M)/(M-1));
D = sqrt(Eb0./(1+Eb0));
%SER for 1 channel
ser = 2*((1-D)./2).*((1+D)./2);
semilogy(Eb0db,ser,'k:') ;
%SER for 2 channel
ser = 2*((1-D)./2).^2.*((3*((1+D)./2).^2));
hold on
semilogy(Eb0db,ser,'k*') ;
%SER for 3 channel
ser = 2*((1-D)./2).^3.*(10.*((1+D)./2).^3);
semilogy(Eb0db,ser,'k+') ;

%simulation SER
ser = zeros(1,length(Eb0));
g= zeros (1,N);
color={'r','g','b'};
for L=1:3
    %rayleigh gain
    n = 1./sqrt(2) .* (randn(1,N) + 1i * randn(1,N));
    %Equivalent gain of MCR with L channel
    g =  sqrt(g.^2+(abs(n) .^ 2) .* (1));
    for l=1:length(Eb0);
        %produce noise
        v1 = sigma (l).* v;
        %add noise and gain
        a1 = a .* g + v1;
        %decision decoding
        s1 = 1*(real(a1)>0&imag(a1)>=0) + 2*(real(a1)<=0&imag(a1)>0) + 3*(real(a1)<0&imag(a1)<=0) + 4*(real(a1)>=0&imag(a1)<0);
        %SER 
        ser(l) = (sum(s1~=s))/N;
    end
semilogy(Eb0db,ser,color{L}) ;
end


xlim([-10 40])
ylim([0.00001 1])

hold off


end
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

function q = qfunR( x)
q = 0.5 * erfc (x/sqrt(2));
end
