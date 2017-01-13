function ser = shi_05_v1
N = 10^6;%Number of symbols
s = randi(4,1,N);%Generate a 2bits/symbol signal
t = [1+1i -1+1i -1-1i 1-1i];
a = t(s);%Transform it in to4QAM random sequence
Eb0db = 0:1:40;% Eb/N0 in db
Eb0 = 10.^(Eb0db/10);% Eb/N0
M=4;

%No fading
Es = sum(abs(a).^2)/N;%The power of signal
N0 = Es./(2*Eb0);
sigma = sqrt( N0 /2);% compute the sigma of noise by SNR
no = randn(1,length(a)) + randn(1,length(a)) * 1i;%produce the basic noise
ser = zeros(1,length(Eb0db));
for l = 1:length(Eb0db)
    no1 = sigma (l).* no;%produce noise
    a1 = a + no1;%add noise and gain
    s1 = 1*(real(a1)>0&imag(a1)>=0) + 2*(real(a1)<=0&imag(a1)>0) + 3*(real(a1)<0&imag(a1)<=0) + 4*(real(a1)>=0&imag(a1)<0);
    % decide the received signal
    ser(l) = (sum(s1~=s))/N;%compute the SER
end
figure('name','assignment5')
semilogy(Eb0db,ser,'b') 
hold on
xlim([0 40])
ylim([10^-5 1])


%Theoretical No fading
r = sqrt(1.5* log2(M)/(M-1));
ser = 4*(1-1/sqrt(M))*qfunR(sqrt(2)*r*sqrt(Eb0));
semilogy(Eb0db,ser,'b*');

%Rayleigh fading K = 0
ser = zeros(1,length(Eb0db));
g = (1/sqrt(2))*(randn(1,length(a))+1i*randn(1,length(a)));
for l = 1:length(Eb0db)
    no1 = sigma (l).* no;%produce noise
    a1 = a.*abs(g) + no1;%add noise and gain
    s1 = 1*(real(a1)>0&imag(a1)>=0) + 2*(real(a1)<=0&imag(a1)>0) + 3*(real(a1)<0&imag(a1)<=0) + 4*(real(a1)>=0&imag(a1)<0);
    % decide the received signal
    ser(l) = (sum(s1~=s))/N;%compute the SER
end
semilogy(Eb0db,ser,'g') 

%Theoretical Rayleigh fading K = 0
d = 0.01:0.01:pi/2;
ser = zeros(1,length(Eb0db));
alpha = sqrt(2)*r*sqrt(Eb0);%alpha of 4qam with EB/N0
for l = 1:length(Eb0db)% Do integral 
    mg = ((1+(alpha(l)^2)./(2*sin(d).^2))).^-1;%Mg
    ser(l) = sum(mg.*0.02)/pi;
end
semilogy(Eb0db,ser,'g*');

%Rician fading with K = 3;
K = 3;
sigemalos = sqrt(K/(K+1));% sigma of LOS
sigemanlos = sqrt(0.5/(K+1));;% sigma of NLOS
ser = zeros(1,length(Eb0db));
g = sigemalos + sigemanlos*(randn(1,length(a))+1i*randn(1,length(a)));
for l = 1:length(Eb0db)
    no1 = sigma (l).* no;%produce noise
    a1 = a.*abs(g) + no1;%add noise and gain
    s1 = 1*(real(a1)>0&imag(a1)>=0) + 2*(real(a1)<=0&imag(a1)>0) + 3*(real(a1)<0&imag(a1)<=0) + 4*(real(a1)>=0&imag(a1)<0);
    % decide the received signal
    ser(l) = (sum(s1~=s))/N;%compute the SER
end
semilogy(Eb0db,ser,'r') 

%Theoretical Rician fading with K = 3
K = 3;
d = 0.01:0.01:pi/2;
ser = zeros(1,length(Eb0db));
for l = 1:length(Eb0db)
    x = alpha(l)^2./(2.*sin(d).^2);% Mg is too complex so i take part of it here.
    mg = (1+K)./(1+K+x).*exp(-(K.*x)./(1+K+x));%Mg
    ser(l) = sum(mg.*0.02)/pi;
end
semilogy(Eb0db,ser,'r*');

%Rician fading with K = 10;
K = 10;
sigemalos = sqrt(K/(K+1));% sigma of LOS
sigemanlos = sqrt(0.5/(K+1));% sigma of NLOS
ser = zeros(1,length(Eb0db));
g = sigemalos + sigemanlos*(randn(1,length(a))+1i*randn(1,length(a)));
for l = 1:length(Eb0db)
    no1 = sigma (l).* no;%produce noise
    a1 = a.*abs(g) + no1;%add noise and gain
    s1 = 1*(real(a1)>0&imag(a1)>=0) + 2*(real(a1)<=0&imag(a1)>0) + 3*(real(a1)<0&imag(a1)<=0) + 4*(real(a1)>=0&imag(a1)<0);
    % decide the received signal
    ser(l) = (sum(s1~=s))/N;%compute the SER
end
semilogy(Eb0db,ser,'k') 
hold on

%Theoretical Rician fading with K = 10
K = 10;
d = 0.01:0.01:pi/2;
ser = zeros(1,length(Eb0db));
for l = 1:length(Eb0db)
    x = alpha(l)^2./(2.*sin(d).^2);% Mg is too complex so i take part of it here.
    mg = (1+K)./(1+K+x).*exp(-(K.*x)./(1+K+x));%Mg
    ser(l) = sum(mg.*0.02)/pi;
end
semilogy(Eb0db,ser,'k*');


xlim([0 40])
ylim([10^-5 1])
legend('No fading','Theoretical No fading','Rayleigh fading K = 0','Theoretical Rayleigh fading K = 0','Rician fading with K = 3','Theoretical Rician fading with K = 3','Rician fading with K = 10;','Theoretical Rician fading with K = 10')
end

function q = qfunR( x)
q = 0.5 * erfc (x/sqrt(2));
end