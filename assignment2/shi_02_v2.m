function shi_02_v1
Q01;
Q02;
Q03;
word= 'Question1: \nfigure(Q01) \n';
fprintf(word);
word= 'Question2: \nfigure(Q02) \n';
fprintf(word);
word= 'Question3: \nfigure(Q031) \nfigure(Q032) \nfigure(Q033) \nfigure(Q034) \n';
fprintf(word);
end
function Q01 %i produce ydb in this queation and change it to y. and get min power by SNR.
% then i do a 1km-40km loop and do the outage prob for each distant by
% mix y with exp loss for each y. and produce the theoretical answer to
% cmpare with it.
K = 1.380650e-23; % Boltzmann's constant
T = 300; %Temperature
c = 2.99792458 * 10^8;%light speed
B = 10^7;%Band width
Pt = 2;%Transmit power
Fc = 10^9;%Carrier Frequency
Gt = 10;%Tansmiter gain
wl = c/Fc;%Wave lenghth
Gr = 1;%Receiver gain
F = 10^0.5;%Noise figure
mu = 0; %value of mu
sigma = 6; %value of sigma
D = 1:1:40;
Pn = F*K*T*B;%noise power
Prmin = 10*Pn;%min power
D0 = 1000;
k = 1/(((4*pi*D0)/wl)^2); % path loss at 1km with free space model
ydb = randn(10000,1) .* sigma + mu;%produce ydb 
y = 10.^(ydb./10);
po1 = zeros(1,40);
po2 = zeros(1,40);
Pr1=0;
Pr2 =0;
for i = 1:1:40
    Pr2 =0;
    for j = 1:1:length(y)
        Pr1 = Pt*Gt*Gr*k*((D0/(D(i)*1000))^3.5)*y(j); % Do the received power for each y
        Pr2 = Pr2+Pr1;
        if Pr1<Prmin
            po1(i) = po1(i)+1; %find the outage
        end
    end
    Pr2=Pr2/length(y);
    po2(i)=Qf((10*log10(Prmin)-10*log10(Pt*Gt*Gr*k*((D0/(D(i)*1000))^3.5)))/sigma);% do the power covrage probility
    po2(i)=1-po2(i);
end
po1=po1./length(y);
figure('name','Q01');
plot(D,po1,D,po2,'b--o');
title('Outage Probility With distant form 1-40km')
xlabel('Distant(Km)')
ylabel('Outage Probility');
end

function Q02 %i do the coverage for (Z,Z) point of the area. and compute the total coverage for those point with in radius R.
% and i do Z point in integral of analytical formula computation.
% i keep the Z small to reduce the compute time.But increase Z will have more accurately result
K = 1.380650e-23; % Boltzmann's constant
T = 300; %Temperature
c = 2.99792458 * 10^8;%light speed
B = 10^7;%Band width
Pt = 2;%Transmit power
Fc = 10^9;%Carrier Frequency
Gt = 10;%Tansmiter gain
wl = c/Fc;%Wave lenghth
Gr = 1;%Receiver gain
F = 10^0.5;%Noise figure
mu = 0; %value of mu
sigma = 6; %value of sigma
D = 3:1:15;
Pn = F*K*T*B;%noise power
Prmin = 10*Pn;%min power
D0 = 1000;
k1 = 1/(((4*pi*D0)/wl)^2);
ydb = randn(1000,1) .* sigma + mu;%produce ydb 
y = 10.^(ydb./10);
pc1=0;
pc2=zeros(1,13);
pc3=zeros(1,13);
R=0;
v=0;
r=0;
r1=0
Z=100;% block number increase this to get more More accurately result
h=0;
Q=0;
x1=0;
y1=0;
for i = 3:1:15;
    r=i*1000;
    v=0;
    for j = 1:1:Z
        for k = 1:1:Z
            x1=j/Z*r;
            y1=k/Z*r;
            R=sqrt(x1^2+y1^2);
            pc1=0;
            if R<=r
                for l = 1:1:length(y)
                    Pr1 = Pt*Gt*Gr*k1*((D0/(R))^3.5)*y(l); % Do the received power for each y within the N*N block and within radius r
                    if Pr1>Prmin
                        pc1 = pc1+1;
                    end
                end
                v=v+1;
                pc2(i-2)=pc2(i-2)+pc1/length(y);%add up the signal coverage for each distant
            end
        end
    end
    pc2(i-2)=pc2(i-2)/v;% do the total coverage 
    h=(r-0)/Z;
    Q=0;
    for j = 0:1:Z-1% do the integral for the Qfunction
        r1=0+j*h+h/2;
        Q=Q+Qf((10*log10(Prmin)-10*log10(Pt*Gt*Gr*k1*((D0/r1)^3.5)))/sigma)*r1;% do the power covrage
    end
    pc3(i-2)=Q*h/r/r*2;
end
figure('name','Q02');
plot(D,pc2,D,pc3,'b--o');
title('Signal coverage With distant form 3-15km')
xlabel('Distant(Km)')
ylabel('Coverage ');
end
function Q03 % i do the PSD of Doppler function and produce the Doppler filter with it. then i filter a white nose with it
%produce the psd and Corresponding correlation
fd=30;
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
h=circshift(h,floor(size(h)/2))
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
m=500;%block number
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
C=circshift(C,ceil(size(C)/2))
T=P/fs*1000;%time table in ms
fx4=T/length(C):T/length(C):T;
figure('name','Q033');
plot(fx4,C)
title('Corresponding correlation function')
xlabel('Time(ms)')
ylabel('correlation magnitude');

end

function Q=Qf(x) %my Qfunction
Q=0.5-0.5*erf(x/sqrt(2))
end

