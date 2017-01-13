function shi_01_v1
Q01;
[a b]=Q02
Q03;
Q04;
word= 'Question1: \nfigure(Q01) \n';
fprintf(word);
word= 'Question2: \nfigure(Q021) \nfigure(Q022) \n';
fprintf(word);
word= 'sigma is %0.3e. total signal power is %0.3e. \n';
fprintf(word,a,b);
word= 'Question3: \nfigure(Q03) \n';
fprintf(word);
word= 'Question4: \nfigure(Q041) \nfigure(Q042) \nfigure(Q043) \nfigure(Q044) \n';
fprintf(word);

end
function Q01% i produce the noise data and build my own hist() in this part and compare the hist histogram of acture data with PDF
mu=4; %value of mu
sigma=2; %value of sigma
data = randn(1000000,1) .* sigma + mu; %produce gussian data with mu and sigma
[y1,x]=histR(data,50);%Do histogram with data to show the distribution. y1 is count. x is x value range
y2 = exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);%compute the pdf of a gussian RV by function
l=x(2)-x(1);% get the length of each distribution range
y11=y1*l; % to the the probility of data in bar. i compute the area of each bar. 
figure('name','Q01');
bar(x,y1/sum(y11));hold on% show the pdf bar chart of histogram
plot(x,y2,'r');hold off% show the pdf 
title('PDF chart of Question 1')
xlabel('Value') ;
ylabel('Probility') ;
end
function [sigma p]=Q02%i build the a blackman window filter with fir1(). 
%i use the blackman window in creat a filter. it can't reach 50db at 12khz but it is best i can find
%then i filter the nose and produe the power density of  filtered noise
%about the 500 block average i should shift the data first.
%but i have hard time in programing my own fft shift.it keep throwing ot "out of memory"(my poor laptop<faceplam>)
%so i just do the average of half power density which cause some trouble in
%frequency<0;
k=1.380000000000000e-23;
f=3.1623;
t=300;
n0=k*t*f;%we have power density of noise here
sigma=sqrt(n0/2);%sigma is the Square root of power density

fs=50000;%sampling frequency
wp=10000/fs*2*pi;% pass band in angle
ws=12000/fs*2*pi;% stop band in angle
wt = ws - wp;%trans band of window
N = ceil(12*pi/wt);% level of window
wc = (ws+wp)/2;% cut off band in angle
h = 100*fir1(N,wc/pi,blackman(N+1));% i use the blackman window in creat a filter. it can't reach 50db at 12khz but it is best i can find
N1=fs/10;%window size of fft with resolution of 10hz
hdft = fft(h,N1);% fouier transform
hdft=ifft([hdft(length(hdft)/2+1:length(hdft)) hdft(1:length(hdft)/2)])
freq = fs/length(hdft):fs/length(hdft):fs% the X axis of filter spectrum

figure('name','Q021');
plot(freq,10*log10(abs(hdft))); 
title('Filter window')
xlabel('Frequency/hz') ;
ylabel('/dB') ;
N1=fs/1;
v = randn(N1,1) .* sigma;% produce the noise
y=filter(h,1,v);% filter the noise
ydft = fft(y);% fouier transform to get 
psdy = (1/(fs*N1)) * abs(ydft).^2;% do the power density
freq = fs/length(y):fs/length(y):fs% the X axis of filter spectrum

figure('name','Q022');
plot(freq,10*log10(psdy))
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)');

p=sum(abs(y).^2);

end
function Q03% i build 2 modle with the equation
k=1.380650e-23; % Boltzmann's constant
t= 300; %Temperature
c=2.99792458 * 10^8;%light speed
Band= 10^7;%Band width
Pt=2;%Transmit power
fc=10^9;%Farrier Frequency
gt=10;%Tansmiter gain
r=c/fc;%Wave lenghth
gr=1;%Receiver gain
F=10^0.5;%Noise figure
d=0;
Pn=F*k*t*Band;
x=1:1:100;
y1=zeros(100,1);
y2=zeros(100,1);
for i=1:1:100
    d=i*10^3;% Get distant form 1km-100km
    pl=((4*pi*d)/r)^2;% get path loss of free space model
    pl1=((4*pi*1000)/r)^2;% get path loss of free space model at 1km
    pl2=pl1*(d/1000)^3.5;% get path loss of exponent model
    Pr=Pt*gt*gr/pl;%Received power of free space model
    Pr2=Pt*gt*gr/pl2;%Received power of exponent model
    y1(i)=10*log10(Pr/Pn);%Gget the SNR into db
    y2(i)=10*log10(Pr2/Pn);%Get the SNR into db
end
figure('name','Q03');
semilogx(x,y1);hold on
semilogx(x,y2)
title('SNR chard of 2 propogation model')
xlabel('Distant(Km)')
ylabel('SNR (dB)');
hold off
end

function Q04%in the qeustion, i compute the r1 r2 r2' first and produce the pathloss with the equation.
%the i plot the snr and path loss form 1-3000m
k=1.380650e-23; % Boltzmann's constant
t= 300; %temperature
c=2.99792458 * 10^8;%light speed
Band= 10^7;%Band width
Pt=2;%Transmit power
fc=10^9;%Farrier Frequency
gt=10;%Tansmiter gain
wl=c/fc;%Wave lenghth
gr=1;%Receiver gain
F=10^0.5;%Noise figure
d=0;
Pn=F*k*t*Band;%Noise power
li=sqrt(-1);%j
h1=30;%Hight of transmitter
h2=10;%Hight of receiver
r1=0;%length from transmitter to receiver
l1=0;%length from bottom of transmitter to reflect point
l2=0;%length from bottom of receiver to reflect point
r21=0;%length from transmitter to reflect point
r22=0;%length from receiver to reflect point
S=0%Distant
d=1:1:3000;%Distant value
loss=zeros(1,3000);%path loss
snr=zeros(1,3000);%SNR
pr=zeros(1,3000);% Power recived 
for i=1:1:3000
    r1=sqrt(i^2+(h1-h2)^2);%get r1
    l1=i/(1+h1/h2);%get l1
    l2=i/(1+h2/h1);%get l2
    r21=sqrt(l1^2+h1^2);%get r21
    r22=sqrt(l2^2+h2^2);%get r22
    S=2*pi*(r21+r22-r1)/wl;%delta
    loss(i)=1/((wl/(4*pi))^2*(abs(1/r1+(-1/(r21+r22))*exp(-li*S)))^2);%compute the path loss at distant i with reflection coefficient of -1
    Pr(i)=Pt*gt*gr/loss(i);%compute the Power recived at distant i
    snr(i)=10*log10(Pr(i)/Pn);%compute the SNR at distant i
end
figure('name','Q041');
semilogx(d,10*log10(loss));
title('Path loss of 2ray model with reflection coefficient of -1')
xlabel('Distant(m)')
ylabel('Path loss (dB)');
figure('name','Q042');
semilogx(d,snr);
title('SNR of 2ray model with reflection coefficient of -1')
xlabel('Distant(m)')
ylabel('SNR (dB)');
for i=1:1:3000
    r1=sqrt(i^2+(h1-h2)^2);%get r1
    l1=i/(1+h1/h2);%get l1
    l2=i/(1+h2/h1);%get l2
    r21=sqrt(l1^2+h1^2);%get r21
    r22=sqrt(l2^2+h2^2);%get r22
    S=2*pi*(r21+r22-r1)/wl;%delta
    loss(i)=1/((wl/(4*pi))^2*(abs(1/r1+(-0.3/(r21+r22))*exp(-li*S)))^2);%compute the path loss at distant i with reflection coefficient of -0.3
    Pr(i)=Pt*gt*gr/loss(i);%compute the Power recived at distant i
    snr(i)=10*log10(Pr(i)/Pn);%compute the SNR at distant i
end
figure('name','Q043');
semilogx(d,10*log10(loss));
title('Path loss of 2ray model with reflection coefficient of -0.3')
xlabel('Distant(m)')
ylabel('Path loss (dB)');
figure('name','Q044');
semilogx(d,snr);
title('SNR of 2ray model with reflection coefficient of -0.3')
xlabel('Distant(m)')
ylabel('SNR (dB)');
end

function [no,xo] = histR(a,b) %replace command hist to get histogram. a is data. b is the total block number.
maxa=0;
mina=0;
maxa=max(a);%get max value
mina=min(a);%get min value
l=0;
l=(maxa-mina)/b;%get the length of each block
xo= zeros(1,b);
for i= 1:1:b
    xo(i)=mina+l*(i-2);%produce the map of block
end
no= zeros(1,b);
t=0;
for j= 1:1:length(a)
    t=round((a(j)-mina)/l+1)+1;%run the data add count to each block
    if t>b
        t=b;
    end
     if t==0
        t=1;
    end
    no(t)=no(t)+1;
end
end

function D = convR(A,B)
N=length(A)+length(B)-1
D=(ifft(fft(A,2*N-1).*fft(B,2*N-1)));
D=D(round((N-1)/2)+1:round((N-1)/2)+N) 
end

