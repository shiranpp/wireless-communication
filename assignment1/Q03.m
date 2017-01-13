function Q03
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
d=0;%Distant
Pn=F*k*t*Band;%Noise power
x=1:1:100;%Distant value
y1=zeros(100,1);%SNR of free space model
y2=zeros(100,1);%SNR of exponent model
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
figure('name','Q031');
semilogx(x,y1);hold on
semilogx(x,y2)
title('SNR chard of 2 propogation model')
xlabel('Distant(Km)')
ylabel('SNR (dB)');
hold off
end