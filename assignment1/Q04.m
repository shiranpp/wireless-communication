function Q04
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