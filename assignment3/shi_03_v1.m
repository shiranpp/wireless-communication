function shi_03_v1
Q01;
Q02;
Q03;
Q04;
word= 'Question1: \nfigure(Q01) \n';
fprintf(word);
word= 'Question2: \nfigure(Q02) \n';
fprintf(word);
word= 'Question3: \nfigure(Q03) \n';
fprintf(word);
word= 'Question4: \nfigure(Q04) \n';
fprintf(word);
end

%Queation 1%
function Q01
clear;
sigma=1; %value of sigma
ji=sqrt(-1);%j
g= randn(1000000,1) .* sigma + ji * randn(1000000,1) .* sigma ;%gain with complex gussian noise
Z=sqrt(real(g).^2+imag(g).^2); %produce z value withe sqrt(real image part of noise) rayleigh
[y1,x1]= histR(Z,200);%produce the histgrom of rayleigh distribution
y1=y1/1000000/diff(x1(1:2));%experimental pdf-s of rayleigh distribution
y2=zeros(1,length(x1));
for i = 1:1:length(x1)
    y2(i)=(x1(i)/(sigma^2))*exp(-(x1(i)^2)/(2*(sigma^2)));% theoretical pdf-s of rayleigh distribution
end
figure('name','Q01');
plot(x1,y1,x1,y2,'.');
title('Rayleigh Distribution Fading gain pdf-s')
xlabel('Value') ;
ylabel('Probility') ;
end

%Queation 2%
function Q02
clear;
sigma = 1; %value of sigma
ji = sqrt(-1);%j
v = 4;%constant
y = randn(1000000,1) .* sigma + ji * randn(1000000,1) .* sigma ;%complex gussian noise
g = v + y;%gain 
Z = abs(g);%produce Z
[y1,x1]= histR(Z,100);%produce the histgrom of raician distribution
y1=y1/1000000/diff(x1(1:2));%experimental pdf-s of raician distribution
y2=zeros(1,length(x1));
for i = 1:1:length(x1)
    y2(i)=(x1(i)/(sigma^2))*exp(-((x1(i)^2)+(v^2))/(2*(sigma^2)))*besseli(0,x1(i).*v./(sigma^2));%theoretical pdf-s of raician distribution
end
figure('name','Q02');
plot(x1,y1,x1,y2,'.');
title('Raician Distribution Fading gain pdf-s')
xlabel('Value') ;
ylabel('Probility') ;
end

%Question 3%
%The experimental value is 8% less then theoretical value%
%i try to slove it but i can't find out where is worng.%
function Q03
clear;
sigma=10;
fd=50;
fs=1000;
ji=sqrt(-1);
fx=-fs/2:1:fs/2;
S=zeros(1,length(fx));
n = randn(1,1000000) .* sigma+ji*randn(1,1000000) .* sigma;%complex white Noise with sigama 10 mu 0
for i=-fs/2:1:fs/2
    if abs(i)<fd
        S(i+fs/2+1)=1/sqrt(1-(i/fd)^2);%Produce the PSD of Doppler function
    else
        S(i+fs/2+1)=0;
    end
end
S1=circshift(S,ceil(size(S)/2));%shift the PSD so it have DC at 0
h=ifft(sqrt(S1))%Doppler filter
h=circshift(h,floor(size(h)/2));%shift impluse of Doppler filter so it can have window at middle
g=filter(h,1,n);% produce the gain by fikter the complex AWGN withe doppler filter
T=length(n)/fs;%total time in second
a=sum(abs(g).^2)/length(g);%power in gain
z=abs(g);%produce Z
x=0:0.1:10;
LR1=zeros(1,101);
LR2=zeros(1,101);
for i=1:1:101
    R=(i-1)/10;
    C=R/sqrt(a);
    LR1(i)=sqrt(2*pi)*fd*C*exp(-C^2);
    for j=fs:1:length(z)-fs
        if z(j)<=R
            if z(j+1)>R
                LR2(i)=LR2(i)+1; % count+1 when the Z cross R
            end
        end
 %try to slove 8% error with wider R but this isn' work well, when the%
 %magitude is close to average%
    end
    LR2(i)=LR2(i)/(T-2);% produce the LCR by count/Time
end
figure('name','Q03');
plot(x,LR1,x,LR2)
title('LCR of fading gain from assignment 2')
ylabel('LCR(Hz)') ;
xlabel('Level R (magitude)') ;
end

%Question 4%
%There is still 8% difference between experimental and theoretical value%%
%I thinke it cause by the LCR error from last queation.%
function Q04
clear;
sigma=10;
fd=30;
fs=1000;
ji=sqrt(-1);
fx=-fs/2:1:fs/2;
S=zeros(1,length(fx));
n = randn(1,1000000) .* sigma+ji*randn(1,1000000) .* sigma;%complex white Noise with sigama 10 mu 0
for i=-fs/2:1:fs/2
    if abs(i)<fd
        S(i+fs/2+1)=1/sqrt(1-(i/fd)^2);%Produce the PSD of Doppler function
    else
        S(i+fs/2+1)=0;
    end
end
S1=circshift(S,ceil(size(S)/2));%shift the PSD so it have DC at 0
h=ifft(sqrt(S1))%Doppler filter
h=circshift(h,floor(size(h)/2));%shift impluse of Doppler filter so it can have window at middle
T=length(n)/fs;%time table in seconds
g=filter(h,1,n);% produce the gain by fikter the complex AWGN withe doppler filter
a=sum(abs(g).^2)/length(g);%power in gain
z=abs(g);%produce Z
x=0:0.1:10;
LR1=zeros(1,101);
t1=zeros(1,101);
LR2=zeros(1,101);
t2=zeros(1,101);
for i=1:1:101
    R=(i-1)/10;
    C=R/sqrt(a);
    LR1(i)=sqrt(2*pi)*fd*C*exp(-C^2);
    t1(i)=(exp(C^2)-1)/(sqrt(2*pi)*fd*C);
    for j=fs:1:length(z)-fs
        if z(j)<=R
            if z(j-1)>R
                LR2(i)=LR2(i)+1; % count+1 when the Z cross R
            end
            t2(i)=t2(i)+1/fs;% time counter + 1ms when under the R
        end
        if (z(j)<R+0.2)&(z(j)>R)
             if (z(j+1)>R+0.2)&(z(j-1)>R+0.2)
                LR2(i)=LR2(i)+1;
           end
        end
    if (z(j)>R-0.2)&(z(j)<R)
       if (z(j+1)<R-0.2)&(z(j-1)<R-0.2)
            LR2(i)=LR2(i)+1;
        end
     end

    end
    t2(i)=t2(i)/LR2(i);%produce the average fading time by total fading time/Level crossed
end
figure('name','Q04');
plot(x,t1,x,t2)
title('average fading time of fading gain from assignment 2')
ylabel('average fading time(s)') ;
xlabel('Level R (magitude)') ;
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