function Q03
clear;
sigma=1;
fd=50;
fs=1000;
ji=sqrt(-1);
fx=-fs/2:1:fs/2;
S=zeros(1,length(fx));
n = randn(1,1000000) .* sigma+ji*randn(1,1000000) .* sigma;%cmplex white Noise with sigama 1 mu 0
for i=-fs/2:1:fs/2
    if abs(i)<fd
        S(i+fs/2+1)=(sigma^2/4/pi/fd)/sqrt(1-(i/fd)^2);%Produce the PSD of Doppler function
    else
        S(i+fs/2+1)=0;
    end
end
S1=circshift(S,ceil(size(S)/2));
h=ifft(sqrt(S1))%Doppler filter
h=circshift(h,floor(size(h)/2));
g=filter(h,1,n);
T=length(n)/fs;%total time in second
a=sum(abs(n).^2)/length(n);%power in gain
z=abs(n);
x=0:0.01:1;
LR1=zeros(1,101);
LR2=zeros(1,101);
for i=1:1:101
    R=(i-1)/100;
    C=R/sqrt(a);
    LR1(i)=sqrt(2*pi)*fd*C*exp(-C^2);
    for j=fs:1:length(z)-fs
        if z(j)<=R
            if z(j+1)>R
                LR2(i)=LR2(i)+1;
            end
        end
        if (z(j)<R+0.02)&(z(j)>R)
            if (z(j+1)>R+0.02)&(z(j-1)>R+0.02)
                LR2(i)=LR2(i)+1;
            end
        end
        if (z(j)>R-0.01)&(z(j)<R)
            if (z(j+1)<R-0.02)&(z(j-1)<R-0.02)
                LR2(i)=LR2(i)+1;
            end
        end


    end
    LR2(i)=LR2(i)/(T-2);
end
plot(x,LR1,x,LR2)
end