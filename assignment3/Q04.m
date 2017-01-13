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
  %      if (z(j)<R+0.02)&(z(j)>R)
   %         if (z(j+1)>R+0.01)&(z(j-1)>R+0.01)
   %             LR2(i)=LR2(i)+1;
   %         end
   %     end
   %     if (z(j)>R-0.02)&(z(j)<R)
   %         if (z(j+1)<R-0.02)&(z(j-1)<R-0.02)
   %             LR2(i)=LR2(i)+1;
   %         end
   %     end

    end
    t2(i)=t2(i)/LR2(i);%produce the average fading time by total fading time/Level crossed
end
figure('name','Q04');
plot(x,t1,x,t2)
title('average fading time of fading gain from assignment 2')
ylabel('average fading time(s)') ;
xlabel('Level R (magitude)') ;
end