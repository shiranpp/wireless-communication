function Q01
clear;
sigma=1; %value of sigma
ji=sqrt(-1);%j
g= randn(1000000,1) .* sigma + ji * randn(1000000,1) .* sigma ;%gain with complex gussian noise
Z=sqrt(real(g).^2+imag(g).^2);
[y1,x1]= histR(Z,200);
l=max(Z)-min(Z);
y1=y1/1000000/diff(x1(1:2));
y2=zeros(1,length(x1));
for i = 1:1:length(x1)
    y2(i)=(x1(i)/(sigma^2))*exp(-(x1(i)^2)/(2*(sigma^2)));
end
plot(x1,y1,x1,y2,'.');
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