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

function Q=Qf(x)
Q=0.5-0.5*erf(x/sqrt(2))
end