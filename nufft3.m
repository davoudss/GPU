%Irregular in X to regular in K


%based on Greengard with x in [a b] and ktilde = (-pi*M/2:pi*M/2-1) to compute
%F(k)=1/N*(sum(f_j*exp(-i*ktilde*xj))) with Matlab FFT
%with spreading Gaussian


clc, clear all, close all

format longe
% I = input('Interval([a b]): ');
% N = input('# of sample points: ');
% M = input('# of oversamples points: ');
% accuracy = input('accuracy (s/d) ','s');
% P = input('# of points to truncate Gaussian (6/12): ');
% a = I(1);
% b = I(2);
tic

a = 2;
b = 7;
L = (b-a)/2;
M = 16;
N = 8;
Mr = 2*M;
h = 2*pi/Mr;
P = 12;
accuracy='d';

l=(-1:1)';
% if(accuracy=='s')
%     Tau = 6/M^2;
%  else
%      Tau = 12/M^2;
%  end
    

Tau = 12/M^2;

k = (-M/2:M/2-1)';
ktilde = pi*k/L;


%initialize vectors
Fexact   = zeros(length(k),1);
ftau     = zeros(Mr,1);
Ftau     = zeros(length(k),1);

%initial values
vec = (1:N)';
yj=a + 2*L*cos(vec).^2;         %yj in [-a b]
xj= yj-(b+a)/2;                 %xj in [-L L]
Xj = pi*xj/L;                   %Xj in [-pi pi]

fj=-1+2*xj;                     %some values


%-------------------------------------------------
%find the exact sum
for kk=1:length(k)
    Fexact(kk) = sum(fj.*exp(-1i*ktilde(kk)*xj))*exp(-1i*ktilde(kk)*(b+a)/2);
end

%Find approximate sum
%--------------------------------------------------

for jj=1:N
    m1 = round(Xj(jj)/h);
    for m=m1-P:m1+P
        cnt = m;
        if (cnt<0)
            cnt = cnt + Mr;
        end
        if (cnt>=Mr)
            cnt = cnt - Mr;
        end
        gtau = sum(exp(-(2*pi*m/Mr-Xj(jj)-2*l*pi).^2/(4*Tau)));
        ftau(cnt+1) = ftau(cnt+1) + fj(jj)*gtau;
    end
end

%manually
% for kk=1:length(k)
%     Ftau(kk) = sum(ftau.*exp(-1i*2*pi*k(kk)*(0:Mr-1)'/Mr))/Mr.*exp(-1i*ktilde*(b+a)/2);    
% end

%with fft
Ftau = fft(ftau)/Mr;
fftshift(fft(ftau))
pause
Ftau = fftshift(Ftau);
Ftau = Ftau(M/2+1:M+M/2);


Fapprox = sqrt(pi/Tau).*exp(k.^2*Tau).*Ftau.*exp(-1i*ktilde*(b+a)/2);

%plot(real(Fexact),imag(Fexact),'.',real(Fapprox),imag(Fapprox),'ro')
%legend('Exact','Appproximation')

disp(['Error:', num2str(abs(norm(Fexact-Fapprox))/norm(Fexact))])

toc