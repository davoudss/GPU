%Irregular in X to regular in K

%based on Greengard with x in [a b] and ktilde = (-pi*M/2:pi*M/2-1) to compute
%F(k)=1/N*(sum(f_j*exp(-i*ktilde*xj))) with Matlab FFT
%with spreading Gaussian


function nufft3d

clc, clear all, close all
tic
a = 2;
b = 7;
L = (b-a)/2;

M = 32;
N = 4;
Mr= 2*M;
h = 2*pi/Mr;
P = 12;

Tau = 12/M^2;

u = (-M/2:M/2-1)';
v = (-M/2:M/2-1)';
w = (-M/2:M/2-1)';

utilde = pi*u/L;
vtilde = pi*v/L;
wtilde = pi*w/L;

%initialize vectors
Fexact   = zeros(length(u),length(v),length(w));
Fapprox   = zeros(length(u),length(v),length(w));
ftau     = zeros(Mr,Mr,Mr);
Ftau     = zeros(length(u),length(v),length(w));

%initial values
vecx = (1:N)';
vecy = (1:N)';
vecz = (1:N)';

x1=a + 2*L*cos(vecx).^2;         %yj in [-a b]
x2= x1-(b+a)/2;                 %xj in [-L L]
x = pi*x2/L;                   %Xj in [-pi pi]

y1=a + 2*L*cos(vecy).^2;         %yj in [-a b]
y2= y1-(b+a)/2;                 %xj in [-L L]
y = pi*y2/L;                   %Xj in [-pi pi]


z1=a + 2*L*cos(vecz).^2;         %yj in [-a b]
z2= z1-(b+a)/2;                 %xj in [-L L]
z = pi*z2/L;                   %Xj in [-pi pi]

[X Y Z]=meshgrid(x1,y1,z1);

%f=-1+2*X+3*Y-2*Z;                     %some values
for n2=1:N
  for n3=1:N
    for n1=1:N
      f(n1,n2,n3)  = -1+2*x2(n1)+3*y2(n2)-2*z2(n3);
    end
  end
end

	


norm(x)
f
pause

%-------------------------------------------------
%find the exact sum
for cu=1:length(u)
    for cv=1:length(v)
        for cw=1:length(w)
            for n1=1:N
                for n2=1:N
                    for n3=1:N
                    Fexact(cu,cv,cw) = Fexact(cu,cv,cw) + ...
                        f(n1,n2,n3)*exp(-1i*utilde(cu)*x2(n1)-1i*vtilde(cv)*y2(n2)-1i*wtilde(cw)*z2(n3)) ...
                        *exp(-1i*(utilde(cu)+vtilde(cv)+wtilde(cw))*(b+a)/2);
                    end
                end
            end
        end
    end
end

%Find approximate sum
%--------------------------------------------------
disp(['Exact: ',num2str(toc)])
tic


for n1=1:N
    cut1 = round(x(n1)/h);
    for n2=1:N
        cut2 = round(y(n2)/h);
        for n3=1:N
            cut3 = round(z(n3)/h);
            for m1=cut1-P:cut1+P
                c1 = m1;
                while (c1<0)
                    c1 = c1 + Mr;
                end
                while (c1>=Mr)
                    c1 = c1 - Mr;
                end
                for m2=cut2-P:cut2+P
                    c2 = m2;
                    while (c2<0)
                        c2 = c2 + Mr;
                    end
                    while (c2>=Mr)
                        c2 = c2 - Mr;
                    end
                    for m3=cut3-P:cut3+P
                        c3 = m3;
                        while (c3<0)
                            c3 = c3 + Mr;
                        end
                        while (c3>=Mr)
                            c3 = c3 - Mr;
                        end
                        g = gtau(x(n1),y(n2),z(n3),c1,c2,c3,Tau,Mr);
                        ftau(c1+1,c2+1,c3+1) = ftau(c1+1,c2+1,c3+1) + f(n1,n2,n3)*g;
                    end
                end
            end
        end
    end
end


%with fft
disp(['Gaussian: ',num2str(toc)])
tic
Ftau = fftn(ftau)/Mr/Mr/Mr;
Ftau = fftshift(fftshift(fftshift(Ftau,1),2),3);
Ftau = Ftau(M/2+1:M+M/2,M/2+1:M+M/2,M/2+1:M+M/2);
disp(['fft+shift: ',num2str(toc)])


for cu=1:numel(u)
    for cv=1:numel(v)
        for cw=1:numel(w)
            Fapprox(cu,cv,cw) = sqrt(pi/Tau)*sqrt(pi/Tau)*sqrt(pi/Tau)*exp((u(cu)^2+v(cv)^2+w(cw)^2)*Tau) ... 
                *Ftau(cu,cv,cw)*exp(-1i*(utilde(cu)+vtilde(cv)+wtilde(cw))*(b+a)/2);
        end
    end
end

disp(['relative error: ',num2str(abs(norm3(Fexact - Fapprox))/norm3(Fexact))])


function g=gtau(x,y,z,ct,cs,cq,Tau,Mr)
g=0;
for l1=-1:1
    for l2=-1:1
        for l3=-1:1
        g=g+exp(-((2*pi*ct/Mr-x-2*l1*pi)^2+(2*pi*cs/Mr-y-2*l2*pi)^2+(2*pi*cq/Mr-z-2*l3*pi)^2)/(4*Tau));
        end
    end
end


function nr = norm3(x)
nr = 0;
if(size(x,3)==1)
    error('not useful')
end
    
for ii=1:size(x,3)
    nr = nr + norm(x(:,:,ii))^2;
end

nr = sqrt(nr);
    
    
    
    
    



