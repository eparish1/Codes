clear all
close all
digits(64)
nx = 50;
dx = 2*pi/nx;
x = linspace(0,2*pi-dx,nx);
k = linspace(-nx/2,nx/2-1,nx);
M = linspace(-nx,nx-1,2*nx);

for i = 1:2*nx
  if (M(i) <= nx/2-1 && M(i) >= -nx/2)
    M(i) = 0;
  end
end
u = sin(x) + cos(20*x);
uhat = myfft(u,nx);

uhat2 = zeros(1,2*nx);
uhat2(nx/2+1:nx/2+nx) = uhat;
v = conv(uhat2,uhat2,'same');
%v = fftshift(cconv(uhat2,uhat2,2*nx));
v2 = 1j*M./2.*v

fin = conv(uhat,v2,'same')

w = conv(uhat2,uhat2,'same');
w2 = 1j.*M./2.*w;
fin2 = conv(w2,uhat,'same');


%% Looping
Gshift = nx + 1;
Fshift = nx/2 + 1;
val = zeros(nx,1);
Farray = linspace(-nx/2,nx/2-1,nx);
Garray = linspace(-nx,-nx/2-1,nx/2);
Garray = [Garray,linspace(nx/2,nx-1,nx/2)];

for k = Farray
    for p = Farray
    for q = Garray
    if (p + q == k)
      convsum = 0.
      for r = Farray
      for s = Farray
      if (r + s == q)
              convsum =  convsum +  uhat(r + Fshift).*uhat(s + Fshift);
      end
      end
      end
      val(k + Fshift) = val(k + Fshift) + uhat(p + Fshift).*1j*q./2.*convsum;
    end
    end
    end
end
val3 = zeros(nx,1);
for k = Farray
    convsum = zeros(2*nx,1);
    for p = Garray
    for q = Farray
    if (p + q == k)
      for r = Farray
      for s = Farray
      if (r + s == p)
              convsum(p + Gshift) = convsum(p + Gshift) +  uhat(r + Fshift).*uhat(s + Fshift);
      end
      end
      end
      val3(k + Fshift) = val3(k + Fshift) + uhat(q + Fshift).*1j*p./2.*convsum(p + Gshift);
    end
    end
    end
end

   
val2 = conv(uhat,1j*M/2.*conv(uhat2,uhat2,'same'),'same');
