%clear all
%close all

N = 9;
uhat =   [0,0,4,1,0,1,4,0,0];
utilde = [4,3,0,0,0,0,0,3,4];
k = linspace(-4,4,9);
c = conv(utilde,uhat,'same');
shift = 5;
du = 1e-4;
dc = zeros(N,N);
for j = -4:4
  i = j + shift;
  if (abs(j) > 2)
    utilde(i) = utilde(i) + du;
    cp = conv(utilde,uhat,'same');
    dc(:,i) =( cp(:) - c(:) ) / du;
    utilde(i) = utilde(i) - du;
  else
    uhat(i) = uhat(i) + du;
    cp = conv(utilde,uhat,'same');
    dc(:,i) =( cp(:) - c(:) ) / du;
    uhat(i) = uhat(i) - du;
  end
end

Rk = -1j.*k/2.*conv(uhat,uhat,'same');
L = ones(N,1);
for i=-4:4
  M = i+shift;
  L = L+Rk(M)*1j.*k(M).*dc(:,M);
end
M = [-4,-3,0,0,0,0,0,3,4];
testConv1 = conv(uhat,uhat,'same');
testConv2 = -1j.*M./2.*testConv1;
testConv3 = -1j.*k.*conv(uhat,testConv2,'same');
testConv3'
%{
clear all
N = 9;
utilde = [14,123,0,0,0,0,0,35,41];
%k = uhat + utilde;
c = conv(utilde,utilde,'same');
shift = 5;
du = 1e2;
dc = zeros(N,N);
for j = -4:4
  i = j + shift;
  if (abs(j) > 2)
    utilde(i) = utilde(i) + du;
    cp = conv(utilde,utilde,'same');
    dc(:,i) =( cp(:,i) - c(:,i) ) / du;
    utilde(i) = utilde(i) - du;
  end
end
%}
