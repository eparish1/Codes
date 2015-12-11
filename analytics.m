clear all
close all
syms um4 um3 um2 um1 u0 up1 up2 up3 up4 nu
k = linspace(-4,4,9)'; 
N = 9;
ka = linspace(-4,4,N);
uhat = [0,0,um2,um1,u0,up1,up2,0,0];
uhatF = [um2;um1;0;up1;up2];
utilde = [um4,um3,0,0,0,0,0,up3,up4];
u = uhat + utilde;
q = k;
q(3:7,1) = 0;
p = k - q;
uhatNum = [0,0,4,1,0,1,4,0,0];
Rk = -1j.*k./2.*symconv(u,u,ka) - k.^2.*u.'.*nu;
Lu = Rk;
PLu = simplify( subs(Lu,[um4,um3,up3,up4],[0,0,0,0]) ); 
QLu = Lu - PLu; 
Jac = sym(zeros(N,N));
for i=1:9
  Jac(:,i) = diff(QLu,u(i));
end
Rk = -1j.*k./2.*symconv(u,u,ka) - k.^2.*u.'.*nu;
RkF= Rk(3:7);
Rk3 = sym(zeros(9,1));
Rk3(3:7,1) = Rk(3:7,1);
shift = 5;
LQLu = Jac*Rk;
PLQLu = subs(LQLu,um4,0);
PLQLu = subs(PLQLu,um3,0);
PLQLu = subs(PLQLu,up3,0);
PLQLu = subs(PLQLu,up4,0);
PLQLunum = subs(PLQLu,[um2,um1,u0,up1,up2],uhatNum(-2+shift:2+shift))
%%%%%%%% Done with zeroth order t-model%%%%%%%%%%%
%% second order model
Jac2 = sym(zeros(N,N));
for i=1:9
  Jac2(:,i) = diff(PLQLu,u(i));
end
LPLQLu = Jac2*Rk;
PLPLQLu = simplify( subs(LPLQLu,[um4,um3,up3,up4],[0,0,0,0]) ); 
PLPLQLunum = subs(PLPLQLu,[um2,um1,u0,up1,up2],[uhatNum(-2+shift:2+shift)])

%%%%%%%% Now generate individual terms
Rkbar = subs(Rk,[um4,um3,up3,up4],[0,0,0,0]);
Rkbar2 = subs(Rk3,[um4,um3,up3,up4],[0,0,0,0]);
Lu = Rk;
PLu = subs(Lu,[um4,um3,up3,up4],[0,0,0,0]);
Jac4 = sym(zeros(N,N));
for i=1:9
  Jac4(:,i) = diff(PLu,u(i));
end
LPLu = Jac4*Rk;
PLPLu = subs(LPLu,[um4,um3,up3,up4],[0,0,0,0]);

JacTemp = sym(zeros(N,N));
for i=1:9
  JacTemp(:,i) = diff(utilde,u(i));
end
Lu_q = JacTemp*Rk;
PLu_q = subs(Lu_q,[um4,um3,up3,up4],[0,0,0,0]);
for i=1:9
  JacTemp(:,i) = diff(PLu_q,u(i));
end
LPLu_q = JacTemp*Rk;
PLPLu_q = subs(LPLu_q,[um4,um3,up3,up4],[0,0,0,0]);

JacTemp = sym(zeros(N,N));
for i=1:9
  JacTemp(:,i) = diff(uhat,u(i));
end
Lu_p = JacTemp*Rk;
PLu_p = subs(Lu_p,[um4,um3,up3,up4],[0,0,0,0]);
for i=1:9
  JacTemp(:,i) = diff(PLu_p,u(i));
end

LPLu_p = JacTemp*Rk;
PLPLu_p = subs(LPLu_p,[um4,um3,up3,up4],[0,0,0,0]);

check = 2.*(-1j.*k./2.*symconv(uhat,PLPLu_q) - 1j.*k./2.*symconv(PLu_p,PLu_q));
t1 = symconv(uhat,PLPLu_q);
t1Check = -1j.*q./2.*symconv(uhat,symconv(2.*uhat , -1j.*q./2.* symconv(uhat,uhat)));



%%% Final Test with no matrix multiplication
PLu_2  = -1j.*k./2.*symconv(uhat,uhat) - nu*k.^2.*uhat.';
PLu_2q = PLu_2; PLu_2q(3:7) = 0.;
PLu_2p = PLu_2 - PLu_2q;
PLPLu_2 = -1j.*k./2.*symconv(2.*uhat,-1j.*p./2.*symconv(uhat,uhat)-nu*p.^2.*uhat.') - p.^2.*nu.'.*Rkbar;
PLPLu_2q = PLPLu_2; PLPLu_2q(3:7) = 0.;
PLPLu_2p = PLPLu_2 - PLPLu_2q;
tmod = 2.*( -1j.*k./2.*symconv(uhat,PLu_2q) + -0.5.*q.^2*nu.*PLu_2q);
tmodNum = subs(tmod,[um2,um1,u0,up1,up2],[uhatNum(-2+shift:2+shift)])
t2mod = 2.*(-1j.*k./2.*symconv(uhat,PLPLu_2q) - 0.5.*q.^2.*nu.*PLPLu_2q  ...
           +  -1j.*k./2.*symconv(PLu_2p,PLu_2q) - 0.5.*p.^2.*nu.*PLu_2q);
t2modNum = subs(t2mod,[um2,um1,u0,up1,up2],[uhatNum(-2+shift:2+shift)])

