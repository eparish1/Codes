clear all
close all
syms um4 um3 um2 um1 u0 up1 up2 up3 up4 nu
k = linspace(-4,4,9)'; 
N = 9;
ka = linspace(-4,4,N);
uhat = [0,0,um2,um1,u0,up1,up2,0,0].';
uhatF = [um2;um1;0;up1;up2];
utilde = [um4,um3,0,0,0,0,0,up3,up4].';
u = uhat + utilde;
p = zeros(N,1);
p(3:7) = k(3:7);
Rk = -1j.*k./2.*symconv(u,u,ka) - k.^2.*u.*nu;
Lu = Rk;
PLu = simplify( subs(Lu,[um4,um3,up3,up4],[0,0,0,0]) ); 
QLu = Lu - PLu; 
Jac = sym(zeros(N,N));
for i=1:9
  Jac(:,i) = diff(QLu,u(i));
end
Jacp = sym(zeros(N,N));
for i=1:9
  Jacp(:,i) = diff(PLu,u(i));
end
LQLu = Jac*Rk;
PLQLu = subs(LQLu,[um4,um3,up3,up4],[0,0,0,0]);
QLQLu =  LQLu - PLQLu;
LPLu = Jacp*Rk;
PLPLu = subs(LPLu,[um4,um3,up3,up4],[0,0,0,0]);

%% Second Order Models
Jac2 = sym(zeros(N,N));
for i=1:9
  Jac2(:,i) = diff(QLQLu,u(i));
end
Jac2p = sym(zeros(N,N));
for i=1:9
  Jac2p(:,i) = diff(PLQLu,u(i));
end
LQLQLu = Jac2*Rk;
PLQLQLu = simplify( subs(LQLQLu,[um4,um3,up3,up4],[0,0,0,0]) ); 
QLQLQLu = LQLQLu - PLQLQLu;

LPLQLu = Jac2p*Rk;
PLPLQLu = simplify( subs(LPLQLu,[um4,um3,up3,up4],[0,0,0,0]) ); 

%% Third Order Models
Jac3 = sym(zeros(N,N));
for i=1:9
  Jac3(:,i) = diff(QLQLQLu,u(i));
end

LQLQLQLu = Jac3*Rk;
PLQLQLQLu = simplify( subs(LQLQLQLu,[um4,um3,up3,up4],[0,0,0,0]) ); 
% Testing
val = -1j.*k./2.*symconv(utilde,utilde);
JacTest = sym(zeros(N,N));
for i=1:9
  JacTest(:,i) = diff(val,u(i));
end
val2 = JacTest*Rk;
Pv = simplify( subs(val2,[um4,um3,up3,up4],[0,0,0,0]) ); 
Qv = val2 - Pv;
JacTest2 = sym(zeros(N,N));
for i=1:9
  JacTest2(:,i) = diff(Qv,u(i));
end
val3 = JacTest2*Rk;
PLQv = simplify( subs(val3,[um4,um3,up3,up4],[0,0,0,0]) ); 
QLQv = val3 - PLQv;
JacTest3 = sym(zeros(N,N));
for i=1:9
  JacTest3(:,i) = diff(QLQv,u(i));
end
val4 = JacTest3*Rk;
PLQLQv = simplify( subs(val4,[um4,um3,up3,up4],[0,0,0,0]) ); 


%% Convolutions
PLu_2 = -1j.*k./2.*symconv(uhat,uhat) - nu*k.^2.*uhat;
PLvisc = -nu.*k.^2.*PLu_2;
PLu_2q = PLu_2; PLu_2q(3:7) = 0.;
PLu_2p = PLu_2 - PLu_2q;

PLPLu_2  = 2.*-1j.*k./2.*symconv(uhat,PLu_2p) - k.^2.*nu.*PLu_2p;
PLPLu_2q = PLPLu_2; PLPLu_2q(3:7) = 0.;
PLPLu_2p = PLPLu_2 - PLPLu_2q;
PLQLu_2  = 2.*-1j.*k./2.*symconv(uhat,PLu_2q) - k.^2.*nu.*PLu_2q;
PLQLu_2q = PLQLu_2; PLQLu_2q(3:7) = 0.;
PLQLu_2p = PLQLu_2 - PLQLu_2q;

PLPLQLu_2 = 2.*-1j.*k./2.*symconv(PLPLu_2q,uhat) + 2.*-1j.*k./2.*symconv(PLu_2q,PLu_2p) - ...
            nu.*k.^2.*PLPLu_2q;
PLQLQLu_2 = 2.*-1j.*k./2.*symconv(PLQLu_2q,uhat) + 2.*-1j.*k./2.*symconv(PLu_2q,PLu_2p) + ...
            2.*-1j.*k./2.*symconv(PLu_2q,PLu_2q) - nu.*k.^2.*PLQLu_2q;
PLQLQLu_2q = PLQLQLu_2; PLQLQLu_2q(3:7) = 0.;
PLQLQLu_2p = PLQLQLu_2 - PLQLQLu_2q;


PLQLQLQLu_2 = 2.*-1j.*k./2.*symconv(uhat,PLQLQLu_2q) + 4.*-1j.*k./2.*symconv(PLu_2p,PLQLu_2q) + ...
              4.*-1j.*k./2.*symconv(PLQLu_2p,PLu_2q) + 2.*-1j.*k./2.*symconv(PLPLu_2p,PLu_2q) + ...
              6.*-1j.*k./2.*symconv(PLQLu_2q,PLu_2q )  + 2.*-1j.*k./2.*symconv(PLPLu_2q,PLu_2q) - ...
              nu.*k.^2.*PLQLQLu_2q;

%PLPLu_2 = -1j.*k./2.*symconv(2.*uhat,-1j.*p./2.*symconv(uhat,uhat)-nu*p.^2.*uhat.') - p.^2.*nu.'.*Rkbar;
%PLu_p = RkFp;
%PLQLu_2 = 2*-1j.*k./2.*symconv(PLu_p,uhat)%; - nu*q.^2.*RkF;
%QLu_2 = -1j.*k.*symconv(uhat,utilde) - 1j.*k./2.*symconv(utilde,utilde) - nu*k.^2.*utilde;
%LQLu_2 = -1j.*k.*symconv(u,Rk) - nu*k.^2.*Rk - (-1j.*k.*symconv(uhat,Rkp) - nu*k.^2.*Rkp);
%QLQLu_2 = -1j.*k.*symconv(utilde,Rk) + 1j.*k.*symconv(utilde,Rkq) - nu.*q.^2.*Rk;
