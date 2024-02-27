function [Po] = update_po(x,Bo,Bc,Pc,a,b,lamda_g,lamda1,lamda2,lamda3,lamda4,c,balloct,guiyihuaxishu_col,Ug,fenzi_col,P,B,G,U)
%x=100000;
cvx_solver Mosek
cvx_begin quiet
variables Po(U,1)
minimize (-a*sum(balloct*Bo./2.*log(1+(x*Po.*fenzi_col)./(lamda_g+guiyihuaxishu_col.*(P-sum(Pc))./Ug)))/log(2)+ ...
     lamda2*(sum(Po)-(1-b)*P) + c/2*(sum(Po)-(1-b)*P)^2)
subject to
Po>=(1-b)*P/U/1.2
%Po>=0
sum(Po)-(1-b)*P<=0
cvx_end