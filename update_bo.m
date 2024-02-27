function [Bo] = update_bo(x,Po,Bc,Pc,a,b,lamda_g,lamda1,lamda2,lamda3,lamda4,c,balloct,guiyihuaxishu_col,Ug,fenzi_col,P,B,G,U)
%x=100000;
cvx_solver Mosek
cvx_begin quiet
variables Bo(G,1)
minimize (-a*sum((balloct*Bo./2).*log(1+(x*Po.*fenzi_col)./(lamda_g+guiyihuaxishu_col.*(P-sum(Pc))./(Ug))))/log(2)+ ...
     (c/2).*((sum(Bo)-(1-b)*B )^2) + lamda4*(sum(Bo)-(1-b)*B))
subject to
Bo>=(1-b)*B/G/1.2
%Bo>=0
sum(Bo)-(1-b)*B <=0
cvx_end