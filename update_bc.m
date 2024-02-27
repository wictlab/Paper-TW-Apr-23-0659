function [Bc] = update_bc(Po,Bo,Pc,a,b,lamda_g,lamda1,lamda2,lamda3,lamda4,c,balloct,guiyihuaxishu_col,Ug,fenzi_col,P,B,G,U)
cvx_solver Mosek 
cvx_begin 
variable Bc(G,1)
minimize (-a*sum(balloct*Bc.*(log(1+fenzi_col.*Pc./(lamda_g+guiyihuaxishu_col.*(P-sum(Po))./(Ug)))/log(2)))+ ...
     (c/2).*(sum(Bc)-b*B )^2 + lamda3*(sum(Bc)-b*B))
subject to
Bc>=b*B/G/1.2;
%Bc>=0
sum(Bc)-b*B <=0;
cvx_end