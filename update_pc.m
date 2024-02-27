function [Pc] = update_pc(Bo,Bc,Po,a,b,lamda_g,lamda1,lamda2,lamda3,lamda4,c,balloct,guiyihuaxishu_col,Ug,fenzi_col,P,B,G,U)
cvx_solver Mosek
cvx_begin quiet
variable Pc(U,1)
minimize (-a*sum(balloct*Bc.*(log(1+fenzi_col.*Pc./(lamda_g+guiyihuaxishu_col.*(P-sum(Po))./(Ug)))/log(2)))+ ...
     lamda1*(sum(Pc)-b*P) + c/2*(sum(Pc)-b*P)^2)
subject to
Pc>=b*P/U/1.2
%Po>=0
sum( Pc)-b*P<=0
cvx_end