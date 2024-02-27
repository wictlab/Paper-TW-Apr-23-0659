function [Bo,Po] = update_bo_po(Bc,Pc,a,b,lamda_g,n,lamda1,lamda2,c,balloct,guiyihuaxishu_col,Ug,fenzi_col,rmin,P,B,G,U)
cvx_solver sedumi
cvx_begin quiet
variables Bo(G,1) Po(U,1) z(U,1)
minimize (a*b*sum(log_sum_exp([log(n)-log(Po.*fenzi_col), log(lamda_g+guiyihuaxishu_col*(P-sum(Pc))./(Ug-1))-log(exp(z)),-log(fenzi_col.*Po)-log(fenzi_col.*Po)]))+ ...
     (c/2).*(sum(Bc + Bo)-B )^2 ...
    + lamda2*(sum(Bo)-B)+lamda1*(sum(Pc)-P) + c/2*(sum(Po + Pc)-P)^2)
subject to
Bo>=B/G/3
Po>=P/U/3
sum(Bc + Bo)-B <=0
sum(Po + Pc)-P<=0
z<=log(Po.*fenzi_col)+log((balloct*Bo))
cvx_end