function [Bc,Pc] = update_bc_pc(Bo,Po,a,b,lamda2_g,n2,lamda1,lamda2,mu,group_information,rmin,P,B,G,U)
cvx_solver sedumi
cvx_begin quiet
variables Bc(U,1) Pc(U,1)
minimize (a*b*sum(rel_entr(Bc,Bc+Pc))/log(2)+ ...
     mu/2*(sum(Bc + Bo)-B)^2 ...
    + lamda2*(sum(Bo)-B)+lamda1*(sum(Pc)-P) + mu/2*(sum(Po + Pc)-P)^2)
subject to
-a*b*rel_entr(Bc,Bc+Pc)/log(2)-rmin*ones(U,1)>=0
Bc>=0
Pc>=0
groupvec*Bc + groupvec*Bo-B*ones(G,1) <=0
sum(Po + Pc)-P<=0
cvx_end