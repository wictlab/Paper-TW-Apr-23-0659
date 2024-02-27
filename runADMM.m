function [totcom_last,totpos_last,Pc,Po] = runADMM(Mt1,v1,balloct,a,b,G,U,P,B,lamda_g,Uset,guiyihuaxishu_col,Ug,fenzi_col)
c = 120; %70 %120 %50 %250  %%penalty parameter of ADMM
x=Mt1^3*v1;
%% initial value
Bo=(B*(1-b)/G).*ones(G,1);
Pc = (P*b/U).*ones(U,1); %mW same for all user 通信初始功率
Po = (P*(1-b)/U).*ones(U,1); %mW same for all user 感知初始功率
Bc=(B*b/G).*ones(G,1);%组间正交 组内相同
lamda1 = 0;%pc
lamda2 = 0;%po
lamda3 = 0;%bc
lamda4 = 0;%bo
%%
IterationNum = 50;
totcom_store = zeros(1,IterationNum); %%contain EE of each iteration
totpos_store= zeros(1,IterationNum);
converged = false;
tol=0.1;
tol_com=1000;
iterationIndex=1;
[totcom,totpos] = cal_EE(Mt1,v1,Ug,P,Bc, Bo, Pc,Po,lamda_g,a,b,balloct,guiyihuaxishu_col,fenzi_col)
while ~converged
    iterationIndex
totcom_last=totcom;
totpos_last=totpos;
[Po] = update_po(x,Bo,Bc,Pc,a,b,lamda_g,lamda1,lamda2,lamda3,lamda4,c,balloct,guiyihuaxishu_col,Ug,fenzi_col,P,B,G,U);
[Bo] = update_bo(x,Po,Bc,Pc,a,b,lamda_g,lamda1,lamda2,lamda3,lamda4,c,balloct,guiyihuaxishu_col,Ug,fenzi_col,P,B,G,U);
[Pc] = update_pc(Bo,Bc,Po,a,b,lamda_g,lamda1,lamda2,lamda3,lamda4,c,balloct,guiyihuaxishu_col,Ug,fenzi_col,P,B,G,U);
 [Bc] = update_bc(Po,Bo,Pc,a,b,lamda_g,lamda1,lamda2,lamda3,lamda4,c,balloct,guiyihuaxishu_col,Ug,fenzi_col,P,B,G,U);   
 
 
 lamda1 = update_lamda1(lamda1,c,Pc,Po,P,b);
    lamda2 = update_lamda2(lamda2,c,Pc,Po,P,b);
    lamda3 = update_lamda3(lamda3,c,Bc,Bo,B,b);
    lamda4 = update_lamda4(lamda4,c,Bc,Bo,B,b);
    
    [totcom,totpos]= cal_EE(Mt1,v1,Ug,P,Bc, Bo, Pc,Po,lamda_g,a,b,balloct,guiyihuaxishu_col,fenzi_col); %calculate EE each iteration. Can use "out_EE" alternatively.
    totcom_store(iterationIndex+1)= totcom;
    totpos_store(iterationIndex+1)= totpos;
    erro_com=norm(totcom - totcom_last);
    erro_pos=norm(totpos - totpos_last);
    if ((erro_pos < tol)||(erro_com<tol_com)) || iterationIndex >= IterationNum
       converged = true;
    end
    totcom
    totpos
    iterationIndex = iterationIndex + 1;
end