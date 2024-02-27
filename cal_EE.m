function [totcom,totpos] = cal_EE(Mt1,v1,Ug,P,Bc, Bo, Pc,Po,lamda_g,a,b,balloct,guiyihuaxishu_col,fenzi_col)

com=a*(balloct*Bc).*log2(1+Pc.*fenzi_col./(lamda_g+guiyihuaxishu_col.*((P-sum(Po))./(Ug))));
pos=a*(balloct*Bo./2).*(log2(1+Mt1^3*v1*Po.*fenzi_col./(lamda_g+guiyihuaxishu_col.*((P-sum(Pc))./(Ug)))));
totcom=sum(com);
totpos=sum(pos);