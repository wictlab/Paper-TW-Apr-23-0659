function out_lamda2 = update_lamda2(lamda2,c,Pc,Po,P,b)
out_lamda2 = lamda2+ c*(P-sum(Po));