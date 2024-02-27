function out_lamda1 = update_lamda1(lamda1,c,Pc,Po,P,b)
out_lamda1 = lamda1+ c*(b*P-sum(Pc));