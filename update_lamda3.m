function out_lamda3 = update_lamda1(lamda3,c,Bc,Bo,B,b)
out_lamda3 = lamda3+ c*(b*B-sum(Bc));