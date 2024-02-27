function out_lamda4 = update_lamda4(lamda4,c,Bc,Bo,B,b)
out_lamda4 = lamda4+ c*((1-b)*B-sum(Bo));