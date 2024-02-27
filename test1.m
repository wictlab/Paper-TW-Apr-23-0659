Ucur = 0;%总用户数
%fenzi=zeros(nums_group,100);
%guiyihuaxishu=zeros(nums_group,100);
Uset=zeros(nums_group,100);
for i=1:nums_group
    Ucur=Ucur+every_group{i,1};%每组用户个数
    for j=1:every_group{i,1}
        Uset(i,j)=every_group{i,3}(j,1);%每组用户编号组
    end
end
Ug = zeros(Ucur,1);%每个用户对应组内用户数
for i=1:nums_group
    for j=1:100
        if Uset(i,j)==0
            break;
        end
        Ug(Uset(i,j))=every_group{i,1};%每组用户编号组
    end
end
lamda_g=zeros(Ucur,1);
fenzi_col=zeros(Ucur,1);
guiyihuaxishu_col=zeros(Ucur,1);
for i=1:nums_group%计算R
    for j=1:every_group{i,1}
        hk=cell1{every_group{i,3}(j,1),2};
        for l=1:every_group{i,1}
            if l==j
                continue
            end
            hl=cell1{every_group{i,3}(l,1),2};
            test=norm(hk*hl')/sqrt(trace1);
            guiyihuaxishu_col(every_group{i,3}(j,1))=guiyihuaxishu_col(every_group{i,3}(j,1))+test;
        end%第i行对应第i个用户的组内其他用户的
        fenzi_col(every_group{i,3}(j,1))=norm(hk*hk'/sqrt(trace1));
        lamda_g(every_group{i,3}(j,1))=cell1{every_group{i,3}(j,1),7};
    end
end
B =100000000;%hz 带宽
a=1/Ucur;%用户权值
P = 10; % maximum transmit power 24dBm
lamda_g=lamda_g*2*P+0.0001;
balloct=zeros(Ucur,nums_group);
for i=1:nums_group
   for j=1: size(Uset(i,:),2)
       if Uset(i,j)==0
           break;
       end
       balloct(Uset(i,j),i)=1;
   end
end
G=nums_group;
U=Ucur;
Btot=rand(2*G,1);
Btot=Btot./sum(Btot).*(B);
Bo=Btot(1:G,1);
Bc=Btot(G+1:end,1);

Ptot=rand(2*U,1);
Ptot=Ptot./sum(Ptot).*(P);
Po=Ptot(1:U,1);
Pc=Ptot(U+1:end,1);
b=0.1;
flag=0;
if flag==0
    [totcom,totpos] = runADMM(balloct,a,b,nums_group,Ucur,P,B,lamda_g,Uset,guiyihuaxishu_col,Ug,fenzi_col);
else
    [totcom,totpos]=cal_EE(Ug,P,Bc, Bo, Pc,Po,lamda_g,a,b,balloct,guiyihuaxishu_col,fenzi_col); 
end