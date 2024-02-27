function [totcom,totpos,Pc1,Po1]=mainADMM(v1,Mt1,Ucur,Uset,P,Bc, Bo, Pc,Po,B,nums_group,every_group,trace1,b,cell1,flag)
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

a=1/Ucur;%用户权值

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

if flag==0
    [totcom,totpos,Pc,Po] = runADMM(Mt1,v1,balloct,a,b,nums_group,Ucur,P,B,lamda_g,Uset,guiyihuaxishu_col,Ug,fenzi_col);
else
    [totcom,totpos]=cal_EE(Mt1,v1,Ug,P,Bc, Bo, Pc,Po,lamda_g,a,b,balloct,guiyihuaxishu_col,fenzi_col); 
end
Pc1=Pc;
Po1=Po;