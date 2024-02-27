run('C:\Users\jeffrey\Desktop\code\Map.m');
%H=generate_mmwave_channel(16, 8, 4, 4);
% Mr=[20 20];%基站天线数
% Mt=[18 14];%用户天线数
Mr=[20];%基站天线数
Mt=[18];%用户天线数
v=pi^2*cos(pi/3)*(Mr-1)/12/2;
P1=[5:10:55];
com1=zeros(size(P1,2),1);
pos1=zeros(size(P1,2),1);
com3=zeros(size(P1,2),1);
pos3=zeros(size(P1,2),1);
com5=zeros(size(P1,2),1);
pos5=zeros(size(P1,2),1);
for nr=1:size(Mr,2)
    f0=60;%GHz中心频率
    nums_ue=36;%用户数量
    nums_bs=3;%基站数
    ax=32.4;%信号衰减指数
    bx=2;%路径损耗指数
    d0=1;%参考距离
    Nx=6.8;%阴影衰落
    c=299792458;%m/s 光速
    
    %%%%%%%%%%%%%%%%%%%%%%%%%系统初始化%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    distance=zeros(nums_ue,nums_bs);%定义一个距离矩阵
    UE= vertcat(u1, u2,u3);
    for i=1:nums_bs
        bs_position=[BS(i,1),BS(i,2)];
        for j=1:nums_ue
            ue_position=[UE(j,1),UE(j,2)];
            distance(j,i)=sqrt((BS(i,1)-UE(j,1))^2+(BS(i,2)-UE(j,2))^2);
        end
    end
    
    Loss=zeros(nums_ue,nums_bs);%路径损耗系数矩阵
    for i=1:nums_bs
        for j =1:nums_ue
            Loss(j,i)=ax+10*bx*log10(distance(j,i)/d0)+20*log10(4*pi*f0*d0*10^6/c)+Nx;
        end
    end
    H0 = cell(nums_ue, nums_bs);%定义一个cell数组放产生的信道矩阵H
    
    for i=1:nums_ue
        for j=1:nums_bs
            H0{i,j}=generate_mmwave_channel(Mr(nr),Mt(nr),randi([5,10]),randi([5,10]));
            %H0{i,j}=generate_mmwave_channel(Mr,Mt,10,10);
        end
    end
    
    H=cell(nums_ue,nums_bs);%路径损耗系数修正后的矩阵
    for i=1:nums_ue
        for j=1:nums_bs
            H{i,j}=H0{i,j}/Loss(i,j);
        end
    end
    %%%%%%%%%%%%%%%%%%%用户分组第一阶段%%%%%%%%%%%%%%%%%%%%%%
    %步骤一
    %T联合干扰矩阵，i,j代表地i个小区内地第j用户的联合干扰矩阵
    T=cell(nums_ue,nums_bs);
    for i=1:nums_ue
        T{i,1}=[H{i,2}',H{i,3}'];
    end
    for i=1:nums_ue
        T{i,2}=[H{i,1}',H{i,3}'];
    end
    for i=1:nums_ue
        T{i,3}=[H{i,1}',H{i,2}'];
    end
    %步骤二
    U=cell(nums_ue,nums_bs);
    S=cell(nums_ue,nums_bs);%对角矩阵
    V=cell(nums_ue,nums_bs);%右奇异矩阵
    
    for i=1:nums_ue
        for j=1:nums_bs
            [U{i,j},S{i,j},V{i,j}]=svd(T{i,j});
        end
    end
    
    lamda=zeros(nums_ue,nums_bs);
    for i=1:nums_ue
        for j=1:nums_bs
            lamda(i,j)=S{i,j}(Mt(nr),Mt(nr));
        end
    end
    
    theta2=zeros(nums_ue,nums_bs);%小区间干扰
    for i=1:nums_ue
        for j=1:nums_bs
            theta2(i,j)=lamda(i,j)*lamda(i,j);
        end
    end
    
    UE_information=zeros(nums_ue,2);
    UE_num=zeros(nums_ue,1);
    UE_Cellnum=zeros(nums_ue,1);
    for  i=1:nums_ue
        UE_num(i,1)=i;
    end
    
    for i=1:nums_ue
        [min_value,min_index]=min(theta2(i,:));
        min_indices=find(theta2(i,:)==min_value);
        UE_Cellnum(i,1)=min_indices;
    end
    UE_information=[UE_num,UE_Cellnum];
    
    W=cell(nums_ue,1);%波束成形矩阵
    for i=1:nums_ue
        W{i,1}=V{i,1}(:,Mt(nr));
    end
    
    rng(1)
    nums_cell1=0;
    for i=1:nums_ue
        if(UE_information(i,2)==3)
            nums_cell1=nums_cell1+1;
        end
    end
    cell1=cell(nums_cell1,9);%第九列放W
    j=0;
    for i=1:nums_ue
        if(UE_information(i,2)==3)
            j=j+1;
            numsrows=size(W{i,1},1);
            cell1{j,9}=W{i,1}(1:numsrows/2, :);%%%%%
            cell1{j,1}=j;
            cell1{j,2}=H{i,3};%第二列是信道矩阵
            cell1{j,7}=theta2(i,3);%%第七列放一下小区间干扰
        end
    end
    for i=1:nums_cell1
        cell1{i,3}=(cell1{i,2}')*cell1{i,2};%第三列所示H'H
    end
    for i=1:nums_cell1
        [V,D]=eig(cell1{i,3});
        cell1{i,4}=D;%第四列就是lamda
        cell1{i,5}=V;%第五列是xi
    end
    for i=1:nums_cell1
        tao=zeros(Mt(nr),1);
        for j=1:Mt(nr)
            tao=cell1{i,4}(j,j)*cell1{i,5}(j,:);
        end
        cell1{i,6}=tao;
    end
    row=zeros(nums_cell1,nums_cell1);%信道相关性矩阵
    for i=1:nums_cell1
        for j=1:nums_cell1
            if i==j
                row(i,j)=0;
            end
            row(i,j)=norm(((cell1{i,6})'*cell1{i,6})'*((cell1{j,6})'*cell1{j,6}))/(norm(((cell1{i,6})'*cell1{i,6}))*norm(((cell1{j,6})'*cell1{j,6})));
        end
    end
    B=100000000;%带宽总和
    for PIDX=1:size(P1,2)
        P=P1(PIDX);
        %P=10;%每个基站传输功率
        n0=-174;
        n=10^(n0/10)/1000;
        disturb=0.1;
        D=40;
        omega=2*pi;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%小区三分组0.2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        result_3=[];
        result_MI_3=[];
        a=0.2;%test 分组阈值
        nums_group=0;
        group_information=cell(100,100);%第一列是组数，第二列是成员信息
        %将第一个用户加入组1
        nums_group=nums_group+1;
        group_information{1,1}=nums_group;
        group1=[];group1(1,1)=1;
        group1_nums=1;
        group_matrix=zeros(nums_cell1,1);%用来标记用户是否进行了分组
        group_matrix(1,1)=1;
        for j=2:nums_cell1
            if row(1,j)<a
                group1_nums=group1_nums+1;
                group1(group1_nums,1)=j;
                group_matrix(j,1)=1;
            end
        end
        group_information{nums_group,2}=group1;
        %%第二阶段步骤3
        while(any(group_matrix==0))
            countzeros=sum(group_matrix==0);
            nums_group=nums_group+1;
            [startindex,colindex]=find(group_matrix==0);
            group_new=[];
            group_new(1,1)=startindex(1,1);
            nums_group_new=1;
            group_matrix(startindex(1))=1;
            j=countzeros;
            [startindex_find,colindex_find]=find(group_matrix==0);%找到另外的还没分组的用户
            for i=1:j-1
                if row(startindex(1,1),startindex_find(i,1))<a
                    nums_group_new=nums_group_new+1;
                    group_new(nums_group_new,1)=startindex_find(i);
                    group_matrix(startindex_find(i,1))=1;
                end
            end
            group_information{nums_group,2}=group_new;
            j=0;
            clear group_new;
        end
        for i=1:nums_group
            group_information{i,1}=i;
        end
        %%第二阶段步骤4
        for i=1:nums_group
            nums_group_i=length(group_information{i,2});
            for j=1:nums_group_i
                nums_less_than_a=0;
                for k=j+1:nums_group_i
                    if k==j
                        continue
                    end
                    if row(group_information{i,2}(j,1),group_information{i,2}(k,1))>a
                        nums_less_than_a=nums_less_than_a+1;
                    end
                end
                if nums_less_than_a>1 %添加新组
                    nums_group=nums_group+1;
                    group_new=[];
                    group_new=group_information{i,2}(j);
                    group_information{nums_group,1}=nums_group;
                    group_information{nums_group,2}=group_new;
                    nums_group_i= nums_group_i-1;%%回溯防止越界 %从该分组中剔除该用户
                    group_information{i,2}(j)=[];
                end
            end
        end
        P_every_ue=P/nums_cell1;
        every_group=cell(nums_group,3);%第一列组内用户数量，第二列组的带宽,第三列用户成员
        for i=1:nums_group
            every_group{i,1}=length(group_information{i,2});
            every_group{i,2}=B*(length(group_information{i,2})/nums_cell1);%组内复用带宽，组间带宽正交
            every_group{i,3}=group_information{i,2};
        end
        H_cell1=[];
        trace1=0;
        for i=1:nums_cell1
            H_cell1=[H_cell1,cell1{i,2}'];
        end
        H_cell1=(H_cell1')*H_cell1;
        for i=1:length(H_cell1(1,:))
            for j=1:length(H_cell1(:,1))
                if i==j
                    trace1=trace1+H_cell1(i,j);
                end
            end
        end
        b=[0.4];
        
        B =100000000;%hz 带宽
        Ucur = 0;%总用户数
        Uset=zeros(nums_group,100);
        for i=1:nums_group
            Ucur=Ucur+every_group{i,1};%每组用户个数
            for j=1:every_group{i,1}
                Uset(i,j)=every_group{i,3}(j,1);%每组用户编号组
            end
        end
        
        for i=1:size(b,2)
            Bo=(B*(1-b(i))/nums_group).*ones(nums_group,1);
            Pc = (P*b(i)/Ucur).*ones(Ucur,1); %mW same for all user 通信初始功率
            Po = (P*(1-b(i))/Ucur).*ones(Ucur,1); %mW same for all user 感知初始功率
            Bc=(B*b(i)/nums_group).*ones(nums_group,1);%组间正交 组内相同
            [totcom,totpos,Pc1,Po1]=mainADMM(v(nr),Mt(nr),Ucur,Uset,P,Bc, Bo, Pc,Po,B,nums_group,every_group,trace1,b(i),cell1,0);%ADMM
            if nr==1
                com1(PIDX)=totcom;
                pos1(PIDX)=totpos;
            end
            if nr==2
                com3(PIDX)=totcom;
                pos3(PIDX)=totpos;
            end
        end
    end
end
figure;
yyaxis right;
semilogy(P1,pos1,'o -- k','linewidth',5);
hold on;
semilogy(P1,pos3,'^ --b','linewidth',5);
hold on;
set(gca,'ycolor',[0,0,0]);
ylabel('\fontname{Times new Roman}\fontsize{40}{total positioning rate/(bits/s)}')
set(gca,'ydir','reverse')
ylim([55000000,69000000]);

yyaxis left;
plot(P1,com1,'o -k','linewidth',5);
hold on;
plot(P1,com3,'^ -b','linewidth',5);
hold on;
ylim([13000000,28000000]);
ylabel('\fontname{Times new Roman}\fontsize{40}{total data rate/(bits/s)}');
xlabel('\fontname{Times new Roman}\fontsize{40}{power/(W)}');
set(gca,'ycolor',[0,0,0]);
legend('\fontname{Times new Roman}\fontsize{28}{\it{R_{com}},M=20,N=18}','\fontname{Times new Roman}\fontsize{28}{\it{R_{com}},M=20,N=14}','\fontname{Times new Roman}\fontsize{28}{\it{R_{pos}},M=20,N=18}','\fontname{Times new Roman}\fontsize{28}{\it{R_{pos}},M=20,N=14}')
