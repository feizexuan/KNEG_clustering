function [idx,T1,T2,T3] = KNEG_cl(data,K)
%idx：Clustering label
%T1：Time taken to construct KNEG
%T2：Time taken to construct subclusters
%T3：Time taken to merge or split
%writen by Zexuan Fei
disp('dyknn start!!');
%data pre-processing
data(:,all(data==0, 1))=[];
n=size(data,1);
data=(data-min(data))./(max(data)-min(data));
%需要调整的部分
shownp_i=0;
shownp_final=0;%展示峰值点
showimage=0;%展示每次合并的过程
show_nparrow=0;%展示np_group之间的指向关系
k=floor(sqrt(n));
kf=0;
%% Construct KNEG
t1=clock;
[knn_neigh,knn_dist] = knnsearch(data,data,'k',k);
knn_neigh0=knn_neigh;%Backup result
knn_oppsite=ones(n,k);%如果等于1说明对方有指向自己的边
for i=2:k
    for j=1:n
        neigh=knn_neigh(j,i);
        if neigh~=0
            loc=find(knn_neigh(neigh,:)==j);
            if isempty(loc) || loc>i %如果对方的k近邻之后才访问到，或者根本不包含。则对方不会连接自己
                %对方结点中对应位置删除
                knn_neigh(neigh,loc)=0;
                knn_dist(neigh,loc)=inf;
                knn_oppsite(j,i)=0;
            end
        end
    end
end

%% obtain peak density point
colordef black;
n_np_list=[];
is_np=[];
is_np2=[];
is_np_list=[];
np_list=[];
% last_flag=0;
for kk=1:k
    if shownp_i && showimage
        figure;
        plot(data(:,1),data(:,2),'linestyle',"none",'marker','.','Color','white','Markersize',10);
        hold on;
    end
    %采用矩阵计算is_np
    is_np=zeros(n,1);%初始化全0
    knn_oppsite2=~knn_oppsite;%取反，对方有连接自己的边的knn_oppsite2=0
    judge=sum(knn_neigh(:,1:kk).*knn_oppsite2(:,1:kk),2);%这行代码逻辑稍微有点复杂
    %judge：自己连接对方knn_neigh>0，如果对方没有连接自己knn_oppsite2>0,此时judge>0,就不是峰值点
    %judge：自己连接对方knn_neigh>0，如果对方也连接自己knn_oppsite2=0,此时judge可能0,可能是峰值点
    %sum就是考虑所有近邻中的情况，必须都是0才是密度峰值点
    is_np(judge==0)=1;%峰值点的judge值一定=0，
    for i=1:n
        if is_np(i)==1 && shownp_i==1 && showimage
            plot(data(i,1),data(i,2),'linestyle',"none",'marker','.','Color','red','Markersize',10);
            hold on;
        end
    end
    is_np_list=[is_np_list,is_np];%存储着密度峰值点的历史记录
    n_np=length(find(is_np==1));
    n_np_list=[n_np_list,n_np];
    if kk~=1
        if n_np_list(end-1)-n_np_list(end)<=n^0.5 %suspensive condition
            kf=kk;
            np_list=(find(is_np==1));%np_list存储峰值点的下标，is_np存储logical数组
%             last_flag=1;
            break;
        end
    end
end
if n_np<K
    is_np=is_np2;
    np_list=(find(is_np==1));
    n_np=length(find(is_np==1));
end
t2=clock;
T1=etime(t2,t1);

%% Get subclusters(saved in group and np_group)
group=zeros(n,1);
group(is_np==1)=1:length(find(is_np==1));
%将较近的密度峰值点归为一类
count=0;
t3=clock;
for i=1:n_np    
    startnp=np_list(i);
    np_around=find(knn_neigh(startnp,1:kf)>0);
    for j=2:size(np_around,2)%遍历周围np，修改数值
        if is_np(knn_neigh(startnp,np_around(j)))==1
            group(knn_neigh(startnp,np_around(j)))=group(startnp);
            count=count+1;
        end
    end
end
%开始计算归属情况
for i=1:n
    if group(i)>0
        continue;
    end
    next_is=i;
    next_i_loc=[];
    loc=find(knn_neigh(i,:)>0);
    for j=1:size(loc,2)
        if knn_oppsite(i,loc(j))==0%走单向箭头
            next_i_loc=loc(j);
            break;
        end
    end
    if isempty(next_i_loc)%如果没找到最近单向箭头，那就直接访问最近的近邻（loc1存储的是自己）
        next_i_loc=loc(2);
    end
    next_i=knn_neigh(i,next_i_loc);
    while group(next_i)==0 
        next_i0=next_i;
        next_is=[next_is,next_i];
        loc=find(knn_neigh(next_i,:)>0);
        for j=2:size(loc,2)
            if knn_oppsite(next_i,loc(j))==0
                next_i_loc=loc(j);
                break;
            end
        end
        if isempty(next_i_loc)
            next_i_loc=loc(2);
        end
        next_i=knn_neigh(next_i,next_i_loc);
        %如果关系中出现重复的
        trace_back=0;
        while ~isempty(find(next_is==next_i,1))%如果next_i本次访问过
            loc=find(knn_neigh(next_i,:)>0);%从这个重复的点开始继续寻找，其他路径
            flag=0;
            for j=2:size(loc,2)
%                 if loc(j)~=next_i_loc %&& knn_oppsite(next_i0,loc(j))==0 
                    if isempty(find(next_is==knn_neigh(next_i,loc(j)),1))
                        next_i_loc=loc(j);
                        flag=1;
                        break;
                    end
%                 end
            end
            if flag==1 %如果找到了本次未访问过的结点，则接着访问
                next_i=knn_neigh(next_i,next_i_loc);
                trace_back=0;
            else
                if size(next_is,2)>trace_back%如果可以回溯
                    next_i=next_is(end-trace_back);%出错了
                    trace_back=trace_back+1;
                else
                    %如果无法回溯，即找遍next_is都无法走出循环
                    %那我们让这些next_is独立成簇
                    is_np(next_is(1))=1;
                    group(next_is(1))=length(find(is_np==1));
                    group(next_is)=group(next_is(1));
                    n_np=n_np+1;
                    np_list=[np_list;next_is(1)];
                    trace_back=0;
                    break;
                end

            end
        end
    end
    group(next_is)=group(next_i);
end

%weight np_group
np_list0=(find(is_np_list(:,kf-1)==1));

is_visited_np=[];
num_delete=0;
for i=1:n_np
    np_cur=np_list(i);
    if isempty(find(is_visited_np==group(np_cur),1))
        is_visited_np=[is_visited_np,group(np_cur)];
    else
        np_list(i)=0;
        num_delete=num_delete+1;
        is_np(np_cur)=0;
        continue;
    end
end
np_list(np_list==0)=[];
n_np=n_np-num_delete;
%确定np_group
max_pec=tabulate(group(:));
max_pec=max(max_pec(:,3));
maxsize=ceil(size(data,1)*max_pec/100)+1;
np_group_mem=zeros(maxsize,n_np);
np_group_pri=zeros(maxsize,n_np);
for i=1:n_np
    np_cur=np_list(i);
    mem=find(group==group(np_cur));
    [~,pri]=ismember(mem,knn_neigh0(np_cur,1:kf));
    fill=zeros(maxsize-size(mem,1),1);
    mem=[mem;fill];
    pri=[pri;fill];
    np_group_mem(:,i)=mem;
    np_group_pri(:,i)=pri;
end
t4=clock;
T2=etime(t4,t3);


%% Calculate the connection between np_group
t5=clock;
relation=zeros(n_np,n_np);%relation(i,j)存i指向j的值
relation_count=zeros(n_np,n_np);
np2np=zeros(n_np,n_np);
is_visited_np=[];
np_group_pri(np_group_pri==0)=kf+0.5;%结点优先级
for i=1:n_np
    if isempty(find(is_visited_np==group(np_list(i))))
        is_visited_np=[is_visited_np,group(np_list(i))];
    else
        continue;
    end
    np_around=np_group_mem(:,i);
    for j=1:find(np_around==0)-1
        np_around_around=knn_neigh(np_around(j),1:floor(2*kf-2));%!!!!!!
        for m=1:size(np_around_around,2)%遍历周围节点不包括自己，
            if np_around_around(m)~=0%取其中有指向关系的
                [row,col]=find(np_group_mem==np_around_around(m));%查看是否在其他np_group中
                for n=1:size(row,1)
                    if(col(n)==i)%访问到自己跳过
                        continue;
                    end                    
                    relation(i,col(n))=relation(i,col(n))+(2*kf+2-np_group_pri(row(n),col(n))-np_group_pri(j,i))/2;
                    relation_count(i,col(n))=relation_count(i,col(n))+1;
                end
            end
        end
    end
end
%绘制密度峰值点指向图，每个小簇之间出边对抗，两个簇之间多的出边指向相对较少的出边
%这是所有np_group之间的关系
% relation(i,j)表示i小簇到j小簇的出边总和
if show_nparrow
    for i=1:n_np
        x=data(np_list(i),:);
        x_row=find(relation(i,:)>0);
        if ~isempty(x_row)
            for j=1:size(x_row,2)
                y=data(np_list(x_row(j)),:);
                sim=relation(x_row(j),i)/relation(i,x_row(j));
                a=0.0001;
                if (sim>=a && sim<=(1/a))
                    if relation(x_row(j),i)>relation(i,x_row(j))
                        drawarrow(y,x,'yellow','arrow');
                    else
                        drawarrow(x,y,'yellow','arrow');
                    end
                end
            end
        end
    end
end
%% np归类
%从所有强连接点的入边开始，根据相似度关系进行连接
relation_norm=data_norm(relation);
in_out_rate=zeros(n_np,n_np);
for i=1:n_np
    x=data(np_list(i),:);
    x_row=find(relation(i,:)>0);
    if ~isempty(x_row)
        for j=1:size(x_row,2)
            y=data(np_list(x_row(j)),:);
            in_out_rate(i,x_row(j))=relation(x_row(j),i)/relation(i,x_row(j));
        end
    end
end
similarity=in_out_rate;
np_similar=zeros(n_np,n_np);
similarity(similarity>1)=0;
np_most_outline=zeros(n_np,n_np);
max_chubian=[];
if showimage
    figure;
    cmap=hsv(max(group));
    gscatter(data(:,1),data(:,2),group,cmap,'.',10);%备注 展示小簇效果
    hold on;
end
relation2=relation;
same_cluster=0.8;
for i=1:n_np
    x=data(np_list(i),:);
    a=similarity(i,:);
    %计算连接关系，指出相似度较小的箭头
    x_row=find((similarity(i,:)<1)&(similarity(i,:)>same_cluster));
    for j=1:size(x_row,2)
        y=data(np_list(x_row(j)),:);
%         drawarrow(x,y,'yellow','arrow');
%         line([x(1),y(1)],[x(2),y(2)],'color','yellow');
        np_similar(i,x_row(j))=1;
%         np_similar(x_row(j),i)=1;
    end
    %备注 计算连接关系，指向出边最多的指向的簇(绿色箭头)，逻辑上形成nps_group
    [x_row_max,index_max]=max(relation(i,:));
    relation2(i,index_max)=0;%删掉该值方便取次最大值
    x_num=length(find(group==group(np_list(i))));
    y_num=length(find(group==group(np_list(index_max))));
    x_chu=relation(i,index_max);
    y=data(np_list(index_max),:);   
    if x_row_max==0
        max_chubian=[max_chubian;i,i,1,x_num,x_num];
        y=x;
    else
        max_chubian=[max_chubian;i,index_max,in_out_rate(i,index_max),x_num,y_num];
    end
    if showimage
        drawarrow(x,y,'green','arrow');
    end
end

%% Merger or split
s=max_chubian(:,1);
t=max_chubian(:,2);
single_chu=max_chubian(~ismember(s,t),:);
single_chu(:,2)=[];
if showimage 
    chu=np_list(single_chu(:,1));
    plot(data(chu,1),data(chu,2),'linestyle',"none",'marker','.','Color','yellow','Markersize',20);
    hold on;
    unnorm=single_chu(:,3)>single_chu(:,4);
    chu=chu(unnorm==1);
    plot(data(chu,1),data(chu,2),'linestyle',"none",'marker','.','Color','red','Markersize',20);
end
for i=1:size(single_chu,1)
    relation(single_chu(i,1),:)=relation(single_chu(i,1),:).*single_chu(i,2);
end
%construct weakly connected components
w=1./(max_chubian(:,3));
G = graph(s,t);
bins = conncomp(G);
finished=0;

if max(bins)<=K %split
    toall=max_chubian(:,end)+max_chubian(:,end-1);
    max_chubian(:,end-2)=max_chubian(:,end-2).*toall;
    np_relation=max_chubian(:,1:3);
    [np_relation,loc]=sortrows(np_relation,3);%np_relation第三列中存储的计算出的最终相似度（绝对从属度*相对从属度）
    %排序之后，选择边进行删除
    for i=1:size(np_relation,1)
        if max(bins)==K
            finished=1;
            break;
        end
        G = rmedge(G,[s(loc(i))],[t(loc(i))]); %图中删除边
        bins = conncomp(G);
    end
end
%分类完毕生成分类结果
np_cluster=1;
np_clusters=[];
np_datanum=[];
i=1;
nps=cell(max(bins),1);
while ~isempty(find(bins==np_cluster,1))
    index=find(bins==np_cluster);
    for j=1:size(index,2)
        np_clusters=[np_clusters,max_chubian(index(j),1)];
    end
    np_clusters=unique(np_clusters);
    nps(np_cluster)={np_clusters};
    sum_datanum=0;
    for j=1:size(np_clusters,2)
        group(group==group(np_list(np_clusters(j))))=np_cluster;%出错了！
    end
    sum_datanum=length(find(group==np_cluster));
    np_datanum=[np_datanum;sum_datanum];
    np_clusters=[];
    np_cluster=np_cluster+1;
end
% merge
nps_relation=[];
if finished==0%如果分裂过程未进行，则开始合并，先计算nps_group之间相似度，存储在nps_relation中
    nps_relation=[];
    for i=1:size(nps,1)-1
        nps_i=nps{i};
        for j=i+1:size(nps,1)
            nps_j=nps{j};
            i2j=relation(nps_i,nps_j);
            j2i=relation(nps_j,nps_i);
            i2j_num=sum(i2j(:));
            j2i_num=sum(j2i(:));
            toall=i2j_num+j2i_num;
            if i2j_num<=j2i_num
                in_out_rate2=i2j_num/j2i_num;
            else
                in_out_rate2=j2i_num/i2j_num;
            end
            if toall~=0
                nps_relation=[nps_relation;[i,j,in_out_rate2,toall]];
            end
        end
    end
end
group_nps=1:size(nps,1);
if ~isempty(nps_relation) && finished==0%计算完nps_relation，开始合并
    nps_relation(:,3)=nps_relation(:,3).*nps_relation(:,4);
    nps_relation_sorted=sortrows(nps_relation,-3);%对于nps_relation排序，第三列存储相似度

    is_incluster=zeros(size(nps,1),1);
    nums_sum=[];
    %合并算法
    u=1:size(np_datanum,1);
    nums_sum=[u',np_datanum];
    if showimage%绘制底图，为下面展示合并顺序做准备
        figure;
        cmap=hsv(max(group));
        gscatter(data(:,1),data(:,2),group,cmap,'.',10);%备注 展示小簇效果
        hold on;
    end
    for i=1:size(nps_relation_sorted,1)
        if showimage%备注 可以在end处打断点观察每次合并的顺序
            nps_x=nps{nps_relation_sorted(i,1)};
            nps_y=nps{nps_relation_sorted(i,2)};
            nodex=np_list(nps_x(1));
            nodey=np_list(nps_y(1));
            x=data(nodex,:);
            y=data(nodey,:);
            line([x(1),y(1)],[x(2),y(2)],'color','white','linewidth',3);
        end
        if ~isempty(nums_sum)
            max_npnum0=max(nums_sum(:,2));
        end
        group_nps0=group_nps;
        if is_incluster(nps_relation_sorted(i,1))==0 && is_incluster(nps_relation_sorted(i,2))==0
            a=group_nps(nps_relation_sorted(i,1));
            b=group_nps(nps_relation_sorted(i,2));
            group_nps(nps_relation_sorted(i,2))=group_nps(nps_relation_sorted(i,1));
            is_incluster(nps_relation_sorted(i,1))=1;
            is_incluster(nps_relation_sorted(i,2))=1;
            %用于记录簇中数据点的数目变化
            b_num=np_datanum(nps_relation_sorted(i,2));
            nums_sum(nums_sum(:,1)==a,2)=nums_sum(nums_sum(:,1)==a,2)+b_num;
            nums_sum(nums_sum(:,1)==b,:)=[];
        elseif is_incluster(nps_relation_sorted(i,1))==1 && is_incluster(nps_relation_sorted(i,2))==1
            if group_nps(nps_relation_sorted(i,1))~=group_nps(nps_relation_sorted(i,2))
                a=group_nps(nps_relation_sorted(i,1));
                b=group_nps(nps_relation_sorted(i,2));
                group_nps(group_nps==group_nps(nps_relation_sorted(i,2)))=group_nps(nps_relation_sorted(i,1));
                %          group_nps(nps_relation_sorted(i,2))=group_nps(nps_relation_sorted(i,1));
                %用于记录簇中数据点的数目变化
                b_num=nums_sum(nums_sum(:,1)==b,2);
                nums_sum(nums_sum(:,1)==a,2)=nums_sum(nums_sum(:,1)==a,2)+b_num;
                nums_sum(nums_sum(:,1)==b,:)=[];
            end
        elseif  is_incluster(nps_relation_sorted(i,1))==0
            a=group_nps(nps_relation_sorted(i,2));
            b=group_nps(nps_relation_sorted(i,1));
            group_nps(nps_relation_sorted(i,1))=group_nps(nps_relation_sorted(i,2));
            is_incluster(nps_relation_sorted(i,1))=1;
            b_num=np_datanum(nps_relation_sorted(i,1));
            nums_sum(nums_sum(:,1)==a,2)=nums_sum(nums_sum(:,1)==a,2)+b_num;
            nums_sum(nums_sum(:,1)==b,:)=[];
        
        elseif  is_incluster(nps_relation_sorted(i,2))==0
            a=group_nps(nps_relation_sorted(i,1));
            b=group_nps(nps_relation_sorted(i,2));
            group_nps(nps_relation_sorted(i,2))=group_nps(nps_relation_sorted(i,1));
            is_incluster(nps_relation_sorted(i,2))=1;
            b_num=np_datanum(nps_relation_sorted(i,2));
            nums_sum(nums_sum(:,1)==a,2)=nums_sum(nums_sum(:,1)==a,2)+b_num;
            nums_sum(nums_sum(:,1)==b,:)=[];
        end
        max_npnum=max(nums_sum(:,2));

        if length(unique(group_nps))==K
            merge_ref=nps_relation_sorted(:,3);
            loc=i;
            break;
        end
    end
else
    merge_ref=0;
    loc=0;
end
 %%最终聚类
 group_f=zeros(size(group,1),1);
 [group_nps_u]=unique(group_nps);
 np_cluster=1;
 for i=1:size(group_nps_u,2)
     group_nps_id=group_nps_u(i);
     onegroup=find(group_nps==group_nps_id);
     for j=1:size(onegroup,2)
         group_f(group==onegroup(j))=np_cluster;
     end
     np_cluster=np_cluster+1;
 end
 t6=clock;
 T3=etime(t6,t5);
% gscatter(data(:,1),data(:,2),group,cmap,'.',10);
idx=group_f;
disp('dyknn end!!');
end