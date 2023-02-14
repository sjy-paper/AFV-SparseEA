function [Population,Dec,Mask,FrontNo,CrowdDis,dis] = EnvironmentalSelection(Population,Dec,Mask,N,len)
% 环境选择


    %% 删除冗余解
    [~,uni] = unique(Population.objs,'rows');
    Population = Population(uni);
    Dec        = Dec(uni,:);
    Mask       = Mask(uni,:);
    N          = min(N,length(Population));

    %% 非支配排序
    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
    Next = FrontNo < MaxFNo;
    
    %% 计算每个解的拥挤距离
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);
    
    %% 根据拥挤距离选择解
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    %% 计算dis
    Objs = Population.objs;
    se = Objs([false(1,len) Next(len+1:end)],:);  %新种群中，选中的个体
    ext = Objs(~Next(1:len),:);     % 从原始种群 剔除的个体

    distance = min(pdist2(se,ext),[],2);
    dis = mean(distance);
    
    
    %% 下一代种群
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
    Dec        = Dec(Next,:);
    Mask       = Mask(Next,:);
end