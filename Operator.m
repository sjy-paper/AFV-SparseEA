function [OffDec,OffMask] = Operator(ParentDec,ParentMask,Fitness,theta,REAL)
    [N,~]       = size(ParentDec);
    Parent1Dec  = ParentDec(1:floor(end/2),:);
    Parent2Dec  = ParentDec(floor(end/2)+1:floor(end/2)*2,:);
    
    Parent1Mask = ParentMask(1:floor(end/2),:);
    Parent2Mask = ParentMask(floor(end/2)+1:floor(end/2)*2,:);
    
     if REAL
        [OffDec,groupIndex,chosengroups] = GLP_OperatorGAhalf(Parent1Dec,Parent2Dec,4); % 4 -- numberofgroups
        OffDec = FuzzyOperator(OffDec,theta);
     end
     
    OffMask = Parent1Mask;
    for i = 1 : N/2
        if rand < 0.5
            index = find(Parent1Mask(i,:)&~Parent2Mask(i,:));
            index = index(TS(-Fitness(index)));
            OffMask(i,index) = 0;
        else
            index = find(~Parent1Mask(i,:)&Parent2Mask(i,:));
            index = index(TS(Fitness(index)));
            OffMask(i,index) = Parent2Mask(i,index);
        end
    end
    
    if REAL
        chosenindex = groupIndex == chosengroups;
        for i = 1 : N/2
            if rand < 0.5
                index = find(OffMask(i,:)&chosenindex(i,:));
                index = index(TS(-Fitness(index)));
                OffMask(i,index) = 0;
            else
                index = find(~OffMask(i,:)&chosenindex(i,:));
                index = index(TS(Fitness(index)));
                OffMask(i,index) = 1;
            end
        end
    else
         OffDec = ones(size(OffMask));
    end
    
    
end

function index = TS(Fitness)
    if isempty(Fitness)
        index = [];
    else
        index = TournamentSelection(2,1,Fitness);
    end
end