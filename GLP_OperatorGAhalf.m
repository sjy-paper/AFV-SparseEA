function [Offspring,outIndexList,chosengroups] = GLP_OperatorGAhalf(Parent1,Parent2, numberOfGroups)

    [proC,disC,~,disM] = deal(1,20,1,20);
    [N,D]   = size(Parent1);
    Problem = PROBLEM.Current();
    
    beta = zeros(N,D);
    mu   = rand(N,D);
    beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
    beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
    beta = beta.*(-1).^randi([0,1],N,D);
    beta(rand(N,D)<0.5) = 1;
    beta(repmat(rand(N,1)>proC,1,D)) = 1;
    Offspring = (Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2;
    Lower = repmat(Problem.lower,N,1);
    Upper = repmat(Problem.upper,N,1);
    [outIndexList,~] = CreateGroups(numberOfGroups,Offspring,D); 
    
    
    chosengroups = randi(numberOfGroups,size(outIndexList,1),1);
    
    
    Site = outIndexList == chosengroups;
    mu    = rand(N,1);
    mu    = repmat(mu,1,D);
    temp  = Site & mu<=0.5;
    Offspring       = min(max(Offspring,Lower),Upper);
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
        (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
    temp = Site & mu>0.5;
    Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
        (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));    

end


function [outIndexArray,numberOfGroupsArray] = CreateGroups(numberOfGroups, xPrime, numberOfVariables)

    outIndexArray = [];
    numberOfGroupsArray = [];
    noOfSolutions = size(xPrime,1);
    for sol = 1:noOfSolutions        
        varsPerGroup = floor(numberOfVariables/numberOfGroups);
        vars = xPrime(sol,:);
        [~,I] = sort(vars);
        outIndexList = ones(1,numberOfVariables);
        for i = 1:numberOfGroups-1
            outIndexList(I(((i-1)*varsPerGroup)+1:i*varsPerGroup)) = i;
        end
        outIndexList(I(((numberOfGroups-1)*varsPerGroup)+1:end)) = numberOfGroups;    
        outIndexArray = [outIndexArray;outIndexList];
        numberOfGroupsArray = [numberOfGroupsArray;numberOfGroups];    
    end
end