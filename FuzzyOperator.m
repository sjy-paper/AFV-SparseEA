function dec = FuzzyOperator(OffDec,theta)
% 变量模糊算子
    Problem = PROBLEM.Current();
    
    R    = Problem.upper-Problem.lower;
    
    gamma{1} = Problem.upper - R*theta^-1.*floor(theta^1*R.^-1.*(Problem.upper - OffDec));
    gamma{2} = Problem.upper - R*theta^-1.*ceil(theta^1*R.^-1.*(Problem.upper - OffDec));


    miu{1}    = 1./(1+ exp(gamma{1} - OffDec));
    miu{2}    = 1./(1+exp(OffDec - gamma{2}));

    logical = miu{1}-miu{2}>0;

    OffDec  = gamma{2};
    OffDec(find(logical)) = gamma{1}(find(logical));
        

    dec = OffDec;
end