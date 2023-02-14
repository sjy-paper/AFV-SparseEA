classdef AFVSparseEA < ALGORITHM
% <multi> <real/binary> <large/none> <constrained/none> <sparse>
% AFV-SparseEA

    methods
        function main(Algorithm,Problem)
           %% 种群初始化
            % 计算每个决策变量的适应度
            TDec    = [];
            TMask   = [];
            TempPop = [];
            Fitness = zeros(1,Problem.D);
            REAL    = ~strcmp(Problem.encoding,'binary');
            for i = 1 : 1+4*REAL
                if REAL
                    Dec = unifrnd(repmat(Problem.lower,Problem.D,1),repmat(Problem.upper,Problem.D,1));
                else
                    Dec = ones(Problem.D,Problem.D);
                end
                Mask       = eye(Problem.D);
                Population = SOLUTION(Dec.*Mask);
                TDec       = [TDec;Dec];
                TMask      = [TMask;Mask];
                TempPop    = [TempPop,Population];
                Fitness    = Fitness + NDSort([Population.objs,Population.cons],inf);
            end
            % 产生初始种群
            if REAL
                Dec = unifrnd(repmat(Problem.lower,Problem.N,1),repmat(Problem.upper,Problem.N,1));
            else
                Dec = ones(Problem.N,Problem.D);
            end
            Mask = zeros(Problem.N,Problem.D);
            for i = 1 : Problem.N
                Mask(i,TournamentSelection(2,ceil(rand*Problem.D),Fitness)) = 1;
            end
            % 实数变量模糊化
            Dec = FuzzyOperator(Dec,ceil(Problem.N*0.1));
            
            Population = SOLUTION(Dec.*Mask);
            
            % 环境选择
            [Population,Dec,Mask,FrontNo,CrowdDis] = EnvironmentalSelection([Population,TempPop],[Dec;TDec],[Mask;TMask],Problem.N,0);

            % 变量初始化
            dis = 0;
            interval = 0;
            rho = 0;
            
            cnt_interval = 0;
            theta = sqrt(Problem.N);
            %% 优化
            while Algorithm.NotTerminated(Population)
                % 锦标赛选择 2N 个父代
                MatingPool       = TournamentSelection(2,2*Problem.N,FrontNo,-CrowdDis);
                % 子代生成 N 个子代
                [OffDec,OffMask] = Operator(Dec(MatingPool,:),Mask(MatingPool,:),...
                                        Fitness , theta , REAL);
                
                Offspring = SOLUTION(OffDec.*OffMask);
                
                predis = dis;
                
                % 环境选择
                [Population,Dec,Mask,FrontNo,CrowdDis,dis] = ...
                    EnvironmentalSelection([Population,Offspring],[Dec;OffDec],[Mask;OffMask],Problem.N,length(Population));

                iter = Problem.FE/Problem.maxFE;
                if iter > 0.6
                    interval = (interval + abs(dis-predis))/2;
                    cnt_interval = cnt_interval + double(abs(interval-abs(dis-predis)) < 1e-2);
                    if cnt_interval >= ceil(sqrt(rho))
                        theta = Problem.N;
                    end
                end
                rho = rho + abs(dis-predis);

            end
        end
    end
end