function [HCluster, CluID, HSort, Clusterspecies] = TwoTrackAggr(data_ev, data_particle, CluID,appear_part,maximalval)
%TWOTRACKAGGR Summary of this function goes here
%   Detailed explanation goes here
%% IF LOOP in case the particle ID is not valid
% If the particles were cut out by the FALSE order and therefore are not
% used in the evaluation, they can not be evaluated via this algorithm. 
[sz,~] = size(data_particle);
HSort = [];
HCluster = [];
Clusterspecies = [];
if sz ~= 0
%% RANGESEARCH for last position of the particle
    % Check last position, time and info of observed particle
    obspos = data_particle(sz,5:6);
    obstime = data_particle(sz,4);
    obspart = data_particle(sz,:);
    % Determination of max particle movement
    t = 30; % in frames 
    t_length_max = sqrt(4*data_particle(1,3)*(t/30)*log(1/(1-0.95))); % MaxMovement of particle in nm to a certain percentage
    r_max = ceil (t_length_max/166); % MaxMovement of particle in pixel
    % Determine Matrix with first appearence of every particle
    AllFades = data_ev(find(diff([data_ev(:,1);max(data_ev(:,1))+1])),:);
    % Which particles are in time window of observed particle?
    TimeMatches = AllFades((AllFades(:,4) > obstime & AllFades(:,4) < obstime+t),:);
    % KDTreeSearch for every particle in r_max:
    MdlKDT_Time = KDTreeSearcher(TimeMatches(:,5:6));
    IdxKDT_r_max = rangesearch(MdlKDT_Time,obspos,r_max);
    lines = sort(IdxKDT_r_max {1,1});
    XYMatches = TimeMatches(sort(lines(1:end)),:);    
    % Cluster out of fading particles, these would be the fading species
    FirstCluster = [obspart;XYMatches(:,:)];
    [fadsz,~] = size (XYMatches);
    [FCsz,~] = size (FirstCluster);
    FdSort = {};
    % If at least two particles are fading in a given time in the area,
    % then the algorithm searches for appearing particles which are bigger
    % and appear in the same area after a given time (same time as in the
    % calculations for the fading particles)
    if FCsz >= 2
        % Searching for particles that appear in the area starting at the
        % observed particle fading +1 frame until the last fading particle
        % frame +30 frames, to cover also this particle in the calculations
        AppearMatches = appear_part((appear_part(:,4) > obstime+1 & appear_part(:,4) < max(FirstCluster(:,4))+t),:);
        SizeMatches = AppearMatches(AppearMatches(:,2) > max(FirstCluster(:,2)),:);
        % KDTreeSearch for every particle in r_max:
        MdlKDT_Appear = KDTreeSearcher(SizeMatches(:,(5:6)));
        IdxKDT_Appear = rangesearch(MdlKDT_Appear,obspos,r_max);
        Appearing = (sort(IdxKDT_Appear {1,1}))';
        AppearList = SizeMatches(Appearing(1:end),:); % List of particle with fitting size and the time frame.
        [appsz,~] = size(AppearList); % how many particles are found
        AppSort = {};
        % If there were particles found during delta t in the searched area,
        % then the algorithm combines both particle lists, generates a sort of
        % the particles (fading or appearing), generates a ClusterID for the
        % particle ensamble and searches via SpeciesDetermine for the
        % corresponding particle species.
        if appsz >=1
            % Now check if particles found fading are fading AFTER the new
            % particle is appearing. (Which is not possible)
            FirstCluster = FirstCluster(FirstCluster(:,4) <= min(AppearList(:,4)),:);
            [FCsz,~] = size (FirstCluster);
            if FCsz >= 2
                CID = [];
                FdSort (1:FCsz,1) =  {'Fading'}; %Define Fading particle list
                AppSort(1:appsz,1) = {'Appearing'}; %Define Appearing particle list
                Sort = [FdSort;AppSort];
                Cluster = [FirstCluster;AppearList]; % Combine Appearing and fading particles
                [Clustersz,~] = size(Cluster);
                CID(1:Clustersz,1) = CluID; % Give the particles a ClusterID
                Cluster = [CID,Cluster]; % Combine CLusterID with Cluster
                FirstCluster = Cluster(1:FCsz,:);
                AppearList = Cluster(FCsz+1:end,:);
                HCluster = [HCluster,Cluster]; % New variable for further calculations
                HSort = [HSort,Sort];
                CluID = CluID+1;
                [Clusterspecies,Code,SpeciesCode] = SpeciesDetermine(data_ev,HCluster,maximalval);
                % Now check if particles belong to different species, for aggregation the
                % species has to be bigger than the later species. If Something is unknown
                % this will be kicked out directly in this step.
                
                % Check the first particle, if this one is unknown kick
                % everything:
                ObsparticleSpecies = cell2mat(SpeciesCode(1,1));
                if ObsparticleSpecies ~= 8
                    Check4Cell = iscell(Code);
                    if Check4Cell == 0
                        NewParticles = SpeciesCode(FCsz+1:end,1);
                        NewClusterspecies = Clusterspecies(FCsz+1:end,1);
                        NewPartCheck = NewParticles(1:end,1) < 8;
                        NewParticleAmount = sum(NewPartCheck);
                        if NewParticleAmount >=1
                            NewParticles = NewParticles(NewPartCheck);
                            NewSort = AppSort(NewPartCheck,1);
                            NewCluster = AppearList(NewPartCheck,:);
                            NewClusterspecies = NewClusterspecies(NewPartCheck,1);
                            % Check for Unknown particles, if so kick them out:
                            StartParticles = SpeciesCode(1:FCsz,1);
                            StartClusterspecies = Clusterspecies(1:FCsz,1);
                            UnknownTest = StartParticles(1:end,1) ~= 8;
                            KnownParticleAmount = sum(UnknownTest);
                            KnownSpecParticles = Code(1:FCsz,1);
                            KnownParticles = KnownSpecParticles(UnknownTest);
                            KnownSort = FdSort(UnknownTest);
                            KnownCluster = FirstCluster(UnknownTest,:);
                            KnownSpecies = StartParticles(UnknownTest);
                            KnownClusterspecies = StartClusterspecies(UnknownTest);
                            % Check for smaller particles, if so kick them out:
                            SpeciesCheck = KnownParticles(1:end,1) >= min(KnownParticles(1:KnownParticleAmount,1));
                            SpeciesAmount = sum(SpeciesCheck);
                            SpeciesParticles = KnownParticles(SpeciesCheck);
                            SpeciesSort = KnownSort(SpeciesCheck);
                            SpeciesCluster = KnownCluster(SpeciesCheck,:);
                            SpeciesSpecies = KnownSpecies(SpeciesCheck);
                            SpeciesClusterspecies = KnownClusterspecies(SpeciesCheck);
                            if SpeciesAmount >= 2 % if there are still two vanishing particles:
                                % Check for acceptable aggregation. This means
                                % the algorithm checks aht results of the
                                % aggregation of the first particle with every
                                % other particle and sees if the resulting
                                % species are found in the appearing particles
                                AggregationSpecies = zeros(SpeciesAmount-1,1);
                                for xx = 2:SpeciesAmount
                                    AggregationSpecies(xx-1,1) = SpeciesSpecies(1,1)+SpeciesSpecies(xx,1);
                                end
                                AggregationParticles = zeros(NewParticleAmount,1);
                                for xx = 1:NewParticleAmount
                                    AggregationCount = 0;
                                    for yy = 1:length(AggregationSpecies)
                                        AggregationTest = NewParticles(xx,1) == AggregationSpecies(yy,1);
                                        AggregationCount = AggregationTest+AggregationCount;
                                    end
                                    if AggregationCount~=0
                                        AggregationParticles(xx,1) = 1;
                                    else
                                        AggregationParticles(xx,1) = 0;
                                    end
                                    AggregationParticles = logical(AggregationParticles);
                                end
                                NewParticles = NewParticles(AggregationParticles,1);
                                NewSort = NewSort(AggregationParticles,1);
                                NewCluster = NewCluster(AggregationParticles,:);
                                NewClusterspecies = NewClusterspecies(AggregationParticles,1);
                                if length(NewParticles) >=1
                                    % Put together the whole results:
                                    Clusterspecies = [SpeciesClusterspecies;NewClusterspecies];
                                    HCluster = [SpeciesCluster(:,:);NewCluster(:,:)];
                                    HSort = [SpeciesSort;NewSort];
                                else 
                                    Clusterspecies = [];
                                    HCluster = [];
                                    HSort = [];
                                    CluID = CluID-1;
                                end
                            else 
                                Clusterspecies = [];
                                HCluster = [];
                                HSort = [];
                                CluID = CluID-1;
                            end
                        else 
                            Clusterspecies = [];
                            HCluster = [];
                            HSort = [];
                            CluID = CluID-1;
                        end
                    else
                        MaxSpec = zeros(length(Code),1);
                        OverallSpec = zeros(length(SpeciesCode),1);
                        for xx = 1:length(Code)
                            CheckRow = isnumeric(cell2mat(Code(xx,1)));
                            CheckRow2 = isnumeric(cell2mat(SpeciesCode(xx,1)));
                            if CheckRow == 0
                                Convert = cell2mat(Code(xx,1));
                                Convert2 = cell2mat(SpeciesCode(xx,1));
                                LengthConvert = length(Convert);
                                LengthConvert2 = length(Convert2);
                                if LengthConvert == 6
                                    MaxSpec(xx,1) = min([str2num(Convert(1,1)),str2num(Convert(1,6))]);
                                else 
                                    MaxSpec(xx,1) = min([str2num(Convert(1:3)),str2num(Convert(7:end))]);
                                end
                                OverallSpec (xx,1) = min([str2num(Convert(1,1)),str2num(Convert(1,6))]);
                            else
                                MaxSpec(xx,1) = cell2mat(Code(xx,1));
                                OverallSpec(xx,1) = cell2mat(SpeciesCode(xx,1)); 
                            end
                        end
                        NewParticles = OverallSpec(FCsz+1:end,1);
                        NewClusterspecies = Clusterspecies(FCsz+1:end,1);
                        NewPartCheck = NewParticles(1:end,1) < 8;
                        NewParticleAmount = sum(NewPartCheck);
                        if NewParticleAmount >=1
                            NewParticles = NewParticles(NewPartCheck);
                            NewSort = AppSort(NewPartCheck,1);
                            NewCluster = AppearList(NewPartCheck,:);
                            NewClusterspecies = NewClusterspecies(NewPartCheck,1);
                            % Check for Unknown particles, if so kick them out:
                            StartParticles = OverallSpec(1:FCsz,1);
                            StartClusterspecies = Clusterspecies(1:FCsz,1);
                            UnknownTest = StartParticles(1:end,1) ~= 8;
                            KnownParticleAmount = sum(UnknownTest);
                            KnownSpecParticles = MaxSpec(1:FCsz,1);
                            KnownParticles = KnownSpecParticles(UnknownTest);
                            KnownSort = FdSort(UnknownTest);
                            KnownCluster = FirstCluster(UnknownTest,:);
                            KnownSpecies = StartParticles(UnknownTest);
                            KnownClusterspecies = StartClusterspecies(UnknownTest);
                            % Check for smaller particles, if so kick them out:
                            SpeciesCheck = KnownParticles(1:end,1) >= min(KnownParticles(1:KnownParticleAmount,1));
                            SpeciesAmount = sum(SpeciesCheck);
                            SpeciesParticles = KnownParticles(SpeciesCheck);
                            SpeciesSort = KnownSort(SpeciesCheck);
                            SpeciesCluster = KnownCluster(SpeciesCheck,:);
                            SpeciesSpecies = KnownSpecies(SpeciesCheck);
                            SpeciesClusterspecies = KnownClusterspecies(SpeciesCheck);
                            if SpeciesAmount >= 2 % if there are still two vanishing particles:
                                % Check for acceptable aggregation. This means
                                % the algorithm checks aht results of the
                                % aggregation of the first particle with every
                                % other particle and sees if the resulting
                                % species are found in the appearing particles
                                AggregationSpecies = zeros(SpeciesAmount-1,1);
                                for xx = 2:SpeciesAmount
                                    AggregationSpecies(xx-1,1) = SpeciesSpecies(1,1)+SpeciesSpecies(xx,1);
                                end
                                AggregationParticles = zeros(NewParticleAmount,1);
                                for xx = 1:NewParticleAmount
                                    AggregationCount = 0;
                                    for yy = 1:length(AggregationSpecies)
                                        AggregationTest = NewParticles(xx,1) == AggregationSpecies(yy,1);
                                        AggregationCount = AggregationTest+AggregationCount;
                                    end
                                    if AggregationCount~=0
                                        AggregationParticles(xx,1) = 1;
                                    else
                                        AggregationParticles(xx,1) = 0;
                                    end
                                    AggregationParticles = logical(AggregationParticles);
                                end
                                NewParticles = NewParticles(AggregationParticles,1);
                                NewSort = NewSort(AggregationParticles,1);
                                NewCluster = NewCluster(AggregationParticles,:);
                                NewClusterspecies = NewClusterspecies(AggregationParticles,1);
                                if length(NewParticles) >=1
                                    % Put together the whole results:
                                    Clusterspecies = [SpeciesClusterspecies;NewClusterspecies];
                                    HCluster = [SpeciesCluster(:,:);NewCluster(:,:)];
                                    HSort = [SpeciesSort;NewSort];
                                else 
                                    Clusterspecies = [];
                                    HCluster = [];
                                    HSort = [];
                                    CluID = CluID-1;
                                end
                            else 
                                Clusterspecies = [];
                                HCluster = [];
                                HSort = [];
                                CluID = CluID-1;
                            end
                        else 
                            Clusterspecies = [];
                            HCluster = [];
                            HSort = [];
                            CluID = CluID-1;
                        end
                    end
                else 
                    Clusterspecies = [];
                    HCluster = [];
                    HSort = [];
                    CluID = CluID-1;
                end
            end
        end
    end
end

end
