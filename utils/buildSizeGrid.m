%
% buildSizeGrid.m
% 
% Takes in all varialoptExpes, and allometric parameters and creates parameter
% arrays to store 
%

function [lp, lz, NP, NZ, bgcRates, idxUptake, idxSat, idxGrazing, idxPredPrey, zeta] = ...
    buildSizeGrid(modeltype,MP,MZ,dx,data_dir, uptakeCoeff, uptakeExp, monodCoeff, grazingCoeff, grazingExp, loptCoeff, loptExp, grazingEff)


% identify index values for bgc rates.
idxUptake = 1;
idxSat = 2;
idxGrazing = 3;
idxPredPrey = 4;
idxTot = 4;


%%%%%%%%%%%%%%%%%%%%% Calculate everything here %%%%%%%%%%%%%%%%%%%%%%%%%
switch modeltype
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 0 % no bgc
        % have a placeholder to read in the files
        lp = [0 0]'; 
        lz = [0 0]';
        zeta = 0;
        NP = 2; 
        NZ = 2;
        
        
        % placeholders
        maxUptakeRate = lp;
        nutrientSaturation = lp;
        maxGrazingRate = lz;
        optimalPredPrey = lz;
        
        % copy values over to a matrix
        bgcRates = zeros(NP,idxTot);
        bgcRates(:,idxUptake) = maxUptakeRate;
        bgcRates(:,idxSat) = nutrientSaturation;
        bgcRates(:,idxGrazing) = maxGrazingRate;
        bgcRates(:,idxPredPrey) = optimalPredPrey;
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    case 1 % size structured npzd model

        %%% Create a vector to hold all of the size classes
        minp = 0.5;
        maxz = 300;
        minz = exp(log(minp./loptCoeff)./loptExp); % micrometers
        maxp = floor(loptCoeff*(maxz).^loptExp);

        if MP == 1
            lp = minp.';
            NP = 1;
        else
            lp = (logspace(log10(minp),log10(maxp),MP)).';
            NP = MP;
        end

        if MZ == 1
            lz = minz.';
            NZ = 1;
        else
            lz = (logspace(log10(minz),log10(maxz),MZ)).';
            NZ = MZ;
        end
        
        % calculate allometric rates
        maxUptakeRate = uptakeCoeff.*(lp).^uptakeExp;
        nutrientSaturation = monodCoeff.*(lp);
        maxGrazingRate = grazingCoeff.*(lz).^grazingExp;
        optimalPredPrey = loptCoeff.*(lz).^loptExp;
        
        % store all variables
        NMAX = max(NP,NZ);
        bgcRates = zeros(NMAX,idxTot);
        bgcRates(1:NP,idxUptake) = maxUptakeRate;
        bgcRates(1:NP,idxSat) = nutrientSaturation;
        bgcRates(1:NZ,idxGrazing) = maxGrazingRate;
        bgcRates(1:NZ,idxPredPrey) = optimalPredPrey;
        
        % calculate the density dependent mortality here
        zeta = grazingEff.*(mean(maxGrazingRate).^2)./(12*mean(maxUptakeRate));
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2 % reduced order ecosystem model
        ssemFile = fullfile(data_dir,'ssemData.mat');
        SF = load(ssemFile);
        
        % size grids
        lp = SF.lp;
        lz = SF.lz;
        
        % load in the locations of the peaks
        plocs = SF.plocs;
        zlocs = SF.zlocs;
        
        % load in the number of peaks to expect
        npvec = SF.ppeakvec;
        nzvec = SF.zpeakvec;
        
        % load in the grazing profile width vector
        dxvec = SF.dxvec;
        
        % variables to calculate the total zeta value
        maxUptakeRate = uptakeCoeff.*(lp).^uptakeExp;
        maxGrazingRate = grazingCoeff.*(lz).^grazingExp;
        zeta = grazingEff.*(mean(maxGrazingRate).^2)./(12*mean(maxUptakeRate));
        
        % find the correct reference run
        [~,ind] = min(abs(dxvec - dx));
        
        lp = plocs(ind,1:npvec(ind)); NP = length(lp);
        lz = zlocs(ind,1:nzvec(ind)); NZ = length(lz);
        
        % calculate allometric rates
        rUptakeRate = uptakeCoeff.*(lp).^uptakeExp;
        rnutrientSaturation = monodCoeff.*(lp);
        rGrazingRate = grazingCoeff.*(lz).^grazingExp;
        roptimalPredPrey = loptCoeff.*(lz).^loptExp;
        
        % store all variables
        NMAX = max(NP,NZ);
        bgcRates = zeros(NMAX,idxTot);
        bgcRates(1:NP,idxUptake) = rUptakeRate;
        bgcRates(1:NP,idxSat) = rnutrientSaturation;
        bgcRates(1:NZ,idxGrazing) = rGrazingRate;
        bgcRates(1:NZ,idxPredPrey) = roptimalPredPrey;
        
    otherwise
        disp('buildSizeGrid: No valid modeltype chosen in creating bgc size vectors')
        
end

% convert to the correct index in C
idxUptake = 0;
idxSat = 1;
idxGrazing = 2;
idxPredPrey = 3;

end