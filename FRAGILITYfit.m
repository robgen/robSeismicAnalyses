function [ fragMean, fragStDev, fragMedian, fragMeanMSAS, fragStDevMSAS, ...
    fragMedianMSAS, powerLaw, rawLSparameters, powerLawMSAS ] = FRAGILITYfit( ...
    EDP, IM, DSthresholds, varargin )
%FRAGILITYfit fits fragility curves based on a multiple stripe or a cloud 
%   analysis (including MSAS sequences)
%   
%   The available fitting methods are:
%   Least squares based on power law
%   Maximum likelihood based on binomial PDF
%   Logistic regression based on logit function
%
%   NOTE: conditional fragility for MSAS sequences can be done only using
%   the energy as EDP

%% Example

% % Mainshock only
%
% EDP = [inf;1.03885714285714;0.253742857142857;0.157697142857143;0.112200000000000;0.180657142857143;inf;inf;1.21914285714286;inf;0.712285714285714;inf;0.385714285714286;inf;0.341028571428571;0.610857142857143;1.35542857142857;0.518000000000000;0.288885714285714;0.238542857142857;0.232628571428571;1.00114285714286;0.873142857142857;inf;0.366285714285714;inf;0.256400000000000;inf;0.121942857142857;0.206771428571429;1.41657142857143;0.538571428571429;0.209628571428571;0.759142857142857;inf;inf;inf;inf;0.444000000000000;0.509428571428571;0.269714285714286;inf;0.704000000000000;0.439885714285714;0.907428571428571;inf;0.283314285714286;0.204514285714286;0.499428571428571;inf;0.716000000000000;inf;0.294314285714286;0.876857142857143;1.10057142857143;0.823142857142857;1.94428571428571;0.313142857142857;0.750571428571429;1.19714285714286;0.272000000000000;0.641428571428572;0.178742857142857;0.144000000000000;0.353142857142857;1.19342857142857;inf;inf;0.0574571428571429;0.765142857142857;0.0690285714285714;0.232714285714286;0.0768571428571429;0.0976428571428571;0.0770457142857143;0.360571428571429;0.354800000000000;inf;0.186200000000000;0.269057142857143;inf;inf;inf;0.491142857142857;1.93828571428571;1.62485714285714;0.372571428571429;0.852000000000000;0.353428571428572;0.833714285714286;0.451714285714286;2.39200000000000;1.44200000000000;0.578857142857143;inf;0.599428571428571;0.802571428571429;1.07228571428571;0.742571428571429;inf;0.481428571428571;0.266028571428571;0.397142857142857;inf;0.254865142857143;0.479142857142857;0.444342857142857;inf;0.694285714285714;0.218257142857143;0.245885714285714;0.370000000000000;1.26714285714286;0.659142857142857;0.503428571428571;0.241200000000000;inf;0.161714285714286;0.358857142857143;1.15685714285714;0.506857142857143;0.611428571428571;inf;0.626285714285714;0.562857142857143;1.83857142857143;inf;0.619142857142857;1.58314285714286;1.91028571428571;1.03857142857143;2.40142857142857;inf;inf;inf;inf;0.611714285714286;0.362000000000000;1.27200000000000;1.59285714285714;2.06457142857143;inf;0.813714285714286;inf;inf;0.692857142857143;0.388571428571429;inf;inf;0.196228571428571];
% IM = [0.486652799930677;0.787319437220646;0.250184099224849;0.136286162442021;0.0979323007935992;0.160409447367736;2.51540077826078;0.641280466625313;0.637928524155712;0.755815662786790;0.373750689615324;0.656488077837695;0.306359259410118;1.08446197521364;0.295354026217986;0.420480455942352;0.531361773276019;0.446562212030412;0.234873290715948;0.159231661829126;0.262337956561675;0.623264584144180;0.497859780059622;0.594799232592503;0.307187117277856;0.674039384654122;0.202491097062811;1.17088898366469;0.0870416759355459;0.130884479952337;0.725952274604069;0.428752652679897;0.207212135212415;0.473540779526088;0.959756392941973;0.858601163618842;1.29194429518169;0.680262457334081;0.373969038805093;0.490942536823474;0.283426468829384;0.372047502013720;0.458676380136675;0.317672773258458;0.466432577564956;0.512864556011263;0.193179405767884;0.175321156684480;0.383402167697341;0.996313986942638;0.358128874361397;0.627415141929770;0.205913788484132;0.580583408793357;0.607669214306288;0.501937287093263;0.820529919517645;0.248726555088684;0.580795779577008;0.616388773076326;0.162746199288536;0.362516860491891;0.144179453194560;0.117667232173954;0.194641562398843;0.556958186846883;0.753174145235171;0.679886041138971;0.0304221904683026;0.368415133726386;0.0610002894860024;0.181041924098290;0.0869496750526241;0.0727065611639131;0.0592407835616732;0.269211188810316;0.423447649792432;1.44879961444343;0.184757738208601;0.201180424394397;1.41715835869259;1.62650573314004;1.16962362362204;0.363396816797663;0.589878841558655;0.549452070750000;0.262038405825126;0.440460456181228;0.300310087965029;0.644685826279394;0.425215276237177;0.779528423185899;0.557367699758368;0.366823700596496;1.80090388460768;0.368410044606816;0.724665712071503;0.503240027314005;0.489751228761089;0.909794943176603;0.283884140849041;0.228788520857757;0.301134360204428;0.535553528376702;0.172104207646056;0.377005658461810;0.334850593444772;1.04331379931016;0.428740595312989;0.158427677052459;0.186719825957502;0.301462926105107;0.559691882815162;0.472795566665489;0.435704386387820;0.217533369087830;0.611275340519240;0.151425997098570;0.293694022623420;0.601167795995101;0.351178784718090;0.429899524299582;0.550550373200590;0.347922355758996;0.370495150121426;0.830054231304741;0.519442554194037;0.527931446162056;0.565281367768579;0.747972874097310;0.681771060656725;0.610743157681256;0.867652351627633;0.908672784262972;0.623457536334091;0.564060317316916;0.456471741177065;0.451313308272209;0.612900152031198;0.550086081977074;0.601400915538297;0.402314767498580;0.635933662754627;1.41393776106299;2.06229597053283;0.442588606508941;0.375373441622377;0.763525869406114;0.937870122831604;0.178417536321610];
% DS_thr = [0.25;0.6;1.25;1.67];
% FRAGILITYfit( EDP, IM, DS_thr,'plot', 0, 0, 0, 'ls' );
% FRAGILITYfit( EDP, IM, DS_thr,'plot', 0, 0, 0, 'ml' );
% FRAGILITYfit( EDP, IM, DS_thr,'plot', 0, 0, 0, 'probit' );

% % MSAS sequences
%
% EDP = [inf;1.03885714285714;0.253742857142857;0.157697142857143;0.112200000000000;0.180657142857143;inf;inf;1.21914285714286;inf;0.712285714285714;inf;0.385714285714286;inf;0.341028571428571;0.610857142857143;1.35542857142857;0.518000000000000;0.288885714285714;0.238542857142857;0.232628571428571;1.00114285714286;0.873142857142857;inf;0.366285714285714;inf;0.256400000000000;inf;0.121942857142857;0.206771428571429;1.41657142857143;0.538571428571429;0.209628571428571;0.759142857142857;inf;inf;inf;inf;0.444000000000000;0.509428571428571;0.269714285714286;inf;0.704000000000000;0.439885714285714;0.907428571428571;inf;0.283314285714286;0.204514285714286;0.499428571428571;inf;0.716000000000000;inf;0.294314285714286;0.876857142857143;1.10057142857143;0.823142857142857;1.94428571428571;0.313142857142857;0.750571428571429;1.19714285714286;0.272000000000000;0.641428571428572;0.178742857142857;0.144000000000000;0.353142857142857;1.19342857142857;inf;inf;0.0574571428571429;0.765142857142857;0.0690285714285714;0.232714285714286;0.0768571428571429;0.0976428571428571;0.0770457142857143;0.360571428571429;0.354800000000000;inf;0.186200000000000;0.269057142857143;inf;inf;inf;0.491142857142857;1.93828571428571;1.62485714285714;0.372571428571429;0.852000000000000;0.353428571428572;0.833714285714286;0.451714285714286;2.39200000000000;1.44200000000000;0.578857142857143;inf;0.599428571428571;0.802571428571429;1.07228571428571;0.742571428571429;inf;0.481428571428571;0.266028571428571;0.397142857142857;inf;0.254865142857143;0.479142857142857;0.444342857142857;inf;0.694285714285714;0.218257142857143;0.245885714285714;0.370000000000000;1.26714285714286;0.659142857142857;0.503428571428571;0.241200000000000;inf;0.161714285714286;0.358857142857143;1.15685714285714;0.506857142857143;0.611428571428571;inf;0.626285714285714;0.562857142857143;1.83857142857143;inf;0.619142857142857;1.58314285714286;1.91028571428571;1.03857142857143;2.40142857142857;inf;inf;inf;inf;0.611714285714286;0.362000000000000;1.27200000000000;1.59285714285714;2.06457142857143;inf;0.813714285714286;inf;inf;0.692857142857143;0.388571428571429;inf;inf;0.196228571428571];
% IM = [0.486652799930677;0.787319437220646;0.250184099224849;0.136286162442021;0.0979323007935992;0.160409447367736;2.51540077826078;0.641280466625313;0.637928524155712;0.755815662786790;0.373750689615324;0.656488077837695;0.306359259410118;1.08446197521364;0.295354026217986;0.420480455942352;0.531361773276019;0.446562212030412;0.234873290715948;0.159231661829126;0.262337956561675;0.623264584144180;0.497859780059622;0.594799232592503;0.307187117277856;0.674039384654122;0.202491097062811;1.17088898366469;0.0870416759355459;0.130884479952337;0.725952274604069;0.428752652679897;0.207212135212415;0.473540779526088;0.959756392941973;0.858601163618842;1.29194429518169;0.680262457334081;0.373969038805093;0.490942536823474;0.283426468829384;0.372047502013720;0.458676380136675;0.317672773258458;0.466432577564956;0.512864556011263;0.193179405767884;0.175321156684480;0.383402167697341;0.996313986942638;0.358128874361397;0.627415141929770;0.205913788484132;0.580583408793357;0.607669214306288;0.501937287093263;0.820529919517645;0.248726555088684;0.580795779577008;0.616388773076326;0.162746199288536;0.362516860491891;0.144179453194560;0.117667232173954;0.194641562398843;0.556958186846883;0.753174145235171;0.679886041138971;0.0304221904683026;0.368415133726386;0.0610002894860024;0.181041924098290;0.0869496750526241;0.0727065611639131;0.0592407835616732;0.269211188810316;0.423447649792432;1.44879961444343;0.184757738208601;0.201180424394397;1.41715835869259;1.62650573314004;1.16962362362204;0.363396816797663;0.589878841558655;0.549452070750000;0.262038405825126;0.440460456181228;0.300310087965029;0.644685826279394;0.425215276237177;0.779528423185899;0.557367699758368;0.366823700596496;1.80090388460768;0.368410044606816;0.724665712071503;0.503240027314005;0.489751228761089;0.909794943176603;0.283884140849041;0.228788520857757;0.301134360204428;0.535553528376702;0.172104207646056;0.377005658461810;0.334850593444772;1.04331379931016;0.428740595312989;0.158427677052459;0.186719825957502;0.301462926105107;0.559691882815162;0.472795566665489;0.435704386387820;0.217533369087830;0.611275340519240;0.151425997098570;0.293694022623420;0.601167795995101;0.351178784718090;0.429899524299582;0.550550373200590;0.347922355758996;0.370495150121426;0.830054231304741;0.519442554194037;0.527931446162056;0.565281367768579;0.747972874097310;0.681771060656725;0.610743157681256;0.867652351627633;0.908672784262972;0.623457536334091;0.564060317316916;0.456471741177065;0.451313308272209;0.612900152031198;0.550086081977074;0.601400915538297;0.402314767498580;0.635933662754627;1.41393776106299;2.06229597053283;0.442588606508941;0.375373441622377;0.763525869406114;0.937870122831604;0.178417536321610];
% EDPas = EDP/2;
% IMas = IM /2;
% DS_thr = [0.5;1.0;1.5;2.0];
% FRAGILITYfit( EDP, IM, DS_thr,'plot', 0, EDPas, IMas, 'ls' );
% FRAGILITYfit( EDP, IM, DS_thr,'plot', 0, EDPas, IMas, 'ml' );
% FRAGILITYfit( EDP, IM, DS_thr,'plot', 0, EDPas, IMas, 'probit' );

%% Optional input

% max number of inputs
NmaxVarargin = 8;
numvarargs = length(varargin);
if numvarargs > NmaxVarargin
    error( 'myfuns:somefun2Alt:TooManyInputs', ...
        sprintf('requires at most %d optional inputs',NmaxVarargin) );
end

% set defaults for optional inputs
optargs = {'plot' 0 0 0 'ls' inf 3 0};

% now put these defaults into the valuesToUse cell array, 
% and overwrite with the ones specified in varargin.
optargs(1:numvarargs) = varargin;
% or ...
% [optargs{1:numvarargs}] = varargin{:};

% Place optional args in memorable variable names
[ plotter, colours, EDPas, IMas, fitMethod, collapseEDPthr, minGMsForFit, yieldEDP ] = optargs{:};

%% Manipulate input

entriesLegendMS = {'DS_{1} Slight damage', ...
                   'DS_{2} Moderate damage', ...
                   'DS_{3} Extensive damage', ...
                   'DS_{4} Complete damage'};
% colours for DS1 to DSn in MS-only condition
GrayScale = fliplr(linspace(0, 0.95, numel(DSthresholds)+1));
if size(colours,2) == 1
    colours = cell(numel(DSthresholds),1);
    for i = 1 : numel(DSthresholds)
        colours{i} = [1 1 1] * GrayScale(i+1);
    end
end

% colours for DS1 to DSn in MS-AS conditions
coloursMSAS = [ 0.0 0.00 1 
                0.0 0.75 0
                1.0 0.50 0
                1.0 0.00 0];

% line types corresponding to the mainshock
% if there are more than 4 DSs, it won't work
lineMSAS    = {'-o', '--x', ':p', '-.^', '-.s'};

%% Fit fragility functions (mainshock only, if MSAS are considered)

switch fitMethod
    case 'ls'
        [ fragMean, fragStDev, fragMedian, powerLaw, rawLSparameters ] = ...
            FRAGILITYls( EDP, IM, DSthresholds, collapseEDPthr, 0, yieldEDP );
    case 'ml'
        [ fragMean, fragStDev, fragMedian ] = FRAGILITYml( EDP, IM, DSthresholds );
        powerLaw = [0 0];
        rawLSparameters = NaN;
    case 'probit'
        [ fragMean, fragStDev, fragMedian ] = FRAGILITYprobit( EDP, IM, DSthresholds );
        powerLaw = [0 0];
        rawLSparameters = NaN;
end

%% Fit conditional fragility functions (if MSAS are considered)

if numel(EDPas) ~= 1 % MSAS sequences
    
    warning(sprintf('MSAS sequences are fit assuming that a monotonic EDP is used. \n Therefore EDP(MS)+EDP(AS) is used to fit aftershock fragilities. \n if you want to use drift, set EDP(AS)=driftAS-driftMS in the input')) %#ok<SPWRN>
    
    % determine damage state after MS and AS (energy-based)
    DS = zeros(numel(EDP), 1);
    for i = 1 : numel(DSthresholds)
        DS(EDP >= DSthresholds(i)) = i;
    end
    
    DSas = zeros(numel(EDP), 1);
    for i = 1 : numel(DSthresholds)
        DSas(EDP+EDPas >= DSthresholds(i)) = i;
    end
    
    % preallocation
    NconditionalDS  = ( numel(DSthresholds)+1 ) * numel(DSthresholds) - sum( 1 : numel(DSthresholds)-1 );
    fragMeanMSAS    = zeros(NconditionalDS, 1);
    fragStDevMSAS   = zeros(NconditionalDS, 1);
    fragMedianMSAS  = zeros(NconditionalDS, 1);
    
    % names of the fragilities (and setting: DSk | DSj )
    nameFrag = cell(1,NconditionalDS);
    settingFrag = zeros(NconditionalDS,2);
    frag = 0;
    for dsAS = 1 : numel(DSthresholds)
        frag = frag + 1;
        nameFrag{frag} = sprintf('DS_{%d}|DS_{0}^{MS}', dsAS);
        settingFrag(frag,:) = [ dsAS 0 ];
    end
    
    for dsMS = 1 : numel(DSthresholds)
        for dsAS = dsMS : numel(DSthresholds)
            frag = frag + 1;
            nameFrag{frag} = sprintf('DS_{%d}|DS_{%d}^{MS}', dsAS, dsMS);
            settingFrag(frag,:) = [ dsAS dsMS ];
        end
    end
    
    switch fitMethod
        case 'ls'
            % fragilities conditional to DS0 in the Mainshock (ONLY ENERGY-BASED)
            ds = 0;
            powerLawMSAS = zeros(numel(DSthresholds)+1,3);
            for dsAS = 1 : numel(DSthresholds)
                ds = ds + 1;
                if sum(DS==0) >= minGMsForFit
                    upwardTransl = 0; %mean(EDP(DS==0));
                    [ fragMeanMSAS(ds,1), fragStDevMSAS(ds,1), ...
                        fragMedianMSAS(ds,1), powerLawMSAS(1,:), ...
                        ~, psdmMSAS(1)] = FRAGILITYls( ...
                        EDP(DS==0)+EDPas(DS==0), IMas(DS==0), ...
                        DSthresholds(dsAS), collapseEDPthr, upwardTransl, yieldEDP );
                else
                    fragMeanMSAS(ds,1) = NaN;
                    fragStDevMSAS(ds,1) = NaN;
                    fragMedianMSAS(ds,1) = NaN;
                    powerLawMSAS(1,:) = [NaN NaN];
                end
            end
            
            % fragilities conditional to DSi in the Mainshock (ONLY ENERGY-BASED)
            for dsMS = 1 : numel(DSthresholds)
                for dsAS = dsMS : numel(DSthresholds)
                    ds = ds + 1;
                    if sum(DS==dsMS) >= minGMsForFit
                        upwardTransl = DSthresholds(dsMS); % mean(EDP(DS==dsMS));
                        [ fragMeanMSAS(ds,1), fragStDevMSAS(ds,1), ...
                            fragMedianMSAS(ds,1), powerLawMSAS(dsMS+1,:),...
                            ~, psdmMSAS(dsMS+1)] = FRAGILITYls( ...
                            EDP(DS==dsMS)+EDPas(DS==dsMS), ...
                            IMas(DS==dsMS), DSthresholds(dsAS), ...
                            collapseEDPthr, upwardTransl, yieldEDP );
                    else
                        fragMeanMSAS(ds,1) = NaN;
                        fragStDevMSAS(ds,1) = NaN;
                        fragMedianMSAS(ds,1) = NaN;
                        powerLawMSAS(dsMS+1,:) = [NaN NaN];
                    end
                end
            end
            
            % control plot
            cols = [0.9 0.9 0.9; coloursMSAS];
            figure; hold on
            for dsMS = 0 : numel(DSthresholds)
                scatter(IMas(DS==dsMS), EDP(DS==dsMS)+EDPas(DS==dsMS), ...
                    20, cols(dsMS+1,:), 'filled', 'MarkerEdgeColor', 0.5*[1 1 1])
            end
            
            IMplot = 0 : 0.005 : 1;
            for dsMS = 0 : numel(DSthresholds)
                plot(IMplot, psdmMSAS(dsMS+1).handle(IMplot), ...
                    'LineWidth', 2, 'Color', cols(dsMS+1,:))
            end
            legend({'DS0', 'DS1', 'DS2', 'DS3', 'DS4'})
            xlabel('IM(AS)')
            ylabel('Eh(MS+AS)')
            set(gca, 'FontSize', 18)
            a = 0;
        case 'ml'
            
            % fragilities conditional to DS0 in the Mainshock (ONLY ENERGY-BASED)
            ds = 0;
            for dsAS = 1 : numel(DSthresholds)
                ds = ds + 1;
                if sum(DS==0) >= minGMsForFit
                    [ fragMeanMSAS(ds,1), fragStDevMSAS(ds,1), ...
                        fragMedianMSAS(ds,1) ] = FRAGILITYml( ...
                        EDP(DS==0)+EDPas(DS==0), IMas(DS==0), DSthresholds(dsAS) );
                else
                    fragMeanMSAS(ds,1) = NaN;
                    fragStDevMSAS(ds,1) = NaN;
                    fragMedianMSAS(ds,1) = NaN;
                end
            end
            
            % fragilities conditional to DSi in the Mainshock (ONLY ENERGY-BASED)
            for dsMS = 1 : numel(DSthresholds)
                for dsAS = dsMS : numel(DSthresholds)
                    ds = ds + 1;
                    if sum(DS==dsMS) >= minGMsForFit
                        [ fragMeanMSAS(ds,1), fragStDevMSAS(ds,1), ...
                            fragMedianMSAS(ds,1) ] = FRAGILITYml( ...
                            EDP(DS==dsMS)+EDPas(DS==dsMS), IMas(DS==dsMS), ...
                            DSthresholds(dsAS) );
                    else
                        fragMeanMSAS(ds,1) = NaN;
                        fragStDevMSAS(ds,1) = NaN;
                        fragMedianMSAS(ds,1) = NaN;
                    end
                end
            end
            
            powerLawMSAS    = 0;
        case 'probit'
            
            % fragilities conditional to DS0 in the Mainshock (ONLY ENERGY-BASED)
            ds = 0;
            for dsAS = 1 : numel(DSthresholds)
                ds = ds + 1;
                if sum(DS==0) >= minGMsForFit
                    [ fragMeanMSAS(ds,1), fragStDevMSAS(ds,1), ...
                        fragMedianMSAS(ds,1) ] = FRAGILITYprobit( ...
                        EDP(DS==0)+EDPas(DS==0), IMas(DS==0), DSthresholds(dsAS) );
                else
                    fragMeanMSAS(ds,1) = NaN;
                    fragStDevMSAS(ds,1) = NaN;
                    fragMedianMSAS(ds,1) = NaN;
                end
            end
            
            % fragilities conditional to DSi in the Mainshock (ONLY ENERGY-BASED)
            for dsMS = 1 : numel(DSthresholds)
                for dsAS = dsMS : numel(DSthresholds)
                    ds = ds + 1;
                    if sum(DS==dsMS) >= minGMsForFit
                        [ fragMeanMSAS(ds,1), fragStDevMSAS(ds,1), ...
                            fragMedianMSAS(ds,1) ] = FRAGILITYprobit( ...
                            EDP(DS==dsMS)+EDPas(DS==dsMS), IMas(DS==dsMS), ...
                            DSthresholds(dsAS) );
                    else
                        fragMeanMSAS(ds,1) = NaN;
                        fragStDevMSAS(ds,1) = NaN;
                        fragMedianMSAS(ds,1) = NaN;
                    end
                end
            end
            
            powerLawMSAS    = 0;
    end
    
    % calculate how many analyses have the same DS in the mainshock
    dsMSpopulation = zeros(numel(DSthresholds), 1);
    for dsMS = 0 : numel(DSthresholds) + 1
        dsMSpopulation(dsMS+1) = sum(DS==dsMS);
    end
    
else
    
    fragMeanMSAS    = 0;
    fragStDevMSAS   = 0;
    fragMedianMSAS  = 0;
    powerLawMSAS    = 0;
end

%% Plot 

if strcmp(plotter, 'plot')
    
    %%% Plot fragilities for the intact structure
    IMfrag = (0:0.005:2.5)';
    PROBfrag  = zeros(numel(IMfrag), numel(DSthresholds));
    figure('Position', [162   495   560   420]);
    hold on
    for ds = numel(DSthresholds) : -1 : 1
        PROBfrag(:,ds) = logncdf(IMfrag, fragMean(ds,1), fragStDev(ds,1));
        p(ds) = plot(IMfrag, PROBfrag(:,ds), ...
            'LineWidth', 3, 'Color', colours{ds});
    end
    
    annotation(gcf,'textbox',...
        [0.751785714285714 0.157142857142857 0.14285714285714 0.0595238095238098],...
        'String', fitMethod, 'FontSize',16, 'FitBoxToText','off', 'EdgeColor','none');

    xlabel('Avg SA'); %'S_{a}(T_{1}) [g]');
    ylabel('Pr[DS>ds_{i}|IM]');
    set(gca,'FontSize',16);
    legend(p, entriesLegendMS, 'Location', 'northwest', 'box', 'off');
    
    
%     % extra plot
%     scatter(IM, EDP, 36, [1 1 1]*0.7, 'filled', 'MarkerEdgeColor', 'k'); hold on
%     plot( min(IM):0.01:max(IM) , powerLaw(1)*(min(IM):0.01:max(IM)).^powerLaw(2), '-', 'Color', [1 1 1]*0.0, 'LineWidth', 2)
%     for ds = 1:numel(DSthresholds)
%         plot( [min(IM) max(IM)] , [DSthresholds(ds) DSthresholds(ds)] , 'LineWidth', 1, 'Color', colours{ds})
%     end
%     xlabel('IM');
%     ylabel('EDP');
%     set(gca,'FontSize',16);
    
    %%% Plot conditional fragilities
    if numel(EDPas) ~= 1 % MSAS sequences
        
        title('Intact building fragilities', 'FontSize', 16)

        
        figure2 = figure('Position', [162     1   560   420]);
        axes2 = axes('Parent',figure2,'YGrid','on','XGrid','on','FontSize',14);
        box(axes2,'on');
        hold(axes2,'all');
        
        IMfragMSAS      = (0:0.005:2.5)';
        PROBfragMSAS    = zeros(numel(IMfragMSAS), NconditionalDS);
        for frag = 1 : NconditionalDS
            if ~isnan(fragMedianMSAS(frag,1)) && fragMedianMSAS(frag,1)>0
                PROBfragMSAS(:,frag)  = logncdf(IMfragMSAS,...
                    fragMeanMSAS(frag,1),fragStDevMSAS(frag,1));
            end
            plot(IMfragMSAS, PROBfragMSAS(:,frag), lineMSAS{settingFrag(frag,2)+1}, ...
                'LineWidth', 2, 'Color', coloursMSAS(settingFrag(frag,1)+1,:), ...
                'MarkerSize', 12, 'MarkerIndices', 1:20:length(IMfragMSAS))
        end
        xlabel('Avg SA'); %'S_{a}(T_{1}) [g]');
        ylabel('Pr[DS>ds_{i}|IM]');
        set(gca,'FontSize',16);
        
%         % add also the mainshock fragilities
%         for dsMS = 1 : numel(DSthresholds)
%             plot(IMfrag, PROBfrag(:,dsMS), ...
%                 'LineWidth', 4, 'Color', coloursMSAS(dsMS+1,:));
%         end
        legend([nameFrag {'DS1^{MS}' 'DS2^{MS}' 'DS3^{MS}' 'DS4^{MS}'}], ...
            'Location', 'SouthEast');
    end
    
end

end


%% Subfunctions

function [ FragMean, FragStDev, FragMedian, powerLaw, rawLSparameters, psdm ] = ...
    FRAGILITYls( EDP, IM, DSthresholds, collapseEDPthr, upwardTranslationPSDM, yieldEDP )
%FRAGILITYls fits fragility curves with least squares given the results
%   of a cloud analysis or a multiple stripe   
% 
%   see Mai2017, section 2.1.2
%   
%   NOTE: 
%   IM vs EDP is EDP = aIM^b
%   a = exp(p(2))
%   b = p(1)
%
%
%%% Example
%
% EDP = [inf;1.03885714285714;0.253742857142857;0.157697142857143;0.112200000000000;0.180657142857143;inf;inf;1.21914285714286;inf;0.712285714285714;inf;0.385714285714286;inf;0.341028571428571;0.610857142857143;1.35542857142857;0.518000000000000;0.288885714285714;0.238542857142857;0.232628571428571;1.00114285714286;0.873142857142857;inf;0.366285714285714;inf;0.256400000000000;inf;0.121942857142857;0.206771428571429;1.41657142857143;0.538571428571429;0.209628571428571;0.759142857142857;inf;inf;inf;inf;0.444000000000000;0.509428571428571;0.269714285714286;inf;0.704000000000000;0.439885714285714;0.907428571428571;inf;0.283314285714286;0.204514285714286;0.499428571428571;inf;0.716000000000000;inf;0.294314285714286;0.876857142857143;1.10057142857143;0.823142857142857;1.94428571428571;0.313142857142857;0.750571428571429;1.19714285714286;0.272000000000000;0.641428571428572;0.178742857142857;0.144000000000000;0.353142857142857;1.19342857142857;inf;inf;0.0574571428571429;0.765142857142857;0.0690285714285714;0.232714285714286;0.0768571428571429;0.0976428571428571;0.0770457142857143;0.360571428571429;0.354800000000000;inf;0.186200000000000;0.269057142857143;inf;inf;inf;0.491142857142857;1.93828571428571;1.62485714285714;0.372571428571429;0.852000000000000;0.353428571428572;0.833714285714286;0.451714285714286;2.39200000000000;1.44200000000000;0.578857142857143;inf;0.599428571428571;0.802571428571429;1.07228571428571;0.742571428571429;inf;0.481428571428571;0.266028571428571;0.397142857142857;inf;0.254865142857143;0.479142857142857;0.444342857142857;inf;0.694285714285714;0.218257142857143;0.245885714285714;0.370000000000000;1.26714285714286;0.659142857142857;0.503428571428571;0.241200000000000;inf;0.161714285714286;0.358857142857143;1.15685714285714;0.506857142857143;0.611428571428571;inf;0.626285714285714;0.562857142857143;1.83857142857143;inf;0.619142857142857;1.58314285714286;1.91028571428571;1.03857142857143;2.40142857142857;inf;inf;inf;inf;0.611714285714286;0.362000000000000;1.27200000000000;1.59285714285714;2.06457142857143;inf;0.813714285714286;inf;inf;0.692857142857143;0.388571428571429;inf;inf;0.196228571428571];
% IM = [0.486652799930677;0.787319437220646;0.250184099224849;0.136286162442021;0.0979323007935992;0.160409447367736;2.51540077826078;0.641280466625313;0.637928524155712;0.755815662786790;0.373750689615324;0.656488077837695;0.306359259410118;1.08446197521364;0.295354026217986;0.420480455942352;0.531361773276019;0.446562212030412;0.234873290715948;0.159231661829126;0.262337956561675;0.623264584144180;0.497859780059622;0.594799232592503;0.307187117277856;0.674039384654122;0.202491097062811;1.17088898366469;0.0870416759355459;0.130884479952337;0.725952274604069;0.428752652679897;0.207212135212415;0.473540779526088;0.959756392941973;0.858601163618842;1.29194429518169;0.680262457334081;0.373969038805093;0.490942536823474;0.283426468829384;0.372047502013720;0.458676380136675;0.317672773258458;0.466432577564956;0.512864556011263;0.193179405767884;0.175321156684480;0.383402167697341;0.996313986942638;0.358128874361397;0.627415141929770;0.205913788484132;0.580583408793357;0.607669214306288;0.501937287093263;0.820529919517645;0.248726555088684;0.580795779577008;0.616388773076326;0.162746199288536;0.362516860491891;0.144179453194560;0.117667232173954;0.194641562398843;0.556958186846883;0.753174145235171;0.679886041138971;0.0304221904683026;0.368415133726386;0.0610002894860024;0.181041924098290;0.0869496750526241;0.0727065611639131;0.0592407835616732;0.269211188810316;0.423447649792432;1.44879961444343;0.184757738208601;0.201180424394397;1.41715835869259;1.62650573314004;1.16962362362204;0.363396816797663;0.589878841558655;0.549452070750000;0.262038405825126;0.440460456181228;0.300310087965029;0.644685826279394;0.425215276237177;0.779528423185899;0.557367699758368;0.366823700596496;1.80090388460768;0.368410044606816;0.724665712071503;0.503240027314005;0.489751228761089;0.909794943176603;0.283884140849041;0.228788520857757;0.301134360204428;0.535553528376702;0.172104207646056;0.377005658461810;0.334850593444772;1.04331379931016;0.428740595312989;0.158427677052459;0.186719825957502;0.301462926105107;0.559691882815162;0.472795566665489;0.435704386387820;0.217533369087830;0.611275340519240;0.151425997098570;0.293694022623420;0.601167795995101;0.351178784718090;0.429899524299582;0.550550373200590;0.347922355758996;0.370495150121426;0.830054231304741;0.519442554194037;0.527931446162056;0.565281367768579;0.747972874097310;0.681771060656725;0.610743157681256;0.867652351627633;0.908672784262972;0.623457536334091;0.564060317316916;0.456471741177065;0.451313308272209;0.612900152031198;0.550086081977074;0.601400915538297;0.402314767498580;0.635933662754627;1.41393776106299;2.06229597053283;0.442588606508941;0.375373441622377;0.763525869406114;0.937870122831604;0.178417536321610];
% DS_thr = [0.5;1.0;1.5;2.0];
% collapseThr = 10;
% [ fragMean, fragStDev, fragMedian, PowerLawParameters ] = FRAGILITYls( EDP, IM, DS_thr, collapseThr );

%% Optional input

if nargin < 6
    yieldEDP = 0;
end

if nargin < 5
    upwardTranslationPSDM = 0;
end

%% Partition

% collapse cases
collapseFLAG = zeros(size(EDP));
collapseFLAG(isinf(EDP) | EDP>=collapseEDPthr) = 1;

% no-collapse cases
noCollapseIM    = IM(~isinf(EDP) & EDP<collapseEDPthr);
noCollapseEDP   = EDP(~isinf(EDP) & EDP<collapseEDPthr);

%% Non-collapse cases: power law regression

[~, gof, powerLaw] = fitPowerLaw(noCollapseIM, ...
    noCollapseEDP-upwardTranslationPSDM, yieldEDP);
a = powerLaw(1); b = powerLaw(2); t = upwardTranslationPSDM;
powerLaw(3) = t;

% EDP = aIM^b + t; ln(EDP-t) = ln(a) + b*ln(IM)
psdm.handle = @(IM) a*IM.^b + t;
psdm.a = a;
psdm.b = b;
psdm.t = t;
psdm.R2 = gof.rsquare;

if b < 0
    warning('Power Law: b is negative. Fragility parameters are NaN')
    FragMean = nan(numel(DSthresholds),1);
    FragStDev = nan(numel(DSthresholds),1);
    FragMedian = nan(numel(DSthresholds),1);
    
    rawLSparameters.noCollapseFragMean = nan(numel(DSthresholds),1);
    rawLSparameters.noCollapseFragStDev = nan(numel(DSthresholds),1);
    rawLSparameters.noCollapseFragMedian = nan(numel(DSthresholds),1);
    
    rawLSparameters.collapseLogit = NaN;
    return
end

% calculate the error of the log(data) to the line
resid = log(noCollapseEDP) - (log(a) + b.*log(noCollapseIM));

noCollapseFragMean(:,1)   = log( DSthresholds ./ a ) ./ b;
noCollapseFragStDev(:,1)  = std( resid ) / abs(b) .* ones(numel(DSthresholds),1);
noCollapseFragMedian(:,1) = exp( noCollapseFragMean ); % intersection of the power law with the EDP threshold line

% functional form of the no-collapse fragilities
lastPower10 = 0;
Npoints = 1000;
PROBfragNoCollapse  = zeros(Npoints, numel(DSthresholds));
while PROBfragNoCollapse(end,1) ~= 1
    lastPower10 = lastPower10 + 3;
    % recursively define IM values to be big enough to reach Frag=1
    IMfrag = logspace(-5, lastPower10, Npoints)';
    for ds = 1:numel(DSthresholds)
        PROBfragNoCollapse(:,ds) = logncdf(IMfrag, noCollapseFragMean(ds,1), noCollapseFragStDev(ds,1));
    end
end

%% Collapse cases: logistic regression

% logistic regression
if sum(collapseFLAG) >= 3
    logisticParam = glmfit(log(IM), collapseFLAG, 'binomial', 'link', 'logit');
    PROBfragCollapse = 1 ./ ( 1 + exp(-logisticParam(1) -logisticParam(2)*log(IMfrag)) );
else
    logisticParam = NaN;
    PROBfragCollapse = zeros(numel(IMfrag), 1);
end

%% Final fragility

epsFactor = 50; % used to avoid equal points in the interpolation

PROBfrag = zeros(numel(IMfrag), numel(DSthresholds));
percentiles = zeros(numel(DSthresholds),3);
for ds = 1 : numel(DSthresholds)
    PROBfrag(:,ds)      = PROBfragNoCollapse(:,ds).*(1-PROBfragCollapse) + PROBfragCollapse;
    percentiles(ds,:)   = interp1(...
        PROBfrag(PROBfrag(:,ds)<1-epsFactor*eps & PROBfrag(:,ds)>epsFactor*eps, ds), ...
        IMfrag(PROBfrag(:,ds)<1-epsFactor*eps & PROBfrag(:,ds)>epsFactor*eps), ...
        [ 0.16 0.5 0.84 ]);
end

FragMedian(:,1) = percentiles(:,2);
FragStDev(:,1)  = ( log(percentiles(:,3)) - log(percentiles(:,1)) ) / 2;
FragMean(:,1)   = log( FragMedian );

% % control plot
% lognorm = logncdf(IMfrag, FragMean(ds,1), FragStDev(ds,1));
% figure; hold on
% plot(IMfrag, PROBfrag(:,ds))
% plot(IMfrag, lognorm)
% legend('original', 'lognormal simplification')

%% Collect raw parameters (noCollapse lognormals + collapse logit)

rawLSparameters.noCollapseFragMean = noCollapseFragMean;
rawLSparameters.noCollapseFragStDev = noCollapseFragStDev;
rawLSparameters.noCollapseFragMedian = noCollapseFragMedian;

rawLSparameters.collapseLogit = logisticParam;

end


function [fitresult, gof, powerLaw] = fitPowerLaw(noCollapseIM, noCollapseEDP, yieldEDP)

if yieldEDP == 0
    % single regression with all the points
    [xData, yData] = prepareCurveData( log(noCollapseIM), log(noCollapseEDP) );

    ft = fittype( 'poly1' );
    [fitresult, gof] = fit( xData, yData, ft );
    a = exp(fitresult.p2);
    b = fitresult.p1;
else
    % line to interpolate the points up to yielding (get a)
    preYield = noCollapseEDP <= yieldEDP;
    
    [xData, yData] = prepareCurveData( ...
        noCollapseIM(preYield), noCollapseEDP(preYield) );
    ft = fittype( {'x'} );
    fitresult = fit( xData, yData, ft );
    a = fitresult.a;
    
    % power law of all the points, given a (get b)
    [xData, yData] = prepareCurveData( log(noCollapseIM), log(noCollapseEDP)-log(a) );
    ft = fittype( {'x'} );
    [fitresult, gof]  = fit( xData, yData, ft );
    b = fitresult.a;
end

powerLaw = [a b];

end


function [ fragMean, fragStDev, fragMedian ] = FRAGILITYprobit( EDP, IM, DSthresholds )
%FRAGILITYlogit fits fragility curves with probit regression given the
%   results of a multiple stripe analysis (based on Jack Baker's code).
%   Note: a cloud can be seen as a multiple stripe in which each stripe has
%   only one GM


% Define collapses for each stripe (or for the cloud)
stripeIM  = unique(IM);

% Ncollapsed 
NgmTotPerStripe         = zeros(numel(stripeIM),1);
NgmCollapsePerStripe    = zeros(numel(stripeIM), numel(DSthresholds));
for consideredIM = 1 : numel(stripeIM)
    NgmTotPerStripe(consideredIM,1) = sum( IM==stripeIM(consideredIM) );
    for ds = 1 : numel(DSthresholds)
        NgmCollapsePerStripe(consideredIM,ds) = sum( EDP(IM==stripeIM(consideredIM)) >= DSthresholds(ds) );
    end
end

fragMedian      = zeros(numel(DSthresholds),1);
fragStDev       = zeros(numel(DSthresholds),1);
for ds = 1 : numel(DSthresholds)
    % probit regression
    Y = [NgmCollapsePerStripe(:,ds) NgmTotPerStripe]; % vector of number of collapses and number of records at each level
    b = glmfit(log(stripeIM), Y, 'binomial', 'link', 'probit');

    % convert probit coefficients to lognormal distribution parameters
    fragMedian(ds,1)  = exp(-b(1)/b(2));
    fragStDev(ds,1)   = 1/b(2);
end

fragMean    = log(fragMedian);

end


function [ fragMean, fragStDev, fragMedian ] = FRAGILITYml( EDP, IM, DSthresholds )
%FRAGILITYml fits fragility curves with maximum likelihood given the results
%   of a multiple stripe analysis (based on Jack Baker's code). Note: a
%   cloud can be seen as a multiple stripe in which each stripe has only
%   one GM
% 
%
% This function fits a lognormal CDF to observed probability of collapse 
% data using optimization on the likelihood function for the data. 
% These calculations are based on equation 11 of the following paper:
%
% Baker, J. W. (2015). “Efficient analytical fragility function fitting 
% using dynamic structural analysis.” Earthquake Spectra, 31(1), 579-599.
%
%
% INPUTS:
% IM            1xn           IM levels of interest
% NgmTot       1x1 or 1xn    number of ground motions used at each IM level
% NgmCollapse 	1xn           number of collapses observed at each IM level
% 
% OUTPUTS:
% theta         1x1           median of fragility function
% beta          1x1           lognormal standard deviation of fragility function


% Define collapses for each stripe (or for the cloud)
stripeIM  = unique(IM);

% Ncollapsed 
NgmTotPerStripe         = zeros(numel(stripeIM),1);
NgmCollapsePerStripe    = zeros(numel(stripeIM), numel(DSthresholds));
for consideredIM = 1 : numel(stripeIM)
    NgmTotPerStripe(consideredIM,1) = sum( IM==stripeIM(consideredIM) );
    for ds = 1 : numel(DSthresholds)
        NgmCollapsePerStripe(consideredIM,ds) = ...
            sum( EDP(IM==stripeIM(consideredIM)) >= DSthresholds(ds) );
    end
end


% maximisation of the likelihood
guessFragParam  = [0.8 0.4];
options         = optimset('MaxFunEvals',1000, 'MaxIter', 10000, 'GradObj', 'off');
fragMedian      = zeros(numel(DSthresholds),1);
fragStDev       = zeros(numel(DSthresholds),1);
for ds = 1 : numel(DSthresholds)
    fragParam           = fminsearch(@getLogLikelihood, ...
        guessFragParam, options, NgmTotPerStripe, ...
        NgmCollapsePerStripe(:,ds), stripeIM) ;
    
    fragMedian(ds,1)	= fragParam(1);
    fragStDev(ds,1)  	= fragParam(2);
end

fragMean = log(fragMedian);

if ~all(fragMedian > 0)
    warning('median of one fragility is negative. DS exceeded for all the analyses? Too-low EDP threshold?')
end

end


function [logLikelihood] = getLogLikelihood(...
    theta, NgmTotPerStripe, NgmCollapsePerStripe, stripeIM)
%getLogLikelihood calculates the log likelihood, based on binomial PDF, for
%   a multiple stripe analysis. Note: a cloud can be seen as a multiple
%   stripe in which each stripe has only one GM

% don't let median of fragility function go below zero
theta(1) = max(theta(1),0);

% estimated probabilities of collapse, given the current fragility functionparameter estimates
p = normcdf(log(stripeIM), log(theta(1)), theta(2)); 

% likelihood of observing NgmCollapse(i) collapses, given NgmTot
% observations, using the current fragility function parameter estimates
likelihood      = binopdf(NgmCollapsePerStripe', NgmTotPerStripe', p'); % 
%likelihood2     = sym( factorial( sym(NgmTotPerStripe) ) ./ ( factorial( sym(NgmCollapsePerStripe) ) .* factorial( sym(NgmTotPerStripe-NgmCollapsePerStripe) ) ) .* p.^NgmCollapsePerStripe .* (1-p).^(NgmTotPerStripe-NgmCollapsePerStripe) );


% sum negative log likelihood (we take the negative value because we want
% the maximum log likelihood, and the function is searching for a minimum)
logLikelihood   = -sum(log(likelihood));

end

