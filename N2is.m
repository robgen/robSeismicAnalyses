function [ IS, CAPretPERIOD, PPdisp, PPacc ] = N2is( ...
    ADRSspectrumLS, pushover, DISPmechanism, DISPcapacityLS, ...
    floorMASSES, PHI1, plotter, varargin )
%N2isv Performs the N2 method according to Eurocode and calculates the IS
%   
%   Firstly, the N2 method is applied for a given ADRS spectrum and a
%   pushover curve, providing the performance point on the curve, if this
%   point exists.
%
%   Then, for a given value of displacement corresponding to a limit state,
%   the Safety index is calculated (as per Circolare PAM) by increasing (or
%   decreasing) the given ADRS spectrum.
%
%   NOTE: If the 5% spectrum related to capacity is just found by scaling
%   the original spectrum, the Safety Index (IS) is equal to the NBS. They
%   would differ only if you select another spectrum based on a different
%   triplet of parameters (ag, Fo, T*) to define the capacity
%
%   Lo spettro deve essere fornito in [metri - g]
%   La massa deve essere definita in maniera coerente con la forza nella
%   pushover. e.g. Pushover in kN, massa in kN/g
%
%
%   INPUT
%
%   ADRSspectrumLS  5% damping, ADRS spectrum for a given intensity (or LS)
%   pushover        pushover curve
%   DISPmechanism   Displacement at the development of the full plastic
%                   mechanism (for bilinearisation) [unit consistent with PO]
%   DISPcapacityLS  displacement capacity for a LS consistent with the
%                   selected ADRS [unit consistent with the pushover]
%   floorMASSES     Mass of each floor [unit consistent with pushover] 
%                   e.g. Pushover in kN, massa in kN/g]
%   PHI1            first modal shape of the system
%   plotter         'plot/noplot'

%% Example

% ADRSspectrumLS= [0,0.356939000000000;0.000188707833496794,0.474636000000000;0.000197305354375608,0.476990000000000;0.000207250202706988,0.479638000000000;0.000217487361514447,0.482287000000000;0.000228019146043329,0.484935000000000;0.000238849186280760,0.487583000000000;0.000252495891044776,0.490820000000000;0.000265297306034601,0.493762000000000;0.000279819960536491,0.496999000000000;0.000296195889566287,0.500530000000000;0.000313135550388310,0.504061000000000;0.000332131634261809,0.507886000000000;0.000351807838722964,0.511711000000000;0.000375369253572282,0.516125000000000;0.000399863067162004,0.520538000000000;0.000425305581112906,0.524952000000000;0.000455304662100784,0.529954000000000;0.000488439968950824,0.535250000000000;0.000524975350849041,0.540841000000000;0.000565195151404757,0.546726000000000;0.000611564580633147,0.553199000000000;0.000635599985137318,0.556436000000000;0.000662475598276850,0.559967000000000;0.000690042553940875,0.563498000000000;0.000718308433028745,0.567029000000000;0.000749726981842801,0.570854000000000;0.000784500343125123,0.574973000000000;0.000820259915965028,0.579093000000000;0.000857014957194100,0.583212000000000;0.000897516412943632,0.587626000000000;0.000942009263002863,0.592333000000000;0.000987856232472527,0.597041000000000;0.00103807206952885,0.602044000000000;0.00109294995354858,0.607340000000000;0.00115281244030097,0.612930000000000;0.00121469499481279,0.618521000000000;0.00128204746626547,0.624406000000000;0.00135524808896455,0.630585000000000;0.00143470344543515,0.637058000000000;0.00152466942151269,0.644120000000000;0.00161812347293649,0.651182000000000;0.00169871010507974,0.657067000000000;0.00178179622463277,0.662952000000000;0.00186741692835183,0.668837000000000;0.00200067562723999,0.677664000000000;0.00209278628501490,0.683549000000000;0.00223594429063482,0.692376000000000;0.00233476425222237,0.698261000000000;0.00248815058299698,0.707088000000000;0.00264783180177572,0.715915000000000;0.00281393023606115,0.724743000000000;0.00298655675398161,0.733570000000000;0.00322709360051803,0.745340000000000;0.00341542169501413,0.754167000000000;0.00367732420679258,0.765937000000000;0.00395181672691765,0.777706000000000;0.00423918990394306,0.789476000000000;0.00461693791863132,0.804188000000000;0.00493430819830476,0.815958000000000;0.00535044062235630,0.830670000000000;0.00587901378799521,0.848325000000000;0.00634447338124399,0.863037000000000;0.00695335985993069,0.873332000000000;0.00742731940969782,0.873332000000000;0.00800002053233310,0.873332000000000;0.00868057783456283,0.873332000000000;0.00903127317907916,0.873332000000000;0.00938891298586315,0.873332000000000;0.00984572839440702,0.873332000000000;0.0102189932412932,0.873332000000000;0.0106953399499649,0.873332000000000;0.0111825373809297,0.873332000000000;0.0117814972515145,0.873332000000000;0.0122925662715244,0.873332000000000;0.0129201720489633,0.873332000000000;0.0135634028665044,0.873332000000000;0.0142222587241477,0.873332000000000;0.0150106722059719,0.873332000000000;0.0158203531034908,0.873332000000000;0.0167717444341588,0.873332000000000;0.0177509136138975,0.873332000000000;0.0187578606427068,0.873332000000000;0.0199238792603345,0.873332000000000;0.0211250542181921,0.873332000000000;0.0226409001225526,0.873332000000000;0.0240645148874209,0.873332000000000;0.0258301444189710,0.873332000000000;0.0276582741109299,0.873332000000000;0.0297092776387913,0.873332000000000;0.0321669662382019,0.873332000000000;0.0347223113382513,0.873332000000000;0.0377364249768574,0.873332000000000;0.0410645585186288,0.873332000000000;0.0449274156550092,0.873332000000000;0.0491702650860977,0.873332000000000;0.0542536114660177,0.873332000000000;0.0600426888238876,0.873332000000000;0.0635270351154456,0.826989000000000;0.0671832848857745,0.781983000000000;0.0714108427183093,0.735690000000000;0.0762095747757297,0.689364000000000;0.0815797030755833,0.643986000000000;0.0878638269291779,0.597927000000000;0.0951763202548472,0.551988000000000;0.103859835458552,0.505837000000000;0.114257286226744,0.459806000000000;0.119970020080725,0.437910000000000;0.126825698041336,0.414240000000000;0.134823562561345,0.389666000000000;0.142821685436619,0.367845000000000;0.151962279914402,0.345719000000000;0.163388093421329,0.321543000000000;0.175956029451730,0.298575000000000;0.190809713646313,0.275333000000000;0.207948541329019,0.252641000000000;0.228514572453488,0.229903000000000;0.253651396082672,0.207120000000000;0.285642594341353,0.183922000000000;0.304411017817532,0.149768000000000;0.304411985612325,0.110475000000000;0.304410438139586,0.0765650000000000;0.304400498531471,0.0122500000000000]; 
% Tc            = 1;
% pushover      = [0, 0; 12.54/1000, 132.5; 90/1000, 132.5]; % capacity curve [m, kN]
% floorMASSES   = 3278.5/8*ones(1,8) / 9.81; % effective mass [Ton]
% PHI1          = (1 : 8) / 8;
% DISPmechanism = pushover(end,1);
% DISPcapacityLS= pushover(end,1);
% DEMretPERIOD  = 1000; % years
% 
% [IS, CAPretPERIOD, dispDEM, accDEM] = N2is( ADRSspectrumLS, pushover, DISPmechanism, DISPcapacityLS, floorMASSES, PHI1, 'plot', DEMretPERIOD, Tc );

%% Optional input

% Maximum number of optional inputs
numvarargs = length(varargin);
if numvarargs > 2
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 2 optional inputs');
end

% set defaults for optional inputs
optargs = {475 1};

% now put these defaults into the valuesToUse cell array,
% and overwrite with the ones specified in varargin.
optargs(1:numvarargs) = varargin;
% or ...
% [optargs{1:numvarargs}] = varargin{:};

% Place optional args in memorable variable names
[ DEMretPERIOD, Tc ] = optargs{:};

%% Calculate main parameters

% Bi-linearise curve
[ PUSHbilin ] = bilinEURO( pushover, DISPmechanism, 'noplot' );

% mass of the SDOF system
MASSsdof = sum(floorMASSES.*PHI1);

% participation factor
gamma = MASSsdof / (sum(floorMASSES.*PHI1.^2));

% period of the SDOF system
Tsdof = 2*pi*( MASSsdof * PUSHbilin(2,1) / PUSHbilin(2,2) )^0.5;

% transform the curve in Capacity spectrum format
CAPspectrum = PUSHbilin ./ gamma;
CAPspectrum(:,2) = CAPspectrum(:,2) / (MASSsdof*9.81);

% yield parameters
ACCyieldSDOF = CAPspectrum(2,2);
DISPyieldSDOF= CAPspectrum(2,1);

%% Perform the N2 method

% intersect demand and capacity
[SdEL, SaEL] = polyxpoly(ADRSspectrumLS(:,1), ADRSspectrumLS(:,2), ...
    CAPspectrum(1:2,1)*1000, CAPspectrum(1:2,2)*1000); % 1000 stands for infinite strength of this system

if ACCyieldSDOF < SaEL % response is nonlinear
    
    q = SaEL / ACCyieldSDOF;
    if Tsdof < Tc % equal energy rule
        % inelastic displacement demand
        PPdisp = SdEL / q * (1 + (q-1)*Tc/Tsdof);
        PPdisp = max(PPdisp, SdEL); % deve essere cmq maggiore dello spost el.
    else % equal displacement rule
        % inelastic displacement demand
        PPdisp = SdEL;
    end
    
    % inelastic acceleration demand
    PPacc = interp1(CAPspectrum(:,1), CAPspectrum(:,2), PPdisp);
else % response is elastic
    PPdisp = SdEL;
    PPacc = SaEL;
end

% for plot
PPdispPLOT = PPdisp;
PPaccPLOT = interp1(CAPspectrum(:,1), CAPspectrum(:,2), PPdispPLOT, 'linear', 'extrap');

% if the structure can't widthstand the demand, PPdisp=NaN
if isnan(PPacc)
    PPdisp = NaN;
end

%% Calculate Safety Index

% convert displacement in sdof displacement
DISPcapLSsdof = DISPcapacityLS / gamma;

% calculate the elastic displacement capacity
if Tsdof < Tc % equal energy rule
    % this derive from inverting Eq. 1.12 in PhD thesis, also using 1.13
    % and the proportion Sa(Y)/Sd(Y) = Sa(EL)/Sd(Y)
    DISPcapLSsdofEL = ( (DISPcapLSsdof/DISPyieldSDOF + Tc/Tsdof -1)*DISPyieldSDOF*Tsdof ) / Tc ;
else % equal displacement rule
    DISPcapLSsdofEL = DISPcapLSsdof;
end

% elastic acceleration capacity (interpolate on the Tsdof line)
ACCcapLSsdofEL = DISPcapLSsdofEL * ACCyieldSDOF / DISPyieldSDOF;

% calculate length of the diagonal on the Tsdof line
diagCAP = ( DISPcapLSsdofEL^2 + ACCcapLSsdofEL^2 )^.5;
diagDEM = ( SdEL^2 + SaEL^2 )^.5;

% calculate Safety Index (conceptually equal to NBS)
IS = diagCAP / diagDEM;

% calculate the capacity return period
etha = 1/0.41;
CAPretPERIOD = DEMretPERIOD * IS^etha; % [years]

%% Plot

if strcmpi(plotter, 'plot')
    
    f = figure;
    hold on
    plot(ADRSspectrumLS(:,1), ADRSspectrumLS(:,2), '--k', 'LineWidth', 2)
    plot(ADRSspectrumLS(:,1)*IS, ADRSspectrumLS(:,2)*IS, 'k', 'LineWidth', 2)
    plot(CAPspectrum(:,1), CAPspectrum(:,2), 'k', 'LineWidth', 2)
    plot(PPdispPLOT, PPaccPLOT, 'o')
    
    xlabel('Displacement [m]')
    ylabel('Acceleration [g]')
    legend('El. demand', 'El. capacity', 'Cap. Spectrum', 'Perf. Point')
    f.CurrentAxes.FontSize = 14;
    
end

end
