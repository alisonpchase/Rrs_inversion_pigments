function pigments_from_Rrs = Rrs_inversion_pigments(Rrs,Rrs_unc,wl,tem,sal)

% Inversion algorithm to estimate phytoplankton pigments from Rrs spectra
%
% Ali Chase, University of Maine, 2017 - 2019
%
% See the following publication for details on the method:
%
%   Chase, A., E. Boss, I. Cetinic, and W. Slade. 2017. "Estimation of Phytoplankton
%   Accessory Pigments from Hyperspectral Reflectance Spectra: Toward a
%   Global Algorithm."
%   Journal of Geophysical Research: Oceans, doi: 10.1002/2017JC012859.
%
% NOTE: This code was developed to estimate phytoplankton pigments from
% hyperspectral remote-sensing reflectance (Rrs) data measured in situ at
% the ocean surface. The Rrs data used in this algorithm development were
% quality controlled, corrected for Raman scattering, normalized to
% eliminate the angular effect of the sun position in the sky relative to
% nadir. Please see above reference for details of these steps.
%
% INPUTS:
%
% Rrs     -  Remote-sensing reflectance measurements, defined as Lw/Ed
% Rrs_unc -  Uncertainy values in Rrs measurements (e.g. could be the standard deviation
%               in Rrs for a given five-minute sample collection), must be on the same 
%               wavlength grid as the Rrs data
% wl      -  wavelengths associated with Rrs (and Rrs_unc) data
% tem     -  water temperature at the time of radiometery data collection
% sal     -  water salinity at the time of radiometery data collection
%
% OUTPUTS:
%
% pigments_from_Rrs - a structure with the following fields:
%               
%   pigments_from_Rrs.est_pigm   - Estimated pigment concentrations 
%   pigments_from_Rrs.pigm_unc   - Uncertainties in estimated pigments,
%               calculated using a Monte Carlo method that in turn uses the 
%               reported uncertainties in the A and B coefficients reported
%               in Chase et al. (2017). 
%   pigments_from_Rrs.vars_units - The names and units of the estimated pigments:
%               chlorophyll a (Chla), chlorophyll b (Chlb), chlorophyll c1
%               +c2 (Chlc12), and photoprotective carotenoids (PPC) defined
%               as abcarotene+zeaxanthin+alloxanthin+diadinoxanthin. All
%               pigments and uncertainties are in mg m^-3.
%
% -------------------------------------------------------------------------
%
% Constants:
g1 = 0.0949; % g1 and g2 are values from Gordon et al., 1988
g2 = 0.0794;
lnot = 400; % reference lambda wavelength (nm)

% The data are cut off at 600 nm to reduce the influence of red wavelengths,
% which have low signal to noise, on the inversion
wlfull = wl;
Iwl = find(wlfull > 400 & wlfull < 600);  
wl = wlfull(Iwl);

Rrs = Rrs(Iwl);
Rrs_unc = Rrs_unc(Iwl);

% Get the absorption and backscattering by water for the temperature and
% salinity measured coincidently with Rrs
[a_sw,bb_sw] = get_water_iops(wl,tem,sal);

% Calculate rrs from Rrs (method from Lee et al., 2002)
rrs = Rrs./(0.52 + 1.7.*Rrs);
rrs_unc = Rrs_unc./(0.52 + 1.7.*Rrs_unc);

% Calculate the positive solution for U using rrs, where U = bb/(bb+a).
% This model and g coefficients are from Gordon et al., 1988
Upos = zeros(length(wl),1);
Uunc = zeros(length(wl),1);
for ii = 1:length(wl)
    Upos(ii) = (-g1 + sqrt(g1.^2 + 4*g2*rrs(ii)))./(2*g2);
    Uunc(ii) = (-g1 + sqrt(g1.^2 + 4*g2*rrs_unc(ii)))./(2*g2);
end

% the U spectra must rows, not columns to pass to the inversion; check for
% this
sz = size(Upos);
if sz(1) > sz(2)
    Upos = Upos';   
end

sz2 = size(Uunc);
if sz2(1) > sz(2)
    Uunc = Uunc';
end

% Define the center peak locations (nm) and widths (nm) of the Gaussian functions
% sig = sigma, where FWHM = sigma*2.355 and FWHM = full width at half max
peaks = [384,413,435,461,464,490,532,583,623,644,655,676];
sig = [23,9,14,11,19,19,20,20,15,12,12,9];

% Define the [lower bound, first guess, upper bound] for each parameter. These will be allowed to vary.
s_nap = [.005 .011 .016]; % slope of non-algal particle absorption (nm^-1)
m_nap = [.0 .005 .05]; % magnitude of non-algal particle absorption (m^-1)
s_cdom = [.005 .0185 .02]; % slope of CDOM absorption (0.0185 first guess from Matsuoka et al., 2013) (nm^-1)
m_cdom = [.01 .1 .8]; % magnitude of CDOM absorption (m^-1)
bbpbp_ratio = [.005 .01 .015];  % ratio of particulate backscattering to total particulate scattering
m_gaus = [.0 .01 0.5]; % magnitude of Gaussian functions representing pigment absorption (m^-1)
cpgam = [.0 1 1.3]; % slope of gamma (cp)
m_cp = [.01 .1 1]; % magnitude of cp spectrum

C_fluor = [.0 .001 .01]; % magnitude of Gaussian representing fluorescence

% First guess array
Amp0 = [s_nap(2), m_nap(2), s_cdom(2), m_cdom(2), bbpbp_ratio(2), cpgam(2), m_cp(2),C_fluor(2),...
    repmat(m_gaus(2),1,length(peaks)),peaks,sig];

% Lower bound array
LB = [s_nap(1), m_nap(1), s_cdom(1), m_cdom(1), bbpbp_ratio(1), cpgam(1), m_cp(1),C_fluor(1),...
    repmat(m_gaus(1),1,length(peaks)),peaks-1,sig-1];

% Upper bound array
UB = [s_nap(3), m_nap(3), s_cdom(3), m_cdom(3), bbpbp_ratio(3), cpgam(3), m_cp(3),C_fluor(3),...
    repmat(m_gaus(3),1,length(peaks)),peaks+1,sig+1];

% Run the inversion using a non-linear least squares inversion function
[Amp2,~,~,~]=...
    lsqnonlin(@lsqnonlin_Amp_gen,Amp0,LB,UB,[],Upos,Uunc,wl,bb_sw,a_sw,lnot);

 %rebuild spectral components with the results of the inversion 
 cdom = exp(-Amp2(3)*(wl-lnot));
 nap = exp(-Amp2(1)*(wl-lnot));
 cp = (wl./lnot).^-Amp2(6);
 
 CDOM = (Amp2(4)*cdom);
 NAP = (Amp2(2)*nap);
 CP = (Amp2(7)*cp);

 peak_locs=Amp2(21:32);
 sigs=Amp2(33:44);
 clear gaus
 % define Gaussian shapes
 for ii=1:max(size(peak_locs))
     for jj=1:max(size(wl))
         gaus(jj,ii) = exp(-0.5.*((wl(jj)-peak_locs(ii))/sigs(ii))^2);
     end
 end
 
 % multiply each Gaussian by the initial guess amplitude
 clear GAUS
 for i=1:length(peak_locs)
     GAUS(:,i) = Amp2(i+8)*gaus(:,i);
 end
 
 % sum all of the Gaussians to get a_phi
 clear APHI
 for j=1:length(wl)
     APHI(j)=0;
     for i=1:length(peak_locs)
         APHI(j) = APHI(j) + sum(GAUS(j,i));
     end
 end
 
 % particulate absorption 
 AP = NAP + APHI';
 
 % particulate backscattering
 BBP = Amp2(5).*(CP-AP);
 
 % Define the fluorescence Gaussian; 685 peak for Chl fluor; 10.6 sigma
 fluor = [];
 for i=1:max(size(wl))
     fluor(i) = exp(-0.5.*((wl(i)-685)/10.6)^2);
 end
 F = Amp2(8)*fluor;
 
 denom = AP' + CDOM' + a_sw' + BBP' + bb_sw;
 numer = BBP' + bb_sw + F;
 Unew = numer./denom;
 rrs_recon = g1.*Unew + g2.*(Unew.^2);      
 
% Uncomment below to plot the measured and rebuilt U spectra to see visually the
% fit of the inversion; can also plot the original Rrs and the
% reconstructed Rrs ("Rrs_recon", just above)

%  figure
%  plot(wl,Upos);hold on
%  plot(wl,Unew)
 
 % Estimate pigment concentrations and their uncertainties with coefficients reported in
 % Chase et al., 2017 (JGR-Oceans) and a Monte Carlo method
 % matrix:   A   A_unc   B   B_unc
 coeffs = [0.048 0.008 0.643 0.068;...
     0.033 0.013 0.327 0.074;...
     0.043 0.009 0.561 0.059;...
     0.079 0.024 0.823 0.105];
 
 % loop through the four pigments that are estimated
 for ii = 1:4
     mc = randn(10000,1) * coeffs(ii,[2 4]);
     As = coeffs(ii,1) + mc(:,1);
     As(As<0)=0; % prevent imaginary pigment concentration values
     Bs = coeffs(ii,3) + mc(:,2);
     pigest = (Amp2(10+ii)./As').^(1./Bs');
     pigmedian(ii) = median(pigest);
     prc = prctile(pigest,[16 84]);
     pigunc(ii) = (prc(2) - prc(1))/2;
 end
 
 % populate structure
 pigments_from_Rrs.est_pigm = pigmedian;
 pigments_from_Rrs.pigm_unc = pigunc;
 pigments_from_Rrs.vars_units = 'Chla, Chlb, Chlc1+c2, PPC; all in mg m^-3';
 
end

function [a_sw,bb_sw] = get_water_iops(wave,T,S)

%function to obtain seawater absorption and backscattering spectra

%pure water absorption from Mason et al 2016 for 250-550, Pope and Frye for
%550-730 nm, and Smith and Baker for 730-800 nm
%salt water backscattering from Zhang et al 2009
%corrected for in situ temperature and salinity conditions Sullivan et al. 2006

wl_water1 = [250   260   270   280   290   300   302   304   306   308   310   312   314   316   318   320   322   324   326   328   330 ...
    332   334   336   338   340   342   344   346   348   350   352   354   356   358   360   362   364   366   368   370   372 ...
    374   376   378   380   382   384   386   388   390   392   394   396   398   400   402   404   406   408   410   412   414 ...
    416   418   420   422   424   426   428   430   432   434   436   438   440   442   444   446   448   450   452   454   456 ...
    458   460   462   464   466   468   470   472   474   476   478   480   482   484   486   488   490   492   494   496   498 ...
    500   502   504   506   508   510   512   514   516   518   520   522   524   526   528   530   532   534   536   538   540 ...
    542   544   546   548   550 ];
wl_water2 = [552.5   555 557.5	560	562.5	565	567.5	570	572.5	575	577.5	580	582.5	585	587.5	590	592.5	595	597.5 ...
    600	602.5	605	607.5	610	612.5	615	617.5	620	622.5	625	627.5	630	632.5	635	637.5	640	642.5 ...
    645	647.5	650	652.5	655	657.5	660	662.5	665	667.5	670	672.5	675	677.5	680	682.5	685	687.5 ...
    690	692.5	695	697.5	700	702.5	705	707.5	710	712.5	715	717.5	720	722.5	725	727.5	730  732.5 ...
    735.0  737.5  740.0  742.5  745.0  747.5  750.0  752.5  755.0  757.5  760.0  762.5  765.0  767.5  770.0  772.5 ...
    775.0  777.5 780.0  782.5  785.0  787.5  790.0  792.5  795.0  797.5  800.0];

wl_water = [wl_water1 wl_water2];

aw1 = [58.7100   51.5000   43.5700   22.3000    9.3900    4.6700    4.3300    3.5600    3.1300    2.7500    2.3600    2.0500 ...
    1.8500    1.7600    1.6300    1.4700    1.3300    1.2400    1.1800    1.1200    1.0700    1.0100    0.9900    0.9500 ...
    0.9100    0.8500    0.8200    0.8100    0.8200    0.8400    0.8900    0.9400    0.9700    0.9800    0.9900    1.0600 ...
    1.1500    1.2000    1.2100    1.2200    1.2400    1.2700    1.2900    1.3300    1.3700    1.4300    1.4700    1.5100 ...
    1.5500    1.6200    1.7000    1.7500    1.8500    1.9600    2.0800    2.2200    2.3700    2.4800    2.5700    2.5900 ...
    2.6600    2.7100    2.8000    2.8800    3.0000    3.1200    3.2200    3.3100    3.4400    3.5800    3.7600    3.9500 ...
    4.1700    4.4200    4.8000    5.2200    5.7400    6.2600    6.9100    7.5100    8.0800    8.4200    8.6300    8.7700 ...
    8.9300    9.0900    9.3300    9.5500    9.7900    9.9900   10.3000   10.6500   11.0000   11.3800   11.7700   12.1400 ...
    12.5400   12.9400   13.3600   13.9100   14.6000   15.4500   16.4800   17.7400   19.2600   20.7300   22.4200   24.2400 ...
    26.6800   29.7100   33.0000   35.6900   37.3800   38.2100   38.7800   39.1700   39.6200   40.1700   40.8800   41.6200 ...
    42.4200   43.3000   44.3600   45.4100   46.4500   47.5400   48.8200   50.4000   52.2400   54.2500   56.2900];

aw2 = [0.0593 0.0596 0.0606 0.0619 0.064...
    0.0642 0.0672 0.0695 0.0733 0.0772 0.0836 0.0896 0.0989 0.11 0.122 0.1351 0.1516 0.1672 0.1925 0.2224 0.247 0.2577...
    0.2629 0.2644 0.2665 0.2678 0.2707 0.2755 0.281 0.2834 0.2904 0.2916 0.2995 0.3012 0.3077 0.3108 0.322 0.325 0.335...
    0.34 0.358 0.371 0.393 0.41 0.424 0.429 0.436 0.439 0.448 0.448 0.461 0.465 0.478 0.486 0.502 0.516 0.538 0.559...
    0.592 0.624	0.663 0.704 0.756 0.827 0.914 1.007 1.119 1.231 1.356 1.489 1.678 1.7845 1.9333 2.0822 2.2311 2.3800...
    2.4025 2.4250 2.4475 2.4700 2.4900 2.5100 2.5300 2.5500 2.5400 2.5300 2.5200 2.5100 2.4725 2.4350 2.3975 2.3600...
    2.3100 2.2600 2.2100 2.1600 2.1375 2.1150 2.0925 2.0700];

a_water = [aw1*1e-3 aw2]; % 1e-3 comes from Mason et al. 2016, Table 2

a_pw = interp1(wl_water,a_water,wave,'linear','extrap');

%temp and salinity correction for water absorption (need to know at what T it was measured):
if S==0 || isnan(S)==1,S=35;end
if T==0 || isnan(T)==1,T=22;end
T_pope=22.0;

% use salt water scattering from Zhang et al 2009
[betasw124, bb_sw, beta90sw, theta] = betasw124_ZHH2009(wave, S, T);

% use the temp and salinity corrections from Sullivan et al. 2006
[psiT,psiS] = tempsal_corr(wave);
%temperature and salinity corrections:
a_sw = (a_pw+psiT*(T-T_pope)+psiS*S);

end

function [psiT,psiS] = tempsal_corr(wavel)
tmp = [400.0000    0.0001         0
    402.0000    0.0001         0
    404.0000    0.0001         0
    406.0000    0.0001         0
    408.0000         0    0.0000
    410.0000         0    0.0000
    412.0000         0    0.0000
    414.0000    0.0001    0.0000
    416.0000    0.0001    0.0000
    418.0000         0    0.0000
    420.0000         0    0.0000
    422.0000         0         0
    424.0000         0         0
    426.0000         0         0
    428.0000         0         0
    430.0000         0         0
    432.0000         0         0
    434.0000         0         0
    436.0000         0         0
    438.0000         0         0
    440.0000         0         0
    442.0000         0         0
    444.0000         0         0
    446.0000         0         0
    448.0000         0         0
    450.0000         0         0
    452.0000         0         0
    454.0000         0         0
    456.0000         0         0
    458.0000         0         0
    460.0000         0         0
    462.0000         0         0
    464.0000         0         0
    466.0000         0         0
    468.0000         0         0
    470.0000         0         0
    472.0000         0         0
    474.0000         0         0
    476.0000         0         0
    478.0000         0         0
    480.0000         0   -0.0000
    482.0000         0   -0.0000
    484.0000         0   -0.0000
    486.0000         0   -0.0000
    488.0000         0   -0.0000
    490.0000         0   -0.0000
    492.0000         0   -0.0000
    494.0000         0   -0.0000
    496.0000         0   -0.0000
    498.0000         0   -0.0000
    500.0000         0   -0.0000
    502.0000         0   -0.0000
    504.0000         0   -0.0000
    506.0000         0   -0.0000
    508.0000    0.0001   -0.0000
    510.0000    0.0001   -0.0000
    512.0000    0.0001   -0.0000
    514.0000    0.0001   -0.0000
    516.0000    0.0001         0
    518.0000    0.0001         0
    520.0000    0.0001         0
    522.0000    0.0001         0
    524.0000    0.0001         0
    526.0000    0.0001         0
    528.0000         0         0
    530.0000         0         0
    532.0000         0         0
    534.0000         0         0
    536.0000         0         0
    538.0000         0         0
    540.0000         0         0
    542.0000         0   -0.0000
    544.0000         0   -0.0000
    546.0000         0         0
    548.0000         0         0
    550.0000         0         0
    552.0000         0         0
    554.0000         0         0
    556.0000         0         0
    558.0000         0         0
    560.0000         0   -0.0000
    562.0000         0   -0.0000
    564.0000         0   -0.0000
    566.0000         0   -0.0000
    568.0000         0   -0.0000
    570.0000         0   -0.0000
    572.0000         0   -0.0000
    574.0000         0   -0.0000
    576.0000    0.0001   -0.0000
    578.0000    0.0001   -0.0000
    580.0000    0.0002   -0.0000
    582.0000    0.0002   -0.0000
    584.0000    0.0003   -0.0000
    586.0000    0.0004   -0.0000
    588.0000    0.0004   -0.0000
    590.0000    0.0005   -0.0000
    592.0000    0.0006   -0.0000
    594.0000    0.0008   -0.0000
    596.0000    0.0009   -0.0000
    598.0000    0.0010   -0.0000
    600.0000    0.0011   -0.0000
    602.0000    0.0011         0
    604.0000    0.0011    0.0000
    606.0000    0.0011    0.0000
    608.0000    0.0011    0.0000
    610.0000    0.0010    0.0000
    612.0000    0.0009    0.0001
    614.0000    0.0008    0.0001
    616.0000    0.0007    0.0001
    618.0000    0.0007    0.0001
    620.0000    0.0006    0.0001
    622.0000    0.0005    0.0001
    624.0000    0.0004    0.0001
    626.0000    0.0003    0.0001
    628.0000    0.0002    0.0001
    630.0000    0.0002    0.0001
    632.0000    0.0001    0.0001
    634.0000         0    0.0001
    636.0000         0    0.0000
    638.0000         0    0.0000
    640.0000   -0.0001    0.0000
    642.0000   -0.0001    0.0000
    644.0000   -0.0001    0.0000
    646.0000   -0.0001    0.0000
    648.0000   -0.0001    0.0000
    650.0000         0    0.0000
    652.0000         0    0.0000
    654.0000    0.0001    0.0000
    656.0000    0.0001    0.0000
    658.0000    0.0002    0.0000
    660.0000    0.0002    0.0000
    662.0000    0.0002    0.0000
    664.0000    0.0002    0.0000
    666.0000    0.0002    0.0000
    668.0000    0.0001    0.0000
    670.0000    0.0001    0.0000
    672.0000         0    0.0000
    674.0000         0    0.0000
    676.0000   -0.0001         0
    678.0000   -0.0001   -0.0000
    680.0000   -0.0001   -0.0000
    682.0000   -0.0002   -0.0000
    684.0000   -0.0002   -0.0000
    686.0000   -0.0002   -0.0001
    688.0000   -0.0001   -0.0001
    690.0000   -0.0001   -0.0001
    692.0000   -0.0001   -0.0001
    694.0000         0   -0.0001
    696.0000    0.0001   -0.0001
    698.0000    0.0002   -0.0001
    700.0000    0.0004   -0.0002
    702.0000    0.0006   -0.0002
    704.0000    0.0009   -0.0002
    706.0000    0.0012   -0.0002
    708.0000    0.0017   -0.0002
    710.0000    0.0022   -0.0002
    712.0000    0.0027   -0.0002
    714.0000    0.0033   -0.0002
    716.0000    0.0040   -0.0002
    718.0000    0.0049   -0.0002
    720.0000    0.0060   -0.0002
    722.0000    0.0071   -0.0003
    724.0000    0.0084   -0.0003
    726.0000    0.0096   -0.0003
    728.0000    0.0108   -0.0002
    730.0000    0.0120   -0.0002
    732.0000    0.0130   -0.0002
    734.0000    0.0139   -0.0001
    736.0000    0.0146    0.0000
    738.0000    0.0150    0.0001
    740.0000    0.0150    0.0003
    742.0000    0.0147    0.0004
    744.0000    0.0142    0.0005
    746.0000    0.0138    0.0007
    748.0000    0.0131    0.0008
    750.0000    0.0124    0.0009];

psiT = interp1(tmp(:,1),tmp(:,2), wavel,'linear',0); % both this line and the next changed from 'spline' on 14 Dec 2015
psiS = interp1(tmp(:,1),tmp(:,3), wavel,'linear',0);

end

function [betasw124, bsw, beta90sw, theta] = betasw124_ZHH2009(lambda, S, Tc, delta)
%
% function [betasw124, bsw, beta90sw, theta]= betasw124_ZHH2009(lambda,S,Tc,delta)
%
%
% Scattering by pure seawater: Effect of salinity
% Xiaodong Zhang, Lianbo Hu, and Ming-Xia He, Optics Express, 2009, accepted
% lambda (nm): wavelength
% Tc: temperauter in degree Celsius, must be a scalar
% S: salinity, must be scalar
% delta: depolarization ratio, if not provided, default = 0.039 will be
% used.
% betasw: volume scattering at angles defined by theta. Its size is [x y],
% where x is the number of angles (x = length(theta)) and y is the number
% of wavelengths in lambda (y = length(lambda))
% beta90sw: volume scattering at 90 degree. Its size is [1 y]
% bw: total scattering coefficient. Its size is [1 y]
% for backscattering coefficients, divide total scattering by 2
%
% Xiaodong Zhang, March 10, 2009
%
%
% MODIFIED on 17/05/2011 to be able to process bbp profiles with coincident T and sal profiles
% MODIFIED on 05 Apr 2013 to use 124 degs instead of 117 degs
%

% make sure that S and Tc are column vectors
S = S(:);
Tc = Tc(:);

% values of the constants
Na = 6.0221417930e23 ;   %  Avogadro's constant
Kbz = 1.3806503e-23 ;    %  Boltzmann constant
Tk = Tc+273.15 ;         %  Absolute tempearture
M0 = 18e-3;              %  Molecular weigth of water in kg/mol

%error(nargchk(3, 4, nargin));
if nargin == 3
    delta = 0.039; % Farinato and Roswell (1976)
end

theta=[0:.01:180];

%if ~isscalar(Tc) || ~isscalar (S)
%    error('Both Tc and S need to be scalar variable');
%end

lambda = lambda(:)'; % a row variable
rad = theta(:)*pi/180; % angle in radian as a colum variable

% nsw: absolute refractive index of seawater
% dnds: partial derivative of seawater refractive index w.r.t. salinity
[nsw dnds] = RInw(lambda,Tc,S);

% isothermal compressibility is from Lepple & Millero (1971,Deep
% Sea-Research), pages 10-11
% The error ~ +/-0.004e-6 bar^-1
IsoComp = BetaT(Tc,S);

% density of water and seawater,unit is Kg/m^3, from UNESCO,38,1981
density_sw = rhou_sw(Tc, S);

% water activity data of seawater is from Millero and Leung (1976,American
% Journal of Science,276,1035-1077). Table 19 was reproduced using
% Eq.(14,22,23,88,107) then were fitted to polynominal equation.
% dlnawds is partial derivative of natural logarithm of water activity
% w.r.t.salinity
dlnawds = dlnasw_ds(Tc, S);

% density derivative of refractive index from PMH model
DFRI = PMH(nsw);  %% PMH model

% volume scattering at 90 degree due to the density fluctuation
beta_df = pi*pi./2.*((lambda*1e-9).^(-4)).*Kbz.*Tk.*IsoComp.*DFRI.^2.*(6+6*delta)./(6-7*delta);

% volume scattering at 90 degree due to the concentration fluctuation
flu_con = S.*M0.*dnds.^2./density_sw./(-dlnawds)./Na;
beta_cf = 2*pi*pi*((lambda*1e-9).^(-4)).*nsw.^2.*(flu_con)*(6+6*delta)/(6-7*delta);

% total volume scattering at 90 degree
beta90sw = beta_df+beta_cf;
bsw = 8*pi/3.*beta90sw.*(2+delta)./(1+delta);

rad124 = find (rad2deg(rad)>=124, 1);

betasw124 = NaN(size(beta90sw));
%for i=1:length(beta90sw)
%    betasw124(i,:) = beta90sw(i) * (  1+((cos(rad(rad124))).^2).*(1-delta)./(1+delta)  );
%end
betasw124 = beta90sw * (  1+((cos(rad(rad124))).^2).*(1-delta)./(1+delta)  );

end

function [nsw dnswds]= RInw(lambda,Tc,S)
% refractive index of air is from Ciddor (1996,Applied Optics)
n_air = 1.0+(5792105.0./(238.0185-1./(lambda/1e3).^2)+167917.0./(57.362-1./(lambda/1e3).^2))/1e8;

% refractive index of seawater is from Quan and Fry (1994, Applied Optics)
n0 = 1.31405; n1 = 1.779e-4 ; n2 = -1.05e-6 ; n3 = 1.6e-8 ; n4 = -2.02e-6 ;
n5 = 15.868; n6 = 0.01155;  n7 = -0.00423;  n8 = -4382 ; n9 = 1.1455e6;

nsw = n0+(n1+n2*Tc+n3*Tc.^2).*S+n4*Tc.^2+(n5+n6*S+n7*Tc)./lambda+n8./lambda.^2+n9./lambda.^3; % pure seawater
nsw = nsw.*n_air;
dnswds = (n1+n2*Tc+n3*Tc.^2+n6./lambda).*n_air;

end

function IsoComp = BetaT(Tc, S)
% pure water secant bulk Millero (1980, Deep-sea Research)
kw = 19652.21+148.4206*Tc-2.327105*Tc.^2+1.360477e-2*Tc.^3-5.155288e-5*Tc.^4;
Btw_cal = 1./kw;

% isothermal compressibility from Kell sound measurement in pure water
% Btw = (50.88630+0.717582*Tc+0.7819867e-3*Tc.^2+31.62214e-6*Tc.^3-0.1323594e-6*Tc.^4+0.634575e-9*Tc.^5)./(1+21.65928e-3*Tc)*1e-6;

% seawater secant bulk
a0 = 54.6746-0.603459*Tc+1.09987e-2*Tc.^2-6.167e-5*Tc.^3;
b0 = 7.944e-2+1.6483e-2*Tc-5.3009e-4*Tc.^2;

Ks =kw + a0.*S + b0.*S.^1.5;

% calculate seawater isothermal compressibility from the secant bulk
IsoComp = 1./Ks*1e-5; % unit is pa
end

function density_sw = rhou_sw(Tc, S)

% density of water and seawater,unit is Kg/m^3, from UNESCO,38,1981
a0 = 8.24493e-1;  a1 = -4.0899e-3; a2 = 7.6438e-5; a3 = -8.2467e-7; a4 = 5.3875e-9;
a5 = -5.72466e-3; a6 = 1.0227e-4;  a7 = -1.6546e-6; a8 = 4.8314e-4;
b0 = 999.842594; b1 = 6.793952e-2; b2 = -9.09529e-3; b3 = 1.001685e-4;
b4 = -1.120083e-6; b5 = 6.536332e-9;

% density for pure water
density_w = b0 + b1.*Tc + b2.*Tc.^2 + b3.*Tc.^3 + b4.*Tc.^4 + b5.*Tc.^5;
% density for pure seawater
density_sw = density_w +((a0 + a1.*Tc + a2.*Tc.^2 + a3.*Tc.^3 + a4.*Tc.^4).*S + (a5 + a6.*Tc + a7.*Tc.^2).*S.^1.5 + a8.*S.^2);

end

function dlnawds = dlnasw_ds(Tc, S)
% water activity data of seawater is from Millero and Leung (1976,American
% Journal of Science,276,1035-1077). Table 19 was reproduced using
% Eqs.(14,22,23,88,107) then were fitted to polynominal equation.
% dlnawds is partial derivative of natural logarithm of water activity
% w.r.t.salinity
% lnaw = (-1.64555e-6-1.34779e-7*Tc+1.85392e-9*Tc.^2-1.40702e-11*Tc.^3)+......
%            (-5.58651e-4+2.40452e-7*Tc-3.12165e-9*Tc.^2+2.40808e-11*Tc.^3).*S+......
%            (1.79613e-5-9.9422e-8*Tc+2.08919e-9*Tc.^2-1.39872e-11*Tc.^3).*S.^1.5+......
%            (-2.31065e-6-1.37674e-9*Tc-1.93316e-11*Tc.^2).*S.^2;

dlnawds = (-5.58651e-4+2.40452e-7*Tc-3.12165e-9*Tc.^2+2.40808e-11*Tc.^3)+...
    1.5*(1.79613e-5-9.9422e-8*Tc+2.08919e-9*Tc.^2-1.39872e-11*Tc.^3).*S.^0.5+...
    2*(-2.31065e-6-1.37674e-9*Tc-1.93316e-11*Tc.^2).*S;

% density derivative of refractive index from PMH model

end

function n_density_derivative=PMH(n_wat)
n_wat2 = n_wat.^2;
n_density_derivative=(n_wat2-1).*(1+2/3*(n_wat2+2).*(n_wat/3-1/3./n_wat).^2);
end


% The following function uses a non-linear least squares solver to minimize
% the difference between the measured Rrs and the modeled Rrs
function spec_min = lsqnonlin_Amp_gen(Amp0,Upos,Uunc,wvns,bb_sw_r,a_sw_r,lnot)

peak_locs = Amp0(21:32); 
sig = Amp0(33:44);

% define cdom and nap functions; both slope and magnitude are allowed to
% vary
cdom = exp(-Amp0(3)*(wvns-lnot));
CDOM = Amp0(4)*cdom;

nap = exp(-Amp0(1)*(wvns-lnot));
NAP = Amp0(2)*nap;

% define Gaussian shapes
for ii=1:max(size(peak_locs))
    for jj=1:max(size(wvns))
        gaus(jj,ii) = exp(-0.5.*((wvns(jj)-peak_locs(ii))/sig(ii))^2);
    end
end

% multiply each Gaussian by the initial guess amplitude
for i=1:length(peak_locs)
    GAUS(:,i) = Amp0(i+8)*gaus(:,i);
end

% Sum all of the Gaussians to get a_phi
for j=1:length(wvns)
    APHI(j)=0;
    for i=1:length(peak_locs)
        APHI(j) = APHI(j) + sum(GAUS(j,i));
    end
end

% define cp (slope and magnitude can vary)
cp = ((wvns./lnot).^-Amp0(6));
CP = (Amp0(7)*cp);

% total particulate absorption
AP = NAP + APHI';

% define bbp in terms of the bbp:bp ratio, ap, and cp following Roesler &
% Boss, 2003
BBP = Amp0(5).*(CP-AP);

% Define the fluorescence Gaussian; 685 peak for Chl fluor; 10.6 sigma
for i=1:max(size(wvns))
    fluor(i) = exp(-0.5.*((wvns(i)-685)/10.6)^2);
end
F = Amp0(8)*fluor;

% modeled U spectrum (U = bb/{a+bb})
denom = APHI + NAP' + CDOM' + a_sw_r' + BBP' + bb_sw_r;
numer = BBP' + bb_sw_r + F;
Unew = numer./denom;

% normalize by uncertainties during the minimization
spec_min = (Upos - Unew)./Uunc;


end
