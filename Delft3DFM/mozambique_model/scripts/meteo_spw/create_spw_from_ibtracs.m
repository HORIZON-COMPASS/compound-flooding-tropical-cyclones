%% Create spiderweb
% v1.0  Nederhoff   Dec-17
% v1.1  Leijnse     Oct-18 Use new wes-script
% v1.2  Leijnse     Jan-2020 - base on Ibtracks data using real RMW (in nautical miles in data)
% v1.3  Leijnse     Feb-2020 - specify parameters using tc instead of spw!
% - made after code change in wes4.m (17.02.2020 10am)
% v1.4  Leijnse     Jan-2021 - select track from ibtracs, and make
% different version using different rainfall relations, and prepare csv of
% tracks to put into CFRSS
%

% steps in script:
% 1) find wanted track
% 2) load data 
% 3) set general TC settings3
%                                                                           
% 4) put data in structure

% 5) set rainfall options per rainfall variation
% 6) make csv for CFRSS 
% 7) make all spiderwebs

clear all
close all
clc

%% 0. Arish changeable user settings
% addpath(genpath('p:\11205281-tcwise\03_currentprojects\MSc_Arish\03_modelling\')); %
    addpath(genpath('p:\11210471-001-compass\02_Models\Delft3DFM\mozambique_model\scripts\meteo_spw\')); %
    addpath(genpath('n:\Deltabox\Bulletin\leijnse\MyBulletin_Github\Matlab\')); %
    % wanted TC
iname = 'IDAI'; % can have different cyclone names in the database
iyear = '2019';  % different years of tropical cyclone (### not working for Harvey 2017 ###)

clear variations_wanted
variations_wanted{1} =  'option_norain';
%variations_wanted{2} =  'option_a';

%% 0. general user settings
foverall                            = 'p:\archivedprojects\11205281-tcwise\'; % to be changed by user  folder of checkout TCWiSE
fmain                               = [foverall,'\00_data\'];
fname                               = [fmain, 'IBTrACS.ALL.v04r00_download261121.nc']; % folder of ibtracks dataset
destout                             = 'p:\11210471-001-compass\02_Models\Delft3DFM\mozambique_model\boundary_conditions\meteo\TC\';

mkdir(destout);
cd(destout);

% wanted source from Ibtracs data:
isource = 'usa';  % #### can it be changed to other countries as well? ####

%% 1. first step
[tc,spw,trackname] = prepare_all(destout, fname, isource, iname, iyear, variations_wanted); 

spw.r35estimate = 1; % new addition from Kees in wes4.m code
% ### what does spw.r35estimate do?
%% 2. parts where you mightwant to change things:

for jj = 1:length(variations_wanted)  % ### what is length varation wanted (is it number of times cyclone captured time interval)
    
    optionname = variations_wanted{jj};
    spw2 = spw; %reload standard spw structure without rainfal selected variables
    %spw2.rain_info.perc = 50; % default  #### Ask for explaination of these parameters####
    spw2.rain_info.asymmetrical = 0;
    spw2.rain_info.random = 0; %random = 1 is mode, random = 1 is stochastic, random = 2 is user specified percentile for BaCla
    
        
    % @Arish, here you could add more options  #### Do you mean different
    % model
%     for ipercentiles=85 %runs different conditions for different percentile senarios
        if strcmp(variations_wanted{jj}, 'option_a')
            optionname = char('ipet_symmetric_mode'); % does option name change makes any change 
            spw2.rain_info.asymmetrical = 0; %added to try if it works        
            spw2.rain_relation = 'ipet';
            spw2.rain_info.rain_relation = spw2.rain_relation;
        elseif strcmp(variations_wanted{jj}, 'option_norain')
            spw2.rainfall = -1;
        end
        
        %% Added check on pressure
        spw2.pn = 1020; %Maarten: better background value
    
        max_pc = nanmax([tc.track.pc]);
    
        if max_pc > spw2.pn
            spw2.pn = max_pc;
          warning('max_pc > spw2.pn')                                                                                                                                                                                                                                                                                                                                                     
        end
    
        %% Actually create spiderweb
    
        trackname_variation = [trackname, '_', optionname];
        mkdir([destout, filesep, trackname])
        
        spwname = [destout, filesep, trackname,'\cyclone_',optionname,'.spw'];
    
%         [tc2,spw3] = wes4_arish_new(tc,'tcstructure',spw2,spwname); 
        %[tc2,spw3] = wes4_new(tc,'tcstructure',spw2,spwname); 
        tc2 = wes4(tc,'tcstructure',spw2,spwname); 

        
        folder = [destout, filesep, trackname,filesep];
        save([folder, 'input_tc_',optionname,'.mat'],'tc2')
        %save([folder,  'input_spw_',optionname,'.mat'],'spw3')
            
       plot_figures(tc2,spwname,folder, trackname, optionname) 
%     end
end
%% Function you don't have to worry about
function [tc,spw, trackname] = prepare_all(destout, fname, isource, iname, iyear, variations_wanted) 
%% 1. find wanted track
% Get TC data from ibtracs:

time     = nc_varget(fname,'time');
t0      = datenum(1858,11,17);
time     = t0+time;

name = nc_varget(fname,'name');

id = 0;
for ii = 1:size(name,1)
    
    nametmp = name(ii,1:length(iname));
    
    timetmp = datestr(time(ii,1),'yyyy');
    
    if strcmp(nametmp,iname) == 1 && strcmp(timetmp,iyear) == 1
        id = ii;
        
        break
    end
end

% check selected track
trackname =   [nametmp, '_',iyear];

disp(['Selected track is: ',trackname,' with genesis at: ',datestr(time(id,1))])

%% 2. Load the data
% data.time = time(id,:);
data.time=datenum(1858,11,17) + nc_varget(fname,'time',[id-1 0],[1 Inf])'; %TL: -1 is needed! otherwise take the wrong timestep

% Read lat/lon
data.x=nc_varget(fname,[isource,'_lon'],[id-1 0],[1 Inf])'; %TL: -1 is needed! otherwise take the wrong timestep
data.y=nc_varget(fname,[isource,'_lat'],[id-1 0],[1 Inf])';
    
% Read Vmax
data.vmax     = nc_varget(fname,[isource,'_wind'],[id-1 0],[1 Inf])'; % in knts

% Read minimum pressure center
data.pc     = nc_varget(fname,[isource,'_pres'],[id-1 0],[1 Inf])';

% Read RMW
data.rmax     = nc_varget(fname,[isource,'_rmw'],[id-1 0],[1 Inf])'; % in nautical mile (not best tracked volgens netcdf-file)

% additional Radial information; availability may depend on Ibtracs data source, see netcdf file
% Read R34 (use as R35)
r35     = nc_varget(fname,[isource,'_r34'],[id-1 0 0],[1 Inf Inf]); % in nautical mile, radius of 34 knot winds, per quadrant
data.r35ne   = (r35(:,1))'; %Watch out with order! > https://groups.google.com/forum/#!searchin/ibtracs-qa/quadrant%7Csort:date/ibtracs-qa/J0wXCeE5PC0/5d2BSykyDQAJ
data.r35se   = (r35(:,2))'; %Watch out with order! > So can you confirm that the 1st is always NE quadrant, 2nd the SE quadrant …
data.r35sw   = (r35(:,3))'; %Watch out with order! > The order is always NE, SE, SW, NW and is true for USA, Reunion and BoM data.
data.r35nw   = (r35(:,4))'; %Watch out with order!
                   
% Read R50
r50     = nc_varget(fname,[isource,'_r50'],[id-1 0 0],[1 Inf Inf]); % in nautical mile
data.r50ne   = r50(:,1)'; %Watch out with order!
data.r50se   = r50(:,2)'; %Watch out with order!
data.r50sw   = r50(:,3)'; %Watch out with order!
data.r50nw   = r50(:,4)'; %Watch out with order!
  
% Read R64
r65     = nc_varget(fname,[isource,'_r64'],[id-1 0 0],[1 Inf Inf]); % in nautical mile
data.r65ne   = r65(:,1)'; %Watch out with order!
data.r65se   = r65(:,2)'; %Watch out with order!
data.r65sw   = r65(:,3)'; %Watch out with order!
data.r65nw   = r65(:,4)'; %Watch out with order!
   
%% 3) set general TC settings 
load 'n:\Deltabox\Bulletin\leijnse\MyBulletin_Github\Matlab\Tropical_Cyclones\tcstructure.mat' % retreived from: p:\11201303-tuvalu-maritime\04_modelsetup\delft3d\TCs\
load 'n:\Deltabox\Bulletin\leijnse\MyBulletin_Github\Matlab\Tropical_Cyclones\spwstucture.mat'

reference_time = datenum(1970,01,01);


%change some values for spiderweb output
spw.cut_off_speed = 0;  %### what are these values ### %
spw.nr_radial_bins = 600;
spw.phi_spiral = 22.6;
spw.asymmetry_option = 'schwerdt1979';
spw.rmax_relation = 'gross2004';
spw.reference_time=datenum(1970,1,1);
spw.merge_frac              = 0.5;

spw.reference_time          = reference_time;
spw.cut_off_speed           = 0;
spw.wind_conversion_factor  = 1; 
spw.radius                  = 900000; % m;   % =900 km

spw.rainfall = 0;

% Length track with data:
%time_active = data.vmax(~isnan(data.vmax));
time_active_ids = find(~isnan(data.vmax));


%% 4) put data in structure
% Conversion 10 to 1 minute
conversion1_10 = 0.93;  % Harper et al., 2010: WMO

% Set dummy values
tc.track = [];
tc.track.time = []; tc.track.rmax = []; 
tc.track.x = []; tc.track.y = [];
tc.track.vmax = []; tc.track.pc = [];
tc.track.r35ne = [] ;    tc.track.r35se =[];     tc.track.r35sw = [];    tc.track.r35nw = [];
tc.track.r50ne = [] ;    tc.track.r50se =[];     tc.track.r50sw = [];    tc.track.r50nw = [];
tc.track.r65ne = [] ;    tc.track.r65se =[];     tc.track.r65sw = [];    tc.track.r65nw = [];
% tc.track.r100ne = [] ;   tc.track.r100se =[];    tc.track.r100sw = [];   tc.track.r100nw = []; % not available for isource = 'usa'

% Create tc
for it = 1:length(time_active_ids) %active length
    time_id = time_active_ids(it)
    
    % Real values
    tc.track.time(it,1)     = data.time(time_id);
    tc.track.x(it,1)        = data.x(time_id);
    tc.track.y(it,1)        = data.y(time_id);
    tc.track.vmax(it,1)     = data.vmax(time_id)*0.514444444*conversion1_10;
    if data.pc(it) > 0
        tc.track.pc(it,1)       = data.pc(time_id);
    else
        tc.track.pc(it,1)       = spw.pn;
    end

    % Determine radius
    tc.track.rmax(it,1)     = data.rmax(time_id);

    % Determine gale force radius as in Nederhof et al. 2019 > in USA we have this data already
    
%     region = 0;
%     probability = 0;
%     [rmax,dr35]  = wind_radii_nederhoff(tc.track.vmax(it,1)/conversion1_10, tc.track.y(it,1), region, probability); % note tcs structure has 10 minute wind, but relationships are derived for 1 minute winds!
%     tc.track.rmax(it,1) = rmax.mode;
%     
%     % as in cyclonestats_write_WES_input.m:
%     if isnan(dr35.mode)
%         r35 = -999;
%     else
%         r35 = rmax.mode + dr35.mode;
%     end
    
    % Starte r35 if possible;
    tc.track.r35ne(it,1) = data.r35ne(time_id);     tc.track.r35se(it,1) = data.r35se(time_id);     tc.track.r35sw(it,1) = data.r35sw(time_id);    tc.track.r35nw(it,1) = data.r35nw(time_id);
    tc.track.r50ne(it,1) = data.r50ne(time_id);     tc.track.r50se(it,1) = data.r50se(time_id);     tc.track.r50sw(it,1) = data.r50sw(time_id);    tc.track.r50nw(it,1) = data.r50nw(time_id);
    tc.track.r65ne(it,1) = data.r65ne(time_id);     tc.track.r65se(it,1) = data.r65se(time_id);     tc.track.r65sw(it,1) = data.r65sw(time_id);    tc.track.r65nw(it,1) = data.r65nw(time_id);
%     tc.track.r100ne(it,1) = data.r100ne(it);   tc.track.r100se(it,1) = data.r100se(it);   tc.track.r100sw(it,1) = data.r100sw(it); tc.track.r100nw(it,1) = data.r100nw(it); % not available for isource = 'usa'
end

%% Optional: do interpolate/extrapolate steps in case there are NaNs in the data
% 
% tc.track.vmax(isnan(tc.track.vmax)) = interp1(tc.track.time(~isnan(tc.track.vmax)),tc.track.vmax(~isnan(tc.track.vmax)),tc.track.time(isnan(tc.track.vmax)),...
%     'nearest','extrap') ;
% 
% tc.track.rmax(isnan(tc.track.rmax)) = interp1(tc.track.time(~isnan(tc.track.rmax)),tc.track.rmax(~isnan(tc.track.rmax)),tc.track.time(isnan(tc.track.rmax)),...
%     'nearest','extrap') ;
% 
% tc.track.r35ne(isnan(tc.track.r35ne)) = interp1(tc.track.time(~isnan(tc.track.r35ne)),tc.track.r35ne(~isnan(tc.track.r35ne)),tc.track.time(isnan(tc.track.r35ne)),...
%     'nearest');%,'extrap') ;
% 
% tc.track.r35se(isnan(tc.track.r35se)) = interp1(tc.track.time(~isnan(tc.track.r35se)),tc.track.r35se(~isnan(tc.track.r35se)),tc.track.time(isnan(tc.track.r35se)),...
%     'nearest');%,'extrap') ;
% 
% tc.track.r35sw(isnan(tc.track.r35sw)) = interp1(tc.track.time(~isnan(tc.track.r35sw)),tc.track.r35sw(~isnan(tc.track.r35sw)),tc.track.time(isnan(tc.track.r35sw)),...
%     'nearest');%,'extrap') ;
% 
% tc.track.r35nw(isnan(tc.track.r35nw)) = interp1(tc.track.time(~isnan(tc.track.r35nw)),tc.track.r35nw(~isnan(tc.track.r35nw)),tc.track.time(isnan(tc.track.r35nw)),...
%     'nearest');%,'extrap') ;

%% last step of input data
tc.wind_speed_unit          = 'ms'; %we already converted from knts to m/s
tc.radius_unit              = 'nm';  %important that this is nautical miles (nm), empirical relationships Kees for TCWiSE are in km
tc.radius_velocity          = [34,50,64,100] * 0.514;
    
%% some figure
cd(destout)
% Check:
addpath(genpath('n:\Deltabox\Bulletin\leijnse\MyBulletin_Github\Matlab\')); %
A4fig; 
subplot(3,1,1)
hold on; box on; grid on;

plot(datetime(tc.track.time, 'convertfrom','datenum'),tc.track.vmax, 'k','linewidth',1)
ylabel('Maximum sustained wind speed [m/s]')

subplot(3,1,2)
hold on; box on; grid on;
plot(datetime(tc.track.time, 'convertfrom','datenum'),tc.track.rmax, 'b','linewidth',1)

plot(datetime(tc.track.time, 'convertfrom','datenum'),tc.track.r35ne, 'r','linewidth',1)
plot(datetime(tc.track.time, 'convertfrom','datenum'),tc.track.r35se, 'r-s','linewidth',1)
plot(datetime(tc.track.time, 'convertfrom','datenum'),tc.track.r35sw, 'r-o','linewidth',1)
plot(datetime(tc.track.time, 'convertfrom','datenum'),tc.track.r35nw, 'r--','linewidth',1)
ylabel('Radius')

legend('rmax','r35ne','r35se','r35sw','r35nw')

subplot(3,1,3)
hold on; box on; grid on;
plot(datetime(tc.track.time, 'convertfrom','datenum'),tc.track.pc, 'r')
ylabel('Pressure drop')

print(['.',filesep,trackname,filesep,'TCinput_variables_',trackname,'.png'],'-dpng','-r500')

[ldbx, ldby] = landboundary('read', 'p:\archivedprojects\11203748-wind_damage\05_TCWiSE\02_application\data\entire_world_coarse.ldb');

%%
A4fig; axis equal; grid on; box on
plot(ldbx, ldby,'k')

plotthick(tc.track.x, tc.track.y)

plot(data.x,data.y,'r.')

scatter(tc.track.x, tc.track.y,[],tc.track.vmax); 
cb = colorbar;
ylabel(cb,'Maximum sustained wind speed [m/s]')

xylim(tc.track.x, tc.track.y)
xylabel_latlon

print(['.',filesep,trackname,filesep,'TCinput_track_',trackname,'.png'],'-dpng','-r500')

%% Make rainfall variations, possible options:
% a) standard ipet
% b) mode bacla, no assymetry, vmax based model
% c) mode bacla, no assymetry, pdef based model 

% not included yet:
% d) percentiles bacla, no assymetry, vmax based model 
% e) percentiles bacla, no assymetry, pdef based model 
% f) random draw bacla, no assymetry, vmax based model 
% g) random draw bacla, no assymetry, pdef based model 
end

%% makes figures of spiderweb

% for jj = 1:length(variations_wanted)
function plot_figures(tc2,spwin, folder,trackname, optionname)

    [ldbx, ldby] = landboundary('read', 'p:\archivedprojects\11203748-wind_damage\05_TCWiSE\02_application\data\entire_world_coarse.ldb');
    
%     spwin = 'P:\11205281-tcwise\03_currentprojects\MSc_Arish\03_modelling\spw_example\HARVEY\cyclone_ipet_symmetric_mode.spw';
%     tc_plot = load('P:\11205281-tcwise\03_currentprojects\MSc_Arish\03_modelling\spw_example\HARVEY\input_tc_ipet_symmetric_mode.mat','tc2');

%     tc_plot = read_spiderweb_file_delft3d(spwin);
    info = asciiwind('open',spwin);

    try
    for id = 1:length(info.Data ) %-20:length(info.Data ) %1:
    [X,Y,UNIT] = asciiwind('grid',info,id ,1:info.Header.n_rows ,1:info.Header.n_cols);

    lon = squeeze(X)';
    lat = squeeze(Y)';
                                                                                                                                                                                                                                                                                                                                                                   
    val = tc2.track(id).precipitation;  %or wind_speed pressure pressure_drop   

    lon = [lon; lon(1,:)];
    lat = [lat; lat(1,:)];
    val = [val; val(end,:)];

    A4fig;axis equal
    plot(ldbx, ldby,'k','linewidth',1)
    pc = pcolor(lon,lat,val); 
    shading flat;
    set(pc, 'facealpha',0.75)
    plot([tc2.track.x], [tc2.track.y],'k','linewidth',2)
    scatter([tc2.track(id).x], [tc2.track(id).y],30,'w','filled')
    
    colormap(jet(100))
    cb = colorbar;
%     caxis([25 85])
    ylabel(cb,'Precipitation rate [mm/hr]')
    xylabel_latlon
    xylim(lon,lat)
%     xlim([-82 -77])
%     ylim([31 34.667])
    title(['Landfall Conditions for ',trackname,'_',optionname,' timestep ',num2str(id,'%02d')],'interpreter','none')
    print([folder,'2D_rain_',trackname,'_',optionname,'_id_',num2str(id,'%02d'),'.png'],'-dpng','-r500')
    end
    end

end
