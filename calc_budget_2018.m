% code to calculate various quantities from both kimura and nasa team data
%
% Paul Holland 22/1/15
%

close all
clear all

tic

% ++++++++++++++++++++
% START MAIN POLE LOOP
% ++++++++++++++++++++

pole=1; %Antarctica

disp(strcat('starting pole=',num2str(pole)))



%disp('bootstrap/kimura: takes 10 mins for both poles')
%load(strcat('kimura_2013_2017_l1_pole',num2str(pole)));
%load(strcat('bootstrap_kimuragrid_2013_2017_nanmean_l1_pole',num2str(pole)));

load(strcat('kimura_2018_l1_pole',num2str(pole)));
load(strcat('amsre_kimuragrid_2018_nanmean_l1_pole',num2str(pole)));

outfname=strcat('budget_siv_smoothed_kimura_2018_l1_pole1_amsre');

% choose months that bracket chosen seasons
%budseas=[1 3; 4 6; 7 9; 10 12;... % JFM AMJ JAS OND
budseas=[2 4; 5 7; 8 10; 11 1;... % FMA MJJ ASO NDJ
         1 12];                   % JFMAMJJASOND

% choose years that bracket period for mean and trends
budyears=[2018 2018];
numyears=budyears(end)-budyears(1)+1;

% use velocity grid
budlat=vellat;
budlon=vellon;

% ------------------------------------------------------
% find times where both velocity and concentration exist
% ------------------------------------------------------

common=[];
for t=1:size(veltime,2)
  index=find(conctime(1,:)==veltime(1,t));
  if (~isempty(index))
    common=[common; veltime(1,t)];
  end
end

budtime=zeros(5,size(common,1));
budconc=zeros(size(icexvel,1),size(icexvel,2),size(common,1));
budxvel=budconc;
budyvel=budconc;

for t=1:size(common,1)

  cindex=find(conctime(1,:)==common(t));
  vindex=find(veltime(1,:)==common(t));

  budtime(:,t)=conctime(:,cindex);
  budconc(:,:,t)=concvelgrid(:,:,cindex);

  budxvel(:,:,t)=icexvel(:,:,vindex);
  budyvel(:,:,t)=iceyvel(:,:,vindex);

end

for t=2:size(budtime,2)
  dt=budtime(1,t)-budtime(1,t-1);
  if (dt>1)
    disp(strcat('WARNING: exceptional budtime dt of length:',num2str(dt))) 
  end
end

% test effect of not blanking land/ocean (note will affect smoothing below)
%budconc(find(isnan(budconc)))=0;
%budxvel(find(isnan(budxvel)))=0;
%budyvel(find(isnan(budyvel)))=0;

% ---------------------------
% smooth input data spatially
% ---------------------------

% smooth over 2*smrad+1 cells

concsmits=1;

concsmrad=0;
%concsmrad=1;
%concsmrad=2;
%concsmrad=3;

velsmits=1;
%velsmits=2;
%velsmits=4;

%velsmrad=0;
%velsmrad=1;
%velsmrad=2;
velsmrad=3;
%velsmrad=4;

disp('WARNING: smoothing input fields spatially - beware variable dx')
%outfname=strcat(outfname,'_csmspat',num2str(concsmrad));
%outfname=strcat(outfname,'_vsmspat',num2str(velsmrad));

for s=1:concsmits
parfor t=1:size(budconc,3)
  budconc(:,:,t)=smooth2a(budconc(:,:,t),concsmrad);
end
end

for s=1:velsmits
parfor t=1:size(budconc,3)
  budxvel(:,:,t)=smooth2a(budxvel(:,:,t),velsmrad);
  budyvel(:,:,t)=smooth2a(budyvel(:,:,t),velsmrad);
end
end

% ----------------------------
% smooth input data temporally
% ----------------------------

% smooth over 2*smrad+1 cells

concsmrad=0;
%concsmrad=1;
%concsmrad=2;
%concsmrad=3;

velsmrad=0;
%velsmrad=1;
%velsmrad=2;
%velsmrad=3;

disp('WARNING: smoothing input fields temporally - beware variable dt')
%outfname=strcat(outfname,'_csmtemp',num2str(concsmrad));
%outfname=strcat(outfname,'_vsmtemp',num2str(velsmrad));

ny=size(budconc,2);
parfor i=1:size(budconc,1)
for j=1:ny
  budconc(i,j,:)=smooth2a(squeeze(budconc(i,j,:)),concsmrad);
end
end

ny=size(budconc,2);
parfor i=1:size(budconc,1)
for j=1:ny
  budxvel(i,j,:)=smooth2a(squeeze(budxvel(i,j,:)),velsmrad);
  budyvel(i,j,:)=smooth2a(squeeze(budyvel(i,j,:)),velsmrad);
end
end

% ------------------------------
% calculate unsteady term (/s)
% (central differencing in time)
% ------------------------------

disp('WARNING: calculating unsteady term - beware variable dt')

buddif=zeros(size(budconc));

buddif(:,:,1)=(budconc(:,:,2)-budconc(:,:,1))/(86400);
for t=2:size(budconc,3)-1
  buddif(:,:,t)=(budconc(:,:,t+1)-budconc(:,:,t-1))/(2*86400);
end
buddif(:,:,end)=(budconc(:,:,end)-budconc(:,:,end-1))/(86400);

% ---------------------------------------------------------
% concentration source from advection (/s)
% (conc and vel are co-located in space, staggered in time)
% ---------------------------------------------------------

budadv=zeros(size(budconc));

xveltmp=zeros(size(budconc));
yveltmp=zeros(size(budconc));

%% get velocities on concentration time-stamps
%disp(strcat('WARNING: nanmean vel to conc time'));
%xveltmp(:,:,1)=budxvel(:,:,1);
%yveltmp(:,:,1)=budyvel(:,:,1);
%for t=2:size(budconc,3)
%  xveltmp(:,:,t)=nanmean(budxvel(:,:,t-1:t),3);
%  yveltmp(:,:,t)=nanmean(budyvel(:,:,t-1:t),3);
%end

% just use daily velocities
disp(strcat('WARNING: ignoring time-stamp difference between conc and vel'));
xveltmp=budxvel;
yveltmp=budyvel;

% A-grid advection terms
% (correct, but would be clearer if i and j were labelled the other way round)
ny=size(budconc,2);
parfor i=2:size(budconc,1)-1
for j=2:ny-1
 budadv(i,j,:)=-xveltmp(i,j,:).*((budconc(i,j+1,:)-budconc(i,j-1,:))/(2*dx))...
               -yveltmp(i,j,:).*((budconc(i+1,j,:)-budconc(i-1,j,:))/(2*dx));
end
end

% ---------------------------------------------------------
% concentration source from divergence (/s)
% (conc and vel are co-located in space, staggered in time)
% ---------------------------------------------------------

buddiv=zeros(size(budconc));

% A-grid divergence terms
% (correct, but would be clearer if i and j were labelled the other way round)
ny=size(budconc,2);
parfor i=2:size(budconc,1)-1
for j=2:ny-1
  buddiv(i,j,:)=-budconc(i,j,:).*((xveltmp(i,j+1,:)-xveltmp(i,j-1,:))/(2*dx))...
                -budconc(i,j,:).*((yveltmp(i+1,j,:)-yveltmp(i-1,j,:))/(2*dx));
end
end

% -----------------------------------------------------------------------------
% smooth advection and divergence over 3 time steps to match centred difference
% -----------------------------------------------------------------------------

disp('WARNING: nanmean adv and div over 3 timesteps - beware variable dt')

budadvtemp=zeros(size(budadv));
buddivtemp=zeros(size(buddiv));

budadvtemp(:,:,1)=nanmean(budadv(:,:,1:2),3);
buddivtemp(:,:,1)=nanmean(buddiv(:,:,1:2),3);
for t=2:size(budadv,3)-1
  budadvtemp(:,:,t)=nanmean(budadv(:,:,t-1:t+1),3);
  buddivtemp(:,:,t)=nanmean(buddiv(:,:,t-1:t+1),3);
end
budadvtemp(:,:,end)=nanmean(budadv(:,:,end-1:end),3);
buddivtemp(:,:,end)=nanmean(buddiv(:,:,end-1:end),3);

budadv=budadvtemp;
buddiv=buddivtemp;

% -----------------------
% calculate residual (/s)
% -----------------------

budres=buddif-budadv-buddiv;

% --------------------------------------------
% find thermodynamic fraction of residual (/s)
% -------------------------------------------

% residual must be thermodynamics if any of these is true
% - divergent flow
% - low concentration
% this is almost all freezing because divergence doesn't seem to be compensated 
% by concentration decrease or advection
% (possibly a sampling issue in divergence or concentration difference)

% in high-concentration, convergent flow
% - positive residual (conc increase) implies freezing exceeds ridging/melting
%   (value of residual is then lower bound estimate of freezing)
% - negative residual implies melting/ridging exceeds freezing
%   (a bit of freezing is missed here)


%% OPTION 1: label residual as thermodynamics where certain on dynamic grounds
%% this results in a lower bound estimate for total freezing
%% remainder is melt plus ridging; this is an upper bound estimate for both
%disp('WARNING: partitioning residual using div and conc only')
%validtmd=find((icedivg>=0)|(budconc<0.85)); % this is definitely thermo
%validrid=find((icedivg<0)&(budconc>=0.85)); % some of this could be thermo
%                                  % (mostly melting) and thus offset the above

%% OPTION 2: as above, but add freezing in convergent case
%disp('WARNING: partitioning residual using div and conc and freez in conv');
%validtmd=find((icedivg>=0)|(budconc<0.85)|(budres>=0));
%validrid=find((icedivg<0)&(budconc>=0.85)&(budres<0));

% OPTION 3: separate all positive and negative parts of residual
% this results in an upper-bound estimate of the total >smoothing but <season
% compensation between freezing and melting/ridging
disp('WARNING: partitioning residual into all positive and negative');
validtmd=find(budres>=0);
validrid=find(budres<0);

budtmd=zeros(size(budconc));
budrid=zeros(size(budconc));

budtmd(validtmd)=budres(validtmd);
budrid(validrid)=budres(validrid);

% ------------------------------------
% other fields for diagnostic purposes
% ------------------------------------

% concentration
budcnc=budconc;

% speed
budspd=sqrt(budxvel.^2+budyvel.^2);

% velocity divergence
budvdv=zeros(size(budconc));
% correct, but would be clearer if i and j were labelled the other way round
ny=size(budconc,2);
parfor i=2:size(budconc,1)-1
for j=2:ny-1
  budvdv(i,j,:)=((budxvel(i,j+1,:)-budxvel(i,j-1,:))/(2*dx))...
                +((budyvel(i+1,j,:)-budyvel(i-1,j,:))/(2*dx));
end
end

% mod of grad concentration
budgrc=zeros(size(budconc));
% correct, but would be clearer if i and j were labelled the other way round
ny=size(budconc,2);
parfor i=2:size(budconc,1)-1
for j=2:ny-1
  budgrc(i,j,:)=sqrt(((budconc(i,j+1,:)-budconc(i,j-1,:))/(2*dx)).^2 ...
                    +((budconc(i+1,j,:)-budconc(i-1,j,:))/(2*dx)).^2);
end
end

% -------------------------------
% mask all fields for consistency
% -------------------------------

budfdf=buddif;  % save a copy of unmasked diff and conc for comparison
budfcn=budcnc;  

% mask with divergence (conservative choice as has smallest footprint)
disp('WARNING: masking all daily terms with divergence')
%outfname=strcat(outfname,'_daydivmask');
valid=find(isnan(buddiv));

%% do not mask (incorrect but pretty)
%disp('WARNING: not masking daily terms')
%%outfname=strcat(outfname,'_daynomask');
%valid=[];

buddif(valid)=NaN;
budadv(valid)=NaN;
buddiv(valid)=NaN;
budres(valid)=NaN;
budtmd(valid)=NaN;
budrid(valid)=NaN;

budedg=budfdf;  % save a copy of masked bits (presumed ice edge) for comparison
budedg(find(~isnan(buddif)))=NaN;

% mask diagnostics
budcnc(valid)=NaN;
budspd(valid)=NaN;
budvdv(valid)=NaN;
budgrc(valid)=NaN;

% --------------------------------------------------
% calculate means for all seasons for all years (/y)
% --------------------------------------------------

budfdfseas=zeros(size(budconc,1),size(budconc,2),numyears,size(budseas,1));
buddifseas=zeros(size(budconc,1),size(budconc,2),numyears,size(budseas,1));
budedgseas=zeros(size(budconc,1),size(budconc,2),numyears,size(budseas,1));
budadvseas=zeros(size(budconc,1),size(budconc,2),numyears,size(budseas,1));
buddivseas=zeros(size(budconc,1),size(budconc,2),numyears,size(budseas,1));
budresseas=zeros(size(budconc,1),size(budconc,2),numyears,size(budseas,1));
budtmdseas=zeros(size(budconc,1),size(budconc,2),numyears,size(budseas,1));
budridseas=zeros(size(budconc,1),size(budconc,2),numyears,size(budseas,1));
budfcnseas=zeros(size(budconc,1),size(budconc,2),numyears,size(budseas,1));
budcncseas=zeros(size(budconc,1),size(budconc,2),numyears,size(budseas,1));
budspdseas=zeros(size(budconc,1),size(budconc,2),numyears,size(budseas,1));
budvdvseas=zeros(size(budconc,1),size(budconc,2),numyears,size(budseas,1));
budgrcseas=zeros(size(budconc,1),size(budconc,2),numyears,size(budseas,1));

for seas=1:size(budseas,1)
for year=1:numyears

  % find valid points (alter logic if DJF)
  if (budseas(seas,1)<budseas(seas,2))
    valid=find((budtime(2,:)==(budyears(1)+year-1))&...
               (budtime(3,:)>=budseas(seas,1))&(budtime(3,:)<=budseas(seas,2)));
  else
    valid=find( ( (budtime(2,:)==(budyears(1)+year-1)) ...
                & (budtime(3,:)>=budseas(seas,1)) ) ...
              | ( (budtime(2,:)==(budyears(1)+year  )) ...
                & (budtime(3,:)<=budseas(seas,2)) ) ...
              );
    if (year==1)
      disp('WARNING: a season is including months from subsequent year')
    end
  end

  % budget: sum up all valid daily values, and convert to per year
  % (nansum means finite output; missing data treated as zero)
  budfdfseas(:,:,year,seas)=nansum(budfdf(:,:,valid),3)*86400*4;
  buddifseas(:,:,year,seas)=nansum(buddif(:,:,valid),3)*86400*4;
  budedgseas(:,:,year,seas)=nansum(budedg(:,:,valid),3)*86400*4;
  budadvseas(:,:,year,seas)=nansum(budadv(:,:,valid),3)*86400*4;
  buddivseas(:,:,year,seas)=nansum(buddiv(:,:,valid),3)*86400*4;
  budresseas(:,:,year,seas)=nansum(budres(:,:,valid),3)*86400*4;
  budtmdseas(:,:,year,seas)=nansum(budtmd(:,:,valid),3)*86400*4;
  budridseas(:,:,year,seas)=nansum(budrid(:,:,valid),3)*86400*4;
  if (seas==1)&(year==1)
    disp('WARNING: seasonmean from nansum with div mask => missing=0')
    %outfname=strcat(outfname,'_seasnansum');
  end

  %% budget: sum up all valid daily values, and convert to per year
  %% (nanmean means finite output; missing data treated as mean; deprecated)
  %numdays=size(valid,2);
  %budfdfseas(:,:,year,seas)=nanmean(budfdf(:,:,valid),3)*numdays*86400*4;
  %buddifseas(:,:,year,seas)=nanmean(buddif(:,:,valid),3)*numdays*86400*4;
  %budedgseas(:,:,year,seas)=nanmean(budedg(:,:,valid),3)*numdays*86400*4;
  %budadvseas(:,:,year,seas)=nanmean(budadv(:,:,valid),3)*numdays*86400*4;
  %buddivseas(:,:,year,seas)=nanmean(buddiv(:,:,valid),3)*numdays*86400*4;
  %budresseas(:,:,year,seas)=nanmean(budres(:,:,valid),3)*numdays*86400*4;
  %budtmdseas(:,:,year,seas)=nanmean(budtmd(:,:,valid),3)*numdays*86400*4;
  %budridseas(:,:,year,seas)=nanmean(budrid(:,:,valid),3)*numdays*86400*4;
  %if (seas==1)&(year==1)
  %  disp('WARNING: seasonmean from nanmean with div mask => missing=mean')
  %  %outfname=strcat(outfname,'_seasnanmean');
  %end

  % diagnostics: sum up all valid values and divide by number of valid days
  % (nansum means finite output; missing data treated as zero)
  numdays=size(valid,2);
  budfcnseas(:,:,year,seas)=nansum(budfcn(:,:,valid),3)/numdays;
  budcncseas(:,:,year,seas)=nansum(budcnc(:,:,valid),3)/numdays;
  budspdseas(:,:,year,seas)=nansum(budspd(:,:,valid),3)/numdays;
  budvdvseas(:,:,year,seas)=nansum(budvdv(:,:,valid),3)/numdays;
  budgrcseas(:,:,year,seas)=nansum(budgrc(:,:,valid),3)/numdays;

end
end

% -----------------------------------------------
% calculate interannual mean for each season (/y)
% -----------------------------------------------

budfdfseasmean=zeros(size(budconc,1),size(budconc,2),size(budseas,1));
buddifseasmean=zeros(size(budconc,1),size(budconc,2),size(budseas,1));
budedgseasmean=zeros(size(budconc,1),size(budconc,2),size(budseas,1));
budadvseasmean=zeros(size(budconc,1),size(budconc,2),size(budseas,1));
buddivseasmean=zeros(size(budconc,1),size(budconc,2),size(budseas,1));
budresseasmean=zeros(size(budconc,1),size(budconc,2),size(budseas,1));
budtmdseasmean=zeros(size(budconc,1),size(budconc,2),size(budseas,1));
budridseasmean=zeros(size(budconc,1),size(budconc,2),size(budseas,1));
budfcnseasmean=zeros(size(budconc,1),size(budconc,2),size(budseas,1));
budcncseasmean=zeros(size(budconc,1),size(budconc,2),size(budseas,1));
budspdseasmean=zeros(size(budconc,1),size(budconc,2),size(budseas,1));
budvdvseasmean=zeros(size(budconc,1),size(budconc,2),size(budseas,1));
budgrcseasmean=zeros(size(budconc,1),size(budconc,2),size(budseas,1));

% calculate interannual mean from all available data
for seas=1:size(budseas,1)

  % find valid points (alter logic if DJF)
  if (budseas(seas,1)<budseas(seas,2))
    valid=find((budtime(3,:)>=budseas(seas,1))&(budtime(3,:)<=budseas(seas,2)));
  else
    valid=find((budtime(3,:)>=budseas(seas,1))|(budtime(3,:)<=budseas(seas,2)));
  end

  % budget: sum up all valid daily values, and convert to per year
  % (nansum means finite output; missing data treated as zero)
  budfdfseasmean(:,:,seas)=nansum(budfdf(:,:,valid),3)*86400*4/numyears;
  buddifseasmean(:,:,seas)=nansum(buddif(:,:,valid),3)*86400*4/numyears;
  budedgseasmean(:,:,seas)=nansum(budedg(:,:,valid),3)*86400*4/numyears;
  budadvseasmean(:,:,seas)=nansum(budadv(:,:,valid),3)*86400*4/numyears;
  buddivseasmean(:,:,seas)=nansum(buddiv(:,:,valid),3)*86400*4/numyears;
  budresseasmean(:,:,seas)=nansum(budres(:,:,valid),3)*86400*4/numyears;
  budtmdseasmean(:,:,seas)=nansum(budtmd(:,:,valid),3)*86400*4/numyears;
  budridseasmean(:,:,seas)=nansum(budrid(:,:,valid),3)*86400*4/numyears;
  if (seas==1)
    disp('WARNING: taking interann from nansum with div mask => missing=0')
    %outfname=strcat(outfname,'_intannnansum');
  end

  % diagnostics: sum up all valid values and divide by number of valid days
  % (nansum means finite output; missing data treated as zero)
  numdays=size(valid,2);
  budfcnseasmean(:,:,seas)=nansum(budfcn(:,:,valid),3)/numdays;
  budcncseasmean(:,:,seas)=nansum(budcnc(:,:,valid),3)/numdays;
  budspdseasmean(:,:,seas)=nansum(budspd(:,:,valid),3)/numdays;
  budvdvseasmean(:,:,seas)=nansum(budvdv(:,:,valid),3)/numdays;
  budgrcseasmean(:,:,seas)=nansum(budgrc(:,:,valid),3)/numdays;

end

% -------------------------------------------------
% calculate means for all months for all years (/y)
% -------------------------------------------------

budfdfmonth=zeros(size(budconc,1),size(budconc,2),numyears,12);
buddifmonth=zeros(size(budconc,1),size(budconc,2),numyears,12);
budedgmonth=zeros(size(budconc,1),size(budconc,2),numyears,12);
budadvmonth=zeros(size(budconc,1),size(budconc,2),numyears,12);
buddivmonth=zeros(size(budconc,1),size(budconc,2),numyears,12);
budresmonth=zeros(size(budconc,1),size(budconc,2),numyears,12);
budtmdmonth=zeros(size(budconc,1),size(budconc,2),numyears,12);
budridmonth=zeros(size(budconc,1),size(budconc,2),numyears,12);
budfcnmonth=zeros(size(budconc,1),size(budconc,2),numyears,12);
budcncmonth=zeros(size(budconc,1),size(budconc,2),numyears,12);
budspdmonth=zeros(size(budconc,1),size(budconc,2),numyears,12);
budvdvmonth=zeros(size(budconc,1),size(budconc,2),numyears,12);
budgrcmonth=zeros(size(budconc,1),size(budconc,2),numyears,12);

for month=1:12
for year=1:numyears

  % find valid points
  valid=find((budtime(2,:)==(budyears(1)+year-1))&...
             (budtime(3,:)==month));

  % budget: sum up all valid daily values, and convert to per year
  % (nansum means finite output; missing data treated as zero)
  budfdfmonth(:,:,year,month)=nansum(budfdf(:,:,valid),3)*86400*12;
  buddifmonth(:,:,year,month)=nansum(buddif(:,:,valid),3)*86400*12;
  budedgmonth(:,:,year,month)=nansum(budedg(:,:,valid),3)*86400*12;
  budadvmonth(:,:,year,month)=nansum(budadv(:,:,valid),3)*86400*12;
  buddivmonth(:,:,year,month)=nansum(buddiv(:,:,valid),3)*86400*12;
  budresmonth(:,:,year,month)=nansum(budres(:,:,valid),3)*86400*12;
  budtmdmonth(:,:,year,month)=nansum(budtmd(:,:,valid),3)*86400*12;
  budridmonth(:,:,year,month)=nansum(budrid(:,:,valid),3)*86400*12;
  if (month==1)&(year==1)
    disp('WARNING: month from nansum with div mask => missing=0')
    %outfname=strcat(outfname,'_monthnansum');
  end

  % diagnostics: sum up all valid values and divide by number of valid days
  % (nansum means finite output; missing data treated as zero)
  numdays=size(valid,2);
  budfcnmonth(:,:,year,month)=nansum(budfcn(:,:,valid),3)/numdays;
  budcncmonth(:,:,year,month)=nansum(budcnc(:,:,valid),3)/numdays;
  budspdmonth(:,:,year,month)=nansum(budspd(:,:,valid),3)/numdays;
  budvdvmonth(:,:,year,month)=nansum(budvdv(:,:,valid),3)/numdays;
  budgrcmonth(:,:,year,month)=nansum(budgrc(:,:,valid),3)/numdays;

end
end

% ----------------------------------------------
% calculate interannual mean for each month (/y)
% ----------------------------------------------

budfdfmonthmean=zeros(size(budconc,1),size(budconc,2),12);
buddifmonthmean=zeros(size(budconc,1),size(budconc,2),12);
budedgmonthmean=zeros(size(budconc,1),size(budconc,2),12);
budadvmonthmean=zeros(size(budconc,1),size(budconc,2),12);
buddivmonthmean=zeros(size(budconc,1),size(budconc,2),12);
budresmonthmean=zeros(size(budconc,1),size(budconc,2),12);
budtmdmonthmean=zeros(size(budconc,1),size(budconc,2),12);
budridmonthmean=zeros(size(budconc,1),size(budconc,2),12);
budfcnmonthmean=zeros(size(budconc,1),size(budconc,2),12);
budcncmonthmean=zeros(size(budconc,1),size(budconc,2),12);
budspdmonthmean=zeros(size(budconc,1),size(budconc,2),12);
budvdvmonthmean=zeros(size(budconc,1),size(budconc,2),12);
budgrcmonthmean=zeros(size(budconc,1),size(budconc,2),12);

% calculate interannual mean from all available data
for month=1:12

  % find valid points
  valid=find(budtime(3,:)==month);

  % budget: sum up all valid daily values, convert to /y
  % (nansum means finite output; missing data treated as zero)
  budfdfmonthmean(:,:,month)=nansum(budfdf(:,:,valid),3)*86400*12/numyears;
  buddifmonthmean(:,:,month)=nansum(buddif(:,:,valid),3)*86400*12/numyears;
  budedgmonthmean(:,:,month)=nansum(budedg(:,:,valid),3)*86400*12/numyears;
  budadvmonthmean(:,:,month)=nansum(budadv(:,:,valid),3)*86400*12/numyears;
  buddivmonthmean(:,:,month)=nansum(buddiv(:,:,valid),3)*86400*12/numyears;
  budresmonthmean(:,:,month)=nansum(budres(:,:,valid),3)*86400*12/numyears;
  budtmdmonthmean(:,:,month)=nansum(budtmd(:,:,valid),3)*86400*12/numyears;
  budridmonthmean(:,:,month)=nansum(budrid(:,:,valid),3)*86400*12/numyears;

  % diagnostics: sum up all valid values and divide by number of valid days
  % (nansum means finite output; missing data treated as zero)
  numdays=size(valid,2);
  budfcnmonthmean(:,:,month)=nansum(budfcn(:,:,valid),3)/numdays;
  budcncmonthmean(:,:,month)=nansum(budcnc(:,:,valid),3)/numdays;
  budspdmonthmean(:,:,month)=nansum(budspd(:,:,valid),3)/numdays;
  budvdvmonthmean(:,:,month)=nansum(budvdv(:,:,valid),3)/numdays;
  budgrcmonthmean(:,:,month)=nansum(budgrc(:,:,valid),3)/numdays;

end

% -------------------------------------------------------
% calculate trend by year for various time periods (/y/y)
% -------------------------------------------------------

budfdftrend=zeros(size(budconc,1),size(budconc,2),size(budseas,1));
budfdfconf=zeros(size(budconc,1),size(budconc,2),size(budseas,1));
budrestrend=zeros(size(budconc,1),size(budconc,2),size(budseas,1));
budresconf=zeros(size(budconc,1),size(budconc,2),size(budseas,1));
budfdftrend(:,:,:)=NaN;
budfdfconf(:,:,:)=NaN;
budrestrend(:,:,:)=NaN;
budresconf(:,:,:)=NaN;

disp('WARNING: removed trend calcs to save time')
%{

% make array of times (in years)
years=[1:numyears]';
years=[ones(size(years))  years];

% for each x,y point for each season, regress valid velocity against years
for seas=1:size(budseas,1)
  nx=size(budconc,1);
  ny=size(budconc,2);
  parfor i=1:nx
  for j=1:ny

    % if desired, only calculate trends when fraction of residual values exist
    threshold=0.0*numyears;

    budresij=zeros(numyears,1);
    budresij(:)=squeeze(budresseas(i,j,:,seas));
    valid=find(~isnan(budresij));

    if (size(valid,1)>=2)&(size(valid,1)>threshold)

      [regcoefs,confints,r,rint,stats]=regress(budresij,years);
      budrestrend(i,j,seas)=regcoefs(2);
      budresconf(i,j,seas)=100*(1-stats(3));

      budfdfij=zeros(numyears,1);
      budfdfij(:)=squeeze(budfdfseas(i,j,:,seas));
      [regcoefs,confints,r,rint,stats]=regress(budfdfij,years);
      budfdftrend(i,j,seas)=regcoefs(2);
      budfdfconf(i,j,seas)=100*(1-stats(3));
    
    end

  end
  end

end

%}

% -----------------
% smooth all fields
% -----------------

% do I want to do this now?

% ---------------
% mask all fields
% ---------------

% mask with divergence field (has smallest stencil)
disp('WARNING: masking final output with divergence')
%outfname=strcat(outfname,'_divmask');
valid=find((buddivseas==0)|(isnan(buddivseas)));

%% mask with difference field (has largest stencil)
%disp('WARNING: masking final output with difference')
%%outfname=strcat(outfname,'_difmask');
%valid=find((buddifseas==0)|(isnan(buddifseas)));

budfdfseas(valid)=NaN;
buddifseas(valid)=NaN;
budedgseas(valid)=NaN;
budadvseas(valid)=NaN;
buddivseas(valid)=NaN;
budresseas(valid)=NaN;
budtmdseas(valid)=NaN;
budridseas(valid)=NaN;
budfcnseas(valid)=NaN;
budcncseas(valid)=NaN;
budspdseas(valid)=NaN;
budvdvseas(valid)=NaN;
budgrcseas(valid)=NaN;

valid=find((buddivseasmean==0)|(isnan(buddivseasmean)));
%valid=find((buddifseasmean==0)|(isnan(buddifseasmean)));

budfdfseasmean(valid)=NaN;
buddifseasmean(valid)=NaN;
budedgseasmean(valid)=NaN;
budadvseasmean(valid)=NaN;
buddivseasmean(valid)=NaN;
budresseasmean(valid)=NaN;
budtmdseasmean(valid)=NaN;
budridseasmean(valid)=NaN;
budfcnseasmean(valid)=NaN;
budcncseasmean(valid)=NaN;
budspdseasmean(valid)=NaN;
budvdvseasmean(valid)=NaN;
budgrcseasmean(valid)=NaN;

budfdftrend(valid)=NaN;
budrestrend(valid)=NaN;
budfdfconf(valid)=NaN;
budresconf(valid)=NaN;

valid=find((buddivmonth==0)|(isnan(buddivmonth)));
%valid=find((buddifmonth==0)|(isnan(buddifmonth)));

budfdfmonth(valid)=NaN;
buddifmonth(valid)=NaN;
budedgmonth(valid)=NaN;
budadvmonth(valid)=NaN;
buddivmonth(valid)=NaN;
budresmonth(valid)=NaN;
budtmdmonth(valid)=NaN;
budridmonth(valid)=NaN;
budfcnmonth(valid)=NaN;
budcncmonth(valid)=NaN;
budspdmonth(valid)=NaN;
budvdvmonth(valid)=NaN;
budgrcmonth(valid)=NaN;

valid=find((buddivmonthmean==0)|(isnan(buddivmonthmean)));
%valid=find((buddifmonthmean==0)|(isnan(buddifmonthmean)));

budfdfmonthmean(valid)=NaN;
buddifmonthmean(valid)=NaN;
budedgmonthmean(valid)=NaN;
budadvmonthmean(valid)=NaN;
buddivmonthmean(valid)=NaN;
budresmonthmean(valid)=NaN;
budtmdmonthmean(valid)=NaN;
budridmonthmean(valid)=NaN;
budfcnmonthmean(valid)=NaN;
budcncmonthmean(valid)=NaN;
budspdmonthmean(valid)=NaN;
budvdvmonthmean(valid)=NaN;
budgrcmonthmean(valid)=NaN;

% ----
% save
% ----

save(outfname,...
  'budlat','budlon','dx', ...
  'budseas','budyears', ...
  'budfdfseas','budfdfseasmean','budfdfmonth','budfdfmonthmean', ...
               'budfdftrend','budfdfconf', ...
  'buddifseas','buddifseasmean','buddifmonth','buddifmonthmean','buddif', ...
  'budedgseas','budedgseasmean','budedgmonth','budedgmonthmean', ...
  'budadvseas','budadvseasmean','budadvmonth','budadvmonthmean','budadv', ...
  'buddivseas','buddivseasmean','buddivmonth','buddivmonthmean','buddiv', ...
  'budresseas','budresseasmean','budresmonth','budresmonthmean','budres', ...
               'budrestrend','budresconf', ...
  'budtmdseas','budtmdseasmean','budtmdmonth','budtmdmonthmean', ...
  'budridseas','budridseasmean','budridmonth','budridmonthmean', ...
  'budfcnseas','budfcnseasmean','budfcnmonth','budfcnmonthmean', ...
  'budcncseas','budcncseasmean','budcncmonth','budcncmonthmean','budcnc',...
  'budspdseas','budspdseasmean','budspdmonth','budspdmonthmean', ...
  'budvdvseas','budvdvseasmean','budvdvmonth','budvdvmonthmean', ...
  'budgrcseas','budgrcseasmean','budgrcmonth','budgrcmonthmean', ...
  '-v7.3');

% ++++++++++++++++++
% END MAIN POLE LOOP
% ++++++++++++++++++


toc

%%
%% save data for David if desired
%% (note need to run above code for pole 1 only)
%%
%
%save('david',...
%  'budlat','budlon','budtime',...
%  'buddif','budadv','buddiv','budres',...
%  'buddifmonth','budadvmonth','buddivmonth','budresmonth',...
%  'buddifseasmean','budadvseasmean','buddivseasmean','budresseasmean',...
%  '-v7.3');

