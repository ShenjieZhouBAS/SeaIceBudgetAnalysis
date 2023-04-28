% code to calculate various quantities from the raw bootstrap data
%
% Paul Holland 29/5/12
%
close all
clear all

% ++++++++++++++++++++
% START MAIN POLE LOOP
% ++++++++++++++++++++
pole=1; %Antarctica



fname=strcat('bootstrap_2013_2017_l1_pole',num2str(pole));
load(fname)

% choose months that bracket chosen seasons
%concseas=[1 3; 4 6; 7 9; 10 12;... % JFM AMJ JAS OND
concseas=[2 4; 5 7; 8 10; 11 1;... % FMA MJJ ASO NDJ
          1 12];                   % JFMAMJJASOND

% choose years that bracket periods for mean and trends
concyears=[2013 2017];

% -----------------------------------------
% calculate means over various time periods
% -----------------------------------------

% make arrays for means
concmean=zeros(size(iceconc,1),size(iceconc,2),...
               size(concseas,1),size(concyears,1));
intensmean=zeros(size(iceconc,1),size(iceconc,2),...
               size(concseas,1),size(concyears,1));

% populate array
tic
for year=1:size(concyears,1)
for seas=1:size(concseas,1)

  % find valid points (alter logic if DJF)
  if (concseas(seas,1)<concseas(seas,2))
    valid=find((conctime(2,:)>=concyears(year,1))&...
               (conctime(2,:)<=concyears(year,2))&...
               (conctime(3,:)>=concseas(seas,1))&...
               (conctime(3,:)<=concseas(seas,2)));
  else
    valid=find((conctime(2,:)>=concyears(year,1))&...
               (conctime(2,:)<=concyears(year,2))&...
               ((conctime(3,:)>=concseas(seas,1))|...
                (conctime(3,:)<=concseas(seas,2))));
  end

  % calculate means
  concmean(:,:,seas,year)=nanmean(iceconc(:,:,valid),3);
  intensmean(:,:,seas,year)=nanmean(iceintens(:,:,valid),3);

end
end
toc

% mask areas with mean less than 0.15
concmean(find(concmean<0.15))=NaN;
intensmean(find(isnan(concmean)))=NaN;

% ------------------------------------------------
% calculate trend by year for various time periods
% ------------------------------------------------

%
% calculate yearly means
%

minyr=min(conctime(2,:));
maxyr=max(conctime(2,:));
numyr=maxyr-minyr+1;

concyearmean=zeros(size(iceconc,1),size(iceconc,2),numyr,size(concseas,1));
intensyearmean=zeros(size(iceconc,1),size(iceconc,2),numyr,size(concseas,1));

% populate array
tic
for year=minyr:maxyr
for seas=1:size(concseas,1)

  % find valid points (alter logic if DJF)
  if (concseas(seas,1)<concseas(seas,2))
    valid=find((conctime(2,:)==year)&...
               (conctime(3,:)>=concseas(seas,1))&...
               (conctime(3,:)<=concseas(seas,2)));
  else
    valid=find((conctime(2,:)==year)&...
               ((conctime(3,:)>=concseas(seas,1))|...
                (conctime(3,:)<=concseas(seas,2))));
  end

  % calculate means
  yindex=year-minyr+1;
  concyearmean(:,:,yindex,seas)=nanmean(iceconc(:,:,valid),3);
  intensyearmean(:,:,yindex,seas)=nanmean(iceintens(:,:,valid),3);

end
end
toc

%
% calculate trend
%

conctrend=zeros(size(iceconc,1),size(iceconc,2),...
                size(concseas,1),size(concyears,1));
concconf=zeros(size(iceconc,1),size(iceconc,2),...
               size(concseas,1),size(concyears,1));
intenstrend=zeros(size(iceconc,1),size(iceconc,2),...
                  size(concseas,1),size(concyears,1));
intensconf=zeros(size(iceconc,1),size(iceconc,2),...
                 size(concseas,1),size(concyears,1));

tic
for year=1:size(concyears,1)

  % make array of times (in years)
  ymindex=concyears(year,1)-minyr+1;
  ymaxdex=concyears(year,2)-minyr+1;
  numyr=concyears(year,2)-concyears(year,1)+1;
  years=[1:numyr]';
  years=[ones(size(years)) years];

  % for each x,y point for each season, regress against years
  for seas=1:size(concseas,1)
    nx=size(conctrend,1);
    ny=size(conctrend,2);
    parfor i=1:nx
    for j=1:ny

      concij=zeros(numyr,1);
      concij(:)=squeeze(concyearmean(i,j,ymindex:ymaxdex,seas));
      valid=find(~isnan(concij));

      if (size(valid,2)>=2)

        [regcoefs,confints,r,rint,stats]=regress(concij,years);
        conctrend(i,j,seas,year)=regcoefs(2);
        concconf(i,j,seas,year)=100*(1-stats(3));

        intensij=zeros(numyr,1);
        intensij(:)=squeeze(intensyearmean(i,j,ymindex:ymaxdex,seas));
        [regcoefs,confints,r,rint,stats]=regress(intensij,years);
        intenstrend(i,j,seas,year)=regcoefs(2);
        intensconf(i,j,seas,year)=100*(1-stats(3));

      end

    end
    end
  end

end
toc

% mask areas with mean concentration masked
conctrend(find(isnan(concmean)))=NaN;
concconf(find(isnan(concmean)))=NaN;
intenstrend(find(isnan(concmean)))=NaN;
intensconf(find(isnan(concmean)))=NaN;

% ----
% save
% ----

fname=strcat('bootstrap_2013_2017_l2_pole',num2str(pole));
save(fname,...
          'conclat','conclon', ...
          'concseas','concyears', ...
          'concmean','intensmean', ...
          'conctrend','intenstrend', ...
          'concconf','intensconf', ...
     '-v7.3');

% ++++++++++++++++++
% END MAIN POLE LOOP
% ++++++++++++++++++



