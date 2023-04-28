% code to calculate various quantities from the raw fowler data
%
% Paul Holland 2/10/13
%

clear all

% ++++++++++++++++++++
% START MAIN POLE LOOP
% ++++++++++++++++++++
for pole=1:2    % 1 for Antarctic, 2 for Arctic
%for pole=1:1
%for pole=2:2

disp(strcat('starting pole=',num2str(pole)))

fname=strcat('fowler_l1_pole',num2str(pole));
load(fname)

% choose months that bracket chosen seasons
%velseas=[1 3; 4 6; 7 9; 10 12;... % JFM AMJ JAS OND
velseas=[2 4; 5 7; 8 10; 11 1;... % FMA MJJ ASO NDJ
         1 12];                   % JFMAMJJASOND

% choose years that bracket periods for mean and trends
velyears=[2003 2010];

% mask all zero data (don't want velocities from areas where no ice exists)
iceuvel(find(iceuvel==0))=NaN;
icevvel(find(icevvel==0))=NaN;
icexvel(find(icexvel==0))=NaN;
iceyvel(find(iceyvel==0))=NaN;

% -----------------------------------------
% calculate means over various time periods
% -----------------------------------------

% make arrays for means
uvelmean=zeros(size(iceuvel,1),size(iceuvel,2),...
               size(velseas,1),size(velyears,1));
vvelmean=uvelmean;

% populate array
tic
for year=1:size(velyears,1)
for seas=1:size(velseas,1)

  % find valid points (alter logic if DJF)
  if (velseas(seas,1)<velseas(seas,2))
    valid=find((veltime(2,:)>=velyears(year,1))&...
               (veltime(2,:)<=velyears(year,2))&...
               (veltime(3,:)>=velseas(seas,1))&...
               (veltime(3,:)<=velseas(seas,2)));
  else
    valid=find((veltime(2,:)>=velyears(year,1))&...
               (veltime(2,:)<=velyears(year,2))&...
               ((veltime(3,:)>=velseas(seas,1))|...
                (veltime(3,:)<=velseas(seas,2))));
  end

  % calculate means
  uvelmean(:,:,seas,year)=nanmean(iceuvel(:,:,valid),3);
  vvelmean(:,:,seas,year)=nanmean(icevvel(:,:,valid),3);

  % mask mean in areas where <30% of possible values exist
  threshold=0.3*size(valid,2);
  for i=1:size(iceuvel,1)
  for j=1:size(iceuvel,2)
    if size(find(~isnan(iceuvel(i,j,valid)))) < threshold
      uvelmean(i,j,seas,year)=NaN;
      vvelmean(i,j,seas,year)=NaN;
    end
  end
  end

end
end
toc

% ------------------------------------------------
% calculate trend by year for various time periods
% ------------------------------------------------

%
% calculate yearly means
%

minyr=min(veltime(2,:));
maxyr=max(veltime(2,:));
numyr=maxyr-minyr+1;

uvelyearmean=zeros(size(iceuvel,1),size(iceuvel,2),numyr,size(velseas,1));
vvelyearmean=uvelyearmean;

% populate array
tic
for year=minyr:maxyr
for seas=1:size(velseas,1)

  % find valid points (alter logic if DJF)
  if (velseas(seas,1)<velseas(seas,2))
    valid=find((veltime(2,:)==year)&...
               (veltime(3,:)>=velseas(seas,1))&...
               (veltime(3,:)<=velseas(seas,2)));
  else
    valid=find((veltime(2,:)==year)&...
               ((veltime(3,:)>=velseas(seas,1))|...
                (veltime(3,:)<=velseas(seas,2))));
  end

  % calculate mean
  yindex=year-minyr+1;
  uvelyearmean(:,:,yindex,seas)=nanmean(iceuvel(:,:,valid),3);
  vvelyearmean(:,:,yindex,seas)=nanmean(icevvel(:,:,valid),3);

end
end
toc

%
% calculate trend
%

uveltrend=zeros(size(iceuvel,1),size(iceuvel,2),...
                size(velseas,1),size(velyears,1));
vveltrend=uveltrend;
uvelconf=uveltrend;
vvelconf=uveltrend;

tic
for year=1:size(velyears,1)

  % make array of times (in years)
  ymindex=velyears(year,1)-minyr+1;
  ymaxdex=velyears(year,2)-minyr+1;
  numyr=velyears(year,2)-velyears(year,1)+1;
  years=[1:numyr]';
  years=[ones(size(years)) years];
  threshold=max(0.90*numyr,2);

  % for each x,y point for each season, regress against years
  for seas=1:size(velseas,1)
    nx=size(uveltrend,1);
    ny=size(uveltrend,2);
    parfor i=1:nx
    %for i=1:nx
    for j=1:ny

      uvelij=zeros(numyr,1);
      uvelij(:)=squeeze(uvelyearmean(i,j,ymindex:ymaxdex,seas));
      valid=find(~isnan(uvelij));

      % only calculate trends where enough data present
      if (size(valid,1)>=threshold)

        [regcoefs,confints,r,rint,stats]=regress(uvelij,years);
        uveltrend(i,j,seas,year)=regcoefs(2);
        uvelconf(i,j,seas,year)=100*(1-stats(3));

        vvelij=zeros(numyr,1);
        vvelij(:)=squeeze(vvelyearmean(i,j,ymindex:ymaxdex,seas));
        [regcoefs,confints,r,rint,stats]=regress(vvelij,years);
        vveltrend(i,j,seas,year)=regcoefs(2);
        vvelconf(i,j,seas,year)=100*(1-stats(3));

      end

    end
    end
  end

end
toc

% mask all zero data
uveltrend(find(uveltrend==0))=NaN;
vveltrend(find(vveltrend==0))=NaN;

% ------------
% lat and long
% ------------

vellat=vellat;
vellon=vellon;

% ----
% save
% ----

fname=strcat('fowler_l2_pole',num2str(pole));
save(fname,...
          'vellat','vellon', ...
          'velseas','velyears', ...
          'uvelmean','vvelmean', ...
          'uveltrend','vveltrend', ...
          'uvelconf','vvelconf', ...
     '-v7.3');

% ++++++++++++++++++
% END MAIN POLE LOOP
% ++++++++++++++++++
end

