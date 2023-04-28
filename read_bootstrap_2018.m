% code to read all daily Bootstrap ice concentrations and save to matlab file
% with lat-lon and time information in the same file
%
% Paul Holland 02/02/16
%

close all
clear all

% ++++++++++++++++++++
% START MAIN POLE LOOP
% ++++++++++++++++++++

pole=1; % analysis in Antarctica

froot=strcat('/Users/Ale/Desktop/dottorato Hobart/data (cruise summer 14-15)/data/SIC_satellite/bootstrap_COMISO_daily/');
fsuff='_v3.1_s.bin';
flatgrid='/Users/Ale/Desktop/dottorato Hobart/data (cruise summer 14-15)/data/SIC_satellite/bootstrap_COMISO_daily/pss25lats_v3.dat';
flongrid='/Users/Ale/Desktop/dottorato Hobart/data (cruise summer 14-15)/data/SIC_satellite/bootstrap_COMISO_daily/pss25lons_v3.dat';
nx=332;
ny=316;
slat=-70;
cmer=0;



% ---------
% read data
% ---------
 
% dysnstart is used to merge the datasets, so all these need to be common to 
% all data reading
yrstart=2018;
yrend=2018;
mnstart=01;
mnend=12;
dystart=01;
%dyend=[31 29 31 30 31 30 31 31 30 31 30 31];  % set below

%
% first determine the number of valid data files so can make output arrays
% (much quicker than using cat to add to arrays)
%

numrec=0;
for yr=yrstart:yrend

  if (mod(yr,4)==0)
    dyend=[31 29 31 30 31 30 31 31 30 31 30 31];
  else
    dyend=[31 28 31 30 31 30 31 31 30 31 30 31];
  end

  for mn=mnstart:mnend
  for dy=dystart:dyend(mn)

    % build filename
    time=1e4*yr+1e2*mn+dy;
    if time < 19870801
       sensor='n07';
    elseif time < 19911203
       sensor='f08';
    elseif time < 19951001
       sensor='f11';
    elseif time < 20080101
       sensor='f13';
    else
       sensor='f17';
    end

    fdate=sprintf('%04d%02d%02d',yr,mn,dy);
    fname=strcat(froot,num2str(yr),'/bt_',fdate,'_',sensor,fsuff);

    % find out if data exists
    fid=fopen(fname,'r');
    if (fid>0)
      numrec=numrec+1;
      fclose(fid);
    end
  end 
  end
end



%
% create counter and output arrays
%
count=0;
dysnstrt=0;
iceconc=zeros(nx,ny,numrec);
conctime=zeros(5,numrec);

for yr=yrstart:yrend

  if (mod(yr,4)==0)
    dyend=[31 29 31 30 31 30 31 31 30 31 30 31];
  else
    dyend=[31 28 31 30 31 30 31 31 30 31 30 31];
  end

  dyoyr=0;


  for mn=mnstart:mnend
  for dy=dystart:dyend(mn)

    dysnstrt=dysnstrt+1;
    dyoyr=dyoyr+1;

    % build filename
    time=1e4*yr+1e2*mn+dy;
    if time < 19870801
       sensor='n07';
    elseif time < 19911203
       sensor='f08';
    elseif time < 19951001
       sensor='f11';
    elseif time < 20080101
       sensor='f13';
    else
       sensor='f17';
    end

    fdate=sprintf('%04d%02d%02d',yr,mn,dy);
    fname=strcat(froot,num2str(yr),'/bt_',fdate,'_',sensor,fsuff);

    % read data if it exists
    fid=fopen(fname,'r');
    if (fid>0)

      count=count+1;

      conc=fread(fid,'int16');
      fclose(fid);

      % shuffle data around into Bedmap-like array
      conc=reshape(conc,ny,nx);
      conc=shiftdim(conc,1);
      conc=conc(nx:-1:1,1:ny);

      % remove non-sea ice data
      conc(find(conc>1000))=NaN;

      % rescale data into concentrations between 0 and 1
      conc=0.001*conc;

      % add field to master array
      iceconc(:,:,count)=conc;
      
      % add times to master array
      conctime(:,count)=[dysnstrt yr mn dy dyoyr];
    end

  end
  end
  
end



% ----------------------
% latitude and longitude
% ----------------------

fid=fopen(flatgrid);
conclat=fread(fid,'int32','ieee-le');
fclose(fid);
conclat=conclat/1e5; % convert to decimal degrees
conclat=reshape(conclat,size(iceconc,2),size(iceconc,1))';
conclat=flipud(conclat);

fid=fopen(flongrid);
conclon=fread(fid,'int32','ieee-le');
fclose(fid);
conclon=conclon/1e5; % convert to decimal degrees
conclon=reshape(conclon,size(iceconc,2),size(iceconc,1))';
conclon=flipud(conclon);


% ------------------------------------
% make daily fields of intensification
% (/y)
% ------------------------------------

iceintens=zeros(size(iceconc));

for t=1:size(iceconc,3)-1
  iceintens(:,:,t)=(iceconc(:,:,t+1)-iceconc(:,:,t))...
                  /((conctime(1,t+1)-conctime(1,t))/365.25);
end

% ----
% save
% ----

fname=strcat('bootstrap_2018_l1_pole',num2str(pole));
save(fname,...
           'conctime','conclat','conclon', ...
                      'iceconc','iceintens', ...
     '-v7.3');

% ++++++++++++++++++
% END MAIN POLE LOOP
% ++++++++++++++++++





%%

lat_min=-80;
lat_max=-60;
lon_min=100;
lon_max=220;

load bedmap2coast.mat
count=0;
for i=1:218

    for j=1:length(coastlon{1})
        
        count=count+1;
        lat_coast(count)=coastlat{1}(j);
        lon_coast(count)=coastlon{1}(j);
        
    end
    
end

clear tmp tmp1 tmp2
tmp = find(lon_coast<=0);
lon_coast(tmp) = lon_coast(tmp)+360;
tmp1 = find(lon_coast==0); 
lon_coast(tmp1) = 0.0001;
tmp2 = find(lon_coast==360); 
lon_coast(tmp2) = 355.9999;

clear tmp tmp1 tmp2
tmp = find(conclon<=0);
conclon(tmp) = conclon(tmp)+360;
tmp1 = find(conclon==0); 
conclon(tmp1) = 0.0001;
tmp2 = find(conclon==360); 
conclon(tmp2) = 355.9999;

figure
m_proj('mercator','longitude',[lon_min lon_max],'latitude',[lat_min lat_max]);
hold on
m_grid
hold on
clear X Y Z
[X Y]=m_ll2xy(conclon,conclat);
Z=iceconc(:,:,180);
hold on
scatter(reshape(X,[332*316,1]),...
    reshape(Y,[332*316,1]),50,reshape(Z,[332*316,1]),'filled')
caxis([0 1])
colormap('jet')
colorbar
hold on
m_plot(lon_coast,lat_coast,'.k')
hold on
m_patch(lon_coast,lat_coast,[0.4 0.4 0.4]) % grounded ice in dark grey
m_grid
set(gca,'XTickLabel',[],'YTickLabel',[])


