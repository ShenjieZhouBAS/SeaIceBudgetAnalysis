% code to read daily Fowler motion vectors and save in lat-lon form
%
% Paul Holland 24/6/15
%
close all
clear all

% ++++++++++++++++++++
% START MAIN POLE LOOP
% ++++++++++++++++++++
for pole=1:1    % 1 for Antarctic, 2 for Arctic
%for pole=1:1
%for pole=2:2

disp(strcat('starting pole=',num2str(pole)))

if pole==1
  froot='/Users/Ale/Desktop/dottorato Hobart/data (cruise summer 14-15)/data/sea_ice_motion_satellite/nsidc/daily/';
  fpref='icemotion_daily_sh_25km_';
  fsuff='_v4.1.nc';
  fgrid=...
      '/Users/Ale/Desktop/dottorato Hobart/data (cruise summer 14-15)/data/sea_ice_motion_satellite/nsidc/daily/coord_sea_ice_motion.txt';
  nx=321;
  slat=-70;
  cmer=0;
  
end

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

numrec=length([datenum('1 Jan 2018'):1:datenum('31 Dec 2018')]);

% create counter and output arrays
%

count=0;
count_1=0;
dysnstrt=0;
icexvel=zeros(nx,nx,numrec);
iceyvel=zeros(nx,nx,numrec);
veltime=zeros(5,numrec);

for yr=yrstart:yrend

  clear a b
  count_2=0;
  if (mod(yr,4)==0)
    dyend=[31 29 31 30 31 30 31 31 30 31 30 31];
  else
    dyend=[31 28 31 30 31 30 31 31 30 31 30 31];
  end

  dyoyr=0;

  disp(strcat('ice motion year-',num2str(yr)))
  
  fname=strcat(froot,fpref,num2str(yr),fsuff);
  a=ncread(fname,'u');
  b=ncread(fname,'v');  
  a=flipdim(a,3);
  b=flipdim(b,3);
  
  for mn=mnstart:mnend
  for dy=dystart:dyend(mn)

    dysnstrt=dysnstrt+1;
    dyoyr=dyoyr+1;
    %fid=fopen(fname,'r','ieee-le');
    
    %if (fid>0)
      count=count+1;
      count_1=count_1+1;
      count_2=count_2+1;

      %fclose(fid);
      %icevel=cat(3,a(:,:,count_2),b(:,:,count_2));
      %icevel=reshape(icevel,[2 nx nx]);
      %icevel=flip(icevel,3);
      %icexvel(:,:,count)=squeeze(icevel(1,:,:))';
      %iceyvel(:,:,count)=squeeze(icevel(2,:,:))';

      icexvel(:,:,count_1)=a(:,:,count_2)';
      iceyvel(:,:,count_1)=b(:,:,count_2)';

      % add times to master array
      veltime(:,count_1)=[dysnstrt yr mn dy dyoyr];

    %end

  end
  end

end



%
% convert data into desired form
%
icexvel(find(icexvel==0))=NaN;
iceyvel(find(iceyvel==0))=NaN;
icexvel=icexvel/100;                        % convert to m/s
iceyvel=iceyvel/100;

% ----------------------
% latitude and longitude
% ----------------------

%vellon=ncread(fname,'longitude');
%vellat=ncread(fname,'latitude');

fid=fopen(fgrid);
latlon=fscanf(fid,'%f',[4 nx*nx]);
fclose(fid);

vellat=reshape(latlon(3,:),size(icexvel,1),size(icexvel,2));
vellon=reshape(latlon(4,:),size(icexvel,1),size(icexvel,2));
vellon=flipud(vellon');        % -180 to +180 version

% ------------------------------------------
% calculate velocities in lat-lon directions
% ------------------------------------------

iceuvel=zeros(size(icexvel));
icevvel=zeros(size(icexvel));
[ivelx,ively]=psll2xy(vellat,vellon,slat,cmer);

[ivelxplusonedeglat,ivelyplusonedeglat]=psll2xy(vellat+1,vellon,slat,cmer);
metperdeglat=sqrt((ivelxplusonedeglat-ivelx).^2+...
                 (ivelyplusonedeglat-ively).^2);

[ivelxplusonedeglon,ivelyplusonedeglon]=psll2xy(vellat,vellon+1,slat,cmer);
metperdeglon=sqrt((ivelxplusonedeglon-ivelx).^2+...
                 (ivelyplusonedeglon-ively).^2);

delt=100;

tic
for t=1:size(icexvel,3)

  % get end points in polar stereo from velocities (m)
  endx=ivelx+icexvel(:,:,t)*delt;
  endy=ively+iceyvel(:,:,t)*delt;

  % convert end points to lat-lon
  [endlat,endlon]=psxy2ll(endx,endy,slat,cmer);

  % calculate lat-lon displacements and convert to metres and then velocities
  iceuvel(:,:,t)=(endlon-vellon).*metperdeglon/delt;
  icevvel(:,:,t)=(endlat-vellat).*metperdeglat/delt;

end

%iceuvel(find(abs(iceuvel)>10^3))=NaN;
%icevvel(find(abs(iceuvel)>10^3))=NaN;
%iceuvel(find(abs(icevvel)>10^3))=NaN;
%icevvel(find(abs(icevvel)>10^3))=NaN;

toc

% ---------------------
% sub-sample if desired
% ---------------------

disp('WARNING: sub-sampling Fowler') 
skip=2

dx=skip*25000;
indices=[ceil(skip/2):skip:nx];
vellat=vellat(indices,indices);
vellon=vellon(indices,indices);
iceuvel=iceuvel(indices,indices,:);
icevvel=icevvel(indices,indices,:);
icexvel=icexvel(indices,indices,:);
iceyvel=iceyvel(indices,indices,:);

% ----
% save
% ----

fname=strcat('fowler_2018_l1_pole',num2str(pole));
save(fname,...
           'veltime','vellat','vellon','dx', ...
                     'iceuvel','icevvel', ...
                     'icexvel','iceyvel', ...
     '-v7.3');

% ++++++++++++++++++
% END MAIN POLE LOOP
% ++++++++++++++++++
end





