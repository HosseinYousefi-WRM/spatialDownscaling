% Contact: ho_yousefi@sbu.ac.ir
% ERA5 Netcdf Downscaling

clear
clc

% % --------------------input----------------------------------------------

% % Specify the directory containing the .nc files
folderPath = 'D:\OMAN\MAIN Folder\ERA5'; % Replace with your directory path
% Start and end dates
% startDate = datetime(1984, 1, 1); % January 1, 2025
% endDate = datetime(2014, 12, 31); % December 31, 2025

observationStations=readtable('observationStations.xlsx');
OBS_lon_lat=table2array(observationStations(:,4:5));

% % -----------------------------------------------------------------------


% PART 1
% %--------------------gridPoints -----------------------------------------
% % Extract the names of the files
ncFiles = dir(fullfile(folderPath, '*.nc'));
fileNames = {ncFiles.name};

% % Extract grid points
% ncdisp(fileNames{1});
lon = ncread(fileNames{1},'longitude')';
lat = ncread(fileNames{1},'latitude')';
c_lon=numel(lon);
c_lat=numel(lat);
NC_gridPoint_lon_lat(1:c_lon,2)=lat(1);
NC_gridPoint_lon_lat(1:c_lon,1)=lon';
i=2;
while i<=c_lat
    x(1:c_lon,2)=lat(i);
    x(1:c_lon,1)=lon';
    NC_gridPoint_lon_lat=vertcat(NC_gridPoint_lon_lat,x);
    i=i+1;
end

save('NC_gridPoint','NC_gridPoint_lon_lat');


% PART 2
% %--------------------sorting files based on year-------------------------
disp("If you are running this code for the first time you should change sort_index to 'y' for sorting files based on year.")
sort_index= 'n';%input('Do you want to sort files based on year (y/n)?  ','s');
if sort_index == 'y'
    for i=1:size(fileNames,2)
        year=ncreadatt(fileNames{i},"/valid_time","units");
        year_num(i)=str2double(year(12:15));
    end
    
    for i=1:size(fileNames,2)
        newName = strcat(num2str(year_num(i)),'-',fileNames{i}); % Replace with the desired file name
        % Rename the file
        status(i) = movefile(fileNames{i}, newName);
    end
end

% % Extract the new names of the files
ncFiles = dir(fullfile(folderPath, '*.nc'));
newFileNames = {ncFiles.name};


% PART 2
% %--------------------spatialDownscaling----------------------------------

% OBS_lon_lat
% NC_gridPoint_lon_lat

for i=1:size(OBS_lon_lat,1)
    for j=1:size(NC_gridPoint_lon_lat,1)
        Station{i,j} = [OBS_lon_lat(i,1),OBS_lon_lat(i,2);NC_gridPoint_lon_lat(j,1),NC_gridPoint_lon_lat(j,2)];
        distance(i,j) = pdist(Station{i,j},'euclidean');
    end
end


for i=1:size(OBS_lon_lat,1)
    Y = sort(distance(i,:));
    val = Y(1:4);
    for gp=1:4
        Sur_gP_ID(i,gp)=find(distance(i,:)==val(gp));
        Sur_gP_Distance(i,gp)=distance(i,Sur_gP_ID(i,gp));
    end
end


IDW_weights=1./Sur_gP_Distance;



% PART 3
% %--------------------Precipitation (output in mm)------------------------
for st=1:size(OBS_lon_lat,1)
    
    for i=1:size(newFileNames,2)
        Precip{i} = ncread(newFileNames{i},'tp');
    end
    
    ID=Sur_gP_ID(st,:);%[517:520,558:561,599:601];% input('Enter the grid point ID (the row number of the nearest grid points to the ground stations, e.g. [558,559,600]):');%559
    
    for id=1:numel(ID)
        if mod(ID(id),c_lon)==0
            column=floor(ID(id)/c_lon);
            row=c_lon;
        else
            column=floor(ID(id)/c_lon)+1;
            row=mod(ID(id),c_lon);
        end
        
        for y=1:size(Precip,2)
            P=Precip{y}(row,column,:);
            PCP{y}=reshape(P,[numel(P),1,1])*1000; % *1000 to convert from m to mm
            clear P;
        end
        
        
        Precipitation=PCP{1};
        i=2;
        while i<=size(Precip,2)
            xx=PCP{i};
            Precipitation=vertcat(Precipitation,xx);
            i=i+1;
        end
    
        % save_name=strcat('PCP_gP_',num2str(ID(id)),'.mat');
        % save(save_name,'Precipitation');
    
        if id==1
            NC_totalPrecipitation{st}(:,1)=Precipitation;
        else
            NC_totalPrecipitation{st}=horzcat(NC_totalPrecipitation{st},Precipitation);
        end
    end
    
    % save('NC_totalPrecipitation','NC_totalPrecipitation');
    W=repmat(IDW_weights(st,:),[size(NC_totalPrecipitation{st},1),1]);
    Sum_W=sum(W,2);

    IDW_newStation(:,st)=sum(NC_totalPrecipitation{st}.*W,2)./Sum_W;
end

%-----------------------------Writing results--------------------------------------

writematrix(IDW_newStation,'ERA5_Downscaled_Precipitation_for_Obs.xlsx');
