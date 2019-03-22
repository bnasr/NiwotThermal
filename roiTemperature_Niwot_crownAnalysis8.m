% Code updated to pull any temperature for a given time (i.e. pulls IRR
% data even if image is missing). Code also gives snow and fog flags based
% on input from Caffe analysis on Niwot4 Phenocam. Default snow and fog
% flag is false (i.e. missing Phenocam data is not noted). ROIs
% automatically shift for the major camera shift experienced in summer
% 2016, after the camera was removed and then reinstalled. Downwelling
% longwave is taken from Sean Burns Niwot climate data (column 28).

clear
clc

codeDir = '/Users/donald/Documents/Research Data/Code/Matlab Code/FLIR Image Processing';
imgDir = {'/Users/donald/Documents/Research Data/FLIR Camera/All Data/FLIR3 - 2015 Master Collection',...
    '/Users/donald/Documents/Research Data/FLIR Camera/All Data/FLIR3 - 2016 Master Collection',...
    '/Users/donald/Documents/Research Data/FLIR Camera/All Data/FLIR3 - 2017 Master Collection'};
yy = [2015; 2016; 2017];
currentYear = 2017; % most recent year that has not had the datetimestamp corrected to local standard time
loggerFile1 = '/Users/donald/Documents/Research Data/Data Loggers/Niwot/CR1000_stats.dat';
loggerFile2 = '/Users/donald/Documents/Research Data/Data Loggers/Niwot/CR800_stats.dat';
loggerFile3 = '/Users/donald/Documents/Research Data/Data Loggers/Niwot/climate_2015_ver.2016.06.09.dat';
loggerFile4 = '/Users/donald/Documents/Research Data/Data Loggers/Niwot/climate_2016_ver.2017.09.15.dat';
loggerFile5 = '/Users/donald/Documents/Research Data/Data Loggers/Niwot/climate_2017_ver.2017.09.25.dat';
snowFile = '/Users/donald/Documents/Research Data/Niwot projects/snow data analysis/snowFlagData.txt';
distanceMapFile = '/Users/donald/Documents/Research Data/Niwot projects/2015-09-15 pointcloud/20150915_niwotPixelDistances.mat';
roiFile = {'niwotflir_F2.mat';'niwotflir_F3.mat';'niwotflir_P2.mat';'niwotflir_P3.mat';...
    'niwotflir_S1.mat';'niwotflir_S3.mat';'niwotflir_referencePlate.mat'};
objectEmissivity = 0.97;
applyCorrections = true;
minSolarElevation = 0;
minAirTemp = 0;
nearestNeighborImages = 7;
shiftFracYear = [0; 2016.743169398907];
shiftX = [0; 5];
shiftY = [0; 3];

startDir = pwd;

for ii = 1:numel(roiFile)
    for jj = 1:numel(shiftFracYear)
        load(roiFile{ii});
        mask{ii,jj} = poly2mask(vertices(:,1)+shiftX(jj),vertices(:,2)+shiftY(jj),480,640);
    end
end

load(distanceMapFile);
pixelDist = pixelDistFinal;
clear pixelDistFinal

fid = fopen(loggerFile1);
loggerOne = textscan(fid,'%s%d%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','Delimiter',',','HeaderLines',4,'TreatAsEmpty','"NAN"');
fclose(fid);
clear fid
fid = fopen(loggerFile2);
loggerTwo = textscan(fid,'%s%d%f%f%f%f%f%f%f%f%f%f%f%f','Delimiter',',','HeaderLines',4,'TreatAsEmpty','"NAN"');
fclose(fid);
clear fid
fid = fopen(snowFile);
snowData = textscan(fid,'%s%s%s%f%s%f%f','Delimiter',',','HeaderLines',1);
fclose(fid);
clear fid
fid = fopen(loggerFile3);
loggerThree = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','Delimiter','\t','HeaderLines',61,'TreatAsEmpty','"NaN"');
fclose(fid);
clear fid
fid = fopen(loggerFile4);
loggerFour = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','Delimiter','\t','HeaderLines',61,'TreatAsEmpty','"NaN"');
fclose(fid);
clear fid
fid = fopen(loggerFile5);
loggerFive = textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','Delimiter','\t','HeaderLines',61,'TreatAsEmpty','"NaN"');
fclose(fid);
clear fid

loggerOne_DT = loggerOne{1};
loggerOne_airT = loggerOne{7};
loggerOne_RH = loggerOne{8};
loggerOne_plateT = loggerOne{21};
loggerOne_skyT = loggerOne{22};
loggerOne_canopyT1 = loggerOne{23};
loggerTwo_DT = loggerTwo{1};
loggerTwo_canopyT2 = loggerTwo{13};
loggerTwo_canopyT3 = loggerTwo{14};
snowData_filename = snowData{2};
snowData_DT = snowData{3};
snowData_solarElevation = snowData{4};
snowData_snowy = snowData{5};
snowData_snowValue = snowData{6};
snowData_treesValue = snowData{7};
loggerThree_year = loggerThree{1};
loggerThree_month = loggerThree{2};
loggerThree_day = loggerThree{3};
loggerThree_hour = loggerThree{4};
loggerThree_minute = loggerThree{5};
loggerThree_second = loggerThree{6};
loggerThree_downwellingLW = loggerThree{28};
loggerFour_year = loggerFour{1};
loggerFour_month = loggerFour{2};
loggerFour_day = loggerFour{3};
loggerFour_hour = loggerFour{4};
loggerFour_minute = loggerFour{5};
loggerFour_second = loggerFour{6};
loggerFour_downwellingLW = loggerFour{28};
loggerFive_year = loggerFive{1};
loggerFive_month = loggerFive{2};
loggerFive_day = loggerFive{3};
loggerFive_hour = loggerFive{4};
loggerFive_minute = loggerFive{5};
loggerFive_second = loggerFive{6};
loggerFive_downwellingLW = loggerFive{28};
loggerUSNR1_year = [loggerThree_year; loggerFour_year; loggerFive_year];
loggerUSNR1_month = [loggerThree_month; loggerFour_month; loggerFive_month];
loggerUSNR1_day = [loggerThree_day; loggerFour_day; loggerFive_day];
loggerUSNR1_hour = [loggerThree_hour; loggerFour_hour; loggerFive_hour];
loggerUSNR1_minute = [loggerThree_minute; loggerFour_minute; loggerFive_minute];
loggerUSNR1_second = [loggerThree_second; loggerFour_second; loggerFive_second];
loggerUSNR1_downwellingLW = [loggerThree_downwellingLW; loggerFour_downwellingLW; loggerFive_downwellingLW];

clear loggerOne loggerTwo snowData loggerThree loggerFour loggerFive loggerThree_year...
    loggerThree_month loggerThree_day loggerThree_hour loggerThree_minute...
    loggerThree_second loggerThree_downwellingLW loggerFour_year loggerFour_month...
    loggerFour_day loggerFour_hour loggerFour_minute loggerFour_second...
    loggerFour_downwellingLW loggerFive_year loggerFive_month loggerFive_day...
    loggerFive_hour loggerFive_minute loggerFive_second loggerFive_downwellingLW

cd(codeDir)

loggerOneFracYear = zeros(numel(loggerOne_DT),1);
loggerOneJD = zeros(numel(loggerOne_DT),1);
for ii = 1:numel(loggerOne_DT)
    datetime = loggerOne_DT{ii};
    datetime = datetime(2:end-1);
    parts = regexp(datetime,' ','split');
    date = parts{1};
    time = parts{2};
    clear parts
    parts = regexp(date,'-','split');
    year = str2double(parts{1});
    month = str2double(parts{2});
    day = str2double(parts{3});
    clear parts
    parts = regexp(time,':','split');
    hour = str2double(parts{1});
    minute = str2double(parts{2});
    second = str2double(parts{3});
    clear parts
    loggerOneFracYear(ii) = date2fracyear(year,month,day,hour,minute,second);
    loggerOneJD(ii) = date2jd(year,month,day,hour,minute,second);
end
clear datetime date time year month day hour minute second
loggerTwoFracYear = zeros(numel(loggerTwo_DT),1);
loggerTwoJD = zeros(numel(loggerTwo_DT),1);
for ii = 1:numel(loggerTwo_DT)
    datetime = loggerTwo_DT{ii};
    datetime = datetime(2:end-1);
    parts = regexp(datetime,' ','split');
    date = parts{1};
    time = parts{2};
    clear parts
    parts = regexp(date,'-','split');
    year = str2double(parts{1});
    month = str2double(parts{2});
    day = str2double(parts{3});
    clear parts
    parts = regexp(time,':','split');
    hour = str2double(parts{1});
    minute = str2double(parts{2});
    second = str2double(parts{3});
    clear parts
    loggerTwoFracYear(ii) = date2fracyear(year,month,day,hour,minute,second);
    loggerTwoJD(ii) = date2jd(year,month,day,hour,minute,second);
end
clear datetime date time year month day hour minute second
snowDataFracYear = zeros(numel(snowData_DT),1);
snowDataJD = zeros(numel(snowData_DT),1);
snowData_snowy1 = false(numel(snowData_DT),1);
for ii = 1:numel(snowData_DT)
    datetime=snowData_DT{ii};
    parts = regexp(datetime,' ','split');
    date = parts{1};
    time = parts{2};
    clear parts
    parts = regexp(date,'-','split');
    year = str2double(parts{1});
    month = str2double(parts{2});
    day = str2double(parts{3});
    clear parts
    parts = regexp(time,':','split');
    hour = str2double(parts{1});
    minute = str2double(parts{2});
    second = str2double(parts{3});
    clear parts
    snowDataFracYear(ii) = date2fracyear(year,month,day,hour,minute,second);
    snowDataJD(ii) = date2jd(year,month,day,hour,minute,second);
    snowData_snowy1(ii,1) = strcmpi(snowData_snowy{ii},'TRUE');
end
snowData_snowy = snowData_snowy1;
clear snowData_snowy1
loggerUSNR1FracYear = zeros(numel(loggerUSNR1_year),1);
loggerUSNR1JD = zeros(numel(loggerUSNR1_year),1);
loggerUSNR1_skyTemp = zeros(numel(loggerUSNR1_year),1);
for ii = 1:numel(loggerUSNR1_year)
    year = loggerUSNR1_year(ii,1);
    month = loggerUSNR1_month(ii,1);
    day = loggerUSNR1_day(ii,1);
    hour = loggerUSNR1_hour(ii,1);
    minute = loggerUSNR1_minute(ii,1);
    second = loggerUSNR1_second(ii,1);
    loggerUSNR1FracYear(ii) = date2fracyear(year,month,day,hour,minute,second);
    loggerUSNR1JD(ii) = date2jd(year,month,day,hour,minute,second);
    loggerUSNR1_skyTemp(ii) = (loggerUSNR1_downwellingLW(ii)/0.0000000567)^(0.25)-273.15;
end
clear datetime date time year month day hour minute second

for ii = 1:numel(yy)
    
    w = waitbar(0,['Processing filenames from directory ',num2str(ii),' of ',num2str(numel(imgDir)),'.']);
    
    processFiles = dir(strcat(imgDir{ii},filesep,'*.seq'));
    
    for jj = 1:numel(processFiles)
        filename = processFiles(jj).name;
        filenameShort = filename(1:end-4);
        parts = regexp(filenameShort,'-','split');
        dt = parts{3};
        clear parts
        parts = regexp(dt,'_','split');
        doy = str2double(parts{1});
        hour = str2double(parts{2});
        minute = str2double(parts{3});
        second = str2double(parts{4})+str2double(parts{5})/1000;
        clear parts
        [~,month,day] = jd2date(yy(ii),doy);
        if (yy(ii) == currentYear)
            jd(jj,1) = date2jd(yy(ii),month,day,hour,minute,second) - (7/24);
        elseif (yy(ii) ~= str2double(datestr(now,'YYYY'))) % everything but current year has had filenames corrected
            jd(jj,1) = date2jd(yy(ii),month,day,hour,minute,second);
        else
            jd(jj,1) = date2jd(yy(ii),month,day,hour,minute,second) - (7/24);
        end
        [~,M,D] = jd2date(yy(ii),floor(jd(jj,1)));
        [hh,mm,ss] = fracDay2time(jd(jj,1)-floor(jd(jj,1)));
        fracYear(jj,1) = date2fracyear(yy(ii),M,D,hh,mm,ss);
        imgDateTime{jj,1} = ['"',num2str(yy(ii)),'-',num2str(M,'%0.2d'),'-',num2str(D,'%0.2d'),...
            ' ',num2str(hh,'%0.2d'),':',num2str(mm,'%0.2d'),':',num2str(ss,'%0.2d'),'"'];
        clear filename filenameShort dt doy hour minute second month day
        
        % Update waitbar every 500 datapoints
        if (mod(jj,500) == 0)
            waitbar(jj/numel(processFiles),w,['Processing filenames from directory ',num2str(ii),' of ',num2str(numel(imgDir)),'.']);
        end
        
    end
    close(w);
    
    loggerOneCurrentYear = loggerOneFracYear((loggerOneFracYear >= yy(ii)) & (loggerOneFracYear <= yy(ii)+1));
    loggerTwoCurrentYear = loggerTwoFracYear((loggerTwoFracYear >= yy(ii)) & (loggerTwoFracYear <= yy(ii)+1));
    loggerUSNR1CurrentYear = loggerUSNR1FracYear((loggerUSNR1FracYear >= yy(ii)) & (loggerUSNR1FracYear <= yy(ii)+1));
    uniqueFracYear = unique([fracYear;loggerOneCurrentYear;loggerTwoCurrentYear]);
    
    w = waitbar(0,['Working on directory ',num2str(ii),' of ',num2str(numel(imgDir)),', datapoint 0']);
    
    for jj = 1:numel(uniqueFracYear)
        
        [~,M,D,hh,mm,ss] = fracyear2date(uniqueFracYear(jj));
        
        dateTime{jj,1} = ['"',num2str(yy(ii)),'-',num2str(M,'%0.2d'),'-',num2str(D,'%0.2d'),...
            ' ',num2str(hh,'%0.2d'),':',num2str(mm,'%0.2d'),':',num2str(ss,'%0.2d'),'"'];
        
        % Pair data from loggers and snow analysis with current image
        [minFracYearImg,imgIdx] = min(abs(fracYear-uniqueFracYear(jj)));
        [minFracYearOne,loggerOneIdx] = min(abs(loggerOneFracYear-uniqueFracYear(jj)));
        [minFracYearTwo,loggerTwoIdx] = min(abs(loggerTwoFracYear-uniqueFracYear(jj)));
        [sortedSnowFracYear,snowFracYearIdx] = sort(abs(snowDataFracYear-uniqueFracYear(jj)),'ascend');
        [minFracYearUSNR1,loggerUSNR1Idx] = min(abs(loggerUSNR1FracYear-uniqueFracYear(jj)));
        filenameStorage{jj,1} = NaN;
        jdStorage(jj,1) = NaN;
        img_DT(jj,1) = {NaN};
        LOne_DT(jj,1) = {NaN};
        airT(jj,1) = NaN;
        RH(jj,1) = NaN;
        plateT(jj,1) = NaN;
        skyT(jj,1) = NaN;
        canopyT1(jj,1) = NaN;
        LTwo_DT(jj,1) = {NaN};
        canopyT2(jj,1) = NaN;
        canopyT3(jj,1) = NaN;
        snowDT(jj,1) = {NaN};
        snowFlag(jj,1) = 0;
        snowValue(jj,1) = NaN;
        treesValue(jj,1) = NaN;
        if (minFracYearImg == 0)
            img_DT(jj,1) = imgDateTime(imgIdx);
            filenameStorage{jj,1} = processFiles(imgIdx).name;
            jdStorage(jj,1) = jd(imgIdx);
        end
        if (minFracYearOne <= 3.794778384942320e-06)  % 2 minute window around image time
            if (minFracYearImg == 0) || (minFracYearImg >= 8.561643835616438e-06) %only process data if image time or no image within 4.5 minutes
                LOne_DT(jj,1) = loggerOne_DT(loggerOneIdx);
                jdStorage(jj,1) = loggerOneJD(loggerOneIdx);
                airT(jj,1) = loggerOne_airT(loggerOneIdx);
                RH(jj,1) = loggerOne_RH(loggerOneIdx);
                plateT(jj,1) = loggerOne_plateT(loggerOneIdx);
                %skyT(jj,1) = loggerOne_skyT(loggerOneIdx);
                canopyT1(jj,1) = loggerOne_canopyT1(loggerOneIdx);
                snowState = snowData_snowy(snowFracYearIdx(sortedSnowFracYear <= 0.002732240437158)); %only use snow data from within 24 hrs of image time
                solarElevation = snowData_solarElevation(snowFracYearIdx(sortedSnowFracYear <= 0.002732240437158));
                snowDateTime = snowData_DT(snowFracYearIdx(sortedSnowFracYear <= 0.002732240437158));
                valueSnow = snowData_snowValue(snowFracYearIdx(sortedSnowFracYear <= 0.002732240437158));
                valueTrees = snowData_treesValue(snowFracYearIdx(sortedSnowFracYear <= 0.002732240437158));
                if ~(isempty(snowState)) %limit data to when sun is above user-defined minimum elevation
                    snowStateRefined = snowState(solarElevation >= minSolarElevation);
                    snowDateTimeRefined = snowDateTime(solarElevation >= minSolarElevation);
                    valueSnowRefined = valueSnow(solarElevation >= minSolarElevation);
                    valueTreesRefined = valueTrees(solarElevation >= minSolarElevation);
                else
                    snowStateRefined = [];
                    snowDateTimeRefined = [];
                    valueSnowRefined = [];
                    valueTreesRefined = [];
                end
                if (~isempty(snowStateRefined) && (numel(snowStateRefined) < nearestNeighborImages) && (airT(jj,1) < minAirTemp))
                    snowDT(jj,1) = {strjoin(snowDateTimeRefined,';')};
                    snowFlag(jj,1) = mode(snowStateRefined);
                    snowValue(jj,1) = mean(valueSnowRefined);
                    treesValue(jj,1) = mean(valueTreesRefined);
                elseif (~isempty(snowStateRefined) && (numel(snowStateRefined) >= nearestNeighborImages) && (airT(jj,1) < minAirTemp))
                    snowDateTimeSubset = snowDateTimeRefined(1:nearestNeighborImages);
                    snowDT(jj,1) = {strjoin(snowDateTimeSubset,';')};
                    snowFlag(jj,1) = mode(snowStateRefined(1:nearestNeighborImages));
                    snowValue(jj,1) = mean(valueSnowRefined(1:nearestNeighborImages));
                    treesValue(jj,1) = mean(valueTreesRefined(1:nearestNeighborImages));
                end
                clear snowState solarElevation snowDateTime valueSnow valueTrees...
                    snowStateRefined snowDateTimeRefined valueSnowRefined...
                    valueTreesRefined snowDateTimeSubset
            end
        end
        if (minFracYearTwo <= 3.794778384942320e-06) % 2 minute window around image time
            if (minFracYearImg == 0) || (minFracYearImg >= 8.561643835616438e-06) %only process data if image time or no image within 4.5 minutes
                LTwo_DT(jj,1) = loggerTwo_DT(loggerTwoIdx);
                jdStorage(jj,1) = loggerTwoJD(loggerTwoIdx);
                canopyT2(jj,1) = loggerTwo_canopyT2(loggerTwoIdx);
                canopyT3(jj,1) = loggerTwo_canopyT3(loggerTwoIdx);
            end
        end
        if (minFracYearUSNR1 <= 5.692167577413479e-05) % 30 minute window around image time
            if (minFracYearImg == 0) || (minFracYearImg >= 8.561643835616438e-06) %only process data if image time or no image within 4.5 minutes
                skyT(jj,1) = loggerUSNR1_skyTemp(loggerUSNR1Idx);
            end
        end
        
        % Open image
        if (minFracYearImg == 0)
            filename = processFiles(imgIdx).name;
            fInfo = irFileOpen(imgDir{ii},filename,'seq',false);
            
            % Process ROIs for uncorrected image
            img = integer2temp(fInfo,false);
            ss = find((shiftFracYear <= uniqueFracYear(jj)),1,'last'); % index for shifted masks ROIs
            if ~isempty(img)
                for kk = 1:numel(roiFile)
                    roiIdx = find(mask{kk,ss});
                    roiPixels = img(roiIdx);
                    roiNumPixels(jj,kk) = numel(roiPixels);
                    roiUncorrectedMean(jj,kk) = mean(roiPixels(:));
                    % roiStd(jj,kk) = std(roiPixels(:));
                    % roiMin(jj,kk) = min(roiPixels(:));
                    % roiMax(jj,kk) = max(roiPixels(:));
                end
            end
            
            % Process ROIs for corrected image
            if applyCorrections && numel(LOne_DT{jj,1} > 1)
                corrFactors = {'distance map',pixelDist,'object emissivity',objectEmissivity,...
                    'air temperature',airT(jj,1),'RH',RH(jj,1),'reflected temperature',skyT(jj,1),...
                    'correction temp units','C'};
                corrFactorsPlate = {'distance map',2,'object emissivity',0.985,...
                    'air temperature',airT(jj,1),'RH',RH(jj,1),'reflected temperature',skyT(jj,1),...
                    'correction temp units','C'};
                img = integer2temp(fInfo,true,corrFactors);
                imgPlate = integer2temp(fInfo,true,corrFactorsPlate);
                if ~isempty(img)
                    for kk = 1:numel(roiFile)
                        if kk == numel(roiFile)
                            roiIdx = find(mask{kk,ss});
                            roiPixels = imgPlate(roiIdx);
                            roiNumPixels(jj,kk) = numel(roiPixels);
                            roiCorrectedMean(jj,kk) = mean(roiPixels(:));
                        else
                            roiIdx = find(mask{kk,ss});
                            roiPixels = img(roiIdx);
                            roiNumPixels(jj,kk) = numel(roiPixels);
                            roiCorrectedMean(jj,kk) = mean(roiPixels(:));
                        end
                        % roiStd(jj,kk) = std(roiPixels(:));
                        % roiMin(jj,kk) = min(roiPixels(:));
                        % roiMax(jj,kk) = max(roiPixels(:));
                    end
                end
            else
                for kk = 1:numel(roiFile)
                    roiCorrectedMean(jj,kk) = NaN;
                    % roiStd(jj,kk) = NaN;
                    % roiMin(jj,kk) = NaN;
                    % roiMax(jj,kk) = NaN;
                end
            end
        else
            for kk = 1:numel(roiFile)
                roiNumPixels(jj,kk) = NaN;
                roiUncorrectedMean(jj,kk) = NaN;
                roiCorrectedMean(jj,kk) = NaN;
                % roiStd(jj,kk) = NaN;
                % roiMin(jj,kk) = NaN;
                % roiMax(jj,kk) = NaN;
            end
        end
        
        clear fInfo img corrFactors roiIdx roiPixels
        
        % Update waitbar every 25 datapoints
        if (mod(jj,100) == 0)
            waitbar(jj/numel(uniqueFracYear),w,['Working on directory ',num2str(ii),' of ',num2str(numel(imgDir)),', datapoint ',num2str(jj),'/',num2str(numel(uniqueFracYear))]);
        end
        
        
    end
    
    % Prepare output structure
    dataStorage{ii}.filenames = filenameStorage;
    dataStorage{ii}.date = [img_DT,LOne_DT,LTwo_DT,snowDT];
    dataStorage{ii}.DOY = jdStorage;
    dataStorage{ii}.fracYear = uniqueFracYear;
    dataStorage{ii}.RH = RH;
    dataStorage{ii}.airT = airT;
    dataStorage{ii}.plateT = plateT;
    dataStorage{ii}.skyT = skyT;
    dataStorage{ii}.snowFlag = snowFlag;
    dataStorage{ii}.snowValue = snowValue;
    dataStorage{ii}.treesValue = treesValue;
    dataStorage{ii}.canopyT1 = canopyT1;
    dataStorage{ii}.canopyT2 = canopyT2;
    dataStorage{ii}.canopyT3 = canopyT3;
    dataStorage{ii}.roiUncorrectedMean = roiUncorrectedMean-273.15;
    if applyCorrections
        dataStorage{ii}.roiCorrectedMean = roiCorrectedMean-273.15;
    end
    
     % Remove points with missing correction, if applying corrections
    if applyCorrections
        deleteIdx = find(isnan(RH));
        dataStorage{ii}.filenames(deleteIdx,:) = [];
        dataStorage{ii}.date(deleteIdx,:) = [];
        dataStorage{ii}.DOY(deleteIdx,:) = [];
        dataStorage{ii}.fracYear(deleteIdx,:) = [];
        dataStorage{ii}.RH(deleteIdx,:) = [];
        dataStorage{ii}.airT(deleteIdx,:) = [];
        dataStorage{ii}.plateT(deleteIdx,:) = [];
        dataStorage{ii}.skyT(deleteIdx,:) = [];
        dataStorage{ii}.snowFlag(deleteIdx,:) = [];
        dataStorage{ii}.snowValue(deleteIdx,:) = [];
        dataStorage{ii}.treesValue(deleteIdx,:) = [];
        dataStorage{ii}.canopyT1(deleteIdx,:) = [];
        dataStorage{ii}.canopyT2(deleteIdx,:) = [];
        dataStorage{ii}.canopyT3(deleteIdx,:) = [];
        dataStorage{ii}.roiUncorrectedMean(deleteIdx,:) = [];
        dataStorage{ii}.roiCorrectedMean(deleteIdx,:) = [];
    end
    
    close(w);
    
    clear jd filenameStorage imgDateTime dateTime loggerOneCurrentYear...
        loggerTwoCurrentYear uniqueFracYear img_DT LOne_DT LTwo_DT jdStorage...
        fracYear RH airT plateT skyT roiNumPixels roiUncorrectedMean roiCorrectedMean...
        snowDT snowFlag snowValue treesValue
    
end

cd(startDir);


% create header and write file
if applyCorrections
    headerMet = {'Filename','DOY','fracYear','RH','airT','plateT','skyT','snowFlag','canopyT1','canopyT2','canopyT3','F2_uncorrected','F3_uncorrected','P2_uncorrected','P3_uncorrected','S1_uncorrected','S3_uncorrected','plateT_camera_uncorrected','F2','F3','P2','P3','S1','S3','plateT_camera'};
    filename = strcat('niwotData_',datestr(now,'yyyymmdd-HHMMSS'),'.txt');
    fid = fopen(filename,'w');
    fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',headerMet{:});
    for ii = 1:numel(yy)
        for jj = 1:size(dataStorage{ii}.filenames,1)
            fprintf(fid,'%s,%f,%f,%f,%f,%f,%f,%i,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n',dataStorage{ii}.filenames{jj,1},...
                dataStorage{ii}.DOY(jj,1),dataStorage{ii}.fracYear(jj,1),dataStorage{ii}.RH(jj,1),...
                dataStorage{ii}.airT(jj,1),dataStorage{ii}.plateT(jj,1),dataStorage{ii}.skyT(jj,1),dataStorage{ii}.snowFlag(jj,1),...
                dataStorage{ii}.canopyT1(jj,1),dataStorage{ii}.canopyT2(jj,1),dataStorage{ii}.canopyT3(jj,1),...
                dataStorage{ii}.roiUncorrectedMean(jj,:),dataStorage{ii}.roiCorrectedMean(jj,:));
        end
    end
    fclose(fid);
else
    headerMet = {'Filename','DOY','fracYear','snowFlag','canopyT1','canopyT2','canopyT3','F2_uncorrected','F3_uncorrected','P2_uncorrected','P3_uncorrected','S1_uncorrected','S3_uncorrected','plateT_camera_uncorrected'};
    filename = strcat('niwotData_',datestr(now,'yyyymmdd-HHMMSS'),'.txt');
    fid = fopen(filename,'w');
    fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',headerMet{:});
    for ii = 1:numel(yy)
        for jj = 1:size(dataStorage{ii}.filenames,1)
            fprintf(fid,'%s,%f,%f,%i,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n',dataStorage{ii}.filenames{jj,1},...
                dataStorage{ii}.DOY(jj,1),dataStorage{ii}.fracYear(jj,1),dataStorage{ii}.snowFlag(jj,1),...
                dataStorage{ii}.canopyT1(jj,1),dataStorage{ii}.canopyT2(jj,1),dataStorage{ii}.canopyT3(jj,1),...
                dataStorage{ii}.roiUncorrectedMean(jj,:));
        end
    end
    fclose(fid);
end