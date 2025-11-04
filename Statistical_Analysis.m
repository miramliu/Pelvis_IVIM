%% Perform basic statistical analyses on images given ROIs
% Input: Patient Identifier, name of MR paremeter (e.g. IVIM f or T1), the full path and name of the excel output
% and the loaded contour (from RT_struct as a mat file) and the parameter map as a matfile
% Output: An added line in the excel sheet of statistical analyses for this particular patient. 

% note that parameter map must be (nx, ny, nz): it can only be of one parameter. 
% Mira 11/4/2025




% assumes have loaded RT struct (with corresponding variable 'contours') and have loaded the parameter map of interest (Image_map)
function Statistical_Analysis(Patient_ID, Parameter_Name, ExcelFileName, contours, Image_Map)

    
    ROI_names = {};
    
    for n = 1:size(contours,2)
        ROI_name =  contours(n).ROIName; %get the nth ROI name
        mask = contours(n).Segmentation; % get the nth contour
        mask = flip(permute(mask,[2, 1, 3]),3);
    
    
        Parameter_Values = Image_Map.*mask;
    
    
        % histogram features
        volumetric_mean = mean(nonzeros(Parameter_Values));
        volumetric_stdev = std(nonzeros(Parameter_Values));
        volumetric_median = median(nonzeros(Parameter_Values));
        volumetric_skew = skewness(nonzeros(Parameter_Values),1,'all');
        volumetric_kurtosis = kurtosis(nonzeros(Parameter_Values),1,'all');
        volumetric_size = size(nonzeros(Parameter_Values),1);

        % accounting for outliers
        volumetric_iqr = iqr(nonzeros(Parameter_Values));
        iqr_25 = volumetric_median - volumetric_iqr;
        iqr_75 = volumetric_median + volumetric_iqr;
        lower_outliers = iqr_25 - 1.5*volumetric_iqr;
        upper_outliers = iqr_75 + 1.5*volumetric_iqr;

        trimmed_Parameter_Values = Parameter_Values(Parameter_Values > lower_outliers);
        trimmed_Parameter_Values = trimmed_Parameter_Values(trimmed_Parameter_Values < upper_outliers);

        trimmed_mean = mean(nonzeros(trimmed_Parameter_Values));
        trimmed_median = median(nonzeros(trimmed_Parameter_Values));

        outlier_size = size(nonzeros(Parameter_Values),1) - size(nonzeros(trimmed_Parameter_Values),1);






   
        Identifying_Info = {Patient_ID, [Patient_ID '_' ROI_name]};

        if exist(ExcelFileName,"file")
            sheets = sheetnames(ExcelFileName); % check if a sheet for this parameter exists already
            if sum(strcmp(sheets,Parameter_Name)) == 0 % if it doesn't exist
            
                disp('saving data in new excel sheet')
                dataarray= {volumetric_mean,volumetric_stdev,volumetric_median, volumetric_skew,volumetric_kurtosis,volumetric_size, volumetric_iqr, trimmed_mean, trimmed_median, outlier_size};

                header_cell = {'Patient ID',	'ROI name',	'mean',	'median',	'stdev',	'skew',	'kurtosis',	'number of data points', 'IQR',	'mean excluding outliers',	'median excluding outliers', 'number of outliers'};
                writecell(header_cell,ExcelFileName,'WriteMode','append','Sheet',Parameter_Name)
                Export_Cell = [Identifying_Info,dataarray];
                writecell(Export_Cell,ExcelFileName,'WriteMode','append','Sheet',Parameter_Name)
            else

                Existing_Data = readcell(ExcelFileName,'Range','A:B','Sheet',Parameter_Name); %read only identifying info that already exists
                MatchFunc = @(A,B)cellfun(@isequal,A,B);
                idx = cellfun(@(Existing_Data)all(MatchFunc(Identifying_Info,Existing_Data)),num2cell(Existing_Data,2));
            
                if sum(idx)==0
                    disp('saving data in excel')
                    dataarray= {volumetric_mean,volumetric_stdev,volumetric_median, volumetric_skew,volumetric_kurtosis,volumetric_size, volumetric_iqr, trimmed_mean, trimmed_median, outlier_size};
                    Export_Cell = [Identifying_Info,dataarray];
                    writecell(Export_Cell,ExcelFileName,'WriteMode','append','Sheet',Parameter_Name)
                else
                    disp('data already existed in excel')
                end
            end
        end

    end


end