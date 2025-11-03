% will run IVIM with bayesian algorithm assuming the folder was sorted by "SortIVIM_Pelvis.m" - Mira Liu Nov 2025
% input folder of interest and name 
function  RunBayesianIVIM_Pelvis(varargin)
    DIR = fullfile(varargin{1},varargin{2});
    namesave = strcat(varargin{2},'_Bayesian.mat'); %number.mat
    SaveDIR = fullfile(varargin{1},namesave);
    
    namesave = strcat(varargin{2},'_SegmentedFit.mat'); %number.mat
    SaveDIR_2 = fullfile(varargin{1},namesave);

    SaveDIR_3 = fullfile(varargin{1},strcat(varargin{2},'_b0map.mat'));
    
    LoadedMask = load(fullfile(varargin{1},'b0_Dicoms','AllMasks.mat'));
    ROI_Mask = LoadedMask.All_Masks;
    if size(dir(fullfile(DIR, '/IM*.dcm')),1)> 0 
        fprintf(strcat('Running: ',DIR,'\n'))
        tic
        disp(datetime)
        [f, Dstar, D, Resid, f_2step, Dstar_2step, D_2step, Resid_2step,b0_map] = IVIM_Reconstruction(DIR,ROI_Mask); %drop need for type, but keep it anyway
        save (SaveDIR,'f','Dstar','D','Resid');
        save (SaveDIR_2,'f_2step','Dstar_2step','D_2step','Resid_2step');
        save (SaveDIR_3,'b0_map');
        fprintf(strcat('Completed: ',DIR,'\n'))

       
        toc
        disp(datetime)
    else
        size(dir(fullfile(DIR, '/*.dcm')),1)
        error('incorrect number of images')
    end
    
    end
    
    
    
    function [IVIM_f, IVIM_Dstar, IVIM_D, IVIM_ADJRsq,   IVIM_f_2step, IVIM_Dstar_2step, IVIM_D_2step, IVIM_ADJRsq_2step, b0_map] = IVIM_Reconstruction(Image_Directory,ROI_Mask)
    
    dat_list = dir(fullfile(Image_Directory,'IM*dcm'));
    datnames = {dat_list.name}; %read them in the correct order
    datnames = natsortfiles(datnames);
    
    
    fname  = fullfile(Image_Directory,dat_list(1).name);
    header = dicominfo(fname);
    nx     = header.Height;
    ny     = header.Width;
     
    [Total_Images, ~] = size(dat_list);
    
    SliceLocation = zeros(1,length(dat_list));
    BVals = zeros(1,length(dat_list));
    
    for j = 1:length(dat_list)
        h=dicominfo(char(fullfile(Image_Directory,dat_list(j).name)));
        SliceLocation(j) = h.SliceLocation;
        bval_name = h.SequenceName;
        bval = extractBetween(bval_name,'*ep_b', 't');
        if isempty(bval)
            bval = extractAfter(bval_name,'*ep_b');
        end
        BVals(j) = string(bval);
    end
    
    
    Bvalues = unique(BVals);
    Total_Slices = length(unique(SliceLocation));
    Num_Bvalues = length(Bvalues);
    Start_Index      = 1;
    Images_Per_Slice = Num_Bvalues;
    Z_Filter         = 0;
    
    
    Images_T         = zeros(Num_Bvalues,nx,ny);
    
    IVIM_f     = zeros(nx,ny,Total_Slices);
    IVIM_Dstar = zeros(nx,ny,Total_Slices);
    IVIM_D     = zeros(nx,ny,Total_Slices);
    IVIM_ADJRsq  = zeros(nx,ny,Total_Slices);
    
    
    IVIM_f_2step     = zeros(nx,ny,Total_Slices);
    IVIM_Dstar_2step = zeros(nx,ny,Total_Slices);
    IVIM_D_2step     = zeros(nx,ny,Total_Slices);
    IVIM_ADJRsq_2step  = zeros(nx,ny,Total_Slices);

    b0_map = zeros(nx,ny,Total_Slices);
    
    if (Z_Filter == 0 )
        Start_Slice = 1;
        End_Slice   = Total_Slices;
    else
        Start_Slice = 3;
        End_Slice   = Total_Slices-3;
    end
    

    %% only process the slices with an ROI
    slices = squeeze(sum(ROI_Mask,[1,2]));
    
    for Islice=min(find(slices)):max(find(slices))%20%Start_Slice:End_Slice %3 - 47
        fprintf(' Working on Slice Number %d \n',Islice);
        
        i1 = Images_Per_Slice*(Islice-1)+Start_Index; %37*(1-1)+28 = 28
        i2 = i1 + Num_Bvalues-1;
        
        jj = 1; %b value
        if (Z_Filter == 0)
            for i= i1:i2
                fname_1 = fullfile(Image_Directory,char(datnames(i)));
                header1 = dicominfo(fname_1);
                Images_T(jj,:,:)= double(dicomread(fname_1)); 
                jj= jj+1;
            end
        else
            w2 = 0.0; %set to zero to see if red dots are continuous
            w1 = 1; %set to zero to see
            w0 = 0.0; 
            for i= i1:i2 % Z-Smoothing across slices? 
                fname_im2 = fullfile(Image_Directory,char(datnames(i-2*(Images_Per_Slice))));
                fname_im1 = fullfile(Image_Directory,char(datnames(i-1*(Images_Per_Slice))));
                fname_im0 = fullfile(Image_Directory,char(datnames(i-0*(Images_Per_Slice)))); % center slice
                fname_ip1 = fullfile(Image_Directory,char(datnames(i+1*(Images_Per_Slice))));
                fname_ip2 = fullfile(Image_Directory,char(datnames(i+2*(Images_Per_Slice))));
                %{
                header_im2 = dicominfo(fname_im2);
                header_im1 = dicominfo(fname_im1);
                header_im0 = dicominfo(fname_im0);
                header_ip1 = dicominfo(fname_ip1);
                header_ip2 = dicominfo(fname_ip2);
                %same slice location, same b value. From slice 3 - 47, all b
                %values, diffusion gradient AVERAGE.
                
                fprintf(' Slice Locs : m2 %3.2f m1 %3.2f m0 %3.2f p1 %3.2f p2 %3.2f ', ...
                          header_im2.SliceLocation                             , ...
                          header_im1.SliceLocation                             , ...
                          header_im0.SliceLocation                             , ...
                          header_ip1.SliceLocation                             , ... 
                          header_ip2.SliceLocation                             );
                
                fprintf(' B values  : m2 %i m1 %i m0 %i P1 %i p2 %i \n', ...
                          header_im2.DiffusionBValue                          , ...
                          header_im1.DiffusionBValue                          , ...
                          header_im0.DiffusionBValue                          , ...
                          header_ip1.DiffusionBValue                          , ... 
                          header_ip2.DiffusionBValue                           );
                
                 Python_Image(Islice,jj,:,:) = w2*double(dicomread(fname_im2)) +           ...   
                                    w1*double(dicomread(fname_im1)) +           ... 
                                    w0*double(dicomread(fname_im0)) +           ... 
                                    w1*double(dicomread(fname_ip1)) +           ... 
                                    w2*double(dicomread(fname_ip2));
                %}
                 Images_T(jj,:,:) = w2*double(dicomread(fname_im2)) +           ...   
                                    w1*double(dicomread(fname_im1)) +           ... 
                                    w0*double(dicomread(fname_im0)) +           ... 
                                    w1*double(dicomread(fname_ip1)) +           ... 
                                    w2*double(dicomread(fname_ip2));
                jj= jj+1;
            end
        end
    
        %imagestack(Images_T,'permute')
        %pause()
        %close()

        %% save b0
        b0_map(:,:,Islice) = squeeze(Images_T(1,:,:));
    
        %% APPLY THE ROIs
        %size(ROI_Mask)
        test = ROI_Mask(:,:,Islice);
        check = repmat(test,[1,1,10]);
        ROI_slice = permute(check,[3,1,2]);
        Images_T = Images_T.*ROI_slice;
    
        
        %% now fit
        f_map = zeros(nx,ny,1);
        D_map=zeros(nx,ny,1);
        Dstar_map=zeros(nx,ny,1);
        Adjrsq_map=zeros(nx,ny,1);
    
        f_map_2step = zeros(nx,ny,1);
        D_map_2step=zeros(nx,ny,1);
        Dstar_map_2step=zeros(nx,ny,1);
        Adjrsq_map_2step=zeros(nx,ny,1);
    
        
        disp(' Fitting...');      
        for i=1:nx
            for j=1:ny
                if (Images_T(1,i,j) >= 1 ) 
                % for normal b values
                    SignalInput = squeeze(double(Images_T(1:Num_Bvalues,i,j)/Images_T(1,i,j))); 
        
                    % Bayesian %% you can also comment out the Bayesian fit with %{ in line 198 if you want ONLY the segmented fit.
                    
                    Output = BayesianIVIM_Pelvis(Bvalues,SignalInput);
                    f_map(i,j,1)=Output.f;
                    D_map(i,j,1)=Output.D;
                    Dstar_map(i,j,1)=Output.Ds;
                    Adjrsq_map(i,j,1)=Output.adj_rsq;
                    %}
                    
                    % Two-step 
                    %% in line 207 below, put %{ to comment out the following block of code meaning it won't be run. Instead now they will all remain zeros.
                    
                    Output = Segmented_Fit_Pelvis(Bvalues, SignalInput);
                    f_map_2step(i,j,1)=Output.f; %% 
                    D_map_2step(i,j,1)=Output.D;
                    Dstar_map_2step(i,j,1)=Output.Ds;
                    Adjrsq_map_2step(i,j,1)=Output.adj_rsq;
                    %} 
                end
            end
        end
        %}
    
         
        IVIM_f(:,:,Islice)     = f_map;
        IVIM_Dstar(:,:,Islice) = Dstar_map;
        IVIM_D(:,:,Islice)     = D_map;
        IVIM_ADJRsq(:,:,Islice)  = Adjrsq_map;
        
        
        
        IVIM_f_2step(:,:,Islice)     = f_map_2step;
        IVIM_Dstar_2step(:,:,Islice) = Dstar_map_2step;
        IVIM_D_2step(:,:,Islice)     = D_map_2step;
        IVIM_ADJRsq_2step(:,:,Islice)  = Adjrsq_map_2step;

        
    
        
    end

end



function Output = BayesianIVIM_Pelvis(b,data)

% https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.28671
% https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.25765
% https://onlinelibrary.wiley.com/doi/abs/10.1002/cmr.a.20043
% make sure they are column arrays
b = double(b(:));
data = double(data(:));
%class(data)
sumData2 = sum(data.^2);

% array of weights in case the NSA at different b-values are different
% weight(i) = NSA for data(i)
weight = ones(size(data));
organ = 'pelvis';


weight = weight(:);

% ADC monoexp fit
%ln_S=log(data);
%p = polyfit(b, ln_S, 1);
%ADC = -(p(1));

% range for the parameters
logDgrid = 120; logD = linspace(-10,-4,logDgrid);
logDsgrid = 120; logDs = linspace(-12,0,logDsgrid); % adjusted Ds to be larger range?
fGrid = 60; f = linspace(0,1,fGrid);

% log-priors
if strcmp(organ,'kidney')
    logDmean = -6.2;
    logDstd = 1;
    logPrD = -0.5*(logD - logDmean).^2/logDstd^2 - 0.5*log(2*pi*logDstd^2);
    %
    logDsmean = -3.5;%init: -3.5 
    logDsStd = 1;
    logPrDs = -0.5*(logDs - logDsmean).^2/logDsStd^2 - 0.5*log(2*pi*logDsStd^2);
elseif strcmp(organ,'liver')
    logDmean = -7.0;
    logDstd = 1;
    logPrD = -0.5*(logD - logDmean).^2/logDstd^2 - 0.5*log(2*pi*logDstd^2);
    %
    logDsmean = -3.5;%init: -3.5 
    logDsStd = 1;
    logPrDs = -0.5*(logDs - logDsmean).^2/logDsStd^2 - 0.5*log(2*pi*logDsStd^2);
elseif strcmp(organ,'brain')
    logDmean = -6.9; %initial Dt = 0.001
    logDstd = 1;
    logPrD = -0.5*(logD - logDmean).^2/logDstd^2 - 0.5*log(2*pi*logDstd^2);
    %
    logDsmean = -4.6;%initial D* = 0.01
    logDsStd = 1;
    logPrDs = -0.5*(logDs - logDsmean).^2/logDsStd^2 - 0.5*log(2*pi*logDsStd^2);

elseif strcmp(organ,'pelvis')
    logDmean = -6; %initial Dt = 0.001 ??????? 
    logDstd = 1;
    logPrD = -0.5*(logD - logDmean).^2/logDstd^2 - 0.5*log(2*pi*logDstd^2);
    %
    logDsmean = -5;%initial D* = 0.05 ??????? 
    logDsStd = 1;
    logPrDs = -0.5*(logDs - logDsmean).^2/logDsStd^2 - 0.5*log(2*pi*logDsStd^2);
end
Pr = (exp(logPrD')*exp(logPrDs)).*(repmat(logD',1,logDsgrid)<repmat(logDs,logDgrid,1));
logPr = log(Pr);

% various modelling arrays
exp_b_logD = repmat(exp(-b*exp(logD)),[1 1 logDsgrid]);
exp_b_logDs = permute(repmat(exp(-b*exp(logDs)),[1 1 logDgrid]),[1 3 2]);

% compute log-likelihood and log posterior by looping over f (loops in D
% and D* are incorporated here using array operations)
[logLh,logPost] = deal(zeros(fGrid,logDgrid,logDsgrid));
for nf = 1:fGrid
    G = f(nf)*exp_b_logDs + (1-f(nf))*exp_b_logD;
    logLh(nf,:,:) = squeeze(-0.5*length(b)*log(sumData2 - sum(repmat(data.*weight,[1 logDgrid logDsgrid]).*G).^2./sum(G.^2.*repmat(weight,[1 logDgrid logDsgrid]))));
    logPost(nf,:,:) = logLh(nf,:,:) + reshape(logPr,[1 logDgrid logDsgrid]);
end

% log-posterior
Post = exp(logPost - max(logPost(:)));
% normalisation
sP = sum(Post(:));


% compute MMSE estimates
Output.D = exp(sum(sum(sum(Post,1),3).*logD)/sP);
Output.Ds = exp(sum(squeeze(sum(sum(Post,1),2)).*logDs')/sP);
Output.f = sum(sum(sum(Post,2),3).*f')/sP;


% compute curve for MMSE estimate
G = Output.f*exp(-b*Output.Ds) + (1-Output.f)*exp(-b*Output.D);
curveFit = G*inv(G'*G)*G'*data;
curveFit1 = curveFit/max(curveFit(:));
data=data/max(data(:));

% compute root square error
yresid = data - curveFit1;
SSresid = sum(yresid.^2);
SStotal = (length(data)-1) * var(data);
rsq = 1 - SSresid/SStotal;

k=3; %number of independent variables
n=length(b); %number of datapoints
adj_rsq = 1- ((1 - rsq)*(n -1)/ (n - k - 1));

Output.SSE = SSresid;
Output.rsq = rsq;
Output.adj_rsq=adj_rsq;


end

%% algorithm 6 is similar to algorithm 2, and carries over f and D. only fits D* to the entire set of b-values
%algorithm 6, segmented fit with f and D held, D* fit in the second run
%but also forced D* upper bound (unlike algo5)
function Output = Segmented_Fit_Pelvis(bvalues, signal,bval_cutoff_idx)


    if ~exist('bval_cutoff_idx','var')
         % third parameter does not exist, so default it to something
          %bval_cutoff_idx = 6;%the threshold of bvalues where b>bvalue(blim) is slow only.
          bval_cutoff_idx = find(bvalues>200,1);
     end
    vec = log(signal/signal(1)); % ln(S/S0) normalize all signal
    vec=double(vec(:)); %force to be column vector
  
    bvalues = double(bvalues(:));
    N_bvalues=length(bvalues);


    %generate D map
    %plot(bvalues(1:N_bvalues),signal/signal(1))
    %ylim([0,1])
    try
        [mono_fitresult, ~] = fit(bvalues(bval_cutoff_idx:N_bvalues),vec(bval_cutoff_idx:N_bvalues),'poly1');
        
        D_fit = -mono_fitresult.p1;
        f_fit = 1-exp(mono_fitresult.p2); % as  ln(Ae^-bx) = ln(a) - bx.
    
    
    
        %generate D* and f map
        ft_bi = fittype('(1-f)*exp(-x*D)+f*exp(-x*(Dstar))','dependent',{'y'},'independent',{'x'},'problem',{'D', 'f'},'coefficients',{'Dstar'});
        fo_bi = 0.005; %startpoints, in alphabetical order, so D, Dstar, f. 
    
        try
            [fitmod_bi,good_bi,~]=fit(bvalues,signal(:)/signal(1),ft_bi,'startpoint',fo_bi, 'upper',0.1, 'problem', {D_fit, f_fit}); %no longer in log space
            Output.D=D_fit;
            Output.Ds=fitmod_bi.Dstar;
            Output.f = f_fit;
            Output.SSE = good_bi.sse;
            Output.rsq = good_bi.rsquare;
            Output.adj_rsq = good_bi.adjrsquare;
        catch
            Output.D=0;
            Output.Ds=0;
            Output.f = 0;
            Output.SSE = NaN;
            Output.rsq = NaN;
            Output.adj_rsq = NaN;
        end
    catch
        Output.D=0;
        Output.Ds=0;
        Output.f = 0;
        Output.SSE = NaN;
        Output.rsq = NaN;
        Output.adj_rsq = NaN;
    end

end
