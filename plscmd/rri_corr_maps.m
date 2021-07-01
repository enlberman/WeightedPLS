function maps=rri_corr_maps(behav,datamat,n,k,cormode, doCCA, varargin);

% creates image-wide correlation map for k scans with behavior vector , or
% whatever is in the behav variable (eg. seed voxel)
%
% written 4.12.98 ARM
% syntax maps=rri_corr_maps(behav,datamat,n,k);

if ~exist('cormode','var'), cormode = 0; end
if ~isempty(varargin)
    weights_dat = varargin{2};
    weights_beh = varargin{4};
end
maps=[];

	for i=1:k
		temp=[];
        if exist('weights_dat','var')
            temp=rri_xcor((behav(1+(n* (i-1) ):n*i,:)),(datamat(1+(n*(i-1)):n*i,:)),cormode, doCCA, 'weights_dat',weights_dat(1+(n*(i-1)):n*i),'weights_beh',weights_beh(1+(n*(i-1)):n*i));
        else
            temp=rri_xcor((behav(1+(n* (i-1) ):n*i,:)),(datamat(1+(n*(i-1)):n*i,:)),cormode, doCCA);
        end

		maps=[maps;temp];

	end

%disp(' ')
%disp('Program complete')
