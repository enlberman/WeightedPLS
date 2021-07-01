% function rri_xcor
%  syntax:  outmat=rri_xcor(design,datamat)
%  Computes crosscorrelation of two matrices
%  Written by ARM 12-15-94


function[outmat]=rri_xcor(design,datamat,cormode,doCCA,varargin)

if ~exist('cormode','var'), cormode=0; end;

if ~isempty(varargin)
    weights_dat = varargin{2};
    weights_beh = varargin{4};
end

[r c]=size(datamat);
[dr dc]=size(design);
  if r ~= dr
    error('Error in rri_xcor: input matrices must have same number of rows');
  end

if doCCA==0
  switch cormode
   case 0

        if ~exist('weights_dat','var')
            avg=mean(datamat);
            stdev=std(datamat);
            checknan=find(stdev==0);

            if (isempty(checknan)==0)
                 datamat(:,checknan)=0;
                 avg(checknan)=0;
                 stdev(checknan)=1;
            end %if

             for i=1:r
              datamat(i,:)=(datamat(i,:)-avg)./stdev;
             end

            davg=mean(design);
            dstdev=std(design);

            checknan=find(dstdev==0);
            if (isempty(checknan)==0)
                 design(:,checknan)=0;
                 davg(checknan)=0;
                 dstdev(checknan)=1;
            end %if

            for i=1:dr
             design(i,:)=(design(i,:)-davg)./dstdev;
            end

            xprod=design'*datamat;
            outmat=xprod./(r-1);
        else
            wavg=mean(datamat.*weights_dat)./sum(weights_dat);
            wstdev=sum(((datamat-wavg).^2).*weights_dat)/sum(weights_dat);
            
            wdavg=mean(design.*weights_beh)./sum(weights_beh);
            wdstdev=sum(((design-wdavg).^2).*weights_beh)/sum(weights_beh);
            
            %for i=1:dr
            % design(i,:)=(design(i,:)-davg)./dstdev;
            %end
            %for i=1:r
            %  datamat(i,:)=(datamat(i,:)-avg)./stdev;
            %end
            cov = (((design-wdavg)./sqrt(weights_beh))'*((datamat-wavg).*sqrt(weights_dat)))/sum(sqrt(weights_dat)./sqrt(weights_beh));
            stdev_mat = wdstdev'*wstdev;            
            outmat = cov./stdev_mat;
        end



   case 2, % covariance

    avg = mean(datamat); % columnwise mean
    davg = mean(design);
    for i=1:r, datamat(i,:)=datamat(i,:)-avg; end    
    for i=1:dr, design(i,:)=design(i,:)-davg; end    
    xprod = design'*datamat;
    outmat = xprod./(r-1);

   case 4, % cosine angle

    stdev=std(datamat);
    checknan = find(stdev==0);     
    if (isempty(checknan)==0)
      datamat(:,checknan)=0;
      stdev(checknan)=1;
    end 
    
    % clean up columns with stdev=0 in design
    dstdev=std(design);    
    checknan = find(dstdev==0);
    if (isempty(checknan)==0)
      design(:,checknan)=0;
      dstdev(checknan)=1;
    end 
    
    for i=1:r, datamat(i,:)= datamat(i,:)./stdev; end    
    for i=1:dr, design(i,:)= design(i,:)./dstdev; end    
    outmat = design'*datamat;

   case 6, % dot product
    outmat = design'*datamat;

  end
  
elseif doCCA==1
    switch cormode
   case 0

        if ~exist('weights_dat','var')
            avg=mean(datamat);
            stdev=std(datamat);
            checknan=find(stdev==0);

            if (isempty(checknan)==0)
                 datamat(:,checknan)=0;
                 avg(checknan)=0;
                 stdev(checknan)=1;
            end %if

             for i=1:r
              datamat(i,:)=(datamat(i,:)-avg)./stdev;
             end

            davg=mean(design);
            dstdev=std(design);

            checknan=find(dstdev==0);
            if (isempty(checknan)==0)
                 design(:,checknan)=0;
                 davg(checknan)=0;
                 dstdev(checknan)=1;
            end %if

            for i=1:dr
             design(i,:)=(design(i,:)-davg)./dstdev;
            end

            xprod=design'*datamat;
            xprod=xprod./(r-1);
            design_cov = (design'*design)./(r-1);
            datamat_cov = (datamat'*datamat)./(r-1);
            design_std = inv(sqrtm(design_cov));
            datamat_std = inv(sqrtm(datamat_cov));
            outmat=design_std*xprod*datamat_std;
        else
            wavg=mean(datamat.*weights_dat)./sum(weights_dat);
            wstdev=sum(((datamat-wavg).^2).*weights_dat)/sum(weights_dat);
            
            wdavg=mean(design.*weights_beh)./sum(weights_beh);
            wdstdev=sum(((design-wdavg).^2).*weights_beh)/sum(weights_beh);
            
            %for i=1:dr
            % design(i,:)=(design(i,:)-davg)./dstdev;
            %end
            %for i=1:r
            %  datamat(i,:)=(datamat(i,:)-avg)./stdev;
            %end
            cov = (((design-wdavg)./sqrt(weights_beh))'*((datamat-wavg).*sqrt(weights_dat)))/sum(sqrt(weights_dat)./sqrt(weights_beh));
            stdev_mat = wdstdev'*wstdev;            
            xprod = cov./stdev_mat;
            
            design_cov = (((design-wdavg)./sqrt(weights_beh))'*((design-wdavg).*sqrt(weights_beh)))/sum(sqrt(weights_beh)./sqrt(weights_beh));
            design_stdev_mat = wdstdev'*wdstdev;            
            design_cov = design_cov./design_stdev_mat;
            datamat_cov = (((datamat-wavg)./sqrt(weights_dat))'*((datamat-wavg).*sqrt(weights_dat)))/sum(sqrt(weights_dat)./sqrt(weights_dat));
            datamat_stdev_mat = wstdev'*wstdev;            
            datamat_cov = datamat_cov./datamat_stdev_mat;
            design_std = inv(sqrtm(design_cov));
            datamat_std = inv(sqrtm(datamat_cov));
            outmat=design_std*xprod*datamat_std;
        end
    end
end

