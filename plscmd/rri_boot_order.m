function [boot_order, new_num_boot] ...
   = rri_boot_order(num_subj_lst, num_cond, num_boot, bscan, incl_seq, boot_type, from_nonstrat)
%
%  USAGE:  [boot_order, new_num_boot] = rri_boot_order(num_subj_lst, num_cond, num_boot [, bscan, incl_seq])
%   
%  Generate bootstrap sample order for all conditions of the groups.  
%  Assume the data rows are stored in the order of subject in condition
%  in group.
%
%  Natasha: add possibility of creating unstratified boot sampling (simply 
%  take any sample with replacment from the rows of the datatmat, no respect
%  for group/cond membership)
%   
%  NOTE: the bootstrap order are applied for the subjects within each
%        group only, i.e. sampling with replacement for the subjects
%        within group.  In order to make sure bootstrap results are
%        reliable, there must have at least 3 subjects in any of the
%        groups.
%
%  Method:  
%     The resampling orders are performed for each group separately.
%     The orders are generated by using MATLAB's 'rand' function. 
%
%  Input:
%         num_subj_lst: a N elements vector, each element specifies the
%			number of subjects of one of groups (subj_group).
%         num_cond:	number of conditions in the experiment
%			(assume the st_datamat has 1 row for 
%			each condition of each group).
%         num_boot:	number of resampling to be performed.
%         bscan:	In Multiblock PLS, you can specify a subset of
%			conditions that are used in multiblock PLS behav
%			block. e.g., bscan=[1 3] for 4 conditions. If it
%			is not specified, it means that all the conditions
%			are selected. Using the above example, bscan is
%			equivalent to [1 2 3 4] for 4 conditions. This 
%			option can be applied to method 4 and 6.
%         incl_seq:	switch to control whether the boot_order includes
%			a sequential order (original un-resampled order).
%			0 (default) means that sequential order should not
%			be included.
%			1 means that sequential order should be included.
%         boot_type:	'strat' (default for standard PLS approach), or
%			'nonstrat' for nonstratified boot samples.
%
%  Output: 
%         boot_order:	an MxN matrix that stores the new order for the N 
%			samples, where M is the total number of rows 
%			of the st_datamat.
%			(i.e M = sum(num_subj_lst) * num_cond)
%         new_num_boot:	revised number of bootstrap resampling
%

   if ~exist('bscan','var')
      bscan = 1:num_cond;
   end

   if ~exist('incl_seq','var')
      incl_seq = 0;
   end

   if ~exist('boot_type','var')
      boot_type = 'strat';
   end

   if isempty(boot_type) | ~ismember(boot_type,{'strat','nonstrat'})
      error('Field "boot_type" should be ''strat'' or ''nonstrat''');
   end

   total_subj = sum(num_subj_lst);
   num_group = length(num_subj_lst);
   num_cond0 = num_cond;
   num_cond = length(bscan);
   total_rows0 = num_cond0 * total_subj;
   total_rows = num_cond * total_subj;
   boot_order = [];
   new_num_boot = num_boot;

   if strcmp(boot_type,'nonstrat')

      from_nonstrat.total_rows0 = total_rows0;
      from_nonstrat.num_group = num_group;
      from_nonstrat.num_subj_lst = num_subj_lst;
      from_nonstrat.num_cond0 = num_cond0;
   
      %  first create order as if all subjects are in one group
      %
      [tmp_order, new_num_boot] = ...
		rri_boot_order(total_subj, num_cond, num_boot, ...
		bscan, 0, 'strat', from_nonstrat);

      %  this way we will respect that data from diff conds belongs to 
      %  the same subject. remap into original subj*cond*group ordering
      %
      boot_order = zeros(size(tmp_order));
      for r = 1:total_rows0
         for p=1:size(tmp_order,2)
            row = tmp_order(r,p);
            ss = mod(row,total_subj)+1; % subject index when all subs are in one group

            %  now find which group he belongs to
            %
            for g=1:num_group
               if ss > sum(num_subj_lst(1:(g-1))) & ss <= sum(num_subj_lst(1:g)), break; end
            end

            s = ss - sum(num_subj_lst(1:(g-1))); % this subject's index within his group (g)

            %  find condition
            %
            cond = ceil(row/total_subj);

            %  now map into row of datatmat
            %
            row1 = sum(num_subj_lst(1:(g-1)))*num_cond0 + num_subj_lst(g)*(cond-1) + s;
            boot_order(r,p) = row1;
         end
      end

      return;
   end

   %  if num_subj <=8 for all groups then create all theoretically possible
   %  boot samples
   %
   [min_subj_per_group, is_boot_samples, boot_samples, new_num_boot] = ...
	rri_boot_check(num_subj_lst, num_cond, num_boot, incl_seq);

   %  determine tmp_boot_order, which is re-ordered subject index matrix
   %
   tmp_boot_order = zeros(total_subj, new_num_boot);

   for p=1:new_num_boot,

      subj_order = cell(1,num_group);
      not_done = 1;
      cnt = 0;

      while (not_done)

         start_subj = 1;

         for g = 1:num_group
            num_subj = num_subj_lst(g);

            %  reorder tasks for the current group.
            all_samples_are_same = 1;

            while (all_samples_are_same)

               if is_boot_samples(g)		% get from boot_samples

                  new_subj_order = boot_samples{g}(p,:);
                  all_samples_are_same = 0;

                  not_done = 0;

               else

                  not_done = 1;

                  new_subj_order = floor(rand(1,num_subj)*num_subj) + 1;

                  % first_subj = new_subj_order(1);
                  % if ~isequal(new_subj_order,repmat(first_subj,1,num_subj))

                  test=length(unique(new_subj_order));

                  %% check to make sure there are more than n/2 people
                  %%
                  %if (test >= num_subj/2)

                  % check to make sure there are at lease min_subj_per_group people
                  %
                  if (test >= min_subj_per_group)
                     all_samples_are_same = 0;
                  end;
               end
            end;

            subj_order{g} = new_subj_order + start_subj - 1; 
            start_subj = start_subj + num_subj;

         end;

         if ~all(is_boot_samples)

            %  make sure the order is not a repeated one
            not_done = 0;
            for i=1:p-1,
               %if isequal(squeeze(tmp_boot_order(:,i)),subj_order)
               if isequal(tmp_boot_order(:,i),[subj_order{:}]')
                  not_done = 1;
                  break;
               end;
            end;

           %  treat sequential order as duplicated one
           %
           if ~incl_seq & isequal([1:total_rows], [subj_order{:}])
              not_done = 1;
           end

            %  discard if all elements are the same
            cnt = cnt+1;
            if (cnt > 500),
               not_done = 0;
               disp('ERROR:  Duplicated bootstrap orders are used!');
            end;

         end;		% if ~all(is_boot_samples)

      end;		% while (not_done)

      tmp_boot_order(:,p) = [subj_order{:}]';

   end;  % for new_num_boot

  %  construct the resampling order matrix for bootstrap
  %
%  row_idx = reshape([1:total_rows],num_cond,total_subj);

  first = 1;
  last = 0;
  row_idx = [];
  for g=1:length(num_subj_lst)
     last = last + num_cond*num_subj_lst(g);
     tmp = reshape([first:last],num_subj_lst(g),num_cond);
     row_idx = [row_idx, tmp'];
     first = last+1;
  end

  boot_order = zeros(num_cond,total_subj, new_num_boot);

  for p=1:new_num_boot
     boot_order(:,:,p) = row_idx(:, tmp_boot_order(:,p));
  end

  b_order = [];
  for g=1:length(num_subj_lst)

     one_group = [];
     for p=1:new_num_boot
        tmp = ...
           boot_order(:,[sum(num_subj_lst(1:(g-1)))+1:sum(num_subj_lst(1:g))],p);
        tmp = reshape(tmp', [num_cond*num_subj_lst(g),1]);
        one_group = [one_group, tmp];
     end

     b_order = [b_order; one_group];

  end

  boot_order = b_order;

   if exist('from_nonstrat','var')
      total_rows0 = from_nonstrat.total_rows0;
      num_group = from_nonstrat.num_group;
      num_subj_lst = from_nonstrat.num_subj_lst;
      num_cond0 = from_nonstrat.num_cond0;
   end

   template = zeros(total_rows0,1);

   for g=1:num_group
      n = num_subj_lst(g);

      for i=bscan
         template(sum(num_cond0*num_subj_lst(1:(g-1)))+(1+(n*(i-1) ):n*i))=1;
      end
   end

   template = find(template);
   boot_order = template(boot_order);

   boot_order0 = repmat([1:total_rows0]', [1 new_num_boot]);
   boot_order0(template,:) = boot_order;
   boot_order = boot_order0;

  return;

