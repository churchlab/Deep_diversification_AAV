function [u_mask, u_mut_pos, u_num_diff, u_full_wt, u_full_mut, u_full_scores, u_sum_scores, ins, sub] = wt_mask_ins(wt, ins, sub, ins_score, sub_score, max_mut, preserve_order)
    % create mask for subs and insertions
    if ~exist('preserve_order', 'var')
        preserve_order = false;
    end
    % first for each mutation, need to generate a new WT and mutant seq of 
    % the correct final length
    if size(ins,1) ~= size(sub,1)
        error('size of ins and seq not equal')
    end
    ins = lower(ins);
    ins_wt = wt; ins_wt(:) = '-'; % all dashes
    
    N = size(ins,1);
    full_wt = cell(N,1);
    full_mut = cell(N,1);
    full_scores = cell(N,1);
    sum_scores = zeros(N,1);
    
    for i = 1:size(ins,1)
        mut_sub = sub(i,:);
        mut_ins = ins(i,:);
        ins_pos = mut_ins ~= '-';
        
        % order of rows is important: ins on top row since insertions come
        % at the position before the substitution
        %%
        pos = [ins_pos; ones(1,numel(wt))];
        wt_nt = [ins_wt;wt];
        mut_nt = [mut_ins;mut_sub];
        nt_scores = [ins_score(i,:); sub_score(i,:);];
        
        % revert some mutations to WT if over max_mut number
        all_muts = find(wt_nt ~= mut_nt);        
        if max_mut < numel(all_muts) %only allow a subset of these mutations
            if preserve_order
                chosen_muts = all_muts(nt_scores(all_muts) <= max_mut);
            else
                chosen_muts = all_muts(randperm(numel(all_muts), max_mut));
            end
            revert = true(size(mut_nt));
            revert(chosen_muts) = false;
            mut_nt(revert) = wt_nt(revert);
            nt_scores(revert) = 0;
            
            % revert the original array
            ins(i,:) = mut_nt(1,:);
            sub(i,:) = mut_nt(2,:);
        end
        mut_ins = mut_nt(1,:); 
        ins_pos = mut_ins ~= '-';
        pos = [ins_pos; ones(1,numel(wt))];
        
        [keep, ix] = sort(pos(:));        
        temp_wt = wt_nt(ix)';
        temp_mut = mut_nt(ix)';
        temp_scores = nt_scores(ix)';
        
        % remove blank insertions
        full_wt{i} = temp_wt(keep'==1);
        full_mut{i} = temp_mut(keep'==1);
        full_scores{i} = temp_scores(keep'==1);
        sum_scores(i) = sum(temp_scores(keep'==1));
    end    
    
    mask = cell(N,1);
    mut_pos = cell(N,1);
    num_diff = zeros(N,1);
    for i = 1:size(ins,1)
        [mask{i}, mut_pos{i}, num_diff(i)] = wt_mask(full_wt{i}, full_mut{i});
    end 
    
    % only choose unique mutants
    [~, ia] = unique(mask);
    N = numel(ia);
    u_mask = cell(N,1);
    u_mut_pos = cell(N,1);
    u_num_diff = zeros(N,1);
    u_full_wt = cell(N,1);
    u_full_mut = cell(N,1);
    u_full_scores = cell(N,1);
    u_sum_scores = zeros(N,1);
    for i = 1:N
        index = ia(i);
        u_mask{i} = mask{index};
        u_mut_pos{i} = mut_pos{index};
        u_num_diff(i) = num_diff(index);
        u_full_wt{i} = full_wt{index};
        u_full_mut{i} = full_mut{index};
        u_full_scores{i} = full_scores{index};
        u_sum_scores(i) = sum_scores(index);
    end
end