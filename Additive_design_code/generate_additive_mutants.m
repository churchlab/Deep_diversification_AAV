function generate_additive_mutants()
% generates and writes the mutant AAV amino acid sequences using data from
% additive models derived from tissue biodistribution data

    data_path = 'input/tissue_fitness'; % raw fitness data folder
    save_data_path ='input/tissue_data.mat'; % data which has already been parsed
    reload = false;
    if reload || ~exist(save_data_path, 'file')
        names = {'liver', 'kidney', 'spleen', 'skin', 'lung', 'heart', 'brain', 'muscle', 'intestine', 'blood'};
        [plasmid_virus_sub, plasmid_virus_ins, wt, ins_wt, aa, ins_aa, pos] = load_tissue_data(data_path, 'plasmid', 'virus');
        [virus_liver_sub, virus_liver_ins] = load_tissue_data(data_path, 'virus', 'liver');
        [virus_spleen_sub, virus_spleen_ins] = load_tissue_data(data_path, 'virus', 'spleen');
        [plasmid_tissue_sub, plasmid_tissue_ins] = load_tissue_data(data_path, 'plasmid', names);
        [virus_tissue_sub, virus_tissue_ins] = load_tissue_data(data_path, 'virus', names);
        save(save_data_path, 'wt', 'ins_wt', 'aa', 'ins_aa', 'pos', 'names', ...
            'plasmid_virus_sub', 'plasmid_virus_ins', 'virus_liver_sub', 'virus_liver_ins', 'virus_spleen_sub', 'virus_spleen_ins', ...
            'plasmid_tissue_sub', 'plasmid_tissue_ins', 'virus_tissue_sub', 'virus_tissue_ins');
    else
        load(save_data_path, 'wt', 'ins_wt', 'aa', 'ins_aa', 'pos', 'names', ...
            'plasmid_virus_sub', 'plasmid_virus_ins', 'virus_liver_sub', 'virus_liver_ins', 'virus_spleen_sub', 'virus_spleen_ins', ...
            'plasmid_tissue_sub', 'plasmid_tissue_ins', 'virus_tissue_sub', 'virus_tissue_ins');
    end
    
    % output
    mut_file = 'output/additive_mutants.txt';
    path_file ='output/additive_paths.txt';
    log_file = 'output/additive_log.txt';
  
    mid = fopen(mut_file, 'w');
    pid = fopen(path_file, 'w');
    lid = fopen(log_file, 'w');

    %% general    
%     plot_entropy_vs_temp(plasmid_tissue_sub); return % use for setting the ranges of temperatures
    highT = 4; % high enough that it's even (random) sampling
    
    %% generates ~10k random mutants
    naive_scatter_mut = 2:6;
    naive_scatter_num = 500;
    naive_paths_num = 100;
    naive_paths_len = 10;    
    
    %% generates ~24k mutational paths
    path_N_less = 4;
    path_N = 12;
    path_N_more = 20;
    
    path_temps_less = linspace(-2, -1, path_N_less);
    path_temps = linspace(-2, 0, path_N);  
    
    path_thresholds = linspace(-1, 2, path_N);
    path_thresholds_more = linspace(-1, 2, path_N_more);
    
    min_path_thresholds_less = linspace(1, 3, path_N_less);
    min_path_thresholds_more = linspace(0, 3, path_N_more);
    
    num_paths_temp = 1;    
    num_paths_threshold = 1;
    num_paths_min_fit = 1;
    
    %% generates ~36k end points (with no connecting paths)
    point_temps_less = path_temps_less;
    point_temps = path_temps;
    point_thresholds = path_thresholds;
    point_thresholds_more = path_thresholds_more;
    min_point_thresholds_less = min_path_thresholds_less;
    min_point_thresholds_more = min_path_thresholds_more;
    
    num_points = 6;
    num_points_temp = num_points;    
    num_points_threshold = num_points;
    num_points_min_fit = num_points;
    
    max_muts = 2:6;
    max_mut_point_temp = max_muts;    
    max_mut_point_threshold = max_muts;    
    max_mut_point_min_fit = max_muts;   
    
    for do_paths = [0 1] % binary variable for if paths to the endpoint should be included
        for t = 1:10 % across all 10 tissue indices
            t
            figure(101);clf;hold on 
            % select for the tissue alone
            make_library(do_paths, names{t}, plasmid_tissue_sub(:,:,t), plasmid_tissue_ins(:,:,t));  
            % detarget liver
            if ~strcmp(names{t}, 'liver')
                figure(102);clf;hold on 
                make_library(do_paths, sprintf('%s-liver', names{t}), {plasmid_tissue_sub(:,:,t), -virus_liver_sub}, {plasmid_tissue_ins(:,:,t), -virus_liver_ins});
            end
            % detarget spleen
            if ~strcmp(names{t}, 'spleen')
                figure(103);clf;hold on 
                make_library(do_paths, sprintf('%s-spleen', names{t}), {plasmid_tissue_sub(:,:,t), -virus_spleen_sub}, {plasmid_tissue_ins(:,:,t), -virus_spleen_ins});
            end
            % detarget liver and spleen
            if ~strcmp(names{t}, 'spleen') && ~strcmp(names{t}, 'liver')
                figure(104);clf;hold on 
                make_library(do_paths, sprintf('%s-liver-spleen', names{t}), {plasmid_tissue_sub(:,:,t), -virus_liver_sub, -virus_spleen_sub}, {plasmid_tissue_ins(:,:,t), -virus_liver_ins, -virus_spleen_ins});
            end
        end
    end
    
    % generate random mutants with no prior info (for clarity called 'naive' since there are several other types of randomness)
    % all ones --> equal probabilities for each mutants
    blank_sub = ones(size(plasmid_virus_sub)); 
    blank_ins = ones(size(plasmid_virus_ins));
    
    % generate random endpoints
    for num_variant = naive_scatter_mut
        generate_mutants('naive', 'random', naive_scatter_num, {'end_point'}, highT, -Inf, -Inf, [0 1 2], {blank_sub}, {blank_ins}, num_variant);
    end
    % generate some longer random paths
    generate_mutants('naive', 'random', naive_paths_num, {'rand_path:1'}, highT, -Inf, -Inf, [0 1 2], {blank_sub}, {blank_ins}, naive_paths_len);

    fclose(mid);     
    fclose(pid);     
    fclose(lid); 
    
    mut = readtable(mut_file, 'HeaderLines', 0);
    number_of_sequences = numel(unique(mut)) % prints the final number of unique amino acid sequences generated
    
    function make_library(do_paths, name, sub, ins)
        multi_property = true; %true for multi-property selection
        if ~iscell(sub)
            sub = {sub};
            multi_property = false;
        end        
        if ~iscell(ins)
            ins = {ins};
        end
        if multi_property
            if do_paths            
                generate_mutants(name, 'temp', num_paths_temp, {'fast_path',}, path_temps_less, -Inf, nan, [0 1], sub, ins);
                generate_mutants(name, 'threshold', num_paths_threshold, {'fast_path'}, highT, path_thresholds_more, nan, [0 1], sub, ins);
                generate_mutants(name, 'min_fit', num_paths_min_fit, {'min_path'}, highT, -Inf, min_path_thresholds_less, [0 1], sub, ins);
            else
                generate_mutants(name, 'temp', num_points_temp, {'end_point'}, point_temps_less, -Inf, nan, [0 1], sub, ins, max_mut_point_temp);
                generate_mutants(name, 'threshold', num_points_threshold, {'end_point'}, highT, point_thresholds_more, nan, [0 1], sub, ins, max_mut_point_threshold);            
                generate_mutants(name, 'min_fit', num_points_min_fit, {'end_point'}, highT, -Inf, min_point_thresholds_less, [0 1], sub, ins, max_mut_point_min_fit);
            end
        else
            if do_paths            
                generate_mutants(name, 'temp', num_paths_temp, {'fast_path', 'rand_path:1'}, path_temps, -Inf, nan, [0 1], sub, ins);
                generate_mutants(name, 'threshold', num_paths_threshold, {'fast_path', 'rand_path:1'}, highT, path_thresholds, nan, [0 2], sub, ins);
                generate_mutants(name, 'min_fit', num_paths_min_fit, {'min_path'}, highT, -Inf, min_path_thresholds_more, [0 2], sub, ins);
            else
                generate_mutants(name, 'temp', num_points_temp, {'end_point'}, point_temps, -Inf, nan, [0 1], sub, ins, max_mut_point_temp);
                generate_mutants(name, 'threshold', num_points_threshold, {'end_point'}, highT, point_thresholds, nan, [0 2], sub, ins, max_mut_point_threshold);            
                generate_mutants(name, 'min_fit', num_points_min_fit, {'end_point'}, highT, -Inf, min_point_thresholds_more, [0 2], sub, ins, max_mut_point_min_fit);
            end
        end
    end
    function mask = mask_else_wt(mask, wt)
        % fix a mask when a given position has no possible aa values by
        % setting the aa choice to wt
        if size(mask, 1) < size(wt,1)
           row_diff = size(wt,1) - size(mask, 1);
           rows = size(wt,1)+ [1:row_diff];
           mask(rows,:) = false;
        end
        for j = 1:size(mask,2)
            if sum(mask(:,j)) == 0
                mask(:,j) = wt(:,j);
            end
        end
    end
    function plot_entropy_vs_temp(data)
        %% scan range of temperatures to find min and mix ranges
        Ts = linspace(-3,2,500);
        % plot diversity for each tissue
        figure(11);clf;
        N = numel(names);
        nr = 5;
        nc = 2;
        for i = 1:N
            for ti = 1:numel(Ts)
                T = 10^Ts(ti);
                p = calc_probs(data(:,:,i), T);
                S(ti,:) = calc_entropy(p);
            end
            subplot(nr,nc,i)
            ph = plot(Ts, 2.^S);
            axis tight
            ylim([1 20]);
            title(names{i});
            set(gca, 'fontsize', 12, 'xtick', -10:10, 'ytick', -0:5:20)
        end
        
        % plotting where we achieve 10aa diversity for each position
        figure(12);clf;
        nr = 4;
        nc = 7;
        for ti = 1:numel(Ts)
            T = 10^Ts(ti);
            p = calc_probs(data(:,:,10), T);
            S(ti,:) = calc_entropy(p);
        end
        for j = 1:size(S,2)
            subplot(nr,nc,j)
            ph = plot(Ts, 2.^S(:,j)); hold on
            x(j) = Ts(find((2.^S(:,j))>10, 1));
            plot([1 1]*x(j), [1 20], 'r')
%             ylim([0 max(S(:))]);
            axis tight
            xlim([-1 1])
            ylim([1 20]);
            title(sprintf('%d, %.2f', j, x(j)));
%             set(gca, 'fontsize', 12, 'xtick', -10:10, 'ytick', -10:10)
            set(gca, 'fontsize', 12, 'xtick', -10:10, 'ytick', -0:5:20)
        end
        % more deleterious positions require higher temperature to achieve
        % diversity
        figure(13);clf;
        avg = mean(data(:,:,1),1);
        plot(avg, x, '-o');
        
%         savefig('figs/S_vs_T.pdf', 'pdf');
    end
    function generate_mutants(name, method, num_variants, types, temps, thresholds, min_fit, enable_insertions, subs, inss, max_muts)
        if ~exist('max_muts', 'var')
            max_muts = Inf;
        end
        if ~iscell(types)
            types = {types};
        end   
        N = numel(subs);
        full_name = sprintf('%s, %s', name, method);
        WT = char('DEEEIRTTNPVATEQYGSVSTNLQRGNR');
        mut_sub = WT; % initially WT 
        mut_ins = WT; mut_ins(:) = '-'; % initially all gaps
        for max_mut = max_muts
            for fmin = min_fit
                for ins_ok = enable_insertions % 0=no insertions, 1 = only in second half of tile, 2 = anywhere ok
                    for T = temps
                        for threshold = thresholds
                            
                            % build sub and allow_sub from all the fitness matrices
                            %
                            allow_sub_temp = true(size(subs{1}));
                            sub = zeros(size(subs{1}));
                            for mi = 1:N
                                allow_sub_temp = allow_sub_temp & subs{mi}>threshold; %enforce a minimum threshold
                                sub = sub+subs{mi}/N; % average across all fitness matrices
                            end
                            allow_sub = mask_else_wt(allow_sub_temp, wt);
                            
                            % same for ins and allow_ins
                            allow_ins_temp = true(size(inss{1}));
                            ins = zeros(size(inss{1}));
                            for mi = 1:N
                                ins = ins+inss{mi}/N;
                            end 
                            if ins_ok
                                for mi = 1:N
                                    allow_ins_temp = allow_ins_temp & inss{mi}>threshold;
                                end
                                if ins_ok < 2
                                    allow_ins_temp(:, 1:14) = false; % disallow insertions in first half of the tile
                                end
                            else
                                allow_ins_temp = false(size(ins)); % disallow all insertions
                            end
                            allow_ins = mask_else_wt(allow_ins_temp, ins_wt); % allow choice of no insertions at each position (making it WT)                            
                            
                            % some plots for debugging purposes
%                             figure(112);clf; hold on
%                             plot(subs{1}(:), subs{2}(:), 'o');
%                             plot(subs{1}(allow_sub(:)), subs{2}(allow_sub(:)), '.r');
%                             axis equal
%                             figure(113);clf;imagesc(cat(1,allow_ins, allow_sub))

                            p_sub = calc_probs(sub, 10^T, allow_sub);
                            p_ins = calc_probs(ins, 10^T, allow_ins);

                            if strcmp(method, 'min_fit')
                                [ins_order, sub_order, ins_indices, sub_indices, ins_fit, sub_fit] = sample_seq_min_fit(fmin, num_variants, ins, sub, allow_ins, allow_sub);
                                min_fit_var.threshold = fmin;
                                min_fit_var.pos_order = {};

                                ins_seqs = ins_aa(ins_indices);
                                sub_seqs = aa(sub_indices);
                                % here we send the indices instead of fitness
                                % scores, in order to keep track of the order
                                [masks, mut_pos, num, full_wt, full_mut, full_order, sum_scores, ins_seqs, sub_seqs] = wt_mask_ins(WT, ins_seqs, sub_seqs, ins_order, sub_order, max_mut, true);                            
                                min_fit_var.pos_order = full_order;

                                % redo to get the scores instead of the order
                                [masks, mut_pos, num, full_wt, full_mut, full_scores, sum_scores] = wt_mask_ins(WT, ins_seqs, sub_seqs, ins_fit, sub_fit, max_mut);
                            else
                                min_fit_var.threshold = nan;                            
                                [ins_indices, ins_fit] = sample_seq(p_ins, num_variants, ins);
                                [sub_indices, sub_fit] = sample_seq(p_sub, num_variants, sub);
                                ins_seqs = ins_aa(ins_indices);
                                sub_seqs = aa(sub_indices);
                                [masks, mut_pos, num, full_wt, full_mut, full_scores, sum_scores] = wt_mask_ins(WT, ins_seqs, sub_seqs, ins_fit, sub_fit, max_mut);
                            end

                            [num_mut, mut_score] = write_paths(full_name, types, full_wt, full_mut, full_scores, T, threshold, min_fit_var, ins_ok, pid, mid, lid);
                            
                            % this will plot the predicted mutant fitness values vs dist from WT                            
                            plot(num_mut, mut_score, 'o');
                            
                            % some plots for debugging purposes
%                             hist(sum_scores);
% 
%                             figure(1);clf;
%                             subplot(2,1,1); imagesc(p_ins); colorbar;
%                             subplot(2,1,2); imagesc(p_sub); colorbar; 
%             
%                             figure(2);clf;
%                             subplot(2,1,1); imagesc(allow_ins); colorbar;
%                             subplot(2,1,2); imagesc(allow_sub); colorbar;
%             
%                             figure(3);clf;
%                             subplot(2,1,1); plot(ins_fit');title('Ins');
%                             subplot(2,1,2); plot(sub_fit');title('Sub');
                        end
                    end
                end
            end
        end
        function [o1, o2, i1, i2, fit1, fit2] = sample_seq_min_fit(fmin, k, f1, f2, f1_allow, f2_allow)
            n = size(f1,2);
            m = size(f1,1); % number of choices for first matrix (insertions)
            
            o1 = zeros(k,n);
            o2 = zeros(k,n);
            i1 = zeros(k,n);
            i2 = zeros(k,n);
            fit1 = zeros(k,n);
            fit2 = zeros(k,n);
            for i = 1:k
                current_f = 0;
                temp_f1 = f1;
                temp_f2 = f2;
                temp_f1(~f1_allow(:)) = nan;
                temp_f2(~f2_allow(:)) = nan;
                
%                 figure(2);clf;imagesc(cat(1,temp_f1,temp_f2));
                for j = 1:(2*n)
                    
                    p1 = zeros(size(f1));
                    p2 = zeros(size(f2));
                    % mark all currently available mutations
                    p1(isnan(temp_f1)) = -1;
                    p2(isnan(temp_f2)) = -1; 
                    % first try to be greater than current fitness
                    p1(~isnan(temp_f1) & current_f+f1 > fmin & f1~=0) = 1;
                    p2(~isnan(temp_f2) & current_f+f2 > fmin & f2~=0) = 1;
                    % if no more choices greater than current, enable equal to current (allows wt)                    
                    if isempty(find([p1(:); p2(:)]==1,1))
                        p1(~isnan(temp_f1) & current_f+f1 >= fmin) = 1;
                        p2(~isnan(temp_f2) & current_f+f2 >= fmin) = 1;
                    end
                    % if no more choices greater or equal to current,
                    % enable equal to WT
                    if isempty(find([p1(:); p2(:)]==1,1))
                        p1(~isnan(temp_f1) & f1 == 0) = 1;
                        p2(~isnan(temp_f2) & f2 == 0) = 1;
                    end
                    % choose a mutation
                    p_all = cat(1,p1,p2);                    
                    choices = find(p_all(:)==1);                    
                    chosen = choices(randperm(numel(choices),1));
                    p_all(chosen) = 2; % marks it as chosen
%                     figure(1);clf;imagesc(p_all);
                    % unstack the matrices
                    p1 = p_all(1:m,:)==2;
                    p2 = p_all(m+1:end,:)==2;
                    
                    r1 = find(sum(p1,2),1);
                    r2 = find(sum(p2,2),1);
                    
                    c1 = find(sum(p1,1),1);
                    c2 = find(sum(p2,1),1);
                    
                    if c1
                        o1(i, c1) = j;
                        i1(i, c1) = r1;
                        fit1(i, c1) = f1(r1,c1);
                        temp_f1(:,c1) = nan;
                        current_f = current_f + f1(r1,c1);
                    else
                        o2(i, c2) = j;
                        i2(i, c2) = r2;
                        fit2(i, c2) = f2(r2,c2);
                        temp_f2(:,c2) = nan;
                        current_f = current_f + f2(r2,c2);
                    end
                end
            end
        end
        function [indices, fit] = sample_seq(p, k, f)
            n = size(p,2);
            indices = zeros(k, n);
            fit = zeros(k, n);
            for j = 1:n
                indices(:,j) = sample_aa(p(:,j), k);
                col = f(:,j);
                fit(:,j) = col(indices(:,j));
            end
        end
        function i = sample_aa(col_p, k)
            if ~exist('k', 'var')
                k = 1;
            end
            i = randsample(numel(col_p), k, true, col_p); 
        end
    end
    function probs = calc_probs(f, T, allow)
        if exist('allow', 'var')
            f(~allow(:)) = -Inf;
        end
        pp = 2.^(f/T);
        z = sum(pp,1);
        Z = repmat(z, size(pp,1),1);
        probs = pp./Z;
        choices = sum(probs>0,1);
        % if no choices due to very low prob (rounding error), set highest f mutants to prob=1
        for c = find(choices==0)
            f_pos = f(:,c);
            probs(:, c) = f_pos == max(f_pos);
        end
    end
    function S = calc_entropy(p)
        s = -p.*log2(p);
        s(p(:)==0) = 0;
        S = sum(s,1);
    end
end