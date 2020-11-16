function [num_muts, mut_scores] = write_paths(name, all_types, full_wt, full_mut, full_scores, T, threshold, min_fit_var, ins_ok, ph, mh, lh)
    
    num_muts = [];
    mut_scores = [];
    for i = 1:size(full_wt,1)
        wt = full_wt{i};
        seq = full_mut{i};
        scores = full_scores{i}; 
        min_fit = min_fit_var.threshold;
        
        for ti = 1:numel(all_types)
            type = all_types{ti};
            if contains(type, 'rand_path') || contains(type, 'sub_sample')>0
               temp = strsplit(type, ':');
               type = temp{1};
               num_paths = str2num(temp{2});
            end
            switch type
                case {'full_sample', 'end_point'}
                    [paths, path_id, mut_id] = generate_paths(type, wt, seq);
                case 'fast_path'
                    [paths, path_id, mut_id] = generate_paths(type, wt, seq, scores);
                case 'min_path'
                    pos_order = min_fit_var.pos_order{i};
                    [paths, path_id, mut_id] = generate_paths(type, wt, seq, pos_order);
                case {'rand_path', 'sub_sample'}
                    [paths, path_id, mut_id] = generate_paths(type, wt, seq, num_paths);
                otherwise
                    type
                    error('incorrect type');
            end

            [~, pos, n] = wt_mask(wt, seq);  
            for pi = 1:numel(paths)
                p = paths(pi);
                mut = wt; % start out with wt sequence
                
                pos_mut = dec2bin(p, max(1,n)); % identify position of mutations to be made
%                 pos_mut
                ii = pos_mut=='1'; % convert to indices
                num_mut(pi) = sum(pos_mut=='1');
%                 pos
%                 ii
                mut(pos(ii)) = seq(pos(ii)); % make the chosen mutations
                mut_mask = wt_mask(wt, mut);
                mut_score(pi) = sum(scores(pos(ii)));
                
                fprintf(ph, sprintf('%s, %8s, %2.2f, %2.2f, %2.2f, %d, %s, %s, %s, %2d, %2d, %2d, %2d, %.4f, %s, %d\n', ...
                    name, type, T, threshold, min_fit, ins_ok, seq, mut, mut_mask, i, path_id(pi), mut_id(pi), num_mut(pi), mut_score(pi), pos_mut, p));
                fprintf(mh, sprintf('%s\n',strrep(mut, '-', '')));
            end
            final_num = numel(unique(paths));
            switch type
                case {'full_sample' , 'sub_sample', 'end_point'}
                    fprintf(lh, '%s, %s, %.4f, %.4f, %d, %d sequences\n', name, type, T, threshold, ins_ok, final_num);
                case {'rand_path', 'fast_path', 'min_path'}
                    num_paths = int32(numel(paths)/(n+1));
                    fprintf(lh, '%s, %s, %.4f, %.4f, %d, %d paths, %d->%d sequences\n', name, type, T, threshold, ins_ok, num_paths, numel(paths), final_num);
                otherwise
                    error('incorrect type');
            end
            num_muts = [num_muts num_mut];
            mut_scores = [mut_scores mut_score];
        end
    end
end