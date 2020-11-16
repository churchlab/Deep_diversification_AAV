function [paths, path_id, mut_id] = generate_paths(type, wt, seq, var1)
    [mask, pos, n] = wt_mask(wt, seq);
    num_mut = 2^n-1;
    switch type
        case 'full_sample'
            paths = 0:num_mut;
            path_id = ones(1,num_mut+1);
            mut_id = [0:num_mut];
        case 'rand_path'
            index = 1;
            path_index = 1;
            num = var1;
            path_pos = [];
            for ii = 1:num % for each path
                pos_order = randperm(n);
                if size(path_pos,1) < size(unique([ path_pos ; pos_order], 'rows'),1) % only if its a new path
                    path_pos(path_index,:) = pos_order;            
                    bin = repmat('0', 1, n); % char array all 0's, length n
                    % add wt
                    paths(index) = 0;
                    path_id(index) = path_index;
                    mut_id(index) = 0;                
                    % add mutants
                    index = index + 1;
                    for jj = 1:numel(pos_order)
                        bin(pos_order(jj)) = '1';
                        p = bin2dec(bin);
                        paths(index) = p;
                        path_id(index) = path_index;
                        mut_id(index) = jj;                    
                        index = index + 1;
                    end
                    path_index = path_index + 1;
                end
            end    
        case 'fast_path'
            index = 1; 
            bin = repmat('0', 1, n); % char array all 0's, length n
            % add wt
            paths(index) = 0;
            path_id(index) = 1;
            mut_id(index) = 0;
            % add mutants
            index = index + 1;
            scores = var1(1,pos);
            [~, pos_order] = sort(scores, 'descend');
            for jj = 1:numel(pos_order)
                bin(pos_order(jj)) = '1';
                p = bin2dec(bin);
                paths(index) = p;
                path_id(index) = 1;
                mut_id(index) = jj;
                index = index + 1;
            end
        case 'min_path'
            index = 1; 
            bin = repmat('0', 1, n); % char array all 0's, length n 
            % add wt
            paths(index) = 0;
            path_id(index) = 1;
            mut_id(index) = 0;
            % add mutants
            index = index + 1;
            [~, pos_order] = sort(var1(1,pos));
            for jj = 1:numel(pos_order)
                bin(pos_order(jj)) = '1';
                p = bin2dec(bin);
                paths(index) = p;
                path_id(index) = 1;
                mut_id(index) = jj;
                index = index + 1;
            end
        case 'sub_sample'
            num = var1;
            paths = sort(randperm(num_mut, num));
            path_id = ones(1,num);
            mut_id = 1:num;            
        case 'end_point'
            paths = num_mut;
            ids{1} = '1';
            path_id = 1;
            mut_id = 1;
        otherwise
            type
            error('incorrect type');
    end
end