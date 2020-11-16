function [f_sub, f_ins, wt, ins_wt, aa, ins_aa, pos] = load_tissue_data(data_path, before, after, subtract_wt)
    if ~exist('subtract_wt', 'var')
        subtract_wt = true;
    end
    %% load tissue data, subtract off wt fitness
    if ~iscell(after) % make after a cell array if not already
        after = {after}; 
    end
    num_files = numel(after);
    % read the first file, get size info
    sub = readtable(sprintf('%s/sub_%s-%s.txt', data_path, 'plasmid', 'virus'));
    cols = sub.abs_pos - min(sub.abs_pos) + 1;
    aa = 'ILAVGMFYWEDQNHCRKSTP*';
    pos = unique(sub.abs_pos);

    num_rows = numel(aa);
    num_cols = numel(pos);
    wt = false(num_rows, num_cols);
    f_sub = zeros(num_rows, num_cols, num_files);
    f_ins = zeros(num_rows, num_cols, num_files);
    % read in each file
    for k = 1:num_files 
        f = after{k};
        [before '-' f]
        sub = readtable(sprintf('%s/sub_%s-%s.txt', data_path, before, f));
        ins = readtable(sprintf('%s/ins_%s-%s.txt', data_path, before, f));
        for i = 1:numel(cols)
            current_aa = char(sub.aa(i));
            row = find(aa == current_aa);
            col = cols(i);
%                 [row col k]
            f_sub(row, col, k) = sub.log2FoldChange(i);
            f_ins(row, col, k) = ins.log2FoldChange(i);
            wt(row, col) = sub.wt(i)==1;
        end
        f_temp = f_sub(:,:,k);
        f_wt = f_temp(find(wt(:),1));
        if subtract_wt
            f_sub(:,:,k) = f_sub(:,:,k) - f_wt;
            f_ins(:,:,k) = f_ins(:,:,k) - f_wt;
        end
    end
    
    % remove stop codons
    good = aa~='*';
    aa = aa(good);
    wt = wt(good,:);
    f_sub = f_sub(good,:,:);
    f_ins = f_ins(good,:,:);
    
    ins_aa = [aa '-']; % dash means no insertion, this is the 'WT'
    f_ins(end+1,:,:) = 0; % add wt fitness at bottom of ins matrix
    ins_wt = false(size(wt)); 
    ins_wt(end+1,:) = true; % last position is wt, set to true
end