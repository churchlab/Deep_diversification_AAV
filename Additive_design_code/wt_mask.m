function [mask, mut_pos, num_diff] = wt_mask(wt, seq)
    mask = seq;
    if size(seq,1) == 1
        mask(seq==wt) = '_';
        mut_pos = find(wt~=seq);
        num_diff = sum(mask~='_');
    else
        for i = 1:size(seq,1)
            mask(i, seq(i,:)==wt) = '_';
            mut_pos{i} = find(wt~=seq(i,:));
            num_diff(i) = sum(mask(i,:)~='_');
        end
    end
end
