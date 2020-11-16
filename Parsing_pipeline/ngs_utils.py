import os, sys, re , csv , gzip 
import subprocess
from collections import defaultdict

def check_directory_exists(path):
    if not os.path.exists(path): 
        print(path)
        os.makedirs(path)
    else:
        print('already exists: ', path)
    
            
def bash_command(cmd):
    subprocess.Popen(['/bin/bash', '-c', cmd])            

def most_common(lst):
    data = Counter(lst)
    return data.most_common(1)[0][0]

def tabulate(d, key, count=1):
    # add 1 to dict value for this key, or add a new key to the dictionary
    if key in d:
        d[key] += count
    else:
        d[key] = count    
    
def tabulate_bar_seq(d, bar, seq, count=1):
    if bar in d.keys():
        if seq in d[bar]:
            d[bar][seq] += count
        else:
            d[bar].update({seq:count})
    else:
        d[bar] = {seq:count}

def write_seq_counts(dict, outputfile, header='seq, count\n'):
    # write sequences and counts to a text file
    print(outputfile)
    with open(outputfile, 'w') as f:
        f.write(header)
        for w in sorted(dict, key=dict.get, reverse=True):
            f.write('{seq}, {count}\n'.format(seq=w, count=str(dict[w])))    

def write_bar_seq_counts(bar_seq_dict, bar_dict, outputfile, header="bar, seq, count\n"):        
    print(outputfile)
    with open(outputfile,'w') as f:
        f.write(header)
        for bar in sorted(bar_dict, key=bar_dict.get, reverse=True):
            if bar in bar_seq_dict:
                seq_dict = bar_seq_dict[bar]
                for seq in sorted(seq_dict, key=seq_dict.get, reverse=True):
                    f.write("{bar}, {seq}, {count}\n".format(bar=bar, seq=seq, count=seq_dict[seq]))

def get_libary_region(seq, scores, fwd, rev):
    # first look for exact match
    fwd_start = seq.find(fwd)
    rev_start = seq.find(rev)
    if fwd_start is -1:
        # print seq 
        # print fwd 
        # print 'no fwd'
        return False, '', []
    if rev_start is -1:
        # print 'no rev'
        return False, '', []

    fwd_spacer = 0
    rev_spacer = 0
    library_start = fwd_start + len(fwd) + fwd_spacer
    library_end = rev_start - rev_spacer
    # print library_start
    # print library_end
    
    if  library_end <= library_start:
        return False, '', []

    if len(seq) < library_end:
        # print 'read too short'Delseq_24_24_L01.txt
        return False, '', []

    library_seq = seq[library_start:library_end]
    library_scores = scores[library_start:library_end]

    return True, library_seq, library_scores     

def get_barcode(seq,scores):
    rev = 'GTGTGGCTGC' ## these priemrs are on forward strand ### 
    fwd = 'ACGCGTAGGA' ## forward strand as well 
    bar_fwd_start = seq.find(fwd)
    bar_rev_start = seq.find(rev)
    if bar_fwd_start is -1:
        # print seq 
        # print fwd 
        # print 'no fwd'
        return False, '', []
    if bar_rev_start is -1:
        # print 'no rev'
        return False, '', []
    barcode_start = bar_fwd_start + len(fwd)
    barcode_end = bar_rev_start
    if  barcode_end <= barcode_start:
        return False, '', []

    if len(seq) < barcode_end:
        # print 'read too short'
        return False, '', []
    barcode_seq = seq[barcode_start:barcode_end]
    barcode_scores = scores[barcode_start:barcode_end]

    return True, barcode_seq, barcode_scores