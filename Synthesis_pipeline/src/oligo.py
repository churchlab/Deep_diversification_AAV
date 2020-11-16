
from Bio.Seq import Seq 
from Bio.Data import CodonTable
from Bio.Restriction import BsaI, HindIII, XbaI, BbsI, BsmBI, EcoRV
import re
import random 
import numpy as np

WT_NT_AAV2_T21= 'GACGAAGAGGAAATCAGGACAACCAATCCCGTGGCTACGGAGCAGTATGGTTCTGTATCTACCAACCTCCAGAGAGGCAACAGA'
WT_AA_AAV2_T21 = 'DEEEIRTTNPVATEQYGSVSTNLQRGNR'
class Oligo():
    def __init__(self,serotype="AAV2",debug=False):
        self.serotype=serotype
        if serotype=="AAV2":
            self.wt_nt=WT_NT_AAV2_T21
            self.wt_aa=WT_AA_AAV2_T21
        else:
            print ("this serotype not supported")
        self.full_sequence="sequence not assembled"
        self.barcode="barcode not set"
        self.variant="variant not set"
        self.AA="variant not set"
        self.debug=debug


    def set_restriction_enzymes_to_avoid(self,re_list):
        """These RE sites should never occur anywhere on the oligo"""
        self.enzymes_to_avoid=augment_with_reverse_complements(re_list)

    def set_targeted_restriction_enzymes(self,re_list):
        """These should only appear where we want them"""
        self.targeted_restriction_enzymes=augment_with_reverse_complements(re_list)

    def check_length(self,sequence,limit=230):
        if len(sequence)>limit:
            if self.debug:
                print ("sequence too large")
            return False
        return True     

    def assemble_upstream(self,primer_forward,enzyme_site,overhang,upstream_cut_sites_num=1):
        upstream=primer_forward+enzyme_site+overhang
        if check_re_sites(upstream,self.enzymes_to_avoid)>0: #no restriction enzyme sites of this type should be allowed
            if self.debug:
               print ("site for enzyme that shouldn't be there")    
            return None
        elif check_re_sites(upstream,self.targeted_restriction_enzymes)!=upstream_cut_sites_num: #exactly one site should be here
            if self.debug:
               print ("Too many sites for targeted re")
            return None
        else:
            self.upstream=upstream
            return upstream

    def assemble_downstream(self,pre_barcode,barcode,post_barcode,downstream_cut_sites_num=4):
        downstream=pre_barcode+barcode+post_barcode
        if check_re_sites(downstream,self.enzymes_to_avoid)>0: #no restriction enzyme sites of this type should be allowed
            if self.debug:
               print ("site for enzyme that shouldn't be there")
               if check_re_sites(post_barcode,self.enzymes_to_avoid)>0:
                  print ("primer issue")
               else:
                  print ("barcode issue")

            return None
        elif check_re_sites(downstream,self.targeted_restriction_enzymes)!=downstream_cut_sites_num: #exactly 4 sites should be here
            if self.debug:
               print ("Too many sites for targeted re")
               print (barcode)
               print(check_re_sites(downstream,self.targeted_restriction_enzymes))
               print (downstream)
               if check_re_sites(post_barcode,self.targeted_restriction_enzymes)!=1:
                  print ("primer issue")
               else:
                  print ("barcode issue")
            return None
        else:
            self.downstream=downstream
            self.barcode=barcode
            return downstream

    def assemble_full_sequence(self,sequence,aa_codon):
        success=False
        all_re_sites=self.targeted_restriction_enzymes+self.enzymes_to_avoid

        if self.downstream == None or self.upstream == None:
            if self.debug:
                print ("set upstream and downstream first")
            return None
        else:
            if check_re_sites(sequence,all_re_sites)>0:
               coding_seq=mutate_away_sites(sequence,all_re_sites,aa_codon) 
               assert(str(Seq(coding_seq).translate()) == str(Seq(sequence).translate()))

            else:
                coding_seq=sequence

            full_sequence=self.upstream+coding_seq+self.downstream

            attempts=0
            MAX_ATTEMPS=10

            while check_re_sites(full_sequence,all_re_sites)!=5 and attempts<MAX_ATTEMPS:
                
                full_oligo = self.upstream + coding_seq + self.downstream
                codon_start = len(self.upstream)
                codon_end = len(full_oligo) - len(self.downstream)
                assert(full_oligo[codon_start:codon_end] == coding_seq)

                new_oligo = mutate_away_sites(full_oligo, all_re_sites, aa_codon, np.arange(codon_start, codon_end, 3))
                assert(str(Seq(full_oligo[codon_start:codon_end]).translate()) == str(Seq(new_oligo[codon_start:codon_end]).translate()))

                full_sequence=new_oligo
                attempts+=1

            if attempts>=MAX_ATTEMPS:
                if self.debug:
                    print ("Could not assemble full sequence with {num} attempts".format(num=MAX_ATTEMPS))    
                    print ( check_re_sites(full_sequence,all_re_sites))
                return None
            else:
                #Do some final checks before returning
                if self.check_length(full_sequence):
                    success=True
                    self.full_sequence=full_sequence

        if success:
            return (self.full_sequence, self.barcode)
        else:
            return None


def augment_with_reverse_complements(re_list):
    for enz_site in re_list:
        rc = str(Seq(enz_site).reverse_complement())
        if rc not in re_list:
           re_list.append(rc)
    return re_list

def check_re_sites(nt_seq, re_list):
    count = 0
    for site in re_list:
        site_loc = re.finditer(site, nt_seq.upper())
        for loc in site_loc:
            count+=1
    return count

def check_re_sites_non_overlapping(nt_seq,re_list):
    count=0
    for site in re_list:
        count+=nt_seq.count(site)
    return count 


def mutate_away_sites(seq, enzyme_list, aa_codon, codon_pos = None):
    if codon_pos is None:
        codon_pos = np.arange(0, len(seq), 3) # all codon positions   
    # print [[p, seq[p:p+3]] for p in codon_pos]
    start_seq = seq
    num_sites = check_re_sites(seq, enzyme_list)
    if num_sites>0:
        sites_in = [site for site in enzyme_list if site in seq]
    else:
        return seq
    seq_len = len(seq)    
    current_seq = seq
    for site in sites_in: # for all enzymes
        for site_loc in re.finditer(site, seq):  # for all cut sites
            codon_iter = 0 
            site_start = site_loc.start()
            site_end = site_loc.end()
            overlapping_codon_pos = codon_pos[(codon_pos+3>site_start) & (codon_pos<site_end)] # codons that overlap the site
            for p in overlapping_codon_pos: # for all overlapping codon positions
                codon = current_seq[p:p+3]
                aa = str(Seq(codon).translate())
                syn_codons = aa_codon[aa]
                random.shuffle(syn_codons) # note that this changes the order in aa_codon
                for new_codon in syn_codons:
                    new_seq = current_seq[:p] + new_codon + current_seq[p+3:]
                    # change the sequence if we ever reduce the number of cut sites
                    if check_re_sites(new_seq, enzyme_list) < check_re_sites(current_seq, enzyme_list):
                        current_seq = new_seq
    return current_seq


def backtranslate(mut_aa_seq, aa_codon_table, wt_nt_seq = '', default_wt = False):
    if default_wt:
        if len(wt_nt_seq)%3 != 0 and len(wt_nt_seq)>0:
            sys.exit(['WT nt sequence not divisible by 3', wt_nt_seq, len(wt_nt_seq)])
        wt_aa_seq = str(Seq(wt_nt_seq).translate())
        # wt_mask = mask_wt(wt_aa_seq, mut_aa_seq)
        wt_codons = [wt_nt_seq[max(i-3,0):i] for i in range(len(wt_nt_seq), 0, -3)][::-1]
        # print wt_aa_seq, wt_mask
        # print wt_nt_seq, wt_codons
    else:
        default_wt = False
    mut_codons = [None]*len(mut_aa_seq)
    
    wt_i = 0
    insertion = False

    for i in range(len(mut_aa_seq)):
        # '-' represents a gap, treat it like an insertion so don't increment the wt counter and don't use the wt codon
        if mut_aa_seq[i] == '-': 
            mut_codons[i] = '' 
            gap = True
            print ('-')
        # lowercase represents an insertion, don't increment the wt counter and don't ever use the wt codon
        elif mut_aa_seq[i].islower(): 
            insertion = True
            gap = False
        else:
            insertion = False
            gap = False

        if not gap:
            if not insertion and default_wt and mut_aa_seq[i] == wt_aa_seq[wt_i]:
                mut_codons[i] = wt_codons[wt_i]
                # print mut_aa_seq[i], wt_codons[wt_i], 'wt'
            else:
                aa = mut_aa_seq[i]
                codons = aa_codon_table[aa.upper()]
                random.shuffle(codons) # note that this changes the order in aa_codon_table
                # print aa, aa_codon_table[aa.upper()]
                mut_codons[i] = codons[0] # choose the first one by default
            if not insertion:
                wt_i += 1
        if insertion:
                mut_codons[i] = mut_codons[i].lower()
    # print mut_codons
    mut_nt = ''.join(mut_codons).upper()
    mut_aa_test = str(Seq(mut_nt).translate())
    # print mask_wt(mut_aa_seq, mut_aa_test)
    # print mut_nt
    return mut_nt


