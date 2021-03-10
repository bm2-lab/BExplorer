# -*- coding: utf-8 -*-

import subprocess, logging, time, re, sys, argparse, os
import pandas as pd, numpy as np
import regex as regex
logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s]:  %(message)s')

def InputParser():
    parser = argparse.ArgumentParser(description = 'Predicting the best gRNA at any location in the human genome for the base editor')
    editor_list = ['BE1', 'BE2', 'BE3', 'HF-BE3', 'BE4', 'BE4max', 'BE4-GAM', 'YE1-BE3 ', 'EE-BE3', 'YE2-BE3',
                   'YEE-BE3', 'VQR-BE3', 'VRER-BE3', 'SaBE3', 'SaBE4', 'SaBE4-Gam', 'SaKKH-BE3', 'Target-AID',
                   'BE-PLUS', 'ABE7.9', 'ABE7.10', 'xABE', 'ABESa', 'VQR-ABE', 'VRER-ABE', 'SaKKH-ABE']
    parser.add_argument('--chr', required=True, choices=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', 
        '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y'] , 
        help='input your target chromosome number(X, Y need to be capitalized)')
    parser.add_argument('--pos', required=True, type=int, help='input your target position number (GRCh38.p12)')
    parser.add_argument('--transtype', required=True, choices=['c2t', 'a2g'], help='c2t or a2t. \n '
                                                                                   '*****c2t(CBE): regard the base at your target position is C, '
                                                                                   'and you want to convert C to T. \n'
                                                                                   '******a2t(ABE):regard the base at your target position is A, and you want to convert A to G')
    parser.add_argument('--editor', required=True, choices= editor_list, help='choose one base editor.\n'
                                                                              '*****CBE: [ BE1, BE2, BE3, HF-BE3, BE4, BE4max, BE4-GAM, YE1-BE3, EE-BE3, YE2-BE3, YEE-BE3, '
                                                                              'VQR-BE3, VRER-BE3, SaBE3, SaBE4, SaBE4-Gam, SaKKH-BE3, Target-AID, BE-PLUS ]. \n'
                                                                              '*****ABE: [ ABE7.9, ABE7.10, xABE, ABESa, VQR-ABE, VRER-ABE, SaKKH-ABE ]')

    args = parser.parse_args()
    return args

def editors_store(args):
    editors_dic = {'BE1': ('NGG', 3, 23, [4,8], 'C'), 'BE2': ('NGG', 3, 23, [4,8], 'C'),
                   'BE3': ('NGG', 3, 23, [4,8], 'C'),
                   'HF-BE3': ('NGG', 3, 23, [4,8], 'C'), 'BE4': ('NGG', 3, 23, [4,8], 'C'),
                   'BE4max': ('NGG', 3, 23, [4,8], 'C'), 'BE4-GAM': ('NGG', 3, 23, [4,8], 'C'),
                   'YE1-BE3': ('NGG', 3, 23, [4,7], 'C'), 'EE-BE3': ('NGG', 3, 23, [5,6], 'C'),
                   'YE2-BE3': ('NGG', 3, 23, [5,6], 'C'), 'YEE-BE3': ('NGG', 3, 23, [5,6], 'C'),
                   'VQR-BE3': ('NGAN', 4, 24, [4,11], 'C'), 'VRER-BE3': ('NGCG', 4, 24, [3,10], 'C'),
                   'SaBE3': ('NNGRRT', 6, 27, [3,12], 'C'), 'SaBE4': ('NNGRRT', 6, 27, [3,12], 'C'),
                   'SaBE4-Gam': ('NNGRRT', 6, 27, [3,12], 'C'), 'SaKKH-BE3': ('NNNRRT', 6, 27, [3,12], 'C'),
                   'Target-AID': ('NGG', 3, 23, [2,8], 'C'), 'BE-PLUS': ('NGG', 3, 23, [4,16], 'C'),
                   'ABE7.9': ('NGG', 3, 23, [4,9], 'A'), 'ABE7.10': ('NGG', 3, 23, [4,7], 'A'),
                   'xABE': ('NG', 2, 22, [4,7], 'A'), 'ABESa': ('NNGRRT', 6, 27, [6,12], 'A'),
                   'VQR-ABE': ('NGA', 3, 23, [4,8], 'A'), 'VRER-ABE': ('NGCG', 4, 24, [4,8], 'A'),
                   'SaKKH-ABE': ('NNNRRT', 6, 27, [8,13], 'A')}
    cbe_editors = ['BE1', 'BE2', 'BE3', 'HF-BE3', 'BE4', 'BE4max', 'BE4-GAM', 'YE1-BE3 ', 'EE-BE3', 'YE2-BE3',
                   'YEE-BE3', 'VQR-BE3', 'VRER-BE3', 'SaBE3', 'SaBE4', 'SaBE4-Gam', 'SaKKH-BE3', 'Target-AID',
                   'BE-PLUS']
    abe_editors = ['ABE7.9', 'ABE7.10', 'xABE', 'ABESa', 'VQR-ABE', 'VRER-ABE', 'SaKKH-ABE']
    editor_info = editors_dic[args.editor]
    return editor_info, cbe_editors, abe_editors

def check_softwares():
    softwares_list = ['vcftools','cas-offinder', 'python2.7','R']
    sw_exist_dict = {}
    for sw in softwares_list:
        rs = subprocess.getoutput('which {}'.format(sw))
        sw_exist_dict[sw] = rs
    for i in sw_exist_dict.keys():
        if sw_exist_dict[i] == '':
            logging.error('{} is not installed correctly or not in environment variable !'.format(i))
    if '' in list(sw_exist_dict.values()):
        sys.exit(1)
    op = open('supportfile/testR_Package.R', 'w')
    op.write('library()')
    op.close()
    Rtest = subprocess.getoutput('Rscript supportfile/testR_Package.R')
    if 'RobustRankAggreg' not in Rtest:
        logging.error("Please install R package 'RobustRankAggreg'!")
        sys.exit(1)
    return

def ckeck_File_Dir():
    if 'supportfile' not in os.listdir():
        logging.error('/BExplore/supportfile/ does not exist! Please check if the programm file is complete!')
        sys.exit(1)

    if 'genomeDATA' in os.listdir() and 'snpDATA' in os.listdir():
        if 'ignoreme' in os.listdir('genomeDATA'):
            os.remove('genomeDATA/ignoreme') 
        if 'ignoreme' in os.listdir('snpDATA'):
            os.remove('snpDATA/ignoreme') 
        if len(os.listdir('genomeDATA')) == 0:
            logging.error('Please put the human genome file (fna format) into /BExplore/genomeDATA/')
            sys.exit(1)
        elif len(os.listdir('genomeDATA')) > 1:
            logging.error('Please make sure there is only one file(human genome file) in /BExplore/genomeDATA/')
            sys.exit(1)

        if len(os.listdir('snpDATA')) == 0:
            logging.error('Please put the human SNP file (vcf format) into /BExplore/snpDATA/')
            sys.exit(1)
        elif len(os.listdir('snpDATA')) > 1:
            logging.error('Please make sure there is only one file(human SNP file) in /BExplore/snpDATA/')
            sys.exit(1)

    if 'genomeDATA' not in os.listdir():
        os.mkdir('genomeDATA')
        logging.error('Please put the human genome file (fna format) into /BExplore/genomeDATA/')
    if 'snpDATA' not in os.listdir():
        os.mkdir('snpDATA')
        logging.error('Please put the human SNP file (vcf format) into /BExplore/snpDATA/')

    primary_genome_exist = 0
    if 'primary_genome' in os.listdir('supportfile'):
        primary_genome_exist = 1
    return primary_genome_exist

def generate_primary_genome_file():
    genome_file = os.listdir('genomeDATA')[0]
    print('')
    print('Genome index has not been detected and is being created... It will take some time...')
    op = open('genomeDATA/{}'.format(genome_file),'r')
    linenum = 0
    numlist, titleseq, primaryseq = [], [], []
    for line in op.readlines():
        line = line.strip()
        linenum = linenum + 1
        numlist.append(linenum)
        if line[0] == '>':
            titleseq.append(linenum) 
        if line.find('Primary Assembly') != -1 and line.find('unlocalized')==-1 and line.find('unplaced') == -1:
            primaryseq.append(linenum) 

    op1 = open('supportfile/primary_genome','w')
    for i in primaryseq:
        op.seek(0)
        linenum = 0
        tt = titleseq.index(i)
        for line1 in op.readlines():
            linenum += 1
            if linenum >= i and linenum < titleseq[tt+1]:
                op1.write(line1)
            if linenum >= titleseq[tt+1]:
                break
        op1.close
    print('Genome index created successfully !')
    return

def find_primary_genome_index():
    grep_rs = subprocess.getoutput('grep -n ">" supportfile/primary_genome').split('\n')
    title_num_list = list(map(lambda x:int(x.split(':>')[0]), grep_rs))
    chr_list = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16',\
            '17', '18', '19', '20', '21', '22', 'X', 'Y']
    tidic = dict(zip(chr_list, title_num_list))
    return tidic

def get_location(args, tidic):
    linenum = 0
    middle, end, front = '', '', ''
    overline, overnum = int(), int()
    seq_double, seq_triple = '', ''
    ti_linenum = -2

    chro = str(args.chr)
    base = args.pos

    ti_linenum = tidic[chro]
    overline = base // 80
    overnum = base % 80
    op_genome.seek(0, 0)

    for line in op_genome.readlines():
        linenum += 1
        if linenum < ti_linenum:
            continue
        if linenum == ti_linenum + overline:
            front = line
        if linenum == ti_linenum + overline + 1:
            middle = line
        if linenum == ti_linenum + overline + 2:
            end = line
        if linenum > ti_linenum + overline + 10 and ti_linenum != -2:
            break
        seq_double = middle + end
        seq_triple = front + middle + end
    seq_double = seq_double.replace('\n', '')
    seq_triple = seq_triple.replace('\n', '')
    seq_triple_noupper = seq_triple  
    seq_triple = seq_triple.upper()
    targetbase = seq_triple[overnum - 1 + 80]
    return overline, overnum, seq_double, seq_triple, targetbase, seq_triple_noupper



def check_pam(granLen, pam):
    seq2 = seq_double[overnum - 1: overnum - 1 + granLen]
    seq2 = seq2.upper()
    re_cmd = {'NGG':'.GG', 'NNGRRT':'..G[AG][AG]T', 'NGCG':'.GCG', 'NNNRRT':'...[AG][AG]T',
              'NGAN':'.GA.', 'NGA':'.GA', 'NG':'.G'}
    if len(re.findall(re_cmd[pam], seq2[2:])) != 0:
        print('PAM test passed！{} exists in sequence downstream'.format(pam))
        pam_whether = 1
    else:
        print('PAM check failed! There is no {} downstream'.format(pam))
        quit()
        pam_whether = 0
    return pam_whether


def get_cangrna(args, pamLen, granLen, actWindow_range, seq_triple):
    cangrna = []  
    cangrna_startbase = {}
    re_cmd = {'NGG': '.GG', 'NNGRRT': '..G[AG][AG]T', 'NGCG': '.GCG', 'NNNRRT': '...[AG][AG]T',
              'NGAN': '.GA.', 'NGA': '.GA', 'NG': '.G'}
    if pam_whether == 1:
        pam_range = [granLen - pamLen + 1, granLen]
        actWindow_range
        baseLoc_in_seqTriple = overnum + 80
        minLocPam_distance = pam_range[0] - actWindow_range[1]
        maxLocPam_distance = pam_range[1] - actWindow_range[0]
        seq_triple = seq_triple[: overnum - 1 + 80] + 'C' + seq_triple[overnum - 1 + 80 + 1 :] 
        re_rs = regex.finditer(re_cmd[pam], seq_triple[baseLoc_in_seqTriple - 1 + minLocPam_distance:
                                                         baseLoc_in_seqTriple - 1 + maxLocPam_distance + 1], overlapped=True)
        pam_loc = []
        for t1 in re_rs:
            pam_loc.append(t1.span()) 
        for t2 in pam_loc:
            cangrna.append(seq_triple[(baseLoc_in_seqTriple - 1 + minLocPam_distance) + t2[1] - granLen:
                                      (baseLoc_in_seqTriple - 1 + minLocPam_distance) + t2[1]])
            cangrna_startbase[seq_triple[(baseLoc_in_seqTriple - 1 + minLocPam_distance) + t2[1] - granLen:
                                         (baseLoc_in_seqTriple - 1 + minLocPam_distance) + t2[1]]] = \
                (baseLoc_in_seqTriple - 1 + minLocPam_distance + t2[1] - granLen) + 1 + (overline - 1) * 80

        return cangrna, cangrna_startbase

def NNN_whether():
    nnn_dic = {}
    for l in cangrna:
        for i in range(len(l) - 6):
            a = 0
            for n in range(7):
                if l[i + n] == l[i]:
                    a += 1
                else:
                    a += 0
            if l not in nnn_dic:
                nnn_dic[l] = 0
            if a == 7:
                nnn_dic[l] = -1
            if a != 7 and nnn_dic[l] != -1:
                nnn_dic[l] = 0
    cangrna_result = []
    for l in cangrna:
        if nnn_dic[l] == 0:
            cangrna_result.append(l)
    return cangrna_result

def get_gc():
    gc_dic = {}
    for l in cangrna:
        gc = 0
        for ll in l:
            if ll == 'G' or ll == 'C':
                gc = gc + 1
        gc_dic[l] = gc / len(l)
    cangrna_result = []
    for l in cangrna:
        if gc_dic[l] <= 0.75 and gc_dic[l] >= 0.3:
            cangrna_result.append(l)
    return cangrna_result, gc_dic


def activity_window(args, actWindow_range):
    sameinwindow = {}
    type_editor_dic = {'c2t':'C', 'a2g':'A'}
    for i in cangrna:
        sameinwindow[i] = i[actWindow_range[0]-1: actWindow_range[1]-1+1].count(type_editor_dic[args.transtype])
    return (sameinwindow)

def gc_concentration():
    gccon_dic = {}
    for l in cangrna:
        gc = 0
        for ll in l:
            if ll == 'G' or ll == 'C':
                gc = gc + 1
        gccon_dic[l] = gc / len(l)
    gccon_dic1 = {}
    for nm in gccon_dic.keys():
        gccon_dic1[nm] = str(gccon_dic[nm] * 100)[:10] + '%'
    return gccon_dic, gccon_dic1

def get_frequency():
    frequency = {}
    for l in cangrna:
        frequency[l] = 0
        op_genome.seek(0, 0)
        a, b = '', ''
        for line in op_genome.readlines():
            line = line.strip('\n')
            line = line.upper()
            b = line
            if l in a + b:
                frequency[l] = frequency[l] + 1
            a = b
    return (frequency)

def count_SNP(args, cangrna_startbase, granLen):
    chro = str(args.chr)
    snprst_dic = {}
    for rna in cangrna:
        op = open('supportfile/snp_output_N2A.txt', 'w')
        cmd = ('vcftools --vcf snpDATA/%s --chr chr%s --from-bp %s --to-bp %s '
               % (os.listdir('snpDATA')[0], chro, cangrna_startbase[rna], cangrna_startbase[rna] + granLen - 1))
        output = subprocess.getoutput(cmd)
        op.write(output)
        op.close()
        result = subprocess.getoutput(" grep 'After filtering' supportfile/snp_output_N2A.txt\
        	| grep 'possible' ")
        result = result[result.find('kept') + 5: result.find('out') - 1]
        snprst_dic[rna] = int(result)

    return snprst_dic

def find_offtarget(pam, pamLen, granLen):
    offscore_dic = {}

    for wtseq in cangrna:
        offnum = 0
        cmd1 = ''
        op1 = open('supportfile/input_find_offtarget_N2A.txt', 'w')
        op1.write('genomeDATA/{}\n'.format(os.listdir('genomeDATA')[0]))
        op1.write('N' * (granLen - pamLen)  + pam + '\n')
        op1.write(str(wtseq) + ' ' + '3')
        op1.close()
        cmd1 = ('cas-offinder supportfile/input_find_offtarget_N2A.txt G \
        	supportfile/output_find_offtarget_N2A.txt')
        cmd1_out = subprocess.getoutput(cmd1)

        if 'No OpenCL devices found' in cmd1_out:
            logging.error('No OpenCL found, cas-offinder cannot run')
            sys.exit(1)

        op2 = open('supportfile/output_find_offtarget_N2A.txt', 'r')
        for line in op2.readlines():
            offnum += 1
            offseq = line.split('\t')[3]
            offseq = offseq.upper()
            cmd2 = ('python2.7 supportfile/CFD_Scoring/cfd-score-calculator.py '
                    '--wt %s --off %s' % (wtseq, offseq))
            CFDscore = subprocess.getoutput(cmd2)
            CFDscore = float(CFDscore.split(':')[1])
            if wtseq in offscore_dic:
                offscore_dic[wtseq] = offscore_dic[wtseq] + CFDscore
            else:
                offscore_dic[wtseq] = CFDscore
        if offnum != 0:
            offscore_dic[wtseq] = offscore_dic[wtseq] / offnum
        else:  
            offscore_dic[wtseq] = 0
    return offscore_dic


def RRA():
    def increase_sortdic(dic):
        dicvalues_list = list(dic.values())
        allsame_whether_list = list(filter(lambda x: dicvalues_list.count(x) == len(dicvalues_list),
                                           dicvalues_list))
        if len(allsame_whether_list) != 0 and len(dicvalues_list) != 1:
            sorted_list = ['None']
        else:
            sorted_tuplelist = sorted(dic.items(), key=lambda x: x[1])
            sorted_list = []
            for i in sorted_tuplelist:
                sorted_list.append(i[0])
        return sorted_list

    def decrease_sortdic(dic):
        sorted_tuplelist = sorted(dic.items(), key=lambda x: x[1], reverse=True)
        sorted_list = []
        for i in sorted_tuplelist:
            sorted_list.append(i[0])
        return sorted_list

    sameinwindow_list = increase_sortdic(sameinwindow_dic)
    gccon_list = increase_sortdic(gccon_dic)
    frequency_list = increase_sortdic(frequency_dic)
    snprst_list = increase_sortdic(snprst_dic)
    offscore_list = increase_sortdic(offscore_dic)

    whole_ranklist = [sameinwindow_list, gccon_list, frequency_list, snprst_list, offscore_list]
    whole_ranklist = list(map(lambda x: str(x).strip('[]'), whole_ranklist))
    whole_ranklist = list(filter(lambda x: x != "'None'", whole_ranklist))

    if len(whole_ranklist) == 0: 
        RRAresult_list, RRAresult_list1, RRAresult_list2 = [], [], []
        RRAresult_list1 = cangrna
        RRAresult_list2 = ['1'] * len(cangrna)
        RRAresult_list.append(','.join(RRAresult_list1))
        RRAresult_list.append(','.join(RRAresult_list2))
        return RRAresult_list

    cmd_str = 'c(' + '),c('.join(whole_ranklist) + ')\n'

    op = open('supportfile/ryuyan_N2A.R', 'w')
    cmd1 = 'library(RobustRankAggreg)\n'
    cmd2 = "glist <- list( %s )\n" % (cmd_str)
    cmd3 = 'aggregateRanks(glist = glist)\n'
    op.write(cmd1)
    op.write(cmd2)
    op.write(cmd3)
    op.close()
    RRAresult = subprocess.getoutput('Rscript supportfile/ryuyan_N2A.R')

    op_writeRRAresult = open('supportfile/RRAresult_N2A', 'w')
    op_writeRRAresult.write(RRAresult)
    op_writeRRAresult.close()
    op_readRRAresult = open('supportfile/RRAresult_N2A', 'r')

    RRAresult_list = []
    RRAresult_list1 = []
    RRAresult_list2 = []
    linecounter = 0
    title_line_num = int(subprocess.getoutput('grep -n "Name Score" "supportfile/RRAresult_N2A"').split(':')[0])
    for line in op_readRRAresult.readlines():
        linecounter += 1
        if linecounter > title_line_num:
            line = re.sub(r' +', ' ', line)
            line = line.replace('\n', '')
            line_list = line.split(' ')
            RRAresult_list1.append(line_list[1])
            RRAresult_list2.append(line_list[2])
    RRAresult_list.append(','.join(RRAresult_list1))
    RRAresult_list.append(','.join(RRAresult_list2))
    return RRAresult_list


def formatting_output(result_dic, RRAresult_list3):
    resultorder_list = RRAresult_list3[0].split(',')
    new_column = []
    for i in resultorder_list:
        new_column.append(str(result_dic[i]))
    column_result = ','.join(new_column)
    return column_result


def keep_decimal(column_str, decimal_number):
    if '.' in column_str:  
        column_list = column_str.split(',')
        new_column_list = []
        for i in column_list:
            i = round(float(i), decimal_number)
            new_column_list.append(str(i))
        new_column = ','.join(new_column_list)
    else:
        new_column = column_str
    return new_column

def show_pleiotropy_rs(args):
    df = pd.read_csv('supportfile/pleiotropy_rs.csv', sep='\t')
    if len(df[df['siteName'] == str(args.chr) + ':' + str(args.pos)]) == 0:
        print(str(args.chr) + ':' + str(args.pos), 'the location is not recorded in pleiotropy prediction results')
        return
    else:
        ple_info = df[df['siteName'] == str(args.chr) + ':' + str(args.pos)]
        chr_base_mutation = (ple_info.iloc[0, 1].split(':')[0] + ':' + ple_info.iloc[0, 1].split(':')[1] + ':'
                             + ple_info.iloc[0, 1].split(':')[2] + '>' + ple_info.iloc[0, 1].split(':')[3])

        num_betaAbove0_p005 = ple_info.iloc[0, 4]
        print(chr_base_mutation, '\t', 'This mutation at the site may significantly promote(beta>0 and p<005) ' \
                                        'the occurrence of the following {} diseases:'.format(num_betaAbove0_p005))
        if num_betaAbove0_p005 == 0:
            print('None', '\n')
        else:
            print('----', ple_info.iloc[0, 8].replace(' | ', '\n----'))
            print('')

        num_betaUnder0_p005 = ple_info.iloc[0, 6]
        print(chr_base_mutation, '\t', 'This mutation at the site may significantly suppress(beta<0 and p<005) ' \
                                 'the occurrence of the following {} diseases:'.format(num_betaUnder0_p005))
        if num_betaUnder0_p005 == 0:
            print('None', '\n')
        else:
            print('----', ple_info.iloc[0, 9].replace(' | ', '\n---- '))
            print('')
        return

if __name__ == '__main__':
    args = InputParser()
    editor_info, cbe_editors, abe_editors = editors_store(args)
    pam, pamLen, granLen, actWindow_range, ori_base = editor_info

    check_softwares()
    primary_genome_exist = ckeck_File_Dir()
    if primary_genome_exist == 0:
        generate_primary_genome_file()
    else:
        print('')
        print('Genome index has been created!')
    tidic = find_primary_genome_index()

    with open('supportfile/primary_genome', 'r') as op_genome, \
            open('output_file_{}:{}.tsv'.format(str(args.chr), str(args.pos)), 'w') as op_output:
        pam_whether = -100

        if (args.transtype == 'c2t' and args.editor in abe_editors) or \
                (args.transtype == 'a2g' and args.editor in cbe_editors):
            print('input error! Conversion type(--transtype) does not match the selected base editor(--editor)')
            quit()

        overline, overnum, seq_double, seq_triple, targetbase, seq_triple_noupper = get_location(args, tidic)
        print('')
        print('Base at the position in reference genome is:', targetbase)

        table_Info = [str(args.chr), str(args.pos), targetbase, args.transtype, args.editor]
        op_output.write('\t'.join(['chr', 'position', 'target_in_ref', 'transtype', 'editor', 'candidate _gRNA_num',
                                   'candidate _gRNA', 'RRAscore', 'identicalBase_num_in_actWindow', 'GC%', 'repeated_gRNA_num_in_genome',
                                   'SNP_num', 'offtarget_score']) + '\n')

        pam_whether = check_pam(granLen, pam)
        if pam_whether == 1 :  
            file_write_whether = 0  

            (cangrna, cangrna_startbase) = get_cangrna(args, pamLen, granLen, actWindow_range, seq_triple)
            print('')
            print('* * * * * * * * 1. CANDIDATE GRNAS ARE AS FOLLOWS * * * * * * * * * * * * * *')
            for grna in cangrna_startbase:
                print(grna, '\t','gRNA start site：', cangrna_startbase[grna], '\t',
                      'target base position on gRNA:', args.pos - cangrna_startbase[grna]+1)

            if len(cangrna) == 0:  
                print('Activity window filter failed')
                op_output.write('\t'.join(table_Info + ['0'] + ['activity_window_error']*7) + '\n')
                file_write_whether = 1
                quit()

            cangrna = NNN_whether()

            if len(cangrna) == 0 and file_write_whether != 1:  
                print('Continuous same bases filter failed')

                op_output.write('\t'.join(table_Info +['0'] + ['continuous_identical_base_error'] * 7) + '\n')
                file_write_whether = 1
                quit()

            (cangrna, gc_dic) = get_gc()

            nopass_gc_list = list(gc_dic.values())
            nopass_gc_list = list(map(lambda x: round(x, 6), nopass_gc_list))
            nopass_gc_list = list(map(str, nopass_gc_list))
            nopass_gc_str = ','.join(nopass_gc_list)
            print('gRNAs that pass all filters：' + ','.join(cangrna))

            if len(cangrna) == 0 and file_write_whether != 1:  
                print('GC% filter failed')
                op_output.write('\t'.join(table_Info + ['0'] + ['GC_ratio_error'] * 3 + [nopass_gc_str] + ['GC_ratio_error'] * 3) + '\n')
                file_write_whether = 1
                quit()

            print('')
            print('* * * * * * * * 2. SCORES OF EACH CRITERIA ARE AS FOLLOWS * * * * * * * * * * ')

            if len(cangrna) != 0:  

                sameinwindow_dic = activity_window(args, actWindow_range)
                print('Number of Bases Identical to Target Site in Activity Window: ', str(sameinwindow_dic).strip('{}').replace("'", ''))

                gccon_dic, gccon_dic1 = gc_concentration()
                print('GC% in the gRNA Sequence:           ', str(gccon_dic1).strip('{}').replace("'", ''))

                frequency_dic = get_frequency()
                print('Number of Tepeated gRNA Sequences in Whole Genome:           ', str(frequency_dic).strip('{}').replace("'", ''))

                snprst_dic = count_SNP(args, cangrna_startbase, granLen)
                print('Number of SNPs in gRNA Sequence:    ', str(snprst_dic).strip('{}').replace("'", ''))

                offscore_dic = find_offtarget(pam, pamLen, granLen)
                print('Potential off-target Score:         ', str(offscore_dic).strip('{}').replace("'", ''))

                RRAresult_list = RRA()
                print('Final Score(RRA):                   ', str(RRAresult_list).strip('[]').replace("'", ''))

                column9 = formatting_output(sameinwindow_dic, RRAresult_list)
                column10 = formatting_output(gccon_dic, RRAresult_list)
                column10 = keep_decimal(column10, 6)
                column11 = formatting_output(frequency_dic, RRAresult_list)
                column12 = formatting_output(snprst_dic, RRAresult_list)
                column13 = formatting_output(offscore_dic, RRAresult_list)
                column13 = keep_decimal(column13, 6)

                op_output.write('\t'.join(table_Info + [str(len(cangrna)), RRAresult_list[0], RRAresult_list[1], column9,
                                                        column10, column11, column12, column13]) + '\n')

                print('')
                print('* * * * * * * * 3. BEST GRNA * * * * * * * * * * * * * * * * * * * * * * * * *')
                print('Best gRNA:', RRAresult_list[0].split(',')[0])

        else:  
            op_output.write('\t'.join(table_Info + ['0'] + ['no_PAM'] * 7) + '\n')
            quit()

    print('')
    print('* * * * * * * * 4. PLEIOTROPY PREDICTION OF THE LOCATION * * * * * * * * * * *')
    show_pleiotropy_rs(args)
    print('')
    print('* * * * * * * * ANALYSIS COMPLETE！ * * * * * * * * ')
