# -*- coding:utf-8 -*-
import os
import sys
import configparser
from datetime import datetime
import argparse
from itertools import islice
import time
parser = argparse.ArgumentParser(description='Put a description of your script here')
parser.add_argument('-s', '--species', type=str, required=True, help='species')
parser.add_argument('-t', '--thres', type=str, required=None, help='choose threshold by yourself')
args = parser.parse_args()
usage = '''Program:  threshold_analysis
Version:  0.0.2
Usage:    python3 threshold_analysis2.py -s species -t 4.3  # fang's threshold
or: python3 threshold_analysis2.py -s species
Contact:  YSQ '''

TimeFormat = '%m/%d/%Y %H:%M:%S'
rootDir = os.path.split(os.path.realpath(__file__))[0]
configFilePath = os.path.join(rootDir, '%s_config.ini' % args.species)
# usage = '''Program:  threshold_analysis
# Version:  0.0.1
# Usage:    python3 threshold_analysis.py
# Contact:  YSQ '''
#
# TimeFormat = '%m/%d/%Y %H:%M:%S'
# rootDir = os.path.split(os.path.realpath(__file__))[0]
# configFilePath = os.path.join(rootDir, '%s_config.ini' % species)

def myprint(string):
    print('[' + datetime.now().strftime(TimeFormat) + ']' + string)
    return


# 峰被enrichment和pvalue过滤
def readpeakfile():
    peak_filtered = os.path.join(out_dir, "%s_peakfilter.txt" % species)
    if not os.path.exists(peak_filtered):
        #cmd = "awk '{if($7>=1&&$8>2&&$9>2)print$1 \"\t\" $2 \"\t\" $3 \"\t\" $7 \"\t\" $8 \"\t\" $9}' %s > %s;echo $?" % (
        #peak_file, peak_filtered)
        cmd = "awk '{if($7>=1&&$8>2)print$1 \"\t\" $2 \"\t\" $3 \"\t\" $7 \"\t\" $8 \"\t\" $9}' %s > %s;echo $?" % (
        peak_file, peak_filtered)
        p = os.popen(cmd)
        a = p.read().split()
        if a[0] == "0" and os.path.exists(peak_filtered) and os.stat(peak_filtered).st_size > 0:
            return peak_filtered
        else:
            myprint("Error in command or input peakfile!")
    else:
        return peak_filtered

# 峰未被enrichment和pvalue过滤
def readpeakfile_unfiltered():
    peak = os.path.join(out_dir, "%s_peak_unfiltered.txt" % species)
    if not os.path.exists(peak):
        cmd = "awk '{print$1 \"\t\" $2 \"\t\" $3 \"\t\" $7 \"\t\" $8 \"\t\" $9}' %s > %s;echo $?" % (
        peak_file, peak)
        p = os.popen(cmd)
        a = p.read().split()
        if a[0] == "0" and os.path.exists(peak) and os.stat(peak).st_size > 0:
            return peak
        else:
            myprint("Error in command or input peakfile!")
    else:
        return peak


# 将6ma点从原始文件中grep出来
def preprocess_modification_file():
    m6A_modifications = os.path.join(out_dir, "%s_m6A_modifications.gff" % species)
    if not os.path.exists(m6A_modifications):
        cmd1 = "grep -v '#' %s |grep m6A >%s"%(modification_file,m6A_modifications)    # ysq mod
        os.system(cmd1)

    tmpfile = os.path.join(out_dir, "%s_tmp.gff" % species)
    if os.path.exists(tmpfile):
        return tmpfile
    else:
        cmd2 = "cut -f 1,3,4,5,6,9 %s|sed 's/coverage=//g'|sed 's/;context=/\t/g'|sed 's/;IPDRatio=/\t/g'|sed 's/;frac=/\t/g'|sed 's/;fracLow=/\tfracflow=/g'|sed 's/;/\t###\t/g'|cut -f 1,2,3,4,5,6,7,8,9>%s"%(m6A_modifications,tmpfile)        # ysq mod
        os.system(cmd2)
        return tmpfile


# 将6ma点进行score和coverage过滤
def readmodificationfile():
    tmpfile = preprocess_modification_file()
    modifications_filtered = os.path.join(out_dir, "%s_modfilter.txt" % species)
    if not os.path.exists(modifications_filtered):
        if os.path.exists(tmpfile):
            cmd = "awk '{if($5>=30&&$6>=50)print$0}' %s > %s;echo $?"%(tmpfile, modifications_filtered)   # ysq mod
        else:
            myprint("Error in cmd or tmpfile!")
        p = os.popen(cmd)
        a = p.read().split()
        if a[0] == "0" and os.path.exists(modifications_filtered) and os.stat(modifications_filtered).st_size > 0:
            return modifications_filtered
        else:
            myprint("Error in command or input modificationfile!")
    else:
        return modifications_filtered


def generate_outfile(outfile, content):
    with open(outfile, 'a') as fo:
        fo.write(content)


def calculate_threshold():
    datasetfile = os.path.join(out_dir, "%s_filtered_dataset.txt" % species)
    if os.path.exists(datasetfile):
        cmd = "Rscript threshold.R %s %s %s;echo $?" % (species, out_dir, datasetfile)
        p = os.popen(cmd)
        a = p.read().split()
        if a[-1] == "0":
            return a[1]
        else:
            myprint("ERROR:%s"%a)
    else:
        myprint("ERROR:there is not dataset.txt, please run construct_dataset function!")
        sys.exit()

# peakfilterfile = readpeakfile()
# modfilterfile = readmodificationfile()
# peakunfilterfile = peak_file
# modunfilterfile = preprocess_modification_file()
# 构建overlap数据集
def construct_overlap_dataset(peakfile,modfile,filter_state):

    info_list = []
    if os.path.exists(peakfile) and os.stat(peakfile).st_size > 0:
        with open(peakfile, 'r') as fp:
            while True:
                lines = list(islice(fp, 10000))
                for line in lines:
                    line = line.strip().split()
                    p_chr = line[0]
                    p_s = line[1]
                    p_e = line[2]
                    p_enrich = line[3]
                    p_pvalue = line[4]
                    p_qvalue = line[5]
                    p_info = [p_chr, p_s, p_e, p_enrich, p_pvalue, p_qvalue]
                    info_list.append(p_info)
                if not lines:
                    break


    if os.path.exists(modfile) and os.stat(modfile).st_size > 0:
        with open(modfile, 'r') as fg:
            n = 0
            dataset = os.path.join(out_dir, "%s_%s_dataset.txt" % (species, filter_state))
            while True:
                ls = list(islice(fg, 10000))
                for l in ls:
                    li = l.strip().split()
                    chr = li[0]
                    m6a_locus = int(li[2])
                    score = li[4]
                    coverage = li[5]
                    context = li[6]
                    ipdratio = li[7]
                    fraction = li[8]
                    for ran in info_list:
                        chr_name = ran[0]
                        start = int(ran[1])
                        end = int(ran[2])
                        enrich = ran[3]
                        if m6a_locus >= start and m6a_locus <= end and chr == chr_name and fraction != "":
                            n += 1
                            t = str(chr_name) + '\t' + str(m6a_locus) + '\t' + str(score) + '\t' + str(coverage) + '\t' + str(context) + '\t' + str(fraction) + '\t' + str(ipdratio) + '\t' + str(enrich) + '\n'
                            generate_outfile(dataset, t)

                if not ls:
                    break
            if filter_state == 'filtered':
                threshold = calculate_threshold()
                print("the quantity of 6ma events in peaks regions: %d" % n)
                print("the threshold of IPDratio: %s" % threshold)
                return threshold
            elif filter_state == 'unfiltered':
                pass



def judge_motif(inputfile,outfile,context_index,motif, motif_loc):
    i = int(motif_loc) - 1
    length = len(motif)
    s = 20 - i
    e = 20 + length - i
    with open(inputfile, 'r') as f:
        m = 0
        for ll in f.readlines():
            li = ll.strip().split()
            context1 = li[context_index][s:e]
            index = -1
            flag = 0
            while index<len(context1)-1:
                index += 1
                if motif[index] in ['A','G','T','C'] and context1[index] == motif[index] or motif[index] not in ['A','G','T','C']:
                    flag = 1    
                    continue
                else:
                    flag = 0
                    break

            if flag == 1:
                m += 1
                c = ll + '\t' + context1        # ysq mod
                generate_outfile(outfile, c)
        return m
       

def find_motifs(motif, motif_loc):
    tmpfile = os.path.join(out_dir, "%s_tmp.gff" % species)
    modfilterfile = os.path.join(out_dir, "%s_modfilter.txt" % species)   # 经score和coverage过滤后的点,*_modfilter.txt
    datasetunfilterd = os.path.join(out_dir, "%s_unfiltered_dataset.txt" % species)
    datasetfilterd = os.path.join(out_dir, "%s_filtered_dataset.txt" % species) # 过滤后的且落在peak里的点
    if os.path.exists(tmpfile) and os.stat(tmpfile).st_size > 0:
        motifsall_unfilterd = os.path.join(out_dir, "%s_%s_motifsall_unfilterd.txt" % (species,motif))
        check_file_exist_rm(motifsall_unfilterd)
        cmd = "touch %s" % motifsall_unfilterd
        os.system(cmd)
        m1 = judge_motif(tmpfile,motifsall_unfilterd,6,motif, motif_loc)
        myprint("the quantity of %s motifsall unfilterd:%d" % (motif, m1))
    else:
        myprint("Don't exist file we need,please modify the value of construct_dataset function!")

    if os.path.exists(modfilterfile) and os.stat(modfilterfile).st_size > 0:
        motifsall_filterd = os.path.join(out_dir, "%s_%s_motifsall_filterd.txt" % (species,motif))
        check_file_exist_rm(motifsall_filterd)
        cmd = "touch %s" % motifsall_filterd
        os.system(cmd)
        m2 = judge_motif(modfilterfile,motifsall_filterd,6,motif, motif_loc)
        myprint("the quantity of %s motifsall filterd:%d" % (motif, m2))
    else:
        myprint("Don't exist file we need,please modify the value of construct_dataset function!")

    if os.path.exists(datasetunfilterd) and os.stat(datasetunfilterd).st_size > 0:
        motifsinpeak_unfilterd = os.path.join(out_dir, "%s_%s_motifsinpeak_unfilterd.txt" % (species,motif))
        check_file_exist_rm(motifsinpeak_unfilterd)
        cmd = "touch %s" % motifsinpeak_unfilterd
        os.system(cmd)
        n1 = judge_motif(datasetunfilterd,motifsinpeak_unfilterd,4,motif, motif_loc)
        myprint("the quantity of %s motif in unfilterd peaks:%d" % (motif, n1))
    else:
        myprint("Don't exist file we need,please modify the value of construct_dataset function!")

    if os.path.exists(datasetfilterd) and os.stat(datasetfilterd).st_size > 0:
        motifsinpeak_filterd = os.path.join(out_dir, "%s_%s_motifsinpeak_filterd.txt" % (species, motif))
        check_file_exist_rm(motifsinpeak_filterd)
        cmd = "touch %s" % motifsinpeak_filterd
        os.system(cmd)
        n2 = judge_motif(datasetfilterd, motifsinpeak_filterd, 4, motif, motif_loc)
        myprint("the quantity of %s motif in filterd peaks:%d" % (motif, n2))
    else:
        myprint("Don't exist file we need,please modify the value of construct_dataset function!")


def analyse(threshold,motif):
    if threshold:
        motifsall_unfilterd = os.path.join(out_dir, "%s_%s_motifsall_unfilterd.txt" % (species, motif))
        motifsall_filterd = os.path.join(out_dir, "%s_%s_motifsall_filterd.txt" % (species, motif))
        motifsinpeak_unfilterd = os.path.join(out_dir, "%s_%s_motifsinpeak_unfilterd.txt" % (species, motif))
        motifsinpeak_filterd = os.path.join(out_dir, "%s_%s_motifsinpeak_filterd.txt" % (species, motif))

        
        #pacbio_unfiltered
        cmd2 = "wc -l %s" % motifsall_unfilterd  #所有的motifs
        p2 = os.popen(cmd2)
        n2 = p2.read().split()[0]

        # pacbio_filtered
        cmd6 = "wc -l %s" % motifsall_filterd  # 经sc过滤后的motifs
        p6 = os.popen(cmd6)
        n6 = p6.read().split()[0]

        # pacbio_unfiltered+MeDIP_unfiltered
        cmd7 = "wc -l %s" % motifsinpeak_unfilterd  # 经peak过滤后的motifs
        p7 = os.popen(cmd7)
        n7 = p7.read().split()[0]

        #pacbio_filtered+MeDIP_filtered
        cmd3 = "wc -l %s" % motifsinpeak_filterd  # 经scep过滤后的motifs
        p3 = os.popen(cmd3)
        n3 = p3.read().split()[0]

        # pacbio_unfiltered+threshold
        cmd8 = "awk '$8>=%f{print}' %s|wc -l" % (float(threshold), motifsall_unfilterd)  # 经阈值过滤后所有的motifs
        p8 = os.popen(cmd8)
        n8 = p8.read().split()[0]
        
        #pacbio_filtered+threshold
        cmd4 = "awk '$8>=%f{print}' %s|wc -l" % (float(threshold), motifsall_filterd)  # 经sc阈值过滤后所有的motifs
        p4 = os.popen(cmd4)
        n4 = p4.read().split()[0]

        # pacbio_unfiltered+MeDIP_unfiltered+threshold
        cmd9 = "awk '$7>=%f{print}' %s|wc -l" % (float(threshold), motifsinpeak_unfilterd)  # 经peak阈值过滤后的motifs
        p9 = os.popen(cmd9)
        n9 = p9.read().split()[0]
        
        #pacbio_filtered+MeDIP_filtered+threshold
        cmd5 = "awk '$7>=%f{print}' %s|wc -l" % (float(threshold), motifsinpeak_filterd)  # 经scep阈值过滤后的motifs
        p5 = os.popen(cmd5)
        n5 = p5.read().split()[0]

        content = "motif:%s" % motif + '\n' + '(pacbio_unfiltered)%s motifs: %s' % (motif, n2) + '\n' + '(pacbio_filtered)%s motifs: %s' % (motif, n6) + '\n' + '(pacbio_unfiltered+MeDIP_unfiltered)%s motifs: %s' % (motif, n7) + '\n' + '(pacbio_filtered+MeDIP_filtered)%s motifs: %s' % (motif, n3) + '\n'+ '(pacbio_unfiltered+threshold)%s motifs: %s' % (motif, n8) + '\n'+'(pacbio_filtered+threshold)%s motifs: %s' % (motif, n4) + '\n'+ '(pacbio_unfiltered+MeDIP_unfiltered+threshold)%s motifs: %s' % (motif, n9) + '\n'+ '(pacbio_filtered+MeDIP_filtered+threshold)%s motifs: %s' % (motif, n5) + '\n'
        return content
    else:
        print("ERROR:there is not threshold")

# 计算未处理前所有6ma点的数量
def caculate_total_modifications():
    tmpfile = preprocess_modification_file()
    cmd = "wc -l %s" % tmpfile
    p = os.popen(cmd)
    n = p.read().split()[0]
    return n

# 计算经score、coverage过滤后的6ma点的数量
def caculate_sc_modifications():
    sc_file = readmodificationfile()
    cmd = "wc -l %s" % sc_file
    p = os.popen(cmd)
    n = p.read().split()[0]
    return n

# 计算所有位点只被阈值过滤的数量
def caculate_filtered_modifications(threshold):
    tmpfile = preprocess_modification_file()
    mod_thres_filter = os.path.join(out_dir, "mod_thres_filter.txt")
    check_file_exist_rm(mod_thres_filter)
    if os.path.exists(tmpfile):
        cmd = "awk '{if($8>=%f)print$0}' %s >%s;echo $?" % (
             float(threshold),tmpfile, mod_thres_filter)  #ysq mod
    else:
        myprint("Error in cmd or tmpfile!")

    p = os.popen(cmd)
    a = p.read().split()[0]
    if a == "0":
        cmd1 = "wc -l %s" % mod_thres_filter  # 阈值过滤后所有的6ma
        p = os.popen(cmd1)
        n = p.read().split()[0]
        return n
    else:
        myprint("Error in command or input modificationfile!")

# 经score、coverage过滤后的点落在经enrichment、pvalue过滤后的峰中的数量
def caculate_peak_modifications_filtered():
    datasetfile = os.path.join(out_dir, "%s_filtered_dataset.txt" % species)
    if os.path.exists(datasetfile):
        cmd = "wc -l %s" % datasetfile
        p = os.popen(cmd)
        n = p.read().split()[0]
        return n
    else:
        print("datasetfile does not exists!")

# 未经过滤后的点落在未经过滤的峰中的数量
def caculate_peak_modifications_unfiltered():
    datasetfile = os.path.join(out_dir, "%s_unfiltered_dataset.txt" % species)
    if os.path.exists(datasetfile):
        cmd = "wc -l %s" % datasetfile
        p = os.popen(cmd)
        n = p.read().split()[0]
        return n
    else:
        print("datasetfile does not exists!")

# 经score、coverage、threshold过滤后的点落在经enrichment、pvalue过滤后的峰中的数量
def caculate_peak_thres_modifications_filtered(threshold):
    datasetfile = os.path.join(out_dir, "%s_filtered_dataset.txt" % species)
    if os.path.exists(datasetfile):
        cmd = "awk '$7>=%f{print}' %s|wc -l" % (float(threshold), datasetfile)
        p = os.popen(cmd)
        n = p.read().split()[0]
        return n
    else:
        print("datasetfile does not exists!")

# 经threshold过滤后的点落在未经过滤的峰中的数量
def caculate_peak_thres_modifications_unfiltered(threshold):
    datasetfile = os.path.join(out_dir, "%s_unfiltered_dataset.txt" % species)
    if os.path.exists(datasetfile):
        cmd = "awk '$7>=%f{print}' %s|wc -l" % (float(threshold), datasetfile)
        p = os.popen(cmd)
        n = p.read().split()[0]
        return n
    else:
        print("datasetfile does not exists!")

# 计算未过滤的峰的数量
def caculate_peak_unfiltered():
    if os.path.exists(peak_file):
        cmd = "wc -l %s" % peak_file
        p = os.popen(cmd)
        n = p.read().split()[0]
        return n
    else:
        print("peak_file does not exists!")

# 计算过滤后的峰的数量
def caculate_peak_filtered():
    peak_filtered = os.path.join(out_dir, "%s_peakfilter.txt" % species)
    if os.path.exists(peak_filtered):
        cmd = "wc -l %s" % peak_filtered
        p = os.popen(cmd)
        n = p.read().split()[0]
        return n
    else:
        print("peak_filtered does not exists!")

def check_file_exist_rm(file):
    if os.path.exists(file):
        os.system("rm -rf %s" % file)
    else:
        pass

def only_check_file_exist(file):
    if os.path.exists(file) and os.stat(file).st_size > 0:
        pass
    else:
        myprint("ERROR: %s is not exist or size is zero!" % file)
        sys.exit()

def only_check_dir_exist(dir):
    if os.path.exists(dir):
        myprint("%s has been exist!" % dir)
    else:
        os.mkdir(dir)

def get_config_values():
    global species
    global peak_file
    global modification_file
    global out_dir
    global con_func
    global find_func
    global analyse_func
    global motifs_list
    global motif_loc_list
    global thres
    config = configparser.ConfigParser()
    config.read(configFilePath)
    species = config.get('species', 'species')
    peak_file = config.get('input', 'narrowpeak')
    modification_file = config.get('input', '6ma_modification')
    out_dir = config.get('out', 'outdir')
    con_func = config.get('function', 'construct_dataset')
    find_func = config.get('function', 'find_motifs')
    analyse_func = config.get('function', 'analyse')
    motifs = config.get('motifs', 'motifs')
    motif_loc = config.get('motifs', 'motif_loc')

    if motifs == "" and motif_loc == "":
        motifs_list = []
        motif_loc_list = []
    else:
        motifs_list = config.get('motifs', 'motifs').split(",")
        motif_loc_list = config.get('motifs', 'motif_loc').split(",")

    if species == "":
        myprint('ERROR:species is required')
        sys.exit()

    elif peak_file == "":
        myprint('ERROR:narrowpeak is required')
        sys.exit()

    elif modification_file == "":
        myprint('ERROR:6ma_modification is required!')
        sys.exit()

    elif out_dir == "":
        myprint('ERROR:out_dir is required!')
        sys.exit()

    elif con_func == "":
        myprint('ERROR:con_function should be TRUE or FALSE!')
        sys.exit()

    elif find_func == "":
        myprint('ERROR:find_function should be TRUE or FALSE!')
        sys.exit()

    elif analyse_func == "":
        myprint('ERROR:verify_function should be TRUE or FALSE!')
        sys.exit()

    elif len(motifs_list) != len(motif_loc_list):
        myprint('ERROR:motif_loc should be correspond to motif')
        sys.exit()


def main():
    if len(sys.argv) > 4:
        myprint(usage)
        sys.exit()

    myprint('program started')
    # 1.读取参数文件
    myprint('reading config file')
    get_config_values()
    # 2.检查是否存在out_dir,不存在就创建
    only_check_dir_exist(out_dir)
    # 3.检查是否存在配置文件,峰文件，原始位点文件,不存在就退出系统
    only_check_file_exist(configFilePath)
    only_check_file_exist(peak_file)
    only_check_file_exist(modification_file)
    # 4.如果构建数据集函数为真，初次产生dataset.txt文件并预测阈值
    content = '[' + datetime.now().strftime(TimeFormat) + ']'
    if con_func == "TRUE":
        # 构建点和峰的重合数据集（点和峰经过score、coverage、pvalue、enrichment过滤），直接产生阈值
        myprint('running tasks: construct_dataset')
        peakfilterfile = readpeakfile()
        modfilterfile = readmodificationfile()
        peakunfilterfile = readpeakfile_unfiltered()
        modunfilterfile = preprocess_modification_file()
        datasetfile = os.path.join(out_dir, "%s_unfiltered_dataset.txt" % species)
        datasetfile1 = os.path.join(out_dir, "%s_filtered_dataset.txt" % species)
        if os.path.exists(datasetfile) and os.stat(datasetfile).st_size > 0:
            s1 = time.time()
            if os.path.exists(datasetfile1) and os.stat(datasetfile1).st_size > 0:
                threshold = calculate_threshold()
            else:
                threshold = construct_overlap_dataset(peakfilterfile, modfilterfile, 'filtered')
            e1 = time.time()
            print(e1 - s1)
        else:
            s=time.time()
            construct_overlap_dataset(peakunfilterfile,modunfilterfile,'unfiltered')
            e=time.time()
            print(e-s)
            s1 = time.time()
            if os.path.exists(datasetfile1) and os.stat(datasetfile1).st_size > 0:
                threshold = calculate_threshold()
            else:
                threshold = construct_overlap_dataset(peakfilterfile, modfilterfile, 'filtered')
            e1 = time.time()
            print(e1 - s1)

        # 未处理前所有6ma点的数量
        n = caculate_total_modifications()
        n4 = caculate_sc_modifications()
        n1 = caculate_peak_modifications_filtered()
        n5 = caculate_peak_modifications_unfiltered()
        n7 = caculate_peak_unfiltered()
        n8 = caculate_peak_filtered()
        if threshold is None:
            content += '%s:%s(threshold)' % (species,threshold) + '\n' + 'Total 6ma sites in genomes before filtered: %s' % n + '\n'+ 'The 6ma sites in genomes after score and coverage filtered: %s' % n4 + '\n' + 'The unfiltered 6ma sites in unfiltered peaks: %s' % n5 +'\n'+'The filtered 6ma sites in filtered peaks: %s' % n1 + '\n'+'The unfiltered peaks: %s' % n7 + '\n'+'The filtered peaks by enrichment and -log(pvalue): %s' % n8 + '\n'
        else:
            n2 = caculate_filtered_modifications(threshold)
            n3 = caculate_peak_thres_modifications_filtered(threshold)
            n6 = caculate_peak_thres_modifications_unfiltered(threshold)
            content += '%s:%s(threshold)' % (species,threshold) + '\n' + 'Total 6ma sites in genomes before filtered: %s' % n + '\n' + 'The 6ma sites in genomes after score and coverage filtered: %s' % n4 + '\n' + 'The unfiltered 6ma sites in unfiltered peaks: %s' % n5 + '\n' + 'The filtered 6ma sites in filtered peaks: %s' % n1 + '\n' + 'The unfiltered peaks: %s' % n7 + '\n' + 'The filtered peaks by enrichment and -log(pvalue): %s' % n8 + '\n'
            content += 'Total 6ma sites in genomes after threshold filtered : %s' % n2 + '\n' + 'The unfiltered 6ma sites in unfiltered peaks after threshold filtered : %s' % n3 + '\n' + 'The filtered 6ma sites in filtered peaks after threshold filtered : %s' % n6 + '\n'
            if motifs_list != []:
                for motif in motifs_list:
                    motif_index = motifs_list.index(motif)
                    motif_loc = motif_loc_list[motif_index]
                    if int(motif_loc) <= 0 or int(motif_loc) > len(motif):
                        myprint('ERROR:motif_loc should over zero and below the length of motif!')
                        sys.exit()
                    else:
                        if find_func == "TRUE" and len(motifs_list) != 0 and len(motif_loc_list) != 0:
                            myprint('running tasks: find_motifs')
                            find_motifs(motif,motif_loc)
                        if analyse_func == "TRUE":
                            myprint('running tasks: analyse')
                            content += analyse(threshold,motif)
        infofile = os.path.join(out_dir, '%s_%s_info.txt' % (species,threshold))
        generate_outfile(infofile, content)
    # 4.如果构建数据集函数为假
    else:
        # 如果传了阈值参数，就用传入的阈值
        if args.thres:
            threshold = args.thres
        # 不然就是已经存在dataset.txt文件，可以预测阈值
        else:
            threshold = calculate_threshold()
        n = caculate_total_modifications()
        n4 = caculate_sc_modifications()
        n1 = caculate_peak_modifications_filtered()
        n5 = caculate_peak_modifications_unfiltered()
        n7 = caculate_peak_unfiltered()
        n8 = caculate_peak_filtered()
        if threshold is None:
            content += '%s:%s(threshold)' % (species,threshold) + '\n' + 'Total 6ma sites in genomes before filtered: %s' % n + '\n' + 'The 6ma sites in genomes after score and coverage filtered: %s' % n4 + '\n' + 'The unfiltered 6ma sites in unfiltered peaks: %s' % n5 + '\n' + 'The filtered 6ma sites in filtered peaks: %s' % n1 + '\n' + 'The unfiltered peaks: %s' % n7 + '\n' + 'The filtered peaks by enrichment and -log(pvalue): %s' % n8 + '\n'
        else:
            n2 = caculate_filtered_modifications(threshold)
            n3 = caculate_peak_thres_modifications_filtered(threshold)
            n6 = caculate_peak_thres_modifications_unfiltered(threshold)
            content += '%s:%s(threshold)' % (species,threshold) + '\n' + 'Total 6ma sites in genomes before filtered: %s' % n + '\n' + 'The 6ma sites in genomes after score and coverage filtered: %s' % n4 + '\n' + 'The unfiltered 6ma sites in unfiltered peaks: %s' % n5 + '\n' + 'The filtered 6ma sites in filtered peaks: %s' % n1 + '\n' + 'The unfiltered peaks: %s' % n7 + '\n' + 'The filtered peaks by enrichment and -log(pvalue): %s' % n8 + '\n'
            content += 'Total 6ma sites in genomes after threshold filtered : %s' % n2 + '\n' + 'The unfiltered 6ma sites in unfiltered peaks after threshold filtered : %s' % n3 + '\n' + 'The filtered 6ma sites in filtered peaks after threshold filtered : %s' % n6 + '\n'
        if motifs_list != []:
            for motif in motifs_list:
                motif_index = motifs_list.index(motif)
                motif_loc = motif_loc_list[motif_index]
                if int(motif_loc) <= 0 or int(motif_loc) > len(motif):
                    myprint('ERROR:motif_loc should over zero and below the length of motif!')
                    sys.exit()
                else:
                    if find_func == "TRUE" and len(motifs_list) != 0 and len(motif_loc_list) != 0:
                        myprint('running tasks: find_motifs')
                        find_motifs(motif,motif_loc)
                    if analyse_func == "TRUE":
                        myprint('running tasks: analyse')
                        content += analyse(threshold,motif)
        infofile = os.path.join(out_dir, '%s_%s_info.txt' % (species,threshold))
        generate_outfile(infofile, content)

if __name__ == '__main__':
    main()
