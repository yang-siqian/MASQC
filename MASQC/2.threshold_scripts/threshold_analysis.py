# -*- coding:utf-8 -*-
import os
import sys
import configparser
from datetime import datetime
import argparse
parser = argparse.ArgumentParser(description='Put a description of your script here')
parser.add_argument('-s', '--species', type=str, required=True, help='species')
args = parser.parse_args()

usage = '''Program:  threshold_analysis 
Version:  0.0.1 
Usage:    python3 threshold_analysis.py
Contact:  YSQ '''

TimeFormat = '%m/%d/%Y %H:%M:%S'
rootDir = os.path.split(os.path.realpath(__file__))[0]
configFilePath = os.path.join(rootDir, '%s_config.ini' % args.species)
usage = '''Program:  threshold_analysis 
Version:  0.0.1 
Usage:    python3 threshold_analysis.py -s species
Contact:  YSQ '''

TimeFormat = '%m/%d/%Y %H:%M:%S'
rootDir = os.path.split(os.path.realpath(__file__))[0]
configFilePath = os.path.join(rootDir, '%s_config.ini' % args.species)


def myprint(string):
    print('[' + datetime.now().strftime(TimeFormat) + ']' + string)
    return


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


def readpeakfile():
    peak_filtered = os.path.join(out_dir, "%s_peakfilter.txt" % species)
    check_file_exist_rm(peak_filtered)
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


def preprocess_modification_file():
    m6A_modifications = os.path.join(out_dir, "%s_m6A_modifications.gff" % species)
    check_file_exist_rm(m6A_modifications)
    cmd = "grep m6A %s >%s"%(modification_file,m6A_modifications)
    os.system(cmd)
    tmpfile = os.path.join(out_dir, "%s_tmp.gff" % species)
    if os.path.exists(tmpfile):
        return tmpfile
    else:
        with open(m6A_modifications, 'r') as f:
            t=""
            for line in f.readlines():
                if line[0] != "#":
                    li = line.strip().split()
                    if "m6A" in li:
                        loc=li.index("m6A")
                        li = li[0] + "\t" + li[loc] + "\t" + li[loc+1] + "\t" + li[loc+2] + "\t" + li[loc+3] + "\t" + li[-1]+ "\n"
                        t += li
            generate_outfile(tmpfile,t)
        return tmpfile


def readmodificationfile():
    tmpfile = preprocess_modification_file()
    modifications_filtered = os.path.join(out_dir, "%s_modfilter.txt" % species)
    check_file_exist_rm(modifications_filtered)
    if os.path.exists(tmpfile):
        cmd = "awk '$2 == \"m6A\"{print $0}' %s | sed 's/coverage=//' | sed 's/;IPDRatio=/\t/' | sed 's/;context=/\t/' | sed 's/;frac=/\t/' | sed 's/;fracLow=/\t/' | sed 's/;fracUp=/\t/' | sed 's/;identificationQv=/\t/' | awk '{if($5>=30&&$6>=50) print $1 \"\t\" $2 \"\t\" $3 \"\t\" $4 \"\t\" $5 \"\t\" $6 \"\t\" $7 \"\t\" $8 \"\t\" $9}'|awk '{if($9>1){$9=\"\t\";print} else{print}}'>%s;echo $?" % (tmpfile, modifications_filtered)
    else:
        myprint("Error in cmd or tmpfile!")
    p = os.popen(cmd)
    a = p.read().split()
    if a[0] == "0" and os.path.exists(modifications_filtered) and os.stat(modifications_filtered).st_size > 0:
        return modifications_filtered
    else:
        myprint("Error in command or input modificationfile!")
    


def generate_outfile(outfile, content):
    with open(outfile, 'a') as fo:
        fo.write(content)


def calculate_threshold():
    datasetfile = os.path.join(out_dir, "%s_dataset.txt" % species)
    if os.path.exists(datasetfile):
        cmd = "Rscript threshold.R %s %s %s;echo $?" % (species, out_dir, '%s_dataset.txt' % species)
        p = os.popen(cmd)
        a = p.read().split()
        if a[-1] == "0":
            return a[1]
        else:
            myprint("ERROR:%s"%a)
    else:
        myprint("ERROR:there is not dataset.txt, please run construct_dataset function!")
        sys.exit()



def construct_overlap_dataset():
    peakfilterfile = readpeakfile()
    modfilterfile = readmodificationfile()
    info_list = []
    if os.path.exists(peakfilterfile) and os.stat(peakfilterfile).st_size > 0:
        with open(peakfilterfile, 'r') as fp:
            for line in fp.readlines():
                line = line.strip().split()
                p_chr = line[0]
                p_s = line[1]
                p_e = line[2]
                p_enrich = line[3]
                p_pvalue = line[4]
                p_qvalue = line[5]
                p_info = [p_chr, p_s, p_e, p_enrich, p_pvalue, p_qvalue]
                info_list.append(p_info)

    if os.path.exists(modfilterfile) and os.stat(modfilterfile).st_size > 0:
        with open(modfilterfile, 'r') as fg:
            n = 0
            t = ""
            for l in fg.readlines():
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
                        t += str(chr_name) + '\t' + str(m6a_locus) + '\t' + str(score) + '\t' + str(coverage) + '\t' + str(context) + '\t' + str(fraction) + '\t' + str(ipdratio) + '\t' + str(enrich) + '\n'
            dataset = os.path.join(out_dir, "%s_dataset.txt" % species)
            generate_outfile(dataset, t)
            threshold = calculate_threshold()
            print("the quantity of 6ma events in peaks regions: %d" % n)
            print("the threshold of IPDratio: %s" % threshold)
            return threshold
    else:
        readmodificationfile()
        construct_overlap_dataset()



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
                generate_outfile(outfile, ll)
        return m
       

def find_motifs(motif, motif_loc):
    modfilterfile = readmodificationfile()
    datasetfile = os.path.join(out_dir, "%s_dataset.txt" % species)
    if os.path.exists(modfilterfile) and os.stat(modfilterfile).st_size > 0:
        motifsall = os.path.join(out_dir, "%s_%s_motifsall.txt" % (species,motif))
        check_file_exist_rm(motifsall)
        cmd = "touch %s" % motifsall
        os.system(cmd)
        m = judge_motif(modfilterfile,motifsall,6,motif, motif_loc)
        myprint("the quantity of %s motifs:%d" % (motif, m))
    else:
        myprint("Don't exist file we need,please modify the value of construct_dataset function!")

    if os.path.exists(datasetfile) and os.stat(datasetfile).st_size > 0:
        motifinpeak = os.path.join(out_dir, "%s_%s_motifsinpeak.txt" % (species,motif))
        check_file_exist_rm(motifinpeak)
        cmd = "touch %s" % motifinpeak
        os.system(cmd)
        n = judge_motif(datasetfile,motifinpeak,4,motif, motif_loc)
        myprint("the quantity of %s motif in peaks:%d" % (motif, n))
    else:
        myprint("Don't exist file we need,please modify the value of construct_dataset function!")


def analyse(threshold,motif):
    if threshold:
        motifinpeak = os.path.join(out_dir, "%s_%s_motifsinpeak.txt" % (species,motif))
        motifall = os.path.join(out_dir, "%s_%s_motifsall.txt" % (species,motif))
        
        #total motif sites in genome
        cmd2 = "wc -l %s" % motifall  #所有的motifs
        p2 = os.popen(cmd2)
        n2 = p2.read().split()[0]
        
        #pacbio+MeDIP
        cmd3 = "wc -l %s" % motifinpeak  # peak内的motifs
        p3 = os.popen(cmd3)
        n3 = p3.read().split()[0]
        
        #pacbio+threshold
        cmd4 = "awk '$8>=%d{print}' %s|wc -l" % (float(threshold), motifall)  # 阈值过滤后所有的motifs
        p4 = os.popen(cmd4)
        n4 = p4.read().split()[0]
        
        #pacbio+MeDIP+threshold
        cmd5 = "awk '$7>=%d{print}' %s|wc -l" % (float(threshold), motifinpeak)  # peak内阈值过滤后的motifs
        p5 = os.popen(cmd5)
        n5 = p5.read().split()[0]

        content = "motif:%s" % motif + '\n' + 'Total %s motifs: %s' % (motif, n2) + '\n' + '(pacbio+MeDIP)motifs in peaks: %s' % n3 + '\n'+ '(pacbio+threshold)motifs filtered by threshold: %s' % n4 + '\n'+ '(pacbio+MeDIP+threshold)motifs in peaks filtered by threshold: %s' % n5 + '\n'
        return content
    else:
        print("ERROR:there is not threshold")


def caculate_total_modifications():
    tmpfile = preprocess_modification_file()
    cmd = "wc -l %s" % tmpfile
    p = os.popen(cmd)
    n = p.read().split()[0]
    return n


def caculate_filtered_modifications(threshold):
    tmpfile = preprocess_modification_file()
    mod_thres_filter = os.path.join(out_dir, "mod_thres_filter.txt")
    check_file_exist_rm(mod_thres_filter)
    if os.path.exists(tmpfile):
        cmd = "awk '$2 == \"m6A\"{print $0}' %s | sed 's/coverage=//' | sed 's/;IPDRatio=/\t/' | sed 's/;context=/\t/' | sed 's/;frac=/\t/' | sed 's/;fracLow=/\t/' | sed 's/;fracUp=/\t/' | sed 's/;identificationQv=/\t/' | awk '{if($8>=%d) print $1 \"\t\" $2 \"\t\" $3 \"\t\" $4 \"\t\" $5 \"\t\" $6 \"\t\" $7 \"\t\" $8 \"\t\" $9}'|awk '{if ($9>1){$9=\"\t\";print} else{print}}'>%s;echo $?" % (
            tmpfile, float(threshold), mod_thres_filter)
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

def caculate_peak_modifications():
    datasetfile = os.path.join(out_dir, "%s_dataset.txt" % species)
    if os.path.exists(datasetfile):
        cmd = "wc -l %s" % datasetfile
        p = os.popen(cmd)
        n = p.read().split()[0]
        return n
    else:
        print("datasetfile does not exists!")

def caculate_peak_filtered_modifications(threshold):
    datasetfile = os.path.join(out_dir, "%s_dataset.txt" % species)
    if os.path.exists(datasetfile):
        cmd = "awk '$7>=%d{print}' %s|wc -l" % (float(threshold), datasetfile)
        p = os.popen(cmd)
        n = p.read().split()[0]
        return n
    else:
        print("datasetfile does not exists!")

def main():
    if len(sys.argv) != 3:
        myprint(usage)
        sys.exit()

    myprint('program started')
    myprint('reading config file')
    get_config_values()
    only_check_dir_exist(out_dir)
    only_check_file_exist(configFilePath)
    only_check_file_exist(peak_file)
    only_check_file_exist(modification_file)
    content = '[' + datetime.now().strftime(TimeFormat) + ']'
    if con_func == "TRUE":
        myprint('running tasks: construct_dataset')
        threshold = construct_overlap_dataset()
        n = caculate_total_modifications()
        n1 = caculate_peak_modifications()
        if threshold is None:
            content += '%s:%s(threshold)' % (species,threshold) + '\n' + 'Total sites in genomes: %s' % n + '\n' + 'Total sites filtered by peaks: %s' % n1 + '\n'
        else:
            n2 = caculate_filtered_modifications(threshold)
            n3 = caculate_peak_filtered_modifications(threshold)
            content += '%s:%s(threshold)' % (species,threshold) + '\n' + 'Total sites in genomes: %s' % n + '\n' + 'Total sites filtered by threshold: %s' % n2 + '\n' + 'Total sites filtered by peaks: %s' % n1 + '\n' + 'Total sites filtered by peaks and threshold: %s' % n3 + '\n'
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
        infofile = os.path.join(out_dir, '%s_info.txt' % species)
        generate_outfile(infofile, content)
    else:
        threshold = calculate_threshold()
        n = caculate_total_modifications()
        n1 = caculate_peak_modifications()
        if threshold is None:
            content += '%s:%s(threshold)' % (species,threshold) + '\n' + 'Total sites in genomes: %s' % n + '\n' + 'Total sites filtered by peaks: %s' % n1 + '\n'
        else:
            n2 = caculate_filtered_modifications(threshold)
            n3 = caculate_peak_filtered_modifications(threshold)
            content += '%s:%s(threshold)' % (species,threshold) + '\n' + 'Total sites in genomes: %s' % n + '\n' + 'Total sites filtered by threshold: %s' % n2 + '\n' + 'Total sites filtered by peaks: %s' % n1 + '\n' + 'Total sites filtered by peaks and threshold: %s' % n3 + '\n'
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
        infofile = os.path.join(out_dir, '%s_info.txt' % species)
        generate_outfile(infofile, content)

if __name__ == '__main__':
    main()
