import argparse
parser = argparse.ArgumentParser(description='Put a description of your script here')
parser.add_argument('-s', '--species', type=str, required=True, help='Path to an input file to be read')
parser.add_argument('-n', '--motif_num', type=str, required=True, help='Path to an input file to be read')
args = parser.parse_args()

def main():
    m = 0
    s = ''
    r1 = ''
    r2 = ''
    l1 = ''
    l2 = ''
    h = ''
    g = ''
    while m < 3:
        m += 1
        #modify the path according to your demand
        file_path = "/xxx/info/info%d/" % m
        file = file_path + "%s_info.txt" % args.species.replace('.','').lower()
        print(args.motif_num)
        dict = extract(file, int(args.motif_num))
        s += str(dict['threshold']) + ','
        r1 += str(dict['increase1']) + ','
        r2 += str(dict['increase2']) + ','
        l1 += str(dict['decrease1']) + ','
        l2 += str(dict['decrease2']) + ','
        h += str(dict['motifs_filtered']) + ','
        g += str(dict['nonmotifs_fitered']) + ','
    c_thres = "%s,%s" % (args.species, s) + '\n'
    generate_outfile("tmp/threshold.csv", c_thres)
    c_increase = "%s,%s%s" % (args.species, r1, r2) + '\n'
    generate_outfile("tmp/increase.csv", c_increase)
    c_decrease = "%s,%s%s" % (args.species, l1, l2) + '\n'
    generate_outfile("tmp/decrease.csv", c_decrease)
    c_motifs = "%s,%s" % (args.species, h) + '\n'
    generate_outfile("tmp/motifs_fitered.csv", c_motifs)
    c_nonmotifs = "%s,%s" % (args.species, g) + '\n'
    generate_outfile("tmp/nonmotifs_fitered.csv", c_nonmotifs)

def generate_outfile(outfile, content):
    with open(outfile, 'a') as fo:
        fo.write(content)


def extract(file,n):
    info = {}
    a = []
    with open(file, 'r') as f:
        for l in f.readlines():
            li = l.strip().split(":")
            a.append(li[-1])
        threshold = a[0].split('(')[0]
        info['threshold'] = threshold
        n_sum = a[1]
        n_peak = a[3]
        n_thres = a[2]
        n_peak_thres = a[4]
        i = 0
        all_motif=0
        all_thres=0
        c_sm = ""
        while i < n:
            i += 1
            motif = a[5*i]
            total = a[5*i+1]
            m_peak = a[5*i+2]
            m_thres = a[5*i+3]
            m_tp = a[5*i+4]
            total_rate = int(total)/int(n_sum)
            peak_rate = int(m_peak)/int(n_peak)
            thres_rate = int(m_thres)/int(n_thres)
            tp_rate = int(m_tp)/int(n_peak_thres)
            c_sm += '%s,%s,%s,%s,%s,%s' % (args.species, motif, total_rate, peak_rate, thres_rate, tp_rate) + '\n'
            all_motif += int(total)
            all_thres += int(m_thres)
        # motifs_filtered = 1 - all_thres/all_motif
        # nonmotifs_fitered = 1 - (int(n_thres)-all_thres)/(int(n_sum)-all_motif)
        motifs_filtered = (all_motif-all_thres)/(int(n_sum)-int(n_thres))
        nonmotifs_fitered = 1 - motifs_filtered
        increase_rate1 = all_motif/int(n_sum)
        increase_rate2 = all_thres/int(n_thres)
        decrease_rate1 = (int(n_sum)-all_motif)/int(n_sum)
        decrease_rate2 = (int(n_thres) - all_thres) / int(n_sum)
        generate_outfile("tmp/single_motif.csv", c_sm)
        info['increase1'] = increase_rate1
        info['increase2'] = increase_rate2
        info['decrease1'] = decrease_rate1
        info['decrease2'] = decrease_rate2
        info['motifs_filtered'] = motifs_filtered
        info['nonmotifs_fitered'] = nonmotifs_fitered
    return info


if __name__ == '__main__':
    main()


