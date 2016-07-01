import numpy as np
import matplotlib.pyplot as plt
import argparse
import csv
import matplotlib
import sys
matplotlib.style.use('ggplot')
font = {'family' : 'Bitstream Vera Sans',
        'weight' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)

def rotate(to_rotate):
    return to_rotate[-1:] + to_rotate[:-1]

def create_bed_color_tree(color_bed):
    """Returns a dictionary for """

def get_args():
    parser = argparse.ArgumentParser(description="Coverage Across a Bed File")
    parser.add_argument("coverage", help="comma list of output from"
                                         " BMFTools depth")
    parser.add_argument("out", help="output file name")
    parser.add_argument("--color_bed", help="bed file to color depth bars",
                        default=None)
    parser.add_argument("--skip_no_cov", help="skips probs/windows with 0 "
                        "coverage", default=None, type=int)
    return parser.parse_args()

def find_gene_ranges(bed, no_cov):
    """gets the range of the bed label from the first input file"""
    curr_gene = None
    gene_range_dict = {}
    with open(bed) as t:
        for ind, line in enumerate(csv.reader(t, delimiter="\t")):
            if "#" in line[0]:
                continue
            if(line[3] == curr_gene):
                gene_range_dict[curr_gene][1] = ind
            else:
                curr_gene = line[3]
                gene_range_dict[curr_gene] = [ind, ind]
    return gene_range_dict


def coverage_plots(doc, outname, color_bed=None, no_cov=None):
    files = doc.split(',')
    #f, axarr = plt.subplots(len(files), sharex=True, figsize=(24,16))
    f, axarr = plt.subplots(len(files), figsize=(24,16))
    gene_ranges = find_gene_ranges(files[0], no_cov)
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off
    for i in range(len(files)):
        uncovered = 0
        coverage = []
        standard_dev = []
        with open(files[i]) as t:
            for line in csv.reader(t, delimiter="\t"):
                if("#" in line [0]):
                    continue
                info = line[4].split("|")[0].split(':')
                cov = float(info[1])
                if no_cov and cov <= no_cov:
                    uncovered += int(line[2]) - int(line[1])
                    continue
                coverage.append(float(info[1]))
                standard_dev.append(float(info[2]))
        ind = np.arange(len(coverage))
        axarr[i].bar(ind, coverage, linewidth=0, color='r')
        tit = files[i].split('/')[-1].split('.')[0]
        axarr[i].set_title(tit)
        axarr[i].set_ylabel("Collapsed Reads")
        heights = [0.7, 0.8, 0.9, 1.0]
        """
        for gene in gene_ranges:
            if(gene == "NamelessInterval"):
                break
            gene_x = np.mean(gene_ranges[gene])
            height = heights[0]
            heights = rotate(heights)
            axarr[i].annotate(gene, xy=(gene_x, max(coverage)*height),
                              xytext=(gene_x, max(coverage)*height),
                              ha='center', va='bottom',
                              bbox=dict(boxstyle="round", fc="w"))
        """
        if(no_cov):
            sys.stdout.write("Total area with %s or less depth of coverage "
                             " for %s: %i \n\n" % (no_cov, tit, uncovered))
    plt.suptitle("AORT Depth of Coverage")
    #plt.xlim(0,len(ind)+1)
    plt.xlabel("Probes")
    plt.savefig(outname, bbox_inches='tight', format="pdf")

    #plt.show()


def main():
    args = get_args()
    coverage_plots(args.coverage, args.out, args.color_bed, args.skip_no_cov)
    return 0


if __name__ == "__main__":
    main()
