

import os
os.chdir("/sc/orga/projects/chdiTrios/NYGC_25WGS_Trio_30x_2")


def main (bam_dir, out_path):
    """This is the top-level function for this module.
    """

    print bam_dir

    #with open(vcf_path + vcf_file, 'r')) as rec_vcf, open(out_path + "test.csv", "w") as w:
    #    for rec in vcf.Reader(rec_vcf):
    #        w.write(rec)
            #rec_list = [rec.CHROM, rec.POS, rec.REF, rec.ALT]


if __name__ == '__main__':
    import argparse
    #import cProfile

    parser = argparse.ArgumentParser()
    parser.add_argument("bam_dir", default = "/sc/orga/projects/chdiTrios/ASE/newRNA/richtf01/subjunc" help = "bam high level directory.")
    args = parser.parse_args()

    #cProfile.run('main (args.fasta, args.pmf, os.getcwd(), args.run_numb)', 'out.stats')  # for profiling only
    # run to find bottlenecks in code
    main (args.bam_dir, os.getcwd())



# FastQC command
#module load fastqc
#fastqc -f bam -o fastqc YOUR_BAMFILE

#Example : 
#fastqc -f bam -o fastqc /sc/orga/projects/chdiTrios/NYGC_25WGS_Trio_30x_2/1-01094/1-01094-02.final.bam