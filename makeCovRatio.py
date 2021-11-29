
import sys

def main():
    if len(sys.argv) != 8:
        sys.stderr.write("Usage:  makeCovRatio.py startScale endScale tumor_trcov.txt exomeMetrics_tumor normal_trcov.txt exomeMetrics_normal cov_ratio.txt\n")
        sys.exit(-1)

    try:
        startScale = float(sys.argv[1])
    except ValueError:
        sys.stderr.write("Error:  startScale is not a number:  %s\n" % sys.argv[1])
        sys.exit(-1)

    try:
        endScale = float(sys.argv[2])
    except ValueError:
        sys.stderr.write("Error:  endScale is not a number:  %s\n" % sys.argv[2])
        sys.exit(-1)

    tumorTrcovFile = sys.argv[3]
    tumorMetricsFile = sys.argv[4]
    normalTrcovFile = sys.argv[5]
    normalMetricsFile = sys.argv[6]
    outputFile = sys.argv[7]

    (regions, tumorMean, tumorCov) = getCoverage(tumorTrcovFile, tumorMetricsFile)
    (dummy, normalMean, normalCov) = getCoverage(normalTrcovFile, normalMetricsFile)

    if len(tumorCov) != len(normalCov):
        sys.stderr.write("Error:  Different number of regions between tumor and normal.\n")
        sys.exit(-1)

    normRatio = []
    for i in range(len(regions)):
        if normalCov[i] < 20:
            normRatio.append(1.0)
        else:
            normRatio.append((tumorCov[i] / normalCov[i]) / (tumorMean / normalMean))

    fp2 = open(outputFile, "w")
    fp2.write("## Bin size = 500000 bp\n")
    fp2.write("# chr\tstart\tend\tnormal\ttumor\ttum_over_norm\tnormalized_tum_over_norm\n")

    binsize = 20000

    chr = ""
    start = 0
    end = 0

    idx = 0
    while idx < len(regions):
        if regions[idx][0] != chr:
            chr = regions[idx][0]
            start = 1
            end = binsize

        if regions[idx][2] < start:
            idx += 1
            continue
        
        if end < regions[idx][1]:
            start += binsize
            end += binsize
            continue

        cnt = 0
        ntotal = 0.0
        ttotal = 0.0
        rtotal = 0.0
        j = idx
        while j < len(regions) and regions[j][0] == chr and regions[j][1] < end and start <= regions[j][2]:
            if normalCov[j] >= 20:
                cnt += 1
                ntotal += normalCov[j]
                ttotal += tumorCov[j]
                rtotal += normRatio[j]
            j += 1

        chr2 = chr
        if not chr2.startswith("chr"):
            chr2 = "chr" + chr
        if chr2 == "chrX":
            chr2 = "chr24"
        elif chr2 == "chrY":
            chr2 = "chr25"
        
        scalenum = int(chr2[3:])
        if scalenum == 24:
            scalenum = 7
        elif scalenum == 25:
            scalenum = 22

        incr = (endScale - startScale) / 22
        scaleFactor = startScale + incr * (scalenum - 1)
        if cnt > 0:
            fp2.write("%s\t%s\t%s\t%.2f\t%.2f\t%.3f\t%.3f\n" % \
                      (chr2, start, end, ntotal / cnt, ttotal / cnt, (ttotal / ntotal) / scaleFactor, (rtotal / cnt) / scaleFactor))

        start += binsize
        end += binsize

    fp2.close()

def getCoverage(trcovFile, metricsFile):
     covMean = 0.0
     regions = []
     cov = []

     fp = open(trcovFile, "r")
     while True:
         line = fp.readline()
         if not line:
             break

         if line.startswith("Chrom"):
             continue

         fields = line.strip().split("\t")
         chr = fields[0]
         if chr.startswith("chr"):
             chr = chr[3:]
         if not (chr.isdigit() or chr in ("X", "Y")):
             continue

         regions.append((fields[0], int(fields[1]), int(fields[2])))
         cov.append(float(fields[3]))
     fp.close()

     fp = open(metricsFile, "r")
     while True:
         line = fp.readline()
         if not line:
             break

         if line.startswith("Median coverage"):
             fields = line.strip().split("\t")
             covMean = float(fields[1])
             break
     fp.close()

     if covMean < 1.0:
         sys.stderr.write("Error:  Unable to find mean coverage in exomeMetrics file:  %s\n" % metricsFile)
         sys.exit(-1)

     return (regions, covMean, cov)

if __name__ == '__main__':
    main()
