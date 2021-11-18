
import sys

def main():
    if len(sys.argv) != 7:
        sys.stderr.write("Usage:  p02callCNV ploidy purity sampleName CBS_calling cnvshort.xls cnvfull.xls\n")
        sys.exit(-1)

    try:
        ploidy = float(sys.argv[1])
        purity = float(sys.argv[2])

        if ploidy == 0.5:
            if purity >= 0.8:
                pass
            elif purity >= 0.4:
                ploidy = 0.4
            else:
                ploidy = 0.32

    except ValueError:
        sys.stderr.write("Error:  Invalid ploidy or purity value:  %s, %s\n" % (sys.argv[1], sys.argv[2]))
        sys.exit(-1)

    sampleName = sys.argv[3]
    cbsfile = sys.argv[4]
    outfile = sys.argv[5]
    fullfile = sys.argv[6]

    regions = []
    fp = open(cbsfile)
    while True:
        line = fp.readline()
        if not line:
            break

        if line.startswith("#chr"):
            continue

        f = line.strip().split("\t")
        regions.append([int(f[0]), int(f[1]), int(f[2]), int(f[3]), int(f[4]), float(f[5]), float(f[6])])
    fp.close()

    outfp = open(outfile, "w")
    fullfp = open(fullfile, "w")

    outfp.write("sampleId\tchr\tstart\tend\tlength\t# markers\tmeanLogR\tgainloss\n")
    fullfp.write("sampleId\tchr\tstart\tend\tlength\t# markers\tmeanLogR\tgainloss\n")

    i = 0
    while i < len(regions):
        step = getStep(ploidy, regions[i][5])
        j = i + 1
        while j < len(regions) and regions[j][0] == regions[i][0]:
            jstep = getStep(ploidy, regions[j][5])
            if jstep != step or step != 2:
                break
            j += 1

        chr = "chr" + str(regions[i][0])
        if chr == "chr24":
            chr = "chrX"

        start = regions[i][1]
        end = regions[j-1][2]
        
        cnt = 0
        total = 0.0
        for x in range(i,j):
            cnt += regions[x][4]
            total += regions[x][4] * regions[x][5]

        avg = total / cnt

        if step > 2:
            gainloss = "gain"
        elif step < 2:
            gainloss = "loss"
        else:
            gainloss = "copy-neutral"

        fullfp.write("%s\t%s\t%d\t%d\t%d\t%d\t%.3f\t%s\n" % \
                     (sampleName, chr, start, end, (end - start + 1), cnt, avg, gainloss))
        if step != 2:
            outfp.write("%s\t%s\t%d\t%d\t%d\t%d\t%.3f\t%s\n" % \
                        (sampleName, chr, start, end, (end - start + 1), cnt, avg, gainloss))

        i = j
    outfp.close()
    fullfp.close()

def getStep(step, value):
    if value < 1.0:
        if value >= 1.0 - (step * 0.6):
            return 2
        elif value >= 1.0 - (step * 1.6):
            return 1
        else:
            return 0
    elif value < 1.0 + step:
        if value < 1.0 + (step * 0.6):
            return 2
        else:
            return 3
    else:
        return 2 + int((value - 1.0) / step + 0.5)

if __name__=='__main__':
    main()
 
