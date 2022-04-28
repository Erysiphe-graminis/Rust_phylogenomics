#Function of script is to read in a GFF and fasta file (optionally a list of genes) and then use these to produce modified fasta files and useful information from a set of pileup files
#TO DO - currently script holds entire CDS of each accession in memory - this is v. inefficient, can hold CDS of each genome in memory one at a time in final step, and only remember SNPs up to that point - then supplement with SNPs before writing, then forget
#TO DO - add functionality for taking alignments (eg a phylip from CD-hit and subsequent alignment of gene families) with no gff3
#TO DO - add printing of a .gff3 for each CDS that is reported
#TO DO - add modes for full genes and mrna, currently just does CDS
import os
import argparse

parser = argparse.ArgumentParser(description='Process a series of pileup files into SNPs for downstream processing.')
parser.add_argument("-i", "--inventory", dest="invfile", help="path to inventory file detailing names and locations of pileup files to be used, Software expects a .tsv file with one header line, formatted as Name,Storage_folder,Assembly_folder,Output_folder and then optionally other arguments synergistically with the --parse-inventory option\nthe accession files are expected to be named Accession_aligned.sorted.bam.pileup2snp etc.",action="store", default="na")
parser.add_argument("-g", "--genome", dest="genomefile", help="path to and name of the .fa file encoding the appropriate genome)", action="store", default="na")
parser.add_argument("-a", "--annotations", dest="gff3file", help="path to and name of the .gff3 file describing the appropriate genome)", action="store", default="na")
parser.add_argument("-b", "--breadth", dest="breadth", help="minimum percent of a gene that must be covered (at the minimum depth) for it to be considered for SNP reporting", action="store", default="0")
parser.add_argument("-d", "--depth", dest="depth", help="minimum depth of coverage for a SNP to be considered present, default 6", action="store", default="6")
parser.add_argument("-c", "--coveragedepth", dest="covdepth", help="minimum depth for a position to be considered \"not an N\" - distinct from SNPs to allow writing of CDS files with more freedom for low coverage accessions, default 3", action="store", default="3")
parser.add_argument("-p", "--percentage", dest="minpcnt", help="minimum percent of reads agreeing with a SNP for it to be considered \"real\", default 30", action="store", default="30")
parser.add_argument("-r", "--report-threshold", dest="report_thres", help="Report a specific number of SNPs, prioritising those with the highest score - reports slightly over this number if there are more SNPs with the exact same score. Incompatible with report-score as the output file is the same. Default is 1000, on by default, set to -1 to report all SNPs", action="store",type=int, default=1000)
parser.add_argument("-o", "--output-file", dest="outfile", help="Prefix for the output files to be created", action="store", default="SNP_process_v10")
parser.add_argument("-s", "--report-score", dest="report_score", help="Report all SNPs >= a specific score (scores are determined by the number of accessions, in general the highest possible score is the number of acceessions - incompatible with report-threshold as the output file is the same). Off by default.", action="store",type=int,default=-1)
parser.add_argument("-x", "--gff-offset-start", dest="gffoffsetstart", help="GFF3 files are supposed to be 1-indexed, but if one is 0 indexed (or something else!) you can provide the modifier and the script will adjust accordingly", action="store",type=int,default=0)
parser.add_argument("-y", "--gff-offset-stop", dest="gffoffsetstop", help="GFF3 files are supposed to be 1-indexed, but if one is 0 indexed (or something else!) you can provide the modifier and the script will adjust accordingly, presumably the start and stop offset are the same, but you never know", action="store",type=int,default=0)
parser.add_argument("-N", "--N-weight", dest="Nweight", help="Weighting against Ns for scoring a SNP. Each accession where a SNP appears adds 1 point to the score, maxing out at the number of accessions if there are two perfectly balanced options. Provide a positive int and it will be subtracted from the score to weight in favour of genes which are expressed in all datasets (default is 6)", action="store",type=int,default=6)
parser.add_argument("--force_genelist", dest="genelist", help="Provide a list of genes (one name per line, no \">\" or other character at beginning) and only SNPs from these genes will be computed", action="store", default="na")
parser.add_argument("--parse-inventory", dest="parseinventory", help="Provide instructions for how to parse a more complex inventory file. Formatted such that each column in the inventory file is represented as a member of a comma seperated list, and the required value is given in the list, no value means any will be accepted. For exampled ,,,,,Reference will only parse samples with \"Reference\" written in the seventh column", action="store", default="na")
parser.add_argument("--pre-coverages", dest="print_coverages", help="extract coverages from genes, and store them in a seperate file. Ignore all other positions", action="store_true", default=False)
parser.add_argument("--no-CDS", dest="noCDS", help="Print CDS files for each accession, defaults to on, toggle to off", action="store_true", default=False)
parser.add_argument("--no-coverage", dest="noCoverage", help="Factoring the coverage of a gene into the output CDS, defaults to on, toggle to off. If off, coverage will still be used to check SNP coverage depth but this can be toggled off by setting the required coverage at -1", action="store_true", default=False)
parser.add_argument("--SNP-list", dest="SNPlist", help="Provide a list of SNPs to use, rather than calculating them. The program currently accepts SNPs in the form of Absolute Values, formatted as Chromosome--Base and will only report them if there is a SNP at that position in the dataset", action="store", default="na")
parser.add_argument("--no-indels", dest="indels", help="Flag to turn off in/del processing. If indels are processed then output CDS files will feature genes of different lengths between accessions. If turned off, the reference genome will be processed for SNPs only and all genes will match the reference at other positions. Indels do not currently count for scores, as they tend not to provide actionable information, however they will be reported on. Representing both genes isoforms in a convenient manner is tricky, currently indels are either accepted or not - heterozygous indels are either dropped or reported as the only variant depending on the coverage thresholds. You can process twice with different thresholds to get both, if you like. Defaults to on", action="store_false", default=True)
(options) = parser.parse_args()


def check_options():
    invfile = options.invfile
    if invfile == "na":
        print("You need to specify an inventory file, listing the accessions and paths to their pileup files in the format\nAccession\t/path/to/file.sorted.bam (I assume the files are named sorted.bam.pileup2snp and sorted.bam.pileup2indel)")
        return(1)
    if os.path.isfile(invfile) == False:
        print("Can't reach the file at:\n" + invfile)
        return(1)
    if invfile.split(".")[-1] == "pileup":
        print("You have selected a specific file, rather than an inventory of files, some things might not work quite right")
    genomefile = options.genomefile
    if genomefile == "na":
        print("You need to specify a genome file, including its extension)")
        return(1)
    if os.path.isfile(genomefile) == False:
        print("Can't reach the file at:\n" + genomefile)
        return(1)
    gff3file = options.gff3file
    if gff3file == "na":
        print("You need to specify a gff3 file, including its extension)")
        return(1)
    if os.path.isfile(gff3file) == False:
        print("Can't reach the file at:\n" + gff3file)
        return(1)
    if options.report_score == -1:
        mode = ["threshold",options.report_thres]
    else:
        mode = ["score",options.report_score]
    if not options.genelist == "na":
        if os.isfile(options.genelist) == False:
            print("Unable to reach the genelist at " + options.genelist)
        else:
            print("Will use the provided gene list at " + options.genelist)
    depth = int(options.depth)
    percentage = float(options.minpcnt)
    offsetstart = int(options.gffoffsetstart)
    offsetstop = int(options.gffoffsetstop)
    print("GFF3 file values will all be offset by " + str(offsetstart) + " and " +str(offsetstop) + " for start and stop, respectively")
    #New code to allow discrimination along species, read, and ratio count for the mixed Pt and Pst samples
    if not options.parseinventory == "na":
        inventoryparse = options.parseinventory.split(",")
    else:
        inventoryparse = "na"
    SNPlist = set()
    if not options.SNPlist == "na":
        print("Parsing SNP file at " + options.SNPlist)
        readfile = open(options.SNPlist, "r")
        contents = readfile.readlines()
        for line in contents:
            sline = line.strip().split("\t")
            SNPlist.add(sline[0])
        print("Found " + str(len(SNPlist)) + " SNPs in file")
    return(mode, depth, percentage, offsetstart, offsetstop,inventoryparse, SNPlist)

nucleotide_pairs = {"A":"T","T":"A","G":"C","C":"G","R":"Y","Y":"R","S":"S","W":"W","K":"M","M":"K","B":"V","D":"H","H":"D","V":"B","N":"N",".":"-","-":"-","a":"t","t":"a","g":"c","c":"g","r":"y","y":"r","s":"s","w":"w","k":"m","m":"k","b":"v","d":"h","h":"d","v":"b","n":"n"}

nucleotides = nucleotide_pairs.keys()
nonsynonymousnucleotides = ["A","T","G","C","R","Y","S","W","K","M","B","V","D","H"]

def reverse_complement(sequence):
    revseq = []
    seqlen = len(sequence)
    count = 1
    while count <= seqlen:
        revseq.append(nucleotide_pairs[sequence[seqlen-count]])
        count += 1
    #Joining a string messes it up - have to make sure the thing to be joined is an actual list
    if len(revseq) > 1:
        revseq = str(''.join(revseq))
    elif len(revseq) == 1:
        revseq = revseq[0]
    else:
        print("Error reversing sequence " + sequence + " in " + gene)
    return(revseq)
    

def read_gff3(infile):
    print("reading gff3 file")
    readfile = open(infile, "r")
    contents = readfile.readlines()
    readfile.close()
    genedict = {}
    geneinv = set()
    chrominv = set()
    chromdict = {}
    #dictionary of genes present, keyed by absolute positions in chromosomes - covers the whole gene not just the CDS which is less efficient but should balance future utility with speed compared to searching the full genelist at every position to count coverage
    positiondict = {}
    positioninv = set()
    for line in contents[1:]:
        if not line[0] == "#":
            sline = line.strip().split("\t")
            chrom = sline[0]
            if not chrom in chrominv:
                chrominv.add(chrom)
                chromdict[chrom] = []
            kind = sline[2]
            start = int(sline[3]) + offsetstart
            stop = int(sline[4]) + offsetstop
            strand = sline[6]
            phase = sline[7]
            if phase in ("0","1","2"):
                phase = int(phase)
            else:
                phase = 0
            desc = sline[8].split(";")
            for entry in desc:
                if len(entry)>0:
                    splentry = entry.split("=")
                    if splentry[0] == "ID":
                        featureid = splentry[1]
                    if splentry[0] == "Parent":
                        parentid = splentry[1]
            if kind == "mRNA":
                #creates a list of all genes, and a set of gene lists seperated by chromosome, and then builds the reference dictionary for all genes - "genes" are actually mRNA's to account for transcripts - in future will have to accomodate multiple overlapping transcripts
                #debug report each gene, remove if uneccesary
                #print(featureid)
                #filters by a genelist if relevant
                if featureid in genelist or genelist == []:
                    geneinv.add(featureid)
                    chromdict[chrom].append(featureid)
                    genedict[featureid] = {"chrom":chrom,"strand":strand,"start":start,"stop":stop,"CDScount":0, "CDS":[],"accessions":{}}
                    for position in range(start,stop+1):
                        absolutepos = chrom + "--" + str(position)
                        if not absolutepos in positioninv:
                            positioninv.add(absolutepos)
                            positiondict[absolutepos] = list()
                        positiondict[absolutepos].append(featureid)
            if kind == "CDS":
                if parentid in genelist or genelist == []:
                    genedict[parentid][kind].append([featureid,phase,start,stop,""])
                    genedict[parentid]["CDScount"] +=1

    print("GFF3 file read, reporting\n" + str(len(geneinv)) + " genes spread over " + str(len(chrominv)) + " chromosomes")
    return(genedict,geneinv,chrominv,chromdict,positioninv, positiondict)
        
def extract_from_genome(infile, genedict, chrominv, chromdict, geneinv):
    print("Reading "+ str(len(geneinv)) + " genes from genome into memory")
    readfile = open(infile,"r")
    contents = readfile.readlines()
    contents.append(">END_OF_FILE")
    readfile.close()
    currentchrom = ""
    chromseq = str("")
    chromcount = 0
    for line in contents:
        sline = line.strip()
        if len(sline) > 0:
            if sline[0] == ">":
                if not currentchrom == "":
                    if currentchrom in chrominv:
                        for gene in chromdict[currentchrom]:
                            cdsgene = genedict[gene]["CDS"]
                            fullcds = ""
                            strand = genedict[gene]["strand"]
                            for cds in cdsgene:
                                #GFF3 format lists the start and stop position on the base, meaning that the start position-1 is needed to count in a 0-indexed list
                                start = cds[2]-1
                                stop = cds[3]
                                cds[4] = chromseq[start:stop]
                                if strand == "+":
                                    fullcds += cds[4]
                                #fullcds is oriented in the +ve direction, to facilitate downstream analysis
                                if strand == "-":
                                    fullcds += reverse_complement(cds[4])
                            for accession in accessions:
                                genedict[gene]["accessions"][accession] = {}
                                genedict[gene]["accessions"][accession]["CDS"] = fullcds
                                fullcdslen = len(fullcds)
                                genedict[gene]["accessions"][accession]["maskedCDS"] = ["N"] * fullcdslen
                                genedict[gene]["accessions"][accession]["totalcoverage"] = 0
                            genedict[gene]["accessions"]["Reference"] = {"CDS":fullcds,"totalcoverage":0,"maskedCDS":fullcds}
                    else:
                        print("No genes detected on chromosome " + currentchrom)
                currentchrom = sline.split("|")[0][1:].split(" ")[0]
                chromseq = str("")
                chromcount +=1
                #TO DO - fix this so it outputs a more user friendly breakdown
                if not currentchrom == "END_OF_FILE":
                    print("Chromosome " + str(chromcount) + " out of " + str(len(chrominv)) + " is " + currentchrom)
            else:
                chromseq += sline
    print("Done!\nDicscrepancies indicate \"chromosome\" with no genes on them, or who's names couldn't be parsed using whitespace and \"|\" as seperators")
    return(genedict)


#code for manually overriding the list of genes to process (IE all of them in the gff file) with a shortlist
def provided_genelist(infile):
    readfile = open(infile,"r")
    contents = readfile.readlines()
    genelist = []
    for line in contents:
        sline = line.strip()
        genelist.append(sline)
    return(genelist)



#Use this section to report on SNP data per accession for quality control. This section also replaces per-accession CDS with the "corrected" SNP CDS
def process_SNPs(infile, genedict, chrominv, chromdict, geneinv, depth, minpcnt,accession,SNPdict,SNPgeneinv,SNPgenelist,mode):
    print("opening SNP2pileup file for accession "+ accession)
    infile = infile  + accession + "_aligned.sorted.bam.pileup2snp"
    readfile = open(infile,"r")
    contents = readfile.readlines()
    readfile.close()
    #currently only reports SNPs in chosen genes
    #posprocset is a correction to not overvalue SNPs in multiple isoforms of a gene
    posprocset = set()
    SNPlist = []
    tempsnpcount = 0
    tempgeneinv = set()
    for line in contents[1:]:
        sline = line.strip().split("\t")
        chrom = sline[0].split("|")[0]
        position = int(sline[1])
        reportbase = sline[3]
        reportcov = str(int(sline[4]) + int(sline[5]))
        reportpcnt = sline[6][0:-1]
        refbase = sline[2]
        absolutepos = chrom + "--" + str(position)
        if chrom in chrominv:
            if float(reportpcnt) > minpcnt and int(reportcov) >= depth:
                if absolutepos in positioninv:
                    for gene in positiondict[absolutepos]:
                        if genedict[gene]["start"] <= position and genedict[gene]["stop"] >= position:
                            #sumcdspos allows the relative position in the gene to be determined from the series of start and stop positions
                            sumcdspos = 0
                            #print(chrom)
                            #print(gene)
                            #print(genedict[gene]["accessions"].keys())
                            for cds in genedict[gene]["CDS"]:
                                if cds[2] <= position and cds [3] >= position:
                                    #temppos describes the position within the gene and is used to amend the CDS (it is 0 indexed) - remember that GFF3 lists -ve strand CDS in reverse order, genepos lists the gene and position within the gene, and is used for the later step of counting "informative" SNPs and reporting them , absolutepos describes the position within the whole genome and is output at the end to allow hard comparisons between iterations of the process. It is NOT 0 indexed.
                                    if genedict[gene]["strand"] == "+":
                                        temppos = position - cds[2]
                                    if genedict[gene]["strand"] == "-":
                                        #Also 0 indexed
                                        temppos = cds[3]  - position
                                    genepos = temppos + sumcdspos
                                    #debug reporting of positions
                                    #print(gene)
                                    #print(genedict[gene]["strand"])
                                    #print(absolutepos)
                                    #print(genepos)
                                    #print(temppos)
                                    #print(genedict[gene]["accessions"])
                                    temp = list(genedict[gene]["accessions"][accession]["CDS"])
                                    #print(len(temp))
                                    #print(sumcdspos)
                                    #print(cds[3] - cds[2] +1 )
                                    #print(cds[2])
                                    #print(cds[3])
                                    if genedict[gene]["strand"] == "+":
                                        temp[genepos] = reportbase
                                    if genedict[gene]["strand"] == "-":
                                        temp[genepos] = reverse_complement(reportbase)
                                    temp = str(''.join(temp))
                                    genedict[gene]["accessions"][accession]["CDS"] = temp
                                    #SNPlist is a temporary list, for the purpose of outputting the report at the end on a per-accession basis
                                    #SNPgeneinv is a permanent dictionary, sorted by gene which lists the number of SNPs in genes, and the accessions which they are from. It is for producing reports structured "by gene", either according to a provided list or sorted by "most informative" genes
                                    #SNPdict is a permanent dictionary of SNP positions, with information about the diversity at that position. It does not list the base at that position, to allow the base to be extracted from the "final" CDS and take coverage and indels into account
                                    SNPlist.append([chrom,gene,refbase,reportbase,reportpcnt,reportcov,genepos,absolutepos])
                                    if not gene in SNPgeneinv:
                                        SNPgeneinv.add(gene)
                                        SNPgenelist[gene] = {"sites":{},"SN":0,"accessions":{},"AN":0}
                                    tempsites = SNPgenelist[gene]["sites"]
                                    tempacc = SNPgenelist[gene]["accessions"]
                                    if not genepos in tempsites.keys():
                                        SNPgenelist[gene]["SN"] += 1
                                        tempsites[genepos] = [[],[]]
                                    tempsites[genepos][0].append(accession)
                                    tempsites[genepos][1].append(reportbase)
                                    if not accession in tempacc.keys():
                                        #Currently does nothing - but space for future implementation
                                        tempacc[accession] = []
                                        SNPgenelist[gene]["AN"] += 1
                                    tempsnpcount +=1
                                    tempgeneinv.add(gene)
                                    if not absolutepos in SNPdict.keys():
                                        SNPdict[absolutepos] = {"chrom":chrom,"gene":gene,"genepos":genepos}
                                        for nucl in nucleotides:
                                            SNPdict[absolutepos][nucl] = 0
                                        SNPdict[absolutepos]["N"] = len(accessions)
                                    if not absolutepos in posprocset:
                                        SNPdict[absolutepos][reportbase] += 1
                                        posprocset.add(absolutepos)
                                    #Ns are taken care of later on
                                    #debugging
                                    #if absolutepos == "scaffold_8--3644033":
                                    #    print(SNPdict[absolutepos])
                                    #    print(gene)
                                #adds the length of the cds to the sumcdspos - remember that the cds is 0 indexed so this needs to be the actual number of bases, if the temppos is going to start at "0" 
                                sumcdspos += (cds[3]-cds[2])+1
    print("Done, identified " + str(tempsnpcount) + " valid SNPs out of a total of " + str(len(contents)-1) +" in the file, across " + str(len(tempgeneinv)) + " genes out of a total of " + str(len(geneinv)))
    writetemp = outfilepaths[accession] + options.outfile + "_" + accession + "_SNP_report.txt"
    writefile = open(writetemp,"w")
    writefile.write("Chromosome\tGene\tReference\tReported\tPercentage_Reported\tTotal_Coverage\tPosition_in_gene\tAbsolute_position\tAccession")
    for snp in SNPlist:
        writefile.write("\n" + snp[0])
        for item in snp[1:]:
            writefile.write("\t" + str(item))
        writefile.write("\t"+accession)
    writefile.close()
    return(genedict, chrominv, chromdict, geneinv, depth, minpcnt,accession,SNPdict,SNPgeneinv, SNPgenelist,mode)



#after all SNPs have been collated, report on them using whatever metrics of SNP and gene frequency you like
def report_SNPs(outfile,SNPdict,accessions,mode,genedict):
    #This section is for total reports, not per-accession reports. Need to list number of total SNPs, SNPs at various accession depth (eg, 50%, 40%, 30% etc. - remember a SNP in 90% of accessions is effectively a SNP in 10% of accessions). Then output the SNPs we want to keep into a single .aln file
    tempsnps = SNPdict.keys()
    scores = []
    maxaccessions = len(accessions)
    halfaccessions = float(maxaccessions)/2
    for snp in tempsnps:
        tempscore = 0
        diversitycount = SNPdict[snp]["N"]
        for nucl in nonsynonymousnucleotides:
            temp = float(SNPdict[snp][nucl])
            diversitycount += temp
            tempscore += halfaccessions-abs(temp-halfaccessions)
        tempscore += abs(halfaccessions - abs((maxaccessions-diversitycount)-halfaccessions))
        #weights against Ns - the per gene process currently ignores this and factors in coverage elsewhere
        nscore = (int(SNPdict[snp]["N"]) * options.Nweight)
        tempscore = tempscore-nscore
        scores.append((snp,tempscore))
        #debugging
        #if snp == "scaffold_8--3644033":
        #    print(SNPdict[snp])
    print(str(len(tempsnps)) + " SNPs inventoried, assigned scores to " + str(len(scores)))
    scores.sort(key=lambda pair: pair[1], reverse = True)
    #Code alternates depending on if given a score or an absolute number
    tempscore = float(mode[1])
    currscore = 0
    if mode[0] == "threshold":
        if mode[1] == -1:
            tempscore = float(scores[-1][1])
        else:
            tempscore = float(scores[mode[1]][1])
    tempcount = 0
    writefile = open(outfilepaths["Reference"] + outfile + "_total_SNP_report.phy","w")
    reportfile = open(outfilepaths["Reference"] + outfile + "_SNP_aln_report.txt","w")
    reportfile.write("Absolute_position\tChromosome\tGene\tPosition\tScore\tA,T,G,C,R,Y,S,W,K,M,B,V,D,H,N")
    print("Writing SNP alignments") 
    accessioncount = 0
    maxaccessions = len(accessions) + 1
    fullaccessions = accessions
    fullaccessions.append("Reference")
    #Create a pair of lists that can be sorted together - the way to do this is creating a pair of tuples each time - will have to do a dictionary then iterate through the keys making the tuples afterwards
    genescoresdict = {}
    genescoresset = set()
    for accession in fullaccessions:
        accessioncount +=1
        print("for accession " + accession + " which is "+ str(accessioncount) + " / " + str(maxaccessions))
        currscore = 0
        fullseq = ""
        #modified this section to permit reporting on a per-SNP or a per-gene basis
        for entry in scores:
            snp = entry[0]
            currscore = entry[1]
            tempgene = SNPdict[snp]["gene"]
            if len(genedict[tempgene]["accessions"]["Reference"]["maskedCDS"])>1000 and accessioncount == 1:
                if not tempgene in genescoresset:
                    genescoresset.add(tempgene)
                    genescoresdict[tempgene] = 0
                genescoresdict[tempgene] += (currscore + SNPdict[snp]["N"] * (options.Nweight - 1))
                #debugging
                #if tempgene == "FUN_006540-T1":
                #    print(currscore)
                #    print(genescoresdict[tempgene])
            if (currscore >= tempscore and options.SNPlist == "na") or (snp in SNPlist):
                tempchrom = SNPdict[snp]["chrom"]
                temppos = SNPdict[snp]["genepos"]
                tempbase = genedict[tempgene]["accessions"][accession]["maskedCDS"][temppos]
                fullseq += tempbase
                if accessioncount == 1:
                    highest = 0
                    SNP = ""
                    for nucl in nonsynonymousnucleotides:
                        SNP += str(SNPdict[snp][nucl])
                        SNP += ","
                    SNP += str(SNPdict[snp]["N"])
                    reportfile.write("\n"+ snp +"\t" + tempchrom + "\t" + tempgene + "\t" + str(temppos) + "\t" + str(currscore) + "\t" + SNP)
        if accessioncount == 1:
            writefile.write(str(len(fullaccessions)) + " " + str(len(fullseq)))
        writefile.write("\n" + accession + " " + fullseq)
    writefile.close()
    reportfile.close()
    print("Writing summary of genes with highest SNP counts")
    #write a file summarising the results, it can be used to grab the desired genes from the CDS files
    genescorespairs = []
    for gene in list(genescoresset):
        genescorespairs.append(tuple(list([gene,genescoresdict[gene]])))
    genescorespairs.sort(key=lambda x:x[1], reverse = True)
    writefile = open(outfilepaths["Reference"] + options.outfile + "_Genes_scored_by_SNPs.txt", "w")
    writefile.write("Gene_ID\tGene_length_(bp)\tSum_of_SNP_scores")
    for gene in genescorespairs:
        writefile.write("\n" + gene[0] + "\t" + str(len(genedict[gene[0]]["accessions"]["Reference"]["maskedCDS"])) + "\t" + str(gene[1]))
    writefile.close()

#alternative to report SNPs where the output is a list of genes with high SNP counts - need to adjust to factor in average coverage (eg, if average coverage, or coverage of too many positions is too low, reject)
#def report_variable_CDS(outfile,SNPdict,accessions,mode,genedict):
    



#can use a modified or alternate version of this function to only report genes of interest
def report_all_CDS(accession):
    writefile = open(outfilepaths[accession] + options.outfile + "_" + accession + "_full_CDS.fa","w")
    print("Writing full .CDS file for accession " + accession)
    for chrom in chrominv:
        for gene in chromdict[chrom]:
            #debug gene lister
            #print(gene)
            #print(chrom)
            writefile.write("\n>"+gene)
            #Can add a flag here to do masking or not
            #if masking == True:
            writefile.write("\n" + ''.join(genedict[gene]["accessions"][accession]["maskedCDS"]))
            #else:
            #writefile.write("\n" + genedict[gene]["accessions"][accession]["CDS"])
    writefile.close()


def read_accessions(infile):
    if infile.split(".")[-1] == "pileup":
        accessions = [infile.split("_aligned")[0]]
        accessionpaths = {accessions[0]:infile}
        outfilepaths = {accessions[0]:accessions[0]}
        print("The only accession is " + accessions[0])
    else:
        print("Reading in accessions") 
        readfile = open(infile, "r")
        contents = readfile.readlines()
        accessions = []
        accessionpaths = {"Reference":""}
        outfilepaths = {"Reference":""}
        for line in contents[1:]:
            maxhits = 0
            sline = line.split("\t")
            sline[-1] = sline[-1].strip()
            hits = 0
            if not inventoryparse == "na":
                for count, element in enumerate(inventoryparse):
                    if not element == "":
                        maxhits += 1
                        if sline[count] == element:
                            hits += 1
            if hits == maxhits:
                accession = sline[0]
                accessions.append(accession)
                accessionpaths[accession] = sline[2] + "/"
                outfilepaths[accession] = sline[3] + "/"
            outfilepaths["Reference"] = sline[3] + "/"
        print("Found " + str(len(accessions)) + " accessions")
    return(accessions,accessionpaths,outfilepaths)

def process_coverage(accession, depth, genedict, chromdict, infile):
    #add functionality to output coverage for genes later - to check for bimodal distribution as observed in melania's paper
    if options.print_coverages == True:
        writefilename = infile  + accession + "_coverages.tsv"
        if len(accessions) == 1:
            writefilename = infile + "_coverages.tsv"
        if len(accessions) > 1:
            infile = infile  + accession + "_aligned.sorted.bam.pileup"
        writefile = open(writefilename, "w")
        writefile.write("Chromosome\tPosition\tCoverage")
        print("Reading pileup files for conversion from accession " + accession)
    elif options.noCoverage == True:
        print("ignoring coverage - all CDS will be reported as if there was coverage for accession " + accession)
        for tempchrom in chrominv:
            for gene in chromdict[tempchrom]:
                tempgene = genedict[gene]["accessions"][accession]
                tempgene["maskedCDS"] = tempgene["CDS"]
        return(genedict, chromdict)
    else:
        infile = infile  + accession + "_coverages.tsv"
        print("Reading coverage files for accession " + accession )
    readfile = open(infile, "r")
    contents = readfile.readlines()
    readfile.close()
    count = 0
    maxcount = len(contents)
    posprocset = set()
    for line in contents[1:]:
        sline = line.strip()
        sline = sline.split("\t")
        tempchrom = sline[0]
        position = int(sline[1])
        if options.print_coverages == True:
            tempcov = sline[3]
            tempbase = sline[2]
        else:
            tempbase = "N"
            tempcov = int(sline[2])
        onehit = False
        absolutepos = tempchrom + "--" + str(position)
        if tempchrom in chrominv:
            if absolutepos in positioninv:
                for gene in positiondict[absolutepos]:
                    tempgene = genedict[gene]
                    if options.print_coverages == True:
                        count +=1
                        writefile.write("\n" + tempchrom + "\t" + str(position) + "\t" + tempcov)
                        genedict[gene]["accessions"][accession]["totalcoverage"] += int(tempcov)
                        #this should marginally improve speed by breaking the for loop once a position has been found once (techinically redundant with the onehit == False statement)
                        #Changed to not do this as it deleted all alternate isoforms - not every gene had one, but if they did it unfairly excluded them. Can add as an option if this is an issue later
                    elif options.print_coverages == False:
                        sumcdspos = 0
                        tempstrand = tempgene["strand"]
                        tempgene = genedict[gene]
                        for cds in tempgene["CDS"]:
                            if cds[2] <= position and cds[3] >= position:
                                genedict[gene]["accessions"][accession]["totalcoverage"] += tempcov
                                #needs to allow a position if coverage is good and remove otherwise - remember that 0 coverage produces nothing oftern in the pileup file (only listed as 0 if between other reads)
                                if tempcov >= depth:
                                    if tempstrand == "+":
                                        temppos = position - cds[2]
                                    if tempstrand == "-":
                                        #Also 0 indexed
                                        temppos = cds[3]  - position
                                    genepos = temppos + sumcdspos
                                    genedict[gene]["accessions"][accession]["maskedCDS"][genepos] = genedict[gene]["accessions"][accession]["CDS"][genepos]
                                    count +=1
                                    #debug reporting of positions
                                    #print(gene)
                                    #print(genedict[gene]["strand"])
                                    #print(absolutepos)
                                    #print(genepos)
                                    #print(temppos)
                                    #print(len(temp))
                                    #print(sumcdspos)
                                    #print(cds[3] - cds[2] +1 )
                                    #print(cds[2])
                                    #print(cds[3])
                                    #print(tempcov)
                                #This should be a set, not a call to a function
                                    if absolutepos in SNPdict.keys():
                                        if not absolutepos in posprocset:
                                            SNPdict[absolutepos]["N"] -=1
                                            posprocset.add(absolutepos)
                            sumcdspos += (cds[3]-cds[2])+1
    if options.print_coverages == True:
         print("Found and reported " + str(count) + " positions in genes out of " + str(maxcount) + " total positions")
    else:
        print("Found and reported " + str(count) + " positions in genes out of " + str(maxcount) + " total positions")
    return(genedict, chromdict)

def process_indels()
    readfile = open(infile + accession + "aligned.sorted.bam.pileup2indel", "r")
    print("Processing indels for accession " + accession)
    contents = readfile.readlines()
    readfile.close()
    for line in contents[1:]:
        sline = line.strip().split("\t")





#Inline code starts here (should change to a main function)

SNPdict = {}
SNPgenelist = {}
SNPgeneinv = set()
mode = ""
depth = options.depth
minpcnt = float(25)
proceed = check_options()
if proceed == 1:
    quit()
else:
    mode, depth, minpcnt, offsetstart, offsetstop, inventoryparse, SNPlist = proceed

if not options.genelist == "na":
    genelist = provided_genelist(options.genelist)
else:
    genelist = list()

accessions, accessionpaths, outfilepaths = read_accessions(options.invfile)

genedict, geneinv, chrominv, chromdict, positioninv, positiondict = read_gff3(options.gff3file)



genedict = extract_from_genome(options.genomefile, genedict, chrominv, chromdict, geneinv)
if options.print_coverages == True:
    for accession in accessions:
        genedict, chromdict = process_coverage(accession, options.covdepth, genedict, chromdict, accessionpaths[accession])
    writefile = open(outfilepaths["Reference"] + options.outfile + "_CDS_coverage_summary.tsv", "w")
    for accession in accessions:
        writefile.write("\t" + accession)
    writefile.write("\nTotal_genes")
    for accession in accessions:
        writefile.write("\t" + str(len(geneinv)))
    writefile.write("\nGenes_with_avg_min_coverage")
    covdict = {}
    print("Assessing coverage for all genes")
    for accession in accessions:
        covdict[accession] = {"Total":0,"genes_with_coverage":0}
        tempcount = 0
        print(accession)
        for chrom in chrominv:
            for gene in chromdict[chrom]:
                totalcoverage = genedict[gene]["accessions"][accession]["totalcoverage"]
                genelen = 0
                for cds in genedict[gene]["CDS"]:
                    genelen += int(cds[3])-int(cds[2]) + 1
                avcoverage = float(float(totalcoverage) / float(genelen))
                covdict[accession]["Total"] += totalcoverage
                covdict[accession][gene] = avcoverage
                covdict[gene] = genelen

                if avcoverage > depth:
                    covdict[accession]["genes_with_coverage"] += 1
        writefile.write("\t" + str(covdict[accession]["genes_with_coverage"]))
    writefile.write("\nTotal_coverage_in_genes")
    for accession in accessions:
        writefile.write("\t" + str(covdict[accession]["Total"]))
    writefile.close()
    print("writing all average coverages")
    writefile = open(outfilepaths["Reference"] + options.outfile + "_CDS_average_coverage_per_gene.tsv","w")
    writefile.write("Gene\tGene_length")
    for accession in accessions:
        writefile.write("\t" + accession)
    for chrom in chrominv:
        for gene in chromdict[chrom]:
            writefile.write("\n" + gene + "\t" + str(covdict[gene]))
            for accession in accessions:
                writefile.write("\t" + str(covdict[accession][gene]))
    writefile.close()


else:
    for accession in accessions:
        genedict, chrominv, chromdict, geneinv, depth, minpcnt,accession,SNPdict,SNPgeneinv, SNPgenelist,mode = process_SNPs(accessionpaths[accession], genedict, chrominv, chromdict, geneinv, depth, minpcnt,accession,SNPdict,SNPgeneinv, SNPgenelist,mode)
        genedict, chromdict = process_coverage(accession, int(options.covdepth), genedict, chromdict, accessionpaths[accession])
        if options.noCDS == False:
            report_all_CDS(accession)
    if options.noCDS == False:
        report_all_CDS("Reference")
    report_SNPs(options.outfile,SNPdict,accessions,mode,genedict)




