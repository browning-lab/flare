# This script uses the FLARE .panels output to build clusters, estimate model parameters, and output a FLARE model file.
# The model file includes a comment line with the lagged cross-correlations.
# This script also outputs a log file that includes some input data summaries.

# Usage: python3 create_model_file.py n_ancestries panels_filename output_prefix 

# "n_ancestries" is a positive integer giving the number of clusters (ancestries) to fit to the data
# "panels_filename" is the name of a file that holds panel probabilities output by FLARE
# "output_prefix" is a string, and the output files will be out_prefix.model and out_prefix.log

import sys
from sklearn.mixture import GaussianMixture
import numpy as np

nanc = int(sys.argv[1])
panelfile = open(sys.argv[2])
outprefix = sys.argv[3]
logfile = open(outprefix+".log","w")
outfile = open(outprefix+".model","w")


T=10 # this is just for the model file, and is an initialization for flare

print("input:",sys.argv[3],file=logfile)

# use header of panels file to get nrefhaps and refpops
nrefhaps = int(panelfile.readline().split('=')[1]) # total number of reference haplotypes across all panels that are being used
refpops = panelfile.readline().split()[4:]
npanel = len(refpops)

# read in the panels file
panelprobs = [] # first dimension is location
rhovals = [] # first dimension is location
panelprobschrom = []
rhovalschrom = []
oldlocindex = -1
hapindices = []
oldchrom = -1
for line in panelfile:
    if line[0]=='#': continue
    bits = line.split()
    locindex = int(bits[0]); chrom = bits[1].split(':')[0]; hapindex = int(bits[2]); thisrho = float(bits[3])
    if chrom != oldchrom:
        oldchrom = chrom
        if len(panelprobschrom)>0:
            panelprobs.extend(panelprobschrom)
            rhovals.extend(rhovalschrom)
        panelprobschrom = []
        rhovalschrom = []
    if locindex != oldlocindex:
        oldlocindex = locindex
        panelprobschrom.append([]) # second dimension is haplotype
        rhovalschrom.append([]) # second dimension is haplotype
        if len(hapindices) > 0:
            if any([x!=y for x,y in zip(hapindices,oldhapindices)]):
                print("haplotype indicies inconsistent across positions",file=sys.stderr)
                raise SystemExit
        oldhapindices = [x for x in hapindices]
        hapindices = []
    hapindices.append(hapindex)
    probs = [float(x) for x in bits[4:]]
    if len(probs)!=npanel:
        print("number of panels varies",file=sys.stderr)
        raise SystemExit
    panelprobschrom[len(panelprobschrom)-1].append(probs) # third dimension is panel
    rhovalschrom[len(panelprobschrom)-1].append(thisrho) # thisrho is a scalar
panelprobs.extend(panelprobschrom)
rhovals.extend(rhovalschrom)
nloc = len(panelprobs)


# summary of panels input printed to log file
print(nloc,"locations,",len(hapindices),"haplotypes,",npanel,"reference panels:",' '.join(refpops),",",nrefhaps,"reference haplotypes",file=logfile)

# do the clustering, including a transformation to avoid having entries sum to 1
def cluster(panelprobs,nanc):
    nhaps = len(panelprobs[0])
    # reshape for clustering
    # combine the locations and haplotypes into a single dimension
    # will first have (trimmed) location 1 for all haps first, etc
    clustinput = np.array([hap for location in panelprobs for hap in location])
    # drop one column of clustinput, which is the one with the highest sum
    sums = np.sum(clustinput,0)
    whichdrop = np.argmax(sums)
    transformed = np.delete(clustinput,whichdrop,1)
    # do the clustering/prediction
    myclust = GaussianMixture(nanc,init_params='k-means++',n_init=10)
    myclust.fit(transformed)
    labels = myclust.predict(transformed)
    # reshape results back to have the location dimension (location first dimension and haps second dimension)
    labelmatrix = [labels[i:i+nhaps] for i in range(0,len(labels),nhaps)]
    # put the removed column back in to the parameter (P) matrix
    # rows are ancestries, columns are ref panels
    extra_column = [1.0-sum(x) for x in myclust.means_]
    param = np.insert(myclust.means_,whichdrop,extra_column,axis=1)
    return (labelmatrix, param)


# compute rho (switch rates)
def make_rho(label_matrix,rho_matrix,nanc):
    # convert label matrix to vector matching the rho vector
    labels = [hap for location in label_matrix for hap in location]
    rho_vec = [hap for location in rho_matrix for hap in location]
    if len(labels)!=len(rho_vec):
        print("trouble in make_rho",len(labels),len(rho_vec),file=sys.stderr)
        raise SystemExit
    rho = []
    for i in range(nanc):
        rhoi_vec = [rho_vec[index] for index in range(len(labels)) if labels[index]==i]
        rhoi = sum(rhoi_vec)/len(rhoi_vec)
        rho.append(rhoi)
    return rho

# calculate correlations at lag 
# labelmatrix has location in first dim, haplotype in second dim
def autoc(labelmatrix,nanc,lag):
    acs=dict()
    maxcrosscor= -1.0
    minautocor=1.0
    if len(labelmatrix)<lag+2: return 0.0
    nhaps = len(labelmatrix[0])
    for i1 in range(nanc):
      for i2 in range(i1+1):
        vec0=[]
        veclag=[]
        for j in range(nhaps):
            vec0+=[1*(x[j]==i1) for x in labelmatrix[:-lag]]
            veclag+=[1*(x[j]==i2) for x in labelmatrix[lag:]]
        ac=np.corrcoef(vec0,veclag)[0][1]
        acs[(i2,i1)]=round(ac,3)
        if i1!=i2 and ac > maxcrosscor: maxcrosscor = ac
        if i1==i2 and ac < minautocor: minautocor = ac
    return acs,maxcrosscor,minautocor


thislabelmatrix,thisparams=cluster(panelprobs,nanc)
thisrho = make_rho(thislabelmatrix,rhovals,nanc)

# print out the model file

print("# model is output from cluster_model_singlewin.py",file=outfile)
lag=1
acs,maxcrosscor,minautocor = autoc(thislabelmatrix,nanc,lag)
print("# lag-",lag,"-window correlations are ",acs,file=outfile,sep="")
print("# max cross-correlation is",maxcrosscor,file=outfile)
print("# min auto-correlation is",minautocor,file=outfile)

anclab="anc"
print("# list of ancestries",file=outfile)
pops=""
for i in range(nanc) :
    pops += anclab+"_"+str(i)
    if i < nanc-1: pops +="\t"
print(pops,file=outfile)
print(file=outfile)

print("# list of reference panels",file=outfile)
print("\t".join(refpops),file=outfile)
print(file=outfile)

print("# T: number of generations since admixture",file=outfile)
print(T,file=outfile)
print(file=outfile)

print("# mu[i]: proportion of target genotypes with ancestry i",file=outfile)
print("\t".join([str(1.0/nanc)]*nanc),file=outfile)
print(file=outfile)

print("# p[i][j]: probability that a model state haplotype is in reference panel j",file=outfile)
print("#          when the model state ancestry is i",file=outfile)
for Prow in thisparams:
    print("\t".join([str(x) for x in Prow]),file=outfile)
print(file=outfile)

print("# theta[i][j]: probability that a model state haplotype and the target",file=outfile)
print("#              haplotype carry different alleles when the model state haplotype",file=outfile)
print("#              is in reference panel j and the model state ancestry is i",file=outfile)
mylambda=1/(np.log(nrefhaps)+0.5)
theta=mylambda/(2*mylambda+2*nrefhaps)
for i in range(nanc):
    print("\t".join([str(theta)]*npanel),file=outfile)
print(file=outfile)

print("# rho[i]: rate of the exponential IBD segment cM-length distribution when the",file=outfile)
print("#         most recent common ancestor precedes admixture and has ancestry i",file=outfile)
print("\t".join([str(x) for x in thisrho]),file=outfile)
