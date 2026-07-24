#!/odinn/groups/ant/software/Python_lib/bin/python2.7

import sys
import os
import argparse
import numpy as np
import gzip
import FileTools
import statTools


'''

shebang lines:
linux: #!/odinn/groups/ant/software/Python_lib/bin/python2.7
windows: #!/usr/bin/env python


Read list of IBD fragments between pairs of individuals and write out coverage grid


ibdFragFile=/mnt/c/projects/aDNA/GWS/impGTs/tmpLongFrags_cM20.txt
outGridFile=/mnt/c/projects/aDNA/GWS/impGTs/tmpLongFrags_cM20_grid.txt

ibdFragFile=/mnt/c/projects/aDNA/GWS/impGTs/tmpLongFrags_cM10.txt.gz
outSumFile=/mnt/c/projects/aDNA/GWS/impGTs/tmpLongFrags_cM10_summary.txt
outGridFile=/mnt/c/projects/aDNA/GWS/impGTs/tmpLongFrags_cM10_grid.txt

ibdFragFile=/mnt/c/projects/aDNA/GWS/impGTs/tmpLongFrags_cM5.txt.gz
outSumFile=/mnt/c/projects/aDNA/GWS/impGTs/tmpLongFrags_cM5_summary.txt
outGridFile=/mnt/c/projects/aDNA/GWS/impGTs/tmpLongFrags_cM5_grid.txt



rangeFile=/mnt/c/projects/locusSet/1000Gfilt/1000G_20200805filt_b38_refAl_chromRange.txt
idFile=/mnt/c/projects/aDNA/GWS/impGTs/frags/allICEomni032021_ids.txt
fragColSTR=0:1:2:4:5:7:8
locCntColSTR=10:11:12:13
useID1=SBT-A-1
cMThresh=7
matchCnt=700
diffCnt=60
missCnt=60
diffProp=0.02
fragFile=/mnt/c/projects/aDNA/GWS/impGTs/frags/allICEomni032021_IBD_cM2_maxDiff60_maf0.03_RC_frags2.txt.gz
outStem=/mnt/c/projects/aDNA/GWS/impGTs/frags/${useID1}_cloverICEomni_cM${cMThresh}_match${matchCnt}_diff${diffCnt}_diffProp${diffProp}_${missCnt}
outSumFile=${outStem}_summary.txt

fragIBDproc.py -f ${fragFile} -fc ${fragColSTR} -cc ${locCntColSTR} -if ${idFile} -i1 ${useID1} -os ${outSumFile} -rf ${rangeFile} -rc 0:3:4 -m ${cMThresh} -n ${matchCnt} -x ${missCnt} -z ${diffCnt} -mh ${diffProp}



ID1     ID2     chrom   phase   startPos        stopPos fragLenPos      startcM stopcM  fragLencM       share1Cnt       share2Cnt       diffGTCnt       missGTCnt
SBT-A-1 381621  chr1    1       26389645        29454709        3065065 48.8925780      51.3157210      2.4231  118     179     3       0
SBT-A-1 365420  chr1    1       178022120       181743154       3721035 191.0993900     194.6009300     3.5015  179     415     6       3

ibdFragFile=/mnt/c/projects/aDNA/GWS/impGTs/frags/KOV-A-2_cloverICEomni_IBD_cM3_maxDiff30_maf0.03_0.005_frags_cM10_RC.txt
outSumFile=/mnt/c/projects/aDNA/GWS/impGTs/frags/KOV-A-2_cloverICEomni_cM10_summary.txt
outGridFile=/mnt/c/projects/aDNA/GWS/impGTs/frags/KOV-A-2_cloverICEomni_cM10_grid.txt

ibdFragFile=/mnt/c/projects/aDNA/GWS/impGTs/frags/HSJ-A-1_cloverICEomni_IBD_cM3_maxDiff30_maf0.03_0.005_frags_cM10_RC.txt
outSumFile=/mnt/c/projects/aDNA/GWS/impGTs/frags/HSJ-A-1_cloverICEomni_cM10_summary.txt
outGridFile=/mnt/c/projects/aDNA/GWS/impGTs/frags/HSJ-A-1_cloverICEomni_cM10_grid.txt

#rangeFile=/mnt/c/projects/locusSet/1000Gfilt/1000G_20200805filt_b38_chromBin.txt
rangeFile=/mnt/c/projects/locusSet/1000Gfilt/1000G_20200805filt_b38_refAl_chromRange.txt



fragIBDproc.py -f ${ibdFragFile} -os ${outSumFile}.2 -og ${outGridFile}.2 -rf ${rangeFile} -rc 0:3:4



fragIBD2grid.py -f ${ibdFragFile} -os ${outSumFile} -og ${outGridFile} -nh -i1 DAV-A-9:HSJ-A-1


fragIBDproc.py -f ${ibdFragFile} -os ${outSumFile} -og ${outGridFile}

fragFile=/mnt/c/projects/aDNA/GWS/impGTs/frags/allICEomni032021_IBD_cM2_maxDiff60_maf0.03_RC_frags2.txt.gz
idFile=/mnt/c/projects/aDNA/GWS/impGTs/frags/allICEomni032021_ids.txt
zcat ${fragFile} | awk '{print $2}' | sort | uniq > ${idFile}



'''


def makeChromDict(autoChromCnt=22):
  chromDict = {}
  for i in xrange(autoChromCnt + 1):
    chromDict["chr"+str(i)] = i
    chromDict[str(i)] = i
    #chromDict[i] = i
  chromDict["chrX"] = autoChromCnt + 1
  chromDict["X"] = autoChromCnt + 1
  chromDict["chrXY"] = autoChromCnt + 2
  chromDict["chrY"] = autoChromCnt + 3
  chromDict["Y"] = autoChromCnt + 3
  chromDict["chrM"] = autoChromCnt + 4
  chromDict["MT"] = autoChromCnt + 4
  return chromDict

 
def openFragFile(fragFile): 
  zip = False
  try:
    # try opening as gzip file
    fragStream = gzip.open(fragFile,'r')
    fragHeadLine = fragStream.readline()
    zip = True
    sys.stdout.write("fragFile: %s is gzipped\n" % (fragFile))
  except:
    fragStream = open(fragFile,'r')
    fragHeadLine = fragStream.readline()
    sys.stdout.write("fragFile: %s is not gzipped\n" % (fragFile))
    
  return fragStream, fragHeadLine, zip




def main(args):
 
  '''
  #HSJ-A-1	74350	0	74350	3	1	81463374	117056427	107.6753900000	127.8259900000	1609	3382	35	3
  ID1:ID2:chrom:pos0:pos1:cM0:cM1
  '''

  # get information about fragFiles
  if args.batch == True:
    fragFileList = FileTools.paramParse(args.fragFile)
  else:
    fragFileList = [args.fragFile]
  fragFileCnt = len(fragFileList)
  colNrList = FileTools.getColNrList(args.fragColSTR)

  sys.stdout.write("###############  Running fragIBDproc.py  ###############\n")
  
  sys.stdout.write("--> %d fragment file(s) were specified\n" % (fragFileCnt))
  fileFoundCnt = 0
  for fIdx, fragFile in enumerate(fragFileList):
    fBool = os.path.isfile(fragFile)
    fileFoundCnt += fBool
    sys.stdout.write("%d [%s] %s\n" % (fIdx, fragFile, fBool))
  if fileFoundCnt < fragFileCnt:
    sys.stderr.write("%d of %d specified fragFiles exist. Aborting!\n" % (fileFoundCnt, fragFileCnt))
    return -9 


  sys.stdout.write("Will retain only IBD fragments with:\n -- cM >%g and <%g\n" % (args.minFragSize, args.maxFragSize))

  if args.locCntColSTR != "":
    # column numbers were provided for four counts of loci - share1Cnt, share2Cnt, diffCnt and missCnt 
    locCntColNrList = FileTools.getColNrList(args.locCntColSTR)
    if len(locCntColNrList) == 4:
      sys.stdout.write(" -- matching loci >%d\n -- mismatching loci <%d\n -- missing GTs <%d\n -- mismatching proportion loci <%g\n" % (args.minMatchLocCnt, args.maxDiffCnt, args.maxMissCnt, args.maxFragDiffProp))
    else:
      sys.stderr.write("Incorrect number of colon-delimted column numbers provided for locus counts. Were %d (%s) but should be 4.\n" % (len(locCntColNrList, args.locCntColSTR)))
      return -8
  else:
    sys.stdout.write("No filtering will be done on number of matching, mismatching or missing loci, because --locCntColSTR value was not provided\n")  

  if args.id1Subset != "":
    id1SubSet = set(FileTools.paramParse(args.id1Subset))
    sys.stdout.write("Seeking fragments for %d individual(s) in ID1 column\n" % (len(id1SubSet)))
  else:  
    id1SubSet = set()
    sys.stdout.write("Seeking fragments for all individuals in ID1 column\n")
    


  # get base ID2 information
  id2Dict = {} # stores and provides idx numbers for ID2
  id2List = []
  tmpList = []
  id2ColNrList = FileTools.getColNrList(args.idColSTR)
  id2ColCnt = len(id2ColNrList)
  if id2ColCnt < 1 or id2ColCnt > 2:
    sys.stderr.write("idColSTR must be either single integer or two integers delimited by a colon\n")
    return -2 
  
  id2Stream = open(args.idFile)
  if args.idNoHeader == False:
    headLine = id2Stream.readline()
    
  for id2Idx, line in enumerate(id2Stream):
    f = line.rstrip("\r\n").split("\t")
    id2 = f[id2ColNrList[0]]
    id2List.append(id2)
    id2Dict[id2] = id2Idx
    if id2ColCnt == 2:
      tmpList.append(f[id2ColNrList[1]])
  id2Stream.close()
  id2Cnt = len(id2List)
  sys.stdout.write("Read ID2 values for %d individuals from idFile [%s]\n" % (id2Cnt, args.idFile))
  # check for duplicate IDs
  
  if id2ColCnt == 1:
    tmpList = ["all" for i in range(id2Cnt)]
  id2GrpIdxArr, id2GrpLabArr = FileTools.recodeArray(tmpList)
  tmpList = None
  id2GrpCnt = id2GrpLabArr.size


  #print "---------///", id2Cnt, id2GrpCnt, id2GrpLabArr


  id1Dict = {} # stores and provides idx numbers for ID1
  id1List = []
  id1Cnt = 0
  
  
  #id2List = []
  #id2Cnt = 0
  chromFragCntArr = np.zeros(args.autoChromCnt + 4, dtype=int)
  chromCnt = 0
  fragList = []

  #print chromRangeMtx
  
  if args.winUnit == "kb":
    winSize = args.winSize * 1000
  else:
    winSize = args.winSize
  
  #print "winSize",winSize   
  chrXIdx = args.autoChromCnt + 1
    
    
  shareCntEvalList = FileTools.getColNrList(args.winCntSTR)
  sharePropEvalList = [float(sharePropSTR) for sharePropSTR in args.winPropSTR.split(":")]
  #print "sharePropEvalList", sharePropEvalList
  
  
  chromDict = makeChromDict(autoChromCnt=args.autoChromCnt)
  chromTypeList = [0 for i in range(args.autoChromCnt + 4)] # autosomes 0, chrX 1, other 2
  chromTypeList[0] = 2
  chromTypeList[chrXIdx] = 1 # chrX
  chromTypeList[args.autoChromCnt + 2] = 2 # chrY
  chromTypeList[args.autoChromCnt + 3] = 2 # chrXY
  

  chromRangeMtx = np.zeros((args.autoChromCnt + 4, 2), dtype=int)
  chromRangeMtx[:,0] = 1e10
  
  if args.rangeFile != "":
    rangeColNrList = FileTools.getColNrList(args.rangeColNrSTR)
    rangeStream = open(args.rangeFile)
    if args.winUnit == "cM":
      sys.stdout.write("Read min and max cM values for chromosomes from file [%s]\n" % (args.rangeFile))
    else:
      sys.stdout.write("Read min and max positions for chromosomes from file [%s]\n" % (args.rangeFile))
    for line in rangeStream:
      f = line.rstrip('\r\n').split('\t')
      chrom = int(f[rangeColNrList[0]])
      chromRangeMtx[chrom, 0] = int(float(f[rangeColNrList[1]]) / winSize) # minWinBin
      chromRangeMtx[chrom, 1] = int(float(f[rangeColNrList[2]]) / winSize) # maxWinBin
      #chromRangeMtx[chrom, 0] = minPos
      #chromRangeMtx[chrom, 1] = maxPos
    rangeStream.close()
  else:
    sys.stdout.write("No min and max positions for cM or positions were provided\n")
  sys.stdout.write("#######################################################\n")
   

  #print chromRangeMtx

  for ffIdx, fragFile in enumerate(fragFileList):
    fragStream, fragHeadLine, zip = openFragFile(fragFile)
    if args.fragNoHeader == True:
      fragStream.seek(0)
      
    if fragFileCnt > 1:
      # in case the same ID occurs in different fragFiles
      fSuffix = "%d_" % (ffIdx)
    else:
      fSuffix = ""
    
    for line in fragStream:
      f = line.rstrip('\r\n').split('\t')
      ID1 = f[colNrList[0]]
      if args.id1Subset == "" or ID1 in id1SubSet:
        ID2 = f[colNrList[1]]
        if ID2 in id2Dict:
          id2Idx = id2Dict[ID2]
          #ID2 = "%d_%s" % (ffIdx, f[colNrList[1]])
          #ID2 = "%s%s" % (fSuffix, f[colNrList[1]])
          chrom = chromDict[f[colNrList[2]]]
          pos0 = int(f[colNrList[3]])
          pos1 = int(f[colNrList[4]])
          cM0 = float(f[colNrList[5]])
          cM1 = float(f[colNrList[6]])
          
          fragLen = cM1 - cM0
          if fragLen >= args.minFragSize and fragLen <= args.maxFragSize:
            ## if args.locCntColSTR != "":  
            ## if locus count columns provided then filter based on that!!!
            #print fragLen, f
            useFragBool = True
            if args.locCntColSTR != "":
              # check for number of loci in fragment - and omit if it fails conditions
              matchCnt = int(f[locCntColNrList[0]]) + int(f[locCntColNrList[1]])
              diffCnt = int(f[locCntColNrList[2]])
              missCnt = int(f[locCntColNrList[3]])
              diffProp = diffCnt / float(matchCnt + diffCnt) 
              if matchCnt < args.minMatchLocCnt or diffCnt > args.maxDiffCnt or missCnt > args.maxMissCnt or diffProp > args.maxFragDiffProp:
                useFragBool = False
              #else:  
              #  print useFragBool, matchCnt, diffCnt, missCnt, diffProp, f

            if useFragBool == True:
              if ID1 not in id1Dict:
                id1Dict[ID1] = id1Cnt
                fragList.append([]) # add list to store fragments for ID1
                id1List.append(ID1)
                id1Cnt += 1
              idxID1 = id1Dict[ID1]
              if args.winUnit == "cM":
                bin0 = int(cM0 / winSize)
                bin1 = int(cM1 / winSize)
              else:
                bin0 = int(pos0 / winSize)
                bin1 = int(pos1 / winSize)
              chromFragCntArr[chrom] += 1
              if bin0 < chromRangeMtx[chrom,0]:
                chromRangeMtx[chrom,0] = bin0
              if bin1 > chromRangeMtx[chrom,1]:
                chromRangeMtx[chrom,1] = bin1

              fragList[idxID1].append([chrom, bin0, bin1, fragLen, id2Idx])
            
    fragStream.close()
  
  #print chromRangeMtx
  
  # parse chromosomes
  #print chromRangeDict
  #print "fragList", len(fragList)
  
  for idxID1, ID1 in enumerate(id1List):
    sys.stdout.write("%d fragments for %s\n" % (len(fragList[idxID1]), ID1))

  if len(fragList) == 0:
    sys.stderr.write("No IBD fragments that met the requirements were found. No output will be generated.\n")
    return -9
  
  
  chromList = np.where(chromFragCntArr > 0)[0]
  chromCnt = len(chromList)
  maxChrom = chromList[-1]
  
  chromTypeFragCntArr = np.zeros(3, dtype=int)
  chromTypeFragCntArr[0] = np.sum(chromFragCntArr[1:args.autoChromCnt+1])
  chromTypeFragCntArr[1] = np.sum(chromFragCntArr[chrXIdx])
  chromTypeFragCntArr[2] = np.sum(chromFragCntArr) - (chromTypeFragCntArr[0] + chromTypeFragCntArr[1])
  
  
  #print chromList
  
  #print id1List
  # populate genome grid 
  binCntArr = [0 for chrom in range(maxChrom + 1)]
  binMtx = []
  for chrom in range(maxChrom + 1):
    if chromFragCntArr[chrom] > 0:
      binCnt = (chromRangeMtx[chrom,1] - chromRangeMtx[chrom,0]) + 1
      #print(binCnt)
      #binChromArr += [chrom] * bi
      #binMtx.append(np.zeros((binCnt, id1Cnt), dtype=int))
      binMtx.append(np.zeros((binCnt, id1Cnt, id2GrpCnt), dtype=int))
      binCntArr[chrom] = binCnt
    else:
      binMtx.append([])
    #print binCntArr[chrom], binMtx[chrom].shape
      
  binCntArr = np.array(binCntArr)
  chromTypeBinCntArr = np.zeros(3, dtype=float)
  chromTypeBinCntArr[0] = np.sum(binCntArr[1:args.autoChromCnt+1])
  if chromTypeFragCntArr[1] > 0:
    chromTypeBinCntArr[1] = np.sum(binCntArr[chrXIdx])
  chromTypeBinCntArr[2] = np.sum(binCntArr) - (chromTypeBinCntArr[0] + chromTypeBinCntArr[1])
  
  #print "chromTypeBinCntArr", chromTypeBinCntArr
  #print "chromTypeFragCntArr", chromTypeFragCntArr
  
  #autoBinCnt = float(np.sum(binCntArr[:args.autoChromCnt+1]))
  #print autoBinCnt
  
  #fragCntArr = np.zeros((maxChrom+1, id1Cnt, id2GrpCnt), dtype=int)
  #fragLenMeanArr = np.zeros((maxChrom+1, id1Cnt, id2GrpCnt), dtype=float)
  #fragLenSDArr = np.zeros((maxChrom+1, id1Cnt, id2GrpCnt), dtype=float)
  fragCntArr = np.zeros((3, id1Cnt, id2GrpCnt), dtype=int)
  fragLenMeanArr = np.zeros((3, id1Cnt, id2GrpCnt), dtype=float)
  fragLenSDArr = np.zeros((3, id1Cnt, id2GrpCnt), dtype=float)
  #autoFragStatList = [[None for j in range(id2GrpCnt)] for i in range(id1Cnt)]
  for idxID1, ID1 in enumerate(id1List):
    #print ID1, len(fragList[idxID1])
    chromGrpArr = [[] for i in range(id2GrpCnt)]
    fragLenGrpArr = [[] for i in range(id2GrpCnt)]
    for chrom, bin0, bin1, fragLen, id2Idx in fragList[idxID1]:
      #chromIdx = chromIdxDict[chrom]
      id2GrpIdx = id2GrpIdxArr[id2Idx]
      b0 = bin0 - chromRangeMtx[chrom,0]
      b1 = (bin1 - chromRangeMtx[chrom,0]) + 1
      binMtx[chrom][b0:b1, idxID1, id2GrpIdx] += 1
      chromType = chromTypeList[chrom]
      chromGrpArr[id2GrpIdx].append(chromType)
      fragLenGrpArr[id2GrpIdx].append(fragLen)
      
      
      #print chrom, chromIdx, np.max(binMtx[chromIdx][b0:b1]), b0, b1, bin0, bin1
    for id2GrpIdx in range(id2GrpCnt):
      #chromTypeArr = np.array(chromGrpArr[id2GrpIdx], dtype=int)
      #fragLenArr = np.array(fragLenGrpArr[id2GrpIdx], dtype=float)
      fragLenStat = statTools.calcGrpSumStats(np.array(chromGrpArr[id2GrpIdx], dtype=int), np.array(fragLenGrpArr[id2GrpIdx], dtype=float))
      
      
      #autoFragStatList.append(statTools.calcSumStats(fragLenArr[chromArr <= args.autoChromCnt]))
      #autoFragStatList[idxID1][id2GrpIdx] = statTools.calcSumStats(fragLenArr[chromArr <= args.autoChromCnt])
      for idx, chromTypeSTR in enumerate(fragLenStat.grpLabArr):
        chromType = int(chromTypeSTR)
        fragCntArr[chromType, idxID1, id2GrpIdx] = fragLenStat.grpCntArr[idx]
        fragLenMeanArr[chromType, idxID1, id2GrpIdx] = fragLenStat.grpMeanArr[idx]
        fragLenSDArr[chromType, idxID1, id2GrpIdx] = fragLenStat.grpSDArr[idx]
        
        
      #  if chromFragCntArr[chrom] > 0:
      #    #chromIdx = chromIdxDict[chrom]
      #    fragCntArr[chrom, idxID1, id2GrpIdx] = fragLenStat.grpCntArr[idx]
      #    fragLenMeanArr[chrom, idxID1, id2GrpIdx] = fragLenStat.grpMeanArr[idx]
      #    fragLenSDArr[chrom, idxID1, id2GrpIdx] = fragLenStat.grpSDArr[idx]

    # add code to calculate for all groups in id2GrpCnt combined


  #print fragLenMeanArr
  ###################################################################
  ### calculate the sharing of ID1 x ID2 - and summarize by ID2Grp
  if args.outIndFile != "":
    outIndStream = open(args.outIndFile, 'w')
    outIndStream.write("ID1\tID2\tID2Grp\tfragCnt_auto\tfragLenSum_auto\tfragCnt_chrX\tfragLenSum_chrX\n")
  #(3, id1Cnt, id2GrpCnt)  
  id2GrpFragSumMtx = [[None for j in range(3)] for i in range(id1Cnt)] 
  for idxID1, ID1 in enumerate(id1List):
    id2FragCntArr = np.zeros((3, id2Cnt), dtype=int) # autosomes and chrX
    id2FragSumArr = np.zeros((3, id2Cnt), dtype=float) # autosomes and chrX
    for chrom, bin0, bin1, fragLen, id2Idx in fragList[idxID1]:
      #idxID2 = id2Dict[ID2]
      chromType = chromTypeList[chrom]
      id2FragCntArr[chromType, id2Idx] += 1
      id2FragSumArr[chromType, id2Idx] += fragLen

    if args.outIndFile != "":
      for idxID2, ID2 in enumerate(id2List):
        #if id2FragCntArr[0, idxID2] + id2FragCntArr[1, idxID2] > 0:
        id2GrpIdx = id2GrpIdxArr[idxID2]
        outIndStream.write("%s\t%s\t%s\t%d\t%.3f\t%d\t%.3f\n" % (ID1, ID2, id2GrpLabArr[id2GrpIdx], id2FragCntArr[0, idxID2], id2FragSumArr[0, idxID2], id2FragCntArr[1, idxID2], id2FragSumArr[1, idxID2]))
    for chromType in range(3):
      if chromTypeBinCntArr[chromType] > 0:
        #print "shape", id2FragSumArr[chromType].shape, id2FragSumArr[chromType]
        id2GrpFragSumMtx[idxID1][chromType] = statTools.calcGrpSumStats(id2GrpIdxArr, id2FragSumArr[chromType, :])
  if args.outIndFile != "":
    outIndStream.close()
  
  
  ## get the number of individuals in each population of ID2
  id2GrpCntArr = np.bincount(id2GrpIdxArr)
  
  #minProp = 1. / float(np.min(id2GrpCntArr))
  #print id2GrpCntArr
  #print "minProp", minProp
  
  #print id2GrpFragSumMtx[0][0].grpCntArr, id2GrpFragSumMtx[0][0].grpMeanArr, id2GrpFragSumMtx[0][0].grpSDArr

  # write outfile
  if args.outGridFile != "":
    outGridStream = open(args.outGridFile, 'w')
    outGridStream.write("ID1\tchrom\tstartPos\tstopPos\t%s\n" % ("\t".join(id2GrpLabArr)))

  for chrom in chromList:
    #chromBinCnt = binMtx[chrom].shape[0]
    #binCntArr[chrom] = binMtx[chrom].shape[0]
    #print chrom, binCntArr[chrom], binCntArr[chrom], binMtx[chrom].shape
    #if chrom <= args.autoChromCnt:
    #  autoBinCnt += chromBinCnt
    #grpSumStat = statTools.calcGrpSumStats(grpArr, valArr, calcCI=False, ciSD=1.96, verbose=False)
    if args.outGridFile != "":
      #for winIdx, cntArr in enumerate(binMtx[chrom]):
      for winIdx in range(binCntArr[chrom]):
        #startPos = (winIdx + chromRangeMtx[chrom,0]) * (args.winSize * 1000) + 1
        #stopPos = (((winIdx + 1) + chromRangeMtx[chrom,0]) * (args.winSize * 1000))
        # startPos and stopPos can be physical position or cM
        #print winIdx, chromRangeMtx[chrom,0]
        if args.winUnit == "cM":
          startPos = ((winIdx + chromRangeMtx[chrom,0]) * winSize)
          stopPos = (((winIdx + 1) + chromRangeMtx[chrom,0]) * winSize)
        else:
          startPos = ((winIdx + chromRangeMtx[chrom,0]) * winSize) + 1
          stopPos = (((winIdx + 1) + chromRangeMtx[chrom,0]) * winSize)
        
        #print chromIdx, chrom, startPos, stopPos, cntArr
        for idxID1, ID1 in enumerate(id1List):
          cntArr = binMtx[chrom][winIdx, idxID1]
          outGridStream.write("%s\t%s\t%.10g\t%.10g\t%s\n" % (ID1, chrom, startPos, stopPos, "\t".join(cntArr.astype("S10"))))
      
    #meanShareArr = np.mean(binMtx[chrom], axis=0)
    #stdShareArr = np.std(binMtx[chrom], axis=0)
    
  if args.outGridFile != "":
    outGridStream.close()


  outSummStream = open(args.outSummFile, 'w')
  outSummStream.write("ID\tID2Grp\tchromType\tID2Cnt\tshareID2Mean\tshareID2SD\tfragCnt\tfragLenMean\tfragLenSD\tshareCntMean\tshareCntSD")
  for shareCntEval in shareCntEvalList:
    outSummStream.write("\tpropWinShareCntGE%d" % (shareCntEval))
  outSummStream.write("\tsharePropMean\tsharePropSD")
  for sharePropEval in sharePropEvalList:
    outSummStream.write("\tpropWinSharePropGE%g" % (sharePropEval))
  outSummStream.write("\n")



  if chromTypeBinCntArr[0] > 0:
    ## there are autosomal bins - deal with autosomes for outSummFile
    autoBinMtx = []
    for chrom in chromList:
      if chrom <= args.autoChromCnt:
        #print chrom, binMtx[chrom].shape
        autoBinMtx.append(binMtx[chrom])
    autoBinMtx = np.concatenate(autoBinMtx, axis=0)
    meanShareCntArr = np.mean(autoBinMtx, axis=0)
    stdShareCntArr = np.std(autoBinMtx, axis=0)
    binShareCntList = []
    for shareCntEval in shareCntEvalList:
      binShareArr = np.sum(autoBinMtx >= shareCntEval, axis=0)
      #binShareCntList.append(binShareArr / autoBinCnt)
      binShareCntList.append(binShareArr /chromTypeBinCntArr[0])

    # evaluate the proportion of sharing in windows
    autoBinPropMtx = autoBinMtx / id2GrpCntArr.astype(dtype=float)
    meanSharePropArr = np.mean(autoBinPropMtx, axis=0)
    stdSharePropArr = np.std(autoBinPropMtx, axis=0)
    binSharePropList = []
    for sharePropEval in sharePropEvalList:
      binShareArr = np.sum(autoBinPropMtx >= (sharePropEval / 100), axis=0)
      #binShareCntList.append(binShareArr / autoBinCnt)
      binSharePropList.append(binShareArr /chromTypeBinCntArr[0])
    #print autoBinMtx.shape
    #print autoBinPropMtx
    
    ## empty large numpy arrays 
    autoBinMtx = None
    autoBinPropMtx = None

    for idxID1, ID1 in enumerate(id1List):
      for id2GrpIdx in range(id2GrpCnt):
        outSummStream.write("%s\t%s\tautosomes\t%d\t%.3f\t%.3f\t%d\t%.3f\t%.3f\t%.4f\t%.4f" % (ID1, id2GrpLabArr[id2GrpIdx], id2GrpCntArr[id2GrpIdx], id2GrpFragSumMtx[idxID1][0].grpMeanArr[id2GrpIdx], id2GrpFragSumMtx[idxID1][0].grpSDArr[id2GrpIdx], fragCntArr[0, idxID1, id2GrpIdx], fragLenMeanArr[0, idxID1, id2GrpIdx], fragLenSDArr[0, idxID1, id2GrpIdx], meanShareCntArr[idxID1, id2GrpIdx], stdShareCntArr[idxID1, id2GrpIdx]))
        for idx, shareCntEval in enumerate(shareCntEvalList):
          outSummStream.write("\t%.5f" % (binShareCntList[idx][idxID1, id2GrpIdx]))
        outSummStream.write("\t%.5f\t%.5f" % (meanSharePropArr[idxID1, id2GrpIdx], stdSharePropArr[idxID1, id2GrpIdx]))
        for idx, sharePropEval in enumerate(sharePropEvalList):
          outSummStream.write("\t%.5f" % (binSharePropList[idx][idxID1, id2GrpIdx]))
        outSummStream.write("\n")
    #print "####", chrom, meanShareArr, binCntArr[chrom], binShareArr, binShareArr / float(binCntArr[chrom]), fragCntArr[chrom], fragLenMeanArr[chrom], fragLenSDArr[chrom]


  if chromTypeBinCntArr[1] > 0:
    ## there are chrX bins - deal with chrX for outSummFile
    meanShareCntArr = np.mean(binMtx[chrXIdx], axis=0)
    stdShareCntArr = np.std(binMtx[chrXIdx], axis=0)
    binShareCntList = []
    for shareCntEval in shareCntEvalList:
      binShareArr = np.sum(binMtx[chrXIdx] >= shareCntEval, axis=0)
      #binShareCntList.append(binShareArr / autoBinCnt)
      binShareCntList.append(binShareArr /chromTypeBinCntArr[1])
      
    propMtx = binMtx[chrXIdx] / id2GrpCntArr.astype(dtype=float)
    meanSharePropArr = np.mean(propMtx, axis=0)
    stdSharePropArr = np.std(propMtx, axis=0)
    binSharePropList = []
    for sharePropEval in sharePropEvalList:
      binShareArr = np.sum(propMtx >= sharePropEval, axis=0)
      binSharePropList.append(binShareArr /chromTypeBinCntArr[1])
      
    
    for idxID1, ID1 in enumerate(id1List):
      for id2GrpIdx in range(id2GrpCnt):
        outSummStream.write("%s\t%s\tchrX\t%d\t%.3f\t%.3f\t%d\t%.3f\t%.3f\t%.4f\t%.4f" % (ID1, id2GrpLabArr[id2GrpIdx], id2GrpCntArr[id2GrpIdx], id2GrpFragSumMtx[idxID1][1].grpMeanArr[id2GrpIdx], id2GrpFragSumMtx[idxID1][1].grpSDArr[id2GrpIdx], fragCntArr[1, idxID1, id2GrpIdx], fragLenMeanArr[1, idxID1, id2GrpIdx], fragLenSDArr[1, idxID1, id2GrpIdx], meanShareCntArr[idxID1, id2GrpIdx], stdShareCntArr[idxID1, id2GrpIdx]))
        for idx, shareCntEval in enumerate(shareCntEvalList):
          outSummStream.write("\t%.5f" % (binShareCntList[idx][idxID1, id2GrpIdx]))
        outSummStream.write("\t%.5f\t%.5f" % (meanSharePropArr[idxID1, id2GrpIdx], stdSharePropArr[idxID1, id2GrpIdx]))
        for idx, sharePropEval in enumerate(sharePropEvalList):
          outSummStream.write("\t%.5f" % (binSharePropList[idx][idxID1, id2GrpIdx]))
        outSummStream.write("\n")
    #print "####", chrom, meanShareArr, binCntArr[chrom], binShareArr, binShareArr / float(binCntArr[chrom]), fragCntArr[chrom], fragLenMeanArr[chrom], fragLenSDArr[chrom]
  
 


  #for idxID1, ID1 in enumerate(id1List):
  #  for id2GrpIdx in range(id2GrpCnt):
  #    if autoFragStatList[idxID1][id2GrpIdx].N != None:
  #      #print ID1, autoFragStatList[idxID1].N, autoFragStatList[idxID1].mean, autoFragStatList[idxID1].sd, meanShareArr[idxID1], stdShareArr[idxID1]
  #      outSummStream.write("%s\t%s\tautosomes\t%d\t%d\t%.3f\t%.3f\t%.4f\t%.4f" % (ID1, id2GrpLabArr[id2GrpIdx], id2GrpCntArr[id2GrpIdx], autoFragStatList[idxID1][id2GrpIdx].N, autoFragStatList[idxID1][id2GrpIdx].mean, autoFragStatList[idxID1][id2GrpIdx].sd, meanShareArr[idxID1, id2GrpIdx], stdShareArr[idxID1, id2GrpIdx]))
  #      for idx, shareCntEval in enumerate(shareCntEvalList):
  #        outSummStream.write("\t%.5f" % (binShareCntList[idx][idxID1, id2GrpIdx]))
  #      outSummStream.write("\n")



  outSummStream.close()


if __name__ == '__main__':

  parser = argparse.ArgumentParser(
    description="Produce genome grid summary of IBD between pairs of individuals based on list of IBD fragments. Summary will be produced for ID1 by chromosome and window\n"
                ".")     

  
  parser.add_argument("-f", "--fragFile",
                      help="path of fragFile containing list of IBD fragments to be processed",
                      required=True,
                      type=str,
                      default="")
                      
  parser.add_argument("-b", "--batch",
                      help="fragFile contains list of paths to actual fragFiles",
                      default=False,
                      action="store_true")
                      
  parser.add_argument("-fc", "--fragColSTR",
                      help="Column numbers of ID1:ID2:chrom:pos0:pos1:cM0:cM1 in fragFile (default: %(default)s)",
                      type=str,
                      default="0:1:2:4:5:7:8")

  parser.add_argument("-fh", "--fragNoHeader",
                      help="No header line in fragFile",
                      default=False,
                      action="store_true")


  parser.add_argument("-if", "--idFile",
                      help="path of containing exact list of IDs used to detect IBD fragments",
                      required=True,
                      type=str,
                      default="")

  parser.add_argument("-ic", "--idColSTR",
                      help="Column number of ID2 in idFile (default: %(default)s). A colon-delimited string can also be provided, where the second digit is the column number of a population label in idFile",
                      type=str,
                      default="0")
                      
  parser.add_argument("-ih", "--idNoHeader",
                      help="No header line in idFile",
                      default=False,
                      action="store_true")

  parser.add_argument("-rf", "--rangeFile",
                      help="A tab-delimited text file with the start and stop positions for each chromosome, for example as created by bimChromRange.py",
                      required=True,
                      type=str,
                      default="")

  parser.add_argument("-rc", "--rangeColNrSTR",
                      help="Colon-delimited column numbers of chrom:start:stop for each chromosome in the same winUnit units",
                      type=str,
                      default="0:1:2")
                      
  parser.add_argument("-os", "--outSummFile",
                      help="Summary for each ID1 against all ID2",
                      required=True,
                      type=str,
                      default="")
                      
  parser.add_argument("-oi", "--outIndFile",
                      help="Summary for each ID1 against each ID2",
                      type=str,
                      default="")
                      
  parser.add_argument("-og", "--outGridFile",
                      help="Output file with counts for genome grid",
                      type=str,
                      default="")

  parser.add_argument("-ws", "--winSize",
                      help="Size of windows in genome grid in units of winUnit (default: %(default)g)",
                      type=float,
                      default=0.05)

  parser.add_argument("-wu", "--winUnit",
                      help="Unit of window size  (default: %(default)s)",
                      choices=['cM', 'kb', 'bp'],
                      type=str,
                      default='cM')


  parser.add_argument("-wc", "--winCntSTR",
                      help="Check the proportion of windows that share with >= this number of individuals (default: %(default)s)",
                      type=str,
                      default="1:2:10:50:100")

  parser.add_argument("-wp", "--winPropSTR",
                      help="Check the proportion of windows that share with >= this percent of individuals (default: %(default)s)",
                      type=str,
                      default="0.3:0.5:3:5")

  parser.add_argument("-i1", "--id1Subset",
                      help="Path of file with subset of IDs to process as ID1 or colon-delimited string with IDs",
                      type=str,
                      default="")

  '''parser.add_argument("-or", "--omitRegFile [This has not been implemented!]",
                      help="File with genomic regions to omit from summary",
                      type=str,
                      default="")
                      
  parser.add_argument("-oc", "--omitRegColNrSTR",
                      help="Colon-delimited column numbers of chrom:startPos:stopPos for omitted regions in omitRegFile [This has not been implemented!]",
                      type=str,
                      default="0:1:2")
  '''         
  parser.add_argument("-m", "--minFragSize",
                      help="threshold for minimum cM map size of IBD fragment (default: %(default)g)",
                      type=float,
                      default=1.0)
                      
  parser.add_argument("-mx", "--maxFragSize",
                      help="threshold for maximum map size of IBD fragment (default: %(default)g)",
                      type=float,
                      default=1e20)
  
  parser.add_argument("-cc", "--locCntColSTR",
                      help="colon-delimited string with column number for the following variables: shareCnt1:shareCnt2:diffCnt:missCnt. If no string is provided, then options -n -x, -z and -mh will do nothing. Be very careful about specifying the values for each of these four threshold parameters - based on the nature of the data set and locus distribution!",
                      type=str,
                      default="")

  parser.add_argument("-n", "--minMatchLocCnt",
                      help="threshold for minimum number of loci with matching alleles in IBD fragment [default: %(default)d]. Only applied if locus count columns are given via parameter -cc.",
                      type=int,
                      default=100)

  parser.add_argument("-x", "--maxMissCnt",
                      help="threshold for maximum number of loci with missing GTs in IBD fragment [default: %(default)d]. Only applied if locus count columns are given via parameter -cc.",
                      type=int,
                      default=100000)

  parser.add_argument("-z", "--maxDiffCnt",
                      help="threshold for maximum number of loci with mismatching GTs in IBD fragment [default: %(default)d]. Only applied if locus count columns are given via parameter -cc.",
                      type=int,
                      default=100000)
  
  parser.add_argument("-mh", "--maxFragDiffProp",
                      help="threshold for maximum proportion of mismatching GTs in fragment [default: %(default)g]. Only applied if locus count columns are given via parameter -cc",
                      type=float,
                      default=1.0)

                      
  parser.add_argument("-a", "--autoChromCnt",
                      help="Number of autosomal chromosomes (default: %(default)d)",
                      type=int,
                      default=22)

  args = parser.parse_args()

  if args.minFragSize > args.maxFragSize:
    parser.error("maxFragSize (%f) must be > minFragSize (%f). Aborting\n" % (args.maxFragSize, args.minFragSize))


  main(args)

