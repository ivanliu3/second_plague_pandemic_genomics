#!/usr/bin/env python3

import scipy.spatial
import numpy as np
import pandas as pd
import sys
import os
from time import time
import argparse

'''
python3 pcaDistCalc.mahalanobis.v2.py -tf /maps/projects/ilab/people/pls394/plague/PCA/reference/europe_more3.projection -tp 4:10 -ti 0 -rf /maps/projects/ilab/people/pls394/plague/PCA/reference/europe_more3.projection -rp 4:10 -rg 1 -ri 0
python3 pcaDistCalc.mahalanobis.v2.py -tf /maps/projects/ilab/people/pls394/plague/results_decode/Lund/pcaProj_ALL_scandice_wSaami.forplot.txt -tp 4:10 -ti 0 -rf /maps/projects/ilab/people/pls394/plague/PCA/reference/europe_more3.replace3toScand.projection  -rp 4:10 -rg 1 -ri 0 
'''

class colNrClass:
    def __init__(self):
      self.valColArr = None # variable used for distance calculation
      self.grpCol = None # coloumn for population
      self.indCol = None # column for individual

def main(args): 
  startTime = time()
  sys.stdout.write("#############################  pcaDistCalc.py  #############################\n")

  ### deal with test data
  ## define column numbers for test data
  testPCColNrSTR = args.testPCColNrSTR
  testStartPCColNr, testStopPCColNr = [int(i) for i in testPCColNrSTR.split(':')]
  #testStartPCColNr=4
  #testStopPCColNr=10

  testCol = colNrClass()
  testCol.valColArr = np.arange(testStartPCColNr, testStopPCColNr + 1)
  #testCol.indCol=np.array([0])
  testCol.indCol = np.array( [int(args.testIDColNrSTR)] )
  ## read testCoordFile
  testCoordMtx = [] # matrix of coordinates
  testIndArr = [] # Array of ID
 
  testCoordStream = open(args.testCoordFile)
  testCoordHeadArr = np.array(testCoordStream.readline().rstrip('\r\n').split('\t'))

  testIDHeadSTR = "\t".join(testCoordHeadArr[testCol.indCol])
  for line in testCoordStream:
    f = np.array(line.rstrip('\r\n').split('\t'))
    testCoordMtx.append((f[testCol.valColArr]).astype(dtype=float))
    testIndArr.append(f[testCol.indCol])
      
  testCoordStream.close()
  testCoordMtx = np.array(testCoordMtx)
  testIndArr = np.array(testIndArr).flatten()


  ## read list of test samples
  testIdList = np.array(['LUN100','LUN101'])
  testIdList = testIndArr
  testMatchIndex = np.where(np.in1d(testIndArr, testIdList))[0]
  testCoordMtx = testCoordMtx[testMatchIndex,]
  testIndArr = testIndArr[testMatchIndex,]
  testIdCnt = testCoordMtx.shape[0]


  ### deal with reference data
  ## define column numbers for referecne data
  refPCColNrSTR = args.refPCColNrSTR 
  refStartPCColNr, refStopPCColNr = [ int(i) for i in refPCColNrSTR.split(':') ]
  #refStartPCColNr=4
  #refStopPCColNr=10

  refCol = colNrClass()
  refCol.valColArr = np.arange(refStartPCColNr, refStopPCColNr + 1)
  refCol.indCol = np.array( [int(args.refIDColNrSTR)] )
  refCol.grpCol = np.array( [int(args.refGrpColNrSTR)] )
  #refCol.indCol = np.array([0])
  #refCol.grpCol = np.array([1])

  if refCol.valColArr.size != testCol.valColArr.size:
      sys.stderr.write("Unequal number of columns for values in the reference [%d] and test [%d] data. Aborting!\n" % (refCol.valColArr.size, testCol.valColArr.size))
      sys.exit()
  else:
      pcCnt = (refStopPCColNr + 1) - refStartPCColNr

  ## read refCoordFile                                                                                                                                                                                    
  refCoordMtx = []
  refIndArr = []
  refGrpArr = []
  refCoordStream = open(args.refCoordFile)
  #refCoordStream = open("/maps/projects/ilab/people/pls394/plague/PCA/reference/europe_more3.projection")
  refCoordHeadList = refCoordStream.readline().rstrip('\r\n').split('\t')
  #print refCoordHeadList                                                                                                                                                                                 
  for line in refCoordStream:
    f = np.array(line.rstrip('\r\n').split('\t'))
    refCoordMtx.append((f[refCol.valColArr]).astype(dtype=float))
    refIndArr.append(f[refCol.indCol])
    refGrpArr.append(f[refCol.grpCol])
      
  refCoordStream.close()

  refCoordMtx = np.array(refCoordMtx)
  refIndArr = np.array(refIndArr).flatten()
  refGrpArr = np.array(refGrpArr).flatten() # array of group names for reference individuals


  ### calculate mahalanobis distances between test data and reference data

  ## check if test ind and reference ind have overlaps
  if not  any(np.isin(testIndArr,refIndArr)) : # if no overlap
    print("There is no overlap between test data and reference data" )
    grpArr, refGrpIdxArr, grpCntArr = np.unique(refGrpArr, return_inverse=True, return_counts=True)
    distMalMtx = np.zeros((testIdCnt, grpArr.size), dtype=float)
    pMalMtx = np.zeros((testIdCnt, grpArr.size), dtype=float)

    ctryArr = []
    grpIdxList = []
    for grpIdx, grp in enumerate(grpArr):
      grpIdxList.append(np.where(refGrpIdxArr == grpIdx)[0])
      tmpList = grp.split("_")
      ctryArr.append(tmpList[0])
      
    for grpIdx in range(grpArr.size):
      coordGrpMeanArr = np.mean(refCoordMtx[grpIdxList[grpIdx], :], axis=0)
      try:
        coordGrpCovInvMtx = np.linalg.inv(np.cov(refCoordMtx[grpIdxList[grpIdx], :].T))
        for testIdx, testCoordArr in enumerate(testCoordMtx):
          malDist = scipy.spatial.distance.mahalanobis(testCoordArr , coordGrpMeanArr, coordGrpCovInvMtx)
          distMalMtx[testIdx, grpIdx] = malDist
          ## df = number of PCs
          pMalMtx[testIdx, grpIdx] = scipy.stats.chi2.sf( malDist ** 2, (testCol.valColArr.size) )

      except:
        for testIdx, testCoordArr in enumerate(testCoordMtx):
          distMalMtx[testIdx, grpIdx] = np.nan
          pMalMtx[testIdx, grpIdx] = np.nan



    Pvalue_df = pd.DataFrame(pMalMtx)
    Pvalue_df.columns = grpArr
    Pvalue_df.index = testIndArr
    Pvalue_df.to_csv('data.pvalue.csv', index=True)

    dist_df = pd.DataFrame(distMalMtx)
    dist_df.columns = grpArr
    dist_df.index = testIndArr
    dist_df.to_csv('data.mahadist.csv',index=True)

  else: # if overlap, do leave-one-out test
    print("There are overlaps between test data and reference data. Assume LOO." )  
    ncol =  len(np.unique(refGrpArr)) # number of unique pop groups in reference
    distMalMtx = np.zeros((testIdCnt, ncol), dtype=float)
    pMalMtx = np.zeros((testIdCnt, ncol), dtype=float)
    for testIdx, testId in enumerate(testIndArr): # outer loop: go through test individual (inner loop: pop group in reference)
      print(testId,end='\r')
      testCoordArr = testCoordMtx[testIdx,:]
      index_LOO = np.where(refIndArr != testId)[0]
      refCoordMtx_LOO = refCoordMtx[index_LOO,:]
      refGrpArr_LOO = refGrpArr[index_LOO]
      refIndArr_LOO = refIndArr[index_LOO]
      
      grpArr, refGrpIdxArr, grpCntArr = np.unique(refGrpArr_LOO, return_inverse=True, return_counts=True)
      ctryArr = []
      grpIdxList = []

      for grpIdx, grp in enumerate(grpArr):
        grpIdxList.append(np.where(refGrpIdxArr == grpIdx)[0])
        tmpList = grp.split("_")
        ctryArr.append(tmpList[0])

      for grpIdx in range(grpArr.size): # inner loop: go through group in reference
        coordGrpMeanArr = np.mean(refCoordMtx_LOO[grpIdxList[grpIdx], :], axis=0)
        try:
          coordGrpCovInvMtx = np.linalg.inv(np.cov(refCoordMtx_LOO[grpIdxList[grpIdx], :].T))
        
          malDist = scipy.spatial.distance.mahalanobis(testCoordArr , coordGrpMeanArr, coordGrpCovInvMtx)
          distMalMtx[testIdx, grpIdx] = malDist
          ## df = number of PCs
          pMalMtx[testIdx, grpIdx] = scipy.stats.chi2.sf( malDist ** 2, (testCol.valColArr.size) )

        except:
        
          distMalMtx[testIdx, grpIdx] = np.nan
          pMalMtx[testIdx, grpIdx] = np.nan
        
    Pvalue_df = pd.DataFrame(pMalMtx)
    Pvalue_df.columns = np.unique(refGrpArr)
    Pvalue_df.index = testIndArr
    Pvalue_df.to_csv('data.pvalue.csv', index=True)

    dist_df = pd.DataFrame(distMalMtx)
    dist_df.columns = grpArr
    dist_df.index = testIndArr
    dist_df.to_csv('data.mahadist.csv',index=True)

    
  sys.stdout.write("Finished in total of %.2f seconds\n" % (time()-startTime))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
    description="Calculate Mahalanobis distances in PCA-space between a set of test individuals and a set of reference individuals with group labels\n"
                "\n")

  parser.add_argument("-tf", "--testCoordFile",
                      help="path of file with PC coordinates for test individuals - must have a column containing IDs",
                      required=True,
                      type=str,
                      default="")

  parser.add_argument("-tp", "--testPCColNrSTR",
                      help="For testCoordFile, colon-delimited string with column numbers of first and last PCs to use for calculating distances (zero-based) [default %(default)s].",
                      type=str,
                      default="4:13")

  parser.add_argument("-ti", "--testIDColNrSTR",
                      help="For testCoordFile, colon-delimited string with column numbers of ID and any other variables to be included in output (zero-based) [default %(default)s].",
                      type=str,
                      default="1")

  parser.add_argument("-rf", "--refCoordFile",
                      help="path of file with PC coordinates for reference individuals - must have a column containing group labels",
                      required=True,
                      type=str,
                      default="")

  parser.add_argument("-rp", "--refPCColNrSTR",
                      help="For refCoordFile, colon-delimited string with column numbers of first and last PCs to use for calculating distances (zero-based) [default %(default)s].",
                      required=True,
                      type=str,
                      default="4:13")

  parser.add_argument("-rg", "--refGrpColNrSTR",
                      help="For refCoordFile, colon-delimited string with column numbers of group and any other variables to be included in output (zero-based)  [default %(default)d].",
                      type=int,
                      default=1)  

  parser.add_argument("-ri", "--refIDColNrSTR",
                      help="For refCoordFile, colon-delimited string with column numbers of ID and any other variables to be included in output (zero-based) [default %(default)s].",
                      type=str,
                      default="1")

  args = parser.parse_args()
  main(args)          
