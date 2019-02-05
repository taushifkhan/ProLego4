import numpy as np
import pandas as pd
import json
import timeit
import os, sys
from operator import itemgetter
from itertools import groupby
sys.path.append('./')
import angleCalc as aC
import structureParser as sP


# param = {'pdbDir':'./',\
#         'dsspDir':'./',\
#         'outDir':'../csJson/','EE':2,'HH':3,'EH':2,'HE':2}

# json.dump(param,open('./param.json','w'),indent=4)
param = json.load(open('/Users/taushif/gpcr/param.json','r'))

cTYPE = {'EE':'E','HH':'H','EH':'C','HE':'C'}
cCONVERT = {'Ea':1,'Ep':2,'Ha':3,'Hp':4,'Ca':5,'Cp':6}

def checkPDBFile(pdbId, pdbDir):
    if pdbDir:
        pdbFl = pdbDir+pdbId+'.pdb'
        if os.path.exists(pdbFl):
            return pdbFl
    else:
        print (pdbId, "pdbFL and directory not found in Path:", pdbDir)
        return 0
    


def checkFiles(pdbId,pdbDir, dsspDir):
    pdbFl = checkPDBFile(pdbId, pdbDir)
    if pdbFl:
        dsspFl = dsspDir+pdbId+'.dssp'
    else:
        print ("Directory from param")
        pdbFl = param['pdbDir']+pdbId+'.pdb'
        dsspFl = param['dsspDir']+pdbId+'.dssp'

    if os.path.exists(pdbFl) and os.path.exists(dsspFl):
        return pdbFl, dsspFl
    elif os.path.exists(pdbFl):
        runStatus = sP.runDSSP(pdbFl,dsspFl)
        if runStatus:
            return pdbFl,dsspFl
        else:
            print ("DSSP can't be generated: Notfound")
            return 0,0
    else:
        print ("Error in PDB file and dssp file")
        return 0,0

def __makeGroup__(data):
    ranges = {}
    for k, g in groupby(enumerate(data), lambda (i, x):i-x):
        group = map(itemgetter(1), g)
        ranges[group[0]] = (group[0], group[-1])
    return ranges

def get_SSEchunks(dD):
    sse = ['H','E']
    sseChunk = []
    sseArray = np.array(dD.sse.values)
    resArray = np.array(dD.pdbResNum.values)

    hIndex = np.where(sseArray=='H')[0]
    eIndex = np.where(sseArray=='E')[0]
    
    groupIni = []
    groupFinal = []
    hgroup = __makeGroup__(hIndex)
    groupIni.extend(hgroup.keys())
    egroup = __makeGroup__(eIndex)
    groupIni.extend(egroup.keys())

    for i in sorted(groupIni):
        if i in hgroup.keys():
            groupFinal.append(['H',hgroup[i][0],hgroup[i][1]])
        else:
            groupFinal.append(['E',egroup[i][0],egroup[i][1]])
    return groupFinal

def pdist(x,y):
    dist = []
    for i in x:
        for j in y:
            dist.append(np.sum(((np.array(i)-np.array(j))**2)))
    return np.min(dist)

def __get_angle__(dD, sse1,sse2):
    ca1 = []
    ca2 = []

    for i,j in zip(np.arange(sse1[0],sse1[1]),np.arange(sse2[0],sse2[1])):
        try:
            ca1.append(dD[dD.pdbResNum==i].catom.values[0])
            ca2.append(dD[dD.pdbResNum==j].catom.values[0])
        except:
            next

    if len(ca1) and len(ca2):
        aCC = aC.angleCalc()
        # angle,distance = aCC.getAngles(cAlpS1,cAlpS2)  # uses angleCal.calc_angles between sses
        angle = aCC.getAngles(ca1,ca2)

        if abs(angle) <=90:
            orientation = 'p'
        else:
            orientation = 'a'

        return angle, orientation
    else:
        return 0,0

def __getResidueContact(pD,dD,ctype,sse1,sse2):
    dist = []
    for res1 in np.arange(sse1[0],sse1[1]+1):
        try:
            x = pD[pD.resNum==res1].coord.values
            x_res = dD[dD.pdbResNum==res1].aa.values[0]
        except:
            next

        for res2 in np.arange(sse2[0],sse2[1]+1):
            try:
                y = pD[pD.resNum==res2].coord.values
                y_res = dD[dD.pdbResNum==res2].aa.values[0]
            except:
                next
            if (len(x) and len(y) and x_res and y_res):
                xy_dist = pdist(x,y)
                if xy_dist <= 25:
                    dist.append([res1,res2,round(np.sqrt(xy_dist),4),x_res,y_res])
            else:
                next
    if len(dist) >= param[ctype]:
        return 1, dist
    else:
        return 0, dist

            # dist.append(np.min(scipy.spatial.distance.cdist(x,y)))


def check_Contact(pD,dD,sseGroup):
    # contMatrix = np.zeros(shape = (len(sseGroup),len(sseGroup)))
    contDetail = {'csCont':[]}#,'cM':contMatrix}
    for i in range(len(sseGroup)):
        i_grp = sseGroup[i][0]
        for j in range(i+1,len(sseGroup)):
            j_grp = sseGroup[j][0]
            ij_typ = "%s%s"%(i_grp,j_grp)
            tmp = {}
            tmp['ssType'] = ij_typ
            tmp['sse1'] = sseGroup[i]
            tmp['sse2'] = sseGroup[j]
            tmp['angle'] = 0
            tmp['orientation'] = 'x'
            tmp['cs_letter'] = '0'
            dstatus, d_detail = __getResidueContact(pD,dD,ij_typ,sseGroup[i][1:],sseGroup[j][1:])
            tmp['distance'] = d_detail
            if dstatus:
                angle,orientation = __get_angle__(dD, sseGroup[i][1:],sseGroup[j][1:])
                # contDetail['cM'][i,j] = cCONVERT["%s%s"%(cTYPE[ij_typ],orientation)]
                tmp['angle'] = round(angle,3)
                tmp['orientation'] = orientation
                tmp['cs_letter'] = "%s%s"%(cTYPE[ij_typ],orientation)
                # import ipdb; ipdb.set_trace();
            else:
                next
            contDetail['csCont'].append(tmp)

    return contDetail #contMatrix,

class ContactString():
    """
    __doc_string__
    Given a pdb Chain and Id, ContactString has the objects related to SSE 
    contact information
    @param
    pdbId : 4 letter alpah numeric code for PDB
    chain : string one letter code for chain
    @return

    """
    def __init__(self,pdbId='', chain='',pdbDir='',dsspDir=''):
        if not pdbId and chain:
            print ("Must provide pdbId and Chain ID")
            exit()
        else:
            self.pdbId = pdbId
            self.chain = chain
            pdbFL, dsspFL = checkFiles(pdbId, pdbDir, dsspDir)
            if pdbFL and dsspFL:
                # get pdb Chain as dataframe of chain
                P = sP.PDB(pdbFL)
                P.to_dataFrame()
                self.pD = P.getChain(chainId=chain)

                # get dssp chain as dataframe of chain
                D = sP.DSSP(dsspFL)
                D.to_dataFrame()
                self.dD = D.getChain(chain)

            else:
                print ("Error in PDB and DSSP files\nExisting")
                exit()
    
    def makecs(self):
        print (self.pD.shape, self.dD.shape)
        startT = timeit.timeit()
        sseGrouped = get_SSEchunks(self.dD) #([<see>,,start_res, end_res])
        sse = []
        for k in sseGrouped:
            sse.append(k[0])
        
        # contMat,
        contactDetail = check_Contact(self.pD,self.dD,sseGrouped)
        contactDetail.update({'pdbId':self.pdbId,'chain':self.chain})
        contactDetail.update({'aa':''.join(self.dD.aa.values),'sse':''.join(self.dD.sse.values)})
        contactDetail.update({'time':startT-timeit.timeit()})
        contactDetail.update({'sseString':''.join(sse)})
        outfl = param['outDir']+self.pdbId+"_"+self.chain+'.json'
        json.dump(contactDetail,open(outfl,'w'))
        print (outfl)



def test_ContactString(pId, chain_id):
    """
    Test module for contactString Class
    """
    cS = ContactString(pdbId=pId,chain= chain_id)
    import ipdb; ipdb.set_trace();
    cS.makecs()

if __name__ == "__main__":
    test_ContactString(sys.argv[1],sys.argv[2])

