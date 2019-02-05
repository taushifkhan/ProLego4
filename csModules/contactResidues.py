import numpy as np  
import pandas as pd
import json
import os, sys
import timeit
import pickle
import structureParser as sP

def __checkPDBFile__(pdbId, pdbDir):
    if pdbDir:
        pdbFl = pdbDir+pdbId+'.pdb'
        if os.path.exists(pdbFl):
            return pdbFl
    else:
        print (pdbId, "pdbFL and directory not found in Path:", pdbDir)
        return 0

vwd_rad = {'C':1.88,'H':1.2,'O':1.42,'N':1.64,'S':1.77,'P':1.8}

def distFunc_2(res1, res2, H):
    d_ = 0
    try:
        for i,row1 in H[H.resNum==res1].iterrows():
            at_1 = vwd_rad[row1['atomType']]
            c_1 = row1['coord']
            for j,row2 in H[H.resNum==res2].iterrows():
                at_2 = vwd_rad[row2['atomType']]
                c_2 = row2['coord']
                d_ = np.sqrt(np.sum((np.array(c_1)-np.array(c_2))**2))
                if (at_1+at_2+0.6) >= d_:
                    # import ipdb; ipdb.set_trace();
                    return round(d_,4)
                else:
                    continue
            # if min(allcont) <= 6:
            #     return min(allcont)
            # else:  
            #     return 0
        return 0
    except:
        import ipdb; ipdb.set_trace();
        return 0            


def distFunc(res1,res2, H, cutoff= 25):
    try:
        v1 =H[H.resNum==res1].coord.values
        v2 = H[H.resNum==res2].coord.values

        for i in v1:
            for j in v2:
                d_ = np.sum((np.array(i)-np.array(j))**2)
                if d_ <=25:
                    return round(np.sqrt(d_),4)
                else:
                    continue
        return 0
    except:
        return 0
    
class contactResidues():
    """
    Calculates all or selective residue distance from a PDB
    @param :
    pdbId --> pdb code
    chain --> specific chain analysis 
    pdbDir --> valid path where PDB Files are saved
    res_set1 --> array of residues 
    res_ser2 --> array of residues
    @return:
    A pickle formatted file
    'Distance file written in ', u'oDir/pdbId_chain.pkl'
    """
    def __init__(self,pdbId='', chain='',pdbDir='',res_set1=[],res_set2=[]):
        if not pdbId and chain:
            print ("Must provide pdbId and Chain ID")
            exit()
        else:
            self.pdbId = pdbId
            self.chain = chain
            pdbFL = __checkPDBFile__(pdbId, pdbDir)
            if pdbFL:
                P = sP.PDB(pdbFL)
                P.to_dataFrame()
                self.pD = P.getChain(chainId=chain)
            else:
                print ("Error in PDB files [Not Found]\nExisting")
                exit()

        self.res_set1 = res_set1
        self.res_set2 = res_set2

    def get_residueDist(self):
        self.residue_contact = []
        if self.res_set1 and self.res_set2:
            for r1 in self.res_set1:
                for r2 in self.res_set2:
                    d = distFunc(r1, r2, self.pD)
                    if d:
                        n1 = self.pD[self.pD.resNum==r1].resName.values[0]
                        n2 = self.pD[self.pD.resNum==r2].resName.values[0]
                        self.residue_contact.append({'resNum1':r1,'resName1':n1,'resNum2':r2,'resName2':n2,'e_dist':d})
                    else:
                        next
            print ("E distance for %d pairs from [%d, %d]"%(len(self.residue_contact),len(self.res_set1),len(self.res_set2)))
        else:
            print ("No residue section provided")
        


    def gen_contactMap(self, minresDist=2):
        all_resi = self.pD.resNum.unique()
        self.all_resiContact = []

        for i in range(len(all_resi)-minresDist):
            r1 = all_resi[i]
            for j in range(i+minresDist,len(all_resi)):
                r2 = all_resi[j]
                d = distFunc_2(r1,r2,self.pD)
                # import ipdb;ipdb.set_trace();
                # d = distFunc(r1,r2,self.pD)
                if d:
                    n1 = self.pD[self.pD.resNum==r1].resName.values[0]
                    n2 = self.pD[self.pD.resNum==r2].resName.values[0]
                    self.all_resiContact.append({'resNum1':r1,'resName1':n1,'resNum2':r2,'resName2':n2,'e_dist':d})
                else:
                    next
            # import ipdb; ipdb.set_trace();

        print ("eDistance computed [=%d] for %d residues [minseqDist = %d]"%(len(self.all_resiContact),len(all_resi),minresDist))


def dump_dictPickl(pId, chain, contMatrix,resNum, aa, elapseTime, oDir):
    outfl = oDir+"%s_%s.pkl"%(pId,chain)
    outDict = {'eDist':contMatrix,'pdbId':pId, 'chain': chain, 'time':elapseTime,'pdbResNum':resNum,'aa':aa}
    pickle.dump(outDict,open(outfl,'w'))
    print ("Distance file written in ", outfl)


def test_contactResidues(pId, chain_id, pdbDir='./',oDir='./'):
    """
    Test contactResidues class
    """
    startT = timeit.timeit()
    cR = contactResidues(pdbId=pId,chain= chain_id,pdbDir=pdbDir)
    # import ipdb; ipdb.set_trace();
    cR.gen_contactMap()
    resNum = cR.pD.resNum.unique()
    aa= [cR.pD[(cR.pD.resNum==i)&(cR.pD.atomName=='CA')].resName.values[0] for i in resNum]
    eT = timeit.timeit()-startT
    dump_dictPickl(pId, chain_id,cR.all_resiContact,resNum, aa, eT, oDir=oDir)

if __name__ == "__main__":
    test_contactResidues(sys.argv[1],sys.argv[2])