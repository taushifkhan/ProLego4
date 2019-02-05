import pandas as pd
import os, re
import subprocess 


def runDSSP(pdbFL,dsspFl):
    dssp = 'dssp -i %s -o %s'%(pdbFL,dsspFl)
    try:
        subprocess.Popen(dssp, shell= True).wait()
        return 1
    except:
        return 0

class DSSP():
    """
    Parse DSSP file and return a DataFrame
    Example:
    D = DSSP(<dsspFile_path>) # initiate dssp instance check for file
    D.to_dataFrame() --> Parsed file in DataFrame
    D.dsspDF --> Dataframe with all residue information
    D.getDF_chain('chain') {parse DSSP dataframe and put in a chain}
    """
    def __init__(self,dsspFL):
        if os.path.isfile(dsspFL):
            self.flName = dsspFL
            self.fl = open(dsspFL,'r').readlines()
        else:
            print ("Error In File Reading")
            exit

    
    def to_dataFrame(self):
        """
        parse dssp file and put information in pandas dataFrame 
        format.
        """
        dssp = []
        start = 0
        for l in self.fl:
            if re.search('#',l):
                start=1
                continue

            if start:
                try:
                    fields = {'chain':l[10:12].strip(),'seqResNum':int(l[:5].strip()),\
                    'pdbResNum':int(l[5:10].strip()),'aa':l[12:14].strip(),\
                    'acc':float(l[34:38].strip()),\
                    'kappa':float(l[91:97].strip()),'alpha':float(l[97:103].strip()),\
                    'phi':float(l[103:109].strip()),'psi':float(l[109:115].strip()),\
                    'catom':(),'NHO_d':[],'NHO_e':[],'OHN_d':[],'OHN_e':[]}

                    sse_tmp = l[14:17].strip()
                    if re.match('[GHITEBS]', sse_tmp):
                        fields['sse'] = sse_tmp
                    else:
                        fields['sse']  = 'C'
                    
                    # calpha x, y and z
                    x = float(l[117:122].strip())
                    y = float(l[124:129].strip())
                    z = float(l[131:].strip())
                    fields['catom'] = (x, y, z)

                    _h1 = l[40:51].strip().split(",")
                    _h2 = l[51:62].strip().split(",")
                    _h3 = l[62:73].strip().split(",")
                    _h4 = l[73:84].strip().split(",")

                    fields['NHO_d'] = [int(_h1[0]),int(_h3[0])]
                    fields['OHN_d'] = [int(_h2[0]),int(_h4[0])]
                    fields['NHO_e'] = [float(_h1[1]),float(_h3[1])]
                    fields['OHN_e'] = [float(_h2[1]),float(_h4[1])]
                    dssp.append(fields)
                except:
                    print ("Warning residue line escaping", l)
                    next
       
        self.dsspDF = pd.DataFrame(dssp)
        assert len(dssp) == self.dsspDF.shape[0]

    def getChain(self,chain):
        """
        @param:
        chain --> protein chain (string)
        @return:
        DFChain --> DataFrame for a certain chain
        """
        self.to_dataFrame()
        return self.dsspDF[self.dsspDF.chain == chain]

class PDB():
    """
    Parse PDB file and return a DataFrame
    Example:
    P = PDB(<PDBFile_path>,vervose=1) # initiate PDB instance check for file
        <pdbFILE_path> = valid file for pdb 
        vervose = if 1; print mode will be ON else run on silent mode
    P.to_dataFrame() --> Parsed file in DataFrame
    P.pdbDF --> Dataframe with all residue information
    P.getDF_chain('chain') {parse PDB dataframe and put in a chain}
    """
    def __init__(self,pdbFL,vervose=1):
        if os.path.exists(pdbFL):
            self.flName = pdbFL
            self.fl = open(pdbFL,'r').readlines()
            self.vervose = vervose
        else:
            print ("Effor in PDB file/Path")
            exit()

    def to_dataFrame(self):
        tmpArray = []
        for l in self.fl:
            if (l[0:6].strip()=='ATOM' or l[0:6].strip()=='HETATM'):
                tmpDict = {}
                atom_line = l.strip()
                tmpDict['type'] = atom_line[0:6].strip()
                tmpDict['atomNum'] = int(l[6:11].strip())
                tmpDict['atomName'] = atom_line[12:16].strip()
                tmpDict['atomType'] = atom_line[76:78].strip()
                tmpDict['resNum'] = int(atom_line[22:26].strip())
                tmpDict['resName'] = atom_line[17:20].strip()
                tmpDict['chain'] = atom_line[21]
                tmpDict['coord'] = (float(atom_line[30:38]),float(atom_line[38:46]),float(atom_line[46:54]))
                tmpDict['occupancy'] = float(atom_line[54:60].strip())
                tmpDict['tmperatureFac'] = float(atom_line[60:66].strip())
                tmpArray.append(tmpDict)
        # put array dict to pandas dataframe            
        self.pdbDF = pd.DataFrame(tmpArray)
        if self.vervose:
            print (self.pdbDF.shape, self.pdbDF.chain.unique())
    
    def getChain(self,chainId):
        if self.pdbDF.shape[0] and chainId:
            pC = self.pdbDF[self.pdbDF.chain == chainId]
            return pC
        else:
            print ("Indicate chainID from one of the : ",self.pdbDF.chain.unique()) 


if __name__ == "__main__":
    import sys
    P = PDB(sys.argv[1])
    P.to_dataFrame()
    print(P.pdbDF.chain.unique())
    P.getChain(sys.argv[1])
    # D = DSSP(sys.argv[1])
    # dChain = D.getDF_chain('A')
    # print dChain.shape[0], dChain.head()
