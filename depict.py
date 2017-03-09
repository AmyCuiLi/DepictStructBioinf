
PROJECT_CASSETTE_DIR = './project_cassette/'


import sys
sys.path.append(PROJECT_CASSETTE_DIR)
from Bio import AlignIO
import Bio.PDB
from Families import *
from Annotations import *
import argparse
import MDAnalysis as MDA
import numpy as np
import os


structDir = './%s/aligned_structures/' %(PROJECT_CASSETTE_DIR)
alnFile = './%s/SequenceAlignment.fasta' %(PROJECT_CASSETTE_DIR)



class HomAlignment:
    def __init__(self):
        #self.alnFile = '../2-alignment/allSequences_aydin_alignment.fasta'
        #self.alnFile = '../2-alignment/all_aligned_aydin.fasta'
        self.alnFile = alnFile
        with open(self.alnFile) as fo:
            self.alnObj = AlignIO.read(fo, 'fasta')

                               
    def getHomResidue(self, fromSeqNum, fromSeqName, toSeqName):
        #print 'getHomResidue called. fromSeqNum:', fromSeqNum
        #toSeqName = 'a3g_ctd'
        fromSeq = self.getSeqByName(fromSeqName)
        toSeq = self.getSeqByName(toSeqName)
        #print 'fromSeq, toSeq', fromSeq.id, toSeq.id
        # So now both sequences have been found
        alignmentCol = self.indexWOGaps2WGaps(fromSeqNum, fromSeq)
        print 'alignmentCol', alignmentCol
        newSeqNum = self.indexWGaps2WOGaps(alignmentCol, toSeq)
        print 'in getHomResidue: fromSeqNum, fromSeqName, toSeqName, newSeqNum', fromSeqNum, fromSeqName, toSeqName, newSeqNum
        if (newSeqNum is None):
            return None
        else:
            return (toSeq[alignmentCol], newSeqNum)
        
    def getSeqByName(self, seqName):
        searchSeq = None
        for sequenceObj in self.alnObj:
            if (sequenceObj.id == seqName):
                if not(searchSeq is None):
                    raise Exception('Two sequences with the same ID found in alignment. fromSeq: %r  sequenceObj: %r' %(searchSeq, sequenceObj))
                searchSeq = sequenceObj
        if searchSeq is None:
            raise Exception('Sequence %s not found in alignment' %(seqName))
        #print "in getSeqByName: seqName searchSeq(to return)", seqName, searchSeq
        return searchSeq
        

    def indexWOGaps2WGaps(self, seqIndex, seqWGaps):
        # NOTE: INDEXES BEGINNING AT 0
        #print 'in indexWOGaps2WGaps: seqIndex, seqWGaps.id', seqIndex, seqWGaps.id
        ngposition = -1
        for gposition, gletter in enumerate(seqWGaps):
            if gletter != '-':
                ngposition += 1
            if ngposition == seqIndex:
                return gposition
        #raise Exception('End of sequence reached without finding corresponsing non-gapped index. index %r seqWGaps %r' %(index, seqWGaps))
        return None


    def indexWGaps2WOGaps(self, alignmentCol, seqWGaps):
        # NOTE: INDEXES BEGINNING AT 0
        ngposition = -1
        #print 'seqWGaps.id', seqWGaps.id
        #print 'in indexWGaps2WOGaps: seqIndex, seqWGaps.id', seqIndex, seqWGaps.id
        for gposition, gletter in enumerate(seqWGaps):
            #print gletter, gletter != '-', ngposition
            if gletter != '-':
                ngposition += 1
            if gposition == alignmentCol:
                if gletter == '-':
                    return None
                else:
                    return ngposition
        #raise Exception('End of sequence reached without finding corresponsing non-gapped index. index %r seqWGaps %r' %(index, seqWGaps))
        return None
    
    def translateAnnotation(self, annDict, newParentStruct):
        oldParent = annDict['parent']
        ## If it's a different parent, translate all of the sequence positions to the new parent's numbering
        newPosnsAAs = [self.getHomResidue(i,oldParent, newParentStruct) for i in annDict['posns']]
        newPosnsAAs = [i for i in newPosnsAAs if not(i is None)]
        print 'newPosnsAAs', newPosnsAAs
        newAnnDict = {}
        newAnnDict['parentStruct'] = newParentStruct
        if oldParent  != newParentStruct:
            newAnnDict['label'] = 'Homologous to:', annDict['label']
        else:
            newAnnDict['label'] = annDict['label']
        newAnnDict['checkAA'] = [i[0] for i in newPosnsAAs]
        newAnnDict['posns'] = [i[1] for i in newPosnsAAs]
        # If none of the residues map onto the new structure, return None
        if len(newPosnsAAs) == 0:
            return None
        return newAnnDict

    
class GeoTools:
    def __init__(self):
        self.structDir = structDir

    def getGeoCenter(self, annDict):
        parentStructFile = self.structDir + annDict['parentStruct'] + '.pdb'
        print annDict['parentStruct']
        U = MDA.Universe(parentStructFile)
        seqPos2ResNum = self.getSeqPos2ResNumFromU(U)
        resnumStr = ' or '.join(['resnum '+str(seqPos2ResNum[posn]) for posn in annDict['posns']])
        sel = U.selectAtoms('( %s ) and name CA' %(resnumStr))
        geoCenter = np.mean(sel.coordinates(), axis=0)
        return geoCenter

    def getPtnGeoCenter(self,annDict):
        parentStructFile = self.structDir + annDict['parentStruct'] + '.pdb'
        U = MDA.Universe(parentStructFile)
        sel = U.selectAtoms('name CA')
        geoCenter = np.mean(sel.coordinates(), axis=0)
        return geoCenter

    def getSeqPos2ResNum(self, structName):
        structFileName = self.structDir + structName + '.pdb'
        U = MDA.Universe(structFileName)
        seqPos2ResNum = self.getSeqPos2ResNumFromU(U)
        return seqPos2ResNum


    def getSeqPos2ResNumFromU(self, U):
        #resnums = U.selectAtoms('protein').residues.resids()
        resnums = U.selectAtoms('protein').residues.resids
        seqPos2ResNum = dict(list(enumerate(resnums)))
        return seqPos2ResNum


    def getSeqPos2ResLet(self, structName):
        structFileName = self.structDir + structName + '.pdb'
        U = MDA.Universe(structFileName)
        seqPos2ResLet = self.getSeqPos2ResLetFromU(U)
        return seqPos2ResLet
    
    def getSeqPos2ResLetFromU(self, U):
        #resnames = [Bio.PDB.Polypeptide.three_to_one(i) for i in U.selectAtoms('protein').residues.resnames()]
        resnames = [Bio.PDB.Polypeptide.three_to_one(i) for i in U.selectAtoms('protein').residues.resnames]
        seqPos2ResLet = dict(list(enumerate(resnames)))
        return seqPos2ResLet
        
class PymolScriptWriter:
    def __init__(self):
        self.scriptText = ''
        self.structPaNames = {}
        self.distNames = []
        self.GT = GeoTools()
        # Record the number of annotations already on each structure, to ensure each annotation is a different color
        self.structAnnCounter = {}
        # A dictionary of unique coloring schemes
        self.colorSchemes = {1:'5',
                             2:'6',
                             3:'9',
                             4:'144',
                             5:'11',
                             6:'13',
                             7:'33',
                             8:'154',}
        self.solidColors = {1:'2',
                            2:'3',
                            3:'6',
                            4:'8',
                            5:'5',
                            6:'13',
                            7:'4',
                            8:'7',
                            9:'9',
                            10:'10'}
        
    def loadStruct(self, structName):
        structFilename = structDir + structName + '.pdb'
        self.scriptText += 'load %s \n' %(structFilename)
        self.structAnnCounter[structName] = 0

    def getNewPANameForStruct(self, structName):
        ## Get an unused pseudoatom name for this structure. 
        nameTaken = True
        index = 0
        while 1:
            candidateName = 'p%i' %(index)
            if not(structName in self.structPaNames.keys()):
                self.structPaNames[structName] = [candidateName]
                return candidateName
            if not(candidateName in self.structPaNames[structName]):
                self.structPaNames[structName] += [candidateName]
                return candidateName
            index += 1
            if index == 100:
                raise Exception('Pseudoatom numbering reached 100 without finding unused atom name.')

    def getNewDistName(self):
        ## Get an unused pseudoatom name for this structure. 
        nameTaken = True
        index = 0
        while 1:
            candidateName = 'D%i' %(index)
            if not(candidateName in self.distNames):
                self.distNames += [candidateName]
                return candidateName
            index += 1
            #if index == 100:
            #    raise Exception('Pseudoatom numbering reached 100 without finding unused atom name.')


    def addAnn(self, fosterAnn, addAnnLines=True, addAnnObj=True):
        parStruct = fosterAnn['parentStruct']
        self.structAnnCounter[parStruct] = self.structAnnCounter.get(parStruct,0) + 1
        if addAnnLines:
            self.addAnnLines(fosterAnn)
        if addAnnObj:
            self.addAnnObj(fosterAnn)
                
    def addAnnLines(self,fosterAnn):
        structName = fosterAnn['parentStruct']
        #seqOffset = self.GT.getStartResFromStructName(structName)
        seqPos2ResNum = self.GT.getSeqPos2ResNum(structName)
        #print structName
        #print seqPos2ResNum
        #print 'In addAnn: structName', structName, 'seqOffset', seqOffset
        #annStructName = fosterAnn['parentStruct'] + '_ann'
        annStructName = structName
        geoCenterStr = '[%f, %f, %f]' %(fosterAnn['labelLoc'][0],
                                        fosterAnn['labelLoc'][1],
                                        fosterAnn['labelLoc'][2])
        labelText = fosterAnn['label']

        # Create lines object to annotation in foster parent structure
        paName = self.getNewPANameForStruct(structName)
        self.scriptText += 'pseudoatom %s, pos=%s, label=%s, name=%s \n' %(annStructName,
                                                                           geoCenterStr,
                                                                           labelText,
                                                                           paName)

        for posn in fosterAnn['posns']:
            if seqPos2ResNum[posn] < 0:
                continue
            #distName = self.getNewDistName()
            distName = annStructName+'_lines'
            #distName = structName
            self.scriptText += 'distance %s, /%s////%s, /%s///`%i/CA \n'%(distName,
                                                                          annStructName,
                                                                          paName,
                                                                          structName,
                                                                          seqPos2ResNum[posn])
            self.scriptText += 'hide labels, %s\n' %(distName)
            
    def addAnnObj(self,fosterAnn):
        # Create new object of homologous residues
        seqPos2ResNum = self.GT.getSeqPos2ResNum(structName)
        surfObjName = fosterAnn['parentStruct']+'_'+fosterAnn['title']
        resnumsStr = '+'.join([str(seqPos2ResNum[i]) for i in fosterAnn['posns']])
        self.scriptText += 'cmd.create("%s", "/%s///`%s")\n' %(surfObjName,
                                                               fosterAnn['parentStruct'],
                                                               resnumsStr)
        #self.scriptText += 'cmd.show("surface","%s")\n' %(surfObjName)
        self.scriptText += 'cmd.show("mesh","%s")\n' %(surfObjName)
        ## Color by atom name
        #colorScheme = self.colorSchemes[self.structAnnCounter[fosterAnn['parentStruct']]]
        #self.scriptText += 'util.cba(%s,"%s")\n' %(colorScheme,surfObjName)

        ## Color in solid color
        solidColor = self.solidColors[self.structAnnCounter[fosterAnn['parentStruct']]]
        self.scriptText += 'cmd.color(%s,"%s")\n' %(solidColor,surfObjName)
                                         
    def finalize(self):
        #self.scriptText += "cmd.bg_color('white')\n"
        self.scriptText += 'cmd.show("cartoon"   ,"all")\n'
        self.scriptText += 'cmd.hide("lines"     ,"all")\n'
        ## Enable orthoscopic
        self.scriptText += 'set orthoscopic, on \n'

        #self.scriptText += "preset.pretty('all',_self=cmd)\n"
        return
    
    def write(self, outFilename):
        with open(outFilename,'wb') as of:
            of.write(self.scriptText)


def checkFosterAnn(fosterAnn):
    # Ensure that this matches up with pdb-based sequence
    GT = GeoTools()
    parentStructFile = structDir + fosterAnn['parentStruct'] + '.pdb'
    U = MDA.Universe(parentStructFile)

    seqPos2ResNum = GT.getSeqPos2ResNumFromU(U)
    seqPos2ResLet = GT.getSeqPos2ResLetFromU(U)
    for checkAAPos, checkAALet in zip(fosterAnn['posns'],fosterAnn['checkAA']):
        pdbResLet = seqPos2ResLet[checkAAPos]
        print 'checkAAPos:',checkAAPos, 'checkAALet:', checkAALet, 'pdbResLet:', pdbResLet
    
    # Ensure that this matches up with alignment-based sequence
    
        
if __name__=='__main__':

    allStructures = []
    for structureList in families.values():
        allStructures += structureList
    # Some may be repeated - Filter to get uniques
    allStructures = list(set(allStructures))
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-s',
                        nargs = '*',
                        default = [],
                        help='Individual structures to load. Options are ' + ', '.join(sorted(allStructures)))
    parser.add_argument('-sf',
                        nargs = '*',
                        default = [],
                        help='Structure families to load. Options are ' + ', '.join(sorted(families.keys())))
    parser.add_argument('-a',
                        nargs = '*',
                        default = [],
                        help='Annotations to load. Options are ' + ', '.join(sorted(annotations.keys())))
    args = parser.parse_args()


    structuresToLoad = []
    annotationsToLoad = []
    ## Process single structure arguments
    for sarg in args.s:
        sarg = sarg.replace(',','')
        if not(sarg in allStructures):
            raise Exception('%s is not a known structure' %(sarg))
        structuresToLoad.append(sarg)

    ## Process structure family arguments
    for sfarg in args.sf:
        sfarg = sfarg.replace(',','')
        if not(sfarg in families.keys()):
            raise Exception('%s is not a known structure family' %(sfarg))
        structuresToLoad += families[sfarg]

    ## Process annotations
    for ann in args.a:
        ann = ann.replace(',','')
        if not(ann in annotations.keys()):
            raise Exception('%s is not a known annotation' %(ann))
        annotationsToLoad.append(ann)

    #print 'structuresToLoad', structuresToLoad
    #print 'annotationsToLoad', annotationsToLoad
    
    homA = HomAlignment()

    ## Translate all orig annotation objects to new ("foster") parents
    fosterAnns = []
    GT = GeoTools()
    for stl in structuresToLoad:
        for atl in annotationsToLoad:
            origAnn = annotations[atl]
            fosterAnn = homA.translateAnnotation(origAnn, stl)
            fosterAnn['title'] = atl
            # Check if none of the annotation residues map onto this structure
            if fosterAnn is None:
                continue
            checkFosterAnn(fosterAnn)
            ## Calculate label xyz positions
            geoCenter = GT.getGeoCenter(fosterAnn)
            fosterAnn['geoCenter'] = geoCenter
            ptnGeoCenter = GT.getPtnGeoCenter(fosterAnn)
            regionOffset = geoCenter - ptnGeoCenter
            regionOffUnit = regionOffset / np.linalg.norm(regionOffset)
            regionOff30 = 30. * regionOffUnit
            labelLoc = regionOff30 + ptnGeoCenter
            fosterAnn['labelLoc'] = labelLoc
            fosterAnns.append(fosterAnn)
            
    
    ## Send structure loading and annotation info to pymolWriter 
    myPSW = PymolScriptWriter()
    for structName in structuresToLoad:
        myPSW.loadStruct(structName)
        for fosterAnn in fosterAnns:
            if fosterAnn['parentStruct'] == structName:
                myPSW.addAnn(fosterAnn)
                                                              
    myPSW.finalize()
    myPSW.write('visualize.pml')

    #os.system('pymol  visualize.pml ')
    
'''
dep = Depict()
print dir(dep.alnObj[0])
sampleResnum = 10
fromSeqName = 'a3a'
toSeqName = 'a3b_ctd'
fromSeq = dep.getSeqByName(fromSeqName)
toSeq = dep.getSeqByName(toSeqName)
print fromSeq.id, fromSeq[sampleResnum], sampleResnum
print "print dep.findHomResidue(%i, '%s','%s')" %(sampleResnum, fromSeqName, toSeqName)
print dep.findHomResidue(sampleResnum, fromSeqName,toSeqName)

#from Annotations import *
for annotation in annotations:
    print annotation
    thisA = annotations[annotation]
    print thisA
    seq = dep.getSeqByName(thisA['parent'])
    for AA, posn in zip(thisA['checkAA'], thisA['posns']):
        gappedPosn = indexWOGaps2WGaps(posn, seq)
        print 'Should be', AA, posn, '(gapped', gappedPosn, ')is',  seq[gappedPosn]

#help(dep.alnObj[0])
#print dep.alnObj[:10,10:20]

#for i in range(len(dep.alnObj)):
#    print dep.alnObj[i][0]
#    print dep.alnObj[i].seq[0]
#    print 
'''
