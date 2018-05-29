#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Laura Pedraza-Gonzalez | May 2018

import os
import sys 
import shutil
import glob # The glob module provides a function for making file lists from directory wildcard searches
import urllib # download files from internet
import tempfile
import re
import numpy as np
import commands 
import time
import os.path as path
##################################################
#Software
##################################################
propkaScript = 'python2.7 /Users/laurapedraza/Documents/Software/propka-3.1/scripts/propka31.py'
pdb2pqr= 'python2.7 /Users/laurapedraza/Documents/Software/apbs-pdb2pqr-master/pdb2pqr/pdb2pqr.py'
putIon= '/Users/laurapedraza/Documents/Software/put_ions_good/./ion.x'
##################################################
#References
##################################################
RCSB = "---> This step uses the https://www.rcsb.org/ web page \n (Berman, H. M. Nucleic Acids Res., (2000), 28, 235–242.)"
vmdCite = "---> This step uses the VMD - Visual Molecular Dynamics software \n (Humphrey, W., Dalke, A. and Schulten, K., J. Molec. Graphics, (1996), 14, 33-38.)"
propkaCite = "---> This step uses the PROPKA3 software. \n (Olsson, M. H., Søndergaard, C. R., Rostkowski, M., and Jensen, J. H., JCTC, (2011), 7(2), 525-537.)"
fpocketCite = "---> This step uses the fpocket software. \n (Le Guilloux, V., Schmidtke, P. and Tuffery, P. BMC Bioinformatics, (2009), 10(1), 168)"
modellerCite = "---> This step uses the MODELLER 9.19 software. \n (Webb, B. and Sali, A. Curr. Protoc. Bioinformatics (2014), 47, 5.6.1–32)"
warning = '\x1b[0;33;31m'+"Warning:"'\x1b[0m'
##################################################
#Global lists
##################################################
rmList =sorted([ "ACE", "HG", "HOH", "ZN", "HTG", "HTO", "MAN", "NAG", "BMA", "SO4", "BNG", "A", "CL", "NA", "BOG", "PLM", "TWT", "OLA", "PEE", "GLC", "GAL", "L2P", "L3P"])
aaList=["ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS","MET","PHE",
        "PRO","SER","THR","TRP","TYR","VAL"]
totalList= aaList+rmList
chargedList= ["C-", "Oco", "N+", "ACE"]
protList = ['ASP', 'HIS', 'LYS', 'GLU', 'ARG']+chargedList
protAA = {'ASP': 'ASH', 'LYS':'LYD', 'GLU':'GLH', 'ARG':'ARN'}
pH = ''
pdbPropkaTemp = ''
##################################################
#This function prints the inital message of the script
##################################################
def initialMessage():
    print '\n', '-' * 100
    print '*' * 16, 'Welcome to the interface for the automatic preparation of the ARM input file', '*' * 16 
    print 'Laura Pedraza-González, Luca De Vico, Federico Melaccio, Prof. Massimo Olivucci. \n Laboratory for Computational Photochemistry and Photobiology (LCPP), University of Siena  - 2018'
    print '-' * 100, '\n'
    pass

##################################################
#This function creates a new folder
##################################################
def createWorkingFolder(newpath):
    newpath = pdbName[:-4]+'_ARM_input' 
    if not os.path.exists(newpath):
        os.makedirs(newpath)
        globals().update({ "workingFolder" : newpath+"/"})

    moveFile(pdbName, newpath+"/"+pdbName)
    os.system("cp *seqmut "+newpath+"/")
    os.chdir(newpath)

##################################################
#This function prints the header of each subroutine
##################################################
def Step(step, developer=' '):
    print '_' *100
    print '-' *15, str('\x1b[3;33;40m'+step+'\x1b[0m').ljust(98,'-'), developer
    print '_' *100
    pass

##################################################
# General Functions
##################################################
def copyFile(FileOrigin, FileCopy):
    shutil.copyfile(FileOrigin, FileCopy)
    pass

def moveFile(FileOrigin, FileMoved):
    shutil.move(FileOrigin, FileMoved)
    pass

def replaceString(oldString, newString):
    regex = re.compile("(%s)" % "|".join(map(re.escape, newString.keys())))
    return regex.sub(lambda x: str(newString[x.string[x.start() :x.end()]]), oldString)

##################################################
# This function replaces a character or string with information contained in a dictionary
##################################################
def replaceLine(oldFile, string1, string2, newString, newFile = "TempFile", mvFile = True):
    with open(oldFile, "r") as oldfile, open(newFile, "w") as newfile:
        oldfile_read = oldfile.readlines()
        for line in oldfile_read:
            line_number = oldfile_read.index(line)        
            if string1 in line and string2 in line:
                oldfile_read[line_number] = replaceString(oldfile_read[line_number],newString)
                newfile.writelines(oldfile_read[line_number])
            else:
                newfile.writelines(oldfile_read[line_number])
        if mvFile == True:
            moveFile(newFile, oldFile)

##################################################
# This function deletes the lines with a determined pattern
##################################################
def deleteLine(oldFile, string1, string2, newString,  newFile = "TempFile"):
    with open(oldFile, "r") as oldfile, open(newFile, "w") as newfile:
        oldfile_read = oldfile.readlines()
        for line in oldfile_read:
            line_number = oldfile_read.index(line)        
            if 'ATOM' in line and not (string1 in line and string2 in line.split()[5]):
                newfile.writelines(oldfile_read[line_number])
    moveFile(newFile, oldFile)

##################################################
#This function ask a yes/no question via raw_input() and return their answer
##################################################
def yes_no(question, default= None):

    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)
    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")

##################################################
#This function is used to download pdb files directly from the RCSB Protein Data Bank web page
##################################################
def DownloadPDB():
    url = "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId="
    question = yes_no('\n Do you want to download the PDB file from the RCSB Protein Data Bank webpage?')
    while(1):
        if question == True:
            while(1):
                pdbID = raw_input('<-> Introduce the pdbID (i.e. 1u19): \t').split(' ')
                for i in pdbID:
                    pdbid = url+str(i)
# Verify if the pdb file requested is found on the RSCB server
                    print '\n', str('Verifying that the requested PDB file is on the RSCB server').rjust(100, '.')
                    if "Not Found" in urllib.urlopen(pdbid).read():
                        question2 = yes_no("\n The requested URL /download/"+'\x1b[0;31;43m'+i+".pdb"+'\x1b[0m'+" was not found on the server.\n Do you want to download another PDB file?")
                        if question2 == True:
                            continue
                        else:
                            print "\n ---> No file has been downloaded."
                            question3 = yes_no("Do you want to continue to the next step?")
                            if question3 == True:
                                return
                            else:
                                print "\n OK. Good bye!"
                                sys.exit()
                    else:
# Download the pdb file
                        print '.'*29, 'Downloading the pdb file'
                        open( i+".pdb", "w" ).write( urllib.urlopen(pdbid).read() )
                        print "\n ---> The", '\x1b[0;31;43m' +i+".pdb"+'\x1b[0m', "file is now in your folder. "

                        with open(i+".pdb") as pdbFile:
                            for line in pdbFile:
                                if 'TITLE' in line:
                                    print "<---", line
                                    return
        else:
            print "\n ---> No file has been downloaded. \n"
            question3 = yes_no("Do you want to continue to the next step?")
            if question3 == True:
                return
            else:
                print "\n OK. Good bye!"
                sys.exit()

##################################################
#This function identify the .ext files in the folder and ask the user to choose the correct one
##################################################
def SearchFileType(ext, message0 = "",  message1 = "", message2 = ""):
    extList = glob.glob('*'+ext)
    ChooseNumOption(extList, file, ext, message0,  message1, message2, True)

##################################################
#This function identify the elements of a list, creates an enumerate list, ask the user to choose the number of 
#the correct one and creates a global variable with that information
#########ç#########################################
def ChooseNumOption(nameList, element, type, message0, message1, message2, pick, dictionary = {}): 
# element=(fileName, residue, chain...), type=(pdb, chromophore, chain), message1=(initial message to show the list), message2=(final message)
    if nameList:
        print message0,  '\x1b[0;31;43m'+str(type)+'\x1b[0m', message1
        for i, element in enumerate(nameList, 1):
            if dictionary:
                print str(i)+str(')'), element, dictionary.get(element)
            else:
                print str(i)+str(')'), element
        
        if pick == True:
            number = PickNumber(len(nameList))
            for i, element in enumerate(nameList, 1):
                if i == number:
                    print '\n ---> The', '\x1b[0;31;43m'+str(element)+'\x1b[0m', message2 # This file is the pdbFile
                    globals().update({type+str("Name") : element})
                    return
        else:
            return
#in case of files
    else:
        if element == file:
            print '\n No '+'\x1b[6;30;42m'+ '.'+type +'\x1b[0m', 'files found. Put a ' +'\x1b[6;30;42m'+ '.'+type +'\x1b[0m', 'file in the folder and start again.  \n Good bye! \n'
            sys.exit()
        if element != file and element != "chromophore" and element != "chain":
            print '\n No '+'\x1b[6;30;42m'+ type +'s'+'\x1b[0m', 'found in the ', '\x1b[0;31;43m'+pdbName+'\x1b[0m', "file"
        else:
            print '\n No '+'\x1b[6;30;42m'+ type +'s'+'\x1b[0m', 'found in the ', '\x1b[0;31;43m'+pdbName+'\x1b[0m', "file. The file is corrupt,  verify your PDB and start again! \n Good bye! \n"
            sys.exit()
##################################################
#This function ask for a number via raw_input() and return the number
##################################################            
def PickNumber(lenList, message = ' To select the correct option pick a number in range ',min = 1, typeInput = int):
    while True:
        try:
            input = typeInput(raw_input('\n'+message+str(min)+'-'+str(lenList)+': \t')) 
        except ValueError: 
            print 'That\'s not a number!'
        else:
            if min <= input <= lenList: 
                return input
            else:
                print 'Number out of range. Try again!'
    
##################################################
#This function identify the possible chains in the pdbName file, ask the user to choose the correct one and stores the chain identifier as the global variable chainName
##################################################
def SearchChain(type, column, pick, group = "ATOM"):
    typeList = []    
    listtype = []    
    with open(pdbName) as pdbInitial:
        for line in pdbInitial:
            if group in line.split()[0] and (line.split()[column]).isalpha():
                listtype.append(line.split()[column])
                typeList = sorted(list(set(listtype)))

        ChooseNumOption(typeList,type, type, '\n The following are identified as possible', 'in the '+str(pdbName)+ ' file:', 'is selected as the '+type, pick, "")        

##################################################
#This function generates a new formatted pdb file, with the information of the selected chain
##################################################
def createpdbARM(pdbName, chainName):          
    string0 = ['mol load pdb %s \n' % pdbName,  
               'set sel [atomselect top \"chain '+chainName+ '\" ]\n',
               '$sel writepdb ' +pdbName[:-3]+'ARM.pdb \n']
    VMDTempFile("centerVMD", string0, "center", False)

    global pdbARM # This is the final *ARM.pdb file
    global pdbARMTemp # This is the temporal .pdb file with the columns separated by tabs
    pdbARM = pdbName[:-3]+'ARM.pdb'
    pdbARMTemp = pdbName[:-3]+chainName+'.temp'

    FormatPDB1(pdbARM, pdbARMTemp) # The *ARM.pdb file is re-written with the columns separated by tabs

##################################################
#This function generates a temporal formatted pdb file,in which the columns are separated by tabs. 
#This separation is fundamental to read the content of the columns in the PDB file 
##################################################
def FormatPDB1(pdb, pdbTemp):
    with open(pdb, "r") as oldfile, open(pdbTemp, "w") as newfile:
        for line in oldfile:
            file = line[:+16]+"\t"+line[16:22]+"\t"+line[22:60]+"\t"+line[60:]
            newfile.writelines(file)

##################################################
#This function identify the residues with different rotamers and ask the user for select the correct one
##################################################
def SearchRotamers(type):
    listtype = []
    typeList = []
    with open(pdbARMTemp) as pdbTemp:
        for line in pdbTemp:
            if 'ATOM' in line and float(line.split()[9]) < 1.0:
                listtype.append(line.split()[3])
                typeList = sorted(list(set(listtype)))
    globals().update({type+str("List") : typeList})

    DicResNameNum(" ", type, typeList, totalList, {}, False)                    

    if rotamerList:
        newString = {}
        print "\n Please check the", pdbName, "file to select the rotamer you think better fits your structure. "
        for i in range(0,len(rotamerList)):
            newString.update({ rotamerList[i] : " "+str(rotamerList[i])[+-3:] })
            question = yes_no('\n <-> Do you want to keep the residue ' + str(rotamerList[i]+" "+rotamerDic.get(rotamerList[i])+"?"))
            if question == True:
                replaceLine(pdbARMTemp, rotamerList[i], rotamerDic.get(rotamerList[i]), newString)
                print "---> The "+'\x1b[0;33;49m'+str(rotamerList[i]+" "+rotamerDic.get(rotamerList[i])+'\x1b[0m'+" rotamer will be KEEPED. Its new label is: "+'\x1b[0;33;49m'+newString.get(rotamerList[i])+" "+rotamerDic.get(rotamerList[i])+'\x1b[0m')
                              
            if question == False:
                deleteLine(pdbARMTemp, rotamerList[i],  rotamerDic.get(rotamerList[i]), newString)
                print "---> The "+'\x1b[0;33;49m'+str(rotamerList[i]+" "+rotamerDic.get(rotamerList[i])+'\x1b[0m'+" rotamer will be DELETED")

##################################################
#This function identifies the residue sequence number and creates a dictionary with the residue name and the residue sequence number
##################################################
def DicResNameNum(group, type, typeList, eraseList, typeDictionary, pick):                    
    typeList = sorted(list(set(typeList).difference(eraseList)))
    globals().update({type+str("List") : typeList})
    with open(pdbARMTemp) as pdbTemp:
        for line in pdbTemp:
            for i in range(0, len(typeList)):
                if "ATOM" in line and typeList[i] in line:
                    typeDictionary.update({ typeList[i] : line.split()[5]} )         

    globals().update({type+str("Dic") : typeDictionary})

    ChooseNumOption(typeList,type, type, '\n The following are identified as possible', 'in the '+str(pdbName)+ ' file:', 'is selected as the '+type, pick, typeDictionary)     
##################################################
#This function identifies the possible chromophores and ask the user to select the correct one
##################################################
def SearchChromophore(type, pick, group = 'ATOM'):
    listtype = []
    typeList = []
    with open(pdbARMTemp) as pdbTemp:
        for line in pdbTemp:
            if group in line:
                listtype.append(line.split()[3])
                typeList = sorted(list(set(listtype)))
    globals().update({type+str("List") : typeList})

    DicResNameNum(" ", type, typeList, totalList, {}, pick)                    

##################################################
#This function identifies the linker atom, the linker amino acid and the main counterion
##################################################
def SearchLinkerAA(diffPar = 1.5):
    numAtomsChr = 0
    with open(pdbARMTemp) as file:
        coordXChr = []
        coordYChr = []
        coordZChr = []
        labelsChr = []
        for line in file:
            numAtomsChr = numAtomsChr+1
            if chromophoreName in line:
                coordXChr.append(float(line.split()[6]))
                coordYChr.append(float(line.split()[7]))
                coordZChr.append(float(line.split()[8]))
                labelsChr.append(line.split()[2])
                globals().update({"coordXChr" : coordXChr})
                globals().update({"coordYChr" : coordYChr})
                globals().update({"coordZChr" : coordZChr})
                globals().update({"labelsChr" : labelsChr})

    with open(pdbARMTemp) as file:
        for line in file:
            for i in range(0, len(aaList)):
                if aaList[i] in line and ("N" in line.split()[2] or "O" in line.split()[2]):
                    for j in range(0, len(coordXChr)):                    
                        diffX = float(abs((float(line.split()[6])) - coordXChr[j]))
                        diffY = float(abs((float(line.split()[7])) - coordYChr[j]))
                        diffZ = float(abs((float(line.split()[8])) - coordZChr[j]))
                        if diffX < diffPar and diffY < diffPar and diffZ < diffPar:
                            linker_aa_ID = line.split()[5]
                            linker_aa = str(line.split()[3]+" "+linker_aa_ID)
                            globals().update({"linker_aa_ID" : linker_aa_ID})
                            globals().update({"linker_aa" : linker_aa})

                            coordXYZChr = [coordXChr[j], coordYChr[j], coordZChr[j]]
                            globals().update({"coordXYZChr" : coordXYZChr})


    while(1):
        counterion_ID = raw_input("\n <-> Type the Residue sequence number of the main counterion: \t")
        question3 = yes_no('---> The residue you selected as main counterion is: '+counterion_ID+'\n Are you sure about your selection?' )

        if question3 == True:
            globals().update({"counterion_ID" : counterion_ID})
            with open(pdbARMTemp) as pdbTemp:
                for line in pdbTemp:
                    if "ATOM" in line and chainName+'{:>4}'.format(str(counterion_ID)) in line:
                        main_counterion = str(line.split()[3]+" "+counterion_ID)
                        globals().update({"main_counterion" : main_counterion})
            break

        if question3 == False:
            continue

    print "\n ---> The linker atom of", chromophoreName, "is "+'\x1b[0;33;49m'+labelsChr[j]+'\x1b[0m'+". The linker amino acid is", '\x1b[0;33;49m'+linker_aa+'\x1b[0m and the main counterion is '+'\x1b[0;33;49m'+main_counterion+'\x1b[0m'

#############################
def VMDTempFile(tempName, content, position, default=False):
    with tempfile.NamedTemporaryFile() as tempNameI:
        tempNameI.writelines(content) 
        tempNameInp = tempNameI.name
        tempNameI.seek(0)

        with tempfile.NamedTemporaryFile() as tempNameO:
            tempNameOut = tempNameO.name
            os.system('vmd -dispdev text -eofexit <' + tempNameInp + '> ' + tempNameOut )

            if default == True:
                for line in tempNameO:
                    if  '{' in line:
                        line = line.replace("{", "").replace("}", "").replace("\n", "")
                        positionXYZ = line.split(" ")
                        for i in range(0,len(positionXYZ)):
                            positionXYZ[i] = float(positionXYZ[i]) * -1 # xyz coordinates
                            globals().update({position+str("XYZ") : positionXYZ}) 
                        
                return positionXYZ

##################################################
#This function generates a PDB file, using the correct format, with the information of chromophore, amino acids and waters
##################################################
def FormatPDB(oldFile, newFile, FinalFile):
    global resNumChr # New assigned residue sequence number of the chromophore

    pdbFormat = "%-6s%5d %s %-4s%3s %s%4d %3s%8.3f%8.3f%8.3f%6.2f%s%6.2f %11s"
    with open(oldFile, "r") as oldfile, open(newFile, "w") as newfile:
        oldfile_read = oldfile.readlines()
        ln = 0 # Atom serial number
        
        for line in oldfile_read:    
            ls = line.split()
            if "ACE" in line:
                ln = ln+1 # Atom serial number
                line = pdbFormat % ("HETATM", int(ln), str(""), ls[2], ls[3],ls[4], int(ls[5]), str(""), float(ls[6]), 
                                    float(ls[7]), float(ls[8]), float(ls[9]),str(""), float(ls[10]), ls[11])
                newfile.writelines(line+"\n")

        for line in oldfile_read:
            for i in range(0, len(aaList)):
                ls = line.split()
                if aaList[i] in line:
                    ln = ln+1
                    line = pdbFormat  % ("ATOM", int(ln), str(""), ls[2], ls[3],ls[4], int(ls[5]), str(""), float(ls[6]), 
                                         float(ls[7]), float(ls[8]), float(ls[9]),str(""),float(ls[10]), ls[11])
                    newfile.writelines(line+"\n")
                    resNum = int(ls[5]) # Residue sequence number
                    globals().update({"resNumID" : resNum}) # Last residue sequence number used of amino acids

        for line in oldfile_read:    
            ls = line.split()
            if "NHH" in line:
                ln = ln+1 # Atom serial number
                line = pdbFormat % ("HETATM", int(ln), str(""), ls[2], ls[3],ls[4], int(ls[5]), str(""), float(ls[6]), 
                                    float(ls[7]), float(ls[8]), float(ls[9]),str(""), float(ls[10]), ls[11])
                newfile.writelines(line+"\n")
                resNum = int(ls[5]) # Residue sequence number
                globals().update({"resNumID" : resNum}) # Last residue sequence number used of amino acids


        newfile.write("TER \n")
        
        for line in oldfile_read:    
            ls = line.split()
            if chromophoreName in line:
                resNumChr = resNum+1
                ln = ln+1
                line = pdbFormat % ("HETATM", int(ln), str(""), ls[2], ls[3],ls[4], int(resNumChr), str(""), float(ls[6]), 
                                    float(ls[7]), float(ls[8]), float(ls[9]),str(""),float(ls[10]), ls[11])
                newfile.writelines(line+"\n")
        newfile.write("TER \n")
        
        resNumHOH = resNumChr # Residue sequence number
        for line in oldfile_read:
            ls = line.split()
            if "HOH" in line:
                ln = ln+1
                resNumHOH = resNumHOH+1
                line = pdbFormat % ("HETATM", int(ln), str(""), ls[2], ls[3],ls[4], int(resNumHOH), str(""), float(ls[6]), 
                                    float(ls[7]), float(ls[8]), float(ls[9]),str(""),float(ls[10]), ls[11])
                newfile.writelines(line+"\n")
        newfile.write("TER \n")
        os.remove(oldFile)
        moveFile(newFile, FinalFile)
    
    globals().update({"counter" : resNumHOH}) # Last residue sequence number used. Employed later for the numeration of the counterions
    globals().update({"lN" : ln}) # Last atom serial number used. Employed later for the numeration of the counterions

##################################################
#This function align the protein along its main rotational axis; center on the center of mass of the chromophore
#This script will center the whole protein+retinal on the center of mass of the whole complex  
##################################################
def alignRotAxis(chainName, pdbName):
    string1 = ['mol load pdb %s \n' % pdbARM,  
               'molinfo 0 get center\n']
    print "\n", str("Centering the protein+"+chromophoreName+" on the center of mass of the complex").rjust(100, ".")
    VMDTempFile("centerVMD", string1, "center", True)

    string2 = ['mol load pdb %s \n' % pdbARM,
               'set sel [atomselect 0 \"all\"]\n',
               'atomselect0 moveby {%.6f %.6f %.6f}\n' %(float(centerXYZ[0]), float(centerXYZ[1]), float(centerXYZ[2])),
               'package require Orient\n',
               'namespace import Orient::orient\n',
               'set sel [atomselect top \"all\"]\n',
               'set I [draw principalaxes $sel]\n',
               'set A [orient $sel [lindex $I 2] {0 0 1}]\n',
               '$sel move $A\n',
               'set I [draw principalaxes $sel]\n',
               'set A [orient $sel [lindex $I 1] {0 1 0}]\n',
               '$sel move $A\n',
               'set I [draw principalaxes $sel]\n',
               '$sel writepdb %s \n' % pdbARM]
    VMDTempFile("orientVMD", string2, "orient", False)

    string3 = ['mol load pdb %s \n' % pdbARM,
               'set sel [atomselect top \"resname '+chromophoreName+'\"]\n',
               '$sel writepdb ' +chromophoreName+"."+chainName+'.pdb \n',
               'mol delete top\n',
               'mol load pdb ' +chromophoreName+"."+chainName+'.pdb \n',
               'molinfo top get center\n' ]
# create file with chromophore coordinates
    VMDTempFile("newcenterVMD", string3, "newcenter", True) 

#center the protein on the retinal center of mass
    string4 = ['mol load pdb %s \n' % pdbARM,
               'set sel [atomselect 0 \"all\"]\n',
               'atomselect0 moveby {%.6f  %.6f  %.6f}\n' %( float(newcenterXYZ[0]), float(newcenterXYZ[1]), float(newcenterXYZ[2])),
               'molinfo top get center\n',
               '$sel writepdb %s\n' % pdbARM] #VMD saves the new pdb file                                                                                               
    VMDTempFile("orientVMD", string4, "orient", True)
    print str("Done!").rjust(100, ".")

    os.remove(chromophoreName+"."+chainName+'.pdb')
    copyFile(pdbARM, pdbARMTemp)
    print "\n ---> The chromophore label has been changed to: ", '\x1b[0;33;49m'+str(chromophoreName)+" "+str(resNumChr)+'\x1b[0m'
#identify the linker atom and the linker amino acid
    SearchLinkerAA()
    print "\n ----> The", '\x1b[0;33;49m'+str(pdbARM)+'\x1b[0m', "ARM input file has been generated."

##################################################
# Luca's proposal to calculate the charge
##################################################
def Charge(Name, NameNum, pKa, pH):
    global resPkaAna
    global The_Charge
    The_Charge = 0.0
    
    if pKa == 99.99: # CYS residue bonded in di-sulfide bond, no charge
        return The_Charge
    #Calculation of the charge based on the calculated pKa and the given pH
    exponent = np.power(10, (pKa - pH))
    The_Charge = exponent / (1.0 + exponent)

    if (Name == 'ASP') or (Name == 'GLU') or (Name == 'C-') or (Name == 'CYS') or (Name == 'TYR') or (Name == 'Oco') or (Name == 'OXT'):
        The_Charge = The_Charge-1

    The_Charge = round(The_Charge)
    resPkaAna= "%3s%4d%5d" % (Name, int(NameNum), The_Charge)
    return The_Charge, resPkaAna    

##################################################
#This function perform the propka analysis and suggest the user a list with possible amino acids which should be protonated
##################################################
def propKa():
    resPkaAnaList = []
    question = True
#Uncomment the next line to ask the user for run the propka analysis
#    question = yes_no('\n <-> Do you want to perform the pKa analysis?')
    if question == True:

        global pdbPropkaTemp, pdbPropka 
        pdbARMFix = pdbARM[:-3]+'fix.pdb'
        pdbPropka = pdbARMFix[:-3]+'pka'
        pdbPropkaTemp = pdbPropka+'.temp'

        os.system(pdb2pqr+" --chain --ff=amber "+str(pdbARM+" "+pdbARMFix))

        print '\n', str('Running PROPKA3.0 analysis for the \x1b[0;33;49m'+str(pdbARM)+'\x1b[0m input file').rjust(100, '.')
        os.system (propkaScript+" "+pdbARMFix+ ">> /dev/null")
        print str("Done! The files \x1b[0;33;49m"+pdbARM[:-3]+"propka_input\x1b[0m and \x1b[0;33;49m"+pdbPropka+"\x1b[0m has been generated. ").rjust(100, '.')

        pH = PickNumber(14.0, '<-> Please write the pH-value (suggested value physiological pH 7.4) in the range ',  0, float)
        globals().update({ "pH"  : str(pH)})

        # A temporal file with the summary of the propka analysis
        with open(pdbPropka) as file, open(pdbPropkaTemp, "w") as file2:
            content = file.read()
            text = re.search(r'SUMMARY OF THIS PREDICTION\n.*?--------------------------------------', content, re.DOTALL).group()
            file2.write(text)

        # Analysis based on Luca's proposal
        print "\n ---> At pH ", pH, "the predicted charge of the residues is: \n", '_'*40 +'\n', '{:^9}'.format('RESIDUE')+'{:^6}'.format('CHARGE')+'{:^6}'.format(' pKa')+'{:^18}'.format(' (pKa - pKa-model)')+ '\n'+ '_'*40

        with open(pdbPropkaTemp) as pkaFile:
            for line in pkaFile:
                for i in range(0, len(protList)):
                    if protList[i] in line.split()[0] and linker_aa_ID not in line.split()[1] and counterion_ID not in line.split()[1]:
                        shift = 2.0
                        pKa_calc = line.split()[3]
                        pKa_model = line.split()[4]
                        diff_pKa = '%.2f' %(abs((float(pKa_calc) - float(pKa_model))))
                        Charge(line.split()[0], line.split()[1], float(pKa_calc), pH)
                        print resPkaAna, '{:^11}'.format(pKa_calc) , '{:^10}'.format(diff_pKa)
                        if (protList[i] == "ASP" or protList[i] == "GLU") and The_Charge != -1:
                            resPkaAnaList.append(resPkaAna)
                        
                        if (protList[i] == "ARG" or protList[i] == "LYS") and The_Charge != 1:
                            resPkaAnaList.append(resPkaAna)
                        # Normal analysis using a shift paramter
                        if protList[i] == "HIS" and float(diff_pKa) > shift:
                            resPkaAnaList.append(resPkaAna)

        ChooseNumOption(resPkaAnaList,"", "ionizable residue", '\n Based on the computed charges, the suggested residues to be protonated are:', '', '', False, "")     
        print warning, " Check carefully if HIS residues must be protonated!"        
        print warning, "The linker amino acid", linker_aa, "and the main counterion", main_counterion, "have been excluded from this analysis"

    else:
        print '---> The pKa analysis will not be performed'

##################################################
#This function perform the change of ionization states for  selected amino acids
##################################################
def protonation(pdbARM):
    protResDictionary = {}
    question = yes_no('\n <-> Based on the computed charges and the experimental information, do you want to change \n the ionization state of any amino acid? (ASP --> ASH, GLU --> GLH, HIS --> (HID or HIE or HIP), LYS --> LYD)')
    if question == True:
        while(1):
            protRes = raw_input('\n <-> Introduce the NUMBER ID of the residues, separated by space \" \" (i.e. 1 2 3): \t').split(' ')
            with open(pdbARM) as pdbTemp:
                for line in pdbTemp:
                    for i in range(0, len(protRes)):
                        if "ATOM" in line and chainName+'{:>4}'.format(protRes[i]) in line:
                            protResDictionary.update({ protRes[i] : line.split()[3] } )
   
            question2 = yes_no('---> The residues you selected are: '+str(protResDictionary)+'\n Are you sure about your selection?' )
            print "\n"
            if question2 == True:

                globals().update({"protResDictionary": protResDictionary })                
                print "*"*80, "\n", warning, "Remember that:", "\n HID: Histidine with hydrogen on the delta nitrogen", "\n HIE: Histidine with hydrogen on the epsilon nitrogen", "\n HIP: Histidine with hydrogens on both nitrogens; this is positively charged. \n", "*"*80, "\n"

                for key in protResDictionary:
                    if protResDictionary[key] == "HIS":
                        print "\n <-> Please select the correct ionization state of the residue \x1b[0;33;49m"+protResDictionary[key]+" "+key+"\x1b[0m:"
                        HISList = ["HIE", "HID", "HIP"]
                        ChooseNumOption(HISList, "HIS", "HIS", '\n Choose the ', 'protonation state.', 'is the new protonation state.', True)
                        replaceLine(pdbARM, "ATOM", chainName+'{:>4}'.format(key), {"HIS" : HISName})
                    
                    else:
                        replaceLine(pdbARM, "ATOM", chainName+'{:>4}'.format(key), protAA)
                        print "\n ---> The protonation state of the residue \x1b[0;33;49m"+protResDictionary[key]+" "+key+"\x1b[0m has been changed to \x1b[0;33;49m"+str(protAA[protResDictionary[key]])+"\x1b[0m \n"
                return
                
            if question2 == False:
                continue
    else:
        print '---> No amino acid requires to change its ionization state'

##################################################
#This function calculates the number of counterions needed to neutralize the system
##################################################
def numberCounterions(pdbARM):
    question = True
#    question = yes_no('\n <-> Do you want to perform the neutralization analysis? ')
    if question == True:
        TotalTop = 0
        TotalBottom = 0
        print "-"*29, "\n Summary \n", "-"*29

        signlinker_aa = '' # Define the outer and inner side of the protein 
        with open(pdbARM) as file:
            for line in file:
                if "ATOM" in line and chainName+'{:>4}'.format(str(linker_aa_ID)) in line and "CA" in line:
                    signlinker_aa = float(line.split()[8])

        def countPosNeg(aa1, aa2, aa3, Charge, List):
            global topCharge, bottomCharge
            topList = []
            bottomList = []
            topCharge = 0
            bottomCharge = 0

#09-05-2018 The count is based on the position of the LA    
#Position of the linker amino acid                                              
            with open(pdbARM) as file:
                for line in file:
#                    if (aa1 in line or aa2 in line or aa3 in line) and "CA" in line and linker_aa_ID not in line.split()[5] and counterion_ID not in line.split()[5]:
                    if (aa1 in line or aa2 in line or aa3 in line) and "CA" in line: #The linker_aa and main_counterion should be conserved in cases where one of them is protonated 
                        if signlinker_aa > 0:
                            resta = (float(line.split()[8]) - float(signlinker_aa))
                            if  (float(line.split()[8]) - float(signlinker_aa)) > 0:
                                topCharge = topCharge + 1
                                topList.append(line.split()[5])
                                print line.split()[3]+"\t"+line.split()[5]+"\t"+"TOP"+"\t"+str(resta)
                            else:
                                print line.split()[3]+"\t"+line.split()[5]+"\t"+"BOTTOM"+"\t"+str(resta)
                                bottomCharge = bottomCharge + 1
                                bottomList.append(line.split()[5])
                        if signlinker_aa < 0:
                            resta = (float(line.split()[8]) - float(signlinker_aa))
                            if  (float(line.split()[8]) - float(signlinker_aa)) < 0:
                                topCharge = topCharge + 1
                                topList.append(line.split()[5])
                                print line.split()[3]+"\t"+line.split()[5]+"\t"+"TOP"+"\t"+str(resta)
                            else:
                                print line.split()[3]+"\t"+line.split()[5]+"\t"+"BOTTOM"+"\t"+str(resta)
                                bottomCharge = bottomCharge + 1
                                bottomList.append(line.split()[5])

            globals().update({ "top"+Charge  : topCharge})
            globals().update({ "top"+List  : topList})
            globals().update({ "bottom"+Charge  : bottomCharge})
            globals().update({ "bottom"+List  : bottomList})

            print "_"*29, "\n top"+Charge+"=", topCharge
            print "bottom"+Charge+"=", bottomCharge, "\n"

        countPosNeg("ARG", "HIS", "LYS", "Positive", "PosList")
        countPosNeg("ASP", "GLU", "ASP", "Negative", "NegList")

        TotalTop = topPositive + (topNegative*-1)
        TotalBottom = bottomPositive + (bottomNegative*-1)
        TotalPos = topPositive + bottomPositive
        TotalNeg = topNegative + bottomNegative
        
        globals().update({ "TotalTop"  : TotalTop})
        globals().update({ "TotalBottom"  : TotalBottom})
        globals().update({ "TotalPos"  : TotalPos})
        globals().update({ "TotalNeg"  : TotalNeg})

        print "Number of positively charged residues: \t", TotalPos
        print "Number of negatively charged residues: \t", TotalNeg
 
        def IonTopBottom(TotalPosition, Position):
            if TotalPosition > 0:
                IonPosition = "CL"
            else:
                IonPosition = "NA"
            globals().update({ "Ion"+Position : IonPosition})
        
        IonTopBottom(TotalTop, "Top")
        IonTopBottom(TotalBottom, "Bottom")
        
        print "-"*29, "\n", "Total charge top: \t", TotalTop, "\t", "|  Suggestion: Add to the top ", abs(TotalTop), IonTop
        print "-"*29, "\n", "Total charge bottom: \t", TotalBottom, "\t",  "|  Suggestion: Add to the bottom", abs(TotalBottom), IonBottom
        print "-"*29
#        print "\n", warning, "The linker amino acid", linker_aa, "and the main counterion", main_counterion, "have been excluded from this analysis"

#Exclude the linker amino acid and the main counterion from the list of target residues

        Step('10. Add the counterions to the inner and outer layer of the proteins', '\n Luca De Vico & Laura Pedraza-González.' )
        addCounterIons(pdbARM)

    else:
        print '---> No counterions will be added'
    
##################################################
#This function adds the counterions needed to neutralice the system
##################################################
def addCounterIons(pdbARM):
    global pH

    if pH:
        pass
    else:
        pH = PickNumber(14.0, '<-> Please write the pH-value (suggested value physiological pH 7.4) in the range ',  0, float)
        globals().update({ "pH"  : str(pH)})

#Uncomment the following two lines and comment the third one to ask the user for the force field to generate the pir file        
#    ForceFieldList = ["amber", "charmm", "parse"]
#    ChooseNumOption(ForceFieldList, "FF", "ForceField", '\n Choose the ', 'to generate the pqr file:', 'is selected as the Force Field to run the PDB2PQR software.', True)
    ForceFieldName = "amber" #Force field by default : amber
        
    global pqrARM
    pqrARM = pdbARM[:-3]+"pqr"

    os.system(pdb2pqr+" --ff="+ForceFieldName+" --with-ph="+str(pH+" "+pdbARM+" "+pqrARM))
    print "The ","\x1b[0;33;49m"+pqrARM+"\x1b[0m", "file has been generated and is ready to be used in the putIon module."

    def PutIon():
        global outIon
        global infoIonFile
        infoIonFile = "infoIon"+pdbARM[:-4]+".sh"
#inner
        with open(infoIonFile, "w") as infoIon:
            infoIon.write("#!/bin/bash \n \n")
            infoIon.write(putIon+" << EOF \n")
            infoIon.write(pqrARM+"\n")
            infoIon.write(str(TotalTop)+"\n")
            if TotalTop > 0:
                infoIon.write(str(topPositive)+"\n")
                for i in range(0,len(topPosList)):
                    infoIon.write(str(topPosList[i])+"\n")
            else:
                infoIon.write(str(topNegative)+"\n")
                for i in range(0,len(topNegList)):
                    infoIon.write(str(topNegList[i])+"\n")
#outer
            infoIon.write(str(TotalBottom)+"\n")
            if TotalBottom > 0:
                infoIon.write(str(bottomPositive)+"\n")
                for i in range(0,len(bottomPosList)):
                    infoIon.write(str(bottomPosList[i])+"\n")
            else:
                infoIon.write(str(bottomNegative)+"\n")
                for i in range(0,len(bottomNegList)):
                    infoIon.write(str(bottomNegList[i])+"\n")
            infoIon.write("EOF")

            outIon = "outputPutIon."+pqrARM[:-4]
            os.system("chmod 777 "+infoIonFile)

    PutIon()
#Execution of the PutIon module
    print ("\n ---> Running PutIon analysis over the \x1b[0;33;49m"+str(pqrARM)+"\x1b[0m file for the placement of the counterions").rjust(100, '.')
    def ExecutePutIon():
        os.system("sh "+infoIonFile+" > "+outIon )
    ExecutePutIon()

#########################                                                                                                                                                                                                                     
    def idenfifyError():
        global topPositive, topNegative, bottomPositive, bottomNegative, ErrorION
        with open(outIon, "r") as outputIon:
            for line in outputIon:
                if ("Error: Charge too small, exiting") in line:
                    ErrorION = True
                    outputIon.seek(0)
                    resError = (outputIon.readlines()[-2]).split()[0]
                    if resError in topPosList:
                        topPosList.remove(resError)
                        topPositive = int(len(topPosList))
                    if resError in topNegList:
                        topNegList.remove(resError)
                        topNegative = int(len(topNegList))
                    if resError in bottomPosList:
                        bottomPosList.remove(resError)
                        bottomPositive = int(len(bottomPosList))
                    if resError in bottomNegList:
                        bottomNegList.remove(resError)
                        bottomNegative = int(len(bottomNegList))
                else:
                    ErrorION = False
        PutIon()
        ExecutePutIon()
    idenfifyError()

    global ErrorION

    while ErrorION == True:
        idenfifyError()

####################
    def writeIonCoord():
        print "---> The following coordinates for the ", IonTop, "and", IonBottom, " counterions have been added to the ", pdbARM, "file: \n"
        with open("ions.pdb", "r") as ions:
            contador = counter+1
            line_num = lN+1
            atomName = []
            resName = []
            x_position = []
            y_position = []
            z_position = []
            for line in ions:
                if ("HETATM") in line:
                    atomName.append(line.split()[2])
                    resName.append(line.split()[3])
                    x_position.append(float(line.split()[6]))
                    y_position.append(float(line.split()[7]))
                    z_position.append(float(line.split()[8]))

        with open(pdbARM, "a") as file2:
            for i in range(0,abs(len(atomName))): 
                cont = contador+i  # continue line number in pdbARM
                line_Num = line_num+i       # continue residue number in pdbARM
                pdbFormat = "%-6s%5d %s %-4s%3s %s%4d %3s%8.3f%8.3f%8.3f"
                IonCoordinates = pdbFormat % ("HETATM", line_Num, "", atomName[i], resName[i], chainName, cont, "", x_position[i], y_position[i], z_position[i])
                print pdbFormat % ("HETATM", line_Num, "", atomName[i], resName[i], chainName, cont, "", x_position[i], y_position[i], z_position[i])
                file2.writelines(IonCoordinates+"\n")
            file2.writelines("END")       
    writeIonCoord()

##################################################
#This function calculates the chromophore cavity using fpocket
##################################################
def fpocket():

    question = yes_no('\n <-> Do you want to calculate the chromophore cavity? ')
    if question == True:
        os.system('fpocket -f '+pdbARM)
        print "\n The folder", pdbARM[:-4]+"_out", "containing the pockets has been generated"

#LMPG 29-05-2018
#This function identifies the pocket which contains the linker AA

        pocketFile= commands.getoutput("grep -lr "+"\""+chainName+'{:>4}'.format(str(linker_aa_ID))+"\""+" "+pdbARM[:-4]+"_out/pockets/*pdb")
        pocketFile= pocketFile.split("\n")[0]

        with open(pocketFile) as cavityFile:
            listcavity = [counterion_ID, linker_aa_ID]
            for line in cavityFile:
                if "ATOM" in line.split()[0]:
                    listcavity.append(line.split()[5])                    
                    cavityList0 = list(set(listcavity))
                    cavityList = sorted(cavityList0, key = int)

        with open("cavity", "w") as cavity:
            for i in range(0, len(cavityList)):
                cavity.writelines(cavityList[i]+"\n")

        print '---> The cavity file has been generated'
    
    else:
        print '---> No cavity file generated'


##################################################
#MUTATIONS
##################################################
def Mutations():

    question = yes_no('\n <-> Do you want to perform mutations of the wild type '+pdbARM+' file?')
    if question == True:

        Mutations_protocol()
    else:
        print "---> No mutations requested."

def Mutations_protocol():

#This step ask the user for select the software for the mutations 
    from mutate import Select_mutation_Software
    Select_mutation_Software(ChooseNumOption, pdbARM, copyFile)

    def Mutations_procedure():
#This step ask the user for the list of mutations. Similar to the seqmut file
        from mutate import Insert_mutations
        Insert_mutations(yes_no, warning, workingFolder, pdbARM, copyFile, mutation_seqmut)
#MODELLER routine
        if mutation_SoftwareName == "Modeller":
            from mutate import Modeller_mutations
            Modeller_mutations(chainName, pdbARM, copyFile, moveFile, FormatPDB1, FormatPDB)
#SCWRL4 routine
        if mutation_SoftwareName == "Scwrl4":
            from mutate import Scwrl4_mutations
            Scwrl4_mutations(resNumID, chainName, FormatPDB1, moveFile, pdbARM, FormatPDB)
# This is the subroutine for change the ionization state and add counterions to the mutations        
        from mutate import protonation_mutant
        protonation_mutant(protonation, numberCounterions, protAA, HISName, protResDictionary, replaceLine, chainName, Step)        

    def Mutations_ID():
        SearchFileType('seqmut', '\n The following', 'files with a list of mutations are found: \n', 'file will be used for the mutations..')  #seqmutName

        mutation_seqmutDic={}
        with open(seqmutName, "r") as file:
            content = file.read().splitlines()
            for line in content:
                if "Mutation" in line:
                    key=line
                    mutation_seqmutDic[key]=[]
                else:
                    mutation_seqmutDic[key].append(line)
        globals().update({"mutation_seqmutDic" : mutation_seqmutDic})

    global wt_pdbARM
    wt_pdbARM = "wt_"+pdbARM
    copyFile(pdbARM, wt_pdbARM)
    globals().update({"wt_pdbARM" : wt_pdbARM})

    Mutations_ID()
# This is the subroutine for the mutations. The mutations are readed from a .seqmut file.      
    for key in mutation_seqmutDic:
        mutation_seqmut = mutation_seqmutDic[key]
        Mutations_procedure()

# This function generates a .pir file for the wild type 
#        from mutate import get_pir_Script
#        get_pir_Script(pdbARM, resNumID, chainName, FormatPDB1, moveFile)
##################################################
#Wild type analysis
##################################################
def wT_analysis():
    wt_workingFolder = "wt_"+workingFolder
    os.system("mkdir "+wt_workingFolder)
    os.system("cp "+pdbARM+" "+wt_workingFolder)
    os.chdir(wt_workingFolder)

##################################################
# CALL FUNCTIONS
##################################################

localtime = time.asctime( time.localtime(time.time()) )
print "\n Local current time :", localtime
t0 = time.time()

initialMessage()

Step('1. Download a crystal structure -PDB file- from the RCSB Protein Data Bank \n',  RCSB)
DownloadPDB()

Step('2. Selection of the PDB file')
SearchFileType('pdb', '\n The following', 'files are found: \n', 'file will be used for preparing the ARM input file.') 
createWorkingFolder("workingFolder")

Step('3. Search and selection of the Chain')
SearchChain('chain', 4, True)
createpdbARM(pdbName, chainName)          

Step('4. Find the residue(s) with different rotamers and select the correct one(s)')
SearchRotamers("rotamer")

Step('5. Search and selection of the Chromophore')
SearchChromophore("chromophore", True)

Step('6. Center the whole protein+chromophore on the center of mass of the whole complex', '\n Luca De Vico & Laura Pedraza-González. \n' + vmdCite)
FormatPDB(pdbARMTemp, pdbARM, pdbARM)
alignRotAxis(chainName, pdbName)
FormatPDB1( pdbARM, pdbARMTemp )
FormatPDB( pdbARMTemp, pdbARM, pdbARM)
wT_analysis()

Step('7. Analysis of the amino acid ionization states:\n pKa analysis using the "PROPKA3.1: A PROTEIN PKA PREDICTOR" software', '\n Luca De Vico & Laura Pedraza-González. \n' + propkaCite)
propKa() 

Step('8. Change the ionization state of selected amino acids', '')
protonation(pdbARM)

Step('9. Calculation of the number of counterions (Cl- and Na+) needed to neutralize the system', '')
numberCounterions(pdbARM)

Step('11. Definiton of the '+chromophoreName+' chromophore cavity using fpocket','\n'+ fpocketCite)
fpocket()
os.chdir("../")


Step('12. Preparation of mutants', '\n'+ modellerCite)
Mutations()

Step('The '+'\x1b[0;33;49m'+pdbARM+' file is ready to be used as input for the ARM protocol! \x1b[0m')

if path.exists(pdbPropkaTemp):
    os.remove(pdbPropkaTemp)

print "Total excecution time (min):", (time.time() - t0)/60


