#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Laura Pedraza-Gonzalez | March 2018

import os
import sys 
import glob # The glob module provides a function for making file lists from directory wildcard searches
import tempfile
import re
import numpy as np
import commands 
import time
import textwrap

modellerScript = "mod9.19"
scrl4Script = "Scwrl4"

threeToOneDic = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                 'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
                 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
                 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'UNK': '*'}
oneToTreeDic = {}

for key in threeToOneDic:
    oneToTreeDic.update({threeToOneDic[key] : key})

#############################
#This step generates the getpir.py script. This script is then executed to obtain the .pir file.
#############################
def get_pir_Script(pdbFile, resNumID, chainName, FormatPDB1, moveFile):

    global sequenceWT
    pdbFileTemp = pdbFile[:-3]+"temp"
    pirFile = pdbFile[:-3]+"pir"
    pirFileTemp = pirFile+".temp"

    getpir = ["import os \n", 
              "import sys \n \n",
              "from modeller import * \n \n",
              "env = environ() \n",
              "aln = alignment(env) \n"
              "mdl = model(env, file='"+pdbFile+"') \n",
              "aln.append_model (mdl, align_codes='wt_"+pdbFile+"',atom_files='"+pdbFile+"') \n",
              "aln.write(file='"+pirFile+"', alignment_format='PIR') \n" ]

    getpirScript = "getpyr.py"
    with open(getpirScript, "w") as getpirFile:
        getpirFile.writelines(getpir)
    os.system(modellerScript +" "+getpirScript  )

#Identify missing residues and write the pir file in correct format
    FormatPDB1(pdbFile, pdbFileTemp)
    realResNumDic = {}
    with open(pdbFileTemp, "r") as file:
        for line in file:
            if "ATOM" in line:
                realResNumDic.update({int(line.split()[5]) : threeToOneDic[str(line.split()[3])]})
        globals().update({ "realResNumDic": realResNumDic})

    sequence = ''
    sequenceWT = ''
    sequenceWTList = []
    with open(pirFile, "r") as pir, open(pirFileTemp, "w") as temp:
        for line in pir:
            if str(pdbFile) in line:
                temp.writelines(line)
        missResList = []
        for i in range(0,resNumID):
            i = i+1
            if i not in realResNumDic:
                missResList.append(i)
                sequence = sequence+"-"
                sequenceWTList.append("-")
            else:
                res = str(realResNumDic[i]).lower()
                sequence = sequence+res
                sequenceWT = sequenceWT+res
                sequenceWTList.append(res)
        temp.writelines('\n'.join(textwrap.wrap(sequence, width=75)))
        temp.writelines('* \n')
    globals().update({ "sequenceWTList": sequenceWTList})

    moveFile(pirFileTemp, pirFile)
    
    print "\n The file "+pirFile+" has been generated using the MODELLER 9.19 software. The missing residues "+str(missResList)+" were considered."

    os.remove(pdbFileTemp)
#############################
#This step generates the mutate_model.py script. This script is then executed to *****
#############################

def mutate_model_Script():
    global mutate_modelScript
    mutate_modelScript = "mutate_model.py"
    
    mutate_model = ["import sys \n",
                    "import os \n",
                    " \n",
                    "from modeller import * \n",
                    "from modeller.optimizers import molecular_dynamics, conjugate_gradients \n",
                    "from modeller.automodel import autosched \n",
                    " \n",
                    "# \n",
                    "#  mutate_model.py \n",
                    "# \n",
                    "#     Usage:   python mutate_model.py modelname respos resname chain > logfile \n",
                    "# \n",
                    "#     Example: python mutate_model.py 1t29 1699 LEU A > 1t29.log \n",
                    "# \n",
                    "# \n",
                    "#  Creates a single in silico point mutation to sidechain type and at residue position \n",
                    "#  input by the user, in the structure whose file is modelname.pdb \n",
                    "#  The conformation of the mutant sidechain is optimized by conjugate gradient and \n",
                    "#  refined using some MD. \n",
                    "# \n",
                    "#  Note: if the model has no chain identifier, specify "" for the chain argument. \n",
                    "# \n",
                    " \n",
                    " \n",
                    "def optimize(atmsel, sched): \n",
                    "    #conjugate gradient \n",
                    "    for step in sched: \n",
                    "        step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001) \n",
                    "    #md \n",
                    "    refine(atmsel) \n",
                    "    cg = conjugate_gradients() \n",
                    "    cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001) \n",
                    " \n",
                    " \n",
                    "#molecular dynamics \n",
                    "def refine(atmsel): \n",
                    "    # at T=1000, max_atom_shift for 4fs is cca 0.15 A. \n",
                    "    md = molecular_dynamics(cap_atom_shift=0.39, md_time_step=4.0, \n",
                    "                            md_return='FINAL') \n",
                    "    init_vel = True \n",
                    "    for (its, equil, temps) in ((200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)), \n",
                    "                                (200, 600, \n",
                    "                                 (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0))): \n",
                    "        for temp in temps: \n",
                    "            md.optimize(atmsel, init_velocities=init_vel, temperature=temp, \n",
                    "                         max_iterations=its, equilibrate=equil) \n",
                    "            init_vel = False \n",
                    " \n",
                    " \n",
                    "#use homologs and dihedral library for dihedral angle restraints \n",
                    "def make_restraints(mdl1, aln): \n",
                    "   rsr = mdl1.restraints \n",
                    "   rsr.clear() \n",
                    "   s = selection(mdl1) \n",
                    "   for typ in ('stereo', 'phi-psi_binormal'): \n",
                    "       rsr.make(s, restraint_type=typ, aln=aln, spline_on_site=True) \n",
                    "   for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'): \n",
                    "       rsr.make(s, restraint_type=typ+'_dihedral', spline_range=4.0, \n",
                    "                spline_dx=0.3, spline_min_points = 5, aln=aln, \n",
                    "                spline_on_site=True) \n",
                    " \n",
                    "#first argument \n",
                    "modelname, respos, restyp, chain, = sys.argv[1:] \n",
                    " \n",
                    " \n",
                    "log.verbose() \n",
                    " \n",
                    "# Set a different value for rand_seed to get a different final model \n",
                    "env = environ(rand_seed=-49837) \n",
                    " \n",
                    "env.io.hetatm = True \n",
                    "#soft sphere potential \n",
                    "env.edat.dynamic_sphere=False \n",
                    "#lennard-jones potential (more accurate) \n",
                    "env.edat.dynamic_lennard=True \n",
                    "env.edat.contact_shell = 4.0 \n",
                    "env.edat.update_dynamic = 0.39 \n",
                    " \n",
                    "# Read customized topology file with phosphoserines (or standardd one) \n",
                    "env.libs.topology.read(file='$(LIB)/top_heav.lib') \n",
                    " \n",
                    "# Read customized CHARMM parameter library with phosphoserines (or standardd one) \n",
                    "env.libs.parameters.read(file='$(LIB)/par.lib') \n",
                    " \n",
                    " \n",
                    "# Read the original PDB file and copy its sequence to the alignment array: \n",
                    "mdl1 = model(env, file=modelname) \n",
                    "ali = alignment(env) \n",
                    "ali.append_model(mdl1, atom_files=modelname, align_codes=modelname) \n",
                    " \n",
                    "#set up the mutate residue selection segment \n",
                    "s = selection(mdl1.chains[chain].residues[respos]) \n",
                    " \n",
                    "#perform the mutate residue operation \n",
                    "s.mutate(residue_type=restyp) \n",
                    "#get two copies of the sequence.  A modeller trick to get things set up \n",
                    "ali.append_model(mdl1, align_codes=modelname) \n",
                    " \n",
                    "# Generate molecular topology for mutant \n",
                    "mdl1.clear_topology() \n",
                    "mdl1.generate_topology(ali[-1]) \n",
                    " \n",
                    " \n",
                    "# Transfer all the coordinates you can from the template native structure \n",
                    "# to the mutant (this works even if the order of atoms in the native PDB \n",
                    "# file is not standardd): \n",
                    "#here we are generating the model by reading the template coordinates \n",
                    "mdl1.transfer_xyz(ali) \n",
                    " \n",
                    "# Build the remaining unknown coordinates \n",
                    "mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES') \n",
                    " \n",
                    "#yes model2 is the same file as model1.  It's a modeller trick. \n",
                    "mdl2 = model(env, file=modelname) \n",
                    " \n",
                    "#required to do a transfer_res_numb \n",
                    "#ali.append_model(mdl2, atom_files=modelname, align_codes=modelname) \n",
                    "#transfers from 'model 2' to 'model 1' \n",
                    "mdl1.res_num_from(mdl2,ali) \n",
                    " \n",
                    "#It is usually necessary to write the mutated sequence out and read it in \n",
                    "#before proceeding, because not all sequence related information about MODEL \n",
                    "#is changed by this command (e.g., internal coordinates, charges, and atom \n",
                    "#types and radii are not updated). \n",
                    " \n",
                    "mdl1.write(file=modelname+restyp+respos+'.tmp') \n",
                    "mdl1.read(file=modelname+restyp+respos+'.tmp') \n",
                    " \n",
                    "#set up restraints before computing energy \n",
                    "#we do this a second time because the model has been written out and read in, \n",
                    "#clearing the previously set restraints \n",
                    "make_restraints(mdl1, ali) \n",
                    " \n",
                    "#a non-bonded pair has to have at least as many selected atoms \n",
                    "mdl1.env.edat.nonbonded_sel_atoms=1 \n",
                    " \n",
                    "sched = autosched.loop.make_for_model(mdl1) \n",
                    " \n",
                    "#only optimize the selected residue (in first pass, just atoms in selected \n",
                    "#residue, in second pass, include nonbonded neighboring atoms) \n",
                    "#set up the mutate residue selection segment \n",
                    "s = selection(mdl1.chains[chain].residues[respos]) \n",
                    " \n",
                    "mdl1.restraints.unpick_all() \n",
                    "mdl1.restraints.pick(s) \n",
                    " \n",
                    "s.energy() \n",
                    " \n",
                    "s.randomize_xyz(deviation=4.0) \n",
                    " \n",
                    "mdl1.env.edat.nonbonded_sel_atoms=2 \n",
                    "optimize(s, sched) \n",
                    " \n",
                    "#feels environment (energy computed on pairs that have at least one member \n",
                    "#in the selected) \n",
                    "mdl1.env.edat.nonbonded_sel_atoms=1 \n",
                    "optimize(s, sched) \n",
                    " \n",
                    "s.energy() \n",
                    " \n",
                    "#give a proper name \n",
                    "mdl1.write(file=modelname+'_'+restyp+respos+'.pdb') \n",
                    " \n",
                    "#delete the temporary file \n",
                    "os.remove(modelname+restyp+respos+'.tmp')"]
    
    with open(mutate_modelScript, "w") as mutate_modelFile:
        mutate_modelFile.writelines(mutate_model)

#############################
#This step ask the user for the list of mutations. Similar to the seqmut file
#############################
def Insert_mutations(yes_no, warning, workingFolder, pdbFile, copyFile, mutation_seqmut):

    for i in range(0, len(mutation_seqmut)):
        mutation_seqmut[i] = mutation_seqmut[i].upper()
        
    globals().update({"mutation_seqmut" : mutation_seqmut})
    mutationsFormat(yes_no, warning)

#Create a new working folder and a new pdb file for the mutation
    mutFile=''
    for i in range(0,len(mutation_seqmut)):
        mutFile = mutFile+mutation_seqmut[i]+"-"
    globals().update({ "mutFile" : mutFile})

    mut_Folder = mutFile+workingFolder 
    os.system("mkdir "+mut_Folder)
    os.system("cp "+pdbFile+" "+mut_Folder)
    os.chdir(mut_Folder)

    global mut_pdbFile, mut_pdbFileTemp, mutation_output
    mut_pdbFile = mutFile+pdbFile
    mut_pdbFileTemp = mut_pdbFile[:-3]+"temp"
    mutation_output = mut_pdbFile[:-3]+"output"
    copyFile(pdbFile, mut_pdbFile) #working mutation File

    NumIDmutations = []
    with open("seqmut"+mut_pdbFile[:-8], "w") as seqmutFile:
        for i in range(0,len(mutation_seqmut)):
            seqmutFile.writelines(mutation_seqmut[i]+"\n")
            NumID = re.findall('\d+', str(mutation_seqmut[i]))[0]
#List with the ResID numbers of the residues to be mutated is stored as NumIDmutationsList
            NumIDmutations.append(NumID)
    globals().update({"NumIDmutationsList" : NumIDmutations})

    print "\n ---> The following mutation(s) will be performed: "
    for i in range(0,len(mutation_seqmut)):
        print str(i+1)+") ",  mutation_seqmut[i]

#############################
#This step recognices non-standard residues in the mutations and unifies the format to 1 letter amino acid 
#############################        
def mutationsFormat(yes_no, warning):
    #Recognizes non-standard residues and ask the user to insert the mutation again
    NonStandarddResID(warning, yes_no)

    #Unifies the format to 1 letter amino acid 
    for i in range(0,len(mutation_seqmut)):
        if (mutation_seqmut[i][0:3]) in threeToOneDic:
            mutation_seqmut[i] = mutation_seqmut[i].replace(mutation_seqmut[i][0:3],threeToOneDic[mutation_seqmut[i][0:3]] )
        if (mutation_seqmut[i][-3:]) in threeToOneDic:
            mutation_seqmut[i] = mutation_seqmut[i].replace(mutation_seqmut[i][-3:],threeToOneDic[mutation_seqmut[i][-3:]] )
        
            globals().update({"mutation_seqmut" : mutation_seqmut})

    #Unifies the format to 3 letter amino acid 
    for i in range(0,len(mutation_seqmut)):
        if (mutation_seqmut[i][0]) in oneToTreeDic:
            mutation_seqmut[i] = mutation_seqmut[i].replace(mutation_seqmut[i][0],oneToTreeDic[mutation_seqmut[i][0]] )
        if (mutation_seqmut[i][-1:]) in oneToTreeDic:
            mutation_seqmut[i] = mutation_seqmut[i].replace(mutation_seqmut[i][-1:],oneToTreeDic[mutation_seqmut[i][-1:]] )
        
            globals().update({"mutation_seqmut" : mutation_seqmut})

#SIMPLIFY!!
def NonStandarddResID(warning, yes_no):
    #Recognizes non-standard residues 3 letter format
    for i in range(0,len(mutation_seqmut)):
        if (mutation_seqmut[i][0:3]).isalpha() == True:
            if (mutation_seqmut[i][0:3]) not in threeToOneDic:
                print "\n", warning, "The following residue is not recognized:", '\x1b[0;33;49m'+(mutation_seqmut[i][0:3])+'\x1b[0m', "\n Try again!"
                Insert_mutations(yes_no, warning, workingFolder, pdbFile, copyFile)
        else:
            if (mutation_seqmut[i][0]).isalpha() == True:
                if (mutation_seqmut[i][0]) not in oneToTreeDic:
                    print "\n", warning, "The following residue is not recognized:", '\x1b[0;33;49m'+(mutation_seqmut[i][0])+'\x1b[0m', "\n Try again!"
                    Insert_mutations(yes_no, warning, workingFolder, pdbFile, copyFile)
            else:
                print "\n", warning, "The following residue is not recognized:", '\x1b[0;33;49m'+(mutation_seqmut[i][0])+'\x1b[0m', "\n Try again!"
                Insert_mutations(yes_no, warning, workingFolder, pdbFile, copyFile)

        if (mutation_seqmut[i][-3:]).isalpha() == True:
            if (mutation_seqmut[i][-3:]) not in threeToOneDic:
                print "\n", warning, "The following residue is not recognized:", '\x1b[0;33;49m'+mutation_seqmut[i][-3:]+'\x1b[0m', "\n Try again!"
                Insert_mutations(yes_no, warning, workingFolder, pdbFile, copyFile)
        else:
            if (mutation_seqmut[i][-1:]).isalpha() == True:
                if (mutation_seqmut[i][-1:]) not in oneToTreeDic:
                    print "\n", warning, "The following residue is not recognized:", '\x1b[0;33;49m'+(mutation_seqmut[i][-1:])+'\x1b[0m', "\n Try again!"
                    Insert_mutations(yes_no, warning, workingFolder, pdbFile, copyFile)
            else:
                print "\n", warning, "The following residue is not recognized:", '\x1b[0;33;49m'+(mutation_seqmut[i][-1:])+'\x1b[0m', "\n Try again!"
                Insert_mutations(yes_no, warning, workingFolder, pdbFile, copyFile)

                
#############################
#This step ask the user for select the software for the mutations
#############################
def Select_mutation_Software(ChooseNumOption, pdbFile, copyFile):
    
    mutationSoftwareList = ["Modeller", "Scwrl4"]
    ChooseNumOption(mutationSoftwareList, "mutation_Software", "mutation_Software", '\n Choose the ', 'to perform the mutations:', 'will be used to perform the mutations.', True)

#############################
#Modeller and SCWR4 complete the missing atoms of the amino acids. 
#To preserve the structure of the wild type is necessary to include the geometry of the substitution in the wild type file. 
#############################
def MutatedToARMFormat(pdbFile, moveFile, FormatPDB1, FormatPDB, NumIDmut, chainName):

#Obtain the geometry of the new residue mutation
    geometry_Mutation=[]
    with open(mutation_output, "r") as oldfile, open(mut_pdbFileTemp, "w") as newfile:
        for line in oldfile:
            if 'ATOM' in line and chainName+'{:>4}'.format(str(NumIDmut)) in line:
                newfile.writelines(line)
                geometry_Mutation.append(line)

    moveFile(mut_pdbFileTemp, mutation_output)

#Calculate the number of atoms of old residue
    i = 0
    numAtomOldRes = ''
    with open(mut_pdbFile, "r") as oldfile:
        for line in oldfile:
            if 'ATOM' in line and chainName+'{:>4}'.format(str(NumIDmut)) in line:
                i = i+1
        numAtomOldRes = i

#Insert the geometry of the new residue mutation in the wild type geometry
    i = 0
    with open(mut_pdbFile, "r") as oldfile, open(mut_pdbFileTemp, "w") as newfile:
        for line in oldfile:
            if 'ATOM' in line and chainName+'{:>4}'.format(str(NumIDmut)) in line:
                i = i+1
                if i == numAtomOldRes:
                    newfile.writelines(geometry_Mutation)
            else:
                newfile.writelines(line)

#Write the new pdb using the correct format
    moveFile(mut_pdbFileTemp, mut_pdbFile)
    FormatPDB1(mut_pdbFile, mut_pdbFileTemp)
    FormatPDB(mut_pdbFileTemp, mut_pdbFile, mut_pdbFile)    

#############################
#Modeller routine                                                                                                               
#############################
def Modeller_mutations(chainName, pdbFile, copyFile, moveFile, FormatPDB1, FormatPDB):

    from mutate import mutate_model_Script
    mutate_model_Script()

    for i in range(0,len(mutation_seqmut)):
        print str("Running MODELLER9.19 for the mutation number " + str(i+1)).rjust(100, '.')
        os.system(modellerScript+" "+mutate_modelScript+" "+mut_pdbFile[:-4]+" "+str(NumIDmutationsList[i])+" "+str(mutation_seqmut[i][-3:])+" "+chainName)  
        moveFile(mut_pdbFile[:-4]+"_"+str(mutation_seqmut[i][-3:])+str(NumIDmutationsList[i])+".pdb", mutation_output)
        NumIDmut = NumIDmutationsList[i]
        print "\n The mutation ", mutation_seqmut[i], "has been succesfully generated!"

        MutatedToARMFormat(pdbFile, moveFile, FormatPDB1, FormatPDB, NumIDmut, chainName)
        os.chdir("../")

#############################
#SCWRL4 routine                                                                                                                 
#############################
def Scwrl4_mutations(resNumID, chainName, FormatPDB1, moveFile, pdbFile, FormatPDB):

    pdbHETATM = mut_pdbFile[:-3]+"HETATM.pdb"
    with open(mut_pdbFile, "r") as pdb, open(pdbHETATM, "w") as hetatm:
        for line in pdb:
            if "HETATM" in line:
                hetatm.writelines(line)

    global sequenceWTList
    seqFileName = mut_pdbFile[:-7]+"seqFileName"

    get_pir_Script(mut_pdbFile, resNumID, chainName, FormatPDB1, moveFile)

    for i in range(0,len(mutation_seqmut)):
        NumIDmut = NumIDmutationsList[i]
        sequenceWTList[int(NumIDmut)-1] = threeToOneDic[mutation_seqmut[i][-3:]]
        
        with open(seqFileName, "w") as seqFile:
            for j in range(0, len(sequenceWTList)):
                if sequenceWTList[j] != "-":
                    seqFile.writelines(sequenceWTList[j])

        print str("Running SCWRL4 for the mutation number " + str(i+1)).rjust(100, '.')
        os.system(scrl4Script+" -i "+mut_pdbFile+" -o "+mutation_output+" -h -f "+pdbHETATM+" -s "+seqFileName+" > scwrl4_mut.log" )
        print (scrl4Script+" -i "+mut_pdbFile+" -o "+mutation_output+" -h -f "+pdbHETATM+" -s "+seqFileName+" > scwrl4_mut.log" )

        print "\n The mutation ", mutation_seqmut[i], "has been succesfully generated!"

        sequenceWTList[int(NumIDmut)-1] = sequenceWTList[int(NumIDmut)-1].lower()

        #Fix the format of the mutation_ouput file
        FormatPDB1(mutation_output, "mut_temp")

        with open("mut_temp", "r") as out, open("mutation_output", "w") as temp:
            for line in out:
                if "ATOM" in line:
                    temp.writelines(line.split()[0]+"\t"+line.split()[1]+"\t"+line.split()[2]+"\t"+line.split()[3]+"\t"+line.split()[4]+"\t"+line.split()[5]+"\t"+line.split()[6]+"\t"+line.split()[7]+"\t"+line.split()[8]+"\t"+line.split()[9]+"\t"+str("0.0")+"\t"+line.split()[10]+"\n")
        FormatPDB("mutation_output", mutation_output, mutation_output)
        os.remove("mut_temp")

        MutatedToARMFormat(pdbFile, moveFile, FormatPDB1, FormatPDB, NumIDmut, chainName)
        os.system("pwd")

#############################
# Propka
#############################
def protonation_mutant(protonation, numberCounterions, protAA, HISName, protResDictionary, replaceLine, chainName, Step):

    def protonation_mut(mut_pdbFile, protAA, protResDictionary,replaceLine, chainName, Step):
        for key in protResDictionary:
            if protResDictionary[key] == "HIS":
                replaceLine(mut_pdbFile, "ATOM", chainName+'{:>4}'.format(key), {"HIS" : HISName})
            else:
                replaceLine(mut_pdbFile, "ATOM", chainName+'{:>4}'.format(key), protAA)

    protonation_mut(mut_pdbFile, protAA, protResDictionary,replaceLine, chainName, Step)
    numberCounterions(mut_pdbFile)
    Step('The '+'\x1b[0;33;49m'+mut_pdbFile+' file is ready to be used as input for the ARM protocol! \x1b[0m')

    os.chdir("../")




