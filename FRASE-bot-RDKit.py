#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
import glob
import shutil
from rdkit import Chem, RDConfig
from rdkit.Chem.rdmolfiles import SDWriter
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Chem import PandasTools
from rdkit.Chem import rdqueries
from rdkit.Chem import Fragments
from rdkit.Geometry import rdGeometry
from itertools import permutations


# In[2]:


def ThreeAlphaCarbon(m):
    """find all alpha carbons"""
    if m is None:
        return
    aaSmarts = '[NX3,NX4+][CX4][CX3](=[OX1])'
    aa = Chem.MolFromSmarts(aaSmarts)
    a = m.GetSubstructMatches(aa)
    aCarbon = []
    if a is not None:
        for i in a:
            for j in i:
                atom = m.GetAtomWithIdx(j)
                ri = atom.GetPDBResidueInfo().GetResidueName()
                if atom.GetSymbol() == 'C' and atom.GetHybridization() == Chem.HybridizationType.SP3 and ri != 'LIG':
                    aCarbon.append(atom.GetIdx())             
    aCarbon.sort()
    
    """Calculate the distance of any two alpha carbons"""
    ll = []
    for i in range(len(aCarbon)-1):
        for j in range(len(aCarbon)):
            if j > i:
                ll.append(aCarbon[i:j+1])

    r1=[]
    r2=[]
    dis=[]
    for i in ll:
        at1 = m.GetConformer().GetAtomPosition(i[0])
        at2 = m.GetConformer().GetAtomPosition(i[-1])
        r12 = at1.Distance(at2)
        if r12 >= 8 and r12 <= 12:
            r1.append(i[0])
            r2.append(i[-1])
            dis.append(r12)
    combo=list(zip(r1,r2,dis))

    """remove duplicates from r1 list and get the number of the 1st alpha-carbon"""
    r1_redu = []
    [r1_redu.append(x) for x in r1 if x not in r1_redu]

    """extract the second elements from the tuple of the combo in order to calculate the distance bewteen them"""
    rlist_empty = []
    a=[]
    for x in r1_redu:
        for i in combo:
            if i[0] == x:            
                rlist_empty.append(i[1])
            if i[0] != x and rlist_empty != [] or i[0] == r1_redu[-1]:
                a.append(rlist_empty)
                rlist_empty = []

    aC2=[ ele for ele in a if ele != []]

    """Calculate the distance of any two alpha carbons on the 2nd search"""
    combo1_tmp=[]
    for x in aC2:
        if len(x) == 1:
            continue
        else:
            ll = []
            for i in range(len(x)-1):
                for j in range(len(x)):
                    if j > i:
                        ll.append(x[i:j+1])                
        r1=[]
        r2=[]
        dis=[]    
        for i in ll:
            at1 = m.GetConformer().GetAtomPosition(i[0])
            at2 = m.GetConformer().GetAtomPosition(i[-1])
            r12 = at1.Distance(at2)
            if r12 >= 8 and r12 <= 12:
                r1.append(i[0])
                r2.append(i[-1])
                dis.append(r12)
        combo1=list(zip(r1,r2,dis))
        combo1_tmp.append(combo1)
    combo2=[ ele for ele in combo1_tmp if ele != []]

    """create triplets and find the resiude through alpha-C atoms"""
    atom1 = []
    atom2 = []
    atom3 = []
    for x in r1_redu:
        for i, ele1 in enumerate(combo2):
            for j, ele2 in enumerate(ele1):
                at1 = m.GetConformer().GetAtomPosition(x)
                at2 = m.GetConformer().GetAtomPosition(combo2[i][j][0])
                at3 = m.GetConformer().GetAtomPosition(combo2[i][j][1])
                r12 = at1.Distance(at2)
                r13 = at1.Distance(at3)
                
                if r12 >=8 and r12 <= 12 and r13 >=8 and r13 <= 12 and x < combo2[i][j][0]:
                    
                    atom1.append(x)
                    atom2.append(combo2[i][j][0])
                    atom3.append(combo2[i][j][1])
    three_atoms_tmp = list(zip(atom1,atom2,atom3))   
    
    three_atoms = []
    [three_atoms.append(x) for x in three_atoms_tmp if x not in three_atoms]
    three_atoms.sort()
    
    return three_atoms


# In[3]:


def fp_database(m,three_atoms):
    lig = []
#    if m is None:
#        return
    for atom in m.GetAtoms():
        ri = atom.GetPDBResidueInfo().GetResidueName()
        num = atom.GetPDBResidueInfo().GetResidueNumber()
        if ri == 'LIG':
            lig.append(ri + '     ' + str(num))
    
    """extract residue name and create triplets"""
    triplet = []
    triplet_list = []
    for ac in three_atoms:
        res = []
        reslist = []
        for i in ac:
            atom = m.GetAtomWithIdx(i)
            ri = atom.GetPDBResidueInfo().GetResidueName()
            num = atom.GetPDBResidueInfo().GetResidueNumber()
            rc = atom.GetPDBResidueInfo().GetChainId()
            #res.append(ri)
            
            if len(str(num)) == 4:
                resId = ri + " " + rc + "" + "%4d"%(num)
                reslist.append(resId)
                reslist = sorted(reslist, key=lambda x: int(x.split()[-1]))
            else:
                resId = ri + " " + rc + " " + "%3d"%(num)
                reslist.append(resId)
                reslist = sorted(reslist, key=lambda x: int(x.split()[-1]))

        res = [i.split()[0] for i in reslist]

        triplet.append(res)
        triplet_list.append(list(set(lig))+reslist)

    """distribute 11-bit fingerprint to each residue in a triplet, 5 physical properties: postive ion; negtive ion; H-bonded donor; H-bonded accptor; aromaticity; hydrophobicity; size"""
    """initalization: #ionpos:='0,';#ionneg:='0,';#hba:='0,';#hbd:='0,';#ar:='0,';#hy:='0,0,0,';#vol:='0,0,0';"""
    """generate the fingerprint encodes for three residues in each triplet"""

    all_fp = []
    for x,ele in enumerate(triplet):
        triplet_fp_six = []
        for i in permutations(ele,3):
            triplet_fp = []
            for res in i:
                if res == 'ARG' or res == 'LYS' or res == 'HIS' or res == 'HIP' or res == 'HID' or res == 'HIE':
                    ionpos = '1,'
                else:
                    ionpos = '0,'
                
                if res == 'ASP' or res == 'GLU':
                    ionneg = '1,'
                else:
                    ionneg = '0,'
                
                if res == 'ASP' or res == 'GLU'  or res == 'SER' or res == 'THR' or res == 'ASN' or res == 'GLN' or res == 'TYR' or res == 'HIS' or res == 'HIP' or res == 'HID' or res == 'HIE':
                    hba = '1,'
                else:
                    hba = '0,'
                
                if res == 'ARG' or res == 'LYS'  or res == 'SER' or res == 'THR' or res == 'ASN' or res == 'GLN' or res == 'TYR' or res == 'HIS' or res == 'HIP' or res == 'HID' or res == 'HIE' or res == 'CYS' or res == 'TRP':
                    hbd = '1,'
                else:
                    hbd = '0,'
                
                if res == 'PHE' or res == 'TYR' or res == 'HIS' or res == 'HIP' or res == 'HID' or res == 'HIE' or res == 'TRP':
                    ar = '1,'
                else:
                    ar = '0,'
                
                if res == 'GLY':
                    hy = '0,0,0,'
                elif res == 'ALA' or res == 'PRO':
                    hy = '1,0,0,'
                elif res == 'VAL' or res == 'ILE' or res == 'LEU' or res == 'THR' or res == 'CYS':
                    hy = '1,1,0,'
                elif res == 'MET' or res == 'PHE' or res == 'TYR' or res == 'TRP':
                    hy = '1,1,1,'
                else:
                    hy = '0,0,0,'
                
                if res == 'GLY':
                    vol = '0,0,0'
                elif res == 'ALA' or res == 'PRO':
                    vol = '1,0,0'
                elif res == 'VAL' or res == 'ILE' or res == 'LEU' or res == 'THR' or res == 'CYS' or res == 'SER' or res == 'ASN' or res == 'ASP' or res == 'HIS' or res == 'HIP' or res == 'HID' or res == 'HIE':
                    vol = '1,1,0'
                elif res == 'MET' or res == 'PHE' or res == 'TYR' or res == 'TRP' or res == 'ARG' or res == 'LYS' or res == 'GLU' or res == 'GLN':
                    vol = '1,1,1'
                else:
                    vol = '0,0,0'
                fp = ionpos + ionneg + hba + hbd + ar + hy + vol
                triplet_fp.append(fp)
                
            triplet_fp_six.append(triplet_fp)
        all_fp.append(triplet_fp_six)
    fp_res = list(zip(all_fp,triplet_list,three_atoms))
    return fp_res


# In[4]:


def fp_target(m,three_atoms):
    """extract residue name and create triplets"""
    triplet = []
    triplet_list = []
#    if m is None:
#        return
    for ac in three_atoms:
        res = []
        reslist = []
        for i in ac:
            atom = m.GetAtomWithIdx(i)
            ri = atom.GetPDBResidueInfo().GetResidueName()
            num = atom.GetPDBResidueInfo().GetResidueNumber()
            rc = atom.GetPDBResidueInfo().GetChainId() 
            res.append(ri)
            
            if len(str(num)) == 4:
                resId = ri + " " + rc + "" + "%4d"%(num)
                reslist.append(resId)
            else:
                resId = ri + " " + rc + " " + "%3d"%(num)
                reslist.append(resId)

        triplet.append(res)
        triplet_list.append(reslist)

    """distribute 11-bit fingerprint to each residue in a triplet, 5 physical properties: postive ion; negtive ion; H-bonded donor; H-bonded accptor; aromaticity; hydrophobicity; size"""
    """initalization: #ionpos:='0,';#ionneg:='0,';#hba:='0,';#hbd:='0,';#ar:='0,';#hy:='0,0,0,';#vol:='0,0,0';"""
    """generate the fingerprint encodes for three residues in each triplet"""

    all_fp = []
    for x,ele in enumerate(triplet):
        triplet_fp = []
        for res in ele:
            if res == 'ARG' or res == 'LYS' or res == 'HIS' or res == 'HIP' or res == 'HID' or res == 'HIE':
                ionpos = '1,'
            else:
                ionpos = '0,'
                
            if res == 'ASP' or res == 'GLU':
                ionneg = '1,'
            else:
                ionneg = '0,'
                
            if res == 'ASP' or res == 'GLU'  or res == 'SER' or res == 'THR' or res == 'ASN' or res == 'GLN' or res == 'TYR' or res == 'HIS' or res == 'HIP' or res == 'HID' or res == 'HIE':
                hba = '1,'
            else:
                hba = '0,'
                
            if res == 'ARG' or res == 'LYS'  or res == 'SER' or res == 'THR' or res == 'ASN' or res == 'GLN' or res == 'TYR' or res == 'HIS' or res == 'HIP' or res == 'HID' or res == 'HIE' or res == 'CYS' or res == 'TRP':
                hbd = '1,'
            else:
                hbd = '0,'
                
            if res == 'PHE' or res == 'TYR' or res == 'HIS' or res == 'HIP' or res == 'HID' or res == 'HIE' or res == 'TRP':
                ar = '1,'
            else:
                ar = '0,'
                
            if res == 'GLY':
                hy = '0,0,0,'
            elif res == 'ALA' or res == 'PRO':
                hy = '1,0,0,'
            elif res == 'VAL' or res == 'ILE' or res == 'LEU' or res == 'THR' or res == 'CYS':
                hy = '1,1,0,'
            elif res == 'MET' or res == 'PHE' or res == 'TYR' or res == 'TRP':
                hy = '1,1,1,'
            else:
                hy = '0,0,0,'
                
            if res == 'GLY':
                vol = '0,0,0'
            elif res == 'ALA' or res == 'PRO':
                vol = '1,0,0'
            elif res == 'VAL' or res == 'ILE' or res == 'LEU' or res == 'THR' or res == 'CYS' or res == 'SER' or res == 'ASN' or res == 'ASP' or res == 'HIS' or res == 'HIP' or res == 'HID' or res == 'HIE':
                vol = '1,1,0'
            elif res == 'MET' or res == 'PHE' or res == 'TYR' or res == 'TRP' or res == 'ARG' or res == 'LYS' or res == 'GLU' or res == 'GLN':
                vol = '1,1,1'
            else:
                vol = '0,0,0'
            fp = ionpos + ionneg + hba + hbd + ar + hy + vol
            triplet_fp.append(fp)
            
        all_fp.append(triplet_fp)
    fp_res = list(zip(all_fp,triplet_list,three_atoms))
    return fp_res


# In[5]:


def WriteTripletTarget(triplet_target,target_file,pdbname,index):
    file_path_target=f'triplet_db/tmp_target_triplet_{pdbname}_{index}_{triplet_target}.pdb'
    for r in triplet_target:
        rfile_target = open(target_file,"r")
        wfile_target = open(file_path_target,"a")
        for line in rfile_target.readlines():
            if r in line:
                wfile_target.write(line)
        rfile_target.close()
        wfile_target.close()
    return

def WriteTripletFrasedb(triplet_db, frasedb_file, pdbname,index):
    file_path_db = f'triplet_db/tmp_FRASEdb_triplet_{pdbname}_{index}_{triplet_db}.pdb'                
    
    if os.path.exists(file_path_db):
        #print(file_path_db, 'exists.')
        return  # Skip the rest of the function
    
    for r in triplet_db:
        rfile_db = open(frasedb_file, "r")                     
        wfile_db = open(file_path_db, "a")
        for line in rfile_db.readlines():
            if r in line:
                wfile_db.write(line)
        wfile_db.close()
        rfile_db.close()
    
    # Add bond connectivity of ligand fragment
    if len(triplet_db) == 4:
        rfile_db = open(frasedb_file, "r")                     
        wfile_db = open(file_path_db, "a")
        for line in rfile_db.readlines():
            if line.split()[0] == 'CONECT':
                wfile_db.write(line)
        wfile_db.close()
        rfile_db.close()


# In[6]:


def mol_with_atom_index(mol):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol

def removecopy(name_m):
    clean = list(set(name_m))
    return clean

# Identify Aromatic Rings
def isRingAromatic(mol, bondRing):
    for id in bondRing:
        if not mol.GetBondWithIdx(id).GetIsAromatic():
            return False
    return True
def AromaticBonds(m):
    ri = m.GetRingInfo()
    #print(ri.BondRings()) 
    bonds_aromatic = []
    for i in ri.BondRings():
        if isRingAromatic(m, i):
            bonds_aromatic.append(isRingAromatic(m, i))
    return bonds_aromatic
def AromaticAtoms(m):
    aromatic_atoms = [i.GetIdx() for i in m.GetAtoms() if i.GetIsAromatic()]
    aromatic_atoms.sort()
    return aromatic_atoms
def AromaticCycles(m):
    ri = m.GetRingInfo()
    AromaticCycles = []
    for i in range(0,len(ri.BondRings())):
        log = isRingAromatic(m, ri.BondRings()[i])
        if log is True:
            AromaticCycles.append(ri.AtomRings()[i])
    return AromaticCycles

# Identify Aromatic Nitrogen Atoms
def AromaticNitrogen(m):
    q_ar = rdqueries.IsAromaticQueryAtom()
    aromatic_nitrogen = [x.GetIdx() for x in m.GetAtomsMatchingQuery(q_ar) if x.GetSymbol() == 'N']
    aromatic_nitrogen.sort()
    return aromatic_nitrogen

# Identify Aromatic Oxygen Atoms
def AromaticOxygen(m):
    q_ar = rdqueries.IsAromaticQueryAtom()
    aromatic_oxygen = [x.GetIdx() for x in m.GetAtomsMatchingQuery(q_ar) if x.GetSymbol() == 'O']
    aromatic_oxygen.sort()
    return aromatic_oxygen

# Identify Aromatic Sulfur Atoms
def AromaticSulfur(m):
    q_ar = rdqueries.IsAromaticQueryAtom()
    aromatic_sulfur = [x.GetIdx() for x in m.GetAtomsMatchingQuery(q_ar) if x.GetSymbol() == 'S']
    aromatic_sulfur.sort()
    return aromatic_sulfur

# Identify Apliphatic Carbon Atoms
"""delete carbonds from alpha-C and C=O group"""
def AliphaticCarbon(m):
    q_ali = rdqueries.IsAliphaticQueryAtom()
    aliphatic_carbon = [x.GetIdx() for x in m.GetAtomsMatchingQuery(q_ali) if not x.IsInRing() and x.GetSymbol() == 'C']
    aliphatic_carbon.sort()
    return aliphatic_carbon

def backboneCarbons(m):
    aaSmarts = '[NX3,NX4+][CX4H][CX3](=[OX1])'
    aa = Chem.MolFromSmarts(aaSmarts)
    a = m.GetSubstructMatches(aa)
    bbCarbon = []
    if a is not None:
        for i in a:
            for j in i:
                atom = m.GetAtomWithIdx(j)
                if atom.GetSymbol() == 'C':
                    bbCarbon.append(atom.GetIdx())
    bbCarbon.sort()
    return bbCarbon

def AliphaticCarbons(m):
    new_list = []
    for elem in AliphaticCarbon(m):
        if elem not in backboneCarbons(m):
            new_list.append(elem)
    return new_list

# Identify halogen
def halogen(m):
    q_atom = rdqueries.AAtomQueryAtom()
    halogen = [x.GetIdx() for x in m.GetAtomsMatchingQuery(q_atom) if x.GetSymbol() == 'Cl' or x.GetSymbol() == 'F' or x.GetSymbol() == 'Br' or x.GetSymbol() == 'I']
    halogen.sort()
    return halogen

# NegativeIonizable
def NegativeIonizable(m):
    q_negatom = rdqueries.FormalChargeLessQueryAtom(0)
    negatom = [x.GetIdx() for x in m.GetAtomsMatchingQuery(q_negatom) if x.GetFormalCharge() < 0]
    negatom.sort()
    return negatom

# PositiveIonizable
def PositiveIonizable(m):
    q_posatom = rdqueries.FormalChargeGreaterQueryAtom(0)
    posatom = [x.GetIdx() for x in m.GetAtomsMatchingQuery(q_posatom) if x.GetFormalCharge() > 0]
    posatom.sort()
    return posatom

# identify NitroOxygen
def NitroOxygen(m):
    #nitroSmarts = '[N;D3](=[O;D1])[O;D1]'
    nitroSmarts = '[N;H0;$(N-[#6]);D3](=[O;D1])~[O;D1]'
    nitro = Chem.MolFromSmarts(nitroSmarts)
    a = m.GetSubstructMatches(nitro)
    nitroOxygen = []
    if a is not None:
        for i in a:
            for j in i:
                atom = m.GetAtomWithIdx(j)
                if atom.GetSymbol() == 'O':
                    nitroOxygen.append(atom.GetIdx())
    nitroOxygen.sort()
    return nitroOxygen

# identify NitroNitrogen
def NitroNitrogen(m):
    #nitroSmarts = '[N;D3](=[O;D1])[O;D1]'
    nitroSmarts = '[N;H0;$(N-[#6]);D3](=[O;D1])~[O;D1]'
    nitro = Chem.MolFromSmarts(nitroSmarts)
    a = m.GetSubstructMatches(nitro)
    nitroNitrogen = []
    if a is not None:
        for i in a:
            for j in i:
                atom = m.GetAtomWithIdx(j)
                if atom.GetSymbol() == 'N':
                    nitroNitrogen.append(atom.GetIdx())
    nitroNitrogen.sort()
    return nitroNitrogen

# identify AmideNitrogen
def AmideNitrogen(m):
    amideSmarts = '[#7]-[#6]=[#8]'
    amide = Chem.MolFromSmarts(amideSmarts)
    a = m.GetSubstructMatches(amide)
    amideNitrogen = []
    if a is not None:
        for i in a:
            for j in i:
                atom = m.GetAtomWithIdx(j)
                if atom.GetSymbol() == 'N' and atom.GetTotalDegree() == 3:
                    amideNitrogen.append(atom.GetIdx())
    amideNitrogen.sort()
    return amideNitrogen

# identify AmineNitrogen
def AmineNitrogen(m):
    #amineNitrogenSmarts = '[N;D1]'
    amineNitrogenSmarts = '[N;$(N-[#6]);!$(N-[!#6;!#1]);!$(N-C=[O,N,S])]'
    amine = Chem.MolFromSmarts(amineNitrogenSmarts)
    a = m.GetSubstructMatches(amine)
    amineNitrogen = []
    if a is not None:
        for i in a:
            for j in i:
                atom = m.GetAtomWithIdx(j)
                if atom.GetTotalDegree() == 3:
                    amineNitrogen.append(atom.GetIdx())
    amineNitrogen.sort()
    return amineNitrogen

# identify CarboxylateOxygen
def CarboxylateOxygen(m):
    #carboxylicSmarts = 'C(=O)[O;D1]'
    carboxylicSmarts = 'C(=O)[O;H,-]'
    carboxylic = Chem.MolFromSmarts(carboxylicSmarts)
    a = m.GetSubstructMatches(carboxylic)
    carboxylateOxygen = []
    if a is not None:
        for i in a:
            for j in i:
                atom = m.GetAtomWithIdx(j)
                if atom.GetSymbol() == 'O':
                    carboxylateOxygen.append(atom.GetIdx())
    carboxylateOxygen.sort()
    return carboxylateOxygen

# identify AlcoholOxygen
def AlcoholOxygen(m):
#    alcoholSmarts = '[O;$(OCC)]'
    alcoholSmarts = '[O;H1;$(O-!@[#6;!$(C=!@[O,N,S])])]'
    alcohol = Chem.MolFromSmarts(alcoholSmarts)
    a = m.GetSubstructMatches(alcohol)
    alcoholOxygen = []
    if a is not None:
        for i in a:
            for j in i:
                atom = m.GetAtomWithIdx(j)
                if atom.GetSymbol() == 'O' and not atom.IsInRing():
                    nbrs = [x for x in atom.GetNeighbors()]
                    for nbr in nbrs:
                        if not nbr.IsInRing() and len(nbr.GetNeighbors()) == 2:
                            alcoholOxygen.append(atom.GetIdx())
    alcoholOxygen.sort()
    return alcoholOxygen

# identify EnolOxygen
def EnolOxygen(m):
    #enolSmarts = '[#6]=[#6]-[#8;H]'
    enolSmarts = '[C;D2]=[C;D2]-[O;D1;H]'
    enol = Chem.MolFromSmarts(enolSmarts)
    a = m.GetSubstructMatches(enol)
    enolOxygen = []
    if a is not None:
        for i in a:
            for j in i:
                atom = m.GetAtomWithIdx(j)
                if atom.GetSymbol() == 'O' and not atom.IsInRing():
                    nbrs = [x for x in atom.GetNeighbors()]
                    for nbr in nbrs:
                        if not nbr.IsInRing() and len(nbr.GetNeighbors()) == 2:
                            enolOxygen.append(atom.GetIdx())
    enolOxygen.sort()
    return enolOxygen

# identify EnamineNitrogen
def EnamineNitrogen(m):
    #enamineSmarts = '[#6]=[#6]-[#7]'
    #enamineSmarts = '[N;D1]-[C;D2]=[C;D2]'
    enamineSmarts = '[N;$(NC=[C])]'
    enamine = Chem.MolFromSmarts(enamineSmarts)
    a = m.GetSubstructMatches(enamine)
    enamineNitrogen = []
    if a is not None:
        for i in a:
            for j in i:
                atom = m.GetAtomWithIdx(j)
                if atom.GetSymbol() == 'N' and not atom.IsInRing() and len(atom.GetNeighbors()) == 3:
                    nbrs = [x for x in atom.GetNeighbors()]
                    for nbr in nbrs:
                        if not nbr.IsInRing() and len(nbr.GetNeighbors()) == 2:
                            enamineNitrogen.append(atom.GetIdx())
    enamineNitrogen.sort()
    return enamineNitrogen

# identify ImineNitrogen
def ImineNitrogen(m):
    imineSmarts = '[N;$(N=C)]'
    imine = Chem.MolFromSmarts(imineSmarts)
    a = m.GetSubstructMatches(imine)
    imineNitrogen = []
    if a is not None:
        for i in a:
            for j in i:
                atom = m.GetAtomWithIdx(j)
                if atom.GetSymbol() == 'N' and not atom.IsInRing() and len(atom.GetNeighbors()) == 2:
                    nbrs = [x for x in atom.GetNeighbors()]
                    for nbr in nbrs:
                        if not nbr.IsInRing() and len(nbr.GetNeighbors()) == 2:
                            imineNitrogen.append(atom.GetIdx())
    imineNitrogen.sort()
    return imineNitrogen

# identify SulfoneSulfur
def SulfoneSulfur(m):
    #sulfoneSmarts = '[S;D4](=O)(=O)'
    sulfoneSmarts = '[$(S-!@[#6])](=O)(=O)'
    sulfone = Chem.MolFromSmarts(sulfoneSmarts)
    a = m.GetSubstructMatches(sulfone)
    sulfoneSulfur = []
    if a is not None:
        for i in a:
            for j in i:
                atom = m.GetAtomWithIdx(j)
                if atom.GetSymbol() == 'S' and not atom.IsInRing() and len(atom.GetNeighbors()) == 4:
                    sulfoneSulfur.append(atom.GetIdx())
    sulfoneSulfur.sort()
    return removecopy(sulfoneSulfur)

# identify SulfoxideSulfur
def SulfoxideSulfur(m):
    sulfoxideSmarts = '[S;D3](=O)'
    sulfoxide = Chem.MolFromSmarts(sulfoxideSmarts)
    a = m.GetSubstructMatches(sulfoxide)
    sulfoxideSulfur = []
    if a is not None:
        for i in a:
            for j in i:
                atom = m.GetAtomWithIdx(j)
                if atom.GetSymbol() == 'S' and not atom.IsInRing() and len(atom.GetNeighbors()) == 3 and atom.GetTotalValence() == 4:
                    nbrs = atom.GetNeighbors()
                    for nbr in nbrs:
                        i = 0
                        if nbr.GetSymbol()== 'O':
                            i = i + 1
                            if i == 2:
                                break
                            sulfoxideSulfur.append(atom.GetIdx())
    sulfoxideSulfur.sort()
    return removecopy(sulfoxideSulfur)

#  Identify Vinyl C=C
def Vinyl(m):
    vinylSmarts = '[#6]=[#6]'
    vinyl = Chem.MolFromSmarts(vinylSmarts)
    a = m.GetSubstructMatches(vinyl)
    vinylgroup = []
    if a is not None:
        for i in a:
            for j in i:
                atom = m.GetAtomWithIdx(j)
                vinylgroup.append(atom.GetIdx())
    return vinylgroup

# Identify HBondDonor and HBondAcceptor
def HBondDA(m):
    atomPharma = {}
    # Define N, O, S
    chalcogen = [7, 8, 16]
    mol = Chem.AddHs(m)
    # Start to search
    for atom in mol.GetAtoms():
        p = []
        if atom.GetAtomicNum() == 1 or atom.GetAtomicNum() not in chalcogen:
            continue
        else:
            # HBD
            if atom.GetFormalCharge() == 0:
                nbrs = [x for x in atom.GetNeighbors()]
                HDflag = False
                for nbr in nbrs:
                    if nbr.GetAtomicNum() == 1:
                        HDflag = True
                if HDflag:
                    p.append('HDonor')
                    
            # HBA
            if atom.GetTotalValence() == 2:
                nbrs = [x for x in atom.GetNeighbors()]
                HAflag_1 = True
                HAflag_2 = True
                if len(nbrs) == 1:
                    nbr = nbrs[0]
                    if nbr.GetAtomicNum() == 7:
                        HAflag_1 = False
                else:
                    for nbr in nbrs:
                        if nbr.GetAtomicNum() == 1:
                            HAflag_2 = False
                if HAflag_1 and HAflag_2:
                    p.append('HAcceptor')
                    
            if atom.GetTotalValence() == 3:
                nbrs = [x for x in atom.GetNeighbors()]
                HAflag_1 = True
                HAflag_2 = True
                if len(nbrs) == 1:
                    nbr = nbrs[0]
                    if nbr.GetAtomicNum() == 8:
                        HAflag_1 = False
                else:
                    for nbr in nbrs:
                        if nbr.GetAtomicNum() == 1:
                            HAflag_2 = False
                if HAflag_1 and HAflag_2:
                    p.append('HAcceptor')
                    
        atomPharma[atom.GetIdx()] = [atom.GetAtomicNum(), ' '.join(p)]
    return atomPharma

def HBAList(m):
    res = HBondDA(m)
    HBA = []
    for k, v in res.items():
        for x in v:        
            if x == 'HAcceptor':
                HBA.append(k)
    HBA.sort()
    return HBA

def HBDList(m):
    res = HBondDA(m)
    HBD = []
    for k, v in res.items():
        for x in v:        
            if x == 'HDonor':
                HBD.append(k)
    HBD.sort()
    return HBD

def GetCoorAromaticCenter(m,aromatic_index):
    """Get geometric center coordinates of aromatic cycles"""
    tot = rdGeometry.Point3D(0.0,0.0,0.0)
    pos = [m.GetConformer().GetAtomPosition(i) for i in aromatic_index]
    for ele in range(0, len(pos)):
        tot = tot + pos[ele]
    aver = tot/len(pos)
    return aver


# In[0]:


# work on a set of molecules
def InteractionFingerprint(pdbname,triplet_target,index):
    
    # print title line
    #sys.stdout = open('FRASE-atomic-convolutions.txt',mode = 'w', encoding = 'utf-8')
    #print('FRASE_ID',end='\t')
    #for i in range(264):
    #    print("i%03d"%(i+1),end='\t')
    #print(end='\n')
    
    SDFFile = str(f'newFragments/frase_in_target_{pdbname}_{index}_{triplet_target}.sdf')
    mol = Chem.SDMolSupplier(SDFFile)
    mols = [x for x in mol if x is not None]
    m = mols[0]
    if m is None:
        return

    #print(m.GetProp('PDB_ID'))    # what we need!
    fraseID = m.GetProp('FRASE_ID')  # what we need!
    output1 = f'tmp_{pdbname}_{index}_{triplet_target}.sdf'
    writer = Chem.SDWriter(output1)
    writer.write(m)
    writer.close()
    
    AtNum = []
    AtName = []
    fp1 = open(output1,"r")
    for line in fp1.readlines():
        if line.startswith('V  '):
            AtNum.append(line.split()[1])
            AtName.append(line.split()[2])
    fp1.close()
    os.remove(output1)
    
    com = list(zip(AtNum,AtName))
    lig_index = []
    res_index = []
    for index,label in com[:]:
        if label == 'LIG':
            lig_index.append(int(index)-1)
        else:
            res_index.append(int(index)-1)

    if AromaticCycles(m) is not None:
        lig_arocycles = [i for i in AromaticCycles(m) if sorted(i)[0] in lig_index]
        res_arocycles = [i for i in AromaticCycles(m) if sorted(i)[0] in res_index]
    if NegativeIonizable(m) is not None:
        lig_neg = [i for i in NegativeIonizable(m) if i in lig_index]
        res_neg = [i for i in NegativeIonizable(m) if i in res_index]
    if PositiveIonizable(m) is not None:
        lig_pos = [i for i in PositiveIonizable(m) if i in lig_index]
        res_pos = [i for i in PositiveIonizable(m) if i in res_index]
    if HBDList(m) is not None:
        lig_hbd = [i for i in HBDList(m) if i in lig_index]
        res_hbd = [i for i in HBDList(m) if i in res_index]
    if HBAList(m) is not None:
        lig_hba = [i for i in HBAList(m) if i in lig_index]
        res_hba = [i for i in HBAList(m) if i in res_index]
    if AmideNitrogen(m) is not None:
        lig_amideN = [i for i in AmideNitrogen(m) if i in lig_index]
        res_amideN = [i for i in AmideNitrogen(m) if i in res_index]
    if AmineNitrogen(m) is not None:
        lig_amineN = [i for i in AmineNitrogen(m) if i in lig_index]
        res_amineN = [i for i in AmineNitrogen(m) if i in res_index]
    if CarboxylateOxygen(m) is not None:
        lig_carO = [i for i in CarboxylateOxygen(m) if i in lig_index]
        res_carO = [i for i in CarboxylateOxygen(m) if i in res_index]
    if AlcoholOxygen(m) is not None:
        lig_alcO = [i for i in AlcoholOxygen(m) if i in lig_index]
        res_alcO = [i for i in AlcoholOxygen(m) if i in res_index]
    if ImineNitrogen(m) is not None:
        lig_imineN = [i for i in ImineNitrogen(m) if i in lig_index]
        res_imineN = [i for i in ImineNitrogen(m) if i in res_index]
    if AromaticNitrogen(m) is not None:
        lig_aroN = [i for i in AromaticNitrogen(m) if i in lig_index]
        res_aroN = [i for i in AromaticNitrogen(m) if i in res_index]
    if AliphaticCarbon(m) is not None:
        lig_aliC = [i for i in AliphaticCarbon(m) if i in lig_index]
        res_aliC = [i for i in AliphaticCarbons(m) if i in res_index]
    if halogen(m) is not None:
        lig_halo = [i for i in halogen(m) if i in lig_index]
        res_halo = [i for i in halogen(m) if i in res_index]
    if Vinyl(m) is not None:
        lig_vinyl = [i for i in Vinyl(m) if i in lig_index]
        res_vinyl = [i for i in Vinyl(m) if i in res_index]
    if NitroOxygen(m) is not None:
        lig_nitroO = [i for i in NitroOxygen(m) if i in lig_index]
        res_nitroO = [i for i in NitroOxygen(m) if i in res_index]
    if NitroNitrogen(m) is not None:
        lig_nitroN = [i for i in NitroNitrogen(m) if i in lig_index]
        res_nitroN = [i for i in NitroNitrogen(m) if i in res_index]
    if SulfoneSulfur(m) is not None:
        lig_sulfoneS = [i for i in SulfoneSulfur(m) if i in lig_index]
        res_sulfoneS = [i for i in SulfoneSulfur(m) if i in res_index]
    if SulfoxideSulfur(m) is not None:
        lig_sulfoxideS = [i for i in SulfoxideSulfur(m) if i in lig_index]
        res_sulfoxideS = [i for i in SulfoxideSulfur(m) if i in res_index]
    if EnolOxygen(m) is not None:
        lig_enolO = [i for i in EnolOxygen(m) if i in lig_index]
        res_enolO = [i for i in EnolOxygen(m) if i in res_index]
    if EnamineNitrogen(m) is not None:
        lig_enamN = [i for i in EnamineNitrogen(m) if i in lig_index]
        res_enamN = [i for i in EnamineNitrogen(m) if i in res_index]
    if AromaticOxygen(m) is not None:
        lig_aroO = [i for i in AromaticOxygen(m) if i in lig_index]
        res_aroO = [i for i in AromaticOxygen(m) if i in res_index]
    if AromaticSulfur(m) is not None:
        lig_aroS = [i for i in AromaticSulfur(m) if i in lig_index]
        res_aroS = [i for i in AromaticSulfur(m) if i in res_index]

    lig_listall = [lig_arocycles] + [lig_neg] + [lig_pos] + [lig_hbd] + [lig_hba] + [lig_halo] + [lig_amideN] + [lig_amineN] + [lig_vinyl] + [lig_carO] + [lig_alcO] + [lig_nitroO] + [lig_nitroN] + [lig_sulfoneS] + [lig_sulfoxideS] + [lig_enolO] + [lig_imineN] + [lig_enamN] + [lig_aroN] + [lig_aroO] + [lig_aroS] + [lig_aliC]
    prot_listall = [res_arocycles] + [res_neg] + [res_pos] + [res_hbd] + [res_hba] + [res_amideN] + [res_amineN] + [res_carO] + [res_alcO] + [res_imineN] + [res_aroN] + [res_aliC]
    
    lig_center_list = lig_arocycles
    prot_center_list = res_arocycles
    
    lig_polar_list = [lig_neg] + [lig_pos] + [lig_hbd] + [lig_hba]
    prot_polar_list = [res_neg] + [res_pos] + [res_hbd] + [res_hba] 

    lig_apolar_list = [lig_halo] + [lig_amideN] + [lig_amineN] + [lig_vinyl] + [lig_carO] + [lig_alcO] + [lig_nitroO] + [lig_nitroN] + [lig_sulfoneS] + [lig_sulfoxideS] + [lig_enolO] + [lig_imineN] + [lig_enamN] + [lig_aroN] + [lig_aroO] + [lig_aroS] + [lig_aliC]
    prot_apolar_list = [res_amideN] + [res_amineN] + [res_carO] + [res_alcO] + [res_imineN] + [res_aroN] + [res_aliC]

    #print(m.GetProp('PDB_ID'),end='\t')
    #print('LIG_polar_length:',[len(i) for i in lig_polar_list], end='\t')
    #print('Prot_polar_length:',[len(i) for i in prot_polar_list])
    
    """Contribution from a single atom pair lig_i/prot_j is equal to a/rij^2, 
    where a=9 for polar interactions (HBD, HBA, NegetiveIons, PositiveIons);
    if nonpolar interaction, then e = a/rij^4""" 
    
    inter_finger = []
    
    # ligand_center to protein_center
    for index_lig, lig in enumerate([lig_center_list]):
        for index_res, res in enumerate([prot_center_list]):
            e = 0
            if len(lig) == 0 or len(res) == 0:
                inter_finger.append(e)
            else:
                for i in lig:
                    for j in res:
                        position_centeri = GetCoorAromaticCenter(m,i)
                        position_centerj = GetCoorAromaticCenter(m,j)
                        rij = position_centeri.Distance(position_centerj)
                        ie = 9 / (rij * rij * rij * rij)
                        e = e + ie
                        #DEBUG3:
                        #print('atom',lig,'atom',res,'rij=',rij,'E(',lig,'-',res,')=',ie,'Etot=',e)
                inter_finger.append(e)
    
    # ligand_center to protein_polar
    for index_lig, lig in enumerate([lig_center_list]):
        for index_res, res in enumerate(prot_polar_list):  
            e = 0
            if len(lig) == 0 or len(res) == 0:
                inter_finger.append(e)
            else:
                for i in lig:
                    for j in res:
                        position_centeri = GetCoorAromaticCenter(m,i)
                        position_j = m.GetConformer().GetAtomPosition(j) 
                        rij = position_centeri.Distance(position_j)
                        ie = 9 / (rij * rij)
                        e = e + ie
                        #print('rij=',rij,'E(',i,'-',j,')=',ie,'Etot=',e)
                inter_finger.append(e)
    
    # ligand_center to protein_nonpolar
    for index_lig, lig in enumerate([lig_center_list]):
        for index_res, res in enumerate(prot_apolar_list):  
            e = 0
            if len(lig) == 0 or len(res) == 0:
                inter_finger.append(e)
            else:
                for i in lig:
                    for j in res:
                        position_centeri = GetCoorAromaticCenter(m,i)
                        position_j = m.GetConformer().GetAtomPosition(j) 
                        rij = position_centeri.Distance(position_j)
                        ie = 9 / (rij * rij * rij * rij)
                        e = e + ie
                        #print('rij=',rij,'E(',i,'-',j,')=',ie,'Etot=',e)
                inter_finger.append(e)
                        
    # ligand_polar to 13 protein types
    for index_lig, lig in enumerate(lig_polar_list):
        for index_res, res in enumerate(prot_listall):  
            e = 0
            if len(lig) == 0 or len(res) == 0:
                inter_finger.append(e)
            else:
                for i in lig:
                    for j in res:
                        # ligand_polar to protein_center
                        if type(j) is tuple:
                            position_i = m.GetConformer().GetAtomPosition(i)
                            position_centerj = GetCoorAromaticCenter(m,j)
                            rij = position_i.Distance(position_centerj)
                            ie = 9 / (rij * rij)
                            e = e + ie
                            #print('rij=',rij,'E(',i,'-',j,')=',ie,'Etot=',e)
                            
                        # ligand_polar to protein_polar and nonpolar
                        if type(j) is not tuple:
                            position_i = m.GetConformer().GetAtomPosition(i)
                            position_j = m.GetConformer().GetAtomPosition(j) 
                            rij = position_i.Distance(position_j)
                            ie = 9 / (rij * rij)
                            e = e + ie
                            #print('rij=',rij,'E(',i,'-',j,')=',ie,'Etot=',e)
                inter_finger.append(e)   
                
    # ligand_nonpolar to 13 protein types
    for index_lig, lig in enumerate(lig_apolar_list):
        for index_res, res in enumerate(prot_listall):  
            e = 0
            if len(lig) == 0 or len(res) == 0:
                inter_finger.append(e)
            else:
                for i in lig:
                    for j in res:
                        # ligand_nonpolar to protein_center
                        if type(j) is tuple:
                            position_i = m.GetConformer().GetAtomPosition(i)
                            position_centerj = GetCoorAromaticCenter(m,j)
                            rij = position_i.Distance(position_centerj)
                            ie = 9 / (rij * rij * rij * rij)
                            e = e + ie
                             
                        # ligand_nonpolar to protein_polar
                        if type(j) is not tuple and j in [x for y in prot_polar_list for x in y]:
                            position_i = m.GetConformer().GetAtomPosition(i)
                            position_j = m.GetConformer().GetAtomPosition(j) 
                            rij = position_i.Distance(position_j)
                            ie = 9 / (rij * rij)
                            e = e + ie
                            
                        # ligand_nonpolar to protein_nonpolar
                        if type(j) is not tuple and j in [x for y in prot_apolar_list for x in y]:
                            position_i = m.GetConformer().GetAtomPosition(i)
                            position_j = m.GetConformer().GetAtomPosition(j) 
                            rij = position_i.Distance(position_j)
                            ie = 9 / (rij * rij *rij *rij)
                            e = e + ie
                inter_finger.append(e)             

    return inter_finger


# In[7]:


def Alignment(triplet_db,triplet_target,pdbname,index):
    file_path_db = f'triplet_db/tmp_FRASEdb_triplet_{pdbname}_{index}_{triplet_db}.pdb'                
    file_path_target=f'triplet_db/tmp_target_triplet_{pdbname}_{index}_{triplet_target}.pdb'
    if not os.path.exists(file_path_db) or not os.path.exists(file_path_target):
        sys.exit(1)  # Exit with an error code
    #try:
    #    mol1 = Chem.MolFromPDBFile(file_path_db)
    #    mol2 = Chem.MolFromPDBFile(file_path_target)
    #    if mol1 is None or mol2 is None:
    #        sys.exit(1)  
   # except Exception as e:
   #     sys.exit(1)  
    mol1 = Chem.MolFromPDBFile(file_path_db)
    mol2 = Chem.MolFromPDBFile(file_path_target)
    if mol1 is None or mol2 is None:
        return  
    """Generate the atomMap for probe (mol1,FRASE_db) and reference (mol2,target_prot) molecules"""
    probe_atomMap = []
    refer_atomMap = []
    if mol1 is not None and mol2 is not None:
        for atom in mol1.GetAtoms():
            ri = atom.GetPDBResidueInfo()
            atname = ri.GetName()
            resname = ri.GetResidueName()
            atn = atname.replace(" ", "")
            #if resname != 'LIG' and atn =='N' or atn =='CA' or atn == 'C' or atn == 'O':
            if resname != 'LIG' and atn =='CA':
                probe_atomMap.append(atom.GetIdx())
    
        for atom in mol2.GetAtoms():
            ri = atom.GetPDBResidueInfo()
            atname = ri.GetName()
            resname = ri.GetResidueName()
            atn = atname.replace(" ", "")
            #if resname != 'LIG' and atn =='N' or atn =='CA' or atn == 'C' or atn == 'O':
            if resname != 'LIG' and atn =='CA' :
                refer_atomMap.append(atom.GetIdx())

    # Align mol1(probe) to mol2(reference)

        atomMap=list(zip(probe_atomMap,refer_atomMap))
        rms=rdMolAlign.AlignMol(mol1, mol2,atomMap=atomMap)

        file_path_align = f'triplet_db/tmp-aligned-FRASEdb_triplet_{pdbname}_{index}_{triplet_target}.pdb'
        Chem.MolToPDBFile(mol1, file_path_align)
        fi = open(file_path_align,"r")
        fragments_aligned=f'triplet_db/aligned-FRASEdb_triplet_{pdbname}_{index}_{triplet_target}.pdb'
        with open(fragments_aligned, 'a') as f:
            for line in fi.readlines():
                if (line.split()[0]) == 'HETATM' and (line.split()[3]) == 'LIG'  or line.split()[0] == 'CONECT':
                    f.write(line)
        f.close()
        fi.close()
        os.remove(file_path_align)
    
        mol = Chem.MolFromPDBFile(fragments_aligned)    
        sdf_file = Chem.SDWriter(f'triplet_db/aligned-FRASEdb_triplet_{pdbname}_{index}_{triplet_target}.sdf')
        if mol is None:
            return  
        sdf_file.write(mol)
        sdf_file.close()
        if os.path.exists(file_path_db):
            os.remove(fragments_aligned)
            os.remove(file_path_db)
            #os.remove(file_path_target)
    return 


# In[8]:


def FindDecoy(triplet_target,pdbname,index):
    clash_threshold = 1.0 # !!! change to 1.0
    sdf_file = f'triplet_db/aligned-FRASEdb_triplet_{pdbname}_{index}_{triplet_target}.sdf'

    if not os.path.exists(sdf_file):
        sys.exit(1)  # Exit with an error code
    try:
        suppl = Chem.SDMolSupplier(sdf_file)
        if suppl is None:
            sys.exit(1)  
    except Exception as e:
        sys.exit(1)  

   # suppl = Chem.SDMolSupplier(sdf_file)
    mol = suppl[0]
    file_path_target=f'triplet_db/tmp_target_triplet_{pdbname}_{index}_{triplet_target}.pdb'
    tmol_triplet = Chem.MolFromPDBFile(file_path_target)
    
    if tmol_triplet is None:
        return
    
    reslist = []
    for atom in mol.GetAtoms():          
        position_i = mol.GetConformer().GetAtomPosition(atom.GetIdx())
        for idx, atom_idx in enumerate(tmol_triplet.GetAtoms()):   
            ri = atom_idx.GetPDBResidueInfo()
            position_j = tmol_triplet.GetConformer().GetAtomPosition(idx)
            rij = position_i.Distance(position_j)
            if rij <= clash_threshold:
                #print('Bond clashes met. Random distance of two atoms between fragment and target:',rij)
                return
            if rij > clash_threshold and rij < 5:
                reslist.append(triplet_target)
                
    os.remove(file_path_target)     
    
    return triplet_target


# In[9]:

def check_and_append_buried_atoms(atom, m, tmol, clash_threshold, reslist):
    position_i = m.GetConformer().GetAtomPosition(atom.GetIdx())

    for idx, atom_idx in enumerate(tmol.GetAtoms()):
        ri = atom_idx.GetPDBResidueInfo()
        position_j = tmol.GetConformer().GetAtomPosition(idx)
        rij = position_i.Distance(position_j)

        if clash_threshold < rij < 5.0:
            residue_number = f"{ri.GetResidueNumber():>4}" if len(str(ri.GetResidueNumber())) == 4 else f"{ri.GetResidueNumber():>3}"
            res = f"{ri.GetResidueName()} {ri.GetChainId()} {residue_number}"
            reslist.append(res)
    return reslist

def Find5ARes(triplet_target,tmol,target,pdbid,pdbname,index):
  
    sdf_file = f"triplet_db/aligned-FRASEdb_triplet_{pdbname}_{index}_{triplet_target}.sdf"
    suppl = Chem.SDMolSupplier(sdf_file)
    m = suppl[0]
    
    if m is None:
        return
    clash_threshold = 1.0
      
    reslist = []
    buried = []

    for atom in m.GetAtoms():          
        position_i = m.GetConformer().GetAtomPosition(atom.GetIdx())
        for idx, atom_idx in enumerate(tmol.GetAtoms()):   
            ri = atom_idx.GetPDBResidueInfo()
            position_j = tmol.GetConformer().GetAtomPosition(idx)
            rij = position_i.Distance(position_j)
            if rij <= clash_threshold:
                return
            if  rij < 5.0 and rij > clash_threshold:              
                buried.append(rij)

    if len(buried) >= 5:
        for atom in m.GetAtoms():
            check_and_append_buried_atoms(atom, m, tmol, clash_threshold, reslist)

    if not reslist:
        return
    # write ligand_fragment pdb file                
    pdb1 = f'newFragments/ligfrag_{pdbname}_{index}_{triplet_target}.pdb'
    Chem.MolToPDBFile(m,pdb1)

    #destination_folder = "newFragments/"
    #shutil.copy(new_sdfname, destination_folder)

    # write residue contact pdb file
    for i in set(reslist):
        output1 = target
        output2 = f'newFragments/surr_res_{pdbname}_{index}_{triplet_target}.pdb'
        fp1 = open(output1,"r")
        fp2 = open(output2,"a")
        for line in fp1.readlines():
            if f"{i}" in line:
                fp2.write(line)
        fp2.close()
        fp1.close()
    # output ligand + surrunding residues
    output3 = f'newFragments/lig_surr_res_{pdbname}_{index}_{triplet_target}.pdb'

    fcom1 = open(pdb1,"r")
    fcom2 = open(output2,"r")
    fcom = open(output3,"a")
    
    smi = Chem.MolToSmiles(m)
    startline = pdbid + '\t' + smi
    fcom.write(startline)
    fcom.write('\r') # start a new line
    for line in fcom1.readlines():
        if "HETATM" in line:
            sline = line.replace('UNL','LIG')
            fcom.write(sline)
    fcom1.close()
    for line in fcom2.readlines():
        if "ATOM" in line:
            fcom.write(line)
    fcom2.close()
    fcom1 = open(pdb1,"r") # pdb1: ligand bond orders
    for line in fcom1.readlines():
        if "CONECT" in line:
            fcom.write(line)
    fcom1.close()
    fcom.close()
    os.remove(output2)
    
    # optimize writing info:
    mol = Chem.MolFromPDBFile(output3)
    if mol is None:
        os.remove(output3)
        return
    hash_value = hash(smi)
    
    atname = []
    riname = []
    for atom in mol.GetAtoms():
        ri = atom.GetPDBResidueInfo()
        if ri.GetResidueName() == 'LIG':
            atname.append(ri.GetResidueName())
        if ri.GetResidueName() != 'LIG':
            atname.append(ri.GetName())
            
    output4 = str(f'newFragments/tmp0_frase_in_target_{pdbname}_{index}_{triplet_target}.sdf')
    writer = Chem.SDWriter(output4)
    mol.SetProp('PDB_ID','%s' %(pdbid))
    mol.SetProp('FRASE_ID', pdbid + '_' + str(hash_value))
    mol.SetProp('_Name',pdbid)
    for i,j in enumerate(atname):
        mol.GetAtomWithIdx(i).SetProp('molFileAlias', 'V' + "%5d"%(i+1) + ' ' + j)
    writer.write(mol)
    writer.close()
    
    output5 = str(f'newFragments/frase_in_target_{pdbname}_{index}_{triplet_target}.sdf')
    fp4 = open(output4,"r")
    fp5 = open(output5,"w")
    for line in fp4.readlines():
        if not (line.startswith('A  ')):
            fp5.write(line)
    fp5.close()
    fp4.close()
    #print(output5)

    # remove intermediate files    
    os.remove(output3)
    os.remove(output4)
                
    return output5, ' was used for generating new FRASE in a target:'

# In [ ]: create directories: newFragments/ and triplet_db/ 
folder_paths = [
    "newFragments",
    "triplet_db"
]

try:
    for folder_path in folder_paths:
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
            print(f"Folder '{folder_path}' created.")
        else:
            print(f"Folder '{folder_path}' already exists. Skipping.")

except FileExistsError:
    print(f"The directory '{folder_path}' already exists.")

# In[ ]:

# Create a temporary file to capture warnings
#temp_stderr = sys.stderr
#with open("temp_stderr.txt", "w") as temp_stderr_file:
#    sys.stderr = temp_stderr_file
    
if len(sys.argv) < 3:
    print("Usage: python FRASE_screen.py frase_db_path pdbid")
    sys.exit(1)

frase_inp = sys.argv[1]
folder_id = sys.argv[2]

frase_db = f'../FRASE_database_PDB{folder_id}/{frase_inp}'
pdb_files = glob.glob(frase_db)
frase_inp_index = frase_db.split('/').index(frase_inp)
pdbname = [i.split('/')[frase_inp_index].split('.')[0] for i in pdb_files]
pmols = [Chem.MolFromPDBFile(x) for x in pdb_files if x is not None]
pdbid = [line.split()[0] for pdb_file_path in pdb_files for line in open(pdb_file_path, "r").readlines()[:1]]
print(len(pdb_files), 'FRASEs were used for matching the triplets in a target')

"""read from target protein"""
target = sys.argv[3]
tmol = Chem.MolFromPDBFile(target)

IF_target = []
pdblabel = []
for x in range(len(pmols)):   
    p3=ThreeAlphaCarbon(pmols[x])
    t3=ThreeAlphaCarbon(tmol)
    if p3 == []:
        print(f'There is no triplet matches found in FRASE {x+1}:',pdb_files[x], 'Skipping...')
        continue
    elif t3 == []:
        print(f'There is no triplet matches found in target, Skipping...')
        continue
    else:

        if pmols[x] is None:
            continue
        print('Triplet matches found in FRASE',x+1,'-',pdb_files[x], 'Processing...')
        for i in fp_target(tmol,t3):
            for fp_db6 in fp_database(pmols[x],p3):
                if i[0] in fp_db6[0]:
                    """start to align""" 
                    ## 1. write triplets from target protein and FRASE database to the pdb files seperately
                    WriteTripletTarget(i[1],target,pdbname[x],x)
                    WriteTripletFrasedb(fp_db6[1],pdb_files[x],pdbname[x],x)
                    ## 2. use the pdb files generated from above step to do alignments, 
                    ##    save the aligned ligand fragment to a single pdb file 
                    Alignment(fp_db6[1],i[1],pdbname[x],x)
                    ## 3. check structural decoy and find the fragment surrounding residues in target protein
                    tri_suv = FindDecoy(i[1],pdbname[x],x)
                    if tri_suv is not None:
                        outp5 = Find5ARes(tri_suv,tmol,target,pdbid[x],pdbname[x],x)
                        ## 4. calculate interaction fingerprint: 264 descriptors each LSE
                        if outp5 is not None:
                            INF = InteractionFingerprint(pdbname[x],tri_suv,x)
                            IF_target.append(INF)
                            na = str(pdbname[x])+'_'+str(x)
                            pdblabel.append(na)
                            #print(INF)
                            SDFFile = str(f'newFragments/frase_in_target_{pdbname[x]}_{x}_{tri_suv}.sdf')
                            os.remove(SDFFile)                           
                    break


if not IF_target:
    print('No interaction fingerprint was calculated. Exit...')
else:
    print('New FRASE was found and interaction fingerprint was shown:')
    for i in range(len(pdblabel)):
        print(f'{pdblabel[i]} :: {IF_target[i]}')

# Restore the original stderr
#sys.stderr = temp_stderr
# Close the temporary stderr file
#temp_stderr_file.close()
# Delete the temporary stderr file
#os.remove("temp_stderr.txt")
# Delete triplet_db folder
#shutil.rmtree(folder_paths[1])
