"""
this file is to calculate the distance of two atoms and calculate angles of the triangle formed by three atoms
caldis1 and calangle1 use the json file formed by bio.py
however
caldis2 and calangle2 use the original pdb file

use numpy to calculate
the prepare function cost the mostly time of all
it cost about 0.025 s to read
if calculate task number is big ,to use prepare and caldis1 calangle1
otherwise,the caldis2 calangle2 are preferred
"""
# import codecs
# import json
import numpy as np
import time
from bio import processAtom
import re
# def prepare():
#     """
#     read the lig_out json file and rep json file and global them for using in other functions
#     this function cost the mostly time of all
#     :return:
#     """
#     global ModelDict
#     global AtomList
#     with codecs.open("lig_out.json", encoding="utf-8") as f:
#         ModelDict = json.loads(f.read())
#
#     with codecs.open("rep.json",encoding="utf-8") as f:
#         AtomList=json.loads(f.read())

def getresult(modelId):
    """
    get the remark vina result of specified model
    use the model dict
    :param modelId:the id of model
    :return:
    """
    Model = ModelDict["MODEL " + str(modelId)]
    return Model["REMARK VINA RESULT"]['x']

# def caldis1(modelId,modelAtomId,AtomId):
#     """
#     use the model dict and atom list to get point info
#     :param modelId: the id of model in the lig_out
#     :param modelAtomId: the id of atom in the specified model
#     :param AtomId: the id of atom in the rep
#     :return: distance of atoms,if wrong id is given ,return -1
#     """
#     try:
#         Atom1=ModelDict["MODEL "+str(modelId)]["ATOMS"][str(modelAtomId)]
#         Atom2=AtomList[str(AtomId)]
#     except KeyError:
#         print("not found")
#         return -1
#
#     dis=calDis(Atom1,Atom2)
#     return dis

def caldis(modelAtomInfo,AtomInfo):
    """
    use the original file to get point info
    because it, the read process need to decide if the atom in the model
    :param modelId: the id of model in the lig_out
    :param modelAtomId: the id of atom in the specified model
    :param AtomId: a tuple of the description of atom in the rep like ("pro2","CA") the
                    the first element to decide which residue and the second is to decide the CA atom
    :return: distance of atoms,if wrong id is given ,return -1
    """
    modelId, modelAtomId=modelAtomInfo
    AtomRes,AtomName=AtomInfo
    AtomName = AtomName.upper()
    AtomResName=AtomRes[0:3].upper()
    AtomResId=str(AtomRes[3])
    modelAtomLine=None
    AtomLine=None      #this two param is to store the line we need
    #get modelatomline
    with open("lig_out.pdbqt") as f:
        inline=False   #to decide atom in model or not
        line=f.readline().strip()
        while(line):
            if line.startswith("MODEL "+str(modelId)):
                inline=True
            if line.startswith("ATOM") and inline:
                #the 6:11 of line if id
                if line[6:11].strip()==str(modelAtomId):
                    modelAtomLine=line
                    break
            line=f.readline().strip()
    #get atomline
    with open("rep.pdb") as f:
        line=f.readline().strip()
        while(line):
            if line.startswith("ATOM"):
                if re.search(AtomName+r"[ ]*"+AtomResName+" A [ ]*"+AtomResId,line):
                    AtomLine = line
                    break
            line=f.readline().strip()
    if modelAtomLine and AtomLine:
        Atom1=processAtom(modelAtomLine)[1]
        Atom2=processAtom(AtomLine)[1]
    else:
        print("Not found")
        return -1

    dis=calDis(Atom1,Atom2)
    return dis

def calDis(Atom1,Atom2):
    """
    use numpy to calculate the distance
    :param Atom1: a dict contains info like points of atom1
    :param Atom2: a dict contains info like points of atom2
    :return: distance
    """
    A = np.array([Atom1['x'], Atom1['y'], Atom1['z']]).astype('float')
    B = np.array([Atom2['x'], Atom2['y'], Atom2['z']]).astype('float')
    d = np.sqrt(np.sum(np.power(A - B, 2)))
    return d

# def calangle1(modelId,modelAtomId,AtomId1,AtomId2):
#     """
#     use modle dict and atom list to get info
#     :param modelId: the id of model in the lig_out
#     :param modelAtomId: the id of atom in the specified model
#     :param AtomId1: the id of atom1 in the rep
#     :param AtomId2: the id of atom2 in the rep
#     :return: a dict form calAngel
#     """
#     Atom1 = ModelDict["MODEL " + str(modelId)]["ATOMS"][str(modelAtomId)]
#     Atom2 = AtomList[str(AtomId1)]
#     Atom3 = AtomList[str(AtomId2)]
#     angle=calAngle(Atom1,Atom2,Atom3)
#     return angle

def calangle(modelAtomInfo,Atom1Info,Atom2Info):
    """
    use original file to get info
    :param modelId: the id of model in the lig_out
    :param modelAtomId: the id of atom in the specified model
    :param Atom1Info: a tuple of the description of atom in the rep like ("pro2","CA") the
                    the first element to decide which residue and the second is to decide the CA atom
    :param Atom2Info: same
    :return: a dict from calAngel,if error ,return -1
    """
    modelId, modelAtomId=modelAtomInfo
    Atom1Res, Atom1Name = Atom1Info
    Atom1Name=Atom1Name.upper()
    Atom1ResName = Atom1Res[0:3].upper()
    Atom1ResId = str(Atom1Res[3])
    Atom2Res, Atom2Name = Atom2Info
    Atom2Name = Atom2Name.upper()
    Atom2ResName = Atom2Res[0:3].upper()
    Atom2ResId = str(Atom2Res[3])
    modelAtomLine=Atom1Line=Atom2Line=None
    with open("lig_out.pdbqt") as f:
        #decide if line is in the model or not
        inline=False
        line=f.readline().strip()
        while(line):
            if line.startswith("MODEL "+str(modelId)):
                inline=True
            if line.startswith("ATOM") and inline:
                if line[6:11].strip()==str(modelAtomId):
                    modelAtomLine=line
                    break
            line=f.readline().strip()

    with open("rep.pdb") as f:
        line=f.readline().strip()
        while(line):
            Atom1found=False
            Atom2found=False    #to decide if we have found the info
            if line.startswith("ATOM"):
                if re.search(Atom1Name+r"[ ]*"+Atom1ResName+" A [ ]*"+Atom1ResId,line):
                    Atom1Line = line
                    Atom1found=True
                if re.search(Atom2Name+r"[ ]*"+Atom2ResName+" A [ ]*"+Atom2ResId,line):
                    Atom2Line=line
                    Atom2found=True
            #both have been found and break
            if Atom1found and Atom2found:
                break
            line=f.readline().strip()
    if modelAtomLine and Atom1Line and Atom2Line:
        Atom1=processAtom(modelAtomLine)[1]
        Atom2=processAtom(Atom1Line)[1]
        Atom3=processAtom(Atom2Line)[1]
        return calAngle(Atom1,Atom2,Atom3)
    else:
        print("Not found")
        return -1



def calAngle(Atom1,Atom2,Atom3):
    """
    use numpy to calculate angle
    x,y,z is three vectors formed by three poins
    :param Atom1: info of atom1
    :param Atom2: same...
    :param Atom3: same...
    :return: a dict like {'modelAtom': angle1 , 'Atom1': angle2, 'Atom2': angle3 }
            the three angel respectively take there points as the pinnacle
    """
    A = np.array([Atom1['x'], Atom1['y'], Atom1['z']]).astype('float')
    B = np.array([Atom2['x'], Atom2['y'], Atom2['z']]).astype('float')
    C = np.array([Atom3['x'], Atom3['y'], Atom3['z']]).astype('float')
    x = B - A
    y = C - A
    z = C - B
    Lx = np.sqrt(x.dot(x))
    Ly = np.sqrt(y.dot(y))
    Lz = np.sqrt(z.dot(z))
    angleA = np.arccos(np.abs(x.dot(y)) / (Lx * Ly)) * 360 / 2 / np.pi
    angleB = np.arccos(np.abs(x.dot(z)) / (Lx * Lz)) * 360 / 2 / np.pi
    angleC = np.arccos(np.abs(y.dot(z)) / (Ly * Lz)) * 360 / 2 / np.pi
    return {"modelAtom": angleA, "Atom1": angleB, "Atom2": angleC}


if __name__=="__main__":
    starttime = time.clock()
    #prepare()
    print(caldis((2,10) ,("pro2","ca")))
    print(calangle((2,10), ("pro2","ca"), ("pro9","cd")))
    # print(getresult(1))
    # print(caldis1(2,10,50))
    # print(calangle1(2,10,50,2))
    endtime=time.clock()
    print(endtime-starttime)



