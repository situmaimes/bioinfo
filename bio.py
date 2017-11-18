"""
read the whole text and tranform to json file
use re module to get infomation
"""


import json
import codecs
import re
import os


def processModel(model):
    """
    split the module and process every line
    :param model: module string
    :return:
    """
    lines=model.split("\n")
    modelDict={}
    #processs the remark
    res=re.match(r"REMARK VINA RESULT:[ ]*([\w\.\-]*)[ ]*([\w\.\-]*)[ ]*([\w\.\-]*)",lines[0])
    modelDict["REMARK VINA RESULT"]=dict(x=res[1],y=res[2],z=res[3])
    res=re.match(r"REMARK[ ]*([\w\.\-]*) (.*):",lines[1])
    num=int(res[1])
    modelDict["REAMARKNUMS"]=res[1]
    modelDict["DESCRIPTION"]=res[2]+" "+re.match(r"REMARK  status: (.*)",lines[2])[1]
    modelDict["TORSDOF"]=re.match(r"TORSDOF ([\w]*)",lines[-2])[1]
    modelDict["REMARK"]={}
    for i in lines[3:(3+num)]:
        result=processRemark(i)
        modelDict["REMARK"][result[0]]=result[1]
    #process the atoms line and put it into the dict
    modelDict["ATOMS"]={}
    for i in lines[(3+num):]:
        if i.startswith("ATOM"):
            result = processAtom(i)
            modelDict["ATOMS"][result[0]] = result[1]
    #process atoms that are in the root and branch
    #the inOrOutDict is to define the atom in the branch or not
    #finally put the info to the dict as before
    inOrOutDict={}
    for i in lines[3+num:]:
        if i=="ROOT" or i.startswith("BRANCH"):
            inOrOutDict[i]=True
            modelDict[i]={}
        elif i=="ENDROOT" or i.startswith("ENDBRANCH"):
            inOrOutDict[i.replace("END","")]=False
        elif i.startswith("ATOM"):
            for j in inOrOutDict:
                if inOrOutDict[j]:
                    result=processAtom(i)
                    modelDict[j][result[0]] = result[1]

    return modelDict

def processRemark(remark):
    """
    use re module to match
    :param remark: remark string
    :return: a tuple the first is the remark id (stirng) and the second is the remark info (dict)
    """
    remarkDict={}
    res=re.match(r"REMARK[ ]*([0-9]*)[ ]*([A-Z]*)[ ]*?between atoms: (.*?)  and  (.*)",remark)
    remarkDict["status"]=res[2]
    remarkDict["betweenLeft"]=res[3]
    remarkDict["betweenRight"]=res[4]

    return (res[1],remarkDict)

def processAtom(atom):
    """

    :param atom: atom string
    :return: a tuple the first is the atom id (string) and the second is the info (dict)
    """
    AtomObj = re.compile(
        r"(ATOM|HETATM)[ ]*([\w]*)[ ]{1,2}([\w\*]{1,4})[ ]*([a-zA-Z\*]{3,4}) ([a-zA-Z])[ ]*([\w\*]{1,5}).[ ]*([\d\.\-]*)[ ]*([\d\.\-]*)[ ]*([\d\.\-]*)[ ]*([\d\.\-]*)[ ]*([\d\.\-]*)[ ]*([\w\.\-]*)[ ]*([A-Z]{1,2})")

    atomDict={}
    res=AtomObj.match(atom)
    atomDict['name']=res[1]
    atomDict["原子名称"]=res[3]
    atomDict["残基名称"]=res[4]
    atomDict['链']=res[5]
    atomDict["残基序列号"]=res[6]
    atomDict["x"]=res[7]
    atomDict['y']=res[8]
    atomDict['z']=res[9]
    atomDict["占有率"]=res[10]
    atomDict['温度']=res[11]
    atomDict['区段标识']=res[12]
    atomDict['元素符号']=res[13]

    return (res[2],atomDict)

def readModel(url):
    """
    this function is to read the module file
    the module file contains some modules starts with module and ends with end module
    split module and use process module to process

    :param url: module filename
    :return:

    """

    with open(url) as f:
        line=f.readline().strip()
        modelsDict={}
        while(line):
            if line.startswith("MODEL"):
                modelName=line
                modelsDict[modelName]=''
            elif line=="ENDMDL":
                pass
            else:
                modelsDict[modelName]=modelsDict[modelName]+line+"\n"
            line=f.readline().strip()
        for i in modelsDict:
            modelsDict[i]=processModel(modelsDict[i])

    #write the dict to lig_out.json
    with codecs.open("lig_out.json","w","utf-8") as f:
        f.write(json.dumps(modelsDict,ensure_ascii=False,indent=4))

def readAtom(url):
    """
    this function is to read the atom file
    the atom file contains atom line starts with ATOM and hetatm line starts with HETATM which is similar to atom
    process every line use processAtom function
    the function can process Atom line and Hetatm line and return a tuple

    :param url: Atom filename
    :return:
    """
    with open(url) as f:
        line=f.readline().strip()
        AtomDict={}
        while(line):
            if line.startswith("ATOM") or line.startswith("HETATM"):
                result = processAtom(line)
                AtomDict[result[0]] = result[1]
            else:
                pass
            line=f.readline().strip()

    #write to a json
    with codecs.open(os.path.splitext(url)[0]+".json","w","utf-8") as f:
        f.write(json.dumps(AtomDict,ensure_ascii=False,indent=4))


if __name__=="__main__":
    modelurl="lig_out.pdbqt"
    readModel(modelurl)

    atomurl="rep.pdb"
    readAtom(atomurl)
