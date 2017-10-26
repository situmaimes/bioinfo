import json
import codecs
import re 
compileObj=re.compile(r"ATOM[ ]*([\w]*)[ ]{1,2}([\w\*]{1,4})[ ]*([a-zA-Z\*]{3,4}) ([a-zA-Z])[ ]*([\w\*]{1,5}).[ ]*([\d\.\-]*)[ ]*([\d\.\-]*)[ ]*([\d\.\-]*)[ ]*([\d\.\-]*)[ ]*([\d\.\-]*)[ ]*([\w\.\-]*)[ ]*([A-Z]{1,2})")
def processModel(model):
    lines=model.split("\n")
    modelDict={}
    res=re.match(r"REMARK VINA RESULT:[ ]*([\w\.\-]*)[ ]*([\w\.\-]*)[ ]*([\w\.\-]*)",lines[0])
    modelDict["REMARK VINA RESULT"]=dict(x=res[1],y=res[2],z=res[3])
    res=re.match(r"REMARK[ ]*([\w\.\-]*) (.*):",lines[1])
    num=int(res[1])
    modelDict["REAMARKNums"]=res[1]
    modelDict["Description"]=res[2]+" "+re.match(r"REMARK  status: (.*)",lines[2])[1]
    modelDict["TORSDOF"]=re.match(r"TORSDOF ([\w]*)",lines[-2])[1]
    modelDict["Remark"]=[]
    for i in lines[3:(3+num)]:
        modelDict["Remark"].append(processRemark(i))
    modelDict["Atoms"]=[]
    for i in lines[(3+num):]:
        if i.startswith("ATOM"):
            modelDict["Atoms"].append(processAtom(i))
    inOrOutDict={}
    for i in lines[3+num:]:
        if i=="ROOT" or i.startswith("BRANCH"):
            inOrOutDict[i]=True
            modelDict[i]=[]
        elif i=="ENDROOT" or i.startswith("ENDBRANCH"):
            inOrOutDict[i.replace("END","")]=False
        elif i.startswith("ATOM"):
            for j in inOrOutDict:
                if inOrOutDict[j]:
                    modelDict[j].append(processAtom(i))
    return modelDict
def processRemark(remark):
    remarkDict={}
    res=re.match(r"REMARK[ ]*([0-9]*)[ ]*([A-Z]*)[ ]*?between atoms: (.*?)  and  (.*)",remark)
    remarkDict["id"]=res[1]
    remarkDict["status"]=res[2]
    remarkDict["betweenLeft"]=res[3]
    remarkDict["betweenRight"]=res[4]
    return remarkDict

def processAtom(atom):
    atomDict={}
    res=compileObj.match(atom)
    atomDict['id']=res[1]
    atomDict["原子名称"]=res[2]
    atomDict["残基名称"]=res[3]
    atomDict['链']=res[4]
    atomDict["残基序列号"]=res[5]
    atomDict["x"]=res[6]
    atomDict['y']=res[7]
    atomDict['z']=res[8]
    atomDict["占有率"]=res[9]
    atomDict['温度']=res[10]
    atomDict['区段标识']=res[11]
    atomDict['元素符号']=res[12]
    return atomDict
def read(url):
    with open(url) as f:
        line=f.readline().strip()
        modelsDict={}
        while(line):
            if line.startswith("MODEL"):
                modelName=line
                modelsDict[modelName]=''
            elif line=="ENDMDL":
                print(line)
            else:
                modelsDict[modelName]=modelsDict[modelName]+line+"\n"
            line=f.readline().strip()
        for i in modelsDict:
            modelsDict[i]=processModel(modelsDict[i])
        with codecs.open("pdb.json","w","utf-8") as f:
            f.write(json.dumps(modelsDict,ensure_ascii=False,indent=4))

if __name__=="__main__":
    url="test.pdb"
    read(url)
