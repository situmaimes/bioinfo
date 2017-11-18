import sys
from calculate import caldis,calangle



def main(argv):
    if len(argv)==1:
        print("Input main.py -h for help")
        return
    if argv[1].lower() == "-h":
        help()
    elif argv[1].lower() == "-d":
        try:
            modelAtomInfo=(argv[2],argv[3])
            atomInfo=(argv[4],argv[5])
            dis=caldis(modelAtomInfo, atomInfo)
            if dis == -1:
                print("please retry")
            else:
                print("the distance of the two atom is " +str(dis))
        except IndexError:
            print("Wrong")
            print("For help")
            help()
    elif argv[1].lower() == "-a":
        try:
            modelAtomInfo = (argv[2], argv[3])
            atom1Info = (argv[4], argv[5])
            atom2Info = (argv[6], argv[7])
            result=calangle(modelAtomInfo, atom1Info,atom2Info)
            if result==-1:
                print("please retry")
            else:
                print("the angles of the three atoms are " + str(calangle(modelAtomInfo, atom1Info,atom2Info)))
        except IndexError:
            print("Wrong")
            print("For help")
            help()
    else:
        print("Input main.py -h for help")


def help():
    print(
        """
        The main.py has two functions
        1)
        calculcate the distance of two atoms
        usage: python main.py -d ModelId AtomId AtomRes AtomSymbol
        ModelId         :in the lig_out.pdbqt, which model you want to use, please give the id of model
        AtomId          :in the model you choosed, which atom you want to use, please give the id of model
        AtomRes         :in the rep.pdb, atom of which residue you want to choose ,please describe the residue
                         by using the residue name and residue id  like "pro2"
        AtomSymbol      :please give the atom symbol to confirm atom you want to use
        finnaly it looks like 
        python main.py -d 2 1 pro2 ca
        Don't care letters capital or not
        
        output is distance
        
        2)
        calculcate the triangle of three atoms
        usage: python main.py -d ModelId AtomId Atom1Res Atom1Symbol Atom2Res Atom2Symbol
        ModelId         :in the lig_out.pdbqt, which model you want to use, please give the id of model
        AtomId          :in the model you choosed, which atom you want to use, please give the id of model
        Atom1Res         :in the rep.pdb, atom of which residue you want to choose ,please describe the residue
                         by using the residue name and residue id  like "pro2"
        Atom1Symbol     :please give the atom symbol to confirm atom you want to use
        Atom2Res        :as the same
        Atom2Symbol     :as the same
        finnaly it looks like 
        python main.py -a 2 1 pro2 ca pro9 cd
        don't care letters capital or not
        
        output is a dict like {'modelAtom': angle1 , 'Atom1': angle2, 'Atom2': angle3 }
        the three atoms and angles respectively take its point as the pinnacle
        
        
        happy using 
        thanks
        """
    )


if __name__=="__main__":
    main(sys.argv)

