import itertools, operator

class getmissing(object):
    def __init__(self):
        self.misres=[]
        self.misatom=[]
        self.ch_nr=None
        self.res_list=[]
        self.misres_sort_ch=[]
    def getmissing(self,pdb):
        missing_res=False
        missing_atm=False
        for line in pdb:
            if "REMARK" in line:
                if "465" in line:
                    missing_res=True
                else:
                    missing_res=False
                if missing_res is True:
                    try:
                        ch=line.split(' ',8)[7].strip()
                        res=line.split(' ',8)[8].strip()
                        if len(res) == 3:
                            tofix=(ch,res)
                            self.misres.append(tofix)
                    except IndexError:
                        pass
        self.misres_sort_ch=sorted((list(group) for key,group in itertools.groupby(self.misres,operator.itemgetter(0))),key=operator.itemgetter(2))
#A complex generator. The first argument argument in sorted groups the misssing residues list by its chain, and creates a sublist. This sublist is then passed to sort to be sorted based on the resid. This is to ensure we create residue ranges for modelling.
# In theory it should be sorted but you cannot be too careful when it comes to pdb files.
        for chain in self.misres_sort_ch:
            for pair in range(0,len(chain)):
                self.res_list.append(chain[pair][1])
        self.ch_nr=len(self.misres_sort_ch)

        return self.misres_sort_ch, self.res_list, self.ch_nr

if __name__ == "__main__":
    print "I am not a stand alone script"
