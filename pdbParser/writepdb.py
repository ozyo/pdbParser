#See COPYING for license 

def writeca(div,fil):
    fmt="ATOM  {:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}"
    with open(fil,'w') as fout:
        for i in div:
            i=[i.item() for i in i]
            fout.write(fmt.format(*i)+'\n')
