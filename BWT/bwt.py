def rotations(t):
    """ Return list of rotations of input string t """
    tt = t * 2
    return [tt[i:i + len(t)] for i in range(0, len(t))]


def bwm(t):
    """ Return lexicographically sorted list of t's rotations, i.e. bw matrix """
    return sorted(rotations(t))


def bwtViaBwm(t):
    """ Given T, return BWT(T) by way of BWM """
    return ''.join(map(lambda x: x[-1], bwm(t)))


print(bwtViaBwm("banana$"))
print(sorted(bwtViaBwm("banana$")))

bwt=list(bwtViaBwm("banana$"))
head=sorted(bwt)

twomer=[str(x[0]+x[1]) for x in zip(bwt,head)]
twomerS=sorted(bwt)
