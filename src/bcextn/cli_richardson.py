from .mpmath import mp, mpf
import matplotlib.pyplot as plt
def main():
    from .fresnelref import load_refpts, load_stepdata
    data_refpts = load_refpts()
    data_steps = load_stepdata()
    if set(data_refpts.keys()) != set(data_steps.keys()):
        raise SystemExit('step data out of sync. Please update with:'
                         ' ./bin/createfresnelrefdata --updatestepdata')
    #import pprint
    for key in data_refpts.keys():

        print(key)
        dref = data_refpts[key]
        #pprint.pprint(dref)
        dsteps = data_steps[key]
        yref = dref['results']['val']
        pts = dsteps['pts']
        sum_n = [mpf(0)]
        for e in pts:
            sum_n.append( sum_n[-1] + e[1] )
        sum_n = sum_n[1:]#remove leading 0 again

        nvals = []
        precvals = []
        for n in range(3,len(sum_n)):
            #if n!=40:
            #   continue
            #if n>80:
            #    break
            r = mp.richardson( sum_n[0:n] )[0] #fixme [1] gives error???
            nvals.append(n)
            precvals.append(float(abs(r-yref)/min(1-yref,yref)))
        plt.plot(nvals,precvals,label=str(key))
            #print(n,float(abs(r-yref)/yref))
    plt.semilogy()
    plt.grid()
    plt.legend()
    plt.show()
