
#import collections
import sys

import matplotlib.pyplot as plt

def main(out, start=None, finish=None):
    first = True
    x = []
    single = []
    multiple = []
    unmapped = []
    for line in sys.stdin:
        if first:
            first = False
            continue
        fields = line.strip('\n').split('\t')
        pos = int(fields[0])
        if (start is None or pos >= start) and (finish is None or pos <= finish):
            x.append(pos)
            total = float(fields[4])
            single.append(int(fields[1]) / total * 100)
            multiple.append(int(fields[2]) / total * 100)
            unmapped.append(int(fields[3]) / total * 100)

    print "Plotting {0} items".format(len(x))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    # at the bottom in red - unmapped
    ax.plot(x, unmapped, label='Unmapped', color='r')
    ax.fill_between(x, unmapped, 0, color='r', alpha=0.5)

    # unmapped + multimapped in blue
    um = [ sum(y) for y in zip(unmapped, multiple) ]
    ax.plot(x, um, label='Multimapped', color='b')
    ax.fill_between(x, um, unmapped, color='b', alpha=0.5)

    # legend
    legend='upper right'
    leg = ax.legend( loc=legend, prop={'size':12})
    leg.get_frame().set_alpha(0.8)
    
    ax.set_ylabel('% references affected')
    ax.set_xlabel('Source reference position')

    # stop rounding off with the x-axis length
    plt.xlim(xmax=max(x), xmin=min(x))
    
    fig.savefig(out, format='pdf', dpi=1000)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Plot bias')
    parser.add_argument('--start', type=int, required=False, help='start position')
    parser.add_argument('--finish', type=int, required=False, help='finish position')
    parser.add_argument('--out', required=False, default="plot.pdf", help='finish position')
    args = parser.parse_args()
    main(args.out, args.start, args.finish)

