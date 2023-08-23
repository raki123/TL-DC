import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

TIMEOUT = 300

def csv2rec(filename):
    return np.recfromtxt(filename, dtype=None, delimiter=',', names=True, encoding='utf-8')

#Switches
EFFICIENCY = True
PROBLOG = True
SMPROBLOG = True
MEU = True
MAP = True
CONCOM = True
WIDTHS = True

#colors
ASPMC_COLOR = "b"
PROBLOG_COLOR = "r"
PITA_COLOR = "m"
CLINGO_COLOR = "g"

#labels
TIME_LABEL = "runtime in seconds"
INSTANCES_LABEL = "number of instances solved"
XWIDTH_LABEL = "X-width"
XDWIDTH_LABEL = "X/D-width"
LABEL_SIZE = 12

files = [
     "results_nauty_600", "results_pd_600", "results_pd_nauty_600", "results_td_600"
]
names = {
    "results_nauty_600" : "Backtracking",
    "results_pd_600" : "Pathwidth",
    "results_pd_nauty_600" : "Pathwidth + Nauty",
    "results_td_600" : "Treewidth"
}
timeout = 300
vbs = [ 300 for i in range(100) ]
for file in files:
    with open(file, "r") as data_file:
        data = []
        for line in data_file:
            line = line.split()
            line = [int(x) for x in line if len(x) > 0]
            if len(line) == 3:
                data.append(line[2])
            else:
                data.append(timeout)
            i = len(data) - 1
            vbs[i] = min(vbs[i], data[i])
        data.sort()
        plt.plot(range(1, 101), data, label=names[file])

vbs.sort()
plt.plot(range(1, 101), vbs, label="Virtual Best")
plt.ylabel(TIME_LABEL, size = LABEL_SIZE)
plt.xlabel(INSTANCES_LABEL, size = LABEL_SIZE)
plt.title("one_pair", size = LABEL_SIZE)
axes = plt.gca()
axes.set_xlim([1,100])
axes.set_ylim([0,timeout])
axes.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.legend(loc="best", prop={'size': LABEL_SIZE})
plt.tight_layout()
plt.savefig()


files = [
     "all_results_pd_600", "all_results_pd_nauty_600", "all_results_td_600"
]
names = {
    "all_results_pd_600" : "Pathwidth",
    "all_results_pd_nauty_600" : "Pathwidth + Nauty",
    "all_results_td_600" : "Treewidth"
}
timeout = 300
vbs = [ 300 for i in range(50) ]
for file in files:
    with open(file, "r") as data_file:
        data = []
        for line in data_file:
            line = line.split()
            line = [int(x) for x in line if len(x) > 0]
            if len(line) == 3:
                data.append(line[2])
            else:
                data.append(timeout)
            i = len(data) - 1
            vbs[i] = min(vbs[i], data[i])
        data.sort()
        plt.plot(range(1, 51), data, label=names[file])

vbs.sort()
plt.plot(range(1, 51), vbs, label="Virtual Best")
plt.ylabel(TIME_LABEL, size = LABEL_SIZE)
plt.xlabel(INSTANCES_LABEL, size = LABEL_SIZE)
plt.title("all_pair", size = LABEL_SIZE)
axes = plt.gca()
axes.set_xlim([1,50])
axes.set_ylim([0,timeout])
axes.xaxis.set_major_locator(MaxNLocator(integer=True))
plt.legend(loc="best", prop={'size': LABEL_SIZE})
plt.tight_layout()
plt.savefig()
