import sys
import os
import subprocess

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def gen_hostfile(P, ppn):
    os.system('rm -f hostfile')
    os.system('touch hostfile')
    f = open('hostfile', "w+")
    f1 = open("nodefile.txt", "r")
    lines = f1.readlines()
    total_groups = len(lines)
    nodes_per_group = [0]

    for line in lines:
        flag = 0
        nodes = line.strip().split(",")
        for node in nodes:
            status = subprocess.call(["ssh", node, "uptime"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if status is 0:
                f.write(node + ":" + str(ppn) + "\n")
                flag = 1
                P = P-1
                nodes_per_group[-1] = nodes_per_group[-1] + 1
            if P is 0:
                f.close()
                f1.close()
                return nodes_per_group
        if flag:
            nodes_per_group.append(0)
    if P is not 0:
        eprint("Required number of nodes not available...terminating")
        os.system('rm -f data.txt')
        exit(-1)

os.system('make clean')
os.system('make')
sys.stdout = open('data.txt', 'w')
for execution in range(10):
    eprint("Execution " + str(execution+1) + " of of 10 starting...")
    for P in [4, 16]:
        for ppn in [1, 8]:
            nodes_per_group = gen_hostfile(P, ppn)
            for D in [16, 256, 2048]:
                eprint("    Executing for P : " + str(P) + ", ppn : " + str(ppn) + ", D : " + str(D))
                string = "mpirun -np " + str(P*ppn) + " -f hostfile ./collectives " + str(D) + " " + str(ppn)
                for nodes_in_this_group in nodes_per_group:
                    string += " " + str(nodes_in_this_group)
                cmd = os.popen(string)
                print(cmd.read())
                cmd.close()

sys.stdout.close()
data = open("data.txt", "r")
lines = data.readlines()
f1 = open('data_Bcast.txt', "w")
f2 = open('data_Reduce.txt', "w")
f3 = open('data_Gather.txt', "w")
f4 = open('data_Alltoallv.txt', "w")

i = 0
for line in lines:
    if i%4 is 0:
        f1.write(line)
    elif i%4 is 1:
        f2.write(line)
    elif i%4 is 2:
        f3.write(line)
    elif i%4 is 3:
        f4.write(line)
    i += 1

f1.close()
f2.close()
f3.close()
f4.close()

os.remove('data.txt')