import os
import sys

os.system('make clean')
os.system('make')
os.system('bash checkhosts.sh')

file = open("hostfile","r") 
data = file.read() 
hosts = data.split("\n")
if len(hosts)<8:
    os.system('make clean')
    print("Minimum number of hosts required to run the program are not available...exiting now")
    exit()

P = [16,36,49,64]
N = [256,1024,4096,16384,65536,262144,1048576]
itrs = 5
for p in range(len(P)):
    sys.stdout = open('data' + str(P[p]) + '.txt', 'w')
    for n in range(len(N)):
        for itr in range(itrs):
            string='mpiexec -np '+ str(P[p]) + ' -ppn '+ str(8) + ' -f hostfile ./halo_exchange '+ str(N[n]) + " 50"
            cmd = os.popen(string)
            print(cmd.read())
            cmd.close()
