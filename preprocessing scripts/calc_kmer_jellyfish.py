#calculate k-mers for every transcript
#the tool jellyfish is utilised here

import os
from subprocess import call

os.chdir("../training_1/sno/")
path = os.getcwd()
files = [f for f in os.listdir(path) if os.path.isfile(f)]

for k in range(2,8):
    #with open (os.path.abspath(os.path.join("Sequences", "jelly"+str(k)+".sh")), "w") as w:
    os.system("pwd")
    if not os.path.isdir(os.path.join("kmer"+str(k))):
        os.makedirs(os.path.join("kmer"+str(k)))
        print (k)
    for i in files:
        os.system("jellyfish-linux count -m %d -s 100M -t 10 -C %s -o kmer%d/%s.jf"
                   %(k, i.rstrip("\n"), k, i.rstrip("\n")))       
        call("jellyfish-linux dump -c kmer"+str(k)+"/"+i.rstrip("\n")+".jf_0 > kmer"+str(k)+"/"+i.rstrip("\n")+".count", shell=True)
        

