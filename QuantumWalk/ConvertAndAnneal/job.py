number = "3"
path = "/workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/"
annealMethods = ["QPU"]
groups=[8]

with (open(f'{path}joblist.txt', 'w') as f):
    for annealMethod in annealMethods:
            f.write("#!/bin/sh\n")
            f.write(f"mkdir -p {path}Results/{annealMethod}{number}; ")
            for group in groups:
                # fname = path+f"Results/{annealMethod}{number}/{annealMethod}_group{group}_BQM"
                # command_str = f"/usr/local/bin/python /workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/ConstructBQM.py {group} > {fname}; "
                # f.write(command_str)
                for i in range(20):
                    fname = path+f"Results/{annealMethod}{number}/{annealMethod}_group{group}_attempt{i+1}"
                    command_str = f"/usr/local/bin/python /workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/anneal.py {annealMethod} > {fname}; "
                    f.write(command_str)