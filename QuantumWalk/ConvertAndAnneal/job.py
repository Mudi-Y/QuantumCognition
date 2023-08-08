number = "5"
path = "/workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/"
annealMethods = ["simulated","QPU"]
groups=[2,3,4,5,6,7,8]

with (open(f'{path}joblist.txt', 'w') as f):
    for annealMethod in annealMethods:
            f.write("#!/bin/sh\n")
            f.write(f"mkdir -p {path}Results/{annealMethod}{number}; ")
            for group in groups:
                fname = path+f"Results/{annealMethod}{number}/{annealMethod}_group{group}_BQM"
                command_str = f"/usr/local/bin/python /workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/ConstructBQM.py {group} > {fname}; "
                f.write(command_str)
                for i in [5,6,7,8,9,10]:
                    fname = path+f"Results/{annealMethod}{number}/{annealMethod}_group{group}_attempt{i+1}"
                    command_str = f"/usr/local/bin/python /workspace/QuantumCognition/QuantumWalk/ConvertAndAnneal/anneal.py {annealMethod} > {fname}; "
                    f.write(command_str)