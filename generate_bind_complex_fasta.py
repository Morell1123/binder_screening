import os
import pandas as pd

# Darian's binders

binders_path = 'ProteinBO/data/data_DBL/Darians_data/Minibinder_fasta'

targets = {}
targets['target_A'] = ('GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTV'
                       'QRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMT'
                       'HHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEP:'
                       'MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPK'
                       'IVKWDRDM:SLLMWITQC')
targets['target_B'] = ('GSHSMRYFFTSVSRPGRGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDGETRKVKAHSQTHRVDLGTLRGYYNQSEAGSHTV'
                       'QRMYGCDVGSDWRFLRGYHQYAYDGKDYIALKEDLRSWTAADMAAQTTKHKWEAAHVAEQLRAYLEGTCVEWLRRYLENGKETLQRTDAPKTHMT'
                       'HHAVSDHEATLRCWALSFYPAEITLTWQRDGEDQTQDTELVETRPAGDGTFQKWAAVVVPSGQEQRYTCHVQHEGLPKPLTLRWEP:'
                       'MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPK'
                       'IVKWDRDM')

f = open('ProteinBO/data/data_DBL/Darians_data/binding_complexes_Darian.fasta', 'w')
for binder_file in os.listdir(binders_path):
    with open(os.path.join(binders_path, binder_file), 'r') as binder:
        lines = binder.readlines()
        binder_name = lines[0][1:-1]
        binder_seq = lines[1][:-1]
        print(binder_seq)
        print(binder_name)
        for target_name, target_seq in targets.items():
            f.write(f">binder_{binder_name}_target_{target_name}\n")
            f.write(f"{binder_seq}\n:{target_seq}\n")
f.close()
hi=2