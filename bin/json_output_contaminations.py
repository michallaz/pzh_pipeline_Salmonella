import json
import sys

# Simple script to preapare json output and save it to a file

pipeline_type = sys.argv[1] # decyduje ktore programy beda wyswietlane np. bacteria_illumina ma wyniki kraken2, metaphlan i kmerfinder, bacterial_nanopore nie ma metaphlan, dla viral jest tylko kraken2
output_name = sys.argv[2] # nazwa pliku z jsonem
initial_status = sys.argv[3]

# json dla kontaminacji ma postac listy

# parsowanie kraken2 pliku o STALEJ nazwie report_kraken2.txt
try:
    f=open("report_kraken2.txt").readlines()
    if len(f) == 0:
        if initial_status == "nie":
            kraken2_json = [{"program_name": "kraken2", \
                    "status": "nie", \
                    "error_message": "Kraken2 module was entered with failed QC status"}]
        elif initial_status == "blad":
            kraken2_json = [{"program_name": "kraken2", \
                    "status": "blad", \
                    "error_message": "Sample failed QC analysis with kraken2"}]
    else:
        genus_dict = {}
        species_dict = {}
        for line in f:
            line = line.split()
            if line[3] == "S":
                species_dict[line[5] + " " + line[6]] = float(line[0])
            elif line[3] == "G":
                genus_dict[line[5]] = float(line[0])
        genus_names_sorted = sorted(genus_dict, key = lambda x: genus_dict[x], reverse = True)
        species_names_sorted = sorted(species_dict, key = lambda x: species_dict[x], reverse = True)    
        kraken2_json =  [{"program_name": "kraken2", \
                "status": "Tak", \
                "main_genus_name": genus_names_sorted[0], \
                "secondary_genus_name": genus_names_sorted[1], \
              "main_species_name" : species_names_sorted[0], \
              "secondary_species_name": species_names_sorted[1], \
              "main_genus_value": f'{genus_dict[genus_names_sorted[0]]:.2f}', \
              "secondary_genus_value" : f'{genus_dict[genus_names_sorted[1]]:.2f}', \
              "main_species_value" : f'{species_dict[species_names_sorted[0]]:.2f}', \
              "secondary_species_value":  f'{species_dict[species_names_sorted[1]]:.2f}'}]
except FileNotFoundError:
    # Nie ma pliku zwracany ustalonego jsonai
    kraken2_json = [{"program_name": "kraken2", \
         "status": "blad", \
         "error_message": "No kraken2 output"}]

### metaphlan
try:
    f1=open("report_metaphlan_genera.txt").readlines()
    f2=open("report_metaphlan_species.txt").readlines()
    if len(f1) == 0 or len(f2) == 0:
        if initial_status == "nie":
            metaphlan_json = [{"program_name": "metaphlan", \
                    "status": "nie", \
                    "error_message": "Metaphlan module was entered with failed QC status"}]
        elif initial_status == "blad":
            metaphlan_json = [{"program_name": "metaphlan", \
                    "status": "blad", \
                    "error_message": "Sample failed QC analysis with metaphlan"}]

    else:
        genus_dict = {}
        species_dict = {}
        for line in f1:
            line=line.split()
            if "g__" in line[0]:
                genus_dict[line[0].split('_')[-1]] = float(line[2])
        for line in f2:
            line=line.split()
            if "s__" in line[0]:
                species_dict[" ".join(line[0].split('_')[-2:])] = float(line[2])

        genus_names_sorted = sorted(genus_dict, key = lambda x: genus_dict[x], reverse = True)
        species_names_sorted = sorted(species_dict, key = lambda x: species_dict[x], reverse = True)
        if len(genus_names_sorted) == 1:
            genus_names_sorted.append('None')
            genus_dict['None'] = 0
        if len(species_names_sorted) == 1:
            species_names_sorted.append('None')
            species_dict['None'] = 0

        metaphlan_json =  [{"program_name": "metaphlan", \
                "status": "Tak", \
              "main_genus_name": genus_names_sorted[0], \
              "secondary_genus_name": genus_names_sorted[1], \
              "main_species_name" : species_names_sorted[0], \
              "secondary_species_name": species_names_sorted[1], \
              "main_genus_value": f'{genus_dict[genus_names_sorted[0]]:.2f}', \
              "secondary_genus_value" : f'{genus_dict[genus_names_sorted[1]]:.2f}', \
              "main_species_value" : f'{species_dict[species_names_sorted[0]]:.2f}', \
              "secondary_species_value":  f'{species_dict[species_names_sorted[1]]:.2f}'}]
except FileNotFoundError:
    # Nie ma pliku zwracany ustalonego jsonai
    metaphlan_json = [{"program_name": "metaphlan", \
         "status": "blad", \
         "error_message": "No metaphlan output"}]

### Kmerfinder
try:
    f=open("results.txt").readlines()
    if len(f) == 0:
        if initial_status == "nie":
            kmerfinder_json = [{"program_name": "kmerfinder", \
                    "status": "nie", \
                    "error_message": "kmerfinder module was entered with failed QC status"}]
        elif initial_status == "blad":
            kmerfinder_json = [{"program_name": "kmerfinder", \
                     "status": "blad", \
                     "error_message": "Sample failed QC analysis with kmerfinder"}]
    else:
        genus_dict = {}
        species_dict = {}
        for line in f:
            line = line.split("\t")
            if "#" in line[0]:
                continue
            rodzaj = line[14].split(' ')[0]
            gatunek = rodzaj + " " + line[14].split(' ')[1]
            pokrycie = float(line[10])
            
            if rodzaj in genus_dict.keys():
                genus_dict[rodzaj] = genus_dict[rodzaj] + pokrycie
            else:
                genus_dict[rodzaj] = pokrycie

            if gatunek in genus_dict.keys():
                species_dict[gatunek] = species_dict[gatunek] + pokrycie
            else:
                species_dict[gatunek] = pokrycie

        genus_names_sorted = sorted(genus_dict, key = lambda x: genus_dict[x], reverse = True)
        species_names_sorted = sorted(species_dict, key = lambda x: species_dict[x], reverse = True)
        if len(genus_names_sorted) == 1:
            genus_names_sorted.append('None')
            genus_dict['None'] = 0
        if len(species_names_sorted) == 1:
            species_names_sorted.append('None')
            species_dict['None'] = 0

        kmerfinder_json =  [{"program_name": "kmerfinder", \
                "status": "Tak", \
              "main_species_name" : species_names_sorted[0], \
              "secondary_species_name": species_names_sorted[1], \
              "main_species_coverage" : f'{species_dict[species_names_sorted[0]]:.2f}', \
              "secondary_species_coverage":  f'{species_dict[species_names_sorted[1]]:.2f}' }]
except FileNotFoundError:
    # Nie ma pliku zwracany ustalonego jsonai
    kmerfinder_json = [{"program_name": "kmerfinder", \
         "status": "blad", \
         "error_message": "No kmerfinder output"}]

if pipeline_type == 'bacterial_illumina':
    full_output = kmerfinder_json + metaphlan_json + kraken2_json


with open(output_name, 'w') as f:
    f.write(json.dumps(full_output))

