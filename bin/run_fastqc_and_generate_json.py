#!/usr/bin/env python3
"""
Skrypt do parsowania wynikow programu fastqc i produkowania json-a zgodnego ze specyfikacja w repozytorium plepiseq_json
"""

import sys
import numpy as np
import click
import json
import subprocess
import re
def run_fastqc(plik, memory, cpu):
    polecenie = (f'fastqc --format fastq --threads ${cpu} --memory {memory} --extract --outdir . {plik}')
    proces = subprocess.Popen(polecenie, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proces.wait()
    outs, errs = proces.communicate()
    status_fastqc = re.findall('Analysis complete', str(outs))
    if len(status_fastqc) > 0:
        status = 'tak'
    else:
        status = 'blad'
    return status
@click.command()
@click.option('-i', '--input_file', help='[INPUT] a path to a file in fastq format ',
              type=click.Path(), default="fastqc_data.txt")
@click.option('-m', '--memory', help='[INPUT] Memory available to fastqc program',
              type=int, default=4048)
@click.option('-c', '--cpu', help='[INPUT] CPUs available to fastqc program',
              type=int, default=4)
@click.option('-s', '--status', help='[INPUT] PREDEFINED status that is transferred to an output json. '
                                     'If this status was either nie or blad fastqc will not run',
              type=click.Choice(['tak', 'nie', 'blad'], case_sensitive=False),  required=True)
@click.option('-e', '--stage', help='[INPUT] Stage on which data is analyzed',
              type=click.Choice(["pre-filtering", "post-filtering"], case_sensitive=False),  required=True)
@click.option('-p', '--publishdir', help='[INPUT] Path with fNEXTFLOW output, required to correctly format json',
              type=str,  required=True)
@click.option('-m', '--error', help='[INPUT] PREDEFINED error message that is put in json. '
                                    'Only used when status was set to nie or blad',
              type=str,  required=False, default="")
@click.option('-o', '--output', help='[Output] Name of a file with json output',
              type=str,  required=False)
def main_program(input_file, memory, cpu, status, stage, publishdir, output, error):
    # Na poczatku obsluzmy sytuacje gdzie predefiniowany status to blad lub nie
    if status == "nie" or status == 'blad':
        json_dict = [{"status": status, "file_name": input_file, "step_name": stage, "error_message": error}]
        with open(output, 'w') as f1:
            f1.write(json.dumps(json_dict))
        return status
    else:
        # upewnijmy sie ze podany plik wogole istnieje
        try:
            open(input_file)
        except FileNotFoundError:
            status = 'blad'
            json_dict = [{"status": {status}, "file_name": input_file, "step_name": stage,
                         "error_message": "Provided file does not exists"}]
            with open(output, 'w') as f1:
                f1.write(json.dumps(json_dict))
            return status
        # to wywoluje polecenie fastqc na pliku, przy okazji obsluzymy wyjatek ze podany plik nie jest
        # fastq wedlug fastqc
        status = run_fastqc(input_file, memory, cpu)
        if status == 'blad':
            json_dict = [{"status": "blad", "file_name": input_file, "step_name": stage,
                          "error_message": "Provided file is not in fastq format"}]
            with open(output, 'w') as f1:
                f1.write(json.dumps(json_dict))
            return status
        # mamy plik z wynikami i fastqc go przeprocesowal
        input_dir = input_file.replace('.fastq.gz', '_fastqc')
        number_of_reads_value = 0  # pole number_of_reads_value
        number_of_bases_value = 0  # pole number_of_bases_value
        reads_median_length_value = 0  # pole reads_median_length_value
        reads_min_length_value = 0  # reads_min_length_value
        reads_max_length_value = 0  # pole reads_max_length_value
        reads_median_quality_value = 0  # pole reads_median_quality_value
        reads_quality_histogram_path = f'{input_dir}_reads_quality_histogram.csv'  # pole reads_quality_histogram_path
        reads_length_histogram_path = f'{input_dir}_reads_length_histogram.csv'  # pole reads_length_histogram_path
        position_quality_plot_path = f'{input_dir}_position_quality_plot.csv'
        gc_content_value = 0  # OPCJONALNE pole gc_content_value

        reads_median_quality_data = []  # pomocnicza lista do obliczenia mediany
        reads_median_length_data = []  # pomocnicza lista do obliczenia mediany
        start = 0
        section_name = ""
        # otwieramy polaczenia do plikow z wymaganymi nazwami
        # tworzymy tez naglowki zgodnie z dokumentacja excelu
        qualty_histogram_file = open(reads_quality_histogram_path, 'w')
        qualty_histogram_file.write(f'#indeks;jakosc_wartosc;liczebnosc\n')

        length_histogram_file = open(reads_length_histogram_path, 'w')
        length_histogram_file.write(f'#indeks;dlugosc_wartosc;liczebnosc\n')

        position_quality_file = open(position_quality_plot_path, "w")
        position_quality_file.write(f'#indeks;pozycja_w_odczycie;mediana_jakosc\n')
        with open(f'{input_dir}/fastqc_data.txt') as f:
            for line in f:
                line = line.split('\t')
                if ">>" in line[0] and start == 0:
                    # Sekcje z rownymi typami wynikow w pliku rozpoczynaja sie od ">>"
                    # a koncza na ">>END_MODULE"
                    start = 1
                    section_name = line[0].replace(">>", "")
                elif ">>" in line[0] and start == 1:
                    start = 0
                    if section_name == "Per sequence quality scores":
                        reads_median_quality_value = np.median(reads_median_quality_data)
                    # koniec sekcji mozna zrzucac wyniki
                    if section_name == "Sequence Length Distribution":
                        reads_median_length_value = np.median(reads_median_length_data)

                elif start == 1:
                    # jestem s sekcji
                    if section_name == "Basic Statistics":
                        if line[0] == "Filename":
                            nazwa_pliku = line[1]
                        if line[0] == "Total Sequences":
                            number_of_reads_value = line[1]
                        if line[0] == "Sequence length":
                            reads_min_length_value, reads_max_length_value = line[1].split("-")
                        if line[0] == "Total Bases":
                            number_of_bases_value = float(line[1].split(" ")[0])
                            unit = line[1].split(" ")[1].rstrip()  # fastqc dodaje okreslecia jak Kbp, Mbp i pewnie Gbp
                            if unit == "Kbp":
                                number_of_bases_value = number_of_bases_value * 1000
                            elif unit == "Mbp":
                                number_of_bases_value = number_of_bases_value * 1000000
                            elif unit == "Gbp":
                                number_of_bases_value = number_of_bases_value * 1000000000
                        if line[0] == "%GC":
                            gc_content_value = line[1]

                    if section_name == "Per sequence quality scores":
                        if "#" in line[0]:
                            indeks = 0
                        else:
                            qualty_histogram_file.write(f'{indeks};{line[0]};{line[1]}')
                            reads_median_quality_data.extend([int(line[0])] * int(float(line[1].rstrip())))
                            indeks += 1
                    if section_name == "Sequence Length Distribution":
                        if "#" in line[0]:
                            indeks = 0
                        else:
                            # fastqc przekazuje dane w postaci od-do np 30-35,
                            length_histogram_file.write(f'{indeks};{line[0].split("-")[0]};{line[1]}')
                            reads_median_length_data.extend([int(line[0].split("-")[0])] * int(float(line[1].rstrip())))
                            indeks += 1
                    if section_name == "Per base sequence quality":
                        if "#" in line[0]:
                            indeks = 0
                        else:
                            if len(line[0].split("-")) == 2:
                                pozycja = line[0].split("-")[1]
                            elif len(line[0].split("-")) == 1:
                                pozycja = line[0].split("-")[0]
                            else:
                                json_dict = [{"status": "error", "file_name": input_file, "step_name": stage,
                                             "error_message": "Error when parsin fastqc file"}]
                                with open(output, 'w') as f1:
                                    f1.write(json.dumps(json_dict))
                                return status

                            position_quality_file.write(f"{indeks};{pozycja};{line[2]}\n")
    json_dict = [{"status": "tak",
                 "file_name": nazwa_pliku.rstrip(),
                 "step_name": stage,
                 "error_message": "",
                 "number_of_reads_value": int(number_of_reads_value),
                 "number_of_bases_value": int(number_of_bases_value),
                 "reads_median_length_value": round(float(reads_median_length_value),2),
                 "reads_min_length_value": int(reads_min_length_value),
                 "reads_max_length_value": int(reads_max_length_value),
                 "reads_median_quality_value": round(float(reads_median_quality_value),2),
                 "reads_quality_histogram_path": f"{publishdir}/{reads_quality_histogram_path}",
                 "reads_length_histogram_path": f"{publishdir}/{reads_length_histogram_path}",
                 "position_quality_plot_path": f"{publishdir}/{position_quality_plot_path}",
                 "gc_content_value": round(float(gc_content_value), 2)}]
    with open(output, 'w') as f1:
        f1.write(json.dumps(json_dict))
    return status


if __name__ == '__main__':
    if len(sys.argv) == 1:
        main_program(['--help'])
    else:
        main_program(sys.argv[1:])
