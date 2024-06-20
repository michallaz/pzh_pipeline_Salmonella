import os
import re
from Bio import SeqIO

def create_profile(profile_file):
    """
    funkcja do tworzenia slownika ST na podstawie pliku profile sciagnietego z repozytorium enterobase
    :param profile_file: String do pliku profiles, rozne profile (7MLST, cgMLSDT, wgMLST) zawiera rozna liczbe kolumn,
    piwerwsza kolumna to ST, koejne kolumny to informacja jaki wariant w danym miejscu genomu wystepuje (w postaci liczby),
    jesli danego wariantu w danym ST nie ma jest "-".  Pierwszy wiersz to naglowki
    :return: slownik slownikow, jako nadrzedny klucz jest sequence type, kazdemu sequence type odpowiada slownik zbudowany
    z kluczy (nazwy genow z pierwszego naglowka), i jedna wartosc (numerek z pliku)
    """
    slownik_profili = {}
    if not os.path.isfile(profile_file):
        raise IOError('Brak sciezki do pliku z profilem')
    else:
        with open(profile_file) as f:
            for line in f:
                if re.search('ST', line):
                    # jestesmy w wierszy naglowkowym tworzymy liste nazw, ktora posluzy jako klucze w slowniku
                    lista_kluczy=line.rsplit()[1:] # ST nie jest kluczem
                else:
                    elementy = line.rsplit()
                    # tworzymy slownik dla danego ST
                    slownik_profili[elementy[0]] = {}
                    for indeks,wartosc in enumerate(elementy[1:]):
                        slownik_profili[elementy[0]][lista_kluczy[indeks]] = wartosc
    return slownik_profili, lista_kluczy

def getST(MLSTout, profile_dict, lista_kluczy):
    """
    Funkcja do porownywania dwoch slownikow. Oba slowniki musza miec identyczne klucze (w tym wypadku odpowiadajace
    loci z danego MLST). MLSTout to zwykly slownik,na tomiast profile_dict to slownik slownikow wynik funkcji
    create_profile. Jeśli nie znajdziemy ST dla slownika MLSToit , to 1. tworzymy nowy ST dla niego (int, +1 w stosunku
    do największego ST obecnego w profile_dict); 2. updatujemy plik profile_file o nowy ST.
    :param MLSTout: slownik
    :param profile_dict: slownik wygenerowany create_profile
    :param profile_file: String do pliku profiles
    :param lista_kluczy: list, lsta posortowanych alleli w danym schemacie
    :return: string, sequence type
    """


    ST = '0'
    slownik_roznic = {} #slownik ktory jako klucze ma nazwe profilu, jako wartosc ilosc roznic do "naszej" probki
    lista_probki = [MLSTout[x] for x in lista_kluczy]
    for ST_dict in profile_dict:
        slownik_roznic[ST_dict] = 0
        lista_ST = [profile_dict[ST_dict][x] for x in lista_kluczy]
        slownik_roznic[ST_dict] = sum(1 for x, y in zip(lista_probki, lista_ST) if x != y)


    if 0 in slownik_roznic.values():
        ST = [x for x, y in slownik_roznic.items() if y == 0][0]
        return ST
    elif 1 in slownik_roznic.values():
        ST = [x for x, y in slownik_roznic.items() if y == 0]
        # ST is a list, we nned to write all possible
        return ST
    else:
        return 0
        # We have smt completly new for 7MLST we better think what to do here ...

    # Na razie nie ma potrzeby tego uzywac
    # dodac oddzielna funkcje jesli znaleziono nowy schemat aby updatowac plik ze schematami ?
    # ale wtedy ten plik musi byc na ZEWNATRZ kontenera ...

    # if ST == 0:
    #
    #     # nowa kombinacja alleli
    #     # sprawdzamy najblizszy ST dla ciekawosci
    #     identity = 0
    #     for klucz in  profile_dict[ST_dict]:
    #         if profile_dict[ST_dict][klucz] == MLSTout[klucz]:
    #             identity +=1
    #     slownik_tmp_identity[ST_dict] = identity
    #
    #     # Na koniec nadajemy nowy St i dopisujemy poprawna kombinacje alleli do pliku profile
    #     ST = max(map(int, profile_dict.keys())) + 1
    #     with open(profile_file, 'a') as f:
    #         msg = "\t".join([str(MLSTout[x]) for x in lista_kluczy])
    #         f.write(f'{ST}\t{msg}\n')
    # #tworzymy plik wsadowy do hiercc
    #
    # with open('to_hiercc.txt', 'w') as f:
    #     msg = [str(ST)] + [str(MLSTout[x]) for x in lista_kluczy]
    #     lista_kluczy = ['ST'] + lista_kluczy
    #     f.write('{naglowek}\n{body}'.format(naglowek = "\t".join(lista_kluczy), body = "\t".join(msg)))
    #
    #
    # return ST, slownik_tmp_identity
def update_MLSTprofile_file():
    pass

def parse_MLST_tsv(file_path, long = True, sep = "\t"):
    """
    Funkcja do tworzenia slownika MLST na podstawie dwulinijkowego pliku, w ktorym naglowek ma klucze (odpowiadajace
    allelom) a drugiw wiersz to odzzielony ("\t", ";", " ", ",") identyfikator allelu
    :param file_path: string, sciezka do pliku z wynikami z enterobase, chewbacca
    :param long: bool, czy plik zawiera jako pierwsze 2 wpisy identyfikator sekwencji i Sequence type (True), czy tylko
    identyfikator (False)
    :param sep: string, separator kolumn w pliku poanych jako argument file_path
    :return: Funkcja zwraca slownik gdzie luczami sa nazwy alleli a wartosciami wariant
    """
    if long:
        zakres = 2
    else:
        zakres = 1
    slownik_alleli= {}
    with open(file_path, 'r') as f:
        for line in f:
            if re.search('ST', line):
                klucze = line.rsplit(sep)[zakres:]
                klucze = [klucz.split('.')[0] for klucz in klucze]
            else:
                elementy = line.rsplit(sep)[zakres:]
                # tworzymy slownik dla danego ST
                for indeks, wartosc in enumerate(elementy):
                    slownik_alleli[klucze[indeks]] = wartosc

    return slownik_alleli

def parse_MLST_fasta(file_path):
    """
    Funkcja do parsowania wynikow etoki zwracanych jako fasta
    :param file_path: string, sciezka do pliku z wynikami etoki, plik w formacie fasta
    :return: slownik, kluczami sa allele a wartosciami ich wariant
    """
    slownik_alleli = {}
    for record in SeqIO.parse(file_path, "fasta"):
        # opis z etoki jest dosc wystandaryzowany wiec zakladam ze pole description[0] zawiera nazwe wariantu
        # pole 2 zawiera id = wartosc
        # pole 6 zwiera identycznosc sekwencyjna miedzy referencja a tym co jest obserwowane w probce
        slownik_alleli[record.description.split(' ')[0]] = record.description.split(' ')[2].split('=')[1]
    return slownik_alleli

def parse_MLST_blastn(plik):
    """
    Funkcja do parsowania outputu mojego skryptu blastowego w ktorym sa 4 kolumny, pierwsza to nazwa allelu, druga numer
    allelu, 3 to identycznosc sekwenyjna miedzy tym allelem a sekwencja w genomie, 4 to info czy bylo wiele mapowan
    :param plik: str, sciezka do pliku z wynikami funkcji run_blastn_verX.sh
    :return: slownik, kluczami sa nazwy alleli w postaci stringow, warosciami numery alleli w postaci int
    """
    slownik_profili = {}
    if not os.path.isfile(plik):
        raise IOError('Brak sciezki do pliku z profilem')
    with open(plik, 'r') as f:
        for line in f:
            line = line.rsplit()
            try:
                int(line[2])
                if int(line[2]) == 100:
                    slownik_profili[line[0]] = int(line[1])
                else:
                    slownik_profili[line[0]] = -1
            except ValueError:
                slownik_profili[line[0]] = -1
    return(slownik_profili)

def compare_2_allel_dict(slownik1, slownik2, file_prefix = 'out'):
    """
    Funkcja do porownywania dwoch slownikow alleli w celu okreslenia ich zgodnosci. Funkcja zwraca ilosc kluczy, ktore
    w obu slownikach przyjmuja takie same wartosci. Ponadtow dw a pliki z unikalnymi kluczami dla obu
    slownikow rowniez zapisywane sa do pliko {file_prefix}_slownik1_uniquekeys.txt i {file_prefix}_slownik2_uniquekeys.txt

    :param slownik1: slownik
    :param slownik2: slownik, o kluczach identycznych z tymi w slowniku 1
    :param file_prefix: string, funkcja zwraca dwa pliki. Kazdy z nich ma 3 kolumny (nazwa alellu, wersja w slowniku1,
    wersja w slowniku 2). Jeden plik dla kluczy o zgodnych wartosciach miedzy slownikami i dtugi dla kluczy niezgodnych.
    Pliki majaca nazwe {file_prefix}_commoncalue_keys.txt, {file_prefix)_differentvalue_keys.txt
    :return: int, liczba kluczy o zgodnych wartosciach miwedzy slownikami, ponadtwo tworzone sa 2 pliki
    """
    #zrobimy ta petla
    common_keys = 0
    with open(f"{file_prefix}_commoncalue_keys.txt", 'w') as f1, open(f"{file_prefix}_differentvalue_keys.txt", 'w') as f2:
        for allel in slownik1:
            if allel in slownik2.keys():
                if slownik1[allel] == slownik2[allel]:
                    common_keys += 1
                    f1.write(f"{allel}\t{slownik1[allel]}\t{slownik2[allel]}\n")
                else:
                    f2.write(f"{allel}\t{slownik1[allel]}\t{slownik2[allel]}\n")
            else:
                #print(f"Allel {allel} nie wystepuje w drugim slowniku !\n")
                pass
    # na koniec wyrzumy jeszcze do plikow unikalne klucze dla slownikow 1 i 2
    slownik1_keys = set(slownik1.keys())
    slownik2_keys = set(slownik2.keys())
    with open(f'{file_prefix}_slownik1_unique_keys.txt', 'w') as f1, \
            open(f'{file_prefix}_slownik2_unique_keys.txt', 'w') as f2:
        f1.write("\n".join(slownik1_keys - slownik2_keys))
        f2.write("\n".join(slownik2_keys - slownik1_keys))

    return common_keys

def sort_profile(profile_file, sorted_keys_file):
    """
    Funkcja do sortowania/uzupelaniania profili o brakujace wpisy, zgodnie z zewnetrzna lista
    :param profile_file: str, sciezka do pliku z profile. Plik ma dwa wiersze, perwszy to lista loci (uwaga niezaleznie
    od nazwy pierwsza kolumna ZAWSZE traktowana jest jako nazwa dla Sequence type), drugi wiersz to int-y z allelami
    :param sorted_keys_file: sciezka do pliku z jednym wierszem gdzie "\t" oddzielono allele. Uwaga NIE podawac kolumny
    identyfikujacej ST
    :return: plik o nazwie basename(profile_file)_sorted_allells.txt
    """

    slownik_alleli = {}
    out_name = f"{os.path.basename(profile_file).split('.')[0]}_sorted_allells.txt"

    sorted_keys = open(sorted_keys_file).readlines()
    sorted_keys = sorted_keys[0].rsplit()

    with open(profile_file) as f:
        for line in f:
            if re.search('ST', line):
                klucze = line.rsplit()
                # czasami jako output jest podawany nazwa.fasta (od nazwy pliku do szukania alleli)
                klucze = [klucz.split('.')[0] for klucz in klucze]
            else:
                elementy = line.rsplit()
                for indeks, wartosc in enumerate(elementy):
                    slownik_alleli[klucze[indeks]] = wartosc


    with open(out_name, 'w') as f:
        out_msg = '800000' #dummy identyfikator ?
        for klucz in sorted_keys:
            try:
                slownik_alleli[klucz]
            except KeyError:
                slownik_alleli[klucz] = '-1'
            out_msg = out_msg + '\t' + slownik_alleli[klucz]


        f.write("ST" + "\t" + "\t".join(sorted_keys) + "\n")
        f.write(out_msg + "\n")

    return True

if __name__ == '__main__':
    pass
