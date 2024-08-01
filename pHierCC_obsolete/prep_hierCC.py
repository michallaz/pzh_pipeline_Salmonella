from all_functions_salmonella import *
moj_profil = parse_MLST_blastn('log.log')
slownik_profili, lista_kluczy = create_profile('/cgMLST2_entero/profiles.list')
ST, slownik_tmp_identity = getST(moj_profil, slownik_profili, 'profiles.list', lista_kluczy)
print(ST)
