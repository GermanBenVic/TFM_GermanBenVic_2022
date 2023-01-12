
import collections
import os
import re
import time
import pandas as pd
from common import *

class duoDragen:

	wd=""
	duos=""
	_filter=""
	
	def __init__(self,wd,duos,_filter): #Define __init___ que es el primer método que utiliza python para crear una instancia de la clase, luego wd, la ruta a los vcf y la ruta a los latest que son el resto de argumentos para definir la instancia self

		self.wd=wd # los argumentos los define como . de self para poder referirse a ellos
		self.duos=duos
		self._filter=_filter
	

	def run(this,caso,parient,sexo):

		if this.Filtrado(caso,parient,sexo):             
			this.Columnas(caso)
			return this.Switch(caso)


	def Switch(this,caso):
	 	#Reordenamiento de columnas

		output_dir=get_output_dir(this.wd,this._filter)

		try:
	 		f_final_file= (output_dir+"/" + str(caso) + "_final.tsv")
	 		f_final = pd.read_csv(f_final_file, sep="\t", header=0)
		except (OSError, IOError) as e:
	 		raise ErrorFinalFile(f_final_file)



		f_finalisimo = (output_dir+"/" + str(caso) + "_finalisimo.tsv")
	
		columns_titles=['annonimous_GENE', 'GEN_OMIM', 'PHENO_OMIM', 'PHENO_DESC', 'INHERITANCE', 'clinvar_id', 'clinvar_CLNSIG', 'clinvar_CLNSIGCONF', 'clinvar_CLNDN', 'annonimous_PREDICTORS', 'CADD_phred', 'occurrence_parents', 'occurrence', 'Locus', 'REF', 'genotype', 'Parient Genotype', 'Origin', 'zigosity', 'zigosity Parient', 'gnomad_AF', 'ALL_sites_2015_08', 'esp6500siv2_all', 'VAF', 'VAF Parient', 'Coverage', 'Parient Depth', 'AlleleRatio', 'AlleleCoverage', 'Parient AlleleCoverage', 'annonimous_ANNOTATION', 'annonimous_Func_refGene', 'annonimous_gNomen', 'annonimous_tCPUNTO', 'reliabilities', 'annonimous_PPUNTO', 'annonimous_EXONS', 'annonimous_INTRONS', 'annonimous_TOTAL', 'distNearestSS', 'nearestSSType', 'dbscSNV_ADA_SCORE', 'dbscSNV_RF_SCORE', 'catalogue_XREF', 'catalogue_VP', 'catalogue_speciality', 'ExClinicoNIM_panel', 'DI', 'SFARI', 'Epilepsia', 'RetNet', 'incidental_findings', 'pseudogenes', 'imprinting', 'PAR_loci', 'dbSNP_id', 'cosmicIds', 'cosmicTissues', 'cosmicSampleCounts', 'ParentsNIMIDs', 'IlluminaProbandsNIMIDs', 'SIFT_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 'LRT_pred', 'MutationTaster_pred', 'MutationAssessor_pred', 'FATHMM_pred', 'MetaSVM_pred', 'gnomad_AF_afr', 'gnomad_AF_amr', 'gnomad_AF_asj', 'gnomad_AF_eas', 'gnomad_AF_sas', 'gnomad_AF_nfe', 'gnomad_AF_fin', 'gnomad_AF_oth', 'gnomad_HetCount_all', 'gnomad_HomCount_all', 'annonimous_HemCount_all', 'gnomad_Filter', 'oe_mis', 'oe_syn', 'pLI', 'oe_lof', 'oe_syn_lower', 'oe_syn_upper', 'oe_mis_lower', 'oe_mis_upper', 'oe_lof_lower', 'oe_lof_upper', 'syn_z', 'mis_z', 'lof_z', 'AFR_sites_2015_08', 'SAS_sites_2015_08', 'EAS_sites_2015_08', 'EUR_sites_2015_08', 'AMR_sites_2015_08', 'esp6500siv2_ea', 'esp6500siv2_aa', 'phyloP100way_vertebrate', 'phyloP20way_mammalian', 'phyloP46way_placental', 'phastCons100way_vertebrate', 'phastCons20way_mammalian', 'APPROVED_NAME', 'SSF_wt', 'SSF_mut', 'MES_wt', 'MES_mut', 'delta_MES', 'delta_SSF', 'SPiCEprobability', 'SPiCEinter_2thr']
	
		f_final = f_final.reindex(columns=columns_titles)
		f_final = f_final.reset_index(drop=True)
		f_final.to_csv(f_finalisimo, sep='\t') 


	def Columnas(this,caso):
		# Abrimos archivos

		path_latest=get_path_latest(caso,this.wd,this._filter)
		my_path_tmp=get_path_tmp(this.wd,this._filter)
		output_dir=get_output_dir(this.wd,this._filter)

		try:
			f = open(path_latest)

		except (OSError, IOError) as e:
			raise ErrorInputFileLatest(path_latest)

		try:
			tmp_file=str(my_path_tmp) + "/duo_intermedio_" + str(caso) + ".tsv"
			f_intermedio = open(tmp_file)
		except (OSError, IOError) as e:
			raise ErrorTempFile(tmp_file)

		f_final = open(output_dir+"/" + str(caso) + "_final.tsv", "w")

		locusIndex = LocusIndex(f)

		f.close()

		f = open(path_latest)

		# AddColummns de intermedio a final (latest + intermedio)


		cabecera = f.readline().split("\n")[0] #Lee la cabecera del latest y se queda con la primera linea
		cabecera_final = cabecera.split("\t") #Parte la cabecera de latest por los tabuladores
		cabecera2 = f_intermedio.readline().split("\t") #Parte la cabecera del intermedio por los tabuladores
		# print (cabecera_final)
		# print (cabecera2)
		# cabecera.remove('\n')
		cabecera2.pop(0) #De la cabecera2 quita la primera posicion
		# print (cabecera)
		# print (cabecera2)

		rows = f.readlines() #Lee el latest tsv linea por linea y lo guarda en la variable rows
		rows2 = f_intermedio.readlines() #Lee el intermedio linea por linea y lo guarda en la variable rows
		rows3 = []

		f.close()
		f_intermedio.close()

		formato = ""
		formato2 = ""

		for i in range(len(cabecera_final)):
			formato2 = formato2 + cabecera_final[i] + "\t" #En formato2 guarda la cabecera final separada por tabuladores
		for j in range(len(cabecera2)):
			if j == len(cabecera2) - 1:
				formato2 = formato2 + cabecera2[j]
			else:
				formato2 = formato2 + cabecera2[j] + "\t"

		rows3.append(formato2) #Añade todas las columnas

		# Añadimos IGV  a las lineas en las que no se han encontrado datos en familiares

		for row in rows:
			prueba = row.split("\t\n") #Parte las lineas por tabulador y salto de línea y lo guarda en la variable prueba
			intermedio = prueba[0].split('\n') #Parte el primer argumento de prueba por salto de linea
			fields = intermedio[0].split('\t') #Parte el primer argumento de intermedio por los tabuladores
			Boo2 = True
			formato = ""
			for row2 in rows2:
				campos = row2.split("\t")
				if campos[0] == fields[locusIndex]:
					campos.pop(0)
					for i in range(len(fields)):
						formato = formato + fields[i] + "\t"
					for j in range(len(campos)):
						if j == len(campos) - 1:
							formato = formato + campos[j]
						else:
							formato = formato + campos[j] + "\t"
					rows3.append(formato)
					Boo2 = False
					break
				else:
					continue
			while Boo2 == True:
				for i in range(len(fields)):
					formato = formato + fields[i] + "\t"
				formato = formato + "IGV\tIGV\tIGV\tIGV\tIGV\n"
				rows3.append(formato)
				Boo2 = False

		for line in rows3:
			f_final.write(line)

		f_final.close()


	def Filtrado(this,caso,parient,sexo):
		# Abrimos archivos

		# Busca el latest en samples

		path_latest=get_path_latest(caso,this.wd,this._filter)
		duo_file=get_path_duo_file(caso,this.duos)
		my_path_tmp=get_path_tmp(this.wd,this._filter)

		try:
			f = open(path_latest) #Guardamos en f el archivo de la ruta latest
		except (OSError, IOError) as e:
			raise ErrorInputFileLatest(path_latest)

		# Busca el duo en análisis

		try:
			f_duo = open(duo_file) #Guardamos en f duo el archivo vcf del duo

		except (OSError, IOError) as e:
			raise ErrorInputFileDuo(duo_file)

		f_intermedio = open(str(my_path_tmp) + "/duo_intermedio_" + str(caso) + ".tsv", "w+") #Guardamos en f_intermedio el archivo intermedio

		locusIndex = LocusIndex(f) #Con la función LocusIndex genero el index del Locus

		locus = []
		duo = []
		lines = f.readlines() #Guardo las lineas del archivo latest en la variable lines

		f.close()

		for line in lines:
			fields = line.split('\t') #Guardo en la variable fields los campos de las lineas del archivo latest
			locus.append(fields[locusIndex]) #Se añaden los campos del locus index a la lista vacia locus

		i = 0

		find_header=False;
		header_proband=None

		while(not find_header): #Indica como header la primera linea del vcf que comienza con CHROM
			header_proband = f_duo.readline()
			if header_proband.startswith("#CHROM"):
				find_header=True


		try:
			probandIndex,parientIndex=getIndexFamily_2(caso,parient,header_proband) #Guarda el indice de cada miembro de la familia
		except NotFoundFamily as e:
			this.__log.error(e)
			return False

		lines_duo = f_duo.readlines() #Guarda en la variable lines_duo las lineas del archivo duo

		f_duo.close()

		for line in lines_duo:
			fields_duo = line.split("\t") #Parte las lineas del duo por los tabuladores generando los campos
			data_proband = fields_duo[probandIndex].split('\n')[0].split(":") #Guarda en la variable data_proband los campos del duo del probando partidos por salto de línea y :
			data_parient = fields_duo[parientIndex].split('\n')[0].split(":") #Lo mismo con el pariente
			
			#Prepara los datos para genotipar a los miembros de la familia
			
			if str(data_parient[0]) == '1' or str(data_parient[0]) == '0' or str(data_parient[0]) == '2':  #Si el genotipo esta representado solo por un numero
				data_parient[0] = str(data_parient[0]) + '/' + str(data_parient[0])

			if str(data_proband[0]) == '1' or str(data_proband[0]) == '0' or str(data_proband[0]) == '2':
				data_proband[0] = str(data_proband[0]) + '/' + str(data_proband[0])


			if '/' in data_proband[0]: #Si el genotipo está separado por un / dividelo en dos variables, una que guarda el primer genotipo y otra el segundo
				genotype_proband = data_proband[0].split("/")
				genotype_proband_0 = genotype_proband[0]
				genotype_proband_1 = genotype_proband[1]
			elif '|' in data_proband[0]:   #Si el genotipo está separado por un | dividelo en dos variables, una que guarda el primer genotipo y otra el segundo
				genotype_proband = data_proband[0].split("|")
				genotype_proband_0 = genotype_proband[0]
				genotype_proband_1 = genotype_proband[1]

			if '/' in data_parient[0]:
				genotype_parient = data_parient[0].split("/")
				genotype_parient_0 = genotype_parient[0]
				genotype_parient_1 = genotype_parient[1]
			elif '|' in data_parient[0]:
				genotype_parient = data_parient[0].split("|")
				genotype_parient_0 = genotype_parient[0]
				genotype_parient_1 = genotype_parient[1]

			#Guarda la cobertura de los miembros de la familia si la tienen

			if len(data_parient) < 3:  
				allele_coverage_parient = '-,-'
			elif len(data_parient[1].split(',')) == 1 and genotype_parient != './.':  
				allele_coverage_parient = data_parient[1] + ',0'
			else:
				allele_coverage_parient = data_parient[1] 

			ref_split = fields_duo[4].split(",")
			lista = range(len(ref_split))

			# HACERLO INDEXANDO VALORES Y QUE BUSQUE EN REF CON EL MISMO INDEX

			for field in locus:
				fields_locus = field.split(":")
				if fields_duo[0] == fields_locus[0] and fields_duo[1] == fields_locus[1]:   #Si el cromosoma y la posición coinciden entre el archivo multivcf y el latest del probando
					for h in lista:
						#print field #TODAS ESTAS RAMAS DE IF AÑADEN LA INFO DE COBERTURA Y GENOTIPO
						
						if genotype_parient[0] == "0": 
							if len(data_parient) < 3:  
								formato = fields_duo[0] + ":" + fields_duo[1] + "\t" + '-' + "\t" + \
										  fields_duo[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
										  fields_duo[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
										  fields_duo[3] + "/"
								break

							else:
								if data_parient[3] == '.':  
									formato = fields_duo[0] + ":" + fields_duo[1] + "\t" + \
											  data_parient[3] + \
											  "\t" + fields_duo[3] + '=' + \
											  allele_coverage_parient.split(',')[
												  0] + ', ' + fields_duo[4] + '=' + allele_coverage_parient.split(',')[
												  1] + '\t' + fields_duo[3] + "/"
								else:   
									formato = fields_duo[0] + ":" + fields_duo[1] + "\t" + \
											  str(int(allele_coverage_parient.split(',')[0]) + int(
												  allele_coverage_parient.split(',')[1])) + \
											  "\t" + fields_duo[3] + '=' + \
											  allele_coverage_parient.split(',')[
												  0] + ', ' + fields_duo[4] + '=' + allele_coverage_parient.split(',')[
												  1] + '\t' + fields_duo[3] + "/"
								break

						elif genotype_parient[0] == str(h + 1):
							if len(data_parient) < 3:
								formato = fields_duo[0] + ":" + fields_duo[1] + "\t" + '-' + "\t" + \
										  fields_duo[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
										  fields_duo[4] + '=' + allele_coverage_parient.split(',')[1] + '\t' + \
										  ref_split[h] + "/"
								break

							else:
								if data_parient[3] == '.':
									formato = fields_duo[0] + ":" + fields_duo[1] + "\t" + \
											  data_parient[3] + \
											  "\t" + fields_duo[3] + '=' + \
											  allele_coverage_parient.split(',')[
												  0] + ', ' + fields_duo[4] + '=' + allele_coverage_parient.split(',')[
												  1] + '\t' + ref_split[h] + "/"
								else:
										formato = fields_duo[0] + ":" + fields_duo[1] + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + \
												  "\t" + fields_duo[3] + '=' + \
												  allele_coverage_parient.split(',')[
													  0] + ', ' + fields_duo[4] + '=' + allele_coverage_parient.split(',')[
													  1] + '\t' + ref_split[h] + "/"
								break

						elif genotype_parient[0] == ".":
							if len(data_parient) < 3:
								formato = fields_duo[0] + ":" + fields_duo[1] + "\t" + '-' + "\t" + \
										  fields_duo[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
										  fields_duo[4] + '=' + allele_coverage_parient.split(',')[1] + '\t' + \
										  genotype_parient[0] + "/"
								break
							else:
								if data_parient[3] == '.':
									formato = fields_duo[0] + ":" + fields_duo[1] + "\t" + \
											  data_parient[3] + \
											  "\t" + fields_duo[3] + '=' + \
											  allele_coverage_parient.split(',')[
												  0] + ', ' + fields_duo[4] + '=' + allele_coverage_parient.split(',')[
												  1] + '\t' + genotype_parient[0] + "/"
								else:
										formato = fields_duo[0] + ":" + fields_duo[1] + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + \
												  "\t" + fields_duo[3] + '=' + \
												  allele_coverage_parient.split(',')[
													  0] + ', ' + fields_duo[4] + '=' + allele_coverage_parient.split(',')[
													  1] + '\t' + genotype_parient[0] + "/"
								break

					for h in lista:

						if genotype_parient[1] == "0":
							formato = formato + fields_duo[3] + "\t"
							break
						elif genotype_parient[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_parient[1] == ".":
							formato = formato + genotype_parient[1] + "\t"
							break

					#Añadimos la zigosidad

					for h in lista:
						if genotype_parient[0] == "0" and genotype_parient[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_parient[0] == "0" and genotype_parient[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t"
							break
						elif genotype_parient[0] == str(h + 1) and genotype_parient[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_parient[0] == "." or genotype_parient[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
							
					#Añadimos VAF del pariente
					
					for h in lista:
						if allele_coverage_parient.split(',') != "0" and '.' in allele_coverage_parient:
							formato = formato + ".,." 
							break
						elif allele_coverage_parient.split(',') != "0" and '-' in allele_coverage_parient:
							formato = formato + ".,." 
							break
						elif allele_coverage_parient.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_parient.split(',')[1])/ (float(allele_coverage_parient.split(',')[0]) + float(allele_coverage_parient.split(',')[1]))) 
							break
						elif allele_coverage_parient.split(',')[0] == "0" and allele_coverage_parient.split(',')[1] == "0":
							formato = formato + "0.0" 
							break
						else:
							formato = formato + "1.0" 
							break

					#Recojo la presencia de la variante en los pacientes a estudio

					if data_proband[0] == data_parient[0]:
						formato = formato + "\t" + "Both"
					else:
						formato = formato + "\t" + "Check"

					duo.append(formato)

		header_latest = "Locus" + "\t" + "Parient Depth" + "\t" + "Parient AlleleCoverage" + "\t" + "Parient Genotype" + "\t" + "zigosity Parient" + "\t" + "VAF Parient" + "\t" + "Origin" + "\n"


		f_intermedio.write(header_latest)

		# Meto los datos de genotipos, profundidades y origen de las variantes en un archivo intermedio

		for line in duo:
			# ~ print line
			f_intermedio.write("%s\n" % line)
		return True
		
		
class duoDragen_2P:

	wd=""
	duos=""
	_filter=""
	
	def __init__(self,wd,duos,_filter): 

		self.wd=wd 
		self.duos=duos
		self._filter=_filter


	

	def run(this,caso,caso2,sexo1,sexo2): #Se introducen los NR y datos de sexo y se corren todas las funciones del script por orden

		if this.Filtrado(caso,caso2,sexo1,sexo2):
			this.Columnas(caso,caso2)
			return this.Switch(caso,caso2)



	def Switch(this,caso,caso2):
	 	#Reordenamiento de columnas

		output_dir=get_output_dir(this.wd,this._filter)

		try:
	 		f_final_file= (output_dir+"/" + str(caso) + '+' + str(caso2) + "_final.tsv")
	 		f_final = pd.read_csv(f_final_file, sep="\t", header=0)
		except (OSError, IOError) as e:
	 		raise ErrorFinalFile(f_final_file)



		f_finalisimo = (output_dir+"/" + str(caso) + '+' + str(caso2) + "_finalisimo.tsv")

		f_final.drop(['genotype', 'zigosity', 'VAF', 'Coverage', 'AlleleRatio'], axis=1, inplace=True)
	
		columns_titles=['annonimous_GENE', 'GEN_OMIM', 'PHENO_OMIM', 'PHENO_DESC', 'INHERITANCE', 'clinvar_id', 'clinvar_CLNSIG', 'clinvar_CLNSIGCONF', 'clinvar_CLNDN', 'annonimous_PREDICTORS', 'CADD_phred', 'occurrence_parents', 'occurrence', 'Locus', 'REF', str(caso) + ' Genotype', str(caso2) + ' Genotype', 'Origin', 'zigosity ' + str(caso), 'zigosity ' + str(caso2), 'gnomad_AF', 'ALL_sites_2015_08', 'esp6500siv2_all', 'VAF ' +  str(caso), 'VAF ' +  str(caso2), str(caso) + ' Coverage', str(caso2) + ' Coverage', str(caso) + ' AlleleCoverage', str(caso2) + ' AlleleCoverage', 'annonimous_ANNOTATION', 'annonimous_Func_refGene', 'annonimous_gNomen', 'annonimous_tCPUNTO', 'reliabilities', 'annonimous_PPUNTO', 'annonimous_EXONS', 'annonimous_INTRONS', 'annonimous_TOTAL', 'distNearestSS', 'nearestSSType', 'dbscSNV_ADA_SCORE', 'dbscSNV_RF_SCORE', 'catalogue_XREF', 'catalogue_VP', 'catalogue_speciality', 'ExClinicoNIM_panel', 'DI', 'SFARI', 'Epilepsia', 'RetNet', 'incidental_findings', 'pseudogenes', 'imprinting', 'PAR_loci', 'dbSNP_id', 'cosmicIds', 'cosmicTissues', 'cosmicSampleCounts', 'ParentsNIMIDs', 'IlluminaProbandsNIMIDs', 'SIFT_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 'LRT_pred', 'MutationTaster_pred', 'MutationAssessor_pred', 'FATHMM_pred', 'MetaSVM_pred', 'gnomad_AF_afr', 'gnomad_AF_amr', 'gnomad_AF_asj', 'gnomad_AF_eas', 'gnomad_AF_sas', 'gnomad_AF_nfe', 'gnomad_AF_fin', 'gnomad_AF_oth', 'gnomad_HetCount_all', 'gnomad_HomCount_all', 'annonimous_HemCount_all', 'gnomad_Filter', 'oe_mis', 'oe_syn', 'pLI', 'oe_lof', 'oe_syn_lower', 'oe_syn_upper', 'oe_mis_lower', 'oe_mis_upper', 'oe_lof_lower', 'oe_lof_upper', 'syn_z', 'mis_z', 'lof_z', 'AFR_sites_2015_08', 'SAS_sites_2015_08', 'EAS_sites_2015_08', 'EUR_sites_2015_08', 'AMR_sites_2015_08', 'esp6500siv2_ea', 'esp6500siv2_aa', 'phyloP100way_vertebrate', 'phyloP20way_mammalian', 'phyloP46way_placental', 'phastCons100way_vertebrate', 'phastCons20way_mammalian', 'APPROVED_NAME', 'SSF_wt', 'SSF_mut', 'MES_wt', 'MES_mut', 'delta_MES', 'delta_SSF', 'SPiCEprobability', 'SPiCEinter_2thr']
	
		f_final = f_final.reindex(columns=columns_titles)
		f_final = f_final.reset_index(drop=True)
		f_final.to_csv(f_finalisimo, sep='\t')


	def Columnas(this,caso,caso2):
		
		# Abrimos archivos

		path_latest=get_path_latest_2P(caso,caso2,this.wd,this._filter)
		output_dir=get_output_dir(this.wd,this._filter)
		my_path_tmp=get_path_tmp(this.wd,this._filter)

		try:
			f = open(path_latest)


		except (OSError, IOError) as e:
			raise ErrorInputFileLatest(path_latest)

		try:
			tmp_file=str(my_path_tmp) + "/duo_intermedio_" + str(caso) + '+' + str(caso2) + ".tsv"
			f_intermedio = open(tmp_file)
		except (OSError, IOError) as e:
			raise ErrorTempFile(tmp_file)

		f_final = open(output_dir+"/" + str(caso) + '+' + str(caso2) + "_final.tsv", "w")

		locusIndex = LocusIndex(f)

		f.close()

		f = open(path_latest)

		# AddColummns de intermedio a final (latest + intermedio)

		cabecera = f.readline().split("\n")[0] #Lee la cabecera del latest y se queda con la primera linea
		cabecera_final = cabecera.split("\t") #Parte la cabecera de latest por los tabuladores
		cabecera2 = f_intermedio.readline().split("\t") #Parte la cabecera del intermedio por los tabuladores
		# print (cabecera_final)
		# print (cabecera2)
		# cabecera.remove('\n')
		cabecera2.pop(0) #De la cabecera2 quita la primera posicion
		# print (cabecera)
		# print (cabecera2)

		rows = f.readlines() #Lee el latest tsv linea por linea
		rows2 = f_intermedio.readlines() #Lee el intermedio linea por linea
		rows3 = []

		f.close()
		f_intermedio.close()

		formato = ""
		formato2 = ""

		for i in range(len(cabecera_final)):
			formato2 = formato2 + cabecera_final[i] + "\t" #En formato2 guarda la cabecera final separada por tabuladores
		for j in range(len(cabecera2)):
			if j == len(cabecera2) - 1:
				formato2 = formato2 + cabecera2[j]
			else:
				formato2 = formato2 + cabecera2[j] + "\t"

		rows3.append(formato2) #Añade todas las columnas

		# Añadimos IGV  a las lineas en las que no se han encontrado datos en los familiares

		for row in rows:
			prueba = row.split("\t\n") #Parte las lineas por tabulador y salto de línea y lo guarda en la variable prueba
			intermedio = prueba[0].split('\n') #Parte el primer argumento de prueba por salto de linea
			fields = intermedio[0].split('\t') #Parte el primer argumento de intermedio por los tabuladores
			Boo2 = True
			formato = ""
			for row2 in rows2:
				campos = row2.split("\t")
				if campos[0] == fields[locusIndex]:
					campos.pop(0)
					for i in range(len(fields)):
						formato = formato + fields[i] + "\t"
					for j in range(len(campos)):
						if j == len(campos) - 1:
							formato = formato + campos[j]
						else:
							formato = formato + campos[j] + "\t"
					rows3.append(formato)
					Boo2 = False
					break
				else:
					continue
			while Boo2 == True:
				for i in range(len(fields)):
					formato = formato + fields[i] + "\t"
				formato = formato + "IGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\n"
				rows3.append(formato)
				Boo2 = False
		
		for line in rows3:
			f_final.write(line)

		f_final.close()


	def Filtrado(this,caso,caso2,sexo1,sexo2):
		# Abrimos archivos
        
		#Generamos el archivo latest común
		
		path_latest_1=get_path_latest(caso,this.wd,this._filter)
		path_latest_2=get_path_latest(caso2,this.wd,this._filter)
		
		latest_1 = pd.read_csv(path_latest_1, sep="\t", header=0)
		latest_2 = pd.read_csv(path_latest_2, sep="\t", header=0)
		
		joint_latest = pd.concat([latest_1,latest_2], ignore_index=True).drop_duplicates(['Locus'])
		
		joint_latest[['chr','position']]=joint_latest.Locus.str.split(':',expand=True)
		joint_latest['chr']=joint_latest.chr.str.replace('chr','')
		joint_latest['chr']=joint_latest.chr.str.replace('X','30')
		joint_latest['chr']=joint_latest.chr.str.replace('Y','31')
		joint_latest['chr']=joint_latest.chr.str.replace('M','32')
		joint_latest['chr']=joint_latest['chr'].astype(int)
		joint_latest['position']=joint_latest['position'].astype(int)

		joint_latest_sort = joint_latest.sort_values(by=['chr','position'])
		joint_latest_sort = joint_latest_sort.drop(['chr','position'], axis=1)
		joint_latest_sort = joint_latest_sort.reset_index(drop=True)
		
		output_dir=get_output_dir(this.wd,this._filter)
		
		file_joint_latest = (output_dir+"/" + str(caso) + '+' + str(caso2) + "_latest.tsv")
		
		joint_latest_sort.to_csv(file_joint_latest, sep='\t')   

		# Busca el latest en samples

		path_latest=get_path_latest_2P(caso,caso2,this.wd,this._filter)
		duo_file=get_path_duo_file_2P(caso,caso2,this.duos)
		my_path_tmp=get_path_tmp(this.wd,this._filter)

		try:
			f = open(path_latest) #Guardamos en f el archivo de la ruta latest
		except (OSError, IOError) as e:
			raise ErrorInputFileLatest(path_latest)

		# Busca el duo en análisis

		try:
			f_duo = open(duo_file) #Guardamos en f duo el archivo vcf duo

		except (OSError, IOError) as e:
			raise ErrorInputFileDuo(duo_file)

		f_intermedio = open(str(my_path_tmp) + "/duo_intermedio_" + str(caso) + "+" + str(caso2) + ".tsv", "w+") #Guardamos en f_intermedio el archivo intermedio

		locusIndex = LocusIndex(f) #Con la función LocusIndex genero el index del Locus

		locus = []
		duo = []
		lines = f.readlines() #Guardo las lineas del archivo latest en la variable lines

		f.close()

		for line in lines:
			fields = line.split('\t') #Guardo en la variable fields los campos de las lineas del archivo latest
			locus.append(fields[locusIndex]) #Se añaden los campos del locus index a la lista vacia locus

		i = 0

		find_header=False;
		header_proband=None

		while(not find_header): #Indica como header la primera linea del vcf que comienza con CHROM
			header_proband = f_duo.readline()
			if header_proband.startswith("#CHROM"):
				find_header=True

		try:
			probandIndex,proband2Index=getIndexFamily_2_2P(caso,caso2,header_proband) #Guarda el indice de cada miembro de la familia
		except NotFoundFamily as e:
			this.__log.error(e)
			return False

		lines_duo = f_duo.readlines() #Guarda en la variable lines_duo las lineas del archivo duo

		f_duo.close()

		for line in lines_duo:
			fields_duo = line.split("\t") #Parte las lineas del duo por los tabuladores generando los campos
			data_proband = fields_duo[probandIndex].split('\n')[0].split(":") #Guarda en la variable data_proband los campos del duo del probando partidos por salto de línea y :
			data_proband2 = fields_duo[proband2Index].split('\n')[0].split(":") #Lo mismo con el pariente
			
			#Prepara los datos para genotipar a los miembros de la familia
			
			if str(data_proband2[0]) == '1' or str(data_proband2[0]) == '0' or str(data_proband2[0]) == '2':  #Si el genotipo esta representado solo por un numero
				data_proband2[0] = str(data_proband2[0]) + '/' + str(data_proband2[0])

			if str(data_proband[0]) == '1' or str(data_proband[0]) == '0' or str(data_proband[0]) == '2':
				data_proband[0] = str(data_proband[0]) + '/' + str(data_proband[0])


			if '/' in data_proband[0]: #Si el genotipo está separado por un / dividelo en dos variables, una que guarda el primer genotipo y otra el segundo
				genotype_proband = data_proband[0].split("/")
				genotype_proband_0 = genotype_proband[0]
				genotype_proband_1 = genotype_proband[1]
			elif '|' in data_proband[0]:   #Si el genotipo está separado por un | dividelo en dos variables, una que guarda el primer genotipo y otra el segundo
				genotype_proband = data_proband[0].split("|")
				genotype_proband_0 = genotype_proband[0]
				genotype_proband_1 = genotype_proband[1]

			if '/' in data_proband2[0]:
				genotype_proband2 = data_proband2[0].split("/")
				genotype_proband2_0 = genotype_proband2[0]
				genotype_proband2_1 = genotype_proband2[1]
			elif '|' in data_proband2[0]:
				genotype_proband2 = data_proband2[0].split("|")
				genotype_proband2_0 = genotype_proband2[0]
				genotype_proband2_1 = genotype_proband2[1]

			#Guarda la cobertura de los miembros de la familia si la tienen
			if len(data_proband2) < 3:
				allele_coverage_proband2 = '-,-'
			elif len(data_proband2[1].split(',')) == 1 and genotype_proband2 != './.':
				allele_coverage_proband2 = data_proband2[1] + ',0'
			else:
				allele_coverage_proband2 = data_proband2[1]

			if len(data_proband) < 3:
				allele_coverage_proband = '-,-'
			elif len(data_proband[1].split(',')) == 1 and genotype_proband != './.':
				allele_coverage_proband = data_proband[1] + ',0'
			else:
				allele_coverage_proband = data_proband[1]

			ref_split = fields_duo[4].split(",")
			lista = range(len(ref_split))

			# HACERLO INDEXANDO VALORES Y QUE BUSQUE EN REF CON EL MISMO INDEX

			for field in locus:
				fields_locus = field.split(":")
				if fields_duo[0] == fields_locus[0] and fields_duo[1] == fields_locus[1]:
					for h in lista:

						if len(data_proband) < 3:   #Si no tienen datos de profundidad en probando 1
							formato = fields_duo[0] + ":" + fields_duo[1] + "\t" + '-' + "\t" + \
									  str(int(allele_coverage_proband2.split(',')[0]) + int(allele_coverage_proband2.split(',')[1])) + "\t" + \
									  fields_duo[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
									  fields_duo[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
									  fields_duo[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
									  fields_duo[4] + '=' + allele_coverage_proband2.split(',')[1]
							break
						elif data_proband[3] == '.':   #Si la profundidad en probando 1 es un .
							formato = fields_duo[0] + ":" + fields_duo[1] + "\t" + data_proband[3] + "\t" + \
									  str(int(allele_coverage_proband2.split(',')[0]) + int(allele_coverage_proband2.split(',')[1])) + "\t" + \
									  fields_duo[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
									  fields_duo[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
									  fields_duo[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
									  fields_duo[4] + '=' + allele_coverage_proband2.split(',')[1] 
						else:   # Si tiene dato de profundidad en el probando 1
							if len(data_proband2) <3:  ##Si el probando 2 no tiene dato de profundidad
								formato = fields_duo[0] + ":" + fields_duo[1] + "\t" + str(int(allele_coverage_proband.split(',')[0]) + \
									  int(allele_coverage_proband.split(',')[1])) + "\t" + '-' + "\t" + \
									  fields_duo[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
									  fields_duo[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
									  fields_duo[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
									  fields_duo[4] + '=' + allele_coverage_proband2.split(',')[1]
								break
							elif data_proband2[3] == '.':
								formato = fields_duo[0] + ":" + fields_duo[1] + "\t" + str(int(allele_coverage_proband.split(',')[0]) + \
									  int(allele_coverage_proband.split(',')[1])) + "\t" + data_proband2[3] + "\t" + \
									  fields_duo[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
									  fields_duo[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
									  fields_duo[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
									  fields_duo[4] + '=' + allele_coverage_proband2.split(',')[1]
							else:   #Si ambos tienen datos reales
								formato = fields_duo[0] + ":" + fields_duo[1] + "\t" + \
										  str(int(allele_coverage_proband.split(',')[0]) + int(
											  allele_coverage_proband.split(',')[1])) + \
										  "\t" + \
										  str(int(allele_coverage_proband2.split(',')[0]) + int(
											  allele_coverage_proband2.split(',')[1])) + \
										  "\t" + fields_duo[3] + '=' + \
										  allele_coverage_proband.split(',')[0] + ', ' + fields_duo[4] + '=' + \
										  allele_coverage_proband.split(',')[1] + "\t" + fields_duo[3] + '=' + \
										  allele_coverage_proband2.split(',')[0] + ', ' + fields_duo[4] + '=' + \
										  allele_coverage_proband2.split(',')[1] 
							break
								


					for h in lista:
						if genotype_proband[0] == "0":
							formato = formato +  "\t" + fields_duo[3] + "/"
							break
						elif genotype_proband[0] == str(h + 1):
							formato = formato +  "\t" + ref_split[h] + "/"
							break
						elif genotype_proband[0] == ".":
							formato = formato +  "\t" + genotype_proband[0] + "/"
							break
					for h in lista:

						if genotype_proband[1] == "0":
							formato = formato + fields_duo[3] + "\t"
							break
						elif genotype_proband[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_proband[1] == ".":
							formato = formato + genotype_proband[1] + "\t"
							break
					for h in lista:
						if genotype_proband2[0] == "0":
							formato = formato + fields_duo[3] + "/"
							break
						elif genotype_proband2[0] == str(h + 1):
							formato = formato + ref_split[h] + "/"
							break
						elif genotype_proband2[0] == ".":
							formato = formato + genotype_proband2[0] + "/"
							break
					for h in lista:
						if genotype_proband2[1] == "0":
							formato = formato + fields_duo[3] + "\t"
							break
						elif genotype_proband2[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_proband2[1] == ".":
							formato = formato + genotype_proband2[1] + "\t"
							break
					
					#Añadimos la zigosidad
					
					for h in lista:
						if genotype_proband[0] == "0" and genotype_proband[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_proband[0] == "0" and genotype_proband[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t" 
							break
						elif genotype_proband[0] == str(h + 1) and genotype_proband[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_proband[0] == "." and genotype_proband[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
							
					for h in lista:
						if genotype_proband2[0] == "0" and genotype_proband2[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_proband2[0] == "0" and genotype_proband2[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t"
							break
						elif genotype_proband2[0] == str(h + 1) and genotype_proband2[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_proband2[0] == "." and genotype_proband2[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
					
					#Añadimos VAF	
					
					for h in lista:
						if allele_coverage_proband.split(',') != "0" and '.' in allele_coverage_proband:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_proband.split(',') != "0" and '-' in allele_coverage_proband:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_proband.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_proband.split(',')[1])/ (float(allele_coverage_proband.split(',')[0]) + float(allele_coverage_proband.split(',')[1]))) + "\t" 
							break
						elif allele_coverage_proband.split(',')[0] == "0" and allele_coverage_proband.split(',')[1] == "0":
							formato = formato + "0.0" + "\t"
							break
						else:
							formato = formato + "1.0" + "\t"
							break
					for h in lista:
						if allele_coverage_proband2.split(',') != "0" and '.' in allele_coverage_proband2:
							formato = formato + ".,." 
							break
						elif allele_coverage_proband2.split(',') != "0" and '-' in allele_coverage_proband2:
							formato = formato + ".,." 
							break
						elif allele_coverage_proband2.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_proband2.split(',')[1])/ (float(allele_coverage_proband2.split(',')[0]) + float(allele_coverage_proband2.split(',')[1]))) 
							break
						elif allele_coverage_proband2.split(',')[0] == "0" and allele_coverage_proband2.split(',')[1] == "0":
							formato = formato + "0.0" 
							break
						else:
							formato = formato + "1.0" 
							break
							
					#Recojo el origen del genotipo
					
					if data_proband[0] == data_proband2[0]:
						formato = formato + "\t" + "Common"
					else:
						formato = formato + "\t" + "Check"


					duo.append(formato)

		header_latest = "Locus" + "\t" +  str(caso) + " Coverage" + "\t" + str(caso2) + " Coverage" + "\t" + str(caso) + " AlleleCoverage" + "\t" + str(caso2) + " AlleleCoverage" + "\t" + str(caso) + " Genotype" + "\t" + str(caso2) + " Genotype" + "\t" + "zigosity " + str(caso) + "\t" + "zigosity " + str(caso2) + "\t" + "VAF " +  str(caso) + "\t" + "VAF " + str(caso2) + "\t" + "Origin" + "\n"

		f_intermedio.write(header_latest)

		# Meto los datos de genotipos, profundidades y origen de las variantes en un archivo intermedio

		for line in duo:
			# ~ print line
			f_intermedio.write("%s\n" % line)
		return True
    
class trioDragen:

	wd=""
	trios=""
	_filter=""
	
	def __init__(self,wd,trios,_filter): #Define __init___ que es el primer método que utiliza python para crear una instancia de trioDragen, luego wd, trios y trio_filter que son el resto de argumentos para definir la instancia self


		self.wd=wd # los argumentos los define como . de self para poder referirse a ellos
		self.trios=trios
		self._filter=_filter


	

	def run(this,proband,father,mother,sexo):


		if this.Filtrado(proband,father,mother,sexo):
			this.Columnas(proband)
			return this.Switch(proband)


	def Switch(this,caso):
	#Reordenamiento de columnas

		output_dir=get_output_dir(this.wd,this._filter)

		try:
	 		f_final_file= (output_dir+"/" + str(caso) + "_final.tsv")
	 		f_final = pd.read_csv(f_final_file, sep="\t", header=0)
		except (OSError, IOError) as e:
	 		raise ErrorFinalFile(f_final_file)



		f_finalisimo = (output_dir+"/" + str(caso) + "_finalisimo.tsv")
	
		columns_titles=['annonimous_GENE', 'GEN_OMIM', 'PHENO_OMIM', 'PHENO_DESC', 'INHERITANCE', 'clinvar_id', 'clinvar_CLNSIG', 'clinvar_CLNSIGCONF', 'clinvar_CLNDN', 'annonimous_PREDICTORS', 'CADD_phred', 'occurrence_parents', 'occurrence', 'Locus', 'REF', 'genotype', 'Mother Genotype', 'Father Genotype', 'Origin', 'zigosity', 'zigosity Mother', 'zigosity Father', 'gnomad_AF', 'ALL_sites_2015_08', 'esp6500siv2_all', 'VAF', 'VAF Mother', 'VAF Father', 'Coverage', 'Mother Depth', 'Father Depth', 'AlleleRatio', 'AlleleCoverage', 'Mother AlleleCoverage', 'Father AlleleCoverage', 'annonimous_ANNOTATION', 'annonimous_Func_refGene', 'annonimous_gNomen', 'annonimous_tCPUNTO', 'reliabilities', 'annonimous_PPUNTO', 'annonimous_EXONS', 'annonimous_INTRONS', 'annonimous_TOTAL', 'distNearestSS', 'nearestSSType', 'dbscSNV_ADA_SCORE', 'dbscSNV_RF_SCORE', 'catalogue_XREF', 'catalogue_VP', 'catalogue_speciality', 'ExClinicoNIM_panel', 'DI', 'SFARI', 'Epilepsia', 'RetNet', 'incidental_findings', 'pseudogenes', 'imprinting', 'PAR_loci', 'dbSNP_id', 'cosmicIds', 'cosmicTissues', 'cosmicSampleCounts', 'ParentsNIMIDs', 'IlluminaProbandsNIMIDs', 'SIFT_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 'LRT_pred', 'MutationTaster_pred', 'MutationAssessor_pred', 'FATHMM_pred', 'MetaSVM_pred', 'gnomad_AF_afr', 'gnomad_AF_amr', 'gnomad_AF_asj', 'gnomad_AF_eas', 'gnomad_AF_sas', 'gnomad_AF_nfe', 'gnomad_AF_fin', 'gnomad_AF_oth', 'gnomad_HetCount_all', 'gnomad_HomCount_all', 'annonimous_HemCount_all', 'gnomad_Filter', 'oe_mis', 'oe_syn', 'pLI', 'oe_lof', 'oe_syn_lower', 'oe_syn_upper', 'oe_mis_lower', 'oe_mis_upper', 'oe_lof_lower', 'oe_lof_upper', 'syn_z', 'mis_z', 'lof_z', 'AFR_sites_2015_08', 'SAS_sites_2015_08', 'EAS_sites_2015_08', 'EUR_sites_2015_08', 'AMR_sites_2015_08', 'esp6500siv2_ea', 'esp6500siv2_aa', 'phyloP100way_vertebrate', 'phyloP20way_mammalian', 'phyloP46way_placental', 'phastCons100way_vertebrate', 'phastCons20way_mammalian', 'APPROVED_NAME', 'SSF_wt', 'SSF_mut', 'MES_wt', 'MES_mut', 'delta_MES', 'delta_SSF', 'SPiCEprobability', 'SPiCEinter_2thr']
	
		f_final = f_final.reindex(columns=columns_titles)
		f_final = f_final.reset_index(drop=True)
		f_final.to_csv(f_finalisimo, sep='\t')  




	def Columnas(this,caso):
		# Abrimos archivos

		path_latest=get_path_latest(caso,this.wd,this._filter)
		my_path_tmp=get_path_tmp(this.wd,this._filter)
		output_dir=get_output_dir(this.wd,this._filter)

		try:
			f = open(path_latest)

		except (OSError, IOError) as e:
			raise ErrorInputFileLatest(path_latest)

		try:
			tmp_file=str(my_path_tmp) + "/trio_intermedio_" + str(caso) + ".tsv"
			f_intermedio = open(tmp_file)
		except (OSError, IOError) as e:
			raise ErrorTempFile(tmp_file)

		f_final = open(output_dir+"/" + str(caso) + "_final.tsv", "w")

		locusIndex = LocusIndex(f)

		f.close()

		f = open(path_latest)

		# AddColummns de intermedio a final (latest + intermedio)

		cabecera = f.readline().split("\n")[0] #Lee la cabecera y se queda con la primera linea
		cabecera_final = cabecera.split("\t") #Parte la cabecera de latest por los tabuladores
		cabecera2 = f_intermedio.readline().split("\t") #Parte la cabecera del intermedio por los tabuladores
		# print (cabecera_final)
		# print (cabecera2)
		# cabecera.remove('\n')
		cabecera2.pop(0) #De la cabecera2 quita la primera posicion
		# print (cabecera)
		# print (cabecera2)

		rows = f.readlines() #Lee el latest tsv linea por linea
		rows2 = f_intermedio.readlines() #Lee el intermedio linea por linea
		rows3 = []

		f.close()
		f_intermedio.close()

		formato = ""
		formato2 = ""

		for i in range(len(cabecera_final)):
			formato2 = formato2 + cabecera_final[i] + "\t" #En formato2 guarda la cabecera final separada por tabuladores
		for j in range(len(cabecera2)):
			if j == len(cabecera2) - 1:
				formato2 = formato2 + cabecera2[j]
			else:
				formato2 = formato2 + cabecera2[j] + "\t"

		rows3.append(formato2) #Fusiona las filas

		# Añadimos IGV  a las lineas en las que no se han encontrado datos en los padres en el archivo trio

		for row in rows:
			prueba = row.split("\t\n") #Parte las lineas por tabulador y salto de línea y lo guarda en la variable prueba
			intermedio = prueba[0].split('\n') #Parte el primer argumento de prueba por salto de linea
			fields = intermedio[0].split('\t') #Parte el primer argumento de intermedio por los tabuladores
			Boo2 = True
			formato = ""
			for row2 in rows2:
				campos = row2.split("\t")
				if campos[0] == fields[locusIndex]:
					campos.pop(0)
					for i in range(len(fields)):
						formato = formato + fields[i] + "\t"
					for j in range(len(campos)):
						if j == len(campos) - 1:
							formato = formato + campos[j]
						else:
							formato = formato + campos[j] + "\t"
					rows3.append(formato)
					Boo2 = False
					break
				else:
					continue
			while Boo2 == True:
				for i in range(len(fields)):
					formato = formato + fields[i] + "\t"
				formato = formato + "IGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\n"
				rows3.append(formato)
				Boo2 = False

		for line in rows3:
			f_final.write(line)

		f_final.close()


	def Filtrado(this,caso, padre, madre, sexo):
		# Abrimos archivos

		# my_path = os.getcwd()

		# Busca el latest en samples

		path_latest=get_path_latest(caso,this.wd,this._filter)
		trio_file=get_path_trio_file(caso,this.trios)
		my_path_tmp=get_path_tmp(this.wd,this._filter)

		try:
			f = open(path_latest) #Guardamos en f erl archivo de la ruta latest
		except (OSError, IOError) as e:
			raise ErrorInputFileLatest(path_latest)

		# Busca el trio en análisis

		try:
			f_trio = open(trio_file) #Guardamos en f trio el archivo trio

		except (OSError, IOError) as e:
			raise ErrorInputFileTrio(trio_file)

		f_intermedio = open(str(my_path_tmp) + "/trio_intermedio_" + str(caso) + ".tsv", "w+") #Guardamos en f_intermedio el archivo intermedio

		locusIndex = LocusIndex(f) #Con la función LocusIndex genero el index del Locus

		locus = []
		trio = []
		lines = f.readlines() #Guardo las lineas del archivo latest en la variable lines

		f.close()

		for line in lines:
			fields = line.split('\t') #Guardo en la variable fields los campos de las lineas del archivo latest
			locus.append(fields[locusIndex]) #Se añade el index del locus a los campos del archivo latest


		i = 0

		find_header=False;
		header_proband=None

		while(not find_header): #Indica como header la primera linea del vcf que comienza con CHROM
			header_proband = f_trio.readline()
			if header_proband.startswith("#CHROM"):
				find_header=True


		try:
			probandIndex,fatherIndex,motherIndex=getIndexFamily(caso,padre,madre,header_proband) #Guarda el indice de cada miembro de la familia
		except NotFoundFamily as e:
			this.__log.error(e)
			return False

		lines_trio = f_trio.readlines() #Guarda en la variable lines_trio las lineas del archivo vcf trio

		f_trio.close()

		for line in lines_trio:
			fields_trio = line.split("\t") #Parte las lineas del trio por los tabuladores generando los campos
			data_proband = fields_trio[probandIndex].split('\n')[0].split(":") #Guarda en la variable data_proband los campos del trio del probando partidos por salto de línea y :
			data_mother = fields_trio[motherIndex].split('\n')[0].split(":") #Lo mismo con la madre
			data_father = fields_trio[fatherIndex].split('\n')[0].split(":") #Lo mismo con el padre
			
			#Prepara los datos para genotipar a los miembros de la familia
			
			if str(data_father[0]) == '1' or str(data_father[0]) == '0' or str(data_father[0]) == '2':
				data_father[0] = str(data_father[0]) + '/' + str(data_father[0])

			if str(data_proband[0]) == '1' or str(data_proband[0]) == '0' or str(data_proband[0]) == '2':
				data_proband[0] = str(data_proband[0]) + '/' + str(data_proband[0])


			if '/' in data_proband[0]:
				genotype_proband = data_proband[0].split("/")
				genotype_proband_0 = genotype_proband[0]
				genotype_proband_1 = genotype_proband[1]
			elif '|' in data_proband[0]:
				genotype_proband = data_proband[0].split("|")
				genotype_proband_0 = genotype_proband[0]
				genotype_proband_1 = genotype_proband[1]

			if '/' in data_mother[0]:
				genotype_mother = data_mother[0].split("/")
				genotype_mother_0 = genotype_mother[0]
				genotype_mother_1 = genotype_mother[1]
			elif '|' in data_mother[0]:
				genotype_mother = data_mother[0].split("|")
				genotype_mother_0 = genotype_proband[0]
				genotype_mother_1 = genotype_proband[1]

			if '/' in data_father[0]:
				genotype_father = data_father[0].split("/")
				genotype_father_0 = genotype_father[0]
				genotype_father_1 = genotype_father[1]
			elif '|' in data_father[0]:
				genotype_father = data_father[0].split("|")
				genotype_father_0 = genotype_father[0]
				genotype_father_1 = genotype_father[1]

			#Guarda la cobertura de los miembros de la familia si la tienen
			if len(data_mother) < 3:
				allele_coverage_mother = '-,-'
			elif len(data_mother[1].split(',')) == 1 and genotype_mother != './.':
				allele_coverage_mother = data_mother[1] + ',0'
			else:
				allele_coverage_mother = data_mother[1]

			if len(data_father) < 3:
				allele_coverage_father = '-,-'
			elif len(data_father[1].split(',')) == 1 and genotype_father != './.':
				allele_coverage_father = data_father[1] + ',0'
			else:
				allele_coverage_father = data_father[1]

			ref_split = fields_trio[4].split(",")
			lista = range(len(ref_split))

			# HACERLO INDEXANDO VALORES Y QUE BUSQUE EN REF CON EL MISMO INDEX

			for field in locus:
				fields_locus = field.split(":")
				if fields_trio[0] == fields_locus[0] and fields_trio[1] == fields_locus[1]:
					for h in lista:
						#print field #TODAS ESTAS RAMAS DE IF AÑADEN LA INFO DE COBERTURA Y GENOTIPO DEL PADRE Y LA MADRE
						if genotype_mother[0] == "0":
							if len(data_mother) < 3:
								if len(data_father) < 3:
									formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + '-' + "\t" + '-' + "\t" + \
											  fields_trio[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
											  fields_trio[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
											  fields_trio[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
											  fields_trio[4] + '=' + allele_coverage_father.split(',')[1] + '\t' + \
											  fields_trio[3] + "/"
									break
								else:
									if data_father[3] == '.':
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + '-' + "\t" + \
												  data_father[3] + \
												  "\t" + fields_trio[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_trio[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_trio[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_trio[4] + '=' + allele_coverage_father.split(',')[1] + '\t' + \
												  fields_trio[3] + "/"
									else:
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + \
												  "\t" + fields_trio[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_trio[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_trio[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_trio[4] + '=' + allele_coverage_father.split(',')[1] + '\t' + \
												  fields_trio[3] + "/"
									break
							else:
								if len(data_father) < 3:
									if data_mother[3] == '.':
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + \
												  data_mother[3] + \
												  "\t" + '-' + "\t" + fields_trio[3] + '=' + \
												  allele_coverage_mother.split(',')[
													  0] + ', ' + fields_trio[4] + '=' + allele_coverage_mother.split(',')[
													  1] + "\t" + fields_trio[3] + '=' + allele_coverage_father.split(',')[
													  0] + ', ' + fields_trio[4] + '=' + allele_coverage_father.split(',')[
													  1] + '\t' + fields_trio[3] + "/"
									else:
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
													  allele_coverage_mother.split(',')[1])) + \
												  "\t" + '-' + "\t" + fields_trio[3] + '=' + \
												  allele_coverage_mother.split(',')[
													  0] + ', ' + fields_trio[4] + '=' + allele_coverage_mother.split(',')[
													  1] + "\t" + fields_trio[3] + '=' + allele_coverage_father.split(',')[
													  0] + ', ' + fields_trio[4] + '=' + allele_coverage_father.split(',')[
													  1] + '\t' + fields_trio[3] + "/"
									break
								else:
									if data_mother[3] == '.' and data_father[3] == '.':
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + \
												  data_mother[3] + \
												  "\t" + \
												  data_father[3] + \
												  "\t" + fields_trio[3] + '=' + \
												  allele_coverage_mother.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_mother.split(',')[1] + "\t" + fields_trio[3] + '=' + \
												  allele_coverage_father.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_father.split(',')[1] + '\t' + fields_trio[3] + "/"
									elif data_father[3] == '.':
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
													  allele_coverage_mother.split(',')[1])) + \
												  "\t" + \
												  data_father[3] + \
												  "\t" + fields_trio[3] + '=' + \
												  allele_coverage_mother.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_mother.split(',')[1] + "\t" + fields_trio[3] + '=' + \
												  allele_coverage_father.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_father.split(',')[1] + '\t' + fields_trio[3] + "/"
									elif data_mother[3] == '.':
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + \
												  data_mother[3] + \
												  "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + \
												  "\t" + fields_trio[3] + '=' + \
												  allele_coverage_mother.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_mother.split(',')[1] + "\t" + fields_trio[3] + '=' + \
												  allele_coverage_father.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_father.split(',')[1] + '\t' + fields_trio[3] + "/"
									else:
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
													  allele_coverage_mother.split(',')[1])) + \
												  "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + \
												  "\t" + fields_trio[3] + '=' + \
												  allele_coverage_mother.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_mother.split(',')[1] + "\t" + fields_trio[3] + '=' + \
												  allele_coverage_father.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_father.split(',')[1] + '\t' + fields_trio[3] + "/"
									break

						elif genotype_mother[0] == str(h + 1):
							if len(data_mother) < 3:
								if len(data_father) < 3:
									formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + '-' + "\t" + '-' + "\t" + \
											  fields_trio[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
											  fields_trio[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
											  fields_trio[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
											  fields_trio[4] + '=' + allele_coverage_father.split(',')[1] + '\t' + \
											  ref_split[h] + "/"
									break
								else:
									if data_father[3] == '.':
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + '-' + "\t" + \
												  data_father[3] + \
												  "\t" + fields_trio[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_trio[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_trio[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_trio[4] + '=' + allele_coverage_father.split(',')[1] + '\t' + \
												  ref_split[h] + "/"
									else:
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + \
												  "\t" + fields_trio[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_trio[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_trio[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_trio[4] + '=' + allele_coverage_father.split(',')[1] + '\t' + \
												  ref_split[h] + "/"
									break
							else:
								if len(data_father) < 3:
									if data_mother[3] == '.':
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + \
												  data_mother[3] + \
												  "\t" + '-' + "\t" + fields_trio[3] + '=' + \
												  allele_coverage_mother.split(',')[
													  0] + ', ' + fields_trio[4] + '=' + allele_coverage_mother.split(',')[
													  1] + "\t" + fields_trio[3] + '=' + allele_coverage_father.split(',')[
													  0] + ', ' + fields_trio[4] + '=' + allele_coverage_father.split(',')[
													  1] + '\t' + ref_split[h] + "/"
									else:
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
													  allele_coverage_mother.split(',')[1])) + \
												  "\t" + '-' + "\t" + fields_trio[3] + '=' + \
												  allele_coverage_mother.split(',')[
													  0] + ', ' + fields_trio[4] + '=' + allele_coverage_mother.split(',')[
													  1] + "\t" + fields_trio[3] + '=' + allele_coverage_father.split(',')[
													  0] + ', ' + fields_trio[4] + '=' + allele_coverage_father.split(',')[
													  1] + '\t' + ref_split[h] + "/"
									break
								else:
									if data_mother[3] == '.' and data_father[3] == '.':
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + \
												  data_mother[3] + \
												  "\t" + \
												  data_father[3] + \
												  "\t" + fields_trio[3] + '=' + \
												  allele_coverage_mother.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_mother.split(',')[1] + "\t" + fields_trio[3] + '=' + \
												  allele_coverage_father.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_father.split(',')[1] + '\t' + ref_split[h] + "/"
									elif data_mother[3] == '.':
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + \
												  data_mother[3] + \
												  "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + \
												  "\t" + fields_trio[3] + '=' + \
												  allele_coverage_mother.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_mother.split(',')[1] + "\t" + fields_trio[3] + '=' + \
												  allele_coverage_father.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_father.split(',')[1] + '\t' + ref_split[h] + "/"
									elif data_father[3] == '.':
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
													  allele_coverage_mother.split(',')[1])) + \
												  "\t" + \
												  data_father[3] + \
												  "\t" + fields_trio[3] + '=' + \
												  allele_coverage_mother.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_mother.split(',')[1] + "\t" + fields_trio[3] + '=' + \
												  allele_coverage_father.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_father.split(',')[1] + '\t' + ref_split[h] + "/"
									else:
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
													  allele_coverage_mother.split(',')[1])) + \
												  "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + \
												  "\t" + fields_trio[3] + '=' + \
												  allele_coverage_mother.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_mother.split(',')[1] + "\t" + fields_trio[3] + '=' + \
												  allele_coverage_father.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_father.split(',')[1] + '\t' + ref_split[h] + "/"
									break

						elif genotype_mother[0] == ".":
							if len(data_mother) < 3:
								if len(data_father) < 3:
									formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + '-' + "\t" + '-' + "\t" + \
											  fields_trio[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
											  fields_trio[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
											  fields_trio[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
											  fields_trio[4] + '=' + allele_coverage_father.split(',')[1] + '\t' + \
											  genotype_mother[0] + "/"
									break
								else:
									if data_father[3] == '.':
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + '-' + "\t" + \
												  data_father[3] + \
												  "\t" + fields_trio[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_trio[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_trio[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_trio[4] + '=' + allele_coverage_father.split(',')[1] + '\t' + \
												  genotype_mother[0] + "/"
									else:
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + \
												  "\t" + fields_trio[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_trio[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_trio[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_trio[4] + '=' + allele_coverage_father.split(',')[1] + '\t' + \
												  genotype_mother[0] + "/"
									break
							else:
								if len(data_father) < 3:
									if data_mother[3] == '.':
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + \
												  data_mother[3] + \
												  "\t" + '-' + "\t" + fields_trio[3] + '=' + \
												  allele_coverage_mother.split(',')[
													  0] + ', ' + fields_trio[4] + '=' + allele_coverage_mother.split(',')[
													  1] + "\t" + fields_trio[3] + '=' + allele_coverage_father.split(',')[
													  0] + ', ' + fields_trio[4] + '=' + allele_coverage_father.split(',')[
													  1] + '\t' + genotype_mother[0] + "/"
									else:
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
													  allele_coverage_mother.split(',')[1])) + \
												  "\t" + '-' + "\t" + fields_trio[3] + '=' + \
												  allele_coverage_mother.split(',')[
													  0] + ', ' + fields_trio[4] + '=' + allele_coverage_mother.split(',')[
													  1] + "\t" + fields_trio[3] + '=' + allele_coverage_father.split(',')[
													  0] + ', ' + fields_trio[4] + '=' + allele_coverage_father.split(',')[
													  1] + '\t' + genotype_mother[0] + "/"
									break
								else:
									if data_mother[3] == '.' and data_father[3] == '.':
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + \
												  data_mother[3] + \
												  "\t" + \
												  data_father[3] + \
												  "\t" + fields_trio[3] + '=' + \
												  allele_coverage_mother.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_mother.split(',')[1] + "\t" + fields_trio[3] + '=' + \
												  allele_coverage_father.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_father.split(',')[1] + '\t' + genotype_mother[0] + "/"
									elif data_mother[3] == '.':
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + \
												  data_mother[3] + \
												  "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + \
												  "\t" + fields_trio[3] + '=' + \
												  allele_coverage_mother.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_mother.split(',')[1] + "\t" + fields_trio[3] + '=' + \
												  allele_coverage_father.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_father.split(',')[1] + '\t' + genotype_mother[0] + "/"
									elif data_father[3] == '.':
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
													  allele_coverage_mother.split(',')[1])) + \
												  "\t" + \
												  data_father[3] + \
												  "\t" + fields_trio[3] + '=' + \
												  allele_coverage_mother.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_mother.split(',')[1] + "\t" + fields_trio[3] + '=' + \
												  allele_coverage_father.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_father.split(',')[1] + '\t' + genotype_mother[0] + "/"
									else:
										formato = fields_trio[0] + ":" + fields_trio[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
													  allele_coverage_mother.split(',')[1])) + \
												  "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + \
												  "\t" + fields_trio[3] + '=' + \
												  allele_coverage_mother.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_mother.split(',')[1] + "\t" + fields_trio[3] + '=' + \
												  allele_coverage_father.split(',')[0] + ', ' + fields_trio[4] + '=' + \
												  allele_coverage_father.split(',')[1] + '\t' + genotype_mother[0] + "/"
									break
					#Añadimos genotipo
					
					for h in lista:

						if genotype_mother[1] == "0":
							formato = formato + fields_trio[3] + "\t"
							break
						elif genotype_mother[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_mother[1] == ".":
							formato = formato + genotype_mother[1] + "\t"
							break
					for h in lista:
						if genotype_father[0] == "0":
							formato = formato + fields_trio[3] + "/"
							break
						elif genotype_father[0] == str(h + 1):
							formato = formato + ref_split[h] + "/"
							break
						elif genotype_father[0] == ".":
							formato = formato + genotype_father[0] + "/"
							break
					for h in lista:
						if genotype_father[1] == "0":
							formato = formato + fields_trio[3] + "\t"
							break
						elif genotype_father[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_father[1] == ".":
							formato = formato + genotype_father[1] + "\t"
							break
							
					#Añadimos la zigosidad
					
					for h in lista:
						if genotype_mother[0] == "0" and genotype_mother[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_mother[0] == "0" and genotype_mother[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t" 
							break
						elif genotype_mother[0] == str(h + 1) and genotype_mother[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_mother[0] == "." or genotype_mother[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
							
					for h in lista:
						if genotype_father[0] == "0" and genotype_father[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_father[0] == "0" and genotype_father[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t" 
							break
						elif genotype_father[0] == str(h + 1) and genotype_father[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_father[0] == "." or genotype_father[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
							
					#Añadimos el VAF
					
					for h in lista:
						if allele_coverage_mother.split(',') != "0" and '.' in allele_coverage_mother:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_mother.split(',') != "0" and '-' in allele_coverage_mother:
							formato = formato + ".,." 
							break
						elif allele_coverage_mother.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_mother.split(',')[1])/ (float(allele_coverage_mother.split(',')[0]) + float(allele_coverage_mother.split(',')[1]))) + "\t" 
							break
						elif allele_coverage_mother.split(',')[0] == "0" and allele_coverage_mother.split(',')[1] == "0":
							formato = formato + "0.0" + "\t"
							break
						else:
							formato = formato + "1.0" + "\t"
							break
					for h in lista:
						if allele_coverage_father.split(',') != "0" and '.' in allele_coverage_father:
							formato = formato + ".,."
							break
						elif allele_coverage_father.split(',') != "0" and '-' in allele_coverage_father:
							formato = formato + ".,." 
							break
						elif allele_coverage_father.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_father.split(',')[1])/ (float(allele_coverage_father.split(',')[0]) + float(allele_coverage_father.split(',')[1])))
							break
						elif allele_coverage_father.split(',')[0] == "0" and allele_coverage_father.split(',')[1] == "0":
							formato = formato + "0.0"
							break
						else:
							formato = formato +"1.0"
							break
							
					#Añadimos el origen

					if data_proband[0] == data_mother[0] and data_proband[0] == data_father[0]:
						formato = formato + "\t" + "Both"
					elif data_mother[0] == data_proband[0]:
						formato = formato + "\t" + "Mother"
					elif data_father[0] == data_proband[0]:
						formato = formato + "\t" + "Father"
					else:
						formato = formato + "\t" + "Check"

					trio.append(formato)

		header_latest = "Locus" + "\t" + "Mother Depth" + "\t" + "Father Depth" + "\t" + "Mother AlleleCoverage" + "\t" + "Father AlleleCoverage" + "\t" + "Mother Genotype" + "\t" + "Father Genotype" + "\t" + "zigosity Mother" + "\t" + "zigosity Father" + "\t" + "VAF Mother" + "\t" + "VAF Father" + "\t" +"Origin" + "\n"


		f_intermedio.write(header_latest)

		# Meto los datos de genotipos, profundidades y origen de las variantes en un archivo intermedio

		for line in trio:
			# ~ print line
			f_intermedio.write("%s\n" % line)
		return True

class cuartetoDragen:


	wd=""
	cuartetos=""
	_filter=""
	
	def __init__(self,wd,cuartetos,_filter): #Define __init___ que es el primer método que utiliza python para crear una instancia de cuartetoDragen, luego wd, cuartetos y cuartetos_filter que son el resto de argumentos para definir la instancia self


		self.wd=wd # los argumentos los define como . de self para poder referirse a ellos
		self.cuartetos=cuartetos
		self._filter=_filter


	

	def run(this,caso,padre,madre,parient,sexo): 

		if this.Filtrado(caso,padre,madre,parient,sexo):
			this.Columnas(caso)
			return this.Switch(caso)


	def Switch(this,caso):
	 	#Reordenamiento de columnas

		output_dir=get_output_dir(this.wd,this._filter)

		try:
	 		f_final_file= (output_dir+"/" + str(caso) + "_final.tsv")
	 		f_final = pd.read_csv(f_final_file, sep="\t", header=0)
		except (OSError, IOError) as e:
	 		raise ErrorFinalFile(f_final_file)



		f_finalisimo = (output_dir+"/" + str(caso) + "_finalisimo.tsv")
	
		columns_titles=['annonimous_GENE', 'GEN_OMIM', 'PHENO_OMIM', 'PHENO_DESC', 'INHERITANCE', 'clinvar_id', 'clinvar_CLNSIG', 'clinvar_CLNSIGCONF', 'clinvar_CLNDN', 'annonimous_PREDICTORS', 'CADD_phred', 'occurrence_parents', 'occurrence', 'Locus', 'REF', 'genotype', 'Mother Genotype', 'Father Genotype', 'Parient Genotype', 'Origin', 'zigosity', 'zigosity Mother', 'zigosity Father', 'zigosity Parient', 'gnomad_AF', 'ALL_sites_2015_08', 'esp6500siv2_all', 'VAF', 'VAF Mother', 'VAF Father', 'VAF Parient', 'Coverage', 'Mother Depth', 'Father Depth', 'Parient Depth', 'AlleleRatio', 'AlleleCoverage', 'Mother AlleleCoverage', 'Father AlleleCoverage', 'Parient AlleleCoverage', 'annonimous_ANNOTATION', 'annonimous_Func_refGene', 'annonimous_gNomen', 'annonimous_tCPUNTO', 'reliabilities', 'annonimous_PPUNTO', 'annonimous_EXONS', 'annonimous_INTRONS', 'annonimous_TOTAL', 'distNearestSS', 'nearestSSType', 'dbscSNV_ADA_SCORE', 'dbscSNV_RF_SCORE', 'catalogue_XREF', 'catalogue_VP', 'catalogue_speciality', 'ExClinicoNIM_panel', 'DI', 'SFARI', 'Epilepsia', 'RetNet', 'incidental_findings', 'pseudogenes', 'imprinting', 'PAR_loci', 'dbSNP_id', 'cosmicIds', 'cosmicTissues', 'cosmicSampleCounts', 'ParentsNIMIDs', 'IlluminaProbandsNIMIDs', 'SIFT_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 'LRT_pred', 'MutationTaster_pred', 'MutationAssessor_pred', 'FATHMM_pred', 'MetaSVM_pred', 'gnomad_AF_afr', 'gnomad_AF_amr', 'gnomad_AF_asj', 'gnomad_AF_eas', 'gnomad_AF_sas', 'gnomad_AF_nfe', 'gnomad_AF_fin', 'gnomad_AF_oth', 'gnomad_HetCount_all', 'gnomad_HomCount_all', 'annonimous_HemCount_all', 'gnomad_Filter', 'oe_mis', 'oe_syn', 'pLI', 'oe_lof', 'oe_syn_lower', 'oe_syn_upper', 'oe_mis_lower', 'oe_mis_upper', 'oe_lof_lower', 'oe_lof_upper', 'syn_z', 'mis_z', 'lof_z', 'AFR_sites_2015_08', 'SAS_sites_2015_08', 'EAS_sites_2015_08', 'EUR_sites_2015_08', 'AMR_sites_2015_08', 'esp6500siv2_ea', 'esp6500siv2_aa', 'phyloP100way_vertebrate', 'phyloP20way_mammalian', 'phyloP46way_placental', 'phastCons100way_vertebrate', 'phastCons20way_mammalian', 'APPROVED_NAME', 'SSF_wt', 'SSF_mut', 'MES_wt', 'MES_mut', 'delta_MES', 'delta_SSF', 'SPiCEprobability', 'SPiCEinter_2thr']
	
		f_final = f_final.reindex(columns=columns_titles)
		f_final = f_final.reset_index(drop=True)
		f_final.to_csv(f_finalisimo, sep='\t')  



	def Columnas(this,caso):
		# Abrimos archivos

		path_latest=get_path_latest(caso,this.wd,this._filter)
		my_path_tmp=get_path_tmp(this.wd,this._filter)
		output_dir=get_output_dir(this.wd,this._filter)

		try:
			f = open(path_latest)

		except (OSError, IOError) as e:
			raise ErrorInputFileLatest(path_latest)

		try:
			tmp_file=str(my_path_tmp) + "/cuarteto_intermedio_" + str(caso) + ".tsv"
			f_intermedio = open(tmp_file)
		except (OSError, IOError) as e:
			raise ErrorTempFile(tmp_file)

		f_final = open(output_dir+"/" + str(caso) + "_final.tsv", "w")

		locusIndex = LocusIndex(f)

		f.close()

		f = open(path_latest)

		# AddColummns de intermedio a final (latest + intermedio)

		cabecera = f.readline().split("\n")[0] #Lee la cabecera y se queda con la primera linea
		cabecera_final = cabecera.split("\t") #Parte la cabecera de latest por los tabuladores
		cabecera2 = f_intermedio.readline().split("\t") #Parte la cabecera del intermedio por los tabuladores
		# print (cabecera_final)
		# print (cabecera2)
		# cabecera.remove('\n')
		cabecera2.pop(0) #De la cabecera2 quita la primera posicion
		# print (cabecera)
		# print (cabecera2)

		rows = f.readlines() #Lee el latest tsv linea por linea
		rows2 = f_intermedio.readlines() #Lee el intermedio linea por linea
		rows3 = []

		f.close()
		f_intermedio.close()

		formato = ""
		formato2 = ""

		for i in range(len(cabecera_final)):
			formato2 = formato2 + cabecera_final[i] + "\t" #En formato2 guarda la cabecera final separada por tabuladores
		for j in range(len(cabecera2)):
			if j == len(cabecera2) - 1:
				formato2 = formato2 + cabecera2[j]
			else:
				formato2 = formato2 + cabecera2[j] + "\t"

		rows3.append(formato2) #Fusiona las filas

		# Añadimos IGV  a las lineas en las que no se han encontrado datos en los padres en el archivo trio

		for row in rows:
			prueba = row.split("\t\n") #Parte las lineas por tabulador y salto de línea y lo guarda en la variable prueba
			intermedio = prueba[0].split('\n') #Parte el primer argumento de prueba por salto de linea
			fields = intermedio[0].split('\t') #Parte el primer argumento de intermedio por los tabuladores
			Boo2 = True
			formato = ""
			for row2 in rows2:
				campos = row2.split("\t")
				if campos[0] == fields[locusIndex]:
					campos.pop(0)
					for i in range(len(fields)):
						formato = formato + fields[i] + "\t"
					for j in range(len(campos)):
						if j == len(campos) - 1:
							formato = formato + campos[j]
						else:
							formato = formato + campos[j] + "\t"
					rows3.append(formato)
					Boo2 = False
					break
				else:
					continue
			while Boo2 == True:
				for i in range(len(fields)):
					formato = formato + fields[i] + "\t"
				formato = formato + "IGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\n"
				rows3.append(formato)
				Boo2 = False

		for line in rows3:
			f_final.write(line)

		f_final.close()


	def Filtrado(this,caso,padre,madre,parient,sexo):
		# Abrimos archivos

		# my_path = os.getcwd()

		# Busca el latest en samples

		path_latest=get_path_latest(caso,this.wd,this._filter)
		cuarteto_file=get_path_cuarteto_file(caso,this.cuartetos)
		my_path_tmp=get_path_tmp(this.wd,this._filter)

		try:
			f = open(path_latest) #Guardamos en f erl archivo de la ruta latest
		except (OSError, IOError) as e:
			raise ErrorInputFileLatest(path_latest)

		# Busca el vcf del cuarteto

		try:
			f_cuarteto = open(cuarteto_file)

		except (OSError, IOError) as e:
			raise ErrorInputFileCuarteto(cuarteto_file)

		f_intermedio = open(str(my_path_tmp) + "/cuarteto_intermedio_" + str(caso) + ".tsv", "w+") #Guardamos en f_intermedio el archivo intermedio

		locusIndex = LocusIndex(f) #Con la función LocusIndex genero el index del Locus

		locus = []
		cuarteto = []
		lines = f.readlines() #Guardo las lineas del archivo latest en la variable lines

		f.close()

		for line in lines:
			fields = line.split('\t') #Guardo en la variable fields los campos de las lineas del archivo latest
			locus.append(fields[locusIndex]) #Se añade el index del locus a los campos del archivo latest

		i = 0

		find_header=False;
		header_proband=None

		while(not find_header): #Indica como header la primera linea del vcf que comienza con CHROM
			header_proband = f_cuarteto.readline()
			if header_proband.startswith("#CHROM"):
				find_header=True


		try:
			probandIndex,fatherIndex,motherIndex,parientIndex=getIndexFamily_4(caso,padre,madre,parient,header_proband) #Guarda el indice de cada miembro de la familia
		except NotFoundFamily as e:
			this.__log.error(e)
			return False

		lines_cuarteto = f_cuarteto.readlines() 

		f_cuarteto.close()

		for line in lines_cuarteto:
			fields_cuarteto = line.split("\t") #Parte las lineas del cuarteto por los tabuladores generano los campos
			data_proband = fields_cuarteto[probandIndex].split('\n')[0].split(":") #Guarda en la variable data_proband campos del cuarteto del probando partidos por salto de línea y :
			data_mother = fields_cuarteto[motherIndex].split('\n')[0].split(":") #Lo mismo con la madre
			data_father = fields_cuarteto[fatherIndex].split('\n')[0].split(":") #Lo mismo con el padre
			data_parient = fields_cuarteto[parientIndex].split('\n')[0].split(":") #Lo mismo con el pariente adicional
			
			#Prepara los datos para genotipar a los miembros de la familia
			
			if str(data_mother[0]) == '1' or str(data_mother[0]) == '0' or str(data_mother[0]) == '2':
				data_mother[0] = str(data_mother[0]) + '/' + str(data_mother[0])
			
			if str(data_father[0]) == '1' or str(data_father[0]) == '0' or str(data_father[0]) == '2':
				data_father[0] = str(data_father[0]) + '/' + str(data_father[0])

			if str(data_proband[0]) == '1' or str(data_proband[0]) == '0' or str(data_proband[0]) == '2':
				data_proband[0] = str(data_proband[0]) + '/' + str(data_proband[0])
			
			if str(data_parient[0]) == '1' or str(data_parient[0]) == '0' or str(data_parient[0]) == '2':
				data_parient[0] = str(data_parient[0]) + '/' + str(data_parient[0])


			if '/' in data_proband[0]:
				genotype_proband = data_proband[0].split("/")
				genotype_proband_0 = genotype_proband[0]
				genotype_proband_1 = genotype_proband[1]
			elif '|' in data_proband[0]:
				genotype_proband = data_proband[0].split("|")
				genotype_proband_0 = genotype_proband[0]
				genotype_proband_1 = genotype_proband[1]

			if '/' in data_mother[0]:
				genotype_mother = data_mother[0].split("/")
				genotype_mother_0 = genotype_mother[0]
				genotype_mother_1 = genotype_mother[1]
			elif '|' in data_mother[0]:
				genotype_mother = data_mother[0].split("|")
				genotype_mother_0 = genotype_proband[0]
				genotype_mother_1 = genotype_proband[1]

			if '/' in data_father[0]:
				genotype_father = data_father[0].split("/")
				genotype_father_0 = genotype_father[0]
				genotype_father_1 = genotype_father[1]
			elif '|' in data_father[0]:
				genotype_father = data_father[0].split("|")
				genotype_father_0 = genotype_father[0]
				genotype_father_1 = genotype_father[1]
				
			if '/' in data_parient[0]:
				genotype_parient = data_parient[0].split("/")
				genotype_parient_0 = genotype_parient[0]
				genotype_parient_1 = genotype_parient[1]
			elif '|' in data_parient[0]:
				genotype_parient = data_parient[0].split("|")
				genotype_parient_0 = genotype_parient[0]
				genotype_parient_1 = genotype_parient[1]

			#Guarda la cobertura de los miembros de la familia si la tienen
			if len(data_mother) < 3: #Si faltan datos, introducirlos
				allele_coverage_mother = '-,-'
			elif len(data_mother[1].split(',')) == 1 and genotype_mother != './.': 
				allele_coverage_mother = data_mother[1] + ',0'
			else:
				allele_coverage_mother = data_mother[1] 

			if len(data_father) < 3:
				allele_coverage_father = '-,-'
			elif len(data_father[1].split(',')) == 1 and genotype_father != './.':
				allele_coverage_father = data_father[1] + ',0'
			else:
				allele_coverage_father = data_father[1]

			if len(data_parient) < 3:
				allele_coverage_parient = '-,-'
			elif len(data_parient[1].split(',')) == 1 and genotype_parient != './.':
				allele_coverage_parient = data_parient[1] + ',0'
			else:
				allele_coverage_parient = data_parient[1]

			ref_split = fields_cuarteto[4].split(",")
			lista = range(len(ref_split))

			# HACERLO INDEXANDO VALORES Y QUE BUSQUE EN REF CON EL MISMO INDEX

			for field in locus:
				fields_locus = field.split(":")
				if fields_cuarteto[0] == fields_locus[0] and fields_cuarteto[1] == fields_locus[1]:
					for h in lista:
						#print field #TODAS ESTAS RAMAS DE IF AÑADEN LA INFO DE COBERTURA Y GENOTIPO
						if len(data_mother) < 3:
							if len(data_father) < 3:
								if len(data_parient) < 3:
									formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + '-' + "\t" + '-' + "\t" + '-' + "\t" + \
											  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
											  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
											  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
											  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
											  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
											  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									break
								else:
									if data_parient[3] == '.':
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + '-' + "\t" +  '-' "\t" + \
												  data_parient[3] + \
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									else:
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + '-' + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + \
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1] 
									break
							else:
								if data_father[3] == '.':
									if len(data_parient) < 3: 
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + '-' + "\t" + \
												  data_father[3] + "\t" + '-' +\
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									
									elif data_parient [3] == '.':
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + '-' + "\t" + \
												  data_father[3] + "\t" + data_parient[3] +\
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									else: 
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + '-' + "\t" + data_father[3] + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + \
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									break
								else:
									if len(data_parient) < 3: 
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  allele_coverage_father.split(',')[1]))  + "\t" + '-' +\
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									
									elif data_parient [3] == '.':
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  allele_coverage_father.split(',')[1])) + "\t" + data_parient[3] + \
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									else: 
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  allele_coverage_father.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + \
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1] 
									break
						else:
							if len(data_father) < 3:
								if data_mother[3] == '.':
									if len(data_parient) < 3: 
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + data_mother[3] + "\t" + \
												  '-' + "\t" + '-' +\
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									
									elif data_parient [3] == '.':
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + data_mother[3] + "\t" + \
												  '-' + "\t" + data_parient[3] +\
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									else: 
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + \
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									break
								else:
									if len(data_parient) < 3: 
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  allele_coverage_mother.split(',')[1])) + "\t" + '-'  + "\t" + '-' +\
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									
									elif data_parient [3] == '.':
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  allele_coverage_mother.split(',')[1])) + "\t" + '-' + "\t" + data_parient[3] + \
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									else: 
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  allele_coverage_mother.split(',')[1])) + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + \
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1] 
									break
							else:
								if data_mother[3] == '.' and data_father[3] == '.':
									if len(data_parient) < 3: 
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + data_mother[3] + "\t" + \
												  data_father[3] + "\t" + '-' +\
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									
									elif data_parient [3] == '.':
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + data_mother[3] + "\t" + \
												  data_father[3] + "\t" + data_parient[3] +\
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									else: 
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + data_mother[3] + "\t" + \
												  data_father[3] + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + \
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									break
								
								elif data_father[3] == '.':
									if len(data_parient) < 3: 
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  allele_coverage_mother.split(',')[1])) + "\t" + data_father[3] + "\t" + '-' +\
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									
									elif data_parient [3] == '.':
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  allele_coverage_mother.split(',')[1])) + "\t" + data_father[3] + "\t" + data_parient[3] + \
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									else: 
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  allele_coverage_mother.split(',')[1])) + "\t" + data_father[3] + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + \
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1] 
									break

								elif data_mother[3] == '.':
								
									if len(data_parient) < 3: 
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  allele_coverage_father.split(',')[1]))  + "\t" + '-' +\
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									
									elif data_parient [3] == '.':
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  allele_coverage_father.split(',')[1])) + "\t" + data_parient[3] + \
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									else: 
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  allele_coverage_father.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + \
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1] 
									break

								else:
									if len(data_parient) < 3: 
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + "\t" + '-' + \
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									
									elif data_parient [3] == '.':
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + "\t" + data_parient[3] +\
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									else: 
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1]))+ "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + \
												  "\t" + fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[
													  0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_parient.split(',')[1]
									break
								break


					for h in lista:

						if genotype_mother[0] == "0":
							formato = formato + "\t" + fields_cuarteto[3] + "/"
							break
						elif genotype_mother[0] == str(h + 1):
							formato = formato + "\t" + ref_split[h] + "/"
							break
						elif genotype_mother[0] == ".":
							formato = formato + "\t" + genotype_mother[0] + "/"
							break					
					for h in lista:

						if genotype_mother[1] == "0":
							formato = formato + fields_cuarteto[3] + "\t"
							break
						elif genotype_mother[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_mother[1] == ".":
							formato = formato + genotype_mother[1] + "\t"
							break
					for h in lista:
						if genotype_father[0] == "0":
							formato = formato + fields_cuarteto[3] + "/"
							break
						elif genotype_father[0] == str(h + 1):
							formato = formato + ref_split[h] + "/"
							break
						elif genotype_father[0] == ".":
							formato = formato + genotype_father[0] + "/"
							break
					for h in lista:
						if genotype_father[1] == "0":
							formato = formato + fields_cuarteto[3] + "\t"
							break
						elif genotype_father[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_father[1] == ".":
							formato = formato + genotype_father[1] + "\t"
							break
					for h in lista:
						if genotype_parient[0] == "0":
							formato = formato + fields_cuarteto[3] + "/"
							break
						elif genotype_parient[0] == str(h + 1):
							formato = formato + ref_split[h] + "/"
							break
						elif genotype_parient[0] == ".":
							formato = formato + genotype_parient[0] + "/"
							break
					for h in lista:
						if genotype_parient[1] == "0":
							formato = formato + fields_cuarteto[3] + "\t"
							break
						elif genotype_parient[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_parient[1] == ".":
							formato = formato + genotype_parient[1] + "\t"
							break
							
					#Añadimos la zigosidad
					
					for h in lista:
						if genotype_mother[0] == "0" and genotype_mother[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_mother[0] == "0" and genotype_mother[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t" 
							break
						elif genotype_mother[0] == str(h + 1) and genotype_mother[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_mother[0] == "." or genotype_mother[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
							
					for h in lista:
						if genotype_father[0] == "0" and genotype_father[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_father[0] == "0" and genotype_father[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t" 
							break
						elif genotype_father[0] == str(h + 1) and genotype_father[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_father[0] == "." or genotype_father[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
							
					for h in lista:
						if genotype_parient[0] == "0" and genotype_parient[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_parient[0] == "0" and genotype_parient[1] == str(h + 1):
							formato = formato + "Heteroz"  + "\t"
							break
						elif genotype_parient[0] == str(h + 1) and genotype_parient[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_parient[0] == "." or genotype_parient[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
							
					#Añadimos el VAF
					
					for h in lista:
						if allele_coverage_mother.split(',') != "0" and '.' in allele_coverage_mother:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_mother.split(',') != "0" and '-' in allele_coverage_mother:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_mother.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_mother.split(',')[1])/ (float(allele_coverage_mother.split(',')[0]) + float(allele_coverage_mother.split(',')[1]))) + "\t"
							break
						elif allele_coverage_mother.split(',')[0] == "0" and allele_coverage_mother.split(',')[1] == "0":
							formato = formato + "0.0" + "\t"
							break
						else:
							formato = formato + "1.0" + "\t"
							break
					for h in lista:
						if allele_coverage_father.split(',') != "0" and '.' in allele_coverage_father:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_father.split(',') != "0" and '-' in allele_coverage_father:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_father.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_father.split(',')[1])/ (float(allele_coverage_father.split(',')[0]) + float(allele_coverage_father.split(',')[1]))) + "\t"
							break
						elif allele_coverage_father.split(',')[0] == "0" and allele_coverage_father.split(',')[1] == "0":
							formato = formato + "0.0" + "\t"
							break
						else:
							formato = formato + "1.0" + "\t"
							break
					for h in lista:
						if allele_coverage_parient.split(',') != "0" and '.' in allele_coverage_parient:
							formato = formato + ".,." 
							break
						elif allele_coverage_parient.split(',') != "0" and '-' in allele_coverage_parient:
							formato = formato + ".,." 
							break
						elif allele_coverage_parient.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_parient.split(',')[1])/ (float(allele_coverage_parient.split(',')[0]) + float(allele_coverage_parient.split(',')[1]))) 
							break
						elif allele_coverage_parient.split(',')[0] == "0" and allele_coverage_parient.split(',')[1] == "0":
							formato = formato + "0.0" 
							break
						else:
							formato = formato + "1.0" 
							break
							
					#Añadimos Origen

					if data_proband[0] == data_mother[0] and data_proband[0] == data_father[0] and data_proband[0] == data_parient[0]:
						formato = formato + "\t" + "All"
					elif data_parient[0] == data_proband[0] and data_mother[0] == data_proband[0]:
						formato = formato + "\t" + "Parient\Mother"
					elif data_parient[0] == data_proband[0] and data_father[0] == data_proband[0]:
						formato = formato + "\t" + "Parient\Father"
					elif data_mother[0] == data_proband[0] and data_father[0] == data_proband[0]:
						formato = formato + "\t" + "Mother\Father"
					elif data_mother[0] == data_proband[0]:
						formato = formato + "\t" + "Mother"
					elif data_father[0] == data_proband[0]:
						formato = formato + "\t" + "Father"
					elif data_parient[0] == data_proband[0]:
						formato = formato + "\t" + "Parient"
					else:
						formato = formato + "\t" + "Check"

					cuarteto.append(formato)

		header_latest = "Locus" + "\t" + "Mother Depth" + "\t" + "Father Depth" + "\t" + "Parient Depth" +"\t" + "Mother AlleleCoverage" + "\t" + "Father AlleleCoverage" + "\t" + "Parient AlleleCoverage" + "\t" + "Mother Genotype" + "\t" + "Father Genotype" + "\t" + "Parient Genotype" + "\t" + "zigosity Mother" + "\t" + "zigosity Father" + "\t" + "zigosity Parient" + "\t" + "VAF Mother" + "\t" + "VAF Father" + "\t" + "VAF Parient" + "\t" + "Origin" + "\n"


		f_intermedio.write(header_latest)

		# Meto los datos de genotipos, profundidades y origen de las variantes en un archivo intermedio

		for line in cuarteto:
			# ~ print line
			f_intermedio.write("%s\n" % line)
		return True
		
class cuartetoDragen_2P:

	wd=""
	cuartetos=""
	_filter=""
	
	def __init__(self,wd,cuartetos,_filter): #Define __init___ que es el primer método que utiliza python para crear una instancia de cuartetoDragen_2P, luego wd, cuartetos y cuartetos_filter que son el resto de argumentos para definir la instancia self


		self.wd=wd # los argumentos los define como . de self para poder referirse a ellos
		self.cuartetos=cuartetos
		self._filter=_filter


	

	def run(this,caso,caso2,father,mother,sexo,sexo2): 

		if this.Filtrado(caso,caso2,father,mother,sexo,sexo2):
			this.Columnas(caso,caso2)
			return this.Switch(caso, caso2)


	def Switch(this,caso,caso2):
	 	#Reordenamiento de columnas

		output_dir=get_output_dir(this.wd,this._filter)

		try:
	 		f_final_file= (output_dir+"/" + str(caso) + '+' + str(caso2) + "_final.tsv")
	 		f_final = pd.read_csv(f_final_file, sep="\t", header=0)
		except (OSError, IOError) as e:
	 		raise ErrorFinalFile(f_final_file)



		f_finalisimo = (output_dir+"/" + str(caso) + '+' + str(caso2) + "_finalisimo.tsv")

		f_final.drop(['genotype', 'zigosity', 'VAF', 'Coverage', 'AlleleRatio'], axis=1, inplace=True)
	
		columns_titles=['annonimous_GENE', 'GEN_OMIM', 'PHENO_OMIM', 'PHENO_DESC', 'INHERITANCE', 'clinvar_id', 'clinvar_CLNSIG', 'clinvar_CLNSIGCONF', 'clinvar_CLNDN', 'annonimous_PREDICTORS', 'CADD_phred', 'occurrence_parents', 'occurrence', 'Locus', 'REF', str(caso) + ' Genotype', str(caso2) + ' Genotype', 'Mother Genotype', 'Father Genotype', 'Origin', 'zigosity ' + str(caso), 'zigosity ' + str(caso2), 'zigosity Mother', 'zigosity Father', 'gnomad_AF', 'ALL_sites_2015_08', 'esp6500siv2_all', 'VAF ' +  str(caso), 'VAF ' +  str(caso2), 'VAF Mother', 'VAF Father', str(caso) + ' Coverage', str(caso2) + ' Coverage', 'Mother Depth', 'Father Depth', str(caso) + ' AlleleCoverage', str(caso2) + ' AlleleCoverage', 'Mother AlleleCoverage', 'Father AlleleCoverage', 'annonimous_ANNOTATION', 'annonimous_Func_refGene', 'annonimous_gNomen', 'annonimous_tCPUNTO', 'reliabilities', 'annonimous_PPUNTO', 'annonimous_EXONS', 'annonimous_INTRONS', 'annonimous_TOTAL', 'distNearestSS', 'nearestSSType', 'dbscSNV_ADA_SCORE', 'dbscSNV_RF_SCORE', 'catalogue_XREF', 'catalogue_VP', 'catalogue_speciality', 'ExClinicoNIM_panel', 'DI', 'SFARI', 'Epilepsia', 'RetNet', 'incidental_findings', 'pseudogenes', 'imprinting', 'PAR_loci', 'dbSNP_id', 'cosmicIds', 'cosmicTissues', 'cosmicSampleCounts', 'ParentsNIMIDs', 'IlluminaProbandsNIMIDs', 'SIFT_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 'LRT_pred', 'MutationTaster_pred', 'MutationAssessor_pred', 'FATHMM_pred', 'MetaSVM_pred', 'gnomad_AF_afr', 'gnomad_AF_amr', 'gnomad_AF_asj', 'gnomad_AF_eas', 'gnomad_AF_sas', 'gnomad_AF_nfe', 'gnomad_AF_fin', 'gnomad_AF_oth', 'gnomad_HetCount_all', 'gnomad_HomCount_all', 'annonimous_HemCount_all', 'gnomad_Filter', 'oe_mis', 'oe_syn', 'pLI', 'oe_lof', 'oe_syn_lower', 'oe_syn_upper', 'oe_mis_lower', 'oe_mis_upper', 'oe_lof_lower', 'oe_lof_upper', 'syn_z', 'mis_z', 'lof_z', 'AFR_sites_2015_08', 'SAS_sites_2015_08', 'EAS_sites_2015_08', 'EUR_sites_2015_08', 'AMR_sites_2015_08', 'esp6500siv2_ea', 'esp6500siv2_aa', 'phyloP100way_vertebrate', 'phyloP20way_mammalian', 'phyloP46way_placental', 'phastCons100way_vertebrate', 'phastCons20way_mammalian', 'APPROVED_NAME', 'SSF_wt', 'SSF_mut', 'MES_wt', 'MES_mut', 'delta_MES', 'delta_SSF', 'SPiCEprobability', 'SPiCEinter_2thr']
	
		f_final = f_final.reindex(columns=columns_titles)
		f_final = f_final.reset_index(drop=True)
		f_final.to_csv(f_finalisimo, sep='\t')  



	def Columnas(this,caso,caso2):
		
		# Abrimos archivos

		path_latest=get_path_latest_2P(caso,caso2,this.wd,this._filter)
		my_path_tmp=get_path_tmp(this.wd,this._filter)
		output_dir=get_output_dir(this.wd,this._filter)

		try:
			f = open(path_latest)

		except (OSError, IOError) as e:
			raise ErrorInputFileLatest(path_latest)

		try:
			tmp_file=str(my_path_tmp) + "/cuarteto_intermedio_" + str(caso) + "+" + str(caso2) + ".tsv"
			f_intermedio = open(tmp_file)
		except (OSError, IOError) as e:
			raise ErrorTempFile(tmp_file)

		f_final = open(output_dir+"/" + str(caso) + '+' + str(caso2) + "_final.tsv", "w")

		locusIndex = LocusIndex(f)

		f.close()

		f = open(path_latest)

		# AddColummns de intermedio a final (latest + intermedio)

		cabecera = f.readline().split("\n")[0] #Lee la cabecera y se queda con la primera linea
		cabecera_final = cabecera.split("\t") #Parte la cabecera de latest por los tabuladores
		cabecera2 = f_intermedio.readline().split("\t") #Parte la cabecera del intermedio por los tabuladores
		# print (cabecera_final)
		# print (cabecera2)
		# cabecera.remove('\n')
		cabecera2.pop(0) #De la cabecera2 quita la primera posicion
		# print (cabecera)
		# print (cabecera2)

		rows = f.readlines() #Lee el latest tsv linea por linea
		rows2 = f_intermedio.readlines() #Lee el intermedio linea por linea
		rows3 = []

		f.close()
		f_intermedio.close()

		formato = ""
		formato2 = ""

		for i in range(len(cabecera_final)):
			formato2 = formato2 + cabecera_final[i] + "\t" #En formato2 guarda la cabecera final separada por tabuladores
		for j in range(len(cabecera2)):
			if j == len(cabecera2) - 1:
				formato2 = formato2 + cabecera2[j]
			else:
				formato2 = formato2 + cabecera2[j] + "\t"

		rows3.append(formato2) #Fusiona las filas

		# Añadimos IGV  a las lineas en las que no se han encontrado datos en los padres en el archivo trio

		for row in rows:
			prueba = row.split("\t\n") #Parte las lineas por tabulador y salto de línea y lo guarda en la variable prueba
			intermedio = prueba[0].split('\n') #Parte el primer argumento de prueba por salto de linea
			fields = intermedio[0].split('\t') #Parte el primer argumento de intermedio por los tabuladores
			Boo2 = True
			formato = ""
			for row2 in rows2:
				campos = row2.split("\t")
				if campos[0] == fields[locusIndex]:
					campos.pop(0)
					for i in range(len(fields)):
						formato = formato + fields[i] + "\t"
					for j in range(len(campos)):
						if j == len(campos) - 1:
							formato = formato + campos[j]
						else:
							formato = formato + campos[j] + "\t"
					rows3.append(formato)
					Boo2 = False
					break
				else:
					continue
			while Boo2 == True:
				for i in range(len(fields)):
					formato = formato + fields[i] + "\t"
				formato = formato + "IGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\n"
				rows3.append(formato)
				Boo2 = False

		for line in rows3:
			f_final.write(line)

		f_final.close()


	def Filtrado(this,caso,caso2,padre,madre,sexo,sexo2):
		# Abrimos archivos

		#Generamos el archivo latest común
		
		path_latest_1=get_path_latest(caso,this.wd,this._filter)
		path_latest_2=get_path_latest(caso2,this.wd,this._filter)
		
		latest_1 = pd.read_csv(path_latest_1, sep="\t", header=0)
		latest_2 = pd.read_csv(path_latest_2, sep="\t", header=0)
		
		joint_latest = pd.concat([latest_1,latest_2], ignore_index=True).drop_duplicates(['Locus'])
		
		joint_latest[['chr','position']]=joint_latest.Locus.str.split(':',expand=True)
		joint_latest['chr']=joint_latest.chr.str.replace('chr','')
		joint_latest['chr']=joint_latest.chr.str.replace('X','30')
		joint_latest['chr']=joint_latest.chr.str.replace('Y','31')
		joint_latest['chr']=joint_latest.chr.str.replace('M','32')
		joint_latest['chr']=joint_latest['chr'].astype(int)
		joint_latest['position']=joint_latest['position'].astype(int)

		joint_latest_sort = joint_latest.sort_values(by=['chr','position'])
		joint_latest_sort = joint_latest_sort.drop(['chr','position'], axis=1)
		joint_latest_sort = joint_latest_sort.reset_index(drop=True)
		
		output_dir=get_output_dir(this.wd,this._filter)
		
		file_joint_latest = (output_dir+"/" + str(caso) + '+' + str(caso2) + "_latest.tsv")
		
		joint_latest_sort.to_csv(file_joint_latest, sep='\t')   

		# Busca el latest 

		path_latest=get_path_latest_2P(caso,caso2,this.wd,this._filter)
		cuarteto_file=get_path_cuarteto_file_2P(caso,caso2,this.cuartetos)
		my_path_tmp=get_path_tmp(this.wd,this._filter)

		try:
			f = open(path_latest) #Guardamos en f erl archivo de la ruta latest
		except (OSError, IOError) as e:
			raise ErrorInputFileLatest(path_latest)

		# Busca el vcf cuarteto

		try:
			f_cuarteto = open(cuarteto_file) #Guardamos en f cuarteto el archivo vcf cuarteto

		except (OSError, IOError) as e:
			raise ErrorInputFileCuarteto(cuarteto_file)

		f_intermedio = open(str(my_path_tmp) + "/cuarteto_intermedio_" + str(caso) + "+" + str(caso2) + ".tsv", "w+") #Guardamos en f_intermedio el archivo intermedio

		locusIndex = LocusIndex(f) #Con la función LocusIndex genero el index del Locus

		locus = []
		cuarteto = []
		lines = f.readlines() #Guardo las lineas del archivo latest en la variable lines

		f.close()

		for line in lines:
			fields = line.split('\t') #Guardo en la variable fields los campos de las lineas del archivo latest
			locus.append(fields[locusIndex]) #Se añade el index del locus a los campos del archivo latest

		i = 0

		find_header=False;
		header_proband=None

		while(not find_header): #Indica como header la primera linea del vcf que comienza con CHROM
			header_proband = f_cuarteto.readline()
			if header_proband.startswith("#CHROM"):
				find_header=True


		try:
			probandIndex,proband2Index,fatherIndex,motherIndex=getIndexFamily_4_2P(caso,caso2,padre,madre,header_proband) #Guarda el indice de cada miembro de la familia
		except NotFoundFamily as e:
			this.__log.error(e)
			return False

		lines_cuarteto = f_cuarteto.readlines() 

		f_cuarteto.close()

		for line in lines_cuarteto:
			fields_cuarteto = line.split("\t") 
			data_proband = fields_cuarteto[probandIndex].split('\n')[0].split(":") 
			data_mother = fields_cuarteto[motherIndex].split('\n')[0].split(":") 
			data_father = fields_cuarteto[fatherIndex].split('\n')[0].split(":") 
			data_proband2 = fields_cuarteto[proband2Index].split('\n')[0].split(":") 
			
			#Prepara los datos para genotipar a los miembros de la familia
			
			if str(data_mother[0]) == '1' or str(data_mother[0]) == '0' or str(data_mother[0]) == '2':
				data_mother[0] = str(data_mother[0]) + '/' + str(data_mother[0])
			
			if str(data_father[0]) == '1' or str(data_father[0]) == '0' or str(data_father[0]) == '2':
				data_father[0] = str(data_father[0]) + '/' + str(data_father[0])

			if str(data_proband[0]) == '1' or str(data_proband[0]) == '0' or str(data_proband[0]) == '2':
				data_proband[0] = str(data_proband[0]) + '/' + str(data_proband[0])
			
			if str(data_proband2[0]) == '1' or str(data_proband2[0]) == '0' or str(data_proband2[0]) == '2':
				data_proband2[0] = str(data_proband2[0]) + '/' + str(data_proband2[0])


			if '/' in data_proband[0]:
				genotype_proband = data_proband[0].split("/")
				genotype_proband_0 = genotype_proband[0]
				genotype_proband_1 = genotype_proband[1]
			elif '|' in data_proband[0]:
				genotype_proband = data_proband[0].split("|")
				genotype_proband_0 = genotype_proband[0]
				genotype_proband_1 = genotype_proband[1]

			if '/' in data_mother[0]:
				genotype_mother = data_mother[0].split("/")
				genotype_mother_0 = genotype_mother[0]
				genotype_mother_1 = genotype_mother[1]
			elif '|' in data_mother[0]:
				genotype_mother = data_mother[0].split("|")
				genotype_mother_0 = genotype_proband[0]
				genotype_mother_1 = genotype_proband[1]

			if '/' in data_father[0]:
				genotype_father = data_father[0].split("/")
				genotype_father_0 = genotype_father[0]
				genotype_father_1 = genotype_father[1]
			elif '|' in data_father[0]:
				genotype_father = data_father[0].split("|")
				genotype_father_0 = genotype_father[0]
				genotype_father_1 = genotype_father[1]
				
			if '/' in data_proband2[0]:
				genotype_proband2 = data_proband2[0].split("/")
				genotype_proband2_0 = genotype_proband2[0]
				genotype_proband2_1 = genotype_proband2[1]
			elif '|' in data_proband2[0]:
				genotype_proband2 = data_proband2[0].split("|")
				genotype_proband2_0 = genotype_proband2[0]
				genotype_proband2_1 = genotype_proband2[1]

			#Guarda la cobertura de los miembros de la familia si la tienen
			if len(data_mother) < 3: #Si faltan datos, introducirlos
				allele_coverage_mother = '-,-'
			elif len(data_mother[1].split(',')) == 1 and genotype_mother != './.': 
				allele_coverage_mother = data_mother[1] + ',0'
			else:
				allele_coverage_mother = data_mother[1] 

			if len(data_father) < 3:
				allele_coverage_father = '-,-'
			elif len(data_father[1].split(',')) == 1 and genotype_father != './.':
				allele_coverage_father = data_father[1] + ',0'
			else:
				allele_coverage_father = data_father[1]

			if len(data_proband) < 3:
				allele_coverage_proband = '-,-'
			elif len(data_proband[1].split(',')) == 1 and genotype_proband != './.':
				allele_coverage_proband = data_proband[1] + ',0'
			else:
				allele_coverage_proband = data_proband[1]
				
			if len(data_proband2) < 3:
				allele_coverage_proband2 = '-,-'
			elif len(data_proband2[1].split(',')) == 1 and genotype_proband2 != './.':
				allele_coverage_proband2 = data_proband2[1] + ',0'
			else:
				allele_coverage_proband2 = data_proband2[1]

			ref_split = fields_cuarteto[4].split(",")
			lista = range(len(ref_split))

			# HACERLO INDEXANDO VALORES Y QUE BUSQUE EN REF CON EL MISMO INDEX

			for field in locus:
				fields_locus = field.split(":")
				if fields_cuarteto[0] == fields_locus[0] and fields_cuarteto[1] == fields_locus[1]:
					for h in lista:
						#print field #TODAS ESTAS RAMAS DE IF AÑADEN LA INFO DE COBERTURA Y GENOTIPO
							if len(data_proband) < 3:
								if len(data_mother) < 3:
									if len(data_father) < 3:
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int( \
												  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
										break
									else:
										if data_father[3] == '.':
											formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int( \
													  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + \
													  data_father[3] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
										else:
											formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int( \
													  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int( \
													  allele_coverage_father.split(',')[1])) + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] 
										break
								else:
									if data_mother[3] == '.':
										if len(data_father) < 3:
											formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int( \
													  allele_coverage_proband2.split(',')[1])) + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
											break
										
										elif data_father[3] == '.':
											formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int( \
													  allele_coverage_proband2.split(',')[1])) + "\t" + data_mother[3] + "\t" + \
													  data_father[3] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
										else: 
											formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int( \
													  allele_coverage_proband2.split(',')[1])) + "\t" + data_mother[3] + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int( \
													  allele_coverage_father.split(',')[1])) + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] 
										break
									else:
										if len(data_father) < 3: 
											formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int( \
													  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int( \
													  allele_coverage_mother.split(',')[1])) + "\t" + '-' + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
										
										elif data_father [3] == '.':
											formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int( \
													  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int( \
													  allele_coverage_mother.split(',')[1])) + "\t" + data_father[3] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
										else: 
											formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int( \
													  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int( \
													  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int( \
													  allele_coverage_father.split(',')[1])) + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] 
										break
							elif data_proband[3] == '.':
								if len(data_mother) < 3:
									if len(data_father) < 3:
										formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int( \
												  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
												  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
										break
									else:
										if data_father[3] == '.':
											formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int( \
													  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + \
													  data_father[3] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
										else:
											formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int( \
													  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int( \
													  allele_coverage_father.split(',')[1])) + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] 
										break
								else:
									if data_mother[3] == '.':
										if len(data_father) < 3:
											formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int( \
													  allele_coverage_proband2.split(',')[1])) + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
											break
										elif data_father[3] == '.':
											formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int( \
													  allele_coverage_proband2.split(',')[1])) + "\t" + data_mother[3] + "\t" + \
													  data_father[3] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
										else: 
											formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int( \
													  allele_coverage_proband2.split(',')[1])) + "\t" + data_mother[3] + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int( \
													  allele_coverage_father.split(',')[1])) + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] 
										break
									else:
										if len(data_father) < 3: 
											formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int( \
													  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int( \
													  allele_coverage_mother.split(',')[1])) + "\t" + '-' + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
										
										elif data_father [3] == '.':
											formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int( \
													  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int( \
													  allele_coverage_mother.split(',')[1])) + "\t" + data_father[3] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
										else: 
											formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int( \
													  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int( \
													  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int( \
													  allele_coverage_father.split(',')[1])) + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] 
										break
							
							else:
								if len(data_proband2) < 3:
									if len(data_mother) < 3:
										if len(data_father) < 3:
											formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int( \
													  allele_coverage_proband.split(',')[1])) + "\t" + '-' + "\t" + '-' + "\t" + '-' + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
											break
										else:
											if data_father[3] == '.':
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + '-' + "\t" + '-' + "\t" + \
														  data_father[3] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
											else:
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
														  '-' + "\t" + '-' + "\t" + \
														  str(int(allele_coverage_father.split(',')[0]) + int( \
														  allele_coverage_father.split(',')[1])) + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
									else:
										if data_mother[3] == '.':
											if len(data_father) < 3:
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
														  '-' + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
												break
											
											elif data_father[3] == '.':
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
													  	  '-' + "\t" + data_mother[3] + "\t" + data_father[3] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
											else: 
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
														  '-' + "\t" + data_mother[3] + "\t" + \
														  str(int(allele_coverage_father.split(',')[0]) + int( \
														  allele_coverage_father.split(',')[1])) + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
										else:
											if len(data_father) < 3: 
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
														  '-' + "\t" + str(int(allele_coverage_mother.split(',')[0]) + int( \
														  allele_coverage_mother.split(',')[1])) + "\t" + '-' + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
											
											elif data_father[3] == '.':
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
														  '-' + "\t" + str(int(allele_coverage_mother.split(',')[0]) + int( \
														  allele_coverage_mother.split(',')[1])) + "\t" + data_father[3] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
											else: 
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
														  '-' + "\t" + str(int(allele_coverage_mother.split(',')[0]) + int( \
														  allele_coverage_mother.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_father.split(',')[0]) + int( \
														  allele_coverage_father.split(',')[1])) + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
								elif data_proband2[3] == '.':
									if len(data_mother) < 3:
										if len(data_father) < 3:
											formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int( \
													  allele_coverage_proband.split(',')[1])) + "\t" + data_proband2[3] + \
													  "\t" + '-' + "\t" + '-' + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
											break
										else:
											if data_father[3] == '.':
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + data_proband2[3] + \
														  "\t" + '-' + "\t" + data_father[3] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
											else:
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
														  data_proband2[3] + "\t" + '-' + "\t" + \
														  str(int(allele_coverage_father.split(',')[0]) + int( \
														  allele_coverage_father.split(',')[1])) + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
									else:
										if data_mother[3] == '.':
											if len(data_father) < 3:
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
														  data_proband2[3] + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
												break
											
											elif data_father[3] == '.':
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
													  	  data_proband2[3] + "\t" + data_mother[3] + "\t" + \
														  data_father[3] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
											else: 
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
														  data_proband2[3] + "\t" + data_mother[3] + "\t" + \
														  str(int(allele_coverage_father.split(',')[0]) + int( \
														  allele_coverage_father.split(',')[1])) + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
										else:
											if len(data_father) < 3: 
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
														  data_proband2[3] + "\t" + str(int(allele_coverage_mother.split(',')[0]) + int( \
														  allele_coverage_mother.split(',')[1])) + "\t" + '-' + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
											
											elif data_father[3] == '.':
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
														  data_proband2[3] + "\t" + str(int(allele_coverage_mother.split(',')[0]) + int( \
														  allele_coverage_mother.split(',')[1])) + "\t" + data_father[3] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
											else: 
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
														  data_proband2[3] + "\t" + str(int(allele_coverage_mother.split(',')[0]) + int( \
														  allele_coverage_mother.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_father.split(',')[0]) + int( \
														  allele_coverage_father.split(',')[1])) + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break

								else:
									if len(data_mother) < 3:
										if len(data_father) < 3:
											formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int( \
													  allele_coverage_proband.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int( \
													  allele_coverage_proband2.split(',')[1])) + \
													  "\t" + '-' + "\t" + '-' + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
													  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
											break
										else:
											if data_father[3] == '.':
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_proband2.split(',')[0]) + int( \
													  	  allele_coverage_proband2.split(',')[1])) + \
														  "\t" + '-' + "\t" + data_father[3] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
											else:
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_proband2.split(',')[0]) + int( \
														  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + \
														  str(int(allele_coverage_father.split(',')[0]) + int( \
														  allele_coverage_father.split(',')[1])) + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
									else:
										if data_mother[3] == '.':
											if len(data_father) < 3:
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_proband2.split(',')[0]) + int( \
														  allele_coverage_proband2.split(',')[1])) + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
												break
											
											elif data_father[3] == '.':
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
													  	  str(int(allele_coverage_proband2.split(',')[0]) + int( \
														  allele_coverage_proband2.split(',')[1])) + "\t" + data_mother[3] + "\t" + \
														  data_father[3] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
											else: 
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_proband2.split(',')[0]) + int( \
														  allele_coverage_proband2.split(',')[1])) + "\t" + data_mother[3] + "\t" + \
														  str(int(allele_coverage_father.split(',')[0]) + int( \
														  allele_coverage_father.split(',')[1])) + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
										else:
											if len(data_father) < 3: 
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_proband2.split(',')[0]) + int( \
														  allele_coverage_proband2.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_mother.split(',')[0]) + int( \
														  allele_coverage_mother.split(',')[1])) + "\t" + '-' + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
											
											elif data_father[3] == '.':
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_proband2.split(',')[0]) + int( \
														  allele_coverage_proband2.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_mother.split(',')[0]) + int( \
														  allele_coverage_mother.split(',')[1])) + "\t" + data_father[3] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1]
											else: 
												formato = fields_cuarteto[0] + ":" + fields_cuarteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int( \
														  allele_coverage_proband.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_proband2.split(',')[0]) + int( \
														  allele_coverage_proband2.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_mother.split(',')[0]) + int( \
														  allele_coverage_mother.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_father.split(',')[0]) + int( \
														  allele_coverage_father.split(',')[1])) + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_mother.split(',')[1]+ "\t" + \
														  fields_cuarteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_cuarteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
								



					for h in lista:
						if genotype_proband[0] == "0":
							formato = formato +  "\t" + fields_cuarteto[3] + "/"
							break
						elif genotype_proband[0] == str(h + 1):
							formato = formato +  "\t" + ref_split[h] + "/"
							break
						elif genotype_proband[0] == ".":
							formato = formato +  "\t" + genotype_proband[0] + "/"
							break
					for h in lista:
						if genotype_proband[1] == "0":
							formato = formato + fields_cuarteto[3] + "\t"
							break
						elif genotype_proband[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_proband[1] == ".":
							formato = formato + genotype_proband[1] + "\t"
							break
					
					for h in lista:
						if genotype_proband2[0] == "0":
							formato = formato + fields_cuarteto[3] + "/"
							break
						elif genotype_proband2[0] == str(h + 1):
							formato = formato + ref_split[h] + "/"
							break
						elif genotype_proband2[0] == ".":
							formato = formato + genotype_proband2[0] + "/"
							break
					for h in lista:
						if genotype_proband2[1] == "0":
							formato = formato + fields_cuarteto[3] + "\t"
							break
						elif genotype_proband2[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_proband2[1] == ".":
							formato = formato + genotype_proband2[1] + "\t"
							break
					
					for h in lista:

						if genotype_mother[0] == "0":
							formato = formato + fields_cuarteto[3] + "/"
							break
						elif genotype_mother[0] == str(h + 1):
							formato = formato + ref_split[h] + "/"
							break
						elif genotype_mother[0] == ".":
							formato = formato + genotype_mother[0] + "/"
							break					
					for h in lista:

						if genotype_mother[1] == "0":
							formato = formato + fields_cuarteto[3] + "\t"
							break
						elif genotype_mother[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_mother[1] == ".":
							formato = formato + genotype_mother[1] + "\t"
							break
					for h in lista:
						if genotype_father[0] == "0":
							formato = formato + fields_cuarteto[3] + "/"
							break
						elif genotype_father[0] == str(h + 1):
							formato = formato + ref_split[h] + "/"
							break
						elif genotype_father[0] == ".":
							formato = formato + genotype_father[0] + "/"
							break
					for h in lista:
						if genotype_father[1] == "0":
							formato = formato + fields_cuarteto[3] + "\t"
							break
						elif genotype_father[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_father[1] == ".":
							formato = formato + genotype_father[1] + "\t"
							break
							
					#Añadimos la zigosidad
					
					for h in lista:
						if genotype_proband[0] == "0" and genotype_proband[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_proband[0] == "0" and genotype_proband[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t" 
							break
						elif genotype_proband[0] == str(h + 1) and genotype_proband[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_proband[0] == "." or genotype_proband[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
							
					for h in lista:
						if genotype_proband2[0] == "0" and genotype_proband2[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_proband2[0] == "0" and genotype_proband2[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t" 
							break
						elif genotype_proband2[0] == str(h + 1) and genotype_proband2[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_proband2[0] == "." or genotype_proband2[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
							
					for h in lista:
						if genotype_mother[0] == "0" and genotype_mother[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_mother[0] == "0" and genotype_mother[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t" 
							break
						elif genotype_mother[0] == str(h + 1) and genotype_mother[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_mother[0] == "." or genotype_mother[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
							
					for h in lista:
						if genotype_father[0] == "0" and genotype_father[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_father[0] == "0" and genotype_father[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t"
							break
						elif genotype_father[0] == str(h + 1) and genotype_father[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_father[0] == "." or genotype_father[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
						
					#Añadimos el VAF
                    
					for h in lista:
						if allele_coverage_proband.split(',') != "0" and '.' in allele_coverage_proband:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_proband.split(',') != "0" and '-' in allele_coverage_proband:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_proband.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_proband.split(',')[1])/ (float(allele_coverage_proband.split(',')[0]) + float(allele_coverage_proband.split(',')[1]))) + "\t" 
							break
						elif allele_coverage_proband.split(',')[0] == "0" and allele_coverage_proband.split(',')[1] == "0":
							formato = formato + "0.0" + "\t"
							break
						else:
							formato = formato + "1.0" + "\t"
							break
					for h in lista:
						if allele_coverage_proband2.split(',') != "0" and '.' in allele_coverage_proband2:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_proband2.split(',') != "0" and '-' in allele_coverage_proband2:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_proband2.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_proband2.split(',')[1])/ (float(allele_coverage_proband2.split(',')[0]) + float(allele_coverage_proband2.split(',')[1]))) + "\t" 
							break
						elif allele_coverage_proband2.split(',')[0] == "0" and allele_coverage_proband2.split(',')[1] == "0":
							formato = formato + "0.0" + "\t"
							break
						else:
							formato = formato + "1.0" + "\t"
							break
							
					for h in lista:
						if allele_coverage_mother.split(',') != "0" and '.' in allele_coverage_mother:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_mother.split(',') != "0" and '-' in allele_coverage_mother:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_mother.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_mother.split(',')[1])/ (float(allele_coverage_mother.split(',')[0]) + float(allele_coverage_mother.split(',')[1]))) + "\t" 
							break
						elif allele_coverage_mother.split(',')[0] == "0" and allele_coverage_mother.split(',')[1] == "0":
							formato = formato + "0.0" + "\t"
							break
						else:
							formato = formato + "1.0" + "\t"
							break
					for h in lista:
						if allele_coverage_father.split(',') != "0" and '.' in allele_coverage_father:
							formato = formato + ".,." 
							break
						elif allele_coverage_father.split(',') != "0" and '-' in allele_coverage_father:
							formato = formato + ".,." 
							break
						elif allele_coverage_father.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_father.split(',')[1])/ (float(allele_coverage_father.split(',')[0]) + float(allele_coverage_father.split(',')[1])))  
							break
						elif allele_coverage_father.split(',')[0] == "0" and allele_coverage_father.split(',')[1] == "0":
							formato = formato + "0.0" 
							break
						else:
							formato = formato + "1.0" 
							break
							
					#Añadimos origen

					if data_proband[0] == data_mother[0] and data_proband[0] == data_father[0] and data_proband[0] == data_proband2[0]:
						formato = formato + "\t" + "All"
					elif data_mother[0] == data_proband[0] and data_mother[0] == data_proband2[0]:
						formato = formato + "\t" + "Mother\Proband1\Proband2"
					elif data_father[0] == data_proband[0] and data_father[0] == data_proband2[0]:
						formato = formato + "\t" + "Father\Proband1\Proband2"
					elif data_father[0] == data_proband[0] and data_mother[0] == data_proband[0]:
						formato = formato + "\t" + "Proband1\Mother\Father"
					elif data_father[0] == data_proband2[0] and data_mother[0] == data_proband2[0]:
						formato = formato + "\t" + "Proband2\Mother\Father"
					elif data_father[0] == data_proband[0]:
						formato = formato + "\t" + "Proband1\Father"
					elif data_mother[0] == data_proband[0]:
						formato = formato + "\t" + "Proband1\Mother"
					elif data_father[0] == data_proband2[0]:
						formato = formato + "\t" + "Proband2\Father"
					elif data_mother[0] == data_proband2[0]:
						formato = formato + "\t" + "Proband2\Mother"
					elif data_proband[0] == data_proband2[0]:
						formato = formato + "\t" + "Proband1\Proband2"
					else:
						formato = formato + "\t" + "Check"

					cuarteto.append(formato)

		header_latest = "Locus" + "\t" +  str(caso) + " Coverage" + "\t" + str(caso2) + " Coverage" + "\t" + "Mother Depth" +"\t" + "Father Depth" + "\t" + str(caso) + " AlleleCoverage" + "\t" + str(caso2) + " AlleleCoverage" + "\t" + "Mother AlleleCoverage" + "\t" + "Father AlleleCoverage" + "\t" + str(caso) + " Genotype" + "\t" + str(caso2) + " Genotype" + "\t" "Mother Genotype" + "\t" + "Father Genotype" + "\t" + "zigosity " + str(caso) + "\t" + "zigosity " + str(caso2) + "\t" "zigosity Mother" + "\t" + "zigosity Father" + "\t" + "VAF " +  str(caso) + "\t" + "VAF " + str(caso2) + "\t" + "VAF Mother" + "\t" + "VAF Father" + "\t" + "Origin" + "\n"


		f_intermedio.write(header_latest)

		# Meto los datos de genotipos, profundidades y origen de las variantes en un archivo intermedio

		for line in cuarteto:
			# ~ print line
			f_intermedio.write("%s\n" % line)
		return True
class quintetoDragen:

	wd=""
	quintetos=""
	_filter=""
	
	def __init__(self,wd,quintetos,_filter): 

		self.wd=wd # los argumentos los define como . de self para poder referirse a ellos
		self.quintetos=quintetos
		self._filter=_filter


	

	def run(this,proband,father,mother,parient,parient_2,sexo): 

		if this.Filtrado(proband,father,mother,parient,parient_2,sexo):
			this.Columnas(proband)
			return this.Switch(proband)



	def Switch(this,caso):
	 	#Reordenamiento de columnas

		output_dir=get_output_dir(this.wd,this._filter)

		try:
	 		f_final_file= (output_dir+"/" + str(caso) + "_final.tsv")
	 		f_final = pd.read_csv(f_final_file, sep="\t", header=0)
		except (OSError, IOError) as e:
	 		raise ErrorFinalFile(f_final_file)



		f_finalisimo = (output_dir+"/" + str(caso) + "_finalisimo.tsv")
	
		columns_titles=['annonimous_GENE', 'GEN_OMIM', 'PHENO_OMIM', 'PHENO_DESC', 'INHERITANCE', 'clinvar_id', 'clinvar_CLNSIG', 'clinvar_CLNSIGCONF', 'clinvar_CLNDN', 'annonimous_PREDICTORS', 'CADD_phred', 'occurrence_parents', 'occurrence', 'Locus', 'REF', 'genotype', 'Mother Genotype', 'Father Genotype', 'Parient Genotype', 'Parient_2 Genotype', 'Origin', 'zigosity', 'zigosity Mother', 'zigosity Father', 'zigosity Parient', 'zigosity Parient_2', 'gnomad_AF', 'ALL_sites_2015_08', 'esp6500siv2_all', 'VAF', 'VAF Mother', 'VAF Father', 'VAF Parient', 'VAF Parient_2', 'Coverage', 'Mother Depth', 'Father Depth', 'Parient Depth', 'Parient_2 Depth','AlleleRatio', 'AlleleCoverage', 'Mother AlleleCoverage', 'Father AlleleCoverage', 'Parient AlleleCoverage', 'Parient_2 AlleleCoverage', 'annonimous_ANNOTATION', 'annonimous_Func_refGene', 'annonimous_gNomen', 'annonimous_tCPUNTO', 'reliabilities', 'annonimous_PPUNTO', 'annonimous_EXONS', 'annonimous_INTRONS', 'annonimous_TOTAL', 'distNearestSS', 'nearestSSType', 'dbscSNV_ADA_SCORE', 'dbscSNV_RF_SCORE', 'catalogue_XREF', 'catalogue_VP', 'catalogue_speciality', 'ExClinicoNIM_panel', 'DI', 'SFARI', 'Epilepsia', 'RetNet', 'incidental_findings', 'pseudogenes', 'imprinting', 'PAR_loci', 'dbSNP_id', 'cosmicIds', 'cosmicTissues', 'cosmicSampleCounts', 'ParentsNIMIDs', 'IlluminaProbandsNIMIDs', 'SIFT_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 'LRT_pred', 'MutationTaster_pred', 'MutationAssessor_pred', 'FATHMM_pred', 'MetaSVM_pred', 'gnomad_AF_afr', 'gnomad_AF_amr', 'gnomad_AF_asj', 'gnomad_AF_eas', 'gnomad_AF_sas', 'gnomad_AF_nfe', 'gnomad_AF_fin', 'gnomad_AF_oth', 'gnomad_HetCount_all', 'gnomad_HomCount_all', 'annonimous_HemCount_all', 'gnomad_Filter', 'oe_mis', 'oe_syn', 'pLI', 'oe_lof', 'oe_syn_lower', 'oe_syn_upper', 'oe_mis_lower', 'oe_mis_upper', 'oe_lof_lower', 'oe_lof_upper', 'syn_z', 'mis_z', 'lof_z', 'AFR_sites_2015_08', 'SAS_sites_2015_08', 'EAS_sites_2015_08', 'EUR_sites_2015_08', 'AMR_sites_2015_08', 'esp6500siv2_ea', 'esp6500siv2_aa', 'phyloP100way_vertebrate', 'phyloP20way_mammalian', 'phyloP46way_placental', 'phastCons100way_vertebrate', 'phastCons20way_mammalian', 'APPROVED_NAME', 'SSF_wt', 'SSF_mut', 'MES_wt', 'MES_mut', 'delta_MES', 'delta_SSF', 'SPiCEprobability', 'SPiCEinter_2thr']
	
		f_final = f_final.reindex(columns=columns_titles)
		f_final = f_final.reset_index(drop=True)
		f_final.to_csv(f_finalisimo, sep='\t') 




	def Columnas(this,caso):
		# Abrimos archivos

		path_latest=get_path_latest(caso,this.wd,this._filter)
		my_path_tmp=get_path_tmp(this.wd,this._filter)
		output_dir=get_output_dir(this.wd,this._filter)

		try:
			f = open(path_latest)

		except (OSError, IOError) as e:
			raise ErrorInputFileLatest(path_latest)

		try:
			tmp_file=str(my_path_tmp) + "/quinteto_intermedio_" + str(caso) + ".tsv"
			f_intermedio = open(tmp_file)
		except (OSError, IOError) as e:
			raise ErrorTempFile(tmp_file)

		f_final = open(output_dir+"/" + str(caso) + "_final.tsv", "w")

		locusIndex = LocusIndex(f)

		f.close()

		f = open(path_latest)

		# AddColummns de intermedio a final (latest + intermedio)

		cabecera = f.readline().split("\n")[0] #Lee la cabecera y se queda con la primera linea
		cabecera_final = cabecera.split("\t") #Parte la cabecera de latest por los tabuladores
		cabecera2 = f_intermedio.readline().split("\t") #Parte la cabecera del intermedio por los tabuladores
		# print (cabecera_final)
		# print (cabecera2)
		# cabecera.remove('\n')
		cabecera2.pop(0) #De la cabecera2 quita la primera posicion
		# print (cabecera)
		# print (cabecera2)

		rows = f.readlines() #Lee el latest tsv linea por linea
		rows2 = f_intermedio.readlines() #Lee el intermedio linea por linea
		rows3 = []

		f.close()
		f_intermedio.close()

		formato = ""
		formato2 = ""

		for i in range(len(cabecera_final)):
			formato2 = formato2 + cabecera_final[i] + "\t" #En formato2 guarda la cabecera final separada por tabuladores
		for j in range(len(cabecera2)):
			if j == len(cabecera2) - 1:
				formato2 = formato2 + cabecera2[j]
			else:
				formato2 = formato2 + cabecera2[j] + "\t"

		rows3.append(formato2) #Fusiona las filas

		# Añadimos IGV  a las lineas en las que no se han encontrado datos 

		for row in rows:
			prueba = row.split("\t\n") #Parte las lineas por tabulador y salto de línea y lo guarda en la variable prueba
			intermedio = prueba[0].split('\n') #Parte el primer argumento de prueba por salto de linea
			fields = intermedio[0].split('\t') #Parte el primer argumento de intermedio por los tabuladores
			Boo2 = True
			formato = ""
			for row2 in rows2:
				campos = row2.split("\t")
				if campos[0] == fields[locusIndex]:
					campos.pop(0)
					for i in range(len(fields)):
						formato = formato + fields[i] + "\t"
					for j in range(len(campos)):
						if j == len(campos) - 1:
							formato = formato + campos[j]
						else:
							formato = formato + campos[j] + "\t"
					rows3.append(formato)
					Boo2 = False
					break
				else:
					continue
			while Boo2 == True:
				for i in range(len(fields)):
					formato = formato + fields[i] + "\t"
				formato = formato + "IGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\n"
				rows3.append(formato)
				Boo2 = False

		for line in rows3:
			f_final.write(line)

		f_final.close()


	def Filtrado(this,caso,padre,madre,parient,parient2,sexo):
		# Abrimos archivos

		# Busca el latest en samples

		path_latest=get_path_latest(caso,this.wd,this._filter)
		quinteto_file=get_path_quinteto_file(caso,this.quintetos)
		my_path_tmp=get_path_tmp(this.wd,this._filter)

		try:
			f = open(path_latest) #Guardamos en f erl archivo de la ruta latest
		except (OSError, IOError) as e:
			raise ErrorInputFileLatest(path_latest)

		# Busca el quinteto en análisis

		try:
			f_quinteto = open(quinteto_file) #Guardamos en f trio el archivo trio

		except (OSError, IOError) as e:
			raise ErrorInputFileQuinteto(quinteto_file)

		f_intermedio = open(str(my_path_tmp) + "/quinteto_intermedio_" + str(caso) + ".tsv", "w+") #Guardamos en f_intermedio el archivo intermedio

		locusIndex = LocusIndex(f) #Con la función LocusIndex genero el index del Locus

		locus = []
		quinteto = []
		lines = f.readlines() #Guardo las lineas del archivo latest en la variable lines

		f.close()

		for line in lines:
			fields = line.split('\t') #Guardo en la variable fields los campos de las lineas del archivo latest
			locus.append(fields[locusIndex]) #Se añade el index del locus a los campos del archivo latest

		i = 0

		find_header=False;
		header_proband=None

		while(not find_header): #Indica como header la primera linea del vcf que comienza con CHROM
			header_proband = f_quinteto.readline()
			if header_proband.startswith("#CHROM"):
				find_header=True

		try:
			probandIndex,fatherIndex,motherIndex,parientIndex,parient2Index=getIndexFamily_5(caso,padre,madre,parient,parient2,header_proband) #Guarda el index de cada miembro de la familia
		except NotFoundFamily as e:
			this.__log.error(e)
			return False

		lines_quinteto = f_quinteto.readlines() 

		f_quinteto.close()

		for line in lines_quinteto:
			fields_quinteto = line.split("\t") 
			data_proband = fields_quinteto[probandIndex].split('\n')[0].split(":") #Guarda en la variable data_proband los campos del quinteto del probando partidos por salto de línea y :
			data_mother = fields_quinteto[motherIndex].split('\n')[0].split(":") #Lo mismo con la madre
			data_father = fields_quinteto[fatherIndex].split('\n')[0].split(":") #Lo mismo con el padre
			data_parient = fields_quinteto[parientIndex].split('\n')[0].split(":") #Lo mismo con el pariente
			data_parient2 = fields_quinteto[parient2Index].split('\n')[0].split(":") #Lo mismo con el segundo pariente
			
			
			#Prepara los datos para genotipar a los miembros de la familia
			
			if str(data_mother[0]) == '1' or str(data_mother[0]) == '0' or str(data_mother[0]) == '2':
				data_mother[0] = str(data_mother[0]) + '/' + str(data_mother[0])
			
			if str(data_father[0]) == '1' or str(data_father[0]) == '0' or str(data_father[0]) == '2':
				data_father[0] = str(data_father[0]) + '/' + str(data_father[0])

			if str(data_proband[0]) == '1' or str(data_proband[0]) == '0' or str(data_proband[0]) == '2':
				data_proband[0] = str(data_proband[0]) + '/' + str(data_proband[0])
				
			if str(data_parient[0]) == '1' or str(data_parient[0]) == '0' or str(data_parient[0]) == '2':
				data_parient[0] = str(data_parient[0]) + '/' + str(data_parient[0])
				
			if str(data_parient2[0]) == '1' or str(data_parient2[0]) == '0' or str(data_parient2[0]) == '2':
				data_parient2[0] = str(data_parient2[0]) + '/' + str(data_parient2[0])


			if '/' in data_proband[0]:
				genotype_proband = data_proband[0].split("/")
				genotype_proband_0 = genotype_proband[0]
				genotype_proband_1 = genotype_proband[1]
			elif '|' in data_proband[0]:
				genotype_proband = data_proband[0].split("|")
				genotype_proband_0 = genotype_proband[0]
				genotype_proband_1 = genotype_proband[1]

			if '/' in data_mother[0]:
				genotype_mother = data_mother[0].split("/")
				genotype_mother_0 = genotype_mother[0]
				genotype_mother_1 = genotype_mother[1]
			elif '|' in data_mother[0]:
				genotype_mother = data_mother[0].split("|")
				genotype_mother_0 = genotype_mother[0]
				genotype_mother_1 = genotype_mother[1]

			if '/' in data_father[0]:
				genotype_father = data_father[0].split("/")
				genotype_father_0 = genotype_father[0]
				genotype_father_1 = genotype_father[1]
			elif '|' in data_father[0]:
				genotype_father = data_father[0].split("|")
				genotype_father_0 = genotype_father[0]
				genotype_father_1 = genotype_father[1]
				
			if '/' in data_parient[0]:
				genotype_parient = data_parient[0].split("/")
				genotype_parient_0 = genotype_parient[0]
				genotype_parient_1 = genotype_parient[1]
			elif '|' in data_parient[0]:
				genotype_parient = data_parient[0].split("|")
				genotype_parient_0 = genotype_parient[0]
				genotype_parient_1 = genotype_parient[1]
				
			if '/' in data_parient2[0]:
				genotype_parient2 = data_parient2[0].split("/")
				genotype_parient2_0 = genotype_parient2[0]
				genotype_parient2_1 = genotype_parient2[1]
			elif '|' in data_parient2[0]:
				genotype_parient2 = data_parient2[0].split("|")
				genotype_parient2_0 = genotype_parient2[0]
				genotype_parient2_1 = genotype_parient2[1]

			#Guarda la cobertura de los miembros de la familia si la tienen
			if len(data_mother) < 3:
				allele_coverage_mother = '-,-'
			elif len(data_mother[1].split(',')) == 1 and genotype_mother != './.':
				allele_coverage_mother = data_mother[1] + ',0'
			else:
				allele_coverage_mother = data_mother[1]

			if len(data_father) < 3:
				allele_coverage_father = '-,-'
			elif len(data_father[1].split(',')) == 1 and genotype_father != './.':
				allele_coverage_father = data_father[1] + ',0'
			else:
				allele_coverage_father = data_father[1]
				
			if len(data_parient) < 3:
				allele_coverage_parient = '-,-'
			elif len(data_parient[1].split(',')) == 1 and genotype_parient != './.':
				allele_coverage_parient = data_parient[1] + ',0'
			else:
				allele_coverage_parient = data_parient[1]
				
			if len(data_parient2) < 3:
				allele_coverage_parient2 = '-,-'
			elif len(data_parient2[1].split(',')) == 1 and genotype_parient2 != './.':
				allele_coverage_parient2 = data_parient2[1] + ',0'
			else:
				allele_coverage_parient2 = data_parient2[1]

			ref_split = fields_quinteto[4].split(",")
			lista = range(len(ref_split))

			# HACERLO INDEXANDO VALORES Y QUE BUSQUE EN REF CON EL MISMO INDEX

			for field in locus:
				fields_locus = field.split(":")
				if fields_quinteto[0] == fields_locus[0] and fields_quinteto[1] == fields_locus[1]:
					for h in lista:
						#print field #TODAS ESTAS RAMAS DE IF AÑADEN LA INFO DE COBERTURA Y GENOTIPO
						if len(data_mother) < 3:
							if len(data_father) < 3:
								if len(data_parient) <3:
									if len(data_parient2) <3:
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + '-' + "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1] 
										break
									else:
										if data_parient2[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + '-' + "\t" + '-' + "\t" + \
													  data_parient2[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + '-' + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_parient2.split(',')[0]) + int(
														  allele_coverage_parient2.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
										break
								else:
									if data_parient[3] == '.':
										if len(data_parient2) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + '-' + "\t" + \
													  data_parient[3] + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
										
										elif data_parient2 [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + '-' + "\t" + \
													  data_parient[3] + "\t" + data_parient2[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + '-' + "\t" + \
													  data_parient[3] + "\t" + \
													  str(int(allele_coverage_parient2.split(',')[0]) + int(
														  allele_coverage_parient2.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
										break
									else:
										if len(data_parient2) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
														  allele_coverage_parient.split(',')[1])) + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
										
										elif data_parient2 [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
														  allele_coverage_parient.split(',')[1])) + "\t" + data_parient2[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
														  allele_coverage_parient.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_parient2.split(',')[0]) + int(
														  allele_coverage_parient2.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1] 
										break
										
							else:
								if data_father[3] == '.':
									if len(data_parient) < 3:
										if len(data_parient2) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + data_father[3] + "\t" + \
													  '-' + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
										
										elif data_parient2 [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + data_father[3] + "\t" + \
													  '-' + "\t" + data_parient2[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + data_father[3] + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_parient2.split(',')[0]) + int(
														  allele_coverage_parient2.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
										break
									else:
										if data_parient[3] == '.':
											if len(data_parient2) < 3: 
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + data_father[3] + "\t" + \
														  data_parient[3] + "\t" + '-' + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
											
											elif data_parient2 [3] == '.':
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + data_father[3] + "\t" + \
														  data_parient[3] + "\t" + data_parient2[3] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
											else: 
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + data_father[3] + "\t" + \
														  data_parient[3] + "\t" + \
														  str(int(allele_coverage_parient2.split(',')[0]) + int(
															  allele_coverage_parient2.split(',')[1])) + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
											break
										else:
											if len(data_parient2) < 3: 
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + data_father[3] + "\t" + \
														  str(int(allele_coverage_parient.split(',')[0]) + int(
															  allele_coverage_parient.split(',')[1])) + "\t" + '-' + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
											
											elif data_parient2 [3] == '.':
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + data_father[3] + "\t" + \
														  str(int(allele_coverage_parient.split(',')[0]) + int(
															  allele_coverage_parient.split(',')[1])) + "\t" + data_parient2[3] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
											else: 
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + data_father[3] + "\t" + \
														  str(int(allele_coverage_parient.split(',')[0]) + int(
															  allele_coverage_parient.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_parient2.split(',')[0]) + int(
															  allele_coverage_parient2.split(',')[1])) + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1] 
											break
								else:
									if len(data_parient) < 3:
										if len(data_parient2) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
														  allele_coverage_father.split(',')[1])) + \
													  "\t" + '-' + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
											
										elif data_parient2 [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
														  allele_coverage_father.split(',')[1])) + \
													  "\t" + '-' + "\t" + data_parient2[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
														  allele_coverage_father.split(',')[1])) + \
													  "\t" + '-' + "\t" + \
													  str(int(allele_coverage_parient2.split(',')[0]) + int(
														  allele_coverage_parient2.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
										break
									else:
										if data_parient[3] == '.':
											if len(data_parient2) < 3: 
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
														  allele_coverage_father.split(',')[1])) + \
													  "\t" + data_parient[3] + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
											
											elif data_parient2 [3] == '.':
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
														  allele_coverage_father.split(',')[1])) + \
													  "\t" + data_parient[3] + "\t" + data_parient2[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
											else: 
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
														  allele_coverage_father.split(',')[1])) + \
													  "\t" + data_parient[3] + "\t" + str(int(allele_coverage_parient2.split(',')[0]) + \
													   int(allele_coverage_parient2.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
											break
										else:
											if len(data_parient2) < 3: 
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
														  allele_coverage_father.split(',')[1])) + "\t" +\
													  str(int(allele_coverage_parient.split(',')[0]) + int(
														  allele_coverage_parient.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
											
											elif data_parient2 [3] == '.':
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
														  allele_coverage_father.split(',')[1])) + "\t" +\
													  str(int(allele_coverage_parient.split(',')[0]) + int(
														  allele_coverage_parient.split(',')[1])) + "\t" + \
													  data_parient2[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
											else: 
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
														  allele_coverage_father.split(',')[1])) + "\t" +\
													  str(int(allele_coverage_parient.split(',')[0]) + int(
														  allele_coverage_parient.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_parient2.split(',')[0]) + int(
														  allele_coverage_parient2.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1] 
											break
			
						elif data_mother[3] == '.':
							if len(data_father) < 3:
								if len(data_parient) <3:
									if len(data_parient2) <3:
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + '-' + "\t" + '-' + \
												  "\t" + '-' + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1] 
										break
									elif data_parient2[3] == '.':
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + '-' + "\t" + '-' + "\t" + \
												  data_parient2[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									else:
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + '-' + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_parient2.split(',')[0]) + int(
													  allele_coverage_parient2.split(',')[1])) + "\t" +\
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									break
								elif data_parient[3] == '.':
									if len(data_parient2) < 3: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
												  data_parient[3] + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
										
									elif data_parient2 [3] == '.':
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
												  data_parient[3] + "\t" + data_parient2[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									else: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
												  data_parient[3] + "\t" + \
												  str(int(allele_coverage_parient2.split(',')[0]) + int(
													  allele_coverage_parient2.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									break
								else:
									if len(data_parient2) < 3: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									
									elif data_parient2 [3] == '.':
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + "\t" + data_parient2[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									else: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_parient2.split(',')[0]) + int(
													  allele_coverage_parient2.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1] 
									break
										
							elif data_father[3] == '.':
								if len(data_parient) < 3:
									if len(data_parient2) < 3: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + \
												  data_father[3] + "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
										
									elif data_parient2 [3] == '.':
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + \
												  data_father[3] + "\t" + '-' + "\t" + data_parient2[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									else: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + \
												  data_father[3] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_parient2.split(',')[0]) + int(
													  allele_coverage_parient2.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									break
								elif data_parient[3] == '.':
									if len(data_parient2) < 3: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + \
												  data_father[3] + "\t" + data_parient[3] + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
											
									elif data_parient2 [3] == '.':
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + \
												  data_father[3] + "\t" + data_parient[3] + "\t" + data_parient2[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									else: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + \
												  data_father[3] + "\t" + data_parient[3] + "\t" + \
												  str(int(allele_coverage_parient2.split(',')[0]) + int(
													  allele_coverage_parient2.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									break
								else:
									if len(data_parient2) < 3: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + data_father[3] + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
											
									elif data_parient2 [3] == '.':
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + data_father[3] + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + "\t" + data_parient2[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									else: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + data_father[3] + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
												 	  allele_coverage_parient.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_parient2.split(',')[0]) + int(
													  allele_coverage_parient2.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1] 
									break
							else:
								if len(data_parient) < 3:
									if len(data_parient2) < 3: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + \
												  "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
											
									elif data_parient2 [3] == '.':
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + \
												  "\t" + '-' + "\t" + data_parient2[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									else: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + \
												  "\t" + '-' + "\t" + \
												  str(int(allele_coverage_parient2.split(',')[0]) + int(
													  allele_coverage_parient2.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									break
								elif data_parient[3] == '.':
									if len(data_parient2) < 3: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + \
												  "\t" + data_parient[3] + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
										
									elif data_parient2 [3] == '.':
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + \
												  "\t" + data_parient[3] + "\t" + data_parient2[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									else: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + \
												  "\t" + data_parient[3] + "\t" + str(int(allele_coverage_parient2.split(',')[0]) + \
												   int(allele_coverage_parient2.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									break
								else:
									if len(data_parient2) < 3: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + "\t" +\
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
										
									elif data_parient2 [3] == '.':
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + "\t" +\
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + "\t" + \
												  data_parient2[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									else: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" +\
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_parient2.split(',')[0]) + int(
													  allele_coverage_parient2.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1] 
									break
							
						else:
							if len(data_father) < 3:
								if len(data_parient) <3:
									if len(data_parient2) <3:
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  '-' + "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1] 
										break
									elif data_parient2[3] == '.':
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  '-' + "\t" + '-' + "\t" + data_parient2[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									else:
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  '-' + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_parient2.split(',')[0]) + int(
													  allele_coverage_parient2.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									break
								elif data_parient[3] == '.':
									if len(data_parient2) < 3: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  '-' + "\t" + data_parient[3] + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
										
									elif data_parient2 [3] == '.':
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  '-' + "\t" + data_parient[3] + "\t" + data_parient2[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									else: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  '-' + "\t" + data_parient[3] + "\t" + \
												  str(int(allele_coverage_parient2.split(',')[0]) + int(
													  allele_coverage_parient2.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									break
								else:
									if len(data_parient2) < 3: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
											  		  allele_coverage_mother.split(',')[1])) + "\t" + \
											  	  '-' + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									
									elif data_parient2 [3] == '.':
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + "\t" + data_parient2[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									else: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_parient2.split(',')[0]) + int(
													  allele_coverage_parient2.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1] 
									break
										
							elif data_father[3] == '.':
								if len(data_parient) < 3:
									if len(data_parient2) < 3: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  data_father[3] + "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
										
									elif data_parient2 [3] == '.':
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  data_father[3] + "\t" + '-' + "\t" + data_parient2[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									else: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  data_father[3] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_parient2.split(',')[0]) + int(
													  allele_coverage_parient2.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									break
								elif data_parient[3] == '.':
									if len(data_parient2) < 3: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  data_father[3] + "\t" + data_parient[3] + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
											
									elif data_parient2 [3] == '.':
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  data_father[3] + "\t" + data_parient[3] + "\t" + data_parient2[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									else: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  data_father[3] + "\t" + data_parient[3] + "\t" + \
												  str(int(allele_coverage_parient2.split(',')[0]) + int(
													  allele_coverage_parient2.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									break
								else:
									if len(data_parient2) < 3: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  data_father[3] + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
											
									elif data_parient2 [3] == '.':
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  data_father[3] + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + "\t" + data_parient2[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									else: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  data_father[3] + "\t" + \
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_parient2.split(',')[0]) + int(
													  allele_coverage_parient2.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1] 
									break
							else:
								if len(data_parient) < 3:
									if len(data_parient2) < 3: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + \
												  "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
											
									elif data_parient2 [3] == '.':
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + \
												  "\t" + '-' + "\t" + data_parient2[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									else: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + \
												  "\t" + '-' + "\t" + \
												  str(int(allele_coverage_parient2.split(',')[0]) + int(
													  allele_coverage_parient2.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									break
								elif data_parient[3] == '.':
									if len(data_parient2) < 3: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + \
												  "\t" + data_parient[3] + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
										
									elif data_parient2 [3] == '.':
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + \
												  "\t" + data_parient[3] + "\t" + data_parient2[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									else: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + \
												  "\t" + data_parient[3] + "\t" + str(int(allele_coverage_parient2.split(',')[0]) + \
												   int(allele_coverage_parient2.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									break
								else:
									if len(data_parient2) < 3: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + "\t" +\
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
										
									elif data_parient2 [3] == '.':
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
													  allele_coverage_father.split(',')[1])) + "\t" +\
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + "\t" + \
												  data_parient2[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1]
									else: 
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  allele_coverage_father.split(',')[1])) + "\t" +\
												  str(int(allele_coverage_parient.split(',')[0]) + int(
													  allele_coverage_parient.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_parient2.split(',')[0]) + int(
													  allele_coverage_parient2.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient2.split(',')[1] 
									break


					for h in lista:

						if genotype_mother[0] == "0":
							formato = formato + "\t" + fields_quinteto[3] + "/"
							break
						elif genotype_mother[0] == str(h + 1):
							formato = formato + "\t" + ref_split[h] + "/"
							break
						elif genotype_mother[0] == ".":
							formato = formato + "\t" + genotype_mother[0] + "/"
							break
					for h in lista:

						if genotype_mother[1] == "0":
							formato = formato + fields_quinteto[3] + "\t"
							break
						elif genotype_mother[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_mother[1] == ".":
							formato = formato + genotype_mother[1] + "\t"
							break
					for h in lista:
						if genotype_father[0] == "0":
							formato = formato + fields_quinteto[3] + "/"
							break
						elif genotype_father[0] == str(h + 1):
							formato = formato + ref_split[h] + "/"
							break
						elif genotype_father[0] == ".":
							formato = formato + genotype_father[0] + "/"
							break
					for h in lista:
						if genotype_father[1] == "0":
							formato = formato + fields_quinteto[3] + "\t"
							break
						elif genotype_father[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_father[1] == ".":
							formato = formato + genotype_father[1] + "\t"
							break		
					for h in lista:
						if genotype_parient[0] == "0":
							formato = formato + fields_quinteto[3] + "/"
							break
						elif genotype_parient[0] == str(h + 1):
							formato = formato + ref_split[h] + "/"
							break
						elif genotype_parient[0] == ".":
							formato = formato + genotype_parient[0] + "/"
							break
					for h in lista:
						if genotype_parient[1] == "0":
							formato = formato + fields_quinteto[3] + "\t"
							break
						elif genotype_parient[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_parient[1] == ".":
							formato = formato + genotype_parient[1] + "\t"
							break
					for h in lista:
						if genotype_parient2[0] == "0":
							formato = formato + fields_quinteto[3] + "/"
							break
						elif genotype_parient2[0] == str(h + 1):
							formato = formato + ref_split[h] + "/"
							break
						elif genotype_parient2[0] == ".":
							formato = formato + genotype_parient2[0] + "/"
							break
					for h in lista:
						if genotype_parient2[1] == "0":
							formato = formato + fields_quinteto[3] + "\t"
							break
						elif genotype_parient2[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_parient2[1] == ".":
							formato = formato + genotype_parient2[1] + "\t"
							break
							
					#Añadimos la zigosidad
					
					for h in lista:
						if genotype_mother[0] == "0" and genotype_mother[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_mother[0] == "0" and genotype_mother[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t" 
							break
						elif genotype_mother[0] == str(h + 1) and genotype_mother[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_mother[0] == "." or genotype_mother[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
							
					for h in lista:
						if genotype_father[0] == "0" and genotype_father[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_father[0] == "0" and genotype_father[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t" 
							break
						elif genotype_father[0] == str(h + 1) and genotype_father[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_father[0] == "." or genotype_father[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
							
					for h in lista:
						if genotype_parient[0] == "0" and genotype_parient[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_parient[0] == "0" and genotype_parient[1] == str(h + 1):
							formato = formato + "Heteroz"  + "\t"
							break
						elif genotype_parient[0] == str(h + 1) and genotype_parient[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_parient[0] == "." or genotype_parient[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
							
					for h in lista:
						if genotype_parient2[0] == "0" and genotype_parient2[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_parient2[0] == "0" and genotype_parient2[1] == str(h + 1):
							formato = formato + "Heteroz"  + "\t"
							break
						elif genotype_parient2[0] == str(h + 1) and genotype_parient2[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_parient2[0] == "." or genotype_parient2[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
							
					#Añadimos el VAF
					
					for h in lista:
						if allele_coverage_mother.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_mother.split(',')[1])// float(allele_coverage_mother.split(',')[0]) + float(allele_coverage_mother.split(',')[1])) + "\t" 
							break
						elif allele_coverage_mother.split(',')[0] == "0" and allele_coverage_mother.split(',')[1] == "0":
							formato = formato + "0.0" + "\t"
							break
						else:
							formato = formato + "1.0" + "\t"
							break
					for h in lista:
						if allele_coverage_father.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_father.split(',')[1])// float(allele_coverage_father.split(',')[0]) + float(allele_coverage_father.split(',')[1])) + "\t"
							break
						elif allele_coverage_father.split(',')[0] == "0" and allele_coverage_father.split(',')[1] == "0":
							formato = formato + "0.0" + "\t"
							break
						else:
							formato = formato +"1.0"+ "\t"
							break
					for h in lista:
						if allele_coverage_parient.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_parient.split(',')[1])// float(allele_coverage_parient.split(',')[0]) + float(allele_coverage_parient.split(',')[1])) + "\t"
							break
						elif allele_coverage_parient.split(',')[0] == "0" and allele_coverage_parient.split(',')[1] == "0":
							formato = formato + "0.0" + "\t"
							break
						else:
							formato = formato + "1.0" + "\t"
							break
					for h in lista:
						if allele_coverage_parient2.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_parient2.split(',')[1])// float(allele_coverage_parient2.split(',')[0]) + float(allele_coverage_parient2.split(',')[1]))
							break
						elif allele_coverage_parient2.split(',')[0] == "0" and allele_coverage_parient2.split(',')[1] == "0":
							formato = formato + "0.0"
							break
						else:
							formato = formato + "1.0"
							break
							
					#Añadimos el origen

					if data_proband[0] == data_mother[0] and data_proband[0] == data_father[0] and data_proband[0] == data_parient[0] and data_proband[0] == data_parient2[0]:
						formato = formato + "\t" + "All"
					elif data_proband[0] == data_mother[0] and data_proband[0] == data_father[0] and data_proband[0] == data_parient[0]:
						formato = formato + "\t" + "Mother/Father/Parient"
					elif data_proband[0] == data_mother[0] and data_proband[0] == data_father[0] and data_proband[0] == data_parient2[0]:
						formato = formato + "\t" + "Mother/Father/Parient2"
					elif data_proband[0] == data_mother[0] and data_proband[0] == data_parient[0] and data_proband[0] == data_parient2[0]:
						formato = formato + "\t" + "Mother/Parient/Parient2"
					elif data_proband[0] == data_father[0] and data_proband[0] == data_parient[0] and data_proband[0] == data_parient2[0]:
						formato = formato + "\t" + "Father/Parient/Parient2"
					elif data_proband[0] == data_mother[0] and data_proband[0] == data_father[0]:
						formato = formato + "\t" + "Mother/Father"
					elif data_proband[0] == data_mother[0] and data_proband[0] == data_parient[0]:
						formato = formato + "\t" + "Mother/Parient"
					elif data_proband[0] == data_mother[0] and data_proband[0] == data_parient2[0]:
						formato = formato + "\t" + "Mother/Parient2"
					elif data_proband[0] == data_father[0] and data_proband[0] == data_parient[0]:
						formato = formato + "\t" + "Father/Parient"
					elif data_proband[0] == data_father[0] and data_proband[0] == data_parient2[0]:
						formato = formato + "\t" + "Father/Parient2"
					elif data_proband[0] == data_parient[0] and data_proband[0] == data_parient2[0]:
						formato = formato + "\t" + "Parient/Parient2"
					else:
						formato = formato + "\t" + "Check"

					quinteto.append(formato)

		header_latest = "Locus" + "\t" + "Mother Depth" + "\t" + "Father Depth" + "\t" + "Parient Depth" + "\t" + "Parient_2 Depth" + "\t" + "Mother AlleleCoverage" + "\t" + "Father AlleleCoverage" + "\t" + "Parient AlleleCoverage" + "\t" + "Parient_2 AlleleCoverage" + "\t" + "Mother Genotype" + "\t" + "Father Genotype" + "\t" + "Parient Genotype" +  "\t" + "Parient_2 Genotype" + "\t" + "zigosity Mother" + "\t" + "zigosity Father" + "\t" + "zigosity Parient" + "\t" + "zigosity Parient_2" + "\t" + "VAF Mother" + "\t" + "VAF Father" + "\t" + "VAF Parient" + "\t" + "VAF Parient_2" + "\t" + "Origin" + "\n"


		f_intermedio.write(header_latest)

		# Meto los datos de genotipos, profundidades y origen de las variantes en un archivo intermedio

		for line in quinteto:
			# ~ print line
			f_intermedio.write("%s\n" % line)
		return True
		
		
class quintetoDragen_2P:
	__log=None


	wd=""
	quintetos=""
	_filter=""
	
	def __init__(self,wd,quintetos,_filter): 

		self.wd=wd # los argumentos los define como . de self para poder referirse a ellos
		self.quintetos=quintetos
		self._filter=_filter


	

	def run(this,proband,proband2,father,mother,parient,sexo,sexo2):


		if this.Filtrado(proband,proband2,father,mother,parient,sexo,sexo2):
			this.Columnas(proband,proband2)
			return this.Switch(proband,proband2)



	def Switch(this,caso,caso2):
	 	#Reordenamiento de columnas

		output_dir=get_output_dir(this.wd,this._filter)

		try:
	 		f_final_file= (output_dir+"/" + str(caso) + '+' + str(caso2) + "_final.tsv")
	 		f_final = pd.read_csv(f_final_file, sep="\t", header=0)
		except (OSError, IOError) as e:
	 		raise ErrorFinalFile(f_final_file)



		f_finalisimo = (output_dir+"/" + str(caso) + '+' + str(caso2) + "_finalisimo.tsv")

		f_final.drop(['genotype', 'zigosity', 'VAF', 'Coverage', 'AlleleRatio'], axis=1, inplace=True)
	
		columns_titles=['annonimous_GENE', 'GEN_OMIM', 'PHENO_OMIM', 'PHENO_DESC', 'INHERITANCE', 'clinvar_id', 'clinvar_CLNSIG', 'clinvar_CLNSIGCONF', 'clinvar_CLNDN', 'annonimous_PREDICTORS', 'CADD_phred', 'occurrence_parents', 'occurrence', 'Locus', 'REF', str(caso) + ' Genotype', str(caso2) + ' Genotype', 'Mother Genotype', 'Father Genotype', 'Parient Genotype', 'Origin', 'zigosity ' + str(caso), 'zigosity ' + str(caso2), 'zigosity Mother', 'zigosity Father', 'zigosity Parient', 'gnomad_AF', 'ALL_sites_2015_08', 'esp6500siv2_all', 'VAF ' +  str(caso), 'VAF ' +  str(caso2), 'VAF Mother', 'VAF Father', 'VAF Parient', str(caso) + ' Coverage', str(caso2) + ' Coverage', 'Mother Depth', 'Father Depth', 'Parient Depth', str(caso) + ' AlleleCoverage', str(caso2) + ' AlleleCoverage', 'Mother AlleleCoverage', 'Father AlleleCoverage', 'Parient AlleleCoverage', 'annonimous_ANNOTATION', 'annonimous_Func_refGene', 'annonimous_gNomen', 'annonimous_tCPUNTO', 'reliabilities', 'annonimous_PPUNTO', 'annonimous_EXONS', 'annonimous_INTRONS', 'annonimous_TOTAL', 'distNearestSS', 'nearestSSType', 'dbscSNV_ADA_SCORE', 'dbscSNV_RF_SCORE', 'catalogue_XREF', 'catalogue_VP', 'catalogue_speciality', 'ExClinicoNIM_panel', 'DI', 'SFARI', 'Epilepsia', 'RetNet', 'incidental_findings', 'pseudogenes', 'imprinting', 'PAR_loci', 'dbSNP_id', 'cosmicIds', 'cosmicTissues', 'cosmicSampleCounts', 'ParentsNIMIDs', 'IlluminaProbandsNIMIDs', 'SIFT_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 'LRT_pred', 'MutationTaster_pred', 'MutationAssessor_pred', 'FATHMM_pred', 'MetaSVM_pred', 'gnomad_AF_afr', 'gnomad_AF_amr', 'gnomad_AF_asj', 'gnomad_AF_eas', 'gnomad_AF_sas', 'gnomad_AF_nfe', 'gnomad_AF_fin', 'gnomad_AF_oth', 'gnomad_HetCount_all', 'gnomad_HomCount_all', 'annonimous_HemCount_all', 'gnomad_Filter', 'oe_mis', 'oe_syn', 'pLI', 'oe_lof', 'oe_syn_lower', 'oe_syn_upper', 'oe_mis_lower', 'oe_mis_upper', 'oe_lof_lower', 'oe_lof_upper', 'syn_z', 'mis_z', 'lof_z', 'AFR_sites_2015_08', 'SAS_sites_2015_08', 'EAS_sites_2015_08', 'EUR_sites_2015_08', 'AMR_sites_2015_08', 'esp6500siv2_ea', 'esp6500siv2_aa', 'phyloP100way_vertebrate', 'phyloP20way_mammalian', 'phyloP46way_placental', 'phastCons100way_vertebrate', 'phastCons20way_mammalian', 'APPROVED_NAME', 'SSF_wt', 'SSF_mut', 'MES_wt', 'MES_mut', 'delta_MES', 'delta_SSF', 'SPiCEprobability', 'SPiCEinter_2thr']
	
		f_final = f_final.reindex(columns=columns_titles)
		f_final = f_final.reset_index(drop=True)
		f_final.to_csv(f_finalisimo, sep='\t')  



	def Columnas(this,caso,caso2):
		
		# Abrimos archivos

		path_latest=get_path_latest_2P(caso,caso2,this.wd,this._filter)
		my_path_tmp=get_path_tmp(this.wd,this._filter)
		output_dir=get_output_dir(this.wd,this._filter)

		try:
			f = open(path_latest)

		except (OSError, IOError) as e:
			raise ErrorInputFileLatest(path_latest)

		try:
			tmp_file=str(my_path_tmp) + "/quinteto_intermedio_" + str(caso) + "+" + str(caso2) + ".tsv"
			f_intermedio = open(tmp_file)
		except (OSError, IOError) as e:
			raise ErrorTempFile(tmp_file)

		f_final = open(output_dir+"/" + str(caso) + "+" + str(caso2) + "_final.tsv", "w")

		locusIndex = LocusIndex(f)

		f.close()

		f = open(path_latest)

		# AddColummns de intermedio a final (latest + intermedio)

		cabecera = f.readline().split("\n")[0] #Lee la cabecera y se queda con la primera linea
		cabecera_final = cabecera.split("\t") #Parte la cabecera de latest por los tabuladores
		cabecera2 = f_intermedio.readline().split("\t") #Parte la cabecera del intermedio por los tabuladores
		# print (cabecera_final)
		# print (cabecera2)
		# cabecera.remove('\n')
		cabecera2.pop(0) #De la cabecera2 quita la primera posicion
		# print (cabecera)
		# print (cabecera2)

		rows = f.readlines() #Lee el latest tsv linea por linea
		rows2 = f_intermedio.readlines() #Lee el intermedio linea por linea
		rows3 = []

		f.close()
		f_intermedio.close()

		formato = ""
		formato2 = ""

		for i in range(len(cabecera_final)):
			formato2 = formato2 + cabecera_final[i] + "\t" #En formato2 guarda la cabecera final separada por tabuladores
		for j in range(len(cabecera2)):
			if j == len(cabecera2) - 1:
				formato2 = formato2 + cabecera2[j]
			else:
				formato2 = formato2 + cabecera2[j] + "\t"

		rows3.append(formato2) #Fusiona las filas

		# Añadimos IGV  a las lineas en las que no se han encontrado datos

		for row in rows:
			prueba = row.split("\t\n") #Parte las lineas por tabulador y salto de línea y lo guarda en la variable prueba
			intermedio = prueba[0].split('\n') #Parte el primer argumento de prueba por salto de linea
			fields = intermedio[0].split('\t') #Parte el primer argumento de intermedio por los tabuladores
			Boo2 = True
			formato = ""
			for row2 in rows2:
				campos = row2.split("\t")
				if campos[0] == fields[locusIndex]:
					campos.pop(0)
					for i in range(len(fields)):
						formato = formato + fields[i] + "\t"
					for j in range(len(campos)):
						if j == len(campos) - 1:
							formato = formato + campos[j]
						else:
							formato = formato + campos[j] + "\t"
					rows3.append(formato)
					Boo2 = False
					break
				else:
					continue
			while Boo2 == True:
				for i in range(len(fields)):
					formato = formato + fields[i] + "\t"
				formato = formato + "IGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\n"
				rows3.append(formato)
				Boo2 = False

		for line in rows3:
			f_final.write(line)

		f_final.close()


	def Filtrado(this,caso,caso2,padre,madre,parient,sexo,sexo2):
		# Abrimos archivos

		#Generamos el archivo latest común
		
		path_latest_1=get_path_latest(caso,this.wd,this._filter)
		path_latest_2=get_path_latest(caso2,this.wd,this._filter)
		
		latest_1 = pd.read_csv(path_latest_1, sep="\t", header=0)
		latest_2 = pd.read_csv(path_latest_2, sep="\t", header=0)
		
		joint_latest = pd.concat([latest_1,latest_2], ignore_index=True).drop_duplicates(['Locus'])
		
		joint_latest[['chr','position']]=joint_latest.Locus.str.split(':',expand=True)
		joint_latest['chr']=joint_latest.chr.str.replace('chr','')
		joint_latest['chr']=joint_latest.chr.str.replace('X','30')
		joint_latest['chr']=joint_latest.chr.str.replace('Y','31')
		joint_latest['chr']=joint_latest.chr.str.replace('M','32')
		joint_latest['chr']=joint_latest['chr'].astype(int)
		joint_latest['position']=joint_latest['position'].astype(int)

		joint_latest_sort = joint_latest.sort_values(by=['chr','position'])
		joint_latest_sort = joint_latest_sort.drop(['chr','position'], axis=1)
		joint_latest_sort = joint_latest_sort.reset_index(drop=True)
		
		output_dir=get_output_dir(this.wd,this._filter)
		
		file_joint_latest = (output_dir+"/" + str(caso) + '+' + str(caso2) + "_latest.tsv")
		
		joint_latest_sort.to_csv(file_joint_latest, sep='\t')  

		# Busca el latest en samples

		path_latest=get_path_latest_2P(caso,caso2,this.wd,this._filter)
		quinteto_file=get_path_quinteto_file_2P(caso,caso2,this.quintetos)
		my_path_tmp=get_path_tmp(this.wd,this._filter)

		try:
			f = open(path_latest) #Guardamos en f el archivo de la ruta latest
		except (OSError, IOError) as e:
			raise ErrorInputFileLatest(path_latest)

		# Busca el quinteto en análisis

		try:
			f_quinteto = open(quinteto_file) #Guardamos en f quinteto el archivo quinteto

		except (OSError, IOError) as e:
			raise ErrorInputFileQuinteto(quinteto_file)

		f_intermedio = open(str(my_path_tmp) + "/quinteto_intermedio_" + str(caso) + "+" + str(caso2) + ".tsv", "w+") #Guardamos en f_intermedio el archivo intermedio

		locusIndex = LocusIndex(f) #Con la función LocusIndex genero el index del Locus

		locus = []
		quinteto = []
		lines = f.readlines() #Guardo las lineas del archivo latest en la variable lines

		f.close()

		for line in lines:
			fields = line.split('\t') #Guardo en la variable fields los campos de las lineas del archivo latest
			locus.append(fields[locusIndex]) #Se añade el index del locus a los campos del archivo latest

		i = 0

		find_header=False;
		header_proband=None

		while(not find_header): #Indica como header la primera linea del vcf que comienza con CHROM
			header_proband = f_quinteto.readline()
			if header_proband.startswith("#CHROM"):
				find_header=True

		try:
			probandIndex,proband2Index,fatherIndex,motherIndex,parientIndex=getIndexFamily_5_2P(caso,caso2,padre,madre,parient,header_proband) #Guarda el indice de cada miembro 
		except NotFoundFamily as e:
			this.__log.error(e)
			return False

		lines_quinteto = f_quinteto.readlines() 

		f_quinteto.close()

		for line in lines_quinteto:
			fields_quinteto = line.split("\t") 
			data_proband = fields_quinteto[probandIndex].split('\n')[0].split(":") #Guarda en la variable data_proband los campos del quinteto del probando partidos por salto de línea y :
			data_proband2 = fields_quinteto[proband2Index].split('\n')[0].split(":") #Lo mismo con el segundo probando2
			data_mother = fields_quinteto[motherIndex].split('\n')[0].split(":") #Lo mismo con la madre
			data_father = fields_quinteto[fatherIndex].split('\n')[0].split(":") #Lo mismo con el padre
			data_parient = fields_quinteto[parientIndex].split('\n')[0].split(":") #Lo mismo con el pariente
			
			
			#Prepara los datos para genotipar a los miembros de la familia
			
			if str(data_mother[0]) == '1' or str(data_mother[0]) == '0' or str(data_mother[0]) == '2':
				data_mother[0] = str(data_mother[0]) + '/' + str(data_mother[0])
			
			if str(data_father[0]) == '1' or str(data_father[0]) == '0' or str(data_father[0]) == '2':
				data_father[0] = str(data_father[0]) + '/' + str(data_father[0])

			if str(data_proband[0]) == '1' or str(data_proband[0]) == '0' or str(data_proband[0]) == '2':
				data_proband[0] = str(data_proband[0]) + '/' + str(data_proband[0])
				
			if str(data_parient[0]) == '1' or str(data_parient[0]) == '0' or str(data_parient[0]) == '2':
				data_parient[0] = str(data_parient[0]) + '/' + str(data_parient[0])
				
			if str(data_proband2[0]) == '1' or str(data_proband2[0]) == '0' or str(data_proband2[0]) == '2':
				data_proband2[0] = str(data_proband2[0]) + '/' + str(data_proband2[0])


			if '/' in data_proband[0]:
				genotype_proband = data_proband[0].split("/")
				genotype_proband_0 = genotype_proband[0]
				genotype_proband_1 = genotype_proband[1]
			elif '|' in data_proband[0]:
				genotype_proband = data_proband[0].split("|")
				genotype_proband_0 = genotype_proband[0]
				genotype_proband_1 = genotype_proband[1]
				
			if '/' in data_proband2[0]:
				genotype_proband2 = data_proband2[0].split("/")
				genotype_proband2_0 = genotype_proband2[0]
				genotype_proband2_1 = genotype_proband2[1]
			elif '|' in data_proband2[0]:
				genotype_proband2 = data_proband2[0].split("|")
				genotype_proband2_0 = genotype_proband2[0]
				genotype_proband2_1 = genotype_proband2[1]

			if '/' in data_mother[0]:
				genotype_mother = data_mother[0].split("/")
				genotype_mother_0 = genotype_mother[0]
				genotype_mother_1 = genotype_mother[1]
			elif '|' in data_mother[0]:
				genotype_mother = data_mother[0].split("|")
				genotype_mother_0 = genotype_mother[0]
				genotype_mother_1 = genotype_mother[1]

			if '/' in data_father[0]:
				genotype_father = data_father[0].split("/")
				genotype_father_0 = genotype_father[0]
				genotype_father_1 = genotype_father[1]
			elif '|' in data_father[0]:
				genotype_father = data_father[0].split("|")
				genotype_father_0 = genotype_father[0]
				genotype_father_1 = genotype_father[1]
				
			if '/' in data_parient[0]:
				genotype_parient = data_parient[0].split("/")
				genotype_parient_0 = genotype_parient[0]
				genotype_parient_1 = genotype_parient[1]
			elif '|' in data_parient[0]:
				genotype_parient = data_parient[0].split("|")
				genotype_parient_0 = genotype_parient[0]
				genotype_parient_1 = genotype_parient[1]


			#Guarda la cobertura de los miembros de la familia si la tienen
			if len(data_proband) < 3:
				allele_coverage_proband = '-,-'
			elif len(data_proband[1].split(',')) == 1 and genotype_proband != './.':
				allele_coverage_proband = data_proband[1] + ',0'
			else:
				allele_coverage_proband = data_proband[1]
			
			if len(data_proband2) < 3:
				allele_coverage_proband2 = '-,-'
			elif len(data_proband2[1].split(',')) == 1 and genotype_proband2 != './.':
				allele_coverage_proband2 = data_proband2[1] + ',0'
			else:
				allele_coverage_proband2 = data_proband2[1]
			
			if len(data_mother) < 3:
				allele_coverage_mother = '-,-'
			elif len(data_mother[1].split(',')) == 1 and genotype_mother != './.':
				allele_coverage_mother = data_mother[1] + ',0'
			else:
				allele_coverage_mother = data_mother[1]

			if len(data_father) < 3:
				allele_coverage_father = '-,-'
			elif len(data_father[1].split(',')) == 1 and genotype_father != './.':
				allele_coverage_father = data_father[1] + ',0'
			else:
				allele_coverage_father = data_father[1]
				
			if len(data_parient) < 3:
				allele_coverage_parient = '-,-'
			elif len(data_parient[1].split(',')) == 1 and genotype_parient != './.':
				allele_coverage_parient = data_parient[1] + ',0'
			else:
				allele_coverage_parient = data_parient[1]


			ref_split = fields_quinteto[4].split(",")
			lista = range(len(ref_split))

			# HACERLO INDEXANDO VALORES Y QUE BUSQUE EN REF CON EL MISMO INDEX

			for field in locus:
				fields_locus = field.split(":")
				if fields_quinteto[0] == fields_locus[0] and fields_quinteto[1] == fields_locus[1]:
					for h in lista:
						#print field #TODAS ESTAS RAMAS DE IF AÑADEN LA INFO DE COBERTURA Y GENOTIPO
						if len(data_proband) < 3:
							if len(data_mother) < 3:
								if len(data_father) <3:
									if len(data_parient) <3:
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
									else:
										if data_parient[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + '-' + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
														  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
								else:
									if data_father[3] == '.':
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + \
													  data_father[3] + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + \
													  data_father[3] + "\t" + data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + data_father[3] + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
														  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
									else:
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										break
										
							elif data_mother[3] == '.':
								if len(data_father) <3:
									if len(data_parient) <3:
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_mother[3] + "\t" + \
												  '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
									else:
										if data_parient[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
														  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
								else:
									if data_father[3] == '.':
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_mother[3] + "\t" + \
													  data_father[3] + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_mother[3] + "\t" + \
													  data_father[3] + "\t" + data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_mother[3] + "\t" + \
													  data_father[3] + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
														  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
									else:
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  data_mother[3] + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  data_mother[3] + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  data_mother[3] + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										break
							
								
							else:
								if len(data_father) <3:
									if len(data_parient) <3:
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
									else:
										if data_parient[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
								else:
									if data_father[3] == '.':
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  data_father[3] + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										
										elif data_parient [3] == '.':
												  formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  data_father[3] + "\t" + data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  data_father[3] + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
									else:
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										break
			
						elif data_proband[3] == '.':
							if len(data_mother) < 3:
								if len(data_father) <3:
									if len(data_parient) <3:
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
									else:
										if data_parient[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + '-' + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
														  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
								else:
									if data_father[3] == '.':
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + \
													  data_father[3] + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + \
													  data_father[3] + "\t" + data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + data_father[3] + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
														  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
									else:
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										break
										
							elif data_mother[3] == '.':
								if len(data_father) <3:
									if len(data_parient) <3:
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_mother[3] + "\t" + \
												  '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
									else:
										if data_parient[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
														  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
								else:
									if data_father[3] == '.':
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_mother[3] + "\t" + \
													  data_father[3] + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_mother[3] + "\t" + \
													  data_father[3] + "\t" + data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_mother[3] + "\t" + \
													  data_father[3] + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
														  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
									else:
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  data_mother[3] + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  data_mother[3] + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  data_mother[3] + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										break
							
								
							else:
								if len(data_father) <3:
									if len(data_parient) <3:
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
									else:
										if data_parient[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
								else:
									if data_father[3] == '.':
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  data_father[3] + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  data_father[3] + "\t" + data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  data_father[3] + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
									else:
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										break
							
						else:
							if len(data_proband2) < 3:
								if len(data_mother) < 3:
									if len(data_father) < 3:
										if len(data_parient) < 3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  '-' + "\t" + \
													  '-' + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
											break
										elif data_parient[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  '-' + "\t" + \
													  '-' + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  '-' + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
									elif data_father[3] == '.':
										if len(data_parient) < 3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  '-' + "\t" + \
													  data_father[3] + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
											
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  '-' + "\t" + \
													  data_father[3] + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  '-' + "\t" + \
													  data_father[3] + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
									else:
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]  
										break
											
								elif data_mother[3] == '.':
									if len(data_father) < 3:
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  data_mother[3] + "\t" + \
													  '-' + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
											
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  data_mother[3] + "\t" + \
													  '-' + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  data_mother[3] + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										break
									elif data_father[3] == '.':
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  data_mother[3] + "\t" + \
													  data_father[3] + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
												
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  data_mother[3] + "\t" + \
													  data_father[3] + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  data_mother[3] + "\t" + \
													  data_father[3] + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										break
									else:
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  data_mother[3] + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
												
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  data_mother[3] + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  data_mother[3] + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										break
								else:
									if len(data_father) < 3:
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
												
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										break
									elif data_father[3] == '.':
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  data_father[3] + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
											
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  data_father[3] + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  data_father[3] + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										break
									else:
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
											
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
										
							elif data_proband2[3] == '.':
								if len(data_mother) < 3:
									if len(data_father) < 3:
										if len(data_parient) < 3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  '-' + "\t" + \
													  '-' + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
											break
										elif data_parient[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  '-' + "\t" + \
													  '-' + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  '-' + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
									elif data_father[3] == '.':
										if len(data_parient) < 3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  '-' + "\t" + \
													  data_father[3] + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
											
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  '-' + "\t" + \
													  data_father[3] + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  '-' + "\t" + \
													  data_father[3] + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
									else:
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]  
										break
											
								elif data_mother[3] == '.':
									if len(data_father) < 3:
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  data_mother[3] + "\t" + \
													  '-' + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
											
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  data_mother[3] + "\t" + \
													  '-' + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  data_mother[3] + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										break
									elif data_father[3] == '.':
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  data_mother[3] + "\t" + \
													  data_father[3] + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
												
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  data_mother[3] + "\t" + \
													  data_father[3] + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  data_mother[3] + "\t" + \
													  data_father[3] + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										break
									else:
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  data_mother[3] + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
												
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  data_mother[3] + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  data_mother[3] + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										break
								else:
									if len(data_father) < 3:
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
												
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										break
									elif data_father[3] == '.':
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  data_father[3] + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
											
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  data_father[3] + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  data_father[3] + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										break
									else:
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
											
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  data_proband2[3] + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break

							else:
								if len(data_mother) < 3:
									if len(data_father) < 3:
										if len(data_parient) < 3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  '-' + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
											break
										elif data_parient[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  '-' + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
									elif data_father[3] == '.':
										if len(data_parient) < 3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  data_father[3] + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
											
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  data_father[3] + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  data_father[3] + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break
									else:
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]  
										break
											
								elif data_mother[3] == '.':
									if len(data_father) < 3:
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  data_mother[3] + "\t" + \
													  '-' + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
											
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  data_mother[3] + "\t" + \
													  '-' + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  data_mother[3] + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										break
									elif data_father[3] == '.':
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  data_mother[3] + "\t" + \
													  data_father[3] + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
												
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  data_mother[3] + "\t" + \
													  data_father[3] + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  data_mother[3] + "\t" + \
													  data_father[3] + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										break
									else:
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  data_mother[3] + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
												
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  data_mother[3] + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  data_mother[3] + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										break
								else:
									if len(data_father) < 3:
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
												
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										break
									elif data_father[3] == '.':
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  data_father[3] + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
											
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  data_father[3] + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  data_father[3] + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										break
									else:
										if len(data_parient) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
											
										elif data_parient [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  data_parient[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												 	  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_parient.split(',')[0]) + int(
													  	  allele_coverage_parient.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_parient.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_parient.split(',')[1] 
										break



					for h in lista:

						if genotype_proband[0] == "0":
							formato = formato + "\t" + fields_quinteto[3] + "/"
							break
						elif genotype_proband[0] == str(h + 1):
							formato = formato + "\t" + ref_split[h] + "/"
							break
						elif genotype_proband[0] == ".":
							formato = formato + "\t" + genotype_proband[0] + "/"
							break
					for h in lista:

						if genotype_proband[1] == "0":
							formato = formato + fields_quinteto[3] + "\t"
							break
						elif genotype_proband[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_proband[1] == ".":
							formato = formato + genotype_proband[1] + "\t"
							break
					for h in lista:

						if genotype_proband2[0] == "0":
							formato = formato + fields_quinteto[3] + "/"
							break
						elif genotype_proband2[0] == str(h + 1):
							formato = formato + ref_split[h] + "/"
							break
						elif genotype_proband2[0] == ".":
							formato = formato + genotype_proband2[0] + "/"
							break
					for h in lista:

						if genotype_proband2[1] == "0":
							formato = formato + fields_quinteto[3] + "\t"
							break
						elif genotype_proband2[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_proband2[1] == ".":
							formato = formato + genotype_proband2[1] + "\t"
							break
					for h in lista:

						if genotype_mother[0] == "0":
							formato = formato + fields_quinteto[3] + "/"
							break
						elif genotype_mother[0] == str(h + 1):
							formato = formato + ref_split[h] + "/"
							break
						elif genotype_mother[0] == ".":
							formato = formato + genotype_mother[0] + "/"
							break
					for h in lista:

						if genotype_mother[1] == "0":
							formato = formato + fields_quinteto[3] + "\t"
							break
						elif genotype_mother[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_mother[1] == ".":
							formato = formato + genotype_mother[1] + "\t"
							break
					for h in lista:
						if genotype_father[0] == "0":
							formato = formato + fields_quinteto[3] + "/"
							break
						elif genotype_father[0] == str(h + 1):
							formato = formato + ref_split[h] + "/"
							break
						elif genotype_father[0] == ".":
							formato = formato + genotype_father[0] + "/"
							break
					for h in lista:
						if genotype_father[1] == "0":
							formato = formato + fields_quinteto[3] + "\t"
							break
						elif genotype_father[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_father[1] == ".":
							formato = formato + genotype_father[1] + "\t"
							break		
					for h in lista:
						if genotype_parient[0] == "0":
							formato = formato + fields_quinteto[3] + "/"
							break
						elif genotype_parient[0] == str(h + 1):
							formato = formato + ref_split[h] + "/"
							break
						elif genotype_parient[0] == ".":
							formato = formato + genotype_parient[0] + "/"
							break
					for h in lista:
						if genotype_parient[1] == "0":
							formato = formato + fields_quinteto[3] + "\t"
							break
						elif genotype_parient[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_parient[1] == ".":
							formato = formato + genotype_parient[1] + "\t"
							break
							
					#Añadimos la zigosidad
					
					for h in lista:
						if genotype_proband[0] == "0" and genotype_proband[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_proband[0] == "0" and genotype_proband[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t" 
							break
						elif genotype_proband[0] == str(h + 1) and genotype_proband[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_proband[0] == "." or genotype_proband[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
							
					for h in lista:
						if genotype_proband2[0] == "0" and genotype_proband2[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_proband2[0] == "0" and genotype_proband2[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t" 
							break
						elif genotype_proband2[0] == str(h + 1) and genotype_proband2[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_proband2[0] == "." or genotype_proband2[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
							
					for h in lista:
						if genotype_mother[0] == "0" and genotype_mother[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_mother[0] == "0" and genotype_mother[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t" 
							break
						elif genotype_mother[0] == str(h + 1) and genotype_mother[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_mother[0] == "." or genotype_mother[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
							
					for h in lista:
						if genotype_father[0] == "0" and genotype_father[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_father[0] == "0" and genotype_father[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t"
							break
						elif genotype_father[0] == str(h + 1) and genotype_father[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_father[0] == "." or genotype_father[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
					
					for h in lista:
						if genotype_parient[0] == "0" and genotype_parient[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_parient[0] == "0" and genotype_parient[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t"
							break
						elif genotype_parient[0] == str(h + 1) and genotype_parient[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_parient[0] == "." or genotype_parient[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
						
					#Añadimos el VAF
					
					for h in lista:
						if allele_coverage_proband.split(',') != "0" and '.' in allele_coverage_proband:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_proband.split(',') != "0" and '-' in allele_coverage_proband:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_proband.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_proband.split(',')[1])/ (float(allele_coverage_proband.split(',')[0]) + float(allele_coverage_proband.split(',')[1]))) + "\t" 
							break
						elif allele_coverage_proband.split(',')[0] == "0" and allele_coverage_proband.split(',')[1] == "0":
							formato = formato + "0.0" + "\t"
							break
						else:
							formato = formato + "1.0" + "\t"
							break
					for h in lista:
						if allele_coverage_proband2.split(',') != "0" and '.' in allele_coverage_proband2:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_proband2.split(',') != "0" and '-' in allele_coverage_proband2:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_proband2.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_proband2.split(',')[1])/ (float(allele_coverage_proband2.split(',')[0]) + float(allele_coverage_proband2.split(',')[1]))) + "\t" 
							break
						elif allele_coverage_proband2.split(',')[0] == "0" and allele_coverage_proband2.split(',')[1] == "0":
							formato = formato + "0.0" + "\t"
							break
						else:
							formato = formato + "1.0" + "\t"
							break
							
					for h in lista:
						if allele_coverage_mother.split(',') != "0" and '.' in allele_coverage_mother:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_mother.split(',') != "0" and '-' in allele_coverage_mother:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_mother.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_mother.split(',')[1])/ (float(allele_coverage_mother.split(',')[0]) + float(allele_coverage_mother.split(',')[1]))) + "\t" 
							break
						elif allele_coverage_mother.split(',')[0] == "0" and allele_coverage_mother.split(',')[1] == "0":
							formato = formato + "0.0" + "\t"
							break
						else:
							formato = formato + "1.0" + "\t"
							break
					for h in lista:
						if allele_coverage_father.split(',') != "0" and '.' in allele_coverage_father:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_father.split(',') != "0" and '-' in allele_coverage_father:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_father.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_father.split(',')[1])/ (float(allele_coverage_father.split(',')[0]) + float(allele_coverage_father.split(',')[1]))) + "\t" 
							break
						elif allele_coverage_father.split(',')[0] == "0" and allele_coverage_father.split(',')[1] == "0":
							formato = formato + "0.0" + "\t"
							break
						else:
							formato = formato + "1.0" + "\t"
							break
					for h in lista:
						if allele_coverage_parient.split(',') != "0" and '.' in allele_coverage_parient:
							formato = formato + ".,." 
							break
						elif allele_coverage_parient.split(',') != "0" and '-' in allele_coverage_parient:
							formato = formato + ".,." 
							break
						elif allele_coverage_parient.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_parient.split(',')[1])/ (float(allele_coverage_parient.split(',')[0]) + float(allele_coverage_parient.split(',')[1]))) 
							break
						elif allele_coverage_parient.split(',')[0] == "0" and allele_coverage_parient.split(',')[1] == "0":
							formato = formato + "0.0" 
							break
						else:
							formato = formato + "1.0" 
							break
							
					#Añadimos el origen

					if data_proband[0] == data_mother[0] and data_proband[0] == data_father[0] and data_proband[0] == data_parient[0] and data_proband[0] == data_proband2[0]:
						formato = formato + "\t" + "All"
					elif data_proband[0] == data_mother[0] and data_proband[0] == data_father[0] and data_proband[0] == data_parient[0]:
						formato = formato + "\t" + "Proband/Mother/Father/Parient"
					elif data_proband2[0] == data_mother[0] and data_proband2[0] == data_father[0] and data_proband2[0] == data_parient[0]:
						formato = formato + "\t" + "Proband2/Mother/Father/Parient"
					elif data_proband[0] == data_mother[0] and data_proband[0] == data_father[0] and data_proband[0] == data_proband2[0]:
						formato = formato + "\t" + "Proband/Proband2/Mother/Father"
					elif data_proband[0] == data_mother[0] and data_proband[0] == data_parient[0] and data_proband[0] == data_proband2[0]:
						formato = formato + "\t" + "Proband/Proband2/Mother/Parient"
					elif data_proband[0] == data_father[0] and data_proband[0] == data_parient[0] and data_proband[0] == data_proband2[0]:
						formato = formato + "\t" + "Proband/Proband2/Father/Parient"
					elif data_proband[0] == data_mother[0] and data_proband2[0] == data_mother[0]:
						formato = formato + "\t" + "Proband/Proband2/Mother"
					elif data_proband[0] == data_father[0] and data_proband2[0] == data_father[0]:
						formato = formato + "\t" + "Proband/Proband2/Father"
					elif data_proband[0] == data_parient[0] and data_proband2[0] == data_parient[0]:
						formato = formato + "\t" + "Proband/Proband2/Parient"
					elif data_proband[0] == data_mother[0] and data_proband[0] == data_father[0]:
						formato = formato + "\t" + "Proband/Mother/Father"
					elif data_proband2[0] == data_mother[0] and data_proband2[0] == data_father[0]:
						formato = formato + "\t" + "Proband2/Mother/Father"
					elif data_proband[0] == data_mother[0] and data_proband[0] == data_parient[0]:
						formato = formato + "\t" + "Proband/Mother/Parient"
					elif data_proband2[0] == data_mother[0] and data_proband2[0] == data_parient[0]:
						formato = formato + "\t" + "Proband2/Mother/Parient"
					elif data_proband[0] == data_father[0] and data_proband[0] == data_parient[0]:
						formato = formato + "\t" + "Proband/Father/Parient"
					elif data_proband2[0] == data_father[0] and data_proband2[0] == data_parient[0]:
						formato = formato + "\t" + "Proband2/Father/Parient"
					elif data_proband[0] == data_mother[0]:
						formato = formato + "\t" + "Proband/Mother"
					elif data_proband2[0] == data_mother[0]:
						formato = formato + "\t" + "Proband2/Mother"
					elif data_proband[0] == data_father[0]:
						formato = formato + "\t" + "Proband/Father"
					elif data_proband2[0] == data_father[0]:
						formato = formato + "\t" + "Proband2/Father"
					elif data_proband[0] == data_parient[0]:
						formato = formato + "\t" + "Proband/Parient"
					elif data_proband2[0] == data_parient[0]:
						formato = formato + "\t" + "Proband2/Parient"
					elif data_proband[0] == data_proband2[0]:
						formato = formato + "\t" + "Proband/Proband2"
					else:
						formato = formato + "\t" + "Check"

					quinteto.append(formato)

		header_latest = "Locus" + "\t" +  str(caso) + " Coverage" + "\t" + str(caso2) + " Coverage" + "\t" + "Mother Depth" +"\t" + "Father Depth" + "\t" + "Parient Depth" + "\t" + str(caso) + " AlleleCoverage" + "\t" + str(caso2) + " AlleleCoverage" + "\t" + "Mother AlleleCoverage" + "\t" + "Father AlleleCoverage" + "\t" + "Parient AlleleCoverage" + "\t" + str(caso) + " Genotype" + "\t" + str(caso2) + " Genotype" + "\t" "Mother Genotype" + "\t" + "Father Genotype" + "\t" + "Parient Genotype" + "\t" + "zigosity " + str(caso) + "\t" + "zigosity " + str(caso2) + "\t" "zigosity Mother" + "\t" + "zigosity Father" + "\t" +"zigosity Parient" + "\t" + "VAF " +  str(caso) + "\t" + "VAF " + str(caso2) + "\t" + "VAF Mother" + "\t" + "VAF Father" + "\t" + "VAF Parient" + "\t" + "Origin" + "\n"


		f_intermedio.write(header_latest)

		# Meto los datos de genotipos, profundidades y origen de las variantes en un archivo intermedio

		for line in quinteto:
			# ~ print line
			f_intermedio.write("%s\n" % line)
		return True

class quintetoDragen_3P:
	__log=None


	wd=""
	quintetos=""
	_filter=""
	
	def __init__(self,wd,quintetos,_filter): #

		self.wd=wd # los argumentos los define como . de self para poder referirse a ellos
		self.quintetos=quintetos
		self._filter=_filter


	

	def run(this,proband,proband2,proband3,father,mother,sexo, sexo2, sexo3): #Obtiene informacion del log


		if this.Filtrado(proband,proband2,proband3,father,mother,sexo,sexo2,sexo3):
			this.Columnas(proband,proband2,proband3)
			return this.Switch(proband, proband2, proband3)



	def Switch(this,caso,caso2,caso3):
	 	#Reordenamiento de columnas

		output_dir=get_output_dir(this.wd,this._filter)

		try:
	 		f_final_file= (output_dir+"/" + str(caso) + "+" + str(caso2) + "+" + str(caso3) + "_final.tsv")
	 		f_final = pd.read_csv(f_final_file, sep="\t", header=0)
		except (OSError, IOError) as e:
	 		raise ErrorFinalFile(f_final_file)



		f_finalisimo = (output_dir+"/" + str(caso) + "+" + str(caso2) + "+" + str(caso3) + "_finalisimo.tsv")

		f_final.drop(['genotype', 'zigosity', 'VAF', 'Coverage', 'AlleleRatio'], axis=1, inplace=True)
	
		columns_titles=['annonimous_GENE', 'GEN_OMIM', 'PHENO_OMIM', 'PHENO_DESC', 'INHERITANCE', 'clinvar_id', 'clinvar_CLNSIG', 'clinvar_CLNSIGCONF', 'clinvar_CLNDN', 'annonimous_PREDICTORS', 'CADD_phred', 'occurrence_parents', 'occurrence', 'Locus', 'REF', str(caso) + ' Genotype', str(caso2) + ' Genotype', str(caso3) + ' Genotype', 'Mother Genotype', 'Father Genotype', 'Origin', 'zigosity ' + str(caso), 'zigosity ' + str(caso2), 'zigosity ' + str(caso3), 'zigosity Mother', 'zigosity Father', 'gnomad_AF', 'ALL_sites_2015_08', 'esp6500siv2_all', 'VAF ' +  str(caso), 'VAF ' +  str(caso2), 'VAF ' +  str(caso3), 'VAF Mother', 'VAF Father', str(caso) + ' Coverage', str(caso2) + ' Coverage', str(caso3) + ' Coverage', 'Mother Depth', 'Father Depth', str(caso) + ' AlleleCoverage', str(caso2) + ' AlleleCoverage', str(caso3) + ' AlleleCoverage', 'Mother AlleleCoverage', 'Father AlleleCoverage', 'annonimous_ANNOTATION', 'annonimous_Func_refGene', 'annonimous_gNomen', 'annonimous_tCPUNTO', 'reliabilities', 'annonimous_PPUNTO', 'annonimous_EXONS', 'annonimous_INTRONS', 'annonimous_TOTAL', 'distNearestSS', 'nearestSSType', 'dbscSNV_ADA_SCORE', 'dbscSNV_RF_SCORE', 'catalogue_XREF', 'catalogue_VP', 'catalogue_speciality', 'ExClinicoNIM_panel', 'DI', 'SFARI', 'Epilepsia', 'RetNet', 'incidental_findings', 'pseudogenes', 'imprinting', 'PAR_loci', 'dbSNP_id', 'cosmicIds', 'cosmicTissues', 'cosmicSampleCounts', 'ParentsNIMIDs', 'IlluminaProbandsNIMIDs', 'SIFT_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 'LRT_pred', 'MutationTaster_pred', 'MutationAssessor_pred', 'FATHMM_pred', 'MetaSVM_pred', 'gnomad_AF_afr', 'gnomad_AF_amr', 'gnomad_AF_asj', 'gnomad_AF_eas', 'gnomad_AF_sas', 'gnomad_AF_nfe', 'gnomad_AF_fin', 'gnomad_AF_oth', 'gnomad_HetCount_all', 'gnomad_HomCount_all', 'annonimous_HemCount_all', 'gnomad_Filter', 'oe_mis', 'oe_syn', 'pLI', 'oe_lof', 'oe_syn_lower', 'oe_syn_upper', 'oe_mis_lower', 'oe_mis_upper', 'oe_lof_lower', 'oe_lof_upper', 'syn_z', 'mis_z', 'lof_z', 'AFR_sites_2015_08', 'SAS_sites_2015_08', 'EAS_sites_2015_08', 'EUR_sites_2015_08', 'AMR_sites_2015_08', 'esp6500siv2_ea', 'esp6500siv2_aa', 'phyloP100way_vertebrate', 'phyloP20way_mammalian', 'phyloP46way_placental', 'phastCons100way_vertebrate', 'phastCons20way_mammalian', 'APPROVED_NAME', 'SSF_wt', 'SSF_mut', 'MES_wt', 'MES_mut', 'delta_MES', 'delta_SSF', 'SPiCEprobability', 'SPiCEinter_2thr']
	
		f_final = f_final.reindex(columns=columns_titles)
		f_final = f_final.reset_index(drop=True)
		f_final.to_csv(f_finalisimo, sep='\t')  

	def Columnas(this,caso,caso2,caso3):
		
		# Abrimos archivos

		path_latest=get_path_latest_3P(caso,caso2,caso3,this.wd,this._filter)
		my_path_tmp=get_path_tmp(this.wd,this._filter)
		output_dir=get_output_dir(this.wd,this._filter)

		try:
			f = open(path_latest)

		except (OSError, IOError) as e:
			raise ErrorInputFileLatest(path_latest)

		try:
			tmp_file=str(my_path_tmp) + "/quinteto_intermedio_" + str(caso) + "+" + str(caso2) + "+" + str(caso3) + ".tsv"
			f_intermedio = open(tmp_file)
		except (OSError, IOError) as e:
			raise ErrorTempFile(tmp_file)

		f_final = open(output_dir+"/" + str(caso) + "+" + str(caso2) + "+" + str(caso3) + "_final.tsv", "w")

		locusIndex = LocusIndex(f)

		f.close()

		f = open(path_latest)

		# AddColummns de intermedio a final (latest + intermedio)

		cabecera = f.readline().split("\n")[0] #Lee la cabecera y se queda con la primera linea
		cabecera_final = cabecera.split("\t") #Parte la cabecera de latest por los tabuladores
		cabecera2 = f_intermedio.readline().split("\t") #Parte la cabecera del intermedio por los tabuladores
		# print (cabecera_final)
		# print (cabecera2)
		# cabecera.remove('\n')
		cabecera2.pop(0) #De la cabecera2 quita la primera posicion
		# print (cabecera)
		# print (cabecera2)

		rows = f.readlines() #Lee el latest tsv linea por linea
		rows2 = f_intermedio.readlines() #Lee el intermedio linea por linea
		rows3 = []

		f.close()
		f_intermedio.close()

		formato = ""
		formato2 = ""

		for i in range(len(cabecera_final)):
			formato2 = formato2 + cabecera_final[i] + "\t" #En formato2 guarda la cabecera final separada por tabuladores
		for j in range(len(cabecera2)):
			if j == len(cabecera2) - 1:
				formato2 = formato2 + cabecera2[j]
			else:
				formato2 = formato2 + cabecera2[j] + "\t"

		rows3.append(formato2) #Fusiona las filas

		# Añadimos IGV  a las lineas en las que no se han encontrado datos

		for row in rows:
			prueba = row.split("\t\n") #Parte las lineas por tabulador y salto de línea y lo guarda en la variable prueba
			intermedio = prueba[0].split('\n') #Parte el primer argumento de prueba por salto de linea
			fields = intermedio[0].split('\t') #Parte el primer argumento de intermedio por los tabuladores
			Boo2 = True
			formato = ""
			for row2 in rows2:
				campos = row2.split("\t")
				if campos[0] == fields[locusIndex]:
					campos.pop(0)
					for i in range(len(fields)):
						formato = formato + fields[i] + "\t"
					for j in range(len(campos)):
						if j == len(campos) - 1:
							formato = formato + campos[j]
						else:
							formato = formato + campos[j] + "\t"
					rows3.append(formato)
					Boo2 = False
					break
				else:
					continue
			while Boo2 == True:
				for i in range(len(fields)):
					formato = formato + fields[i] + "\t"
				formato = formato + "IGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\tIGV\n"
				rows3.append(formato)
				Boo2 = False

		for line in rows3:
			f_final.write(line)

		f_final.close()


	def Filtrado(this,caso,caso2,caso3,padre,madre,sexo,sexo2,sexo3):

		#Generamos el archivo latest común
		
		path_latest_1=get_path_latest(caso,this.wd,this._filter)
		path_latest_2=get_path_latest(caso2,this.wd,this._filter)
		path_latest_3=get_path_latest(caso3,this.wd,this._filter)
		
		latest_1 = pd.read_csv(path_latest_1, sep="\t", header=0)
		latest_2 = pd.read_csv(path_latest_2, sep="\t", header=0)
		latest_3 = pd.read_csv(path_latest_3, sep="\t", header=0)
		
		joint_latest = pd.concat([latest_1,latest_2,latest_3], ignore_index=True).drop_duplicates(['Locus'])
		
		joint_latest[['chr','position']]=joint_latest.Locus.str.split(':',expand=True)
		joint_latest['chr']=joint_latest.chr.str.replace('chr','')
		joint_latest['chr']=joint_latest.chr.str.replace('X','30')
		joint_latest['chr']=joint_latest.chr.str.replace('Y','31')
		joint_latest['chr']=joint_latest.chr.str.replace('M','32')
		joint_latest['chr']=joint_latest['chr'].astype(int)
		joint_latest['position']=joint_latest['position'].astype(int)

		joint_latest_sort = joint_latest.sort_values(by=['chr','position'])
		joint_latest_sort = joint_latest_sort.drop(['chr','position','reliabilities'], axis=1)
		joint_latest_sort = joint_latest_sort.reset_index(drop=True)
		
		output_dir=get_output_dir(this.wd,this._filter)
		
		file_joint_latest = (output_dir+"/"+ str(caso) + "+" + str(caso2) + "+" + str(caso3) + "_latest.tsv")
		
		joint_latest_sort.to_csv(file_joint_latest, sep='\t')  

		# Busca el latest en samples

		path_latest=get_path_latest_3P(caso,caso2,caso3,this.wd,this._filter)
		quinteto_file=get_path_quinteto_file_3P(caso,caso2,caso3,this.quintetos)
		my_path_tmp=get_path_tmp(this.wd,this._filter)

		try:
			f = open(path_latest) #Guardamos en f el archivo de la ruta latest
		except (OSError, IOError) as e:
			raise ErrorInputFileLatest(path_latest)

		# Busca el quinteto en análisis

		try:
			f_quinteto = open(quinteto_file) #Guardamos en f quinteto el archivo quinteto

		except (OSError, IOError) as e:
			raise ErrorInputFileQuinteto(quinteto_file)

		f_intermedio = open(str(my_path_tmp) + "/quinteto_intermedio_"+ str(caso) + "+" + str(caso2) + "+" + str(caso3) + ".tsv", "w+") #Guardamos en f_intermedio el archivo intermedio

		locusIndex = LocusIndex(f) #Con la función LocusIndex genero el index del Locus

		locus = []
		quinteto = []
		lines = f.readlines() #Guardo las lineas del archivo latest en la variable lines

		f.close()

		for line in lines:
			fields = line.split('\t') #Guardo en la variable fields los campos de las lineas del archivo latest
			locus.append(fields[locusIndex]) #Se añade el index del locus a los campos del archivo latest


		i = 0

		find_header=False;
		header_proband=None

		while(not find_header): #Indica como header la primera linea del vcf que comienza con CHROM
			header_proband = f_quinteto.readline()
			if header_proband.startswith("#CHROM"):
				find_header=True

		try:
			probandIndex,proband2Index,proband3Index,fatherIndex,motherIndex=getIndexFamily_5_3P(caso,caso2,caso3,padre,madre,header_proband) #Guarda el indice de cada miembro 
		except NotFoundFamily as e:
			this.__log.error(e)
			return False

		lines_quinteto = f_quinteto.readlines() 

		f_quinteto.close()

		for line in lines_quinteto:
			fields_quinteto = line.split("\t") 
			data_proband = fields_quinteto[probandIndex].split('\n')[0].split(":") #Guarda en la variable data_proband los campos del quinteto del probando partidos por salto de línea y :
			data_proband2 = fields_quinteto[proband2Index].split('\n')[0].split(":") #Lo mismo con el segundo probando2
			data_proband3 = fields_quinteto[proband3Index].split('\n')[0].split(":") #Lo mismo con el segundo probando3
			data_mother = fields_quinteto[motherIndex].split('\n')[0].split(":") #Lo mismo con la madre
			data_father = fields_quinteto[fatherIndex].split('\n')[0].split(":") #Lo mismo con el padre

			
			#Prepara los datos para genotipar a los miembros de la familia
			
			if str(data_mother[0]) == '1' or str(data_mother[0]) == '0' or str(data_mother[0]) == '2':
				data_mother[0] = str(data_mother[0]) + '/' + str(data_mother[0])
			
			if str(data_father[0]) == '1' or str(data_father[0]) == '0' or str(data_father[0]) == '2':
				data_father[0] = str(data_father[0]) + '/' + str(data_father[0])

			if str(data_proband[0]) == '1' or str(data_proband[0]) == '0' or str(data_proband[0]) == '2':
				data_proband[0] = str(data_proband[0]) + '/' + str(data_proband[0])
				
			if str(data_proband2[0]) == '1' or str(data_proband2[0]) == '0' or str(data_proband2[0]) == '2':
				data_proband2[0] = str(data_proband2[0]) + '/' + str(data_proband2[0])
				
			if str(data_proband3[0]) == '1' or str(data_proband3[0]) == '0' or str(data_proband3[0]) == '2':
				data_proband3[0] = str(data_proband3[0]) + '/' + str(data_proband3[0])


			if '/' in data_proband[0]:
				genotype_proband = data_proband[0].split("/")
				genotype_proband_0 = genotype_proband[0]
				genotype_proband_1 = genotype_proband[1]
			elif '|' in data_proband[0]:
				genotype_proband = data_proband[0].split("|")
				genotype_proband_0 = genotype_proband[0]
				genotype_proband_1 = genotype_proband[1]
				
			if '/' in data_proband2[0]:
				genotype_proband2 = data_proband2[0].split("/")
				genotype_proband2_0 = genotype_proband2[0]
				genotype_proband2_1 = genotype_proband2[1]
			elif '|' in data_proband2[0]:
				genotype_proband2 = data_proband2[0].split("|")
				genotype_proband2_0 = genotype_proband2[0]
				genotype_proband2_1 = genotype_proband2[1]
				
			if '/' in data_proband3[0]:
				genotype_proband3 = data_proband3[0].split("/")
				genotype_proband3_0 = genotype_proband3[0]
				genotype_proband3_1 = genotype_proband3[1]
			elif '|' in data_proband3[0]:
				genotype_proband3 = data_proband3[0].split("|")
				genotype_proband3_0 = genotype_proband3[0]
				genotype_proband3_1 = genotype_proband3[1]

			if '/' in data_mother[0]:
				genotype_mother = data_mother[0].split("/")
				genotype_mother_0 = genotype_mother[0]
				genotype_mother_1 = genotype_mother[1]
			elif '|' in data_mother[0]:
				genotype_mother = data_mother[0].split("|")
				genotype_mother_0 = genotype_mother[0]
				genotype_mother_1 = genotype_mother[1]

			if '/' in data_father[0]:
				genotype_father = data_father[0].split("/")
				genotype_father_0 = genotype_father[0]
				genotype_father_1 = genotype_father[1]
			elif '|' in data_father[0]:
				genotype_father = data_father[0].split("|")
				genotype_father_0 = genotype_father[0]
				genotype_father_1 = genotype_father[1]


			#Guarda la cobertura de los miembros de la familia si la tienen
			if len(data_proband) < 3:
				allele_coverage_proband = '-,-'
			elif len(data_proband[1].split(',')) == 1 and genotype_proband != './.':
				allele_coverage_proband = data_proband[1] + ',0'
			else:
				allele_coverage_proband = data_proband[1]
			
			if len(data_proband2) < 3:
				allele_coverage_proband2 = '-,-'
			elif len(data_proband2[1].split(',')) == 1 and genotype_proband2 != './.':
				allele_coverage_proband2 = data_proband2[1] + ',0'
			else:
				allele_coverage_proband2 = data_proband2[1]
				
			if len(data_proband3) < 3:
				allele_coverage_proband3 = '-,-'
			elif len(data_proband3[1].split(',')) == 1 and genotype_proband3 != './.':
				allele_coverage_proband3 = data_proband3[1] + ',0'
			else:
				allele_coverage_proband3 = data_proband3[1]
			
			if len(data_mother) < 3:
				allele_coverage_mother = '-,-'
			elif len(data_mother[1].split(',')) == 1 and genotype_mother != './.':
				allele_coverage_mother = data_mother[1] + ',0'
			else:
				allele_coverage_mother = data_mother[1]

			if len(data_father) < 3:
				allele_coverage_father = '-,-'
			elif len(data_father[1].split(',')) == 1 and genotype_father != './.':
				allele_coverage_father = data_father[1] + ',0'
			else:
				allele_coverage_father = data_father[1]
				


			ref_split = fields_quinteto[4].split(",")
			lista = range(len(ref_split))

			# HACERLO INDEXANDO VALORES Y QUE BUSQUE EN REF CON EL MISMO INDEX

			for field in locus:
				fields_locus = field.split(":")
				if fields_quinteto[0] == fields_locus[0] and fields_quinteto[1] == fields_locus[1]:
					for h in lista:
						#print field #TODAS ESTAS RAMAS DE IF AÑADEN LA INFO DE COBERTURA Y GENOTIPO
						if len(data_proband) < 3:
							if len(data_proband2) <3:
								if len(data_mother) < 3:
									if len(data_father) <3:
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										break
									else:
										if data_father[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + '-' + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										break
								else:
									if data_mother[3] == '.':
										if len(data_father) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
										
										elif data_father [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + data_mother[3] + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										break
									else:
										if len(data_father) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + "-" + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										
										elif data_father [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
										break
										
							elif data_proband2[3] == '.':
								if len(data_mother) <3:
									if len(data_father) <3:
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										break
									else:
										if data_father[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + '-' + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										break
								else:
									if data_mother[3] == '.':
										if len(data_father) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
										
										elif data_father [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + data_mother[3] + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										break
									else:
										if len(data_father) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + "-" + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										
										elif data_father [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
										break
									
							else:
								if len(data_proband3) <3:
									if len(data_mother) <3:
										if len(data_father) <3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
										else:
											if data_father[3] == '.':
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + '-' + "\t" + \
													  data_father[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											else:
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											break
									else:
										if data_mother[3] == '.':
											if len(data_father) <3:
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + \
													  data_mother[3] + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
												break
											else:
												if data_father[3] == '.':
													formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
														  str(int(allele_coverage_proband2.split(',')[0]) + int(
														  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + \
														  data_mother[3] + "\t" + data_father[3] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
												else:
													formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
														  str(int(allele_coverage_proband2.split(',')[0]) + int(
														  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + \
														  data_mother[3] + "\t" + \
														  str(int(allele_coverage_father.split(',')[0]) + int(
														  	  allele_coverage_father.split(',')[1])) + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]   
												break
										else:
											if len(data_father) < 3: 
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											
											elif data_father [3] == '.':
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  data_father[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
											else: 
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
											break
			
								elif data_proband3[3] == '.':
									if len(data_mother) < 3:
										if len(data_father) <3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_proband3[3] + "\t" + \
												  '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
										elif data_father[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_proband3[3] + "\t" + \
												  '-' + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_proband3[3] + "\t" + \
												  '-' + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
									elif data_mother[3] == '.':
										if len(data_father) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_proband3[3] + "\t" + \
												  data_mother[3] + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
												
										elif data_father [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_proband3[3] + "\t" + \
												  data_mother[3] + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_proband3[3] + "\t" + \
												  data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
									else:
										if len(data_father) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_proband3[3] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
												
										elif data_father [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_proband3[3] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_proband3[3] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
											break
												
								else:
									if len(data_mother) <3:
										if len(data_father) <3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											break
										elif data_father[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  '-' + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
									elif data_mother[3] == '.':
										if len(data_father) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  data_mother[3] + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
												
										elif data_father [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  data_mother[3] + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
									else:
										if len(data_father) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
												
										elif data_father [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
											break
									
						elif data_proband[3] == '.':
							if len(data_proband2) <3:
								if len(data_mother) < 3:
									if len(data_father) <3:
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										break
									else:
										if data_father[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + '-' + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										break
								else:
									if data_mother[3] == '.':
										if len(data_father) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
										
										elif data_father [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + data_mother[3] + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										break
									else:
										if len(data_father) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + "-" + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										
										elif data_father [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
										break
										
							elif data_proband2[3] == '.':
								if len(data_mother) <3:
									if len(data_father) <3:
										formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										break
									else:
										if data_father[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + '-' + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										break
								else:
									if data_mother[3] == '.':
										if len(data_father) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
										
										elif data_father [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + data_mother[3] + "\t" + \
												  data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										break
									else:
										if len(data_father) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + "-" + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										
										elif data_father [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
										break
									
							else:
								if len(data_proband3) <3:
									if len(data_mother) <3:
										if len(data_father) <3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
										else:
											if data_father[3] == '.':
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + '-' + "\t" + \
													  data_father[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											else:
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											break
									else:
										if data_mother[3] == '.':
											if len(data_father) <3:
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + \
													  data_mother[3] + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
												break
											else:
												if data_father[3] == '.':
													formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
														  str(int(allele_coverage_proband2.split(',')[0]) + int(
														  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + \
														  data_mother[3] + "\t" + data_father[3] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
												else:
													formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
														  str(int(allele_coverage_proband2.split(',')[0]) + int(
														  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + \
														  data_mother[3] + "\t" + \
														  str(int(allele_coverage_father.split(',')[0]) + int(
														  	  allele_coverage_father.split(',')[1])) + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]   
												break
										else:
											if len(data_father) < 3: 
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											
											elif data_father [3] == '.':
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  data_father[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
											else: 
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
											break
			
								elif data_proband3[3] == '.':
									if len(data_mother) < 3:
										if len(data_father) <3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_proband3[3] + "\t" + \
												  '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
										elif data_father[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_proband3[3] + "\t" + \
												  '-' + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_proband3[3] + "\t" + \
												  '-' + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
									elif data_mother[3] == '.':
										if len(data_father) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + '-' + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_proband3[3] + "\t" + \
												  data_mother[3] + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
												
										elif data_father [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_proband3[3] + "\t" + \
												  data_mother[3] + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_proband3[3] + "\t" + \
												  data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
									else:
										if len(data_father) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_proband3[3] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
												
										elif data_father [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_proband3[3] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + data_proband3[3] + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
											break
												
								else:
									if len(data_mother) <3:
										if len(data_father) <3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											break
										elif data_father[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  '-' + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
									elif data_mother[3] == '.':
										if len(data_father) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  data_mother[3] + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
												
										elif data_father [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  data_mother[3] + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
									else:
										if len(data_father) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
												
										elif data_father [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + data_proband[3] + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
											break

						else:
							if len(data_proband2) < 3:
								if len(data_proband3) <3:
									if len(data_mother) <3:
										if len(data_father) <3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  '-' + "\t" + '-' + "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
										else:
											if data_father[3] == '.':
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  '-' + "\t" + '-' + "\t" + '-' + "\t" + data_father[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											else:
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  '-' + "\t" + '-' + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											break
									else:
										if data_mother[3] == '.':
											if len(data_father) <3:
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  '-' + "\t" + '-' + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
												break
											else:
												if data_father[3] == '.':
													formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int(
														  	  allele_coverage_proband.split(',')[1])) + "\t" + \
														  '-' + "\t" + '-' + "\t" + data_mother[3] + "\t" + data_father[3] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
												else:
													formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int(
														  	  allele_coverage_proband.split(',')[1])) + "\t" + \
														  '-' + "\t" + '-' + "\t" + data_mother[3] + "\t" + \
														  str(int(allele_coverage_father.split(',')[0]) + int(
														  	  allele_coverage_father.split(',')[1])) + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]   
												break
										else:
											if len(data_father) < 3: 
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  '-' + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											
											elif data_father [3] == '.':
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  '-' + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  data_father[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
											else: 
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  '-' + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
											break
			
								elif data_proband3[3] == '.':
									if len(data_mother) <3:
										if len(data_father) <3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  '-' + "\t" + data_proband3[3] + "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
										else:
											if data_father[3] == '.':
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  '-' + "\t" + data_proband3[3] + "\t" + '-' + "\t" + data_father[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											else:
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  '-' + "\t" + data_proband3[3] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											break
									else:
										if data_mother[3] == '.':
											if len(data_father) <3:
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  '-' + "\t" + data_proband3[3] + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
												break
											else:
												if data_father[3] == '.':
													formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int(
														  	  allele_coverage_proband.split(',')[1])) + "\t" + \
														  '-' + "\t" + data_proband3[3] + "\t" + data_mother[3] + "\t" + data_father[3] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
												else:
													formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int(
														  	  allele_coverage_proband.split(',')[1])) + "\t" + \
														  '-' + "\t" + data_proband3[3] + "\t" + data_mother[3] + "\t" + \
														  str(int(allele_coverage_father.split(',')[0]) + int(
														  	  allele_coverage_father.split(',')[1])) + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]   
												break
										else:
											if len(data_father) < 3: 
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  '-' + "\t" + data_proband3[3] + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											
											elif data_father [3] == '.':
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  '-' + "\t" + data_proband3[3] + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  data_father[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
											else: 
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  '-' + "\t" + data_proband3[3] + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
											break
												
								else:
									if len(data_mother) <3:
										if len(data_father) <3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											break
										elif data_father[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  '-' + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]   
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											break
									elif data_mother[3] == '.':
										if len(data_father) <3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  data_mother[3] + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											break
										elif data_father[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  data_mother[3] + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]   
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											break
									else:
										if len(data_father) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
												
										elif data_father [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break

							elif data_proband2[3] == '.':
								if len(data_proband3) <3:
									if len(data_mother) <3:
										if len(data_father) <3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  data_proband2[3] + "\t" + '-' + "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
										else:
											if data_father[3] == '.':
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  data_proband2[3] + "\t" + '-' + "\t" + '-' + "\t" + data_father[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											else:
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  data_proband2[3] + "\t" + '-' + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											break
									else:
										if data_mother[3] == '.':
											if len(data_father) <3:
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  data_proband2[3] + "\t" + '-' + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
												break
											else:
												if data_father[3] == '.':
													formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int(
														  	  allele_coverage_proband.split(',')[1])) + "\t" + \
														  data_proband2[3] + "\t" + '-' + "\t" + data_mother[3] + "\t" + data_father[3] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
												else:
													formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int(
														  	  allele_coverage_proband.split(',')[1])) + "\t" + \
														  data_proband2[3] + "\t" + '-' + "\t" + data_mother[3] + "\t" + \
														  str(int(allele_coverage_father.split(',')[0]) + int(
														  	  allele_coverage_father.split(',')[1])) + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]   
												break
										else:
											if len(data_father) < 3: 
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  data_proband2[3] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											
											elif data_father [3] == '.':
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  data_proband2[3] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  data_father[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
											else: 
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  data_proband2[3] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
											break
			
								elif data_proband3[3] == '.':
									if len(data_mother) <3:
										if len(data_father) <3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  data_proband2[3] + "\t" + data_proband3[3] + "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
										else:
											if data_father[3] == '.':
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  data_proband2[3] + "\t" + data_proband3[3] + "\t" + '-' + "\t" + data_father[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											else:
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  data_proband2[3] + "\t" + data_proband3[3] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											break
									else:
										if data_mother[3] == '.':
											if len(data_father) <3:
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  data_proband2[3] + "\t" + data_proband3[3] + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
												break
											else:
												if data_father[3] == '.':
													formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int(
														  	  allele_coverage_proband.split(',')[1])) + "\t" + \
														  data_proband2[3] + "\t" + data_proband3[3] + "\t" + \
														  data_mother[3] + "\t" + data_father[3] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
												else:
													formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int(
														  	  allele_coverage_proband.split(',')[1])) + "\t" + \
														  data_proband2[3] + "\t" + data_proband3[3] + "\t" + data_mother[3] + "\t" + \
														  str(int(allele_coverage_father.split(',')[0]) + int(
														  	  allele_coverage_father.split(',')[1])) + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]   
												break
										else:
											if len(data_father) < 3: 
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  data_proband2[3] + "\t" + data_proband3[3] + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											
											elif data_father [3] == '.':
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  data_proband2[3] + "\t" + data_proband3[3] + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  data_father[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
											else: 
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  data_proband2[3] + "\t" + data_proband3[3] + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
											break
												
								else:
									if len(data_mother) <3:
										if len(data_father) <3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											break
										elif data_father[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  '-' + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]   
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											break
									elif data_mother[3] == '.':
										if len(data_father) <3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  data_mother[3] + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											break
										elif data_father[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  data_mother[3] + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]   
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											break
									else:
										if len(data_father) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
												
										elif data_father [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  data_proband2[3] + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
							else:
								if len(data_proband3) <3:
									if len(data_mother) <3:
										if len(data_father) <3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  '-' + "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
										else:
											if data_father[3] == '.':
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  '-' + "\t" + '-' + "\t" + data_father[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]   
											else:
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  '-' + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]   
											break
									else:
										if data_mother[3] == '.':
											if len(data_father) <3:
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  '-' + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
												break
											else:
												if data_father[3] == '.':
													formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int(
														  	  allele_coverage_proband.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_proband2.split(',')[0]) + int(
														  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
														  '-' + "\t" + data_mother[3] + "\t" + data_father[3] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
												else:
													formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int(
														  	  allele_coverage_proband.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_proband2.split(',')[0]) + int(
														  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
														  '-' + "\t" + data_mother[3] + "\t" + \
														  str(int(allele_coverage_father.split(',')[0]) + int(
														  	  allele_coverage_father.split(',')[1])) + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]   
												break
										else:
											if len(data_father) <3:
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  '-' + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
												break
											else:
												if data_father[3] == '.':
													formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int(
														  	  allele_coverage_proband.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_proband2.split(',')[0]) + int(
														  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
														  '-' + "\t" + \
														  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  	  allele_coverage_mother.split(',')[1])) + "\t" + data_father[3] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
												else:
													formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int(
														  	  allele_coverage_proband.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_proband2.split(',')[0]) + int(
														  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
														  '-' + "\t" + \
														  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  	  allele_coverage_mother.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_father.split(',')[0]) + int(
														  	  allele_coverage_father.split(',')[1])) + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]   
												break
			
								elif data_proband3[3] == '.':
									if len(data_mother) <3:
										if len(data_father) <3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  data_proband3[3] + "\t" + '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break
										else:
											if data_father[3] == '.':
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  data_proband3[3] + "\t" + '-' + "\t" + data_father[3] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]   
											else:
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  data_proband3[3] + "\t" + '-' + "\t" + \
													  str(int(allele_coverage_father.split(',')[0]) + int(
													  	  allele_coverage_father.split(',')[1])) + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]   
											break
									else:
										if data_mother[3] == '.':
											if len(data_father) <3:
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  data_proband3[3] + "\t" + data_mother[3] + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
												break
											else:
												if data_father[3] == '.':
													formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int(
														  	  allele_coverage_proband.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_proband2.split(',')[0]) + int(
														  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
														  data_proband3[3] + "\t" + data_mother[3] + "\t" + data_father[3] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
												else:
													formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int(
														  	  allele_coverage_proband.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_proband2.split(',')[0]) + int(
														  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
														  data_proband3[3] + "\t" + data_mother[3] + "\t" + \
														  str(int(allele_coverage_father.split(',')[0]) + int(
														  	  allele_coverage_father.split(',')[1])) + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]   
												break
										else:
											if len(data_father) <3:
												formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
													  str(int(allele_coverage_proband.split(',')[0]) + int(
													  	  allele_coverage_proband.split(',')[1])) + "\t" + \
													  str(int(allele_coverage_proband2.split(',')[0]) + int(
													  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
													  data_proband3[3] + "\t" + \
													  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  allele_coverage_mother.split(',')[1])) + "\t" + '-' + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
													  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
													  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
												break
											else:
												if data_father[3] == '.':
													formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int(
														  	  allele_coverage_proband.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_proband2.split(',')[0]) + int(
														  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
														  data_proband3[3] + "\t" + \
														  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  	  allele_coverage_mother.split(',')[1])) + "\t" + data_father[3] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]
												else:
													formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
														  str(int(allele_coverage_proband.split(',')[0]) + int(
														  	  allele_coverage_proband.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_proband2.split(',')[0]) + int(
														  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
														  data_proband3[3] + "\t" + \
														  str(int(allele_coverage_mother.split(',')[0]) + int(
													  	  	  allele_coverage_mother.split(',')[1])) + "\t" + \
														  str(int(allele_coverage_father.split(',')[0]) + int(
														  	  allele_coverage_father.split(',')[1])) + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
														  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
														  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]   
												break
												
								else:
									if len(data_mother) <3:
										if len(data_father) <3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  '-' + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											break
										elif data_father[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  '-' + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]   
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											break
									elif data_mother[3] == '.':
										if len(data_father) <3:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  data_mother[3] + "\t" + '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											break
										elif data_father[3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  data_mother[3] + "\t" + data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]   
										else:
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  data_mother[3] + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
											break
									else:
										if len(data_father) < 3: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  '-' + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
												
										elif data_father [3] == '.':
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  data_father[3] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1]  
										else: 
											formato = fields_quinteto[0] + ":" + fields_quinteto[1] + "\t" + \
												  str(int(allele_coverage_proband.split(',')[0]) + int(
												  	  allele_coverage_proband.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband2.split(',')[0]) + int(
												  	  allele_coverage_proband2.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_proband3.split(',')[0]) + int(
												  	  allele_coverage_proband3.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_mother.split(',')[0]) + int(
												  	  allele_coverage_mother.split(',')[1])) + "\t" + \
												  str(int(allele_coverage_father.split(',')[0]) + int(
												  	  allele_coverage_father.split(',')[1])) + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband2.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband2.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_proband3.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_proband3.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_mother.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_mother.split(',')[1] + "\t" + \
												  fields_quinteto[3] + '=' + allele_coverage_father.split(',')[0] + ', ' + \
												  fields_quinteto[4] + '=' + allele_coverage_father.split(',')[1] 
											break				
					for h in lista:

						if genotype_proband[0] == "0":
							formato = formato + "\t" + fields_quinteto[3] + "/"
							break
						elif genotype_proband[0] == str(h + 1):
							formato = formato + "\t" + ref_split[h] + "/"
							break
						elif genotype_proband[0] == ".":
							formato = formato + "\t" + genotype_proband[0] + "/"
							break
					for h in lista:

						if genotype_proband[1] == "0":
							formato = formato + fields_quinteto[3] + "\t"
							break
						elif genotype_proband[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_proband[1] == ".":
							formato = formato + genotype_proband[1] + "\t"
							break
					for h in lista:

						if genotype_proband2[0] == "0":
							formato = formato + fields_quinteto[3] + "/"
							break
						elif genotype_proband2[0] == str(h + 1):
							formato = formato + ref_split[h] + "/"
							break
						elif genotype_proband2[0] == ".":
							formato = formato + genotype_proband2[0] + "/"
							break
					for h in lista:

						if genotype_proband2[1] == "0":
							formato = formato + fields_quinteto[3] + "\t"
							break
						elif genotype_proband2[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_proband2[1] == ".":
							formato = formato + genotype_proband2[1] + "\t"
							break
					for h in lista:

						if genotype_proband3[0] == "0":
							formato = formato + fields_quinteto[3] + "/"
							break
						elif genotype_proband3[0] == str(h + 1):
							formato = formato + ref_split[h] + "/"
							break
						elif genotype_proband3[0] == ".":
							formato = formato + genotype_proband3[0] + "/"
							break
					for h in lista:

						if genotype_proband3[1] == "0":
							formato = formato + fields_quinteto[3] + "\t"
							break
						elif genotype_proband3[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_proband3[1] == ".":
							formato = formato + genotype_proband3[1] + "\t"
							break
					for h in lista:

						if genotype_mother[0] == "0":
							formato = formato + fields_quinteto[3] + "/"
							break
						elif genotype_mother[0] == str(h + 1):
							formato = formato + ref_split[h] + "/"
							break
						elif genotype_mother[0] == ".":
							formato = formato + genotype_mother[0] + "/"
							break
					for h in lista:

						if genotype_mother[1] == "0":
							formato = formato + fields_quinteto[3] + "\t"
							break
						elif genotype_mother[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_mother[1] == ".":
							formato = formato + genotype_mother[1] + "\t"
							break
					for h in lista:
						if genotype_father[0] == "0":
							formato = formato + fields_quinteto[3] + "/"
							break
						elif genotype_father[0] == str(h + 1):
							formato = formato + ref_split[h] + "/"
							break
						elif genotype_father[0] == ".":
							formato = formato + genotype_father[0] + "/"
							break
					for h in lista:
						if genotype_father[1] == "0":
							formato = formato + fields_quinteto[3] + "\t"
							break
						elif genotype_father[1] == str(h + 1):
							formato = formato + ref_split[h] + "\t"
							break
						elif genotype_father[1] == ".":
							formato = formato + genotype_father[1] + "\t"
							break		
							
					#Añadimos la zigosidad
					
					for h in lista:
						if genotype_proband[0] == "0" and genotype_proband[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_proband[0] == "0" and genotype_proband[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t" 
							break
						elif genotype_proband[0] == str(h + 1) and genotype_proband[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_proband[0] == "." or genotype_proband[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
							
					for h in lista:
						if genotype_proband2[0] == "0" and genotype_proband2[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_proband2[0] == "0" and genotype_proband2[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t" 
							break
						elif genotype_proband2[0] == str(h + 1) and genotype_proband2[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_proband2[0] == "." or genotype_proband2[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
							
					for h in lista:
						if genotype_proband3[0] == "0" and genotype_proband3[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_proband3[0] == "0" and genotype_proband3[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t" 
							break
						elif genotype_proband3[0] == str(h + 1) and genotype_proband3[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_proband3[0] == "." or genotype_proband3[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
							
					for h in lista:
						if genotype_mother[0] == "0" and genotype_mother[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_mother[0] == "0" and genotype_mother[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t" 
							break
						elif genotype_mother[0] == str(h + 1) and genotype_mother[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_mother[0] == "." or genotype_mother[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
							
					for h in lista:
						if genotype_father[0] == "0" and genotype_father[1] == '0':
							formato = formato +  "WT" + "\t"
							break
						elif genotype_father[0] == "0" and genotype_father[1] == str(h + 1):
							formato = formato + "Heteroz" + "\t"
							break
						elif genotype_father[0] == str(h + 1) and genotype_father[1] == str(h + 1):
							formato = formato + "Homoz" + "\t"
							break
						elif genotype_father[0] == "." or genotype_father[1] == ".":
							formato = formato +  "NA" + "\t"
							break
						else:
							formato = formato + "NA" + "\t"
							break
						
					#Añadimos el VAF
					
					for h in lista:
						if allele_coverage_proband.split(',') != "0" and '.' in allele_coverage_proband:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_proband.split(',') != "0" and '-' in allele_coverage_proband:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_proband.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_proband.split(',')[1])/ (float(allele_coverage_proband.split(',')[0]) + float(allele_coverage_proband.split(',')[1]))) + "\t" 
							break
						elif allele_coverage_proband.split(',')[0] == "0" and allele_coverage_proband.split(',')[1] == "0":
							formato = formato + "0.0" + "\t"
							break
						else:
							formato = formato + "1.0" + "\t"
							break
					for h in lista:
						if allele_coverage_proband2.split(',') != "0" and '.' in allele_coverage_proband2:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_proband2.split(',') != "0" and '-' in allele_coverage_proband2:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_proband2.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_proband2.split(',')[1])/ (float(allele_coverage_proband2.split(',')[0]) + float(allele_coverage_proband2.split(',')[1]))) + "\t" 
							break
						elif allele_coverage_proband2.split(',')[0] == "0" and allele_coverage_proband2.split(',')[1] == "0":
							formato = formato + "0.0" + "\t"
							break
						else:
							formato = formato + "1.0" + "\t"
							break
							
					for h in lista:
						if allele_coverage_proband3.split(',') != "0" and '.' in allele_coverage_proband3:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_proband3.split(',') != "0" and '-' in allele_coverage_proband3:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_proband3.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_proband3.split(',')[1])/ (float(allele_coverage_proband3.split(',')[0]) + float(allele_coverage_proband3.split(',')[1]))) + "\t" 
							break
						elif allele_coverage_proband3.split(',')[0] == "0" and allele_coverage_proband3.split(',')[1] == "0":
							formato = formato + "0.0" + "\t"
							break
						else:
							formato = formato + "1.0" + "\t"
							break
					for h in lista:
						if allele_coverage_mother.split(',') != "0" and '.' in allele_coverage_mother:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_mother.split(',') != "0" and '-' in allele_coverage_mother:
							formato = formato + ".,." + "\t"
							break
						elif allele_coverage_mother.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_mother.split(',')[1])/ (float(allele_coverage_mother.split(',')[0]) + float(allele_coverage_mother.split(',')[1]))) + "\t" 
							break
						elif allele_coverage_mother.split(',')[0] == "0" and allele_coverage_mother.split(',')[1] == "0":
							formato = formato + "0.0" + "\t"
							break
						else:
							formato = formato + "1.0" + "\t"
							break	
					for h in lista:
						if allele_coverage_father.split(',') != "0" and '.' in allele_coverage_father:
							formato = formato + ".,." 
							break
						elif allele_coverage_father.split(',') != "0" and '-' in allele_coverage_father:
							formato = formato + ".,." 
							break
						elif allele_coverage_father.split(',')[0] != "0":
							formato = formato + str(float(allele_coverage_father.split(',')[1])/ (float(allele_coverage_father.split(',')[0]) + float(allele_coverage_father.split(',')[1]))) 
							break
						elif allele_coverage_father.split(',')[0] == "0" and allele_coverage_father.split(',')[1] == "0":
							formato = formato + "0.0" 
							break
						else:
							formato = formato + "1.0" 
							break
							
					#Añadimos el origen

					if data_proband[0] == data_mother[0] and data_proband[0] == data_father[0] and data_proband[0] == data_proband3[0] and data_proband[0] == data_proband2[0]:
						formato = formato + "\t" + "All"
					elif data_proband[0] == data_mother[0] and data_proband2[0] == data_mother[0] and data_proband3[0] == data_mother[0]:
						formato = formato + "\t" + "Proband/Proband2/Proband3/Mother"
					elif data_proband[0] == data_father[0] and data_proband2[0] == data_father[0] and data_proband3[0] == data_father[0]:
						formato = formato + "\t" + "Proband/Proband2/Proband3/Father"
					elif data_proband[0] == data_mother[0] and data_proband[0] == data_father[0] and data_proband[0] == data_proband2[0]:
						formato = formato + "\t" + "Proband/Proband2/Mother/Father"
					elif data_proband[0] == data_mother[0] and data_proband[0] == data_father[0] and data_proband[0] == data_proband3[0]:
						formato = formato + "\t" + "Proband/Proband3/Mother/Father"
					elif data_proband2[0] == data_mother[0] and data_proband2[0] == data_father[0] and data_proband2[0] == data_proband3[0]:
						formato = formato + "\t" + "Proband2/Proband3/Mother/Father"
					elif data_proband[0] == data_mother[0] and data_proband2[0] == data_mother[0]:
						formato = formato + "\t" + "Proband/Proband2/Mother"
					elif data_proband[0] == data_mother[0] and data_proband3[0] == data_mother[0]:
						formato = formato + "\t" + "Proband/Proband3/Mother"
					elif data_proband2[0] == data_mother[0] and data_proband3[0] == data_mother[0]:
						formato = formato + "\t" + "Proband2/Proband3/Mother"
					elif data_proband[0] == data_father[0] and data_proband2[0] == data_father[0]:
						formato = formato + "\t" + "Proband/Proband2/Father"
					elif data_proband[0] == data_father[0] and data_proband3[0] == data_father[0]:
						formato = formato + "\t" + "Proband/Proband3/Father"
					elif data_proband2[0] == data_father[0] and data_proband3[0] == data_father[0]:
						formato = formato + "\t" + "Proband2/Proband3/Father"
					elif data_proband[0] == data_mother[0]:
						formato = formato + "\t" + "Proband/Mother"
					elif data_proband2[0] == data_mother[0]:
						formato = formato + "\t" + "Proband2/Mother"
					elif data_proband3[0] == data_mother[0]:
						formato = formato + "\t" + "Proband3/Mother"
					elif data_proband[0] == data_father[0]:
						formato = formato + "\t" + "Proband/Father"
					elif data_proband2[0] == data_father[0]:
						formato = formato + "\t" + "Proband2/Father"
					elif data_proband3[0] == data_father[0]:
						formato = formato + "\t" + "Proband3/Father"
					elif data_proband[0] == data_proband2[0] and data_proband[0] == data_proband[3]:
						formato = formato + "\t" + "Proband/Proband2/Proband3"
					elif data_proband2[0] == data_proband2[0]:
						formato = formato + "\t" + "Proband/Proband2"
					elif data_proband[0] == data_proband3[0]:
						formato = formato + "\t" + "Proband/Proband3"
					elif data_proband2[0] == data_proband3[0]:
						formato = formato + "\t" + "Proband2/Proband3"
					else:
						formato = formato + "\t" + "Check"

					quinteto.append(formato)

		header_latest = "Locus" + "\t" +  str(caso) + " Coverage" + "\t" + str(caso2) + " Coverage" + "\t" +  str(caso3) + " Coverage" + "\t" + "Mother Depth" +"\t" + "Father Depth" + "\t" + str(caso) + " AlleleCoverage" + "\t" + str(caso2) + " AlleleCoverage" + "\t" + str(caso3) + " AlleleCoverage" + "\t" + "Mother AlleleCoverage" + "\t" + "Father AlleleCoverage" + "\t" + str(caso) + " Genotype" + "\t" + str(caso2) + " Genotype" + "\t" + str(caso3) + " Genotype" + "\t" + "Mother Genotype" + "\t" + "Father Genotype" + "\t" + "zigosity " + str(caso) + "\t" + "zigosity " + str(caso2) + "\t" + "zigosity " + str(caso3) + "\t" + "zigosity Mother" + "\t" + "zigosity Father" + "\t" + "VAF " +  str(caso) + "\t" + "VAF " + str(caso2) + "\t" + "VAF " +  str(caso3)+ "\t" + "VAF Mother" + "\t" + "VAF Father" + "\t" + "Origin" + "\n"


		f_intermedio.write(header_latest)

		# Meto los datos de genotipos, profundidades y origen de las variantes en un archivo intermedio

		for line in quinteto:
			# ~ print line
			f_intermedio.write("%s\n" % line)
		return True

#Ejemplo de ejecución

#ejemplo = duoDragen_2P("WD", "ruta_VCFs", "ruta_latest")

#ejemplo.run("probando1", "probando2","sexo1","sexo2")
