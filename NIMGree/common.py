def get_info_family_duo(file_):
	with open(file_, 'r') as fp:  #Abre el archivo file_duo

		line = fp.readline()  # Devuelve la primera linea del archivo fp
		while(line):
			line=line.replace("\n","") #Remplaza los saltos de linea por espacios
			if(line.startswith('#CHROM')):  
				line=line.split("\t")[-2:]
				# line.sort(key = str)
				return tuple(line) #Devuelve una tupla
			line = fp.readline()             
	return None

def get_info_family(file_trio):
	with open(file_trio, 'r') as fp:  #Abre el archivo file_trio

		line = fp.readline()  # Devuelve la primera linea del archivo fp
		while(line):
			line=line.replace("\n","") #Remplaza los saltos de linea por espacios
			if(line.startswith('#CHROM')):  #Si la linea comienza por #CHROM parte los 3 últimos (que son la familia)
				line=line.split("\t")[-3:]
				# line.sort(key = str)
				return tuple(line) #Devuelve una tupla
			line = fp.readline()             
	return None
	
def get_info_family_cuarteto(file_cuarteto):
	with open(file_cuarteto, 'r') as fp:  #Abre el archivo file_cuarteto

		line = fp.readline()  # Devuelve la primera linea del archivo fp
		while(line):
			line=line.replace("\n","") #Remplaza los saltos de linea por espacios
			if(line.startswith('#CHROM')):  #Si la linea comienza por #CHROM parte los 4 últimos (que son la familia)
				line=line.split("\t")[-4:]
				# line.sort(key = str)
				return tuple(line) #Devuelve una tupla
			line = fp.readline()             
	return None

def get_info_family_5(file_quinteto):
	with open(file_quinteto, 'r') as fp:  #Abre el archivo file_quinteto

		line = fp.readline()  # Devuelve la primera linea del archivo fp
		while(line):
			line=line.replace("\n","") #Remplaza los saltos de linea por espacios
			if(line.startswith('#CHROM')):  #Si la linea comienza por #CHROM parte los 5 últimos (que son la familia)
				line=line.split("\t")[-5:]
				# line.sort(key = str)
				return tuple(line) #Devuelve una tupla
			line = fp.readline()             
	return None

def get_path_duo_file(case,path_dir_duo):
	return path_dir_duo+"/duo_"+case+".vcf"
	
def get_path_duo_file_2P(case,case2,path_dir_duo):
	return path_dir_duo+"/duo_"+case+"+"+case2+".vcf"

def get_path_trio_file(case,path_dir_trio):
	return path_dir_trio+"/trio_"+case+".vcf"

def get_path_cuarteto_file(case,path_dir_cuarteto):
	return path_dir_cuarteto+"/cuarteto_"+case+".vcf"
	
def get_path_cuarteto_file_2P(case,case2,path_dir_cuarteto):
	return path_dir_cuarteto+"/cuarteto_"+case+"+"+case2+".vcf"
	
def get_path_quinteto_file(case,path_dir_quinteto):
	return path_dir_quinteto+"/quinteto_"+case+".vcf"
	
def get_path_quinteto_file_2P(case,case2,path_dir_quinteto):
	return path_dir_quinteto+"/quinteto_"+case+"+"+case2+".vcf"
	
def get_path_quinteto_file_3P(case,case2,case3,path_dir_quinteto):
	return path_dir_quinteto+"/quinteto_"+case+"+"+case2+"+"+case3+".vcf"

def get_path_latest(case,wd,filter_apply):
	return wd+"/%s/%s_latest.tsv"%(filter_apply,case)
	
def get_path_latest_2P(case,case2,wd,filter_apply):
	return wd+"/%s/%s+%s_latest.tsv"%(filter_apply,case,case2)

def get_path_latest_3P(case,case2,case3,wd,filter_apply):
	return wd+"/%s/%s+%s+%s_latest.tsv"%(filter_apply,case,case2,case3)

def get_path_tmp(wd,filter_apply):
	return wd+"/"+filter_apply

def get_output_dir(wd,filter_apply):
	return wd+"/"+filter_apply



def LocusIndex(f):
	# Saco la columna de locus con el valor del index de la cabecera. Locus en formato: chr11:54654 

	header_proband = f.readline() #Lee f y guardalo en headerproband

	header_probandSplit = header_proband.split('\t') #parte la linea de header_proband por los tabuladores y guardalo como header_probandSplit

	return header_probandSplit.index('Locus') #devuelve el index de headerproband de locus

	# for field in header_probandSplit:
	# 	if field == "Locus":
	# 		locusIndex = header_probandSplit.index(field)

	# return locusIndex


def getIndexFamily_2(caso,parient,header_proband):
	header_probandSplit = header_proband.split('\t') #parte a header proband por los tabuladores
	
	parientIndex=-1
	probandIndex=-1

	for field in header_probandSplit:
		# print (field)
		campo = field.split("\n") #parte el campo por los saltos de linea
		if campo[0] == str(parient):
			# print (field)
			parientIndex = header_probandSplit.index(field) #Si el primer campo es el pariente, indexalo
		elif campo[0] == str(caso):
			# print (field)
			probandIndex = header_probandSplit.index(field) #si es el probando indexalo

	if parientIndex==-1:
		raise NotFoundFamily('parient',parient)

	if probandIndex==-1:
		raise NotFoundFamily('proband',caso)

	return probandIndex,parientIndex
	
def getIndexFamily_2_2P(caso,caso2,header_proband):
	header_probandSplit = header_proband.split('\t') #parte a header proband por los tabuladores
	
	proband2Index=-1
	probandIndex=-1

	for field in header_probandSplit:
		# print (field)
		campo = field.split("\n") #parte el campo por los saltos de linea
		if campo[0] == str(caso2):
			# print (field)
			proband2Index = header_probandSplit.index(field) #Si el primer campo es el probando2, indexalo
		elif campo[0] == str(caso):
			# print (field)
			probandIndex = header_probandSplit.index(field) #si es el probando indexalo

	if proband2Index==-1:
		raise NotFoundFamily('proband2',caso2)

	if probandIndex==-1:
		raise NotFoundFamily('proband',caso)

	return probandIndex,proband2Index

def getIndexFamily(caso,padre,madre,header_proband):
	header_probandSplit = header_proband.split('\t') #parte a header proband por los tabuladores
	
	fatherIndex=-1
	motherIndex=-1
	probandIndex=-1

	for field in header_probandSplit:
		# print (field)
		campo = field.split("\n") #parte el campo por los saltos de linea
		if campo[0] == str(padre):
			# print (field)
			fatherIndex = header_probandSplit.index(field) #Si el primer campo es padre, indexalo
		elif campo[0] == str(madre):
			# print (field)
			motherIndex = header_probandSplit.index(field) #Si es madre indexalo
		elif campo[0] == str(caso):
			# print (field)
			probandIndex = header_probandSplit.index(field) #si es el probando indexalo

	if fatherIndex==-1:
		raise NotFoundFamily('father',padre)

	if motherIndex==-1:
		raise NotFoundFamily('mother',madre)

	if probandIndex==-1:
		raise NotFoundFamily('proband',caso)

	return probandIndex,fatherIndex,motherIndex
	
def getIndexFamily_4(caso,padre,madre,parient,header_proband):
	header_probandSplit = header_proband.split('\t') #parte a header proband por los tabuladores
	
	fatherIndex=-1
	motherIndex=-1
	probandIndex=-1
	parientIndex=-1
	

	for field in header_probandSplit:
		# print (field)
		campo = field.split("\n") #parte el campo por los saltos de linea
		if campo[0] == str(padre):
			# print (field)
			fatherIndex = header_probandSplit.index(field) #Si el primer campo es padre, indexalo
		elif campo[0] == str(madre):
			# print (field)
			motherIndex = header_probandSplit.index(field) #Si es madre indexalo
		elif campo[0] == str(caso):
			# print (field)
			probandIndex = header_probandSplit.index(field) #si es el probando indexalo
		elif campo[0] == str(parient):
			# print (field)
			parientIndex = header_probandSplit.index(field) #si es el pariente indexalo

	if fatherIndex==-1:
		raise NotFoundFamily('father',padre)

	if motherIndex==-1:
		raise NotFoundFamily('mother',madre)

	if probandIndex==-1:
		raise NotFoundFamily('proband',caso)
		
	if parientIndex==-1:
		raise NotFoundFamily('parient',parient)

	return probandIndex,fatherIndex,motherIndex,parientIndex
	
def getIndexFamily_4_2P(caso,caso2,padre,madre,header_proband):
	header_probandSplit = header_proband.split('\t') #parte a header proband por los tabuladores
	
	fatherIndex=-1
	motherIndex=-1
	probandIndex=-1
	proband2Index=-1
	

	for field in header_probandSplit:
		# print (field)
		campo = field.split("\n") #parte el campo por los saltos de linea
		if campo[0] == str(padre):
			# print (field)
			fatherIndex = header_probandSplit.index(field) #Si el primer campo es padre, indexalo
		elif campo[0] == str(madre):
			# print (field)
			motherIndex = header_probandSplit.index(field) #Si es madre indexalo
		elif campo[0] == str(caso):
			# print (field)
			probandIndex = header_probandSplit.index(field) #si es el probando indexalo
		elif campo[0] == str(caso2):
			# print (field)
			proband2Index = header_probandSplit.index(field) #si es el probando2 indexalo

	if fatherIndex==-1:
		raise NotFoundFamily('father',padre)

	if motherIndex==-1:
		raise NotFoundFamily('mother',madre)

	if probandIndex==-1:
		raise NotFoundFamily('proband',caso)
		
	if proband2Index==-1:
		raise NotFoundFamily('proband2',caso2)

	return probandIndex,proband2Index,fatherIndex,motherIndex

def getIndexFamily_5(caso,padre,madre,parient,parient2,header_proband):
	header_probandSplit = header_proband.split('\t') #parte a header proband por los tabuladores
	
	fatherIndex=-1
	motherIndex=-1
	probandIndex=-1
	parientIndex=-1
	parient2Index=-1

	for field in header_probandSplit:
		# print (field)
		campo = field.split("\n") #parte el campo por los saltos de linea
		if campo[0] == str(padre):
			# print (field)
			fatherIndex = header_probandSplit.index(field) #Si el primer campo es padre, indexalo
		elif campo[0] == str(madre):
			# print (field)
			motherIndex = header_probandSplit.index(field) #Si es madre indexalo
		elif campo[0] == str(parient):
			# print (field)
			parientIndex = header_probandSplit.index(field) #Si es pariente1 indexalo
		elif campo[0] == str(parient2):
			# print (field)
			parient2Index = header_probandSplit.index(field) #Si es pariente2 indexalo
		elif campo[0] == str(caso):
			# print (field)
			probandIndex = header_probandSplit.index(field) #si es el probando indexalo

	if fatherIndex==-1:
		raise NotFoundFamily('father',padre)

	if motherIndex==-1:
		raise NotFoundFamily('mother',madre)

	if probandIndex==-1:
		raise NotFoundFamily('proband',caso)
	
	if parientIndex==-1:
		raise NotFoundFamily('parient',parient)
		
	if parient2Index==-1:
		raise NotFoundFamily('parient2',parient2)

	return probandIndex,fatherIndex,motherIndex,parientIndex,parient2Index
	
def getIndexFamily_5_2P(caso,caso2,padre,madre,parient,header_proband):
	header_probandSplit = header_proband.split('\t') #parte a header proband por los tabuladores
	
	fatherIndex=-1
	motherIndex=-1
	parientIndex=-1
	probandIndex=-1
	proband2Index=-1
	

	for field in header_probandSplit:
		# print (field)
		campo = field.split("\n") #parte el campo por los saltos de linea
		if campo[0] == str(padre):
			# print (field)
			fatherIndex = header_probandSplit.index(field) #Si el primer campo es padre, indexalo
		elif campo[0] == str(madre):
			# print (field)
			motherIndex = header_probandSplit.index(field) #Si es madre indexalo
		elif campo[0] == str(parient):
			# print (field)
			parientIndex = header_probandSplit.index(field) #Si es pariente indexalo
		elif campo[0] == str(caso):
			# print (field)
			probandIndex = header_probandSplit.index(field) #si es el probando indexalo
		elif campo[0] == str(caso2):
			# print (field)
			proband2Index = header_probandSplit.index(field) #si es el probando2 indexalo

	if fatherIndex==-1:
		raise NotFoundFamily('father',padre)

	if motherIndex==-1:
		raise NotFoundFamily('mother',madre)
		
	if parientIndex==-1:
		raise NotFoundFamily('parient',parient)

	if probandIndex==-1:
		raise NotFoundFamily('proband',caso)
		
	if proband2Index==-1:
		raise NotFoundFamily('proband2',caso2)

	return probandIndex,proband2Index,fatherIndex,motherIndex,parientIndex
	
def getIndexFamily_5_3P(caso,caso2,caso3,padre,madre,header_proband):
	header_probandSplit = header_proband.split('\t') #parte a header proband por los tabuladores
	
	fatherIndex=-1
	motherIndex=-1
	probandIndex=-1
	proband2Index=-1
	proband3Index=-1
	

	for field in header_probandSplit:
		# print (field)
		campo = field.split("\n") #parte el campo por los saltos de linea
		if campo[0] == str(padre):
			# print (field)
			fatherIndex = header_probandSplit.index(field) #Si el primer campo es padre, indexalo
		elif campo[0] == str(madre):
			# print (field)
			motherIndex = header_probandSplit.index(field) #Si es madre indexalo
		elif campo[0] == str(caso):
			# print (field)
			probandIndex = header_probandSplit.index(field) #si es el probando indexalo
		elif campo[0] == str(caso2):
			# print (field)
			proband2Index = header_probandSplit.index(field) #si es el probando2 indexalo
		elif campo[0] == str(caso3):
			# print (field)
			proband3Index = header_probandSplit.index(field) #si es el probando3 indexalo

	if fatherIndex==-1:
		raise NotFoundFamily('father',padre)

	if motherIndex==-1:
		raise NotFoundFamily('mother',madre)

	if probandIndex==-1:
		raise NotFoundFamily('proband',caso)
		
	if proband2Index==-1:
		raise NotFoundFamily('proband2',caso2)
	
	if proband3Index==-1:
		raise NotFoundFamily('proband2',caso3)

	return probandIndex,proband2Index,proband3Index,fatherIndex,motherIndex

class NotFoundFamily(Exception):
	def __init__(self,family,idfamily):
		message = 'Not found %s=%s in file '%(family,idfamily)
		super(NotFoundFamily, self).__init__(message)

class ErrorInputFileLatest(Exception):
	def __init__(self,file):
		message = 'Latest remodeled not found: (%s)'%(file)
		super(ErrorInputFileLatest, self).__init__(message)

class ErrorInputFileDuo(Exception):
	def __init__(self,file):
		message = 'Duo merge not found: (%s)'%(file)
		super(ErrorInputFileTrio, self).__init__(message)
		
class ErrorInputFileTrio(Exception):
	def __init__(self,file):
		message = 'Trio merge not found: (%s)'%(file)
		super(ErrorInputFileTrio, self).__init__(message)
		
class ErrorInputFileCuarteto(Exception):
	def __init__(self,file):
		message = 'Cuarteto merge not found: (%s)'%(file)
		super(ErrorInputFileCuarteto, self).__init__(message)

class ErrorInputFileQuinteto(Exception):
	def __init__(self,file):
		message = 'Quinteto merge not found: (%s)'%(file)
		super(ErrorInputFileQuinteto, self).__init__(message)
		
class ErrorFinalFile(Exception):
	def __init__(self,file):
		message = 'Final file not found: (%s)'%(file)
		super(ErrorFinalFile, self).__init__(message)

class ErrorTempFile(Exception):
	def __init__(self,file):
		message = 'Temp file not found: (%s)'%(file)
		super(ErrorTempFile, self).__init__(message)
