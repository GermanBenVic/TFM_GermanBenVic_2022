##SCRIPT PARA OBTENER UN UNICO ARCHIVO MULTISAMPLE GVCF

shopt -s expand_aliases
alias gatk='/media/hdd-10T/homes/gbvicente/TFM_German/GATK-4.2.6.1/gatk'

##Combinando los archivos gVCFs en un multisample gVCF, se puede cambiar el número de muestras a combinar dependiendo del número de gVCFs

if [ "$#" -eq 5 ]
then
	echo "Combinando los archivos gVCFs de la cohorte"
	sampleid1=$1
	sampleid2=$2
	sampleid3=$3
	genomefasta=$4
	output=$5

	gatk CombineGVCFs \
		-R /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/Ref_Genome/${genomefasta}.fa \
		--variant /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/gVCFs/${sampleid1}.g.vcf.gz \
		--variant /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/gVCFs/${sampleid2}.g.vcf.gz \
		--variant /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/gVCFs/${sampleid3}.g.vcf.gz \
		-O /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/MultiGVCF/${output}.g.vcf.gz
else
	echo "Error: los argumentos introducidos no son correctos. Introduzca en primer lugar los nombres de los archivos gVCFs, en segundo lugar el nombre del genoma de referencia y en tercer lugar el nombre que desee que tenga el output, todos sin extensión"
fi
