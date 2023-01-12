##SCRIPT PARA OBTENER EL ARCHIVO MULTISAMPLE VCF

shopt -s expand_aliases
alias gatk='/media/hdd-10T/homes/gbvicente/TFM_German/GATK-4.2.6.1/gatk'

#Se establece como primeros argumentos el nombre de las muestras sin la extensión, como segundo el genoma de referencia sin la extensión y como tercero el nombre del output

if [ "$#" -eq 5 ]
then
	echo "Preparando los archivos para el llamado de variantes"
	
	sampleid=$1
	sampleid2=$2
	sampleid3=$3
	genomefasta=$4
	output=$5
	
	if [ -e /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/Ref_Genome/${genomefasta}.fa.fai ] && [ -e /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/Ref_Genome/${genomefasta}.dict ]
	then
		echo "Ya existen los archivos .fai y .dict del genoma de referencia, se procede con el llamado de variantes"
	else
		#Se indexa el genoma de referencia
		
		samtools faidx /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/Ref_Genome/${genomefasta}.fa
		
		#Se obtiene el archivo .dict necesario para el llamado de variantes

		java -jar picard.jar CreateSequenceDictionary R=/media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/Ref_Genome/${genomefasta}.fa O=/media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/Ref_Genome/${genomefasta}.dict
	
	fi
	
	#Se ejecuta HaplotypeCaller para obtener el archivo single sample gVCF de la muestra

	echo "Comenzando el llamado de variantes..."
	
	gatk --java-options "-Xmx16g -XX:ParallelGCThreads=4" HaplotypeCaller -R /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/Ref_Genome/${genomefasta}.fa -I /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/BAMs/${genomefasta}/${sampleid}.bam -I /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/BAMs/${genomefasta}/${sampleid2}.bam -I /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/BAMs/${genomefasta}/${sampleid3}.bam -O /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_JC/${output}.vcf.gz
else
	echo "Error: los argumentos introducidos no son correctos. Introduzca en primer lugar el nombre de las muestras, en segundo lugar el nombre del genoma de referencia y en tercer lugar el nombre del output, todos sin extensión .fa y .bam"
fi


