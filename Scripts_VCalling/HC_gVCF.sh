##SCRIPT PARA OBTENER EL ARCHIVO GVCF DE CADA MUESTRA

shopt -s expand_aliases
alias gatk='/media/hdd-10T/homes/gbvicente/TFM_German/GATK-4.2.6.1/gatk'

#Se establece como primer argumento el nombre de la muestra sin la extensión y como segundo el genoma de referencia sin la extensión

if [ "$#" -eq 2 ]
then
	echo "Preparando los archivos para el llamado de variantes"
	
	sampleid=$1
	genomefasta=$2
	
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
	
	gatk --java-options "-Xmx8g -XX:+UseParallelGC -XX:ParallelGCThreads=2" HaplotypeCaller -R /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/Ref_Genome/${genomefasta}.fa -I /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/BAMs/${genomefasta}/${sampleid}.bam -O /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/gVCFs_2/${sampleid}.g.vcf.gz -ERC GVCF

else
	echo "Error: los argumentos introducidos no son correctos. Introduzca en primer lugar el nombre de la muestra y en segundo lugar el nombre del genoma de referencia, los dos sin extensión .fa y .bam"
fi


