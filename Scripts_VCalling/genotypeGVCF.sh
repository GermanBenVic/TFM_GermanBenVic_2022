##SCRIPT PARA REALIZAR EL VARIANT CALLING DE LOS MULTIGVCFS

shopt -s expand_aliases
alias gatk='/media/hdd-10T/homes/gbvicente/TFM_German/GATK-4.2.6.1/gatk'

##Llamando las variantes del archivo gVCF combinado

if [ "$#" -eq 3 ]
then
	echo "Preparando los archivos para el variant calling"

	multiGVCF=$1
	genomefasta=$2
	output=$3

	gatk --java-options "-Xmx4g -XX:ParallelGCThreads=2" GenotypeGVCFs \
		-R /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/Ref_Genome/${genomefasta}.fa \
		-V /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/MultiGVCF/${multiGVCF}.g.vcf.gz \
		-O /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/${output}.vcf.gz
else
	echo "Error: los argumentos introducidos no son correctos. Introduzca en primer lugar el nombre del gVCF combinado, en segundo lugar el nombre del genoma de referencia y en tercer lugar el nombre que tendrá el output, todos sin la extensión"
fi
