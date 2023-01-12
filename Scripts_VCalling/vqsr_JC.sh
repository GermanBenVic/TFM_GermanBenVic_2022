##SCRIPT PARA OBTENER EL ARCHIVO MULTISAMPLE VCF FILTRADO

shopt -s expand_aliases
alias gatk='/media/hdd-10T/homes/gbvicente/TFM_German/GATK-4.2.6.1/gatk'

#Se establecen argumentos para llevar a cabo el programa

##IMPORTANTE: TEN EN CUENTA QUE HAY QUE CAMBIAR LOS ARCHIVOS DE FILTRADO A LA VERSIÓN DEL GENOMA CORRESPONDIENTE

if [ "$#" -eq 2 ]
then
	echo "Preparando los archivos para el llamado de variantes"

	sampleid=$1
	genomefasta=$2
	
	mkdir -p /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_recal/${sampleid}_JC

	##CREANDO EL MODELO DE RECALIBRADO PARA SNPs
	
	echo "Creando el modelo para el recalibrado de SNVs"
	
	gatk --java-options "-Xmx12g -Xms13g" VariantRecalibrator \
		-R /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/Ref_Genome/${genomefasta}.fa\
		-V /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_JC/${sampleid}.vcf.gz \
		-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
		-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
		-mode SNP \
		--max-gaussians 6 \
		--resource:hapmap,known=false,training=true,truth=true,prior=15 /media/hdd-10T/homes/gbvicente/TFM_German/vqsr/hapmap_3.3.hg19.vcf.gz \
		--resource:omni,known=false,training=true,truth=true,prior=12 /media/hdd-10T/homes/gbvicente/TFM_German/vqsr/1000G_omni2.5.hg19.vcf.gz \
		--resource:1000G,known=false,training=true,truth=false,prior=10 /media/hdd-10T/homes/gbvicente/TFM_German/vqsr/1000G_phase1.snps.high_confidence.hg19.vcf.gz \
		--resource:dbsnp,known=true,training=false,truth=false,prior=7 /media/hdd-10T/homes/gbvicente/TFM_German/vqsr/Homo_sapiens_assembly19.dbsnp138.vcf \
		-O /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_recal/${sampleid}_JC/${sampleid}_JC.SNP.recal \
		--output-model /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_recal/${sampleid}_JC/${sampleid}_JC.SNP.model\
		--tranches-file /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_recal/${sampleid}_JC/${sampleid}_JC.SNP.tranches
	
	##CREANDO EL MODELO DE RECALIBRADO PARA INDELs
	
	echo "Creando el modelo para el recalibrado de INDELs"
	
	gatk --java-options "-Xmx12g -Xms12g" VariantRecalibrator \
		-R /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/Ref_Genome/${genomefasta}.fa\
		-V /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_JC/${sampleid}.vcf.gz \
		-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
		-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
		-mode INDEL \
		--max-gaussians 4 \
		--resource:mills,known=false,training=true,truth=true,prior=12 /media/hdd-10T/homes/gbvicente/TFM_German/vqsr/Mills_and_1000G_gold_standard.indels.hg19.vcf.gz \
		--resource:axiomPoly,known=false,training=true,truth=false,prior=10 /media/hdd-10T/homes/gbvicente/TFM_German/vqsr/Axiom_Exome_Plus.genotypes.all_populations.poly.hg19.vcf.gz \
		--resource:dbsnp,known=true,training=false,truth=false,prior=7 /media/hdd-10T/homes/gbvicente/TFM_German/vqsr/Homo_sapiens_assembly19.dbsnp138.vcf \
		-O /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_recal/${sampleid}_JC/${sampleid}_JC.INDEL.recal \
		--output-model /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_recal/${sampleid}_JC/${sampleid}_JC.INDEL.model\
		--tranches-file /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_recal/${sampleid}_JC/${sampleid}_JC.INDEL.tranches
	
	##APLICANDO EL MODELO DE RECALIBRADO PARA SNPs
	
	echo "Aplicando el modelo para el recalibrado de SNPs"
	
	mkdir -p /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_filtered/${sampleid}_JC

	gatk ApplyVQSR \
  		-R /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/Ref_Genome/${genomefasta}.fa \
   		-V /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_JC/${sampleid}.vcf.gz \
   		-O /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_filtered/${sampleid}_JC/${sampleid}_JC.SNP-filtered.vcf.gz \
   		--truth-sensitivity-filter-level 99.9 \
   		--tranches-file /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_recal/${sampleid}_JC/${sampleid}_JC.SNP.tranches \
   		--recal-file /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_recal/${sampleid}_JC/${sampleid}_JC.SNP.recal \
   		-mode SNP
	
	##APLICANDO EL MODELO DE RECALIBRADO PARA INDELs
	
	echo "Aplicando el modelo para el recalibrado de INDELs"

	gatk ApplyVQSR \
  		-R /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/Ref_Genome/${genomefasta}.fa \
   		-V /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_JC/${sampleid}.vcf.gz \
   		-O /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_filtered/${sampleid}_JC/${sampleid}_JC.INDEL-filtered.vcf.gz \
   		--truth-sensitivity-filter-level 99.9 \
   		--tranches-file /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_recal/${sampleid}_JC/${sampleid}_JC.INDEL.tranches \
   		--recal-file /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_recal/${sampleid}_JC/${sampleid}_JC.INDEL.recal \
   		-mode INDEL
   	
	##FILTRANDO SOLO LAS VARIANTES QUE PASAN EL FILTRO
	#SNP

	gatk --java-options '-Xmx12g' SelectVariants \
		-R /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/Ref_Genome/${genomefasta}.fa \
		-V /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_filtered/${sampleid}_JC/${sampleid}_JC.SNP-filtered.vcf.gz \
		-O /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_filtered/${sampleid}_JC/${sampleid}_JC.SNP-filtered.PASS.vcf.gz \
		-L /media/hdd-10T/homes/gbvicente/TFM_German/vqsr/Twist_Exome_RefSeq_targets_${genomefasta}.bed \
		--exclude-filtered \
		--select-type-to-exclude INDEL
	
	#INDELS

	gatk --java-options '-Xmx12g' SelectVariants \
		-R /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/Ref_Genome/${genomefasta}.fa \
		-V /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_filtered/${sampleid}_JC/${sampleid}_JC.INDEL-filtered.vcf.gz \
		-O /media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_filtered/${sampleid}_JC/${sampleid}_JC.INDEL-filtered.PASS.vcf.gz \
		-L /media/hdd-10T/homes/gbvicente/TFM_German/vqsr/Twist_Exome_RefSeq_targets_${genomefasta}.bed \
		--exclude-filtered \
		--select-type-to-exclude SNP



   	##HACIENDO EL MERGE DE LOS VCFs FILTRADOS
   	
   	echo "Uniendo los VCFs filtrados"

	java -jar picard.jar SortVcf \
		I=/media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_filtered/${sampleid}_JC/${sampleid}_JC.SNP-filtered.PASS.vcf.gz \
		I=/media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_filtered/${sampleid}_JC/${sampleid}_JC.INDEL-filtered.PASS.vcf.gz \
		O=/media/hdd-10T/homes/gbvicente/TFM_German/TFM_NIM/vcfs/vcfs_filtered/${sampleid}_JC/${sampleid}_JC.filtered.vcf.gz

else

	echo "Los argumentos introducidos no son correctos en primer lugar introduzca el nombre del archivo vcf, en segundo lugar el genoma de referencia"
fi
