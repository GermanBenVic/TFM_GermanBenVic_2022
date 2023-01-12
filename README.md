# TFM_GermanBenVic_2022
In this repository you can find the code used in the development of my master thesis and also a summarised version of the inputs and outputs of NIMGree.

# This repository is divided into following directories:
```shell
├── NIMGree
│   ├── Archivos_ejemplo
│   │   ├── ejemplo_input_latest_P1.tsv
│   │   ├── ejemplo_input_latest_P2.tsv
│   │   ├── ejemplo_input_multivcf.vcf
│   │   └── ejemplo_output_NIMGree.tsv
│   ├── common.py
│   └── NIMGree.py 
└── Scripts_VCalling
    ├── combinegVCF.sh
    ├── genotypeGVCF.sh
    ├── HC_gVCF.sh
    ├── HC_JointCalling.sh
    ├── vqsr_JC.sh
    └── vqsr.sh
```

# The following files are located in the NIMGree directory:
  - common.py : in it are the functions that NIMGree.py imports in order to execute its code.
  - NIMGree. py: it contains the code divided into classes to make the script more easily executable. This code performs the crosswalk between the multiVCF     file and the annotation file to add the segregation information. It is divided in the following classes:
      - duoDragen: for duo exome analysis of a family member and a proband.
      - duoDragen_2P: for duo exome analysis of two probands.
      - trioDragen: for trio exome analysis of a proband and its progenitors.
      - cuartetoDragen: for quad exome analysis of one proband and three relatives
      - cuartetoDragen_2P: for quad exome analysis of two probands and its progenitors.
      - quintetoDragen: for exome analysis of one proband and four relatives.
      - quintetoDragen_2P: for exome analysis of two probands and three relatives.
      - quintetoDragen_3P: for exome analysis of three probands and its progenitors.
    
   - Also you can find a summarised example of the inputs and output of NIMGree.py:
      - ejemplo_input_latest_P1.tsv
      - ejemplo_input_latest_P2.tsv
      - ejemplo_input_multivcf.vcf
      - ejemplo_output_NIMGree.tsv
     
  This code should be executed in the following way:
    
   ```shell
   variable = class("WD", "path to VCF","path to annotation file")
   variable.run ("sample1", "sample2" ..., "sex1", "sex2" ...)
   ```
   
# The following files are located in the Scripts_VCalling directory:
  - combinegVCF.sh : script with code from GATK that combine the GVCF files.
  - genotypeGVCF.sh : script with code from GATK that performs the variant calling of the multisampleGVCF.
  - HC_gVCF.sh : script with code from GATK that processes the bam files to obtain the GVCF files.
  - HC_JointCalling.sh : script with code from GATK that processes the bam files to obtain the multiVCF files.
  - vqsr_JC.sh : script that performs variant filtering on files generated with JC.
  - vqsr.sh : script that performs variant filtering on files generated with HC.
 
 
If you have any doubts or questions please contact me at german.ben.vic@gmail.com
