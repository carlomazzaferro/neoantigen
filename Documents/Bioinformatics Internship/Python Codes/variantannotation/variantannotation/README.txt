Required library: myvariant. The other ones should be installed natively (all of them are in the CCBB server)

Download and open main.py. Change the paths so that variantAnnotation can added to your python path (with the command
python.sys.append(). Change then file paths to the vcf and csv file (if any).

Fo method A, a processed csv file from annovar is required, and it will provide the user with integrated data from
annovar and myvariant.info. For methods B, the VCF file is enough, and the functions used will create a list of
dictionaries with information retrieved from myvariant.info query service.

ANNOVAR standard script:
[sudo] perl /database/annovar/table_annovar.pl (/filepath) /database/annovar/humandb/ -buildver hg19 -out (/output filepath and name) -remove -protocol knownGene,tfbsConsSites,cytoBand,targetScanS,genomicSuperDups,gwasCatalog,esp6500siv2_all,1000g2015aug_all,snp138,ljb26_all,cg46,cg69,popfreq_all,clinvar_20140929,cosmic70,nci60 -operation g,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -csvout



