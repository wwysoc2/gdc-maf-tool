# GDC MAF Concatenation Tool

Easy GDC MAF file concatenation

The GDC MAF Concatenation Tool automatically downloads a set of open-access XXXX MAFs from the Genomic Data Commons and concatenates them.  Each subset of MAFs can be specified using a manifest downloaded from the GDC Data Portal (https://portal.gdc.cancer.gov) or using project IDs.  

Either a manifest (-m) or project (-p) must be specified. More than one project can be specified by passing a comma-delimited string with no spaces; i.e. TCGA-LUAD,TCGA-LUSC

Each file downloaded is checked automatically for integrity by comparing the observed (file) and expected (API) md5sum.  


usage: python gdc-maf-cat.py <-m MANIFEST or -p PROJECT_ID>

----GDC MAF Concatenation Tool v1.0----

optional arguments:

  -h, --help            show this help message and exit

  -m MANIFEST, --manifest MANIFEST
                        Specify MAF files with GDC Manifest

  -p PROJECT, --project PROJECT
                        Specify MAF files by project

  -o FILE_PREFIX, --output FILE_PREFIX
                        Designates a name for the output file
