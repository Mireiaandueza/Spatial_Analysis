#!/bin/bash

#SLIDE 1
spaceranger count --id=ST1-PER226-189 \
                   --transcriptome=./reference/refdata-gex-mm10-2020-A \
                   --fastqs=./data/ \
                   --sample=ST1-PER226-189 \
                   --image=./images/ST2-A1.tif \
                   --slide=V12N18-044 \
		   --probe-set=./reference/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv \
                   --area=A1 

spaceranger count --id=ST2-PER182-185 \
                   --transcriptome=./refdata-gex-mm10-2020-A \
                   --fastqs=./data/ \
                   --sample=ST2-PER182-185 \
                   --image=./images/ST2-B1.tif \
                   --slide=V12N18-044 \
		   --probe-set=./reference/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv \
                   --area=B1 

spaceranger count --id=ST3-PER158-157-155 \
                   --transcriptome=./reference/refdata-gex-mm10-2020-A \
                   --fastqs=./data/ \
                   --sample=ST3-PER158-157-155  \
                   --image=./images/ST2-C1.tif \
                   --slide=V12N18-044 \
                   --probe-set=./reference/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv \
                   --area=C1

spaceranger count --id=ST4-PER174-152 \
                   --transcriptome=./reference/refdata-gex-mm10-2020-A \
                   --fastqs=./data/ \
                   --sample=ST4-PER174-152 \
                   --image=./images/ST2-D1.tif \
                   --slide=V12N18-044 \
                   --probe-set=./reference/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv \
                   --area=D1


#SLIDE 2
spaceranger count --id=ST5-PER104-111 \
                   --transcriptome=./reference/refdata-gex-mm10-2020-A \
                   --fastqs=./data/ \
                   --sample=ST5-PER104-111 \
                   --image=./images/ST1-A1.tif \
                   --slide=V12N18-045 \
                   --probe-set=./reference/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv \
                   --area=A1

spaceranger count --id=ST6-PER124-106 \
                   --transcriptome=./reference/refdata-gex-mm10-2020-A \
                   --fastqs=./data/ \
                   --sample=ST6-PER124-106 \
                   --image=./images/ST1-B1.tif \
                   --slide=V12N18-045 \
                   --probe-set=./reference/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv \
                   --area=B1

spaceranger count --id=ST8-PER197-222 \
                   --transcriptome=./reference/refdata-gex-mm10-2020-A \
                   --fastqs=./data/ \
                   --sample=ST8-PER197-222 \
                   --image=./images/ST1-C1.tif \
                   --slide=V12N18-045 \
                   --probe-set=./reference/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv \
                   --area=C1

spaceranger count --id=ST7-PER099-112 \
                   --transcriptome=./reference/refdata-gex-mm10-2020-A \
                   --fastqs=./data/ \
                   --sample=ST7-PER099-112 \
                   --image=./images/ST1-D1.tif \
                   --slide=V12N18-045 \
                   --probe-set=./reference/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv \
                   --area=D1

#SLIDE 3
spaceranger count --id=ST13_PEO680-630 \
                   --transcriptome=./reference/refdata-gex-mm10-2020-A \
                   --fastqs=./data/ \
                   --sample=cb-ST13_PEO680-630 \
                   --image=./images/ST3_A1.tiff \
                   --slide=V12N18-046 \
                   --probe-set=./reference/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv \
                   --area=A1

spaceranger count --id=ST14_PEO425-661 \
                   --transcriptome=./reference/refdata-gex-mm10-2020-A \
                   --fastqs=./data/ \
                   --sample=cb-ST14_PEO425-661 \
                   --image=./images/ST3_B1.tiff \
                   --slide=V12N18-046 \
                   --probe-set=./reference/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv \
                   --area=B1

spaceranger count --id=ST15_PEO384-360-696 \
                   --transcriptome=./reference/refdata-gex-mm10-2020-A \
                   --fastqs=./data/ \
                   --sample=cb-ST15_PEO384-360-696  \
                   --image=./images/ST3_C1.tiff \
                   --slide=V12N18-046 \
                   --probe-set=./reference/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv \
                   --area=C1

spaceranger count --id=ST16_PEO372-298 \
                   --transcriptome=./reference/refdata-gex-mm10-2020-A \
                   --fastqs=./data/ \
                   --sample=cb-ST16_PEO372-298 \
                   --image=./images/ST3_D1.tiff \
                   --slide=V12N18-046 \
		   --probe-set=./reference/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv \
                   --area=D1



# SLIDE 4
spaceranger count --id=ST09_PEO688-PER135 \
                   --transcriptome=./reference/refdata-gex-mm10-2020-A \
                   --fastqs=./data/ \
                   --sample=cb-ST09_PEO688-PER135 \
                   --image=./images/ST4_A1.tiff \
                   --slide=V12N18-047 \
                   --probe-set=./reference/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv \
                   --area=A1

spaceranger count --id=ST10_PER156-160 \
                   --transcriptome=./reference/refdata-gex-mm10-2020-A \
                   --fastqs=./data/ \
                   --sample=cb-ST10_PER156-160 \
                   --image=./images/ST4_B1.tiff \
                   --slide=V12N18-047 \
                   --probe-set=./reference/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv \
                   --area=B1

spaceranger count --id=ST11_PER226-189 \
                   --transcriptome=./reference/refdata-gex-mm10-2020-A \
                   --fastqs=./data/ \
                   --sample=cb-ST11_PER226-189 \
                   --image=./images/ST4_C1.tiff \
                   --slide=V12N18-047 \
                   --probe-set=./reference/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv \
                   --area=C1

spaceranger count --id=ST12_PIR615-636 \
                   --transcriptome=./reference/refdata-gex-mm10-2020-A \
                   --fastqs=./data/ \
                   --sample=cb-ST12_PIR615-636 \
                   --image=./images/ST4_D1.tiff \
                   --slide=V12N18-047 \
                   --probe-set=./reference/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv \
                   --area=D1


spaceranger count --id=ST12_PIR615-636 \
                   --transcriptome=./reference/refdata-gex-mm10-2020-A \
                   --fastqs=/storage/scratch01/users/mandueza/pancreas_st/raw/data/ \
                   --sample=cb-ST12_PIR615-636 \
                   --image=./images/ST4_D1.tiff \
                   --slide=V12N18-047 \
                   --probe-set=./reference/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv \
                   --area=D1


