#!/bin/bash

# replace *...* with respective value

cellranger count --id=*sampleName* \
        --transcriptome=*path to reference* \
        --fastqs=*path to folder with all fastqs* \
        --sample=*sampleName with all variations of flowcell* \
        --include-introns true \
	--disable-ui