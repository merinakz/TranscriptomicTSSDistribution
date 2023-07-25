# Transcriptomic TSS Distribution

This code analyses the output file from TSSPredator and distributes TSS coordinates according to their enrichment scores. 

This code will only work if you have 3 different conditions. 

Case 1: Enriched for every condition. 
Case 2: Enriched for NDC only. 
Case 3: Enriched for condition 1 only
Case 4: Enriched for condition 2 only
Case 5: Enriched for both conditions 1 & 2

Transcriptome annotation is performed based on the start and end position of each gene. 5'UTR and 3'UTR were assigned 150 nt from each side. 
