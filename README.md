**Poplar Microbiome Project**

## **Bioinformatics project for the Big Data in Agriculture Internship**

Student: Morgan Miller

Mentor: Ricardo Alcala Briseno

PI: Jared LeBoldus

## **Objective**

Use ITS metabarcoding to characterize the fungal wood microbiome of diseased (infected with Sphaerulina musiva) and healthy Populus trichocarpa trees. Comparison of microbiome makeup and diversity between diseased and healthy trees could illuminate how fungal endophytes contribute to phytopathogen resistance.

See (**Busby LeBoldus) NSF RAPID Proposal.pdf** in the **literature** directory for more information

## **Background**

In 2018, an outbreak of the phytopathogen Sphaerulina musiva was reported in a Populus trichocarpa (Black Cottonwood) plantation in northeastern Oregon. Introduction of S. musiva to the Pacific Northwest threatens both regional hybrid-poplar plantations and native riparian ecosystems where cottonwoods are instrumental in floodplain succession and ecological functioning.

Recent studies have demonstrated that fungal endophytes (fungal microbes living within plant tissues) can modify the disease resistance and susceptibility of their plant hosts. Despite a growing interest in fungal endophyte research, the specifics of these interactions remain poorly understood. Fungal endophytes are seemingly capable of modifying plant-disease in a variety of ways including:

-   Directly antagonizing phytopathogens

-   Indirectly competing with phytopathogens

-   Activating host-resistance

-   Suppressing host-resistance.

The diversity and complexity of fungal endophyte communities makes characterizing these interactions challenging. Quantitative molecular approaches such as ITS metabarcoding will be critical techniques for developing our understanding of fungal endophyte communities and their effect on disease resistance and susceptibility.

## **Acknowledging Biases Associated with ITS Amplification for Determining Species Composition.**

In a perfect world, every ITS region of the DNA fragments present in the original samples would be replicated evenly and proportionately before sequencing, giving us a dataset that represents the biological reality of the samples. Unfortunately, biases and errors introduced during the PCR process can skew sequencing results by overrepresenting or underrepresenting certain taxa or by producing erroneous DNA fragments. There are many factors that contribute to PCR bias and error within the ITS metabarcoding process:

-   Taxonomic Bias: Certain combinations of ITS primers may favor certain taxonomic groups over others leading to disproportionate amplification and representation.

-   Length Bias: ITS region length will vary between taxonomic groups. Because shorter sequences are more easily amplified than longer sequences when PCR is performed on sequences with varied length simultaneously, inherent diversity in ITS region length could contribute to disproportionate representation of certain taxonomic groups.

-   Mismatch Bias: Certain combinations of ITS primers respond differently to PCR conditions that would allow or disallow \'mismatches\' where the primer sequence differed from the binding region by a few base pairs. This means that under strict PCR conditions, primer \'preference\' might be more severe.

Clearly, the primers used for amplification considerably contribute to the potential for PCR bias and error. Careful consideration of primer strengths and challenges is critical for minimizing amplification error prior to sequencing. For this study, ITS3_KY01 (forward) and ITS4 (reverse) were used to sequence the ITS2 region. . . **WIPWIPWIPWIP**

## **ASVs and OTUs**

Sequencing technology introduces error. In order for any analysis to meaningfully represent the biological reality of a microbiome, sequencing error must be addressed. Otherwise, every sequence with an incorrect base call would be treated as if it represented a distinct organism in the community. A standardized technique for denoising sequence data is the construction of Operational Taxonomic Units (OTUs). OTUs are created by grouping sequences together based on similarity. In essence OTUs resolve sequence error by \'blurring\' similar sequences together so that minimal error won\'t lead to erroneously identified taxa. More recently, algorithms have been developed that can produce Amplicon Sequence Variants (ASVs) by inferring exact sequences and separating them from those likely caused by error with high confidence. . . **WIPWIPWIPWIPW**

For this study, ASVs are particularly advantageous. Because fungal wood microbiomes are relatively under-studied, ASVs allow for identification of fungal endophytes that might not be documented in a reference database. Open Reference Clustering techniques, and the OTUs they would produce, would likely fail to identify undocumented members of the microbial community due to reference bias\*. Additionally, an ASV approach allows these results to be meaningfully compared to other ongoing and future research.

\*Of course ASVs do not eliminate reference bias altogether as they still need to be compared to a database for taxonomic assignment.

### **Video Resource**

This page contains an excellent video comparing ASVs and OTUs. Produced by Michael Weinstein at Zymo Research.[ https://www.zymoresearch.com/blogs/blog/microbiome-informatics-otu-vs-asv](https://www.zymoresearch.com/blogs/blog/microbiome-informatics-otu-vs-asv).

## **Data Cleaning and Denoising Procedure**

### **Raw Data**

528 wood samples were collected from the infected Populus trichocarpa plantation. Samples were processed and ITS PCR performed to amplify the region of interest. Samples were multiplexed and combined before being sequenced (Illumina). . . **WIPWIPWIPW**

### **Remove Primers With Cutadapt**

### **Denoise With DADA2**

\
\
