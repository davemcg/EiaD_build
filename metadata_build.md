# Find Study
Usually you can find a GEO accession in a paper (GSE - and a bunch of digits). 

[https://ncbi.nlm.nih.gov/geo](https://ncbi.nlm.nih.gov/geo) -> enter in GEO ID (e.g. `GSE129479`)

# Get SRA xml
This contains a lot of sample related info

Click "SRA" in "Relations"

<img width="719" alt="Screen Shot 2022-10-19 at 14 25 05" src="https://user-images.githubusercontent.com/10225430/196773831-0a8569c3-ac06-4f8b-89fe-f6a70d865712.png">

Click "Send to", "File", and export as "Full XML"
<img width="816" alt="Screen Shot 2022-10-19 at 14 25 45" src="https://user-images.githubusercontent.com/10225430/196773951-887e7ca0-f0e0-4318-b473-c8b8f88ef729.png">

# Get SRA run info

The xml metadata above has just about everything - except the SRR run level ID necessary to get the raw fastq file

Click "SRA Metadata" on the GEO page (see first image). Then click Download -> Metadata

<img width="895" alt="Screen Shot 2022-10-19 at 14 28 27" src="https://user-images.githubusercontent.com/10225430/196774602-5f39f414-a5a6-449d-8b1a-b320d82ed6db.png">

# Go to R 
Example script: [scripts/metadata_builder/SRP191440.R](https://github.com/davemcg/EiaD_build/blob/recount3/scripts/metadata_custom/SRP191440.R)
