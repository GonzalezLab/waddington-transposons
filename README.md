<a name="readme-top"></a>

# Standing genetic variation explains trait evolution by Waddington genetic assimilation
Scripts for analyzing TE insertions the paper "Standing genetic variation explains trait evolution by Waddington genetic assimilation".

## Folders
- **0. Data processing**:  retrieving fastq files for each sample and run `fastp`
- **2. PoPoolationTE2**: scripts for preparing inputs, running `PoPoolationTE2` and process results
- **3. TEMP2**: scripts for preparing inputs, running `TEMP2` and process results
- **4. Common insertions**: find common TE insertions between `PoPoolationTE2` and `TEMP2` using `bedtools`
- **5. Post analysis**: obtaining `TEMP2` frequency for *de novo* and segregating TE insertions, perform Fisher's exact test in `R`

## Citation
Sabarís (2024) Standing genetic variation explains trait evolution by Waddington genetic assimilation. BioRxiv

## Contact

Project: [https://github.com/GonzalezLab/waddington-transposons](https://github.com/GonzalezLab/waddington-transposons)

<p align="right">(<a href="#readme-top">back to top</a>)</p>
