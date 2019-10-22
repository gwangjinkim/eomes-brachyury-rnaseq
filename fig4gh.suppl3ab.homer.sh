####################################
# analysis for Fig 4 g+h using HOMER
####################################


wt_atac_peaks="/media/josephus/Elements/galaxy_atac_data/atac-wt-merge-86.bed"
dko_atac_peaks="/media/josephus/Elements/galaxy_atac_data/atac-dko-merge-92.bed"

mkdir -p ${wt_atac_out}
mkdir -p ${dko_atac_out}

findMotifsGenome.pl ${wt_atac_peaks} mm10 ${wt_atac_out}/ -bg ${dko_atac_peaks} -size 200 -p 6 && \
findMotifsGenome.pl ${dko_atac_peaks} mm10 ${dko_atac_out} -bg ${wt_atac_peaks} -size 200 -p 6

rsync -av /media/josephus/Elements/ChIPresults/ /media/josephus/24AA554AAA551A1E/ChIPresults/


################################################################################

####################################
# analysis for Suppl Fig 3 a+b
####################################

path_to_tosicj190315="/media/josephus/My Book/1TB_backup/tosicj190315"
out="/media/josephus/Elements/ChIPresults/homer/ChIPseq2019"
mkdir -p ${out}

eo_final_peak_fpath="${path_to_tosicj190315}/Peaks_ChIP_EoGFP_pileup120.bed"
br_final_peak_fpath="${path_to_tosicj190315}/Peaks_ChIP_BraGFP_pileup40.bed"

out_eo_bg_br="${out_eo}/eo_final_bg_br_final"
out_br_bg_eo="${out_br}/br_final_bg_eo_final"

mkdir -p ${out_eo_bg_br}
mkdir -p ${out_br_bg_eo}

findMotifsGenome.pl "${eo_final_peak_fpath}" mm10 "${out_eo_bg_br}"/ -bg "${br_final_peak_fpath}" -size 200 -p 6 && \
findMotifsGenome.pl "${br_final_peak_fpath}" mm10 "${out_br_bg_eo}"/ -bg "${eo_final_peak_fpath}" -size 200 -p 6



#####################################
# count of intervals
#####################################

wc -l "/media/josephus/Elements/galaxy_atac_data/atac-wt-merge-86.bed"
## 56089 /media/josephus/Elements/galaxy_atac_data/atac-wt-merge-86.bed
wc -l "/media/josephus/Elements/galaxy_atac_data/atac-dko-merge-92.bed"
## 53224 /media/josephus/Elements/galaxy_atac_data/atac-dko-merge-92.bed
wc -l "/media/josephus/My Book/1TB_backup/tosicj190315/Peaks_ChIP_EoGFP_pileup120.bed"
## 30144 /media/josephus/My Book/1TB_backup/tosicj190315/Peaks_ChIP_EoGFP_pileup120.bed
wc -l "/media/josephus/My Book/1TB_backup/tosicj190315/Peaks_ChIP_BraGFP_pileup40.bed"
## 27598 /media/josephus/My Book/1TB_backup/tosicj190315/Peaks_ChIP_BraGFP_pileup40.bed






