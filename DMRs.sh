while getopts "u:f:g:m:n:c:o:t:h:z" opt; do
  case $opt in
    u) UPD="$OPTARG" ;;     
    f) fatherHap1="$OPTARG" ;;   
    g) fatherHap2="$OPTARG" ;;  
    m) motherHap1="$OPTARG" ;;
    n) motherHap2="$OPTARG" ;;
    c) childHaps="$OPTARG" ;;
    o) out="$OPTARG" ;;
    t) thresh="$OPTARG" ;;
    z) code="$OPTARG" ;;
    h)                        # -h for help
       echo "Usage: $0 -u <upd.bed> -f <fatherHap1.bed> -g <fatherHap2.bed> -m <motherHap1.bed> -n <motherHap2.bed> -c <childHaps.bed> -t <threshold> -o <outDirName> -z <codeDir>"
       exit 0
       ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

F1=$out/temp/F1.bed #father Hap1 Bed
F2=$out/temp/F2.bed #father Hap2 Bed
M1=$out/temp/M1.bed #mother Hap1 Bed
M2=$out/temp/M2.bed #mother Hap2 Bed
Child=$out/temp/child.bed #child Bed

F12=$out/temp/F12.txt #DMRs between Father Hap1 and Father Hap2
F1C=$out/temp/F1C.txt #DMRs between Father Hap1 and Child
F2C=$out/temp/F2C.txt #DMRs between Father Hap2 and Child
M12=$out/temp/M12.txt #DMRs between Mother Hap1 and Mother Hap2
M1C=$out/temp/M1C.txt #DMRs between Mother Hap1 and Child
M2C=$out/temp/M2C.txt #DMRs between Mother Hap2 and Child

chroms=$(cut -f1 $UPD | sort -u)

mkdir $out
mkdir $out/temp

awk '{print $1, $2, $6, $7}' $fatherHap1 > $F1
awk '{print $1, $2, $6, $7}' $fatherHap2 > $F2
awk '{print $1, $2, $6, $7}' $motherHap1 > $M1
awk '{print $1, $2, $6, $7}' $motherHap2 > $M2
> $Child

for chrom in $chroms; do
    echo "Processing $chrom"
    awk -v chrom=$chrom '$1==chrom {print $1, $2, $6, $7}' $childHaps >> $Child
done

echo "running DMR"
echo "F1 and F2"
Rscript $code/DMRs.R  --hap1 $F1 --hap2 $F2 --output $F12 --threshold $thresh
echo "F1 and Child"
Rscript $code/DMRs.R  --hap1 $F1 --hap2 $Child --output $F1C --threshold $thresh
echo "F2 and Child"
Rscript $code/DMRs.R  --hap1 $F2 --hap2 $Child --output $F2C --threshold $thresh
echo "M1 and M2"
Rscript $code/DMRs.R  --hap1 $M1 --hap2 $M2 --output $M12 --threshold $thresh
echo "M1 and Child"
Rscript $code/DMRs.R  --hap1 $M1 --hap2 $Child --output $M1C --threshold $thresh
echo "M2 and Child"
Rscript $code/DMRs.R  --hap1 $M2 --hap2 $Child --output $M2C --threshold 0.3


awk '{print $1"\t"$2"\t"$3}' $F12 | sed '1d' | bedtools sort -i - > $out/temp/F12.bed
awk '{print $1"\t"$2"\t"$3}' $F1C | sed '1d' | bedtools sort -i - > $out/temp/F1C.bed
awk '{print $1"\t"$2"\t"$3}' $F2C | sed '1d' | bedtools sort -i - > $out/temp/F2C.bed
awk '{print $1"\t"$2"\t"$3}' $M12 | sed '1d' | bedtools sort -i - > $out/temp/M12.bed
awk '{print $1"\t"$2"\t"$3}' $M1C | sed '1d' | bedtools sort -i - > $out/temp/M1C.bed
awk '{print $1"\t"$2"\t"$3}' $M2C | sed '1d' | bedtools sort -i - > $out/temp/M2C.bed

cat $out/temp/F1C.bed $out/temp/F2C.bed | bedtools sort -i - | bedtools merge -i - | bedtools intersect -a - -b $out/temp/F12.bed > $out/temp/fatherROI.bed
cat $out/temp/M1C.bed $out/temp/M2C.bed | bedtools sort -i - | bedtools merge -i - | bedtools intersect -a - -b $out/temp/M12.bed > $out/temp/motherROI.bed

cat $out/temp/fatherROI.bed $out/temp/motherROI.bed | bedtools sort -i - | bedtools merge -i - | bedtools intersect -a - -b $UPD | bedtools sort -i - | bedtools merge -i - > $out/UPD_DMRs.bed



