# *Oryza sativa* Nipponbare 倍数因子测试


[TOC levels=1-3]: # ""

- [*Oryza sativa* Nipponbare 倍数因子测试](#oryza-sativa-nipponbare-倍数因子测试)
  - [基本信息](#基本信息)
  - [Symlink](#symlink)
  - [Trim and cutoff](#trim-and-cutoff)
  - [`kat hist` and `kat gcp`](#kat-hist-and-kat-gcp)
  - [Mapping](#mapping)
  - [Depth](#depth)
  - [Merge all results](#merge-all-results)
  - [Remove intermediate files](#remove-intermediate-files)

## 基本信息

`cutoff = 倍数因子 * 覆盖深度`

+ 因子值 0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64
+ 测试文件 [SRR545231](https://www.ncbi.nlm.nih.gov/sra/SRX179254)
+ 覆盖度 17,220,721,594 / 373,870,564 = 46

```shell script
BASES=$(
    faops n50 -N 0 -S \
        ~/data/plastid/ena/SRR545231_1.fastq.gz \
        ~/data/plastid/ena/SRR545231_2.fastq.gz |
    cut -f 2
)

GENOME=$(
    faops size ~/data/plastid/genome/nip/genome.fa |
        cut -f 2 |
        perl -nl -e '
            $sum += $_;
            END {print $sum}
        '
)

echo ${BASES}
# 17220721594

echo ${GENOME}
# 373870564

bc <<< "${BASES} / ${GENOME}"

```

## Symlink

```shell script
mkdir -p ~/data/plastid/evaluation/nip
cd ~/data/plastid/evaluation/nip

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR545231_${NAME}
    
    mkdir -p ${BASE_NAME}/1_genome
    pushd ${BASE_NAME}/1_genome
    
    ln -fs ../../../../genome/nip/genome.fa genome.fa
    popd
    
    mkdir -p ${BASE_NAME}/2_illumina
    pushd ${BASE_NAME}/2_illumina
    
    ln -fs ../../../../ena/SRR545231_1.fastq.gz R1.fq.gz
    ln -fs ../../../../ena/SRR545231_2.fastq.gz R2.fq.gz
    popd

done

```

## Trim and cutoff

```shell script
cd ~/data/plastid/evaluation/nip

# 倍数因子::cutoff
ARRAY=(
    '0::0'
    '0.25::11'
    '0.5::23'
    '1::46'
    '2::92'
    '4::184'
    '8::368'
    '16::736'
    '32::1472'
    '64::2944'
)

for item in "${ARRAY[@]}" ; do
    NAME="${item%%::*}"
    CUTOFF="${item##*::}"

    echo 1>&2 "==> ${item}"
    
    BASE_NAME=SRR545231_${NAME}
    pushd ${BASE_NAME}
    
    rm *.sh
        
    if [[ "${NAME}" == "0" ]]; then
        anchr template \
            --genome 373870564 \
            --parallel 24 \
            --xmx 80g \
            \
            --fastqc \
            --insertsize \
            --kat \
            \
            --trim "--dedupe" \
            --qual "25" \
            --len "60" \
            --filter "adapter artifact"

        bsub -q mpi -n 24 -J "${BASE_NAME}" "
            bash 2_fastqc.sh
            bash 2_insert_size.sh
            bash 2_kat.sh
            bash 2_trim.sh
            bash 9_stat_reads.sh
        "        
    else
        anchr template \
            --genome 373870564 \
            --parallel 24 \
            --xmx 80g \
            \
            --fastqc \
            --insertsize \
            --kat \
            \
            --trim "--dedupe --cutoff ${CUTOFF} --cutk 31" \
            --qual "25" \
            --len "60" \
            --filter "adapter artifact"

        bsub -q mpi -n 24 -J "${BASE_NAME}" "
            bash 2_trim.sh
            bash 9_stat_reads.sh
        "
    fi

    popd

done

```

## `kat hist` and `kat gcp`

```shell script
cd ~/data/plastid/evaluation/nip

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR545231_${NAME}
    
    mkdir -p ${BASE_NAME}/kat
    pushd ${BASE_NAME}/kat
    
    bsub -q mpi -n 24 -J "${BASE_NAME}-kat" "
        kat hist \
            -t 24 -m 31 \
            ../2_illumina/trim/Q25L60/R1.fq.gz \
            ../2_illumina/trim/Q25L60/R2.fq.gz \
            ../2_illumina/trim/Q25L60/Rs.fq.gz \
            -o R-hist-31

        kat gcp \
            -t 24 -m 31 \
            ../2_illumina/trim/Q25L60/R1.fq.gz \
            ../2_illumina/trim/Q25L60/R2.fq.gz \
            ../2_illumina/trim/Q25L60/Rs.fq.gz \
            -o R-gcp-31
    "
    
    popd

done

```

## Mapping

```shell script
cd ~/data/plastid/evaluation/nip

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR545231_${NAME}
    
    mkdir -p ${BASE_NAME}/mapping
    pushd ${BASE_NAME}/mapping
    
    bsub -q mpi -n 24 -J "${BASE_NAME}-mapping" '

        # Pipe all reads together as we do not need mate info
        gzip -dcf \
            ../2_illumina/trim/Q25L60/R1.fq.gz \
            ../2_illumina/trim/Q25L60/R2.fq.gz \
            ../2_illumina/trim/Q25L60/Rs.fq.gz |
            bwa mem -M -t 20 \
                ../../../../genome/nip/genome.fa \
                /dev/stdin |
            pigz -p 4 \
            > R.sam.gz
        
        picard CleanSam \
            --INPUT R.sam.gz \
            --OUTPUT R.clean.bam \
            --VALIDATION_STRINGENCY LENIENT
    
        picard SortSam \
            --INPUT R.clean.bam \
            --OUTPUT R.sort.bam \
            --SORT_ORDER coordinate \
            --VALIDATION_STRINGENCY LENIENT
    
        picard BuildBamIndex \
            --INPUT R.sort.bam \
            --VALIDATION_STRINGENCY LENIENT
            
        find . -name "R.sam.gz" | xargs rm
        find . -name "R.clean.bam" | xargs rm
    '
    
    popd

done

```


## Depth

* Depth via `mosdepth`

```shell script
cd ~/data/plastid/evaluation/nip

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR545231_${NAME}
    
    echo 1>&2 "==> ${BASE_NAME}"
    
    mkdir -p ${BASE_NAME}/depth
    pushd ${BASE_NAME}/depth

    mosdepth R ../mapping/R.sort.bam
    
    popd

done

```

* Covered regions via `spanr`

```shell script
cd ~/data/plastid/evaluation/nip

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR545231_${NAME}
    
    echo 1>&2 "==> ${BASE_NAME}"
    
    mkdir -p ${BASE_NAME}/depth
    pushd ${BASE_NAME}/depth

    gzip -dcf R.per-base.bed.gz |
        perl -nla -F"\t" -e '
            $F[3] == 0 and next;
            $start = $F[1] + 1;
            if ($start == $F[2]) {
                print qq($F[0]:$start);
            }
            else {
                print qq($F[0]:$start-$F[2]);
            }
        ' |
        spanr cover stdin -o covered.yml
        
    spanr stat ../../../../genome/nip/chr.sizes covered.yml -o stdout |
        grep -v "^all" |
        sed 's/^chr/chrom/' |
        sed 's/,size/,covLength/' |
        sed 's/,coverage/,covRate/' |
        sed 's/,/\t/g' \
        > coverage.tsv
    
    popd

done

```

* Combine results

```shell script
cd ~/data/plastid/evaluation/nip

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR545231_${NAME}
    
    echo 1>&2 "==> ${BASE_NAME}"
    
    mkdir -p ${BASE_NAME}/depth
    pushd ${BASE_NAME}/depth

    cat coverage.tsv |
        tsv-join -H --filter-file R.mosdepth.summary.txt \
            --key-fields chrom --append-fields 3-6 \
        > join.tsv

    cat join.tsv |
        grep -v "^Mt" |
        grep -v "^Pt" |
        tsv-summarize -H --sum chrLength,covLength,bases --min min --max max |
        sed '1d' |
        perl -e '
            my $line = <>;
            chomp $line;
            my ($chrLength, $covLength, $bases, $min, $max, ) = split qq(\t), $line;
            my $covRate = sprintf qq(%.4f), $covLength / $chrLength;
            my $mean = sprintf qq(%.2f), $bases / $chrLength;
            print join qq(\t), (
                "Nc", $chrLength, $covLength, $covRate, $bases, $mean, $min, $max, 
            ); 
            print qq(\n);
        ' |
        (cat join.tsv | sed '2,13d' && cat) \
        > combine.tsv

    popd

done

```

## Merge all results

```shell script
cd ~/data/plastid/evaluation/nip

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR545231_${NAME}
    
    echo 1>&2 "==> ${BASE_NAME}"
    
    mkdir -p ${BASE_NAME}/depth
    pushd ${BASE_NAME}/depth > /dev/null

    echo -e "Fold\tchrom\n${NAME}\tNc\n${NAME}\tMt\n${NAME}\tPt" |
        tsv-join -H --filter-file combine.tsv --key-fields chrom --append-fields 2-8

    popd > /dev/null

done |
    tsv-uniq \
    > SRR545231_folds.tsv

for PART in Nc Mt Pt; do
    cat SRR545231_folds.tsv |
        tsv-filter -H --str-eq chrom:${PART} |
        mlr --itsv --omd cat
    echo
    echo
done

```

| Fold | chrom | chrLength | covLength | covRate | bases       | mean  | min | max   |
|:-----|:------|:----------|:----------|:--------|:------------|:------|:----|:------|
| 0    | Nc    | 373245519 | 368880800 | 0.9883  | 13576809466 | 36.38 | 0   | 26417 |
| 0.25 | Nc    | 373245519 | 366981283 | 0.9832  | 13553005337 | 36.31 | 0   | 26377 |
| 0.5  | Nc    | 373245519 | 356497402 | 0.9551  | 13218187244 | 35.41 | 0   | 26367 |
| 1    | Nc    | 373245519 | 174643554 | 0.4679  | 4695599969  | 12.58 | 0   | 26299 |
| 2    | Nc    | 373245519 | 121463116 | 0.3254  | 3331082797  | 8.92  | 0   | 26323 |
| 4    | Nc    | 373245519 | 99399478  | 0.2663  | 2712470322  | 7.27  | 0   | 26400 |
| 8    | Nc    | 373245519 | 83231703  | 0.2230  | 2266009858  | 6.07  | 0   | 26287 |
| 16   | Nc    | 373245519 | 68138024  | 0.1826  | 1859787055  | 4.98  | 0   | 26399 |
| 32   | Nc    | 373245519 | 51746684  | 0.1386  | 1366345135  | 3.66  | 0   | 26215 |
| 64   | Nc    | 373245519 | 34772643  | 0.0932  | 868522731   | 2.33  | 0   | 26060 |


| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:----|
| 0    | Mt    | 490520    | 488208    | 0.9953  | 62199415 | 126.80 | 0   | 519 |
| 0.25 | Mt    | 490520    | 488202    | 0.9953  | 62208977 | 126.82 | 0   | 522 |
| 0.5  | Mt    | 490520    | 488196    | 0.9953  | 62229153 | 126.86 | 0   | 530 |
| 1    | Mt    | 490520    | 488120    | 0.9951  | 62192323 | 126.79 | 0   | 530 |
| 2    | Mt    | 490520    | 488144    | 0.9952  | 62179239 | 126.76 | 0   | 503 |
| 4    | Mt    | 490520    | 434859    | 0.8865  | 46184786 | 94.15  | 0   | 525 |
| 8    | Mt    | 490520    | 107245    | 0.2186  | 8023602  | 16.36  | 0   | 528 |
| 16   | Mt    | 490520    | 59875     | 0.1221  | 5022299  | 10.24  | 0   | 510 |
| 32   | Mt    | 490520    | 43957     | 0.0896  | 4272313  | 8.71   | 0   | 512 |
| 64   | Mt    | 490520    | 10379     | 0.0212  | 533498   | 1.09   | 0   | 363 |


| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max  |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:-----|
| 0    | Pt    | 134525    | 132976    | 0.9885  | 49878823 | 370.78 | 0   | 1885 |
| 0.25 | Pt    | 134525    | 132904    | 0.9880  | 50004519 | 371.71 | 0   | 1885 |
| 0.5  | Pt    | 134525    | 133032    | 0.9889  | 49917098 | 371.06 | 0   | 1885 |
| 1    | Pt    | 134525    | 133023    | 0.9888  | 49937784 | 371.22 | 0   | 1885 |
| 2    | Pt    | 134525    | 132974    | 0.9885  | 49971095 | 371.46 | 0   | 1885 |
| 4    | Pt    | 134525    | 133102    | 0.9894  | 49893353 | 370.89 | 0   | 1885 |
| 8    | Pt    | 134525    | 132955    | 0.9883  | 49873101 | 370.73 | 0   | 1885 |
| 16   | Pt    | 134525    | 132966    | 0.9884  | 49919246 | 371.08 | 0   | 1885 |
| 32   | Pt    | 134525    | 116544    | 0.8663  | 26829492 | 199.44 | 0   | 1882 |
| 64   | Pt    | 134525    | 19865     | 0.1477  | 1997824  | 14.85  | 0   | 419  |


## Remove intermediate files

```shell script
cd ~/data/plastid/evaluation/nip

find . -type d -name "trim" | xargs rm -fr
find . -type d -name "mapping" | xargs rm -fr

find . -type f -name "*.tadpole.contig.*" | xargs rm

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm

```

