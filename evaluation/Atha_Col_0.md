# *Arabidopsis thaliana* Col-0 倍数因子测试


[TOC levels=1-3]: # ""

- [*Arabidopsis thaliana* Col-0 倍数因子测试](#arabidopsis-thaliana-col-0-倍数因子测试)
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

+ 因子值 0, 0.25, 0.5, 1, 2, 4, 8, 16, 32
+ 测试文件 SRR616966
+ 覆盖度 4,970,359,200 / 119,667,750 = 41

```shell script
BASES=$(
    faops n50 -N 0 -S \
        ~/data/organelles/ena/SRR616966_1.fastq.gz \
        ~/data/organelles/ena/SRR616966_2.fastq.gz |
    cut -f 2
)

GENOME=$(
    faops size ~/data/organelles/genome/col_0/genome.fa |
        cut -f 2 |
        perl -nl -e '
            $sum += $_;
            END {print $sum}
        '
)

echo ${BASES}
# 4970359200

echo ${GENOME}
# 119667750

bc <<< "${BASES} / ${GENOME}"

```

* 末被包含在 `anchr` 里的工具

```shell script
# 将 reads 比对到参考序列上
brew install bwa

# 计算深度
brew install mosdepth

# 转换得到 bigwig, 使用 IGV 查看比对图
pip install deeptools

```

## Symlink

```shell script
mkdir -p ~/data/organelles/evaluation/col_0
cd ~/data/organelles/evaluation/col_0

for NAME in 0 0.25 0.5 1 2 4 8 16 32; do
    BASE_NAME=SRR616966_${NAME}
    
    mkdir -p ${BASE_NAME}/1_genome
    pushd ${BASE_NAME}/1_genome
    
    ln -fs ../../../../genome/col_0/genome.fa genome.fa
    popd
    
    mkdir -p ${BASE_NAME}/2_illumina
    pushd ${BASE_NAME}/2_illumina
    
    ln -fs ../../../../ena/SRR616966_1.fastq.gz R1.fq.gz
    ln -fs ../../../../ena/SRR616966_2.fastq.gz R2.fq.gz
    popd

done

```

## Trim and cutoff

```shell script
cd ~/data/organelles/evaluation/col_0

# 倍数因子::cutoff
ARRAY=(
    '0::0'
    '0.25::10'
    '0.5::20'
    '1::40'
    '2::80'
    '4::160'
    '8::320'
    '16::640'
    '32::1280'
)

for item in "${ARRAY[@]}" ; do
    NAME="${item%%::*}"
    CUTOFF="${item##*::}"

    echo 1>&2 "==> ${item}"
    
    BASE_NAME=SRR616966_${NAME}
    pushd ${BASE_NAME}
    
    rm *.sh
        
    if [[ "${NAME}" == "0" ]]; then
        anchr template \
            --genome 119667750 \
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
            --genome 119667750 \
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
cd ~/data/organelles/evaluation/col_0

for NAME in 0 0.25 0.5 1 2 4 8 16 32; do
    BASE_NAME=SRR616966_${NAME}
    
    mkdir -p ${BASE_NAME}/kat
    pushd ${BASE_NAME}/kat
    
    bsub -q mpi -n 24 -J "${BASE_NAME}" "
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
cd ~/data/organelles/evaluation/col_0

for NAME in 0 0.25 0.5 1 2 4 8 16 32; do
    BASE_NAME=SRR616966_${NAME}
    
    mkdir -p ${BASE_NAME}/mapping
    pushd ${BASE_NAME}/mapping
    
    bsub -q mpi -n 24 -J "${BASE_NAME}-mapping" '

        # Pipe all reads together as we do not need mate info
        gzip -dcf \
            ../2_illumina/trim/Q25L60/R1.fq.gz \
            ../2_illumina/trim/Q25L60/R2.fq.gz \
            ../2_illumina/trim/Q25L60/Rs.fq.gz |
            bwa mem -M -t 20 \
                ../../../../genome/col_0/genome.fa \
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
cd ~/data/organelles/evaluation/col_0

for NAME in 0 0.25 0.5 1 2 4 8 16 32; do
    BASE_NAME=SRR616966_${NAME}
    
    echo 1>&2 "==> ${BASE_NAME}"
    
    mkdir -p ${BASE_NAME}/depth
    pushd ${BASE_NAME}/depth

    mosdepth R ../mapping/R.sort.bam
    
    popd

done

```

* Covered regions via `spanr`

```shell script
cd ~/data/organelles/evaluation/col_0

for NAME in 0 0.25 0.5 1 2 4 8 16 32; do
    BASE_NAME=SRR616966_${NAME}
    
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
        
    spanr stat ../../../../genome/col_0/chr.sizes covered.yml -o stdout |
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
cd ~/data/organelles/evaluation/col_0

for NAME in 0 0.25 0.5 1 2 4 8 16 32; do
    BASE_NAME=SRR616966_${NAME}
    
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
        (cat join.tsv | sed '2,6d' && cat) \
        > combine.tsv

    popd

done

```

## Merge all results

```shell script
cd ~/data/organelles/evaluation/col_0

for NAME in 0 0.25 0.5 1 2 4 8 16 32; do
    BASE_NAME=SRR616966_${NAME}
    
    echo 1>&2 "==> ${BASE_NAME}"
    
    mkdir -p ${BASE_NAME}/depth
    pushd ${BASE_NAME}/depth > /dev/null


    echo -e "Fold\tchrom\n${NAME}\tNc\n${NAME}\tMt\n${NAME}\tPt" |
        tsv-join -H --filter-file combine.tsv --key-fields chrom --append-fields 2-8

    popd > /dev/null

done |
    tsv-uniq \
    > SRR616966_folds.tsv

for PART in Nc Mt Pt; do
    cat SRR616966_folds.tsv |
        tsv-filter -H --str-eq chrom:${PART} |
        mlr --itsv --omd cat
    echo
    echo
done

```

| Fold | chrom | chrLength | covLength | covRate | bases      | mean  | min | max    |
|:-----|:------|:----------|:----------|:--------|:-----------|:------|:----|:-------|
| 0    | Nc    | 119146348 | 118886605 | 0.9978  | 2866737944 | 24.06 | 0   | 110403 |
| 0.25 | Nc    | 119146348 | 118741171 | 0.9966  | 2848440949 | 23.91 | 0   | 110853 |
| 0.5  | Nc    | 119146348 | 78550298  | 0.6593  | 1438090100 | 12.07 | 0   | 110836 |
| 1    | Nc    | 119146348 | 14305617  | 0.1201  | 619831007  | 5.20  | 0   | 110930 |
| 2    | Nc    | 119146348 | 9520335   | 0.0799  | 557199560  | 4.68  | 0   | 110690 |
| 4    | Nc    | 119146348 | 6788261   | 0.0570  | 515526753  | 4.33  | 0   | 110740 |
| 8    | Nc    | 119146348 | 4073269   | 0.0342  | 452934718  | 3.80  | 0   | 110423 |
| 16   | Nc    | 119146348 | 2822024   | 0.0237  | 427482435  | 3.59  | 0   | 110335 |
| 32   | Nc    | 119146348 | 1814569   | 0.0152  | 405006611  | 3.40  | 0   | 109534 |


| Fold | chrom | chrLength | covLength | covRate | bases    | mean   | min | max  |
|:-----|:------|:----------|:----------|:--------|:---------|:-------|:----|:-----|
| 0    | Mt    | 366924    | 363892    | 0.9917  | 58975468 | 160.73 | 0   | 1928 |
| 0.25 | Mt    | 366924    | 363896    | 0.9917  | 58830939 | 160.34 | 0   | 1928 |
| 0.5  | Mt    | 366924    | 363895    | 0.9917  | 58763345 | 160.15 | 0   | 1918 |
| 1    | Mt    | 366924    | 363874    | 0.9917  | 58758453 | 160.14 | 0   | 1900 |
| 2    | Mt    | 366924    | 363627    | 0.9910  | 58687612 | 159.94 | 0   | 1928 |
| 4    | Mt    | 366924    | 358240    | 0.9763  | 55603705 | 151.54 | 0   | 1902 |
| 8    | Mt    | 366924    | 98428     | 0.2683  | 7947061  | 21.66  | 0   | 1915 |
| 16   | Mt    | 366924    | 31051     | 0.0846  | 1355422  | 3.69   | 0   | 1879 |
| 32   | Mt    | 366924    | 28738     | 0.0783  | 1292139  | 3.52   | 0   | 1911 |


| Fold | chrom | chrLength | covLength | covRate | bases     | mean    | min | max  |
|:-----|:------|:----------|:----------|:--------|:----------|:--------|:----|:-----|
| 0    | Pt    | 154478    | 154478    | 1.0000  | 483364405 | 3129.02 | 4   | 4507 |
| 0.25 | Pt    | 154478    | 154478    | 1.0000  | 482256217 | 3121.84 | 4   | 4504 |
| 0.5  | Pt    | 154478    | 154478    | 1.0000  | 482127125 | 3121.01 | 4   | 4504 |
| 1    | Pt    | 154478    | 154478    | 1.0000  | 482063003 | 3120.59 | 4   | 4504 |
| 2    | Pt    | 154478    | 154478    | 1.0000  | 481992148 | 3120.13 | 4   | 4504 |
| 4    | Pt    | 154478    | 154478    | 1.0000  | 481933052 | 3119.75 | 4   | 4504 |
| 8    | Pt    | 154478    | 154478    | 1.0000  | 481921900 | 3119.68 | 4   | 4504 |
| 16   | Pt    | 154478    | 154478    | 1.0000  | 481928876 | 3119.72 | 4   | 4504 |
| 32   | Pt    | 154478    | 154478    | 1.0000  | 481751895 | 3118.58 | 2   | 4504 |


## Remove intermediate files

```shell script
cd ~/data/organelles/evaluation/col_0

find . -type d -name "trim" | xargs rm -fr
find . -type d -name "mapping" | xargs rm -fr

find . -type f -name "core.*" | xargs rm
find . -type f -name "output.*" | xargs rm

```
