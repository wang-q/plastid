# 🌾 *Oryza sativa* 50 accessions

[TOC levels=1-3]: # ""

- [🌾 *Oryza sativa* 50 accessions](#-oryza-sativa-50-accessions)
  - [Basic info](#basic-info)
  - [Project](#project)
  - [Other Projects](#other-projects)
  - [Download](#download)
    - [Reference](#reference)
    - [Illumina](#illumina)
  - [Symlink](#symlink)
  - [Run](#run)
  - [Pack and clean](#pack-and-clean)
  - [VCF](#vcf)


## Basic info

* Genome: GCF_001433935.1, IRGSP-1.0, 374.422 Mb
* Chloroplast: [NC_001320](https://www.ncbi.nlm.nih.gov/nuccore/NC_001320), **Japonica**, 134525 bp
* Mitochondrion: [NC_011033](https://www.ncbi.nlm.nih.gov/nuccore/NC_011033), **Japonica**, 490520
  bp

* Chloroplast: [NC_008155](https://www.ncbi.nlm.nih.gov/nuccore/NC_008155), **Indica**, 134496 bp
* Mitochondrion: [NC_007886](https://www.ncbi.nlm.nih.gov/nuccore/NC_007886), **Indica**, 491515 bp

Japonica 的叶绿体可能有些问题, 换成 Indica 的

## Project

* [SRP003189](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP003189)

* <https://www.nature.com/articles/nbt.2050>

| GROUP    | ABBREVIATION | EXPLAIN            | 中文名   |
|:---------|:-------------|:-------------------|:-------|
| Japonica | TRJ          | tropical japonica  | 热带粳稻 |
| Japonica | TEJ          | temperate japonica | 温带粳稻 |
| Japonica | ARO          | aromatic           | 香稻    |
| Indica   | AUS          | aus                | aus 稻  |
| Indica   | IND          | indica             | 籼稻    |
| III      | III          | deepwater rices    | 深水水稻 |
| IV       | IV           | deepwater rices    | 深水水稻 |

III 和 IV 是在孟加拉国和东北亚地区发现, 是一种叫做深水稻 (deepwater rices) 的类型. 目前大部分的栽培深水稻是 Indica
类型, 但是也有发现 Japonica 类型 (参考 wiki-Deepwater rice). III 只在孟加拉国和印度曼尼普尔邦发现.
它由特殊的水稻组成, 周期短, 光周期不敏感, 适应深水条件. IV 相当于孟加拉国的 Rayada rices. 这些是非常特殊的水稻, 在 11
月到 12 月播种, 12 个月后收获, 早期耐寒, 对光周期敏感, 能经受 12 天的洪水, 并能将其伸长调整到 6 米深.

## Other Projects

+ [3000水稻基因组计划 PRJEB6180](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJEB6180)
+ [《The 3,000 rice genomes project》](https://gigascience.biomedcentral.com/articles/10.1186/2047-217X-3-7)


PRJNA554986 GWAS and domesticated selection

PRJNA522896 Resequencing data of 147 rice varieties

PRJNA522923 Resequencing data of 120 rice RILs

查询 IRGC 编号: https://gringlobal.irri.org/gringlobal/search.aspx

查询种质资源信息: https://www.genesys-pgr.org/

## Download

### Reference

```shell script
mkdir -p ~/data/plastid/Osat_50/genome
cd ~/data/plastid/Osat_50/genome

for ACCESSION in "NC_008155" "NC_007886"; do
    URL=$(printf "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=%s&id=%s&retmode=text" "fasta" "${ACCESSION}");
    curl $URL -o ${ACCESSION}.fa
done

TAB=$'\t'
cat <<EOF > replace.tsv
NC_008155${TAB}Pt
NC_007886${TAB}Mt
EOF

cat NC_008155.fa NC_007886.fa |
    faops filter -s stdin stdout |
    faops replace stdin replace.tsv stdout |
    faops order stdin <(echo Pt; echo Mt) genome.fa

```

### Illumina

* Download `Metadata` from NCBI SRA Run Selector via a web browser
  * <https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP003189>
  * Save it to `SraRunTable.txt`

* Extract "Supplementary Table 1" from P37 of
  <https://static-content.springer.com/esm/art%3A10.1038%2Fnbt.2050/MediaObjects/41587_2012_BFnbt2050_MOESM7_ESM.pdf>

```shell script
mkdir -p ~/data/plastid/Osat_50/ena
cd ~/data/plastid/Osat_50/ena

cat SraRunTable.txt |
    mlr --icsv --otsv cat \
    > SraRunTable.tsv

cat <<EOF > rice.csv
Sample No.,Accession name,Sample Name,Status,Origin,Variety group
1,Mehr,IRGC12883,Landrace,Iran,AUS
2,Kalamkati,IRGC45975,Landrace,India,AUS
3,Jhona 349,IRGC6307,Landrace,India,AUS
4,DZ78,IRGC8555,Landrace,Bangladesh,AUS
5,Binulawan,IRGC26872,Landrace,Philippines,TRJ
6,Leung Pratew,IRGC27762,Landrace,Thailand,IND
7,IR 36,IRGC30416,Improve,Brazil,IND
8,Popot 165,IRGC43545,Landrace,Indonesia (E. Kalimantan),IND
9,Ai-Chiao-Hong,IRGC51250,Landrace,China,IND
10,Guan-Yin-Tsan,IRGC51300,Landrace,China,IND
11,Gie 57,IRGC8231,Landrace,Vietnam,IND
12,TD2,IRGC9148,Elite,Thailand,IND
13,JC91,IRGC9177,Elite,India,IND
14,Ta Hung Ku,IRGC1107,Landrace,China,TEJ
15,Haginomae Mochi,IRGC2540,Elite,Japan,TEJ
16,Darmali,IRGC27630,Landrace,Nepal,TEJ
17,Phudugey,IRGC32399,Landrace,Bhutan,TEJ
18,Norin 20,IRGC418,Landrace,Japan,TEJ
19,Chodongji,IRGC55471,Landrace,South Korea,TEJ
20,Mansaku,IRGC8191,Landrace,Japan,TEJ
21,Nipponbare,NP,Elite,Japan,TEJ
22,Maintmolotsy,IRGC11010,Elite,Madagascar,TRJ
23,Jambu,IRGC17757,Landrace,Indonesia,TRJ
24,Miriti,IRGC25901,Landrace,Bangladesh,IND
25,AZUCENA,IRGC328,Landrace,Philippines,TRJ
26,NPE 844,IRGC38698,Landrace,Pakistan,TRJ
27,Arias,IRGC43325,Landrace,Indonesia (West Java),TRJ
28,Gotak Gatik,IRGC43397,Landrace,Indonesia (C. Java),TRJ
29,Trembese,IRGC43675,Landrace,Indonesia (East Java),TRJ
30,Canella De Ferro,IRGC50448,Elite,Brazil,TRJ
31,Lemont,IRGC66756,Elite,"TX,USA",TRJ
32,Davao,IRGC8244,Landrace,Philippines,TRJ
33,Kitrana 508,IRGC12793,Elite,Madagascar,ARO
34,Bico Branco,IRGC38994,Elite,Brazil,ARO
35,JC101,IRGC9060,Elite,India,ARO
36,JC111,IRGC9062,Elite,India,ARO
37,Firooz,RA4952,Landrace,Iran,ARO
38,KUI SALI,IRGC31856,Landrace,Bangladesh,ARO
39,HAISHA CAMAN,IRGC60542,Landrace,Bangladesh,IV
40,BADAL 89,IRGC6513,Landrace,Bangladesh,III
41,042/87/34,IRGC105327,wild rice,"Dhoni, India",nivara
42,MV 89-80,IRGC106105,wild rice,"Medinipur, India",nivara
43,L 89-12,IRGC106154,wild rice,"Vientiane, Laos",nivara
44,HK 47,IRGC80470,wild rice,"Madhya Pradesh, India",nivara
45,CA 97-053,IRGC89215,wild rice,"Sopoir Tep, Cambodia",nivara
46,PADI PADIAN,IRGC105958,wild rice,"Kromat Watu, Indonesia",rufipogon
47,DAL DHAN,IRGC105960,wild rice,"Chakaria, Bangladesh",rufipogon
48,VOC4,VOC4,wild rice,Nepal,rufipogon
49,P46,P46,wild rice,"Hainan, China",rufipogon
50,Yuan 3-9,Yuan3-9,wild rice,"Yunnan, China",rufipogon
EOF

cat rice.csv |
    mlr --icsv --omd cat

cat rice.csv |
    mlr --icsv --otsv cat |
    tsv-join -H --key-fields "Sample\ Name" \
        -f <(
            cat SraRunTable.tsv | tsv-filter -H --ge AvgSpotLen:100
            ) \
        --append-fields AvgSpotLen,Bases,Experiment |
    tsv-select -H -f Experiment,"Sample\ Name","Variety\ group" |
    mlr --itsv --ocsv cat |
    sed 1d \
    > source.csv

anchr ena info | perl - -v source.csv > ena_info.yml
anchr ena prep | perl - ena_info.yml

mlr --icsv --omd cat ena_info.csv

aria2c -j 4 -x 4 -s 2 --file-allocation=none -c -i ena_info.ftp.txt

md5sum --check ena_info.md5.txt

```

| Sample No. | Accession name   | Sample Name | Status    | Origin                    | Variety group |
|:-----------|:-----------------|:------------|:----------|:--------------------------|:--------------|
| 1          | Mehr             | IRGC12883   | Landrace  | Iran                      | AUS           |
| 2          | Kalamkati        | IRGC45975   | Landrace  | India                     | AUS           |
| 3          | Jhona 349        | IRGC6307    | Landrace  | India                     | AUS           |
| 4          | DZ78             | IRGC8555    | Landrace  | Bangladesh                | AUS           |
| 5          | Binulawan        | IRGC26872   | Landrace  | Philippines               | TRJ           |
| 6          | Leung Pratew     | IRGC27762   | Landrace  | Thailand                  | IND           |
| 7          | IR 36            | IRGC30416   | Improve   | Brazil                    | IND           |
| 8          | Popot 165        | IRGC43545   | Landrace  | Indonesia (E. Kalimantan) | IND           |
| 9          | Ai-Chiao-Hong    | IRGC51250   | Landrace  | China                     | IND           |
| 10         | Guan-Yin-Tsan    | IRGC51300   | Landrace  | China                     | IND           |
| 11         | Gie 57           | IRGC8231    | Landrace  | Vietnam                   | IND           |
| 12         | TD2              | IRGC9148    | Elite     | Thailand                  | IND           |
| 13         | JC91             | IRGC9177    | Elite     | India                     | IND           |
| 14         | Ta Hung Ku       | IRGC1107    | Landrace  | China                     | TEJ           |
| 15         | Haginomae Mochi  | IRGC2540    | Elite     | Japan                     | TEJ           |
| 16         | Darmali          | IRGC27630   | Landrace  | Nepal                     | TEJ           |
| 17         | Phudugey         | IRGC32399   | Landrace  | Bhutan                    | TEJ           |
| 18         | Norin 20         | IRGC418     | Landrace  | Japan                     | TEJ           |
| 19         | Chodongji        | IRGC55471   | Landrace  | South Korea               | TEJ           |
| 20         | Mansaku          | IRGC8191    | Landrace  | Japan                     | TEJ           |
| 21         | Nipponbare       | NP          | Elite     | Japan                     | TEJ           |
| 22         | Maintmolotsy     | IRGC11010   | Elite     | Madagascar                | TRJ           |
| 23         | Jambu            | IRGC17757   | Landrace  | Indonesia                 | TRJ           |
| 24         | Miriti           | IRGC25901   | Landrace  | Bangladesh                | IND           |
| 25         | AZUCENA          | IRGC328     | Landrace  | Philippines               | TRJ           |
| 26         | NPE 844          | IRGC38698   | Landrace  | Pakistan                  | TRJ           |
| 27         | Arias            | IRGC43325   | Landrace  | Indonesia (West Java)     | TRJ           |
| 28         | Gotak Gatik      | IRGC43397   | Landrace  | Indonesia (C. Java)       | TRJ           |
| 29         | Trembese         | IRGC43675   | Landrace  | Indonesia (East Java)     | TRJ           |
| 30         | Canella De Ferro | IRGC50448   | Elite     | Brazil                    | TRJ           |
| 31         | Lemont           | IRGC66756   | Elite     | TX,USA                    | TRJ           |
| 32         | Davao            | IRGC8244    | Landrace  | Philippines               | TRJ           |
| 33         | Kitrana 508      | IRGC12793   | Elite     | Madagascar                | ARO           |
| 34         | Bico Branco      | IRGC38994   | Elite     | Brazil                    | ARO           |
| 35         | JC101            | IRGC9060    | Elite     | India                     | ARO           |
| 36         | JC111            | IRGC9062    | Elite     | India                     | ARO           |
| 37         | Firooz           | RA4952      | Landrace  | Iran                      | ARO           |
| 38         | KUI SALI         | IRGC31856   | Landrace  | Bangladesh                | ARO           |
| 39         | HAISHA CAMAN     | IRGC60542   | Landrace  | Bangladesh                | IV            |
| 40         | BADAL 89         | IRGC6513    | Landrace  | Bangladesh                | III           |
| 41         | 042/87/34        | IRGC105327  | wild rice | Dhoni, India              | nivara        |
| 42         | MV 89-80         | IRGC106105  | wild rice | Medinipur, India          | nivara        |
| 43         | L 89-12          | IRGC106154  | wild rice | Vientiane, Laos           | nivara        |
| 44         | HK 47            | IRGC80470   | wild rice | Madhya Pradesh, India     | nivara        |
| 45         | CA 97-053        | IRGC89215   | wild rice | Sopoir Tep, Cambodia      | nivara        |
| 46         | PADI PADIAN      | IRGC105958  | wild rice | Kromat Watu, Indonesia    | rufipogon     |
| 47         | DAL DHAN         | IRGC105960  | wild rice | Chakaria, Bangladesh      | rufipogon     |
| 48         | VOC4             | VOC4        | wild rice | Nepal                     | rufipogon     |
| 49         | P46              | P46         | wild rice | Hainan, China             | rufipogon     |
| 50         | Yuan 3-9         | Yuan3-9     | wild rice | Yunnan, China             | rufipogon     |


| name       | srx       | platform | layout | ilength | srr       | spots    | bases |
|:-----------|:----------|:---------|:-------|:--------|:----------|:---------|:------|
| IRGC105327 | SRX025244 | ILLUMINA | PAIRED | 210     | SRR063622 | 33451859 | 6.23G |
| IRGC105958 | SRX025231 | ILLUMINA | PAIRED | 207     | SRR063609 | 32866744 | 6.12G |
| IRGC105960 | SRX025242 | ILLUMINA | PAIRED | 180     | SRR063620 | 33384518 | 6.22G |
| IRGC106105 | SRX025243 | ILLUMINA | PAIRED | 188     | SRR063621 | 31910374 | 5.94G |
| IRGC106154 | SRX025232 | ILLUMINA | PAIRED | 216     | SRR063610 | 30653350 | 5.71G |
| IRGC11010  | SRX025212 | ILLUMINA | PAIRED | 211     | SRR063590 | 32691730 | 4.57G |
| IRGC1107   | SRX025259 | ILLUMINA | PAIRED | 479     | SRR063637 | 31125089 | 5.8G  |
| IRGC12793  | SRX025228 | ILLUMINA | PAIRED | 209     | SRR063606 | 36129012 | 6.73G |
| IRGC12883  | SRX025217 | ILLUMINA | PAIRED | 184     | SRR063595 | 34542675 | 6.43G |
| IRGC17757  | SRX025255 | ILLUMINA | PAIRED | 494     | SRR063633 | 31586148 | 5.88G |
| IRGC2540   | SRX025225 | ILLUMINA | PAIRED | 192     | SRR063603 | 35360096 | 6.59G |
| IRGC25901  | SRX025256 | ILLUMINA | PAIRED | 481     | SRR063634 | 30739199 | 5.73G |
| IRGC26872  | SRX025247 | ILLUMINA | PAIRED | 480     | SRR063625 | 35401204 | 6.59G |
| IRGC27630  | SRX025213 | ILLUMINA | PAIRED | 186     | SRR063591 | 31322347 | 4.38G |
| IRGC27762  | SRX025237 | ILLUMINA | PAIRED | 202     | SRR063615 | 31875248 | 5.94G |
| IRGC30416  | SRX025238 | ILLUMINA | PAIRED | 209     | SRR063616 | 33949026 | 6.32G |
| IRGC31856  | SRX025240 | ILLUMINA | PAIRED | 212     | SRR063618 | 31806012 | 5.92G |
| IRGC32399  | SRX025218 | ILLUMINA | PAIRED | 184     | SRR063596 | 31073220 | 5.79G |
| IRGC328    | SRX025253 | ILLUMINA | PAIRED | 479     | SRR063631 | 33861047 | 6.31G |
| IRGC38698  | SRX025220 | ILLUMINA | PAIRED | 211     | SRR063598 | 35350187 | 6.58G |
| IRGC38994  | SRX025227 | ILLUMINA | PAIRED | 182     | SRR063605 | 33452203 | 6.23G |
| IRGC418    | SRX025251 | ILLUMINA | PAIRED | 486     | SRR063629 | 33833808 | 6.3G  |
| IRGC43325  | SRX025261 | ILLUMINA | PAIRED | 465     | SRR063639 | 32123287 | 5.98G |
| IRGC43397  | SRX025257 | ILLUMINA | PAIRED | 496     | SRR063635 | 31483704 | 5.86G |
| IRGC43545  | SRX025248 | ILLUMINA | PAIRED | 480     | SRR063626 | 33663316 | 6.27G |
| IRGC43675  | SRX025258 | ILLUMINA | PAIRED | 475     | SRR063636 | 31226985 | 5.82G |
| IRGC45975  | SRX025226 | ILLUMINA | PAIRED | 182     | SRR063604 | 32928990 | 6.13G |
| IRGC50448  | SRX025229 | ILLUMINA | PAIRED | 205     | SRR063607 | 35397687 | 6.59G |
| IRGC51250  | SRX025249 | ILLUMINA | PAIRED | 480     | SRR063627 | 34581348 | 6.44G |
| IRGC51300  | SRX025250 | ILLUMINA | PAIRED | 468     | SRR063628 | 34005135 | 6.33G |
| IRGC55471  | SRX025222 | ILLUMINA | PAIRED | 212     | SRR063600 | 32770369 | 6.1G  |
| IRGC60542  | SRX025236 | ILLUMINA | PAIRED | 200     | SRR063614 | 30512958 | 5.68G |
| IRGC6307   | SRX025214 | ILLUMINA | PAIRED | 199     | SRR063592 | 25577558 | 4.76G |
| IRGC6513   | SRX025219 | ILLUMINA | PAIRED | 195     | SRR063597 | 32458498 | 6.05G |
| IRGC66756  | SRX025235 | ILLUMINA | PAIRED | 176     | SRR063613 | 29913714 | 5.57G |
| IRGC80470  | SRX025230 | ILLUMINA | PAIRED | 235     | SRR063608 | 29176451 | 5.43G |
| IRGC8191   | SRX025252 | ILLUMINA | PAIRED | 471     | SRR063630 | 33751895 | 6.29G |
| IRGC8231   | SRX025223 | ILLUMINA | PAIRED | 209     | SRR063601 | 32425786 | 6.04G |
| IRGC8244   | SRX025254 | ILLUMINA | PAIRED | 477     | SRR063632 | 32209551 | 6G    |
| IRGC8555   | SRX025221 | ILLUMINA | PAIRED | 224     | SRR063599 | 33404831 | 6.22G |
| IRGC89215  | SRX025233 | ILLUMINA | PAIRED | 212     | SRR063611 | 28881277 | 5.38G |
| IRGC9060   | SRX025234 | ILLUMINA | PAIRED | 205     | SRR063612 | 29478266 | 5.49G |
| IRGC9062   | SRX025216 | ILLUMINA | PAIRED | 191     | SRR063594 | 32191480 | 6G    |
| IRGC9148   | SRX025215 | ILLUMINA | PAIRED | 279     | SRR063593 | 22284854 | 4.15G |
| IRGC9177   | SRX025239 | ILLUMINA | PAIRED | 187     | SRR063617 | 34117844 | 6.35G |
| NP         | SRX025260 | ILLUMINA | PAIRED | 463     | SRR063638 | 29784011 | 5.55G |
| P46        | SRX025245 | ILLUMINA | PAIRED | 185     | SRR063623 | 27365410 | 5.1G  |
| RA4952     | SRX025224 | ILLUMINA | PAIRED | 201     | SRR063602 | 34059182 | 6.34G |
| VOC4       | SRX025241 | ILLUMINA | PAIRED | 195     | SRR063619 | 34890729 | 6.5G  |
| Yuan3-9    | SRX025246 | ILLUMINA | PAIRED | 230     | SRR063624 | 23867062 | 4.45G |


## Symlink

* 采用的倍数因子值: `2`

* Rsync to hpcc

```shell script
rsync -avP \
    ~/data/plastid/Osat_50/ \
    wangq@202.119.37.251:data/plastid/Osat_50

# rsync -avP wangq@202.119.37.251:data/plastid/Osat_50/ ~/data/plastid/Osat_50

```

```shell script
cd ~/data/plastid/Osat_50/

export FOLD=2
export GENOME_SIZE=$(
    cat ~/data/plastid/genome/nip/chr.sizes |
        tsv-summarize --sum 2
)

cat ena/ena_info.csv |
    mlr --icsv --otsv cat |
    tsv-select -H -f name,srr,bases |
    perl -nla -F'\t' -e '
        BEGIN { our %seen }
        /^name/ and next;
        $seen{$F[0]} and next;
        my $bases = $F[2];
        $bases =~ s/G$//;
        my $cutoff = $bases * 1000 * 1000 * 1000 / $ENV{GENOME_SIZE} * $ENV{FOLD};
        $cutoff = int $cutoff;
        print join qq(\t), ($F[0], $F[1], $cutoff);
        $seen{$F[0]}++;
    ' \
    > opts.tsv

cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ ! -f ena/{2}_1.fastq.gz ]; then
            exit;
        fi
        if [ ! -f ena/{2}_2.fastq.gz ]; then
            exit;
        fi

        if [ -f {1}.tar.gz ]; then
            exit;
        fi

        mkdir -p {1}/1_genome
        pushd {1}/1_genome

        cp ../../genome/genome.fa genome.fa
        popd > /dev/null

        mkdir -p {1}/2_illumina
        pushd {1}/2_illumina

        ln -fs ../../ena/{2}_1.fastq.gz R1.fq.gz
        ln -fs ../../ena/{2}_2.fastq.gz R2.fq.gz
        popd > /dev/null
    '

```

## Run

```shell script
cd ~/data/plastid/Osat_50/

cat opts.tsv | #head -n 20 | #tail -n 10 |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ -f {1}.tar.gz ]; then
            exit;
        fi

        if [ ! -d {1} ]; then
            exit;
        fi

        if [ ! -e {1}/2_illumina/R1.fq.gz ]; then
            exit;
        fi
        if [ ! -e {1}/2_illumina/R2.fq.gz ]; then
            exit;
        fi

        if bjobs -w | tr -s " " | cut -d " " -f 7 | grep -w "^{1}$"; then
            echo Job {1} exists
            exit;
        fi

        cd {1}

        echo {1}

        rm *.sh
        anchr template \
            --genome 1000000 \
            --parallel 24 \
            --xmx 80g \
            \
            --fastqc \
            --insertsize \
            --kat \
            \
            --trim "--dedupe --cutoff {3} --cutk 31" \
            --qual "25" \
            --len "60" \
            --filter "adapter artifact" \
            \
            --bwa Q25L60 \
            --gatk

        bsub -q mpi -n 24 -J "{1}" "
            bash 2_fastqc.sh
            bash 2_insert_size.sh
            bash 2_kat.sh
            bash 2_trim.sh
            bash 9_stat_reads.sh
            bash 3_bwa.sh
            bash 3_gatk.sh
        "
    '

```

## Pack and clean

```shell script
cd ~/data/plastid/Osat_50/

cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
        if [ -f {1}.tar.gz ]; then
            echo "==> {1} .tar.gz"
            exit;
        fi

        if [ ! -f {1}/3_gatk/R.filtered.vcf ]; then
            echo "==> {1} 3_gatk"
            exit;
        fi

#        if [ ! -f {1}/7_merge_anchors/anchor.merge.fasta ]; then
#            echo "==> {1} 7_merge_anchors"
#            exit;
#        fi
#
#        if [ ! -d {1}/9_quast ]; then
#            echo "==> {1} 9_quast"
#            exit;
#        fi

        echo "==> Clean {1}"
        bash {1}/0_cleanup.sh

        echo "==> Create {1}.tar.gz"

        tar -czvf {1}.tar.gz \
            {1}/1_genome/genome.fa \
            {1}/2_illumina/fastqc \
            {1}/2_illumina/insert_size \
            {1}/2_illumina/kat \
            {1}/3_bwa/join.tsv \
            {1}/3_bwa/R.dedup.metrics \
            {1}/3_bwa/R.wgs.metrics \
            {1}/3_bwa/R.sort.bam \
            {1}/3_bwa/R.sort.bai \
            {1}/3_gatk \
            {1}/7_merge_anchors/anchor.merge.fasta \
            {1}/7_merge_anchors/others.non-contained.fasta \
            {1}/8_megahit/anchor/anchor.fasta \
            {1}/8_megahit/megahit.non-contained.fasta \
            {1}/8_mr_megahit/anchor/anchor.fasta \
            {1}/8_mr_megahit/megahit.non-contained.fasta \
            {1}/8_spades/anchor/anchor.fasta \
            {1}/8_spades/spades.non-contained.fasta \
            {1}/8_mr_spades/anchor/anchor.fasta \
            {1}/8_mr_spades/spades.non-contained.fasta \
            {1}/8_platanus/anchor/anchor.fasta \
            {1}/8_platanus/platanus.non-contained.fasta \
            {1}/9_quast \
            {1}/*.md

        echo
    '

cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ -f {1}.tar.gz ]; then
            if [ -d {1} ]; then
                echo "==> Remove {1}/"
                rm -fr {1}
            fi
        fi
    '

```

* Remove processed files

```shell script
cd ~/data/plastid/Osat_50/

cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ ! -f {1}.tar.gz ]; then
            exit;
        fi

        find ena -type f -name "*{2}*"
    ' |
    xargs rm

```

* Unpack

```shell script
cd ~/data/plastid/Osat_50/

cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        if [ -d {1} ]; then
            echo "==> {1} exists"
            exit;
        fi

        if [ ! -f {1}.tar.gz ]; then
            echo "==> {1}.tar.gz not exists"
            exit;
        fi

        tar -xzvf {1}.tar.gz
        rm {1}.tar.gz
    '

```

## VCF

```shell script
cd ~/data/plastid/Osat_50/

mkdir -p vcf

cat opts.tsv |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 4 '
        if [ ! -f {1}.tar.gz ]; then
            echo "==> {1}.tar.gz not exists"
            exit;
        fi

        tar -xOzvf {1}.tar.gz {1}/3_gatk/R.filtered.vcf |
            bcftools reheader --samples <(echo {1}) |
            bcftools view \
                --apply-filters PASS --types snps --max-alleles 2 --targets Pt -Oz |
            bcftools view --include "AF>0.01" -Oz -o vcf/{1}.vcf.gz

        bcftools index -f vcf/{1}.vcf.gz
    '

bcftools merge --merge all -l <(
        cat opts.tsv |
            cut -f 1 |
            parallel -k -j 1 ' [ -f vcf/{}.vcf.gz ] && echo "vcf/{}.vcf.gz" '
    ) \
    > Osat_50.vcf

rm -fr vcf

```

