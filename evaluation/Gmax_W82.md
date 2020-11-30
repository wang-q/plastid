# *Glycine max* Williams 82 倍数因子测试


[TOC levels=1-3]: # ""

- [*Glycine max* Williams 82 倍数因子测试](#glycine-max-williams-82-倍数因子测试)
  - [基本信息](#基本信息)
  - [Symlink](#symlink)
  - [Trim and cutoff](#trim-and-cutoff)

## 基本信息

`cutoff = 倍数因子 * 覆盖深度`

+ 因子值 0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64
+ 测试文件 [SRR10296600](https://www.ncbi.nlm.nih.gov/sra/SRX7009428)
+ 覆盖度 48,633,106,500 / 949,738,161 = 51

```shell script
BASES=$(
    faops n50 -N 0 -S \
        ~/data/plastid/ena/SRR10296600_1.fastq.gz \
        ~/data/plastid/ena/SRR10296600_2.fastq.gz |
    cut -f 2
)

GENOME=$(
    faops size ~/data/plastid/genome/w82/genome.fa |
        cut -f 2 |
        perl -nl -e '
            $sum += $_;
            END {print $sum}
        '
)

echo ${BASES}
# 48633106500

echo ${GENOME}
# 949738161

bc <<< "${BASES} / ${GENOME}"

```

## Symlink

```shell script
mkdir -p ~/data/plastid/evaluation/w82
cd ~/data/plastid/evaluation/w82

for NAME in 0 0.25 0.5 1 2 4 8 16 32 64; do
    BASE_NAME=SRR10296600_${NAME}
    
    mkdir -p ${BASE_NAME}/1_genome
    pushd ${BASE_NAME}/1_genome
    
    ln -fs ../../../../genome/w82/genome.fa genome.fa
    popd
    
    mkdir -p ${BASE_NAME}/2_illumina
    pushd ${BASE_NAME}/2_illumina
    
    ln -fs ../../../../ena/SRR10296600_1.fastq.gz R1.fq.gz
    ln -fs ../../../../ena/SRR10296600_2.fastq.gz R2.fq.gz
    popd

done

```
