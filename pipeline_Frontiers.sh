[TOC]

    # 名称 Name: 易扩增子(EasyAmplicon)
    # 作者 Authors: 刘永鑫(Yong-Xin Liu), 陈同(Tong Chen)，周欣(Xin Zhou)等
    # 版本 Version: v1.10
    # 更新 Update: 2020-12-04

    # 设置流程EasyAmplicon(ea)和工作目录(work directory, wd)，添加环境变量并进入wd
    # **每次打开Rstudio必须运行下面4行**
    ea=/d/gmu/EasyAmplicon
    wd=/d/gmu/lt1/
    PATH=$PATH:${ea}/win/
    cd ${wd}
# 22、扩增子分析流程 Analysis pipeline for 16S amplicon 

    # 系统要求 System requirement: Windows 10 / Mac OS 10.12+ / Ubuntu 18.04+
    # 依赖软件和数据 Sofware and database dependencies: 复制public目录到C盘
    # gitforwidnows 2.28.0 http://gitforwindows.org/(Windows only)
    # R 4.0.2 https://www.r-project.org/
    # Rstudio 1.3.1056 https://www.rstudio.com/products/rstudio/download/#download
    # vsearch v2.15.0 https://github.com/torognes/vsearch/releases
    # usearch v10.0.240 https://www.drive5.com/usearch/download.html

    # 运行前准备
    # 1. 将EasyAmplicon目录下发 /复制到指定位置，如Widnows的C盘(C:/) 或 Mac/Linux服务器家目录(~/)
    # 2. 按EasyAmplicon中Readme.md中说明安装依赖软件和数据库
    # 3. 准备流程脚本pipeline.sh、样本元数据metadata.txt和测序数据seq/*.fq.gz
 
## 1. 准备输入数据(测试数据是流程正对照)

    #1. 分析流程pipeline.sh，每个项目复制一份，可修改实现分析过程全记录(可重复分析)
    #2. 样本元数据/实验设计metadata.txt
    # wget -c http://210.75.224.110/github/EasyAmplicon/data/metadata.txt
    head -n2 metadata.txt
    #3. 原始测序数据保存于seq目录，通常以`.fq.gz`结尾，每个样品一对文件
    # 下载测试数据，按元数据中GSA的CRA(批次)和CRR(样品)编号从公共数据库下载
    # 示例下载单个文件并改名
    # wget -c ftp://download.big.ac.cn/gsa/CRA002352/CRR117575/CRR117575_f1.fq.gz -O seq/KO1_1.fq.gz
    # 按实验设计编号批量下载并改名，速度单位10M/s，家500K/s
    # mkdir -p seq
    # awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$5"/"$6"/"$6"_f1.fq.gz -O seq/"$1"_1.fq.gz")}' \
    #     <(tail -n+2 metadata.txt)
    # awk '{system("wget -c ftp://download.big.ac.cn/gsa/"$5"/"$6"/"$6"_r2.fq.gz -O seq/"$1"_2.fq.gz")}' \
    #     <(tail -n+2 metadata.txt)
    ls -sh seq
    #4. 创建临时和结果目录，临时目录分析结束可删除
    mkdir -p temp result

### 1.1. metadata.txt 实验设计文件

    #cat查看前3行，-A显示符号
    cat -A metadata.txt | head -n3
    #windows用户如果结尾有^M，运行sed命令去除，并cat -A检查结果
    sed -i 's/\r/\n/' metadata.txt
    cat -A metadata.txt | head -n3

### 1.2. seq/*.fq.gz 原始测序数据

    # 公司返回的测序结果，通常为一个样品一对fq/fastq.gz格式压缩文件
    # 文件名与样品名务必对应：不一致时手工修改，批量改名见"常见问题6"
    # 测序数据是.gz的压缩文件，可以使用zcat/zless查看，连用head查看几行
    # zless/less按页查看，空格翻页、q退出；head查看前10行，-n指定行
    # zcat seq/KO1_1.fq.gz|head -n4
    # # gz压缩文件可以使用gunzip解压作为USEARCH输入，vsearch可直接读取gz文件作为输入
    # gunzip seq/*.fq.gz
    # # gzip seq/*.fastq
    # ls -sh seq/
    # # # 每行太长，指定查看每行的1-60个字符
    # cut -c 1-60 seq/LI1.fastq | head -n4

### 1.3. pipeline.sh 流程依赖数据库

    # 数据库第一次使用需要解压，解压过可跳过此段

    # usearchs可用16S/18S/ITS数据库：RDP, SILVA和UNITE，本地文件位置 ${ea}/usearch/
    # usearch数据库database下载页: http://www.drive5.com/sintax
    # 解压rdp用于物种注释和silva用于去嵌合体，*代表任意字符，选择以fa.gz结尾的文件
    # time gunzip ${ea}/usearch/*.fa.gz  # 51s
    # # # 解压缩数据库
    # gunzip ${ea}/usearch/rdp_16s_v16.fa.gz
    # gzip ${ea}/usearch/rdp_16s_v18.fa.gz
    # # # 查看注释数据库的格式
    # head -n2 ${ea}/usearch/rdp_16s_v16_sp.fa
    # # 
    # # # QIIME greengene 13_8有参数据库用于功能注释: ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_8_otus.tar.gz
    # # gunzip ${ea}/gg/97_otus.fasta.gz
    # # 不用时压缩结省空间
    # gzip ${ea}/gg/97_otus.fasta


## 2. 合并双端序列并按样品重命名 Merge paired reads and label samples

    # 以WT1单样品合并为例
    # vsearch --fastq_mergepairs seq/WN2_2.fastq.gz \
    #   --reverse seq/WN2_1.fastq.gz \
    #   --fastqout temp/WN2.fq \
    #   --relabel WN2.

    #依照实验设计批处理并合并。
    #tail -n+2去表头，cut -f 1取第一列，即获得样本列表。5万对序列的18个样本合并耗时约2分钟。
    #Windows下复制命令Ctrl+C在Linux下终止，为防止异常中断，结尾添加&使命令转后台。 

    for i in `tail -n+2 metadata.txt | cut -f 1`;do
       vsearch --fastq_mergepairs seq/${i}_1.fastq.gz --reverse seq/${i}_2.fastq.gz --fastqout temp/${i}.merged.fq --relabel ${i}.
    done &
    
    # usearch -fastx_relabel seq/LE1.fastq -fastqout temp/LE1.merged.fq -prefix LE1.
    # # # 己经合并样本直接改名，接入分析流程；替换上面vsearch --fastq_mergepairs命令为
    # for i in `tail -n+2 metadata.txt | cut -f 1`;do
    #    usearch -fastx_relabel seq/${i}.fastq -fastqout temp/${i}.merged.fq -prefix ${i}.
    # done
    # 另一种方法参考“常见问题2”

    # #查看文件序列名
    head -n1 seq/L_1_1.fastq.gz
    #序列按样本命名，并输出到新文件夹
    for i in `tail -n+2 metadata2.txt | cut -f 1`;do
    vsearch --fastq_convert seq/${i}.fastq.gz --fastqout temp/${i}.merged.fq --relabel ${i}.
    done
    #查看转换后33编码格式，质量值多为大写字母
    head -n1 FAQ/relabel/L_1_1.fq
    

# ## 3. 引物切除和质量控制 Cut primers and quality filter

    # 采用等长的方式切除引物，引物外侧如果有标签，标签的长度需要计算在内。
    # 如本例的结果为左端标签10 bp + 正向引物V5 19 bp共为29 bp，反向引物V7 18 bp。
    # Cut barcode 10bp + V5 19bp in left and V7 18bp in right
    # 注意：务必清楚实验设计中引物和标签长度，如果引物已经去除，可在下方参数处填0表示无需去除。76万条序列37s  --fastq_maxlen 400 --fastq_minlen 150 
    
    for i in `tail -n+2 metadata.txt | cut -f 1`;do
    vsearch --fastq_filter temp/${i}.merged.fq --fastq_stripleft 23 --fastq_stripright 27 -fastq_maxee_rate 0.01 --fastaout temp/${i}.filtered.fa --fastq_qmax 42
    done
    # # 
    # 
   # #合并所有样品至同一文件
   cat temp/*.filtered.fa > temp/filtered.fa
   cat temp/NC/*.filtered.fa > temp/Wall.filtered.fa
   #查看文件大小634M，软件不同版本结果略有差异
   ls -lsh temp/filtered.fa
   #查看序列名.之前是否为样本名，样本名绝不允许有.
   head -n 6 temp/filtered.fa|cut -c1-60
   
    # 查看文件了解fa文件格式
    head temp/filtered.fa
    less temp/filtered.fa|grep '>' -c ## 查看序列数   # usearch -fastx_relabel seq/${i}.fq -fastqout temp/${i}.merged.fq -prefix ${i}.
    # less temp/Wall.filtered.fa|grep '>' -c
    # 另一种方法参考“常见问题2”

## 4. 序列去冗余并挑选代表序列(OTU/ASV) Dereplicate and cluster/denoise

### 4.1 序列去冗余 Dereplication
 vsearch --derep_fulllength temp/filtered.fa \
      --output temp/uniques.fa \
      --relabel Uni --minuniquesize 10 --sizeout
    #高丰度非冗余序列非常小(<2Mb),名称后有size和频率
    ls -lsh temp/uniques.fa
    head -n 2 temp/uniques.fa
    less temp/uniques.fa|grep '>' -c
### 4.2 聚类OTU/去噪ASV Cluster OTUs / denoise ASV

    # 去冗余数据量起码降低1个数量级，减小下游分析工作量，也更适合基于丰度鉴定真实特征序列，如可操作分类单元(OTUs)或扩增序列变化(ASVs)。
    # 参数--miniuniqusize设置使用序列的最小出现次数，默认为8，此处设置为10，推荐最小为总数据量的百万分之一，可实现去除低丰度噪音并增加计算速度。
    # -sizeout输出丰度, --relabel添加序列前缀。
   
    # OTU和ASV统称为特征(Feature)，它们的区别是：
    # OTU通常按97%聚类后挑选最高丰度或中心的代表性序列；
    # ASV是基于序列进行去噪(排除或校正错误序列，并挑选丰度较高的可信序列)作为代表性序列

    #有两种方法：推荐unoise3去噪获得单碱基精度ASV，传统的97%聚类OTU (属水平精度)供备选
    #usearch两种特征挑选方法均自带de novo去嵌合体

    #方法1. 97%聚类OTU，适合大数据/ASV规律不明显/reviewer要求
    #结果耗时6s, 产生878 OTUs, 去除320 chimeras
    # usearch -cluster_otus temp/uniques.fa \
    #  -otus temp/otus.fa \
    #  -relabel OTU_

    #方法2. ASV去噪 Denoise: predict biological sequences and filter chimeras
    #59s, 2920 good, 227 chimeras
    usearch -unoise3 temp/uniques.fa \
      -zotus temp/zotus.fa
    #修改序列名：Zotu为改为ASV方便识别
    sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/otus.fa
    head -n 2 temp/otus.fa
    # ## 去singleton序列
    # vsearch --sortbysize otus1.fa --output otus.fa --minsize 2
    
    #方法3. 数据过大无法使用usearch时，备选vsearch方法见"常见问题3"

### 4.3 基于参考去嵌合 Reference-based chimera detect

    # 不推荐，容易引起假阴性，因为参考数据库无丰度信息，
    # 而de novo去嵌合时要求亲本丰度为嵌合体16倍以上防止假阴性
    # 因为已知序列不会被去除，数据库选择越大越合理，假阴性率最低
    mkdir -p result/raw

    # 方法1. vsearch+rdp去嵌合(快但容易假阴性)，或
    # silva去嵌合(silva_16s_v123.fa)，推荐(慢，耗时15m ~ 3h，但更好)
    vsearch --uchime_ref temp/otus.fa \
    -db ${ea}/usearch/rdp_16s_v16_sp.fa \
    --nonchimeras result/raw/otus.fa
    # RDP: 51s, 250 (8.6%) chimeras; SILVA：10m, 255 (8.7%) chimeras
    # Win vsearch结果添加了windows换行符^M需删除，mac不要执行此命令
    sed -i 's/\r//g' result/raw/otus.fa

    # 方法2. 不去嵌合
    cp -f temp/otus.fa result/raw/otus.fa
    less result/raw/otus.fa|grep '>' -c
    
# ## 5. 特征表生成和筛选 Feature table create & filter
# # 删除污染率高于0.99的序列，按需设置
vsearch --usearch_global result/raw/otus.fa \
  --db temp/Wall.filtered.fa \
  --id 0.99 --query_cov 0.99 \
  --strand both --biomout result/raw/otus.txt --alnout result/raw/otus.aln \
  --blast6out result/raw/otus_blast.xls --fastapairs result/raw/otus_pairs.fa \
  --notmatched result/raw/otus_notmatched.fa --userfields query+target+id+qcov+tcov \
  --userout result/raw/otus_stat.xls
less result/raw/otus.fa|grep '>' -c
less result/raw/otus_notmatched.fa|grep '>' -c
less temp/filtered.fa|grep '>' -c
# # # ### 5.1 生成特征表 Creat Feature table

    # 方法1. usearch生成特征表，小样本(<30)快；但大样本受限且多线程效率低，84.1%, 4核1m
    usearch -otutab temp/filtered.fa -otus result/raw/otus_notmatched.fa \
        -otutabout result/raw/otutab.txt -threads 12

    # # 方法2. vsearch生成特征表
    vsearch --usearch_global temp/filtered.fa --db result/raw/otus_notmatched.fa \
      --otutabout result/raw/otutab.txt --id 0.97 --threads 12
      
    # # # # # #660953 of 761432 (86.80%)可比对，耗时8m
    # windows用户删除换行符^M
    sed -i 's/\r//' result/raw/otutab.txt
    head -n3 result/raw/otutab.txt |cat -A
    wc -l result/raw/otutab.txt
    less result/raw/otus.fa|grep '>' -c
    less result/raw/otus_notmatched.fa|grep '>' -c
    
### 5.2 物种注释-去除质体和非细菌/古菌并统计比例(可选) Remove plastid and non-Bacteria

    # RDP物种注释(rdp_16s_v16_sp)数据库小，分类速度极快，但缺少完整真核来源数据，耗时15s;
    # SILVA数据库(silva_16s_v123.fa)更好注释真核、质体序列，3h
    # --sintax_cutoff设置分类的可信度阈值，常用0.6/0.8，越大注释比例越低。
   
    vsearch --sintax result/raw/otus_notmatched.fa --db ${ea}/usearch/silva_16s_v123.fa.gz \
    --tabbedout result/raw/otus.sintax --sintax_cutoff 0.6
    vsearch --sintax result/raw/otus_notmatched.fa --db ${ea}/usearch/ltp_16s_v123.fa.gz \
      --tabbedout result/raw/ltp/otus.sintax --sintax_cutoff 0.6
    vsearch --sintax result/raw/otus_notmatched.fa --db ${ea}/usearch/rdp_16s_v18_sp.fa.gz \
       --tabbedout result/raw/rdp/otus.sintax --sintax_cutoff 0.6
    less result/raw/otus_notmatched.fa|grep '>' -c

    
    # # 原始特征表行数
    # wc -l result/raw/otutab.txt
    ## 9. 有参比对——功能预测，如Greengene，可用于picurst, bugbase分析
# # 
#      mkdir -p result/gg/
#     #与GG所有97% OTUs比对，用于功能预测
# 
#     #  #方法1. usearch比对更快，但文件超限报错选方法2
#     #  usearch -otutab temp/filtered.fa -otus ${ea}/gg/97_otus.fasta \
#     #  	-otutabout result/gg/otutab.txt -threads 12
#     # # #79.9%, 4核时8m；12核5m
#     #  head -n3 result/gg/otutab.txt
# 
#     # #方法2. vsearch比对，更准更慢，但并行更强
#      vsearch --usearch_global temp/filtered.fa --db ${ea}/gg/97_otus.fasta \
#        --otutabout result/gg/otutab.txt --id 0.97 --threads 12
#     #80.9%, 12cores 20m, 1core 1h, 594Mb
# #
# #     #统计
#      usearch -otutab_stats result/gg/otutab.txt  \
#        -output result/gg/otutab.stat
#      cat result/gg/otutab.stat
    
     # usearch -otutab_stats result/otutab.txt \
     #   -output result/otutab.stat
     # cat result/otutab.stat
     # 
    # 为去除16S rDNA测序中的非特异性扩增和按需分配污染，
    #细菌和古菌(原核生物)、以及去除叶绿体和线粒体并统计比例，输出筛选并按丰度排序的OTU表。
    #输入为特征表(result/raw/otutab.txt)和物种注释(result/raw/otus.sintax)，
    #输出筛选并排序的特征表(result/otutab.txt)、统计污染比例文件(result/raw/otutab_nonBac.txt)和过滤细节(otus.sintax.discard)。
    #注：真菌ITS数据，请改用otutab_filter_nonFungi.R脚本，只筛选注释为真菌的序列。
    # 查看脚本帮助，请运行Rscript ${ea}/script/otutab_filter_nonBac.R -h
    Rscript ${ea}/script/otutab_filter_nonBac.R \
      --input result/raw/otutab.txt \
      --taxonomy result/raw/otus.sintax \
      --output result/otutab.txt\
      --stat result/raw/otutab_nonBac.stat \
      --discard result/raw/otus.sintax.discard
    # 筛选后特征表行数，2921-2904，删除17个ASV
    wc -l result/otutab.txt
    
    # #过滤特征表对应序列
    cut -f 1 result/otutab.txt | tail -n+2 > result/otutab.id
    usearch -fastx_getseqs result/raw/otus.fa -labels result/otutab.id -fastaout result/otus.fa
    # # #过滤特征表对应序列注释
    awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}'\
        result/raw/otus.sintax result/otutab.id \
        > result/otus.sintax
    # #补齐末尾列
    sed -i 's/\t$/\td:Unassigned/' result/otus.sintax
    head -n2 result/otus.sintax
    less result/otus.fa|grep '>' -c
    wc -l result/otutab.txt
    wc -l result/otutab.id
    # wc -l result/otutab_rare.txt
    # wc -l result/taxonomy_rare.txt
    # # 方法2. 觉得筛选不合理可以不筛选
    # cp result/raw/otu* result/

    #可选统计方法：OTU表简单统计 Summary OTUs table
    usearch -otutab_stats result/otutab.txt \
      -output result/otutab.stat
    cat result/otutab.stat
    #注意最小值、分位数，或查看result/raw/otutab_nonBac.stat中样本详细数据量，用于重采样

### 5.3 等量抽样标准化 normlize by subsample

    #使用vegan包进行等量重抽样，输入reads count格式Feature表result/otutab.txt
    #可指定输入文件、抽样量和随机数，输出抽平表result/otutab_rare.txt和多样性alpha/vegan.txt
    mkdir -p result/alpha/
    Rscript ${ea}/script/otutab_rare.R --input result/otutab.txt \
      --depth 45310 --seed 123 \
      --normalize result/otutab_rare.txt \
      --output result/alpha/vegan.txt
    usearch -otutab_stats result/otutab_rare.txt \
      -output result/otutab_rare.stat
    cat result/otutab_rare.stat

    cut -f 1-2 metadata.txt > temp/group.txt

## 6. Alpha多样性 Alpha diversity

### 6.1. 计算多样性指数 Calculate alpha diversity index
    #Calculate all alpha diversity index(Chao1有错误勿用)
    #details in http://www.drive5.com/usearch/manual/alpha_metrics.html
    usearch -alpha_div result/otutab.txt \
      -output result/alpha/alpha.txt

### 6.2. 计算稀释过程的丰富度变化 Rarefaction
    #稀释曲线：取1%-100%的序列中OTUs数量，20s
    #Rarefaction from 1%, 2% .. 100% in richness (observed OTUs)-method fast / with_replacement / without_replacement https://drive5.com/usearch/manual/cmd_otutab_subsample.html
    usearch -alpha_div_rare result/otutab_rare.txt \
      -output result/alpha/alpha_rare.txt \
      -method without_replacement

### 6.3. 筛选各组高丰度菌用于比较

    #计算各特征的均值，有组再求分组均值，需根据实验设计metadata.txt修改组列名
    #输入文件为feautre表result/otutab.txt，实验设计metadata.txt
    #输出为特征表按组的均值-一个实验可能有多种分组方式
    Rscript ${ea}/script/otu_mean.R --input result/otutab_rare.txt \
      --design metadata.txt \
      --group Group --thre 1 \
      --output result/otutab_mean.txt
    head -n3 result/otutab_mean.txt

    #如以平均丰度频率高于千一(0.1%)为筛选标准，可选千五或万五，得到每个组的OTU组合
    awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=2;i<=NF;i++) a[i]=$i;} \
        else {for(i=2;i<=NF;i++) if($i>0.1) print $1, a[i];}}' \
        result/otutab_mean.txt > result/alpha/otu_group_exist.txt
    head result/alpha/otu_group_exist.txt
    # 结果可以直接在http://www.ehbio.com/ImageGP绘制Venn、upSetView和Sanky
    # 可在http://ehbio.com/test/venn/中绘图并显示各组共有和特有元素
   
## 7. Beta多样性 Beta diversity

    #结果有多个文件，需要目录
    mkdir -p result/beta/
    #基于OTU构建进化树 Make OTU tree, 30s
    usearch -cluster_agg result/otus.fa -treeout result/otus.tree
    #生成5种距离矩阵：bray_curtis, euclidean, jaccard, manhatten, unifrac, 3s
    usearch -beta_div result/otutab_rare.txt -tree result/otus.tree \
      -filename_prefix result/beta/ # 1s
    usearch -beta_div result/otutab_rare.txt -filename_prefix result/beta/
    
## 8. 物种注释及分类汇总 Taxonomy summary

     #OTU对应物种注释2列格式：去除sintax中置信值，只保留物种注释，替换:为_，删除引号
    cut -f 1,4 result/otus.sintax \
      |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' \
      > result/taxonomy2.txt
    head -n3 result/taxonomy2.txt

    #OTU对应物种8列格式：注意注释是非整齐
    #生成物种表格OTU/ASV中空白补齐为Unassigned
    awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
      split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
      result/taxonomy2.txt > temp/otus.tax
    sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | \
      sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
      > result/taxonomy.txt
    head -n3 result/taxonomy.txt

    # 统计门纲目科属，使用 rank参数 p c o f g，为phylum, class, order, family, genus缩写
    mkdir -p result/tax
    # for i in s;do
    for i in p c o f g;do
    # for i in f g s;do
      usearch -sintax_summary result/otus.sintax \
      -otutabin result/otutab_rare.txt -rank ${i} \
      -output result/tax/sum_${i}.txt
    done
    sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' result/tax/sum_*.txt
    # 列出所有文件
    ls -sh result/tax/sum_*.txt
    head -n3 result/tax/sum_g.txt


## 9. 有参比对——功能预测，如Greengene，可用于picurst, bugbase分析

    mkdir -p result/gg/
    #与GG所有97% OTUs比对，用于功能预测
    
    #方法1. usearch比对更快，但文件超限报错选方法2
    usearch -otutab temp/filtered.fa -otus ${ea}/gg/97_otus.fasta.gz \
    	-otutabout result/gg/otutab.txt -threads 12
    #79.9%, 4核时8m；12核5m
    head -n3 result/gg/otutab.txt

    # #方法2. vsearch比对，更准更慢，但并行更强
    # vsearch --usearch_global temp/filtered.fa --db ${ea}/gg/97_otus.fasta \
    #   --otutabout result/gg/otutab.txt --id 0.97 --threads 12
    #80.9%, 12cores 20m, 1core 1h, 594Mb

    #统计
    usearch -otutab_stats result/gg/otutab.txt  \
      -output result/gg/otutab.stat
    cat result/gg/otutab.stat


## 10. 空间清理及数据提交

  
    rm -rf temp/*.fq
    
    #短期不用数库压缩节省空间
    gzip ${ea}/usearch/*.fa
    gzip ${ea}/gg/*.fasta
    gzip seq/*.fq
    gzip seq/*.fastq
    # 原始序列统计md5值，用于数据提交
    cd seq
    md5sum *_1.fq.gz> md5sum1.txt
    md5sum *_2.fq.gz > md5sum2.txt
    paste md5sum1.txt md5sum2.txt | awk '{print $2"\t"$1"\t"$4"\t"$3}' | \
      sed 's/*//g' > ../result/md5sum.txt
    rm md5sum*
    cd ..
    cat result/md5sum.txt



# 23、R语言多样性和物种分析


## 1. Alpha多样性

### 1.1 Alpha多样性箱线图
    # 查看帮助(ANOVA+LSD计算，符合正态分布和方差齐性才能进行计算)
    Rscript ${ea}/script/alpha_boxplot.R -h
    # 完整参数，多样性指数可选richness chao1 ACE shannon simpson invsimpson
    Rscript ${ea}/script/alpha_boxplot.R --alpha_index richness \
      --input result/alpha/vegan.txt --design metadata.txt \
      --group Group --output result/alpha/ \
      --width 89 --height 59
    # 使用循环绘制6种常用指数
    mkdir -p result/alpha/
     for i in `head -n1 result/alpha/vegan.txt|cut -f 2-`;do
      Rscript ${ea}/script/alpha_boxplot.R --alpha_index ${i} \
        --input result/alpha/vegan.txt --design metadata.txt \
        --group Group --output result/alpha/ \
        --width 139 --height 109;done
      for i in `head -n1 result/alpha/alpha.txt|cut -f 2-`;do
      Rscript ${ea}/script/alpha_boxplot.R --alpha_index ${i} \
        --input result/alpha/alpha.txt --design metadata.txt \
        --group Group --output result/alpha/ \
        --width 139 --height 109;done

    # Alpha多样性柱状图+标准差
    for i in `head -n1 result/alpha/alpha.txt|cut -f 2-`;do
    Rscript ${ea}/script/alpha_barplot.R --alpha_index ${i} \
      --input result/alpha/alpha.txt --design metadata.txt \
      --group Group --output result/alpha/ \
      --width 89 --height 59;done
    for i in `head -n1 result/alpha/vegan.txt|cut -f 2-`;do
    Rscript ${ea}/script/alpha_barplot.R --alpha_index ${i} \
      --input result/alpha/vegan.txt --design metadata.txt \
      --group Group --output result/alpha/ \
      --width 89 --height 59;done
   
### 1.2 稀释曲线
    Rscript ${ea}/script/alpha_rare_curve.R \
      --input result/alpha/alpha_rare.txt --design metadata.txt \
      --group Group --output result/alpha/ \
      --width 139 --height 109

### 1.3 多样性维恩图
    # 三组比较:-f输入文件,-a/b/c/d/g分组名,-w/u为宽高英寸,-p输出文件名后缀
    bash ${ea}/script/sp_vennDiagram.sh \
      -f result/alpha/otu_group_exist.txt \
      -a WT -b KO -c OE \
      -w 3 -u 3 \
      -p WT_KO_OE
    # 四组比较，图和代码见输入文件目录，运行目录为当前项目根目录
    bash ${ea}/script/sp_vennDiagram.sh \
      -f result/alpha/otu_group_exist.txt \
      -a WT -b KO -c OE -d All \
      -w 3 -u 3 \
      -p WT_KO_OE_All

## 2. Beta多样性

### 2.1 距离矩阵热图pheatmap
    # 添加分组注释，如2，4列的基因型和地点
    cut -f 1-2 metadata.txt > temp/group.txt
    # 以bray_curtis为例，-f输入文件,-h是否聚类TRUE/FALSE,-u/v为宽高英寸
    for i in weighted_unifrac unifrac bray_curtis jaccard;do
    bash ${ea}/script/sp_pheatmap.sh \
      -f result/beta/${i}.txt \
      -H 'TRUE' -u 5 -v 5
    # -P添加行注释文件，-Q添加列注释
    bash ${ea}/script/sp_pheatmap.sh \
      -f result/beta/${i}.txt \
      -H 'TRUE' -u 8.9 -v 5.6 \
      -P temp/group.txt -Q temp/group.txt;done
    # 距离矩阵与相关类似，可尝试corrplot或ggcorrplot绘制更多样式
    # - [绘图相关系数矩阵corrplot](http://mp.weixin.qq.com/s/H4_2_vb2w_njxPziDzV4HQ)
    # - [相关矩阵可视化ggcorrplot](http://mp.weixin.qq.com/s/AEfPqWO3S0mRnDZ_Ws9fnw)

### 2.2 主坐标分析PCoA
    # 输入文件，选择分组，输出文件，图片尺寸mm，统计见beta_pcoa_stat.txt
    for i in weighted_unifrac unifrac bray_curtis;do
    Rscript ${ea}/script/beta_pcoa.R \
      --input result/beta/${i}.txt --design metadata.txt \
      --group Group --output result/beta/${i}.pcoa.pdf \
      --width 139 --height 109;done
    
### 2.3 限制性主坐标分析CPCoA
    for i in wweighted_unifrac unifrac bray_curtis;do
    Rscript ${ea}/script/beta_cpcoa.R \
      --input result/beta/${i}.txt --design metadata.txt \
      --group Group --output result/beta/${i}.cpcoa.pdf \
      --width 139 --height 109;done

## 3. 物种组成Taxonomy

### 3.1 堆叠柱状图Stackplot
    # 以门(p)水平为例，结果包括output.sample/group.pdf两个文件
    Rscript ${ea}/script/tax_stackplot.R \
      --input result/tax/sum_p.txt --design metadata.txt \
      --group Group --output result/tax/sum_p.stackplot \
      --legend 5 --width 89 --height 59
    # 批量绘制输入包括p/c/o/f/g共5级
    for i in k p c o f g s; do
    Rscript ${ea}/script/tax_stackplot.R \
      --input result/tax/sum_${i}.txt --design metadata.txt \
      --group Group --output result/tax/sum_${i}.stackplot \
      --legend 838 --width 909 --height 209; done

## 24-2. STAMP输入文件准备

### 2.1 命令行生成输入文件
    Rscript ${ea}/script/format2stamp.R -h
    mkdir -p result/stamp
    Rscript ${ea}/script/format2stamp.R --input result/otutab_rare.txt \
      --taxonomy result/taxonomy.txt --threshold 0 \
      --output result/stamp/tax
    

### 2.2 Rmd生成输入文件
    #1. 24compare/stamp目录中准备otutab.txt和taxonomy.txt文件；
    #2. Rstudio打开format2stamp.Rmd，设置参数；
    #3. 点击Knit在当前目录生成stamp输入文件和可重复计算网页。


## 24-3. LEfSe输入文件准备

### 3.1. 命令行生成文件

    # 可选命令行生成输入文件
    Rscript ${ea}/script/format2lefse.R -h
    mkdir -p result/lefse
    Rscript ${ea}/script/format2lefse.R --input result/otutab_rare.txt \
      --taxonomy result/taxonomy.txt --design metadata.txt \
      --group Group --threshold 0 \
      --output result/lefse/LEfSe

   
### 3.2 Rmd生成输入文件
    #1. 24Compare/LEfSe目录中准备otutab.txt, metadata.txt, taxonomy.txt三个文件；
    #2. Rstudio打开format2lefse.Rmd并Knit生成输入文件和可重复计算网页；

### 3.3 LEfSe分析
    #方法1. 打开LEfSe.txt并在线提交 http://www.ehbio.com/ImageGP/index.php/Home/Index/LEFSe.html
    #方法2. LEfSe本地分析(限Linux服务器、选学)，参考代码见附录
    #方法3. LEfSe官网在线使用



# 25、QIIME 2分析流程
    # 代码详见 25QIIME2/pipeline_qiime2.sh
    #外部导入特征表和代表序列(常用)

	# # 上传其他流程生成的OTU表otutab.txt和代表序列otus.fa
	# wget -c http://210.75.224.110/github/MicrobiomeProtocol/e2.QIIME2/otutab.txt
	# wget -c http://210.75.224.110/github/MicrobiomeProtocol/e2.QIIME2/otus.fa
	# # 转换文本为Biom1.0，注意biom --version 2.1.5/8可以，2.1.7报错
	biom convert -i otutab.txt -o otutab.biom --table-type="OTU table" --to-json
	# 导入特征表，9s
	qiime tools import --input-path otutab.biom --type 'FeatureTable[Frequency]' --input-format BIOMV100Format --output-path table.qza
	# 导入代表序列，8s
	qiime tools import --input-path otus.fa --type 'FeatureData[Sequence]' --output-path rep-seqs.qza

  # 统计特征表。
  qiime feature-table summarize \
    --i-table table.qza \
    --o-visualization table.qzv \
    --m-sample-metadata-file metadata.txt
  
  # 统计代表序列。
  qiime feature-table tabulate-seqs \
    --i-data rep-seqs.qza \
    --o-visualization rep-seqs.qzv
  
  ## 1. Alpha和beta多样性分析

## 构建进化树用于多样性分析 53s
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza


### 计算核心多样性
	# 13s，采样深度通常选择最小值，来自table.qzv##
	qiime diversity core-metrics-phylogenetic \
	  --i-phylogeny rooted-tree.qza \
	  --i-table feature-frequency-filtered-table.qza \
	  --p-sampling-depth 45310 \
	  --m-metadata-file metadata.txt \
	  --output-dir core-metrics-results

### Alpha多样性组间显著性分析和可视化
	# 7s, 可选的alpha指数有 faith_pd、shannon、observed_features、evenness
	index=shannon
	qiime diversity alpha-group-significance \
	  --i-alpha-diversity core-metrics-results/${index}_vector.qza \
	  --m-metadata-file metadata.txt \
	  --o-visualization core-metrics-results/${index}-group-significance.qzv

	index=faith_pd
	qiime diversity alpha-group-significance \
	  --i-alpha-diversity core-metrics-results/${index}_vector.qza \
	  --m-metadata-file metadata.txt \
	  --o-visualization core-metrics-results/${index}-group-significance.qzv

	index=observed_features
	qiime diversity alpha-group-significance \
	  --i-alpha-diversity core-metrics-results/${index}_vector.qza \
	  --m-metadata-file metadata.txt \
	  --o-visualization core-metrics-results/${index}-group-significance.qzv

	index=evenness
	qiime diversity alpha-group-significance \
	  --i-alpha-diversity core-metrics-results/${index}_vector.qza \
	  --m-metadata-file metadata.txt \
	  --o-visualization core-metrics-results/${index}-group-significance.qzv

### Beta多样性组间显著性分析和可视化
	# 可选的beta指数有 unweighted_unifrac、bray_curtis、weighted_unifrac和jaccard
	# 7s, 指定分组是减少计算量，置换检验较耗时

	distance=weighted_unifrac
	column=Group
	qiime diversity beta-group-significance \
	  --i-distance-matrix core-metrics-results/${distance}_distance_matrix.qza \
	  --m-metadata-file metadata.txt \
	  --m-metadata-column ${column} \
	  --o-visualization core-metrics-results/${distance}-${column}-significance.qzv \
	  --p-pairwise

distance=unweighted_unifrac
	column=Group
	qiime diversity beta-group-significance \
	  --i-distance-matrix core-metrics-results/${distance}_distance_matrix.qza \
	  --m-metadata-file metadata.txt \
	  --m-metadata-column ${column} \
	  --o-visualization core-metrics-results/${distance}-${column}-significance.qzv \
	  --p-pairwise

distance=bray_curtis
	column=Group
	qiime diversity beta-group-significance \
	  --i-distance-matrix core-metrics-results/${distance}_distance_matrix.qza \
	  --m-metadata-file metadata.txt \
	  --m-metadata-column ${column} \
	  --o-visualization core-metrics-results/${distance}-${column}-significance.qzv \
	  --p-pairwise

 #导出table
    qiime tools export \
      --input-path core-metrics-results/weighted_unifrac_distance_matrix.qza \
      --output-path weighted_unifrac_distance_matrix


