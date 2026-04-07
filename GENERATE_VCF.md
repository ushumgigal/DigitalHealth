generate_vcf.py -- Sentetik Genomik Pipeline
=============================================

Bu script, sentetik bir hasta için FASTQ -> BAM -> VCF pipeline'ının tamamını uçtan uca çalıştırır. Rastgele bir referans genom üretir, üzerine SNP ve küçük indel varyantları yerleştirerek iki haplotipi olan (diploid) sentetik bir birey oluşturur, bu bireyden paired-end FASTQ read'leri simüle eder, bunları referansa hizalar ve varyant çağrısı (variant calling) yaparak geçerli bir VCF dosyası üretir.

Script tek bir hasta için kullanılabileceği gibi, `--auto` parametresiyle 3 klinik x 10 hasta olmak üzere toplam 30 hastanın VCF'lerini otomatik olarak üretebilir. `--joint` parametresi ise aynı işlemi yapmanın yanı sıra her kliniğin 10 hastasını tek bir çoklu-örnek (multi-sample) VCF dosyasında birleştirir.

Tüm varyantlar **heterozigot (0/1)** olarak üretilir. Yani her varyant yalnızca bir haplotipte (hap2) bulunur; diğer haplotip (hap1) referans alelini taşır.


Gerekli Araçlar
----------------

Script, harici biyoinformatik araçlarını `subprocess` üzerinden çağırmaktadır. Çalıştırmadan önce aşağıdaki araçların kurulu olması gerekmektedir:

| Araç         | Görevi                            | Minimum Versiyon |
|--------------|-----------------------------------|------------------|
| **bwa**      | Read hizalama (alignment)         | 0.7.x            |
| **samtools** | SAM/BAM dosya işleme              | 1.10+            |
| **bcftools** | Varyant çağrısı (variant calling) | 1.10+            |

### macOS Kurulumu

```bash
brew install bwa samtools bcftools
```

### Ubuntu / Debian Kurulumu

```bash
sudo apt-get install bwa samtools bcftools
```

Kurulumu doğrulamak için:

```bash
bwa 2>&1 | head -3
samtools --version
bcftools --version
```


Kullanım
--------

### Tek Hasta (Varsayılan)

```bash
python3 scripts/generate_vcf.py
```

Bu komut `pipeline_output/` dizinine tek bir hastanın tüm ara dosyaları ve nihai VCF dosyasını üretir.

### Tek Hasta (Özelleştirilmiş)

```bash
python3 scripts/generate_vcf.py \
    --outdir klinik1/raw_variants/patient1 \
    --sample-name patient1 \
    --num-variants 20 \
    --seed 100
```

### Toplu Üretim (Auto Mod)

```bash
python3 scripts/generate_vcf.py --auto
```

Bu komut aşağıdaki yapıyı oluşturur:

```
klinik1/raw_variants/patient1.vcf  ... patient10.vcf
klinik2/raw_variants/patient1.vcf  ... patient10.vcf
klinik3/raw_variants/patient1.vcf  ... patient10.vcf
```

Auto modda ara dosyalar (FASTQ, BAM, referans) geçici dizinlerde oluşturulup işlem bittikten sonra otomatik olarak silinir. Yalnızca nihai VCF dosyaları ilgili klinik dizinine kopyalanır.

Auto mod diğer parametrelerle birleştirilebilir:

```bash
python3 scripts/generate_vcf.py --auto --num-variants 30 --num-reads 20000
```

### Toplu Üretim + Birleştirme (Joint Mod)

```bash
python3 scripts/generate_vcf.py --joint
```

`--auto` ile aynı işlemi yapar (3 klinik x 10 hasta), ek olarak her klinikteki 10 hastanın VCF'lerini `bcftools merge` ile birleştirerek tek bir çoklu-örnek VCF oluşturur:

```
klinik1/
  klinik1.vcf                          <-- 10 hastanın birleşik VCF'i
  raw_variants/
    patient1.vcf ... patient10.vcf     <-- bireysel VCF'ler (korunur)

klinik2/
  klinik2.vcf
  raw_variants/
    patient1.vcf ... patient10.vcf

klinik3/
  klinik3.vcf
  raw_variants/
    patient1.vcf ... patient10.vcf
```

Joint modda aynı klinikteki tüm hastalar **ortak bir referans genom** kullanır. Bu, `bcftools merge` için zorunlu bir koşuldur. Farklı klinikler farklı referanslara sahip olabilir.

Joint mod da diğer parametrelerle birleştirilebilir:

```bash
python3 scripts/generate_vcf.py --joint --num-variants 25 --seed 99
```


Parametreler
------------

| Parametre          | Varsayılan            | Açıklama                                                                                                  |
|--------------------|-----------------------|-----------------------------------------------------------------------------------------------------------|
| `--outdir`         | `pipeline_output`     | Çıktı dizini. Tek hasta modunda tüm dosyalar (FASTQ, BAM, VCF) bu dizine yazılır.                        |
| `--ref-length`     | `50000`               | Sentetik referans genomun uzunluğu (baz çifti). Daha büyük değerler daha gerçekçi ama daha yavaş üretim.  |
| `--num-variants`   | `15`                  | Yerleştirilecek toplam varyant sayısı (SNP + insertion + deletion karışımı).                               |
| `--num-reads`      | `10000`               | Simüle edilecek paired-end read çifti sayısı. Daha fazla read = daha yüksek kapsam (coverage).            |
| `--read-length`    | `150`                 | Her bir read'in uzunluğu (baz çifti). Illumina standardı olan 150 bp varsayılandır.                       |
| `--fragment-size`  | `400`                 | Ortalama fragment boyutu (baz çifti). Gerçekte Illumina kütüphanelerinde ~300-500 bp arasındadır.         |
| `--error-rate`     | `0.005`               | Baz başına sekanslama hatası olasılığı (%0.5). Illumina platformu için gerçekçi bir değerdir.             |
| `--seed`           | `42`                  | Rastgelelik tohumu. Aynı seed aynı çıktıyı üretir (reproducibility). Auto modda her hasta için türetilir. |
| `--sample-name`    | `SYNTH_PATIENT_01`    | BAM header ve VCF örnek kolonunda görünecek isim. Auto modda otomatik olarak `patient1`..`patient10` atanır.|
| `--auto`           | *(kapalı)*            | Aktif edildiğinde 3 klinik x 10 hasta = 30 VCF üretir. `--outdir` ve `--sample-name` yok sayılır.        |
| `--joint`          | *(kapalı)*            | `--auto` ile aynı + her klinik için birleşik çoklu-örnek VCF üretir (`klinikN.vcf`). Aynı klinikteki hastalar ortak referans genom kullanır. |


Pipeline Aşamaları
------------------

Script aşağıdaki aşamaları sırasıyla çalıştırır:

### 1. Referans Genom Üretimi

Belirtilen uzunlukta (`--ref-length`) rastgele bir DNA dizisi (A, C, G, T) üretilir ve `reference.fa` dosyasına FASTA formatında yazılır. Tek bir kontig (`chr_synth`) içerir. Bu referans, hizalama ve varyant çağrısı için kullanılır.

### 2. Varyant Yerleştirilmesi ve Haplotip Oluşturma

Referans genom üzerine rastgele konumlara varyantlar yerleştirilir:

- **SNP** (Single Nucleotide Polymorphism): Tek baz değişimi (örn. A -> G)
- **Insertion**: 1-3 baz eklenmesi (örn. C -> CGA)
- **Deletion**: 1-3 baz silinmesi (örn. TCA -> T)

Tüm varyantlar **0/1 (heterozigot)** olarak atanır: varyant yalnızca hap2'de bulunur, hap1 referans alelini taşır. Bu sayede read'lerin yaklaşık yarısı referans, yarısı alternatif aleli taşır.

Varyant dağılımı yaklaşık olarak %55 SNP, %25 insertion, %20 deletion şeklindedir.

### 3. Paired-End Read Simülasyonu

Diploid genomdan (hap1 + hap2) paired-end read'ler üretilir:

- Her read çifti iki haplotipin birinden %50 olasılıkla seçilir
- Fragment boyutu normal dağılımdan örneklenir (ortalama: `--fragment-size`, standart sapma: %15)
- R1 (forward read) ve R2 (reverse-complement read) olarak iki FASTQ dosyasına yazılır
- Her bazda `--error-rate` olasılığıyla sekanslama hatası eklenir
- Kalite skorları (Phred) gerçekçi biçimde üretilir: read'in 3' ucuna doğru kalite düşmesi Illumina platformuna özgü bir özelliktir

### 4. Hizalama (Alignment)

Üretilen FASTQ dosyaları `bwa mem` ile referans genoma hizalanır, ardından `samtools sort` ile koordinata göre sıralanır ve `samtools index` ile indekslenir:

```
bwa mem reference.fa reads_R1.fastq reads_R2.fastq | samtools sort -o aligned.sorted.bam
samtools index aligned.sorted.bam
```

BAM dosyasına `@RG` (Read Group) bilgisi eklenir; böylece örnek adı BAM'dan VCF'e taşınır.

### 5. Varyant Çağrısı (Variant Calling)

`bcftools mpileup` ile hizalanmış read'lerden pileup oluşturulur, ardından `bcftools call` ile varyantlar çağrılır:

```
bcftools mpileup -f reference.fa aligned.sorted.bam | bcftools call -mv --ploidy 2 -Oz -o variants.vcf.gz
```

Çıktı VCF dosyası sıkıştırılmış halden açık metin formatına dönüştürülür.

### 6. Özet (Yalnızca Tek Hasta Modunda)

Tek hasta modunda pipeline sonunda bir özet yazdırılır:

```
==================================================
  PIPELINE SUMMARY
==================================================
  Total variant records : 15
  SNPs                  : 7
  Indels                : 8
  Genotypes             : {'0/1': 15, '1/1': 0, '0/0': 0}
==================================================
```


Çıktı Dosyaları
----------------

### Tek Hasta Modu

`--outdir` ile belirtilen dizinde aşağıdaki dosyalar oluşturulur:

| Dosya                    | Format | Açıklama                                                                         |
|--------------------------|--------|----------------------------------------------------------------------------------|
| `reference.fa`           | FASTA  | Sentetik referans genom. Tek kontig (`chr_synth`).                               |
| `reference.fa.fai`       | FAI    | samtools faidx indeksi.                                                          |
| `reference.fa.bwt` vb.  | BWA    | bwa indeks dosyaları (`.bwt`, `.pac`, `.ann`, `.amb`, `.sa`).                    |
| `reads_R1.fastq`        | FASTQ  | Forward (R1) read'ler. Her kayıt: isim, dizi, kalite skoru.                     |
| `reads_R2.fastq`        | FASTQ  | Reverse (R2) read'ler. R1 ile eşleştirilmiş (paired).                           |
| `aligned.sorted.bam`    | BAM    | Koordinata göre sıralanmış hizalama dosyası. Read group bilgisi içerir.          |
| `aligned.sorted.bam.bai`| BAI    | BAM indeks dosyası. IGV gibi görüntüleyiciler için gereklidir.                   |
| `variants.vcf`          | VCF    | Nihai varyant çağrısı dosyası. SNP ve indel'leri genotip bilgisiyle içerir.      |
| `variants.vcf.gz`       | VCF.GZ | Sıkıştırılmış VCF (bcftools çıktısı).                                            |
| `variants.vcf.gz.csi`   | CSI    | Sıkıştırılmış VCF indeksi.                                                       |

### Auto Modu

Auto modda yalnızca VCF dosyaları kalıcı olarak saklanır:

```
klinik1/raw_variants/patient1.vcf
klinik1/raw_variants/patient2.vcf
...
klinik1/raw_variants/patient10.vcf
klinik2/raw_variants/patient1.vcf
...
klinik3/raw_variants/patient10.vcf
```

Ara dosyalar (FASTQ, BAM, referans) geçici dizinlerde oluşturulup işlendi̇kten sonra otomatik olarak temizlenir.

### Joint Modu

Joint modda yukarıdaki bireysel VCF'lere ek olarak her klinik dizininde birleşik bir VCF oluşturulur:

```
klinik1/klinik1.vcf   <-- 10 hastanın tüm varyantları, tek dosyada
klinik2/klinik2.vcf
klinik3/klinik3.vcf
```

Birleşik VCF'de her varyant satırı tüm 10 hastanın genotipini içerir. Bir hastada varyant yoksa o hasta için `0/0` (homozigot-referans) veya `./.` (veri yok) görünür.


VCF Çıktısı Nasıl Okunur
-------------------------

Üretilen VCF dosyası standart VCF 4.2 formatındadır. Bir örnek satır:

```
#CHROM     POS    ID  REF  ALT  QUAL   FILTER  INFO                 FORMAT           patient1
chr_synth  5318   .   C    G    222.2  .       DP=62;...;AC=1;AN=2  GT:PL:DP:AD      0/1:255,0,247:62:24,38
```

### Sütunlar

| Sütun      | Açıklama                                                                                    |
|------------|---------------------------------------------------------------------------------------------|
| `CHROM`    | Kromozom adı. Bu pipeline'da her zaman `chr_synth`.                                         |
| `POS`      | Varyant pozisyonu (1-tabanlı).                                                              |
| `ID`       | Varyant tanımlayıcısı (bu pipeline'da `.` yani tanımsız).                                   |
| `REF`      | Referans alel. SNP için tek baz, deletion için birden fazla baz.                            |
| `ALT`      | Alternatif alel. SNP için tek baz, insertion için birden fazla baz.                         |
| `QUAL`     | Varyant kalite skoru (Phred ölçeğinde). Yüksek = daha güvenilir.                            |
| `FILTER`   | Filtre durumu. `.` = filtre uygulanmamış.                                                   |
| `INFO`     | Varyant hakkında ek bilgiler (aşağıya bakınız).                                             |
| `FORMAT`   | Örnek alanlarının formatı.                                                                  |
| Örnek      | Hasta için genotip ve diğer bilgiler.                                                       |

### Önemli INFO Alanları

| Alan   | Açıklama                                                              |
|--------|-----------------------------------------------------------------------|
| `DP`   | Toplam okuma derinliği (read depth).                                  |
| `AC`   | Alternatif alel sayısı. Tüm varyantlar heterozigot olduğu için her zaman 1. |
| `AN`   | Toplam alel sayısı (diploid için her zaman 2).                        |
| `MQ`   | Ortalama hizalama kalitesi (mapping quality).                         |
| `INDEL`| Bu alanın varlığı kaydın bir indel olduğunu belirtir.                 |

### Önemli FORMAT Alanları

| Alan  | Açıklama                                                                                     |
|-------|----------------------------------------------------------------------------------------------|
| `GT`  | Genotip. `0/0` = homozigot-referans, `0/1` = heterozigot, `1/1` = homozigot-alternatif.     |
| `PL`  | Phred ölçekli genotip olasılıkları (0/0, 0/1, 1/1 sırasıyla). Düşük = daha olası.          |
| `DP`  | Bu örnekteki okuma derinliği.                                                                |
| `AD`  | Alel derinliği. Örn: `24,38` = referans alel 24 read, alternatif alel 38 read.              |


Tekrarlanabilirlik (Reproducibility)
-------------------------------------

Pipeline tamamen deterministiktir. Aynı `--seed` değeri aynı çıktıyı üretir. Auto modda her hasta için benzersiz bir seed hesaplanır:

```
seed = base_seed + (klinik_indeksi * 1000) + (hasta_numarasi * 100)
```

Örneğin `--seed 42` ile:

| Hasta                | Seed |
|----------------------|------|
| klinik1 / patient1   | 142  |
| klinik1 / patient2   | 242  |
| klinik2 / patient1   | 1142 |
| klinik3 / patient10  | 3042 |


Teknik Notlar
-------------

- **Diploid simülasyon**: İki ayrı haplotip dizisi oluşturulur. Tüm varyantlar heterozigot (`0/1`) olduğu için yalnızca hap2'ye uygulanır; hap1 referans ile aynıdır. Read'ler %50-%50 olasılıkla iki haplotipin birinden örneklenir.
- **Fragment boyutu varyasyonu**: Fragment boyutu sabit değil, normal dağılımdan (sd = %15) örneklenir. Bu, gerçekçi Illumina kütüphanelerini taklit etmenin yanı sıra, bwa'nın indel içeren read'leri "proper pair" olarak işaretlemesi için kritiktir. Sabit fragment boyutunda, indel'in neden olduğu 1-3 bp'lik insert boyutu kayması, bwa tarafından "improper pair" olarak değerlendirilir ve bcftools bu read'leri yok sayar.
- **Kalite skoru profili**: Read'lerin 3' ucuna doğru kalite düşmesi uygulanır (Illumina'ya özgü). Hata eklenmiş bazlara düşük kalite skoru atanır.
- **Varyant aralığı**: Varyantlar birbirinden en az 10 bp uzakta yerleştirilir; referansın her iki ucunda 200 bp'lik bir marj bırakılır.