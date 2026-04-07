ÖN BİLGİ
========
Genomik veriler (ör. VCF dosyaları) birey hakkında son derece hassas bilgiler içerir: hastalık yatkınlıkları, kalıtsal özellikler, hatta akrabalık ilişkileri. Bu nedenle bu verilerin merkezi bir sunucuda toplanması hem etik hem de hukuki (KVKK, GDPR benzeri regülasyonlar) riskler doğurur.

Klasik makine öğrenmesi yaklaşımında veri merkezileştirilir; ancak sağlık sektöründe veri genellikle:
- Dağıtık (farklı hastaneler/laboratuvarlar)
- Heterojen (farklı cihazlar, farklı popülasyonlar)
- Paylaşımı kısıtlıdır (hukuki/etik sebeplerden ötürü.)

Federated Learning (FL), bu problemi çözmek için geliştirilmiş bir paradigmadır. Model eğitimi veri üzerinde değil, veri bulunduğu yerde yapılır. Böylece ham veri paylaşılmadan kolektif öğrenme sağlanır.

AMAÇ
====
Bu çalışmanın amacı, dağıtık hasta verileri üzerinde gizliliği koruyarak model eğitimi yapabilen bir Federated Learning prototipi geliştirmektir.

Alt hedefler:
- VCF veya anotasyonlanmış varyantlardan anlamlı feature çıkarımı
- Birden fazla “client” (örn. hastane) üzerinde lokal model eğitimi
- Dağıtık eğitim sonucu elde edilen gradient değerleri kullanılarak merkezi sunucuda model ağırlıklarının birleştirilmesi (aggregation)
- Ham veriyi paylaşmadan performanslı bir global model elde edilmesi

SENTETİK VERİ
=============

Bu prototip, gerçek hasta genomik verisi yerine **sentetik veri** ile çalışmaktadır. Sentetik veri; rastgele üretilmiş bir referans genom üzerine yapay varyantlar (SNP, insertion, deletion) yerleştirilerek, gerçekçi bir FASTQ -> BAM -> VCF pipeline'ından geçirilerek elde edilir. Üretilen VCF dosyaları format ve içerik olarak gerçek sekanslama verisiyle aynı yapıdadır; ancak herhangi bir gerçek bireye ait bilgi içermez.

Sentetik veri kullanımının gerekçeleri:
- Gerçek genomik veriye erişim hukuki ve etik kısıtlamalar nedeniyle prototip aşamasında mümkün değildir
- Pipeline'ın ve federated learning mimarisinin doğrulanması için gerçek veri şart değildir
- Tekrarlanabilir (reproducible) ve kontrollü deneyler yapılabilir
- Varyant sayısı, hasta sayısı, klinik sayısı gibi parametreler serbestçe ayarlanabilir

Sentetik Veri Üretimi
---------------------

Sentetik veri üretimi `scripts/generate_vcf.py` scripti ile sağlanır. Detaylı kullanım bilgisi için bkz. [GENERATE_VCF.md](GENERATE_VCF.md).

VERİ ORGANİZASYONU
==================

Dağıtık Birimler

----------------
Her dağıtık birim (örn. klinik) kendi lokal veri dizinine sahiptir ve veri **ham (VCF)** ve **işlenmiş (annotated)** olarak iki seviyede tutulur.

Dizin yapısı:
```
klinik1/
├──── klinik1.vcf              <-- Birleşik (joint) multi-sample VCF
├──── raw_variants/
|   └── patient1.vcf
|   └── patient2.vcf
|   └── ...
│  
└──── annotated_variants/
    └── patient1.csv
    └── patient2.csv
    └── ...
```

### Veri katmanları

#### `klinikN.vcf`
- Kliniğe ait tüm hastaların varyantlarını içeren birleşik (joint) multi-sample VCF
- Her satır bir varyant pozisyonunu, her sütun bir hastanın genotipini temsil eder
- Annotation ve eğitim aşamalarında girdi olarak bu dosya kullanılır
- Bu veri **asla client dışına çıkmaz**

#### 1. `raw_variants/`
- Hasta bazında bireysel ham genomik veri
- Format: VCF (Variant Call Format)
- Bu veri **asla client dışına çıkmaz**

#### 2. `annotated_variants/`
- Annotation sonrası üretilmiş veri
- Format: CSV (veya benzeri tablo formatı)
- İçerik (örnek):

! Bu veriler **asla dağıtık birim dışına çıkmaz**.

Global Birim
------------

Federated Learning mimarisinde merkezi birim (`global/`), modelin versiyonlarını ve client’lardan gelen güncellemeleri (gradient / weight delta) yönetir.

Dizin yapısı:
```
global/
└──── gradients/
|   └── v1_0_0_klinik1.delta
|   └── v1_0_1_klinik2.delta
|   └── ...
│  
└──── versions/
    └── v1_0_0.weights
    └── v1_0_1.weights
    └── ...
```

---

### 1. `gradients/`

- Client’lardan gelen model güncellemelerini içerir
- Her dosya:
  - belirli bir model versiyonuna ait
  - belirli bir client (klinik) tarafından üretilmiş
şekilde isimlendirilir. Örneğin 1.0.0 versiyonu üzerine uygulanacak, klinik1 tarafından sağlanan gradient ismi: `v1_0_0_klinik1.delta`
- Model weight farkları (ΔW) veya gradient’ler
- Ham veri içermez
- Sadece modelin nasıl değiştiğini temsil eder

### 2. `versions/`

- Global modelin farklı versiyonlarını saklar
- Her dosya tam model ağırlıklarını içerir (full weights)
- Versiyon ismi `semantic versioning` usulune uygun olarak oluşturulur.

İŞLEYİŞ
=======

Veri Toplama
------------
Her dağıtık birim (hastane/laboratuvar) kendi lokal verisini işler.  

- Bireysel hasta VCF dosyaları lokal olarak saklanır.  
- Bu VCF’ler birleştirilerek kliniğe ait tek bir joint multi-sample VCF (`klinikN.vcf`) oluşturulur.
- Joint VCF üzerinden annotation araçları (VEP, ANNOVAR vb.) ile anotasyon yapılır.  
- Anotasyonlanmış CSV veya feature vektörleri oluşturulur.  
- Ham veri **asla dışarı çıkmaz**; sadece eğitim için gerekli feature’lar lokal olarak hazırlanır.

Eğitim
-------
### Dağıtık Birimler
- Her birim, merkezi global modelin güncel kopyasını alır.  
- Joint VCF'den (`klinikN.vcf`) üretilen anotasyonlanmış veri üzerinde lokal eğitim yapar.  
- Eğitim sırasında **sadece model ağırlıklarındaki değişiklikler veya gradientler** hesaplanır. Bu işlem prototip için oluşturulan `train` Python script'i ile sağlanır.  
- Bu güncellemeler merkezi global birime iletilir, ham veri paylaşılmaz.

### Global Birim
- Gelen güncellemeler mevcut model ile birleştirilir. Bu işlem prototip için oluşturulan `merge` Python script'i ile sağlanır.
- Oluşturulan yeni model bir sonraki eğitim turu için dağıtık birimlere gönderilir.  

Kullanım
--------
- Dağıtık birimler, merkezi global modelin en güncel sürümünü çeker ve kendi lokal verisi üzerinde tahmin veya analiz yapabilir.  
- Global model, veri paylaşımı olmadan tüm birimlerden öğrenilmiş bilgiyi taşır.  
- Sistem, sürekli güncellenebilir ve performans/ doğruluk iyileştirmeleri için yeni eğitim turları eklenebilir.

NOTLAR
======

- Dağıtık birimlerin ve global birimin farklı klasörlerde tutulması, dahili veri saklama ve yazılım çağrıları prototipin tanıtımı amacı ile temsilen oluşturulup, pilot uygulamalarda ve nihai saha faaliyetlerinde verilerin farklı cihazlarda ve fiziksel merkezlerde saklanması amaçlanmıştır. Benzer şekilde, dağıtık birimler dahilindeki veri organizasyonu, dosya formatları, kullanılacak yazılım vb işlemlerin homojen olması beklenmemektedir. Prototipin uygulamaya dökümündeki koşullar ve ön kabuller:
  - Ham verinin asla bir dağıtık birim dışında çıkmaması,
  - Dağıtık birimlerin ortak bir model için lokal olarak eğitim işlemlerini gerçekletirmeleri,
  - Bir dağıtık birimden global birime sadece lokal eğitimler ile elde edilen `gradient` setlerin iletilmesinden ibarettir.