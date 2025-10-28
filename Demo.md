

---

# **本地端 BLAST+ 完整實戰教學：從零到完成分析**

## **一、 前言**

本教學旨在引導初學者在 Linux 環境下，從頭開始安裝 NCBI BLAST+ 工具，並完成一次完整的序列比對分析。我們將學習如何：

1. 下載並設定 BLAST+ 執行環境。
    
2. 從 NCBI 下載序列檔案。
    
3. 建立一個小型的客製化本地資料庫。
    
4. 執行 `blastn` 進行核酸序列比對並解讀結果。
    

本教學所有指令皆可直接複製貼上至您的終端機 (Terminal) 執行。

---

## **二、 Part 1: 安裝與環境設定**

首先，我們需要準備好 BLAST+ 這個分析工具。

#### **步驟 1.1：建立工作目錄**

為了保持整潔，我們先建立一個專案資料夾並進入。

猛擊

```
mkdir my_blast_project
cd my_blast_project
```

#### **步驟 1.2：下載 BLAST+ 程式**

從 NCBI 官方 FTP 伺服器下載最新版的 BLAST+ for Linux。

猛擊

```
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-x64-linux.tar.gz
```

#### **步驟 1.3：解壓縮檔案**

將下載回來的壓縮檔解開。

猛擊

```
tar -zxvf ncbi-blast-2.17.0+-x64-linux.tar.gz
```

#### **步驟 1.4：設定環境變數 (關鍵步驟)**

為了讓系統能在任何路徑下都能找到 `blastn` 等指令，我們需要將 BLAST+ 的路徑加入到系統的 `PATH` 環境變數中。

**注意**：請根據你**當前所在的路徑**修改指令。以下指令假設你在 `/workspaces/BLAST-` 或類似的目錄下。你可以用 `pwd` 指令查看當前完整路徑。

猛擊
將底下這行指令中的 "/workspaces/BLAST-" 換成你用 pwd 指令看到的真實路徑

```
export PATH="/workspaces/BLAST-/my_blast_project/ncbi-blast-2.17.0+/bin:$PATH"
```

- **說明**：這個指令告訴系統：「除了原本預設的路徑，也請到我指定的這個 `bin` 資料夾裡尋找指令」。這個設定只在當前的終端機視窗有效，關閉後就會消失，非常適合臨時教學使用。
    

#### **步驟 1.5：驗證安裝**

執行以下指令，檢查是否安裝成功。

猛擊

```
blastn -version
```

如果畫面成功顯示版本號 (例如 `blastn: 2.17.0+`)，恭喜你，環境已經準備就緒！

---

## **三、 Part 2: 準備分析序列**

工具備妥後，我們需要準備「要查詢的序列」以及「用來比對的資料庫序列」。

#### **步驟 2.1：下載查詢序列 (Query)**

我們下載大腸桿菌 (_E. coli_) 的 16S rRNA 基因序列，作為我們這次搜尋的目標。

猛擊

```
wget -O e_coli_16s.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NR_024570.1&rettype=fasta&retmode=text"
```

#### **步驟 2.2：下載資料庫用的序列**

為了展示效果，我們一次性從 NCBI 下載 12 種不同微生物的 16S rRNA 序列，用它們來建立我們的本地資料庫。

猛擊

```
wget -O more_sequences.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NR_024570.1,NR_112115.1,NR_026078.1,NR_118997.1,NR_117723.1,NR_113278.1,NR_115629.1,NR_113958.1,NR_074246.1,NR_114421.1,NR_029202.1,NR_075073.1&rettype=fasta&retmode=text"
```

---

## **四、 Part 3: 建立客製化 BLAST 資料庫**

直接對 FASTA 檔案進行搜尋效率很低。我們需要用 `makeblastdb` 將序列檔案格式化成 BLAST 專用的索引檔，這能極大化提升搜尋速度。

#### **步驟 3.1：執行`makeblastdb`**

猛擊

```
makeblastdb -in more_sequences.fasta -dbtype nucl -out demo_db_medium
```

- `-in`: 指定輸入的 FASTA 檔案。
- `-dbtype nucl`: 指定序列類型為核酸 (`nucl`)。
- `-out`: 為你的資料庫取一個名字。

#### **步驟 3.2：檢查資料庫檔案**

執行 `ls`，你會發現多了幾個以 `demo_db_medium` 為開頭的檔案 (`.nhr`, `.nin`, `.nsq` 等)，這代表資料庫已成功建立。

---

## **五、 Part 4: 執行 BLAST 比對與分析**

萬事俱備，現在就來執行最終的比對！

#### **步驟 4.1：執行`blastn`**

我們用大腸桿菌的序列 (`e_coli_16s.fasta`)，去搜尋剛剛建立好的中型資料庫 (`demo_db_medium`)。

猛擊

```
blastn -query e_coli_16s.fasta -db demo_db_medium -out medium_blast_results.txt
```

- `-query`: 指定查詢序列。
- `-db`: 指定要搜尋的本地資料庫名稱。
- `-out`: 將結果儲存到指定的檔案。

#### **步驟 4.2：查看結果**

使用 `cat` 或 `less` 指令來查看純文字的結果報告。

猛擊

```
cat medium_blast_results.txt
```

#### **步驟 4.3：解讀結果**

在結果報告中，你會看到一個排序列表。

- **第一名** 應該是 _E. coli_ 自己，因為我們的資料庫裡也包含了它，所以會是 100% 匹配。
    
- **接下來** 會是與 _E. coli_ 親緣關係很近的物種，如志賀氏菌 (_Shigella_)，它們的比對分數 (Score) 和相似度 (Ident) 都會非常高。
    
- **排名較後** 的則是親緣關係較遠的物種，分數和相似度會依次降低。
    

**至此，你已經完成了一次從無到有的完整本地端 BLAST 分析！**
