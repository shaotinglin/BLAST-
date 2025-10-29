

---

# **本地端 BLAST+ 完整實戰教學**

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




下一步！從序列資料建立演化樹!

要將之前下載的 12 條 16S rRNA 序列轉換成演化樹，我們需要使用到一些不同的工具。BLAST 和 MMseqs2 主要用於序列比對和搜尋，而建立演化樹則需要經過以下兩個主要步驟：

1.  **多序列比對 (Multiple Sequence Alignment, MSA)**：將所有序列進行比對，找出它們之間保守和變異的區域。這一步是構建演化樹的基礎。
2.  **演化樹建構 (Phylogenetic Tree Construction)**：根據多序列比對的結果，利用演化模型推斷序列之間的親緣關係，並繪製成樹狀圖。

-----

# **演化樹建構教學：從多序列比對到樹狀圖**

#### **事前準備：你需要額外安裝的工具**

請注意，以下兩個工具你需要先安裝。它們通常在大多數 Linux 發行版的軟體倉庫中都可以找到，或者你可以手動下載安裝（類似於 BLAST+ 的手動安裝方式）。

1.  **多序列比對工具：MAFFT** (或 Clustal Omega)

      * **MAFFT** 以其速度和準確性而聞名。
      * **安裝指令 (Ubuntu/Debian 系統)：**
        ```bash
        sudo apt update
        sudo apt install mafft
        ```
      * **安裝指令 (CentOS/Fedora 系統)：**
        ```bash
        sudo yum install mafft # 或 sudo dnf install mafft
        ```
      * **驗證安裝：** `mafft --version`

2.  **演化樹建構工具：FastTree** (或 IQ-TREE, RAxML)

      * **FastTree** 速度極快，適合較大的數據集。
      * **安裝指令 (Ubuntu/Debian 系統)：**
        ```bash
        sudo apt install fasttree
        ```
      * **安裝指令 (CentOS/Fedora 系統)：**
        ```bash
        sudo yum install fasttree # 或 sudo dnf install fasttree
        ```
      * **驗證安裝：** `fasttree -h`

-----

### **步驟一：執行多序列比對 (使用 MAFFT)**

我們將使用之前下載的 `more_sequences.fasta` 檔案，它包含了 12 條 16S rRNA 序列。

```bash
# 執行 MAFFT 多序列比對
mafft --auto more_sequences.fasta > aligned_sequences.fasta
```

  * `mafft`: 呼叫 MAFFT 程式。
  * `--auto`: 讓 MAFFT 自動選擇最適合比對策略。
  * `more_sequences.fasta`: 輸入檔案 (包含所有 12 條序列)。
  * `>`: 將 MAFFT 的輸出導向到一個新的檔案。
  * `aligned_sequences.fasta`: 輸出檔案，包含了所有已經比對好的序列。

你可以使用 `cat aligned_sequences.fasta` 查看比對後的序列，你會看到序列中多了許多 `-` (gap)，表示插入或刪除。

-----

### **步驟二：建構演化樹 (使用 FastTree)**

接下來，我們將使用 MAFFT 輸出比對好的序列檔案 (`aligned_sequences.fasta`)，透過 FastTree 來建構演化樹。

```bash
# 執行 FastTree 建構演化樹
FastTree -nt -gtr -gamma -out phylogenetic_tree.nwk aligned_sequences.fasta
```

  * `FastTree`: 呼叫 FastTree 程式。
  * `-nt`: 指定輸入序列是核酸 (nucleotide) 類型。
  * `-gtr`: 使用廣義時間可逆 (General Time Reversible) 模型，這是核酸序列演化常用的模型。
  * `-gamma`: 考慮不同位點演化速率的變異 (Gamma 分佈)。
  * `-out phylogenetic_tree.nwk`: 將建構好的演化樹輸出到名為 `phylogenetic_tree.nwk` 的檔案。
      * `.nwk` 副檔名代表 Newick 格式，這是演化樹最常見的文字格式。

-----

### **步驟三：查看演化樹檔案**

現在，你的資料夾中應該有一個 `phylogenetic_tree.nwk` 的檔案。它是一個純文字檔，但內容是 Newick 格式，看起來像一串用括號和逗號組成的字串。

```bash
cat phylogenetic_tree.nwk
```

**範例 Newick 格式 (這不是實際結果，僅為示意)：**

```
((A:0.1,B:0.2):0.05,(C:0.3,D:0.4):0.01);
```

這串文字就代表了序列之間的關係和分支長度。

-----

### **步驟四：演化樹視覺化 (重要！)**

雖然我們已經生成了演化樹的數據 (`.nwk` 檔案)，但要把它變成直觀的樹狀圖，通常需要一個**圖形介面 (GUI)** 的軟體。這一步無法在命令列中直接完成，你需要將 `phylogenetic_tree.nwk` 檔案下載到你的個人電腦上，然後使用以下軟體開啟：

1.  **FigTree** (推薦，免費且常用)：

      * 下載網址：[http://tree.bio.ed.ac.uk/software/figtree/](http://tree.bio.ed.ac.uk/software/figtree/)
      * 這是一個跨平台軟體，支援 Windows, macOS, Linux。你可以載入 `.nwk` 檔案，然後調整樹的樣式、標籤、顏色等。

2.  **iTOL (Interactive Tree Of Life)** (線上工具)：

      * 網址：[https://itol.embl.de/](https://itol.embl.de/)
      * 這是一個強大的線上工具，你可以上傳你的 `.nwk` 檔案，它會幫你繪製並提供豐富的註釋功能。

-----

### **總結流程與指令懶人包**

以下是從頭到尾的所有指令，假設你已經安裝了 MAFFT 和 FastTree。

```bash
# --- Part 0: 環境設定與工具安裝 (請根據你的系統安裝 MAFFT 和 FastTree) ---
# For Ubuntu/Debian:
# sudo apt update
# sudo apt install mafft fasttree

# For CentOS/Fedora:
# sudo yum install mafft fasttree # 或 sudo dnf install mafft fasttree

# (假設你已經在 /my_blast_project/ 資料夾下，並已下載 more_sequences.fasta)
# --- Part 1: 多序列比對 (Multiple Sequence Alignment) ---
echo "--- 正在執行 MAFFT 多序列比對 ---"
mafft --auto more_sequences.fasta > aligned_sequences.fasta
echo "比對完成，結果儲存至 aligned_sequences.fasta"

# --- Part 2: 建構演化樹 (Phylogenetic Tree Construction) ---
echo "--- 正在執行 FastTree 建構演化樹 ---"
FastTree -nt -gtr -gamma -out phylogenetic_tree.nwk aligned_sequences.fasta
echo "演化樹建構完成，結果儲存至 phylogenetic_tree.nwk (Newick 格式)"

# --- Part 3: 查看演化樹檔案 (文字格式) ---
echo "--- 演化樹檔案內容 (Newick 格式) ---"
cat phylogenetic_tree.nwk

# --- Part 4: 視覺化提醒 ---
echo "--- 重要：請將 phylogenetic_tree.nwk 下載到你的電腦 ---"
echo "然後使用 FigTree (軟體) 或 iTOL (線上工具) 進行視覺化，才能看到樹狀圖。"
```

透過這些步驟，你就可以將那 12 條序列成功轉換為一個視覺化的演化樹，展示它們之間的親緣關係了！

