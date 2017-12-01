## Laboratoria data mining - (AGH IET) Marcin Piechota, PhD

## Spis treści
1. [Lab 1](#lab1) 10.11.2017
2. [Lab 2](#lab2) 17.11.2017
3. [Lab 3](#lab3) 1.12.2017

## Lab1
Celem zadania jest identyfikacja populacji komórek w korze mózgowej na podstawie profilu ekspresji genów (patrz wykład) oraz wybranie genów markerowych dla każdej populacji. Na zadanie będą poświęcone trzy laboratoria.

1. Feature selection, redukcja wymiarowości
2. Hierachiczna klasteryzacja, identyfikacja markerów
3. Wizualizacja i omówienie wyników

Należy poklasyfikować komórki według i wzorca ekspresji. Zbiór danych dotyczy >1800 komórek nerwowych i ich trankryptomów (> 20000 genów każdy).

Dane znajdują się pod adresem https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71585.
Należy skorzystać z pliku GSE71585_RefSeq_RPKM.csv.gz.
Publikacja znajduje się pod adresem https://www.ncbi.nlm.nih.gov/pubmed/26727548.

### 1
```
otworzyć rstudio i utworzyć nowy projekt
```

### 2
```
pobrać plik
```

### 3
```
podejrzeć plik
```

Można wykorzystać `cat`, `head`, `cut`, `gzip -d`, `gunzip`

### 4
```
wczytać plik
```

R. Funkcja `read.csv`

### 5
```
policzyć średnią z każdego wiersza
```

Wykorzystać `apply` zamiast `for`. W domu przeczytać dlaczego w R korzystamy z `apply`.

### 6
```
policzyć decyle, średnią, medianę, kwartyle z każdego wiersza
```

Można wykorzystać funckję `summary` dla wektorów liczbowych oraz funckję `quantile`

### 7
```
obejrzec rozkłady średnich, median, różnic między kwartylami
```

Wykorzystać funckję `hist`. Rozkłady skośne można oglądać wykorzystując transformację `log2(1 + x)`

### 8
```
wykonać redukcję wymiarów przy pomocy PCA dla pierwszych 300 genów
```

Wykorzystać funkcję `prcomp`

### 9
```
przedstawić na wykresie typu scatterplot rozłożenie próbek na dwóch pierwszych osiach PCA1 i PCA2
```

Funckja `plot`

### 10
```
przefiltrować features na podstawie rozkładu i powtórzyć kroki 8 i 9 aż do uzyskania rozdzielonych "chmurek"
```

Można wykorzystać funcję `which`

## Lab2
Kod do Lab 1
```
url <- "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE71585&format=file&file=GSE71585%5FRefSeq%5FRPKM%2Ecsv%2Egz"
file <- "rpkm.csv.gz"
download.file(url,file)
data <- read.csv(
  gzfile(file),
  sep=",",
  header=TRUE,
  stringsAsFactors=FALSE,
  row.names=1)

colnames(data) <- 1:ncol(data)
which.columns <- 1:ncol(data)

selected.data <- data[,which.columns]
selected.desc <- data.frame(t(apply(selected.data, 1, summary)),
                   t(apply(selected.data, 1, quantile, prob = seq(0, 1, length = 11))))
which.rows <- which(selected.desc$Min < 0.5 & selected.desc$Max < 500 & selected.desc$X20. < 1 & selected.desc$X80. > 30)
pca <- prcomp(t(selected.data[which.rows,]))
plot(pca$x[,"PC1"], pca$x[,"PC2"])
```
### 1
```
Utworzyć macierz odległości euklidesowych pomiędzy "komórkami" w przestrzeni rozpiętej na osiach PCA1 i PCA2
```
Wykorzystać funkcję `dist`

### 2
```
Wykonać hierarchiczną klasteryzację na utworzonej macierzy odległości. Wykorzystać metodę składania 'average'
```
Funckja `hclust`

### 3
```
Wykonać dendrogram 
```
Funckja `as.dendrogram`

### 4
```
Narysować dendrogram
```

### 5
```
Przedstawić macierz odległości za pomocą mapy cieplnej. Do mapy cieplnej wykorzystać utworzony dendrogram.
```
Funckja `heatmap`

### 6
```
Wykonać kroki 2 do 5 (bez 4) na innych metodach składania. Wnioski?
```

### 7
```
Podzielić drzewo na dwie gałęzie odchodzące od korzenia.
```
Funkcja `cutree`

### 8
```
Na wszystkich genach poszukać najlepszego markera rozdzielającego dwie gałęzie.
```
Czy istnieje pojedynczy marker rozdzielający elementy z czułością i specyficznością większą od 98%. Jaka jest nazwa tego genu? Czy istnieje marker mniej licznej gałęzi? Czy istnieje marker bardziej licznej gałęzi?

### 8.1
```
Narysować krzywą ROC dla predyktora
```
pakiet `pROC`. Funkcja `roc`

### 9
```
Usunąć komórki źle sklasyfikowane przez marker - można przypisać im klasę "0". Pozostałym przypisać klasy "1" i "2". Można także nazwać je od markera.
```

### 10
```
Na wydzielonych klasach powtórzyć iteracyjnie procedurę, aż do momentu, w którym nie da się znaleźc markera dla żadnej z gałęzi.
```
Opisać komórki klasami

Kod
```
url <- "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE71585&format=file&file=GSE71585%5FRefSeq%5FRPKM%2Ecsv%2Egz"
file <- "rpkm.csv.gz"
download.file(url,file)
data <- read.csv(
  gzfile(file),
  sep=",",
  header=TRUE,
  stringsAsFactors=FALSE,
  row.names=1)
colnames(data) <- 1:ncol(data)

MIN_CLUSTER_COUNT = 10

mergeVectors <- function(x) {
  apply(x, 1, function(y) {
    wh <- which(y != "-")
    if(length(wh) == 0) {return("-")}
    return(y[wh])
  })
}

score <- function(x, cutree, clusternum = 1) {
  x <- as.numeric(x)
  division <- (table(cutree))[3 - clusternum]
  threshold <- sort(x)[division]
  threshold
  sum(x[cutree == clusternum] > threshold) + sum(x[cutree == (3 - clusternum)] <= threshold)
}


getCellIndexes <- function(x, cutree, clusternum = 1) {
  out <- cutree
  x <- as.numeric(x)
  division <- (table(cutree))[3 - clusternum]
  threshold <- sort(x)[division]
  which(cutree == clusternum & x > threshold)
}

analyse <- function(which.columns, name = "-") {
  selected.data <- data[,which.columns]
  selected.desc <- data.frame(
                      t(apply(selected.data, 1, summary)),
                      t(apply(selected.data, 1, quantile, prob = seq(0, 1, length = 11)))
                      )
  which.rows <- which(selected.desc$Min < 0.5 & selected.desc$Max < 500 & selected.desc$X20. < 1 & selected.desc$X80. > 30)
  print("Computing PCA...")
  pca <- prcomp(t(selected.data[which.rows,]))
  print("Plotting PCA...")
  plot(pca$x[,1:2])
  print("Clustering...")
  dist <- dist(pca$x[,1:2])
  hc <- hclust(dist, method = 'average')
  dend <- as.dendrogram(hc)
  print("Plotting heatmap...")
  heatmap(as.matrix(dist), Rowv = dend, Colv = dend)
  cutree <- (cutree(hc, k = 2))
  if(min(table(cutree)) < MIN_CLUSTER_COUNT) {
    out <- rep("-", ncol(data))
    indexOfBiggerCluster = order(table(cutree))[2]
    out[which.columns[cutree == indexOfBiggerCluster]] <- name
    return(out)
  }
  print("Scoring genes...")
  scores_one <- apply(selected.data, 1, score, cutree = cutree, clusternum = 1)
  print("Scoring genes...")
  scores_two <- apply(selected.data, 1, score, cutree = cutree, clusternum = 2)
  print("Sorting...")
  gene1id = order(scores_one, decreasing = T)[1]
  gene1name = rownames(data)[gene1id]
  gene2id = order(scores_two, decreasing = T)[1]
  gene2name = rownames(data)[gene2id]
  plot(1:ncol(selected.data),selected.data[gene1id,hc$order])
  plot(1:ncol(selected.data),selected.data[gene2id,hc$order])
  out <- rep("-", ncol(data))
  out[which.columns[getCellIndexes(selected.data[gene1id,], cutree, 1)]] <- gene1name
  out[which.columns[getCellIndexes(selected.data[gene2id,], cutree, 2)]] <- gene2name
  out
}

markers <- matrix(data = "-", nrow = 25, ncol = ncol(data))
markers[1,] <- analyse(1:ncol(data))
for(i in 1:25) {
  genes <- names(table(markers[i - 1,]))[table(markers[i - 1,]) > 10]
  genes <- genes[genes != "-"]
  markers[i,] <- mergeVectors(
    sapply(genes, function(x) {analyse(which(markers[i - 1,] == x), name = x)})
  )
}

```

## Lab3
Kod do Lab 3
```
url <- "http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE71585&format=file&file=GSE71585%5FRefSeq%5FRPKM%2Ecsv%2Egz"
file <- "rpkm.csv.gz"
download.file(url,file)
data <- read.csv(
  gzfile(file),
  sep=",",
  header=TRUE,
  stringsAsFactors=FALSE,
  row.names=1)
colnames(data) <- 1:ncol(data)
load(url("https://github.com/marpiech/bioinfo-recipes/raw/master/markers.RData"))
```
### 1
Zapoznać się z macierzą `markers`
### 2
Do dalszej pracy wybieramy 6 markerów
```
which.genes <- c("Snap25", "Slc6a1", "Synpr", "Satb2", "Qk", "Sept4")`
```
### 3
Wybrać odpowiednie numery wierszy w macierzy `data` na podstawie nazw markerów

### 4
Obliczyć macierz odległości na podstawie korelacji. Wykorzystać funkcję `cor`. Miejsca zawierające wartości `NA` uzupełnić jako zerowa korelacja.

Pytanie pomocnicze: czy wysoka korelacja oznacza małą czy dużą odległość pomiędzy wektorami?

### 5

Posortować kolumny w macierzy według hierarchicznej klasteryzacji

Uwaga: jeśli nie zmienimy na tym etapie nazw kolumn zostaną one posortowane alfanumerycznie i podczas wizualizacji ggplot nasze sortowanie zostanie odwrócone.

### 6

Macierz przetransformować na postać normalną do `data.frame`. Funkcja `melt`

### 7

Poziomy kespresji wybranych genów przedstawić na wykresie słupkowym.
Funkcja `ggplot`. Rodzaj wykresu `geom_bar`. Każdy gen powinien znaleźć się na osobnym panelu. Wykorzystać do tego `facet_wrap`

### 8

Procedurę wykonać dla innych zestawów markerów

### 9
Na wykresie pca z Lab1 przedstawić dowolny marker
