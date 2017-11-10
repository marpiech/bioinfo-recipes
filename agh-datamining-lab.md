### Lab 1
Celem zadania jest identyfikacja populacji komórek w korze mózgowej na podstawie profilu ekspresji genów (patrz wykład) oraz wybranie genów markerowych dla każdej populacji. Na zadanie będą poświęcone trzy laboratoria.

1. Feature selection, redukcja wymiarowości
2. Hierachiczna klasteryzacja, identyfikacja markerów
3. Wizualizacja i omówienie wyników

Należy poklasyfikować komórki według i wzorca ekspresji. Zbiór danych dotyczy >1800 komórek nerwowych i ich trankryptomów (> 20000 genów każdy).

Dane znajdują się pod adresem https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71585.
Należy skorzystać z pliku GSE71585_RefSeq_RPKM.csv.gz.
Publikacja znajduje się pod adresem https://www.ncbi.nlm.nih.gov/pubmed/26727548.

## 1
`otworzyć rstudio`

# 2
`pobrać plik`

# 3
`podejrzeć plik`
Można wykorzystać cat, head, cut, gzip -d, gunzip

# 4
`wczytać plik`
R funkcja read.csv

# 5
`policzyć średnią z każdego wiersza`
Wykorzystać apply zamiast for. W domu przeczytać dlaczego w R korzystamy z apply.

# 6
`policzyć decyle, średnią, medianę, kwartyle z każdego wiersza`
Można wykorzystać funckje summary dla wektorów liczbowych oraz quantile

# 7
`obejrzec rozkłady średnich, median, różnic między kwartylami`
Wykorzystać funckję hist. Rozkłady skośne można oglądać wykorzystując transformację log(1 + x)

# 8
`wykonać redukcję wymiarów przy pomocy PCA dla pierwszych 300 genów`
wykorzystać funkcję prcomp

# 9
`przedstawić na wykresie typu scatterplot rozłożenie próbek na dwóch pierwszych osiach PCA1 i PCA2`
funckja plot

# 10
`przefiltrować features na podstawie rozkładu i powtórzyć kroki 8 i 9 aż do uzyskania rozdzielonych "chmurek"`
można wykorzystać funcję which
