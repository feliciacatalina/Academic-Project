library(readxl)
data = read_excel("D:/data pengangguran.xlsx", sheet = "Sheet1")
tahun2019 = data$`Tingkat Pengangguran Agustus 2019`
tahun2020 = data$`Tingkat Pengangguran Agustus 2020`

thn2019 = as.numeric(tahun2019)
thn2020 = as.numeric(tahun2020)
da = data.frame(thn2019,thn2020)

# sari numerik
summary(da)

# kruskal wallis test
kruskal.test(list(thn2019, thn2020))

# wilcoxon rank sum test
wilcox.test(thn2019,thn2020,paired = FALSE, data = data)

# boxplot
boxplot(da$thn2019, ylab = 'Tingkat Pengangguran (%)', main = 'Tingkat Pengangguran di Indonesia pada Agustus 2019')
boxplot(da$thn2020, ylab = 'Tingkat Pengangguran (%)', main = 'Tingkat Pengangguran di Indonesia pada Agustus 2020')

barplot(da[, 'thn2019'], main = 'Tingkat Pengangguran di Indonesia pada Agustus 2019 (%)')
barplot(da[, 'thn2020'], main = 'Tingkat Pengangguran di Indonesia pada Agustus 2020 (%)')

dens1 = density(da$thn2019)
hist(da$thn2019, freq = FALSE, col = 'blue', main = 'Tingkat Pengangguran di Indonesia pada Agustus 2019 (%)')
polygon(dens1, border='red')

dens2 = density(da$thn2020)
hist(da$thn2020, freq = FALSE, col = 'yellow', main = 'Tingkat Pengangguran di Indonesia pada Agustus 2020 (%)')
polygon(dens2, border='darkblue')

cor.test(x=da$thn2019, y=da$thn2020, method = 'spearman')
