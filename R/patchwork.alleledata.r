patchwork.alleledata <- function (Pileup, vcf, hg_build)
{
    packagepath = system.file(package = "patchwork")
    system(paste("cp -r ", packagepath, "/python .python", sep = ""))
    system(paste("python .python/mpile2alleles.py ", vcf, " >",
        getwd(), "/pile.alleles", sep = ""))
    system("rm -r .python")
    alf = read.csv("pile.alleles", sep = "\t", header = F)[,
        1:6]
    colnames(alf) = c("achr", "apos", "atype", "aqual", "atot",
        "amut")
    alf$aref = alf$atot - alf$amut
    alf$amax = apply(alf[, 6:7], 1, max)
    alf$amin = alf$atot - alf$amax
    if (length(grep("M", alf$achr)) != 0) {
        alf = alf[-grep("M", alf$achr), ]
    }
    alf$achr = as.character(alf$achr)
    x_x = strsplit(alf$achr[1], "chr")
    if (length(x_x[[1]]) == 1) {
        for (i in 1:22) {
            alf$achr[alf$achr == as.character(i)] = paste("chr",
                as.character(i), sep = "")
        }
        i = "X"
        alf$achr[alf$achr == as.character(i)] = paste("chr",
            as.character(i), sep = "")
        i = "Y"
        alf$achr[alf$achr == as.character(i)] = paste("chr",
            as.character(i), sep = "")
    }
    if (hg_build == "hg18") {
        data(commonSnpsHG18, package = "patchworkData")
    }
    else if (hg_build == "hg19") {
        data(commonSnps132, package = "patchworkData")
    }
    dbSnp = dbSnp[, c(1, 3)]
    dbSnp$dbSnp = T
    alf = merge(alf, dbSnp, all.x = T, all.y = F, by = 1:2)
    alf$dbSnp[is.na(alf$dbSnp)] = F
    return(alf)
}

