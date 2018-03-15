patchwork.plot <- function (Tumor.bam, Tumor.vcf = NULL, Normal.bam = NULL, Normal.vcf = NULL,
    Reference = NULL, Alpha = 1e-04, SD = 1, hg.build = "hg19")
{
    name_s = strsplit(basename(Tumor.bam), "\\.")
    name = paste(name_s[[1]][1:length(name_s[[1]]) - 1], collapse = ".")
    alf <- normalalf <- NULL
    try(load(paste(name, "_copynumbers.Rdata", sep = "")), silent = TRUE)
    if (is.null(alf)) {
        try(load(paste(name, "_pile.alleles.Rdata", sep = "")),
            silent = TRUE)
        if (is.null(alf)) {
            normalalf <- patchwork.alleledata(vcf = Normal.vcf,
                hg_build = hg.build)
            alf = patchwork.alleledata(vcf = Tumor.vcf, hg_build = hg.build)
            cat("Allele Data Generation Complete \n")
            save(alf, normalalf, file = paste(name, "_pile.alleles.Rdata",
                sep = ""))
        }
        if (!is.null(normalalf)) {
            normaltemp = data.frame(achr = normalalf$achr, apos = normalalf$apos,
                normal = T)
            somatic <- merge(normaltemp, alf, by = 1:2, all.x = F,
                all.y = T)
            somatic <- somatic[is.na(somatic$normal), ]
            somatic <- somatic[!somatic$amut == 0, ]
            save(somatic, file = paste(name, "_somatic.Rdata",
                sep = ""))
            normalalf <- normalalf[normalalf$amin/normalalf$atot >
                0.2, ]
            alf <- merge(normalalf[, 1:2], alf, by = 1:2, all = F)
        }
        else {
            alf <- alf[alf$dbSnp == T, ]
            alf <- alf[alf$amin > 1 & alf$aref > 2, ]
        }
        data(ideogram, package = "patchworkData")
        chroms = as.character(unique(ideogram$chr))
        data <- normaldata <- NULL
        try(load(paste(name, "_data.Rdata", sep = "")), silent = TRUE)
        if (is.null(data)) {
            cat("Initiating Read Chromosomal Coverage \n")
            data = patchwork.readChroms(Tumor.bam, chroms)
            cat("Read Chromosomal Coverage Complete \n")
            save(data, file = paste(name, "_data.Rdata", sep = ""))
        }
        if (length(data) < 7) {
            cat("Initiating GC Content Normalization \n")
            data = patchwork.GCNorm(data)
            cat("GC Content Normalization Complete \n")
            save(data, file = paste(name, "_data.Rdata", sep = ""))
        }
        if (!is.null(Normal.bam) & is.null(Reference)) {
            normaldata = NULL
            try(load(paste(name, "_normaldata.Rdata", sep = "")),
                silent = TRUE)
            if (is.null(normaldata)) {
                cat("Initiating Read Chromosomal Coverage (matched normal) \n")
                normaldata = patchwork.readChroms(Normal.bam,
                  chroms)
                cat("Read Chromosomal Coverage Complete (matched normal) \n")
            }
            if (length(normaldata) != 3) {
                cat("Initiating GC Content Normalization (matched normal) \n")
                normaldata = patchwork.GCNorm(normaldata[, 1:6])
                cat("GC Content Normalization Complete (matched normal) \n")
                normaldata <- data.frame(chr = normaldata$chr,
                  pos = normaldata$pos, normal = normaldata$norm)
                save(normaldata, file = paste(name, "_normaldata.Rdata",
                  sep = ""))
            }
        }
        kbsegs = NULL
        try(load(paste(name, "_smoothed.Rdata", sep = "")), silent = TRUE)
        if (length(kbsegs) == 0) {
            cat("Initiating Smoothing \n")
            kbsegs = patchwork.smoothing(data, normaldata, Reference,
                chroms)
            save(kbsegs, file = paste(name, "_smoothed.Rdata",
                sep = ""))
            cat("Smoothing Complete \n")
        }
        segs = NULL
        try(load(paste(name, "_Segments.Rdata", sep = "")), silent = TRUE)
        if (length(segs) == 0) {
            cat("Initiating Segmentation \n")
            cat("Note: If segmentation fails to initiate the probable reason is that you have not ")
            cat("installed the R package DNAcopy. See the homepage, http://patchwork.r-forge.r-project.org/ ,\n\t\t\tor ?patchwork.readme for installation instructions. \n")
            segs = patchwork.segment(kbsegs, chroms, Alpha, SD)
            save(segs, file = paste(name, "_Segments.Rdata",
                sep = ""))
            cat("Segmentation Complete \n")
        }
        if (length(segs) == 6) {
            cat("Initiating Segment data extraction (Medians and AI) \n")
            segs = patchwork.Medians_n_AI(segs, kbsegs, alf)
            save(segs, file = paste(name, "_Segments.Rdata",
                sep = ""))
            cat("Segment data extraction Complete \n \n \n")
        }
        cat(paste("Saving information objects needed for patchwork.copynumbers() in ",
            name, "_copynumbers.Rdata \n \n \n", sep = ""))
        save(segs, alf, kbsegs, file = paste(name, "_copynumbers.Rdata",
            sep = ""))
    }
    cat("Initiating Plotting \n")
    karyotype(as.character(segs$chr), segs$start, segs$end, segs$median,
        segs$ai, as.character(kbsegs$chr), kbsegs$pos, kbsegs$ratio,
        alf$achr, alf$apos, (1 - alf$amin/alf$amax), name = as.character(name))
    karyotype_chroms(as.character(segs$chr), segs$start, segs$end,
        segs$median, segs$ai, as.character(kbsegs$chr), kbsegs$pos,
        kbsegs$ratio, alf$achr, alf$apos, (1 - alf$amin/alf$amax),
        name = as.character(name))
    cat("Plotting Complete \n")
    cat("patchwork.plot Complete.\n")
    cat("Below you may see some warning messages, you can read about these on our homepage. They are either nothing to be worried about (\"Tried to load file, it didn't exist.\") or something you should send us an email about. \n")
}

