signatures <- list(
  microglia = list(
    `Cain A. et al (2022)` = list(
      `Homeostatic` = c("FRMD4A", "P2RY12", "SYNDIG1", "CX3CR1"),
      `Stress response/anti-inflammatory` = c("TMEM163", "SPP1", "DDIT4", "SLC2A3"),
      `Enhanced redox` = c("RPL19", "RPS6", "FTL", "FTH1", "APOE","C1QA","C1QB","C1QC", "RPS19"),
      `Interferon Response` = c("IFI44L", "IFI6","IFIT3","PTPRG"),
      `Proliferating` = c("ARHGAP11B" ,"TOP2A","HMMR")
    ),
    `Keren-Shaul et al (2017)` = list(
      `DAM1 down` = c("P2RY12","P2RY13","SERINC3","CX3CR1","TGFBR1","TMEM119","TXNIP","GLUL"),
      `DAM1` = c("APOE","B2M","CSTB","TYROBP","TIMP2","HLA-A","HLA-B","HLA-C","FTH1","LYZ","CTSB","CTSD"),
      `DAM2` = c("ANK1","SPP1","AXL","CSF1","CST7","CD9","CADM1","CLEC7A","CCL6","ITGAX","CD63","CD68","CTSA","LPL","GUSB","SERPINE2","CTSZ","CD52","CTSL","HIF1A")
    ),
    `Gerrits E. et al (2021)` = list(
      Homeostatic = c("SH3RF3","TMEM163", "P2RY12","CX3CR1","SOX5"),
      AD1.c7 = c("DSCAM","IPCEF1","SOCS6","ADAMTS17"),
      AD1.c9 = c("ITGAX","MYO1E","GLDN","DTNA","SPATS2L","TPRG1","PTPRG","LPL"),
      AD1.c10 = c("P2RY12","CX3CR1", "STARD13","CPM","EYA2","GPNMB","APOE"),
      AD2 = c("GRID2","ADGRB3","DPP10","DDX17","NAV2","FOXP2"),
      `Inflammation` = c("GPNMB", "CD83", "IL1B", "NFKB1","CCL3","CCL4","SLC2A3","SRGN","GPR183"),
      `Stress` = c("FOS","JUNB", "HSPA1A", "HSPA1B","HSP90AA1","SLC2A3","SRGN","GPR183","HIST1H2AC"),
      `Proliferating` = c("TOP2A", "BRIP1", "MKI67", "FANCI")
    )
  ),
  astrocytes = list(
    `Cain A. et al (2022)` = list(
      `Homeostatic protoplasmic-like` = c("SLC1A2", "DLC1", "WIF1"),
      `Non-homeostatic` =c("GFAP","SERPINA3","OSMR","TPST1","IGFBP7","CHI3L1"),
      `Fibrous-like` = c("GFAP","DPP10","ID3","SLC38A1","FOS", "CD44", "AQP4"),
      `Ast.4` = c("S100A6","MT1G","MT1F","SNCG","MT1A","COL5A3","PDE4DIP"),
      `Interferon-Response` = c("IFI44L", "IFI6","IFIT1","B2M")
    ),
    # Astrocytes and oligodendrocytes undergo subtypespecific transcriptional changes in Alzheimer’s disease
    `Sandick J.S. et al (2022)` = list(
      `Ast.0` = c("TRPM3","ARHGAP24","ZNF98","CABLES1","HPSE2","EGFR", "LRRC4C","EPHB1"),
      `Ast.1` = c("CST3","APOE","CLU","FTH1","HEPN1","PSAP", "MT-CO1","MT-ND1","MT-ND3","ITM2B","ITM2C"),
      `Ast.2` = c("DPP10","LINC00609","ADAMTSL3","L3MBTL4","AC012405.1"),
      `Ast.3` = c("RHGEF3","SAMD4A","TPST1","SERPINA3","CHI3L1","C3","OSMR"),
      `Ast.4` = c("FOS","LINC01411","DPP6","SLC38A1","DCLK1","NTNG1","GRIA1","GRIK4","SHISA6"),
      `Ast.5` = c("KAZN","LINC01088","ADAMTSL3","FBN1","SORBS1","SPIRE1"),
      `Ast.6` = c("DDIT4","ID3","JUN","HSPA1A"),
      `Ast.7` = c("AL353138.1","ZMAT3","MIR34AHG","KCNMB2−AS1","PTCHD4"),
      `Ast.8` = c("CREB5","SGCD","LTBP1","DACH1","MOXD1","EPHA4","AKAP12","NLGN4X")
    ),
    `Sandick J.S et al (2022) - Integrated` = list(
      Int.0 = c("NF98","ARHGAP24","HPSE2","CABLES1","SGCD"),
      Int.1 = c("HEPN1","FTH1","APOE","ITM2C","CST3"),
      Int.2 = c("DPP10","L3MBTL4","ADAMTSL3","CACNA2D3","AC012405.1"),
      Int.3 = c("FOS","DPP6","LINC01411","SLC38A1","DCLK1"),
      Int.4 = c("SAMD4A","ARHGEF3","TPST1","SERPINA3","CHI3L1"),
      Int.5 = c("KAZN","LINC00609","LINC01088"),
      Int.6 = c("HS3ST3A1","DPP10−AS3","ID3"),
      Int.7 = c("DST","PITPNC1","AC105052.4","MACF1","SYNE1"),
      Int.8 = c("EDIL3","ELMO1","IL1RAPL1","ANK3","SLC24A2")
    ),
    `Habib N. et al (2020)` = list(
      `Homeostatic + GFAP-low` = c("LUZP2","MGAT4C","TRPM3","HSALNG0103764","GPC5","GRIA2","LSAMP","SLC7A10","ST6GALNAC5","UPF2"), # Ast 1-2
      `low-DAA` = c("APOE","CLU","GPC5","NNAT","TRPM3","GRM5","NCAM2"), # Ast 3
      DAA = c("APOE","GFAP","CTSB","VIM","OSMR","SERPINA3","GSN","GGTA","CTSD","CTSL","CHI3L1"), # Ast 4
      `low-high` = c("GRIN3A","MYOC","AQP4","ID3","SNHG11","GPM6B"), # Ast 5
      `GFAP-high` = c("GFAP","ID3","AQP4","MYOC","ID1","FABP3","SLC38A1") # Ast 6
    )
  ),
  oligodendroglia = list(
    `Cain A. et al (2022)` = list(
      Oli.1 = c("MOG","SVEP1","FCHSD2","GRIN2A","CFTR","NAALADL2","DMD","SLC7A14","PLPPR1","DOCK2","NRXN3","PPP2R2B","LTBP1"),
      Oli.2 = c("QDPR","CPQ","S100A6","SLC38A2","PIM3","IGF1R","TPD52L1","CHSY3","CXCR4","FBXO2","ETV5","PLEKHH2","CTNNA2"),
      Oli.3 = c("MOG","KCTD8","RBFOX1","PLEKHG1","ACTN2","SLC5A11","RASGRF1","MACROD2","NAV3","DTNA","NEGR1","RASGRF2","PDE4D"),
      Oli.4 = c("QDPR","CLU","UBL5","COX6A1","MT2A","HINT1","UBB","TMSB4X","RPS23","RPL32","SOD1","CALM1","CKB")
    ),
    
    # Astrocytes and oligodendrocytes undergo subtypespecific transcriptional changes in Alzheimer’s disease
    `Sadick J.S. et al (2022)` = list(
      Oli.0 = c("HIP1","ANKRD18A", "SLC5A11","LURAP1L−AS1"),
      Oli.1 = c("SELENOP","DBNDD2","CRYAB","PLP1","FTH1","CD9","CNP","SEPTIN4","SERINC3"),
      Oli.2 = c("COL18A1","AFF3","ACSBG1","RASGRF1","RBFOX1"," MSMO1","FDFT1","LSS"),
      Oli.3 = c("DPP10","SLC1A2","GPM6A","ADGRV1","NRG3"),
      Oli.4 = c("NAV2","CAMK2D","BIRC3","IFIT2","ZC3HAV1","PSMB1","B2M","HLA-A")
    ),
    # # Integrated with Mathys, Grubman and Zhou
    `Sadick J.S. et al (2022) - Integrated` = list(
      Int.0 = c("DYSF","PLXDC2","LINC01608","SVEP1"),
      Int.1 = c("CRYAB","QDPR","ST3GAL6","CNDP1","CTNNA2","FTH1", "PLP1", "APOD", "DBNDD2"),
      Int.2 = c("MT−ATP6","MT−CO2","MT−ND3","MT−ND4","FP236383.3"),
      Int.3 = c("ANKRD18A","LINC00609","RASGRF1","SLC5A11","ACTN2"),
      Int.4 = c("FRY","KCNIP4","CNTN1","MDGA2","SGCZ","FTH1", "PLP1", "APOD", "DBNDD2"),
      Int.5 = c("COL18A1","ACSBG1","AFF3","RBFOX1"),
      Int.6 = c("NEAT1","CAMK2D","NAV2","LUCAT1","NRP2")
    ),
    # Marques S. et al., (2016) Oligodendrocyte heterogeneity in the mouse juvenile and adult central nervous system
    `Marques` = list(
      OPC = c("PTPRZ1","PDGFRA","SERPINE2","CSPG5","VCAN","CSPG4"),
      COP = c("CD9", "NEU1", "BMP4","GPR17","VCAN"),
      NFOL = c("ARPC1B","TMEM2","CCHN2","MPZL1","FRMD4A","MOBP","DDR1","TSPAN2"),
      MFOL = c("CTPS","TMEM141","OPALIN","MAL","PRGDS","EVI2A","EVI2B"),
      MOL = c("APOD","SEPP1","S100B","FOSB","DUP1","DNAJB1","ANXA5","KLK6","MGST3","CAR2","CNTN2","GAD2","SERPINB1A","NEAT1","CYP51","DHCR24","PDLIM2","IL33","APOE","PTGDS")
    ),
    # Pandey S. et al, (2022) Disease-associated oligodendrocyte responses across neurodegenerative diseases
    `Pandey S. et al (2022)` = list(
      MOL56 = c("PTGDS", "OPALIN", "QDPR", "TAGLN2", "GADD45B", "VIM", "SYT4", "SCN1B", "ENO2", "BTBD16", "GADD45G", "CNTN1"),
      MOL56_DA1 = c("PTGDS", "OPALIN", "QDPR", "TAGLN2", "GADD45B", "VIM", "SYT4", "SCN1B", "ENO2", "BTBD16", "GADD45G", "CNTN1", "C4B","SERPINA3", "ANXA2", "PLVAP", "THBS3", "STEAP3", "EMP3", "PARVB", "S100A10", "TNFRSF1A", "COL6A1", "SEMA4F"),
      MOL56_DA2 = c("PTGDS", "OPALIN", "QDPR", "TAGLN2", "GADD45B", "VIM", "SYT4", "SCN1B", "ENO2", "BTBD16", "GADD45G", "CNTN1", "CDKN1A", "BAX", "DDIT3", "FOS", "ATF4", "EGR1", "CCND1", "TNFRSF12A", "BTG1", "ERG2", "KLF4", "FGF7", "RRAD", "GDF15"),
      MOL56_IFN = c("PTGDS", "OPALIN", "QDPR", "STAT1", "BST2", "IRGM", "PSMB8", "IFIT1", "IRF7", "PSME1", "OASL", "HLA-A", "TAP1", "IFIT2")
    )
  ),
  endo = list(
    `Cain A. et al. (2022)` = list(
      `End.1` = c("TMEM132C", "SLC5A6","GPCPD1"),
      `End.2` = c("RPL34","FTH1","MT-ND4","RGS5"),
      `End.3` = c("SAMD4A","KLF6", "TNFRSF12A","KLF4"),
      `End.4` = c("TSHZ2","CACNA1C","ACKR1","SNTG2"),
      `End.5` = c("ARL15","MGP","COL8A1","CXCL2"),
      `End.6` = c("MT1X","MT1M","DEPP1")
    ),
    # A human brain vascular atlas reveals diverse cell mediators of Alzheimer’s disease risk
    `Yang A.C. et al (2022)` = list(
      `Arterial`  = c("VEGFC", "ALPL"),
      `Capillary` = c("MFSD2A", "SLC7A5"),
      `Venous`    = c("IL1R1" ,"NR2F2"),
      `Tip Cells` = c("PLAUR", "LAMB1"),
      `T-Pericytes` = c("ABCC9", "PTN","SLC6A1", "SLC1A3", "SLC12A7","SLC6A12","SLC6A13","SLC20A2", "PTPRK"),
      `M-Pericytes` = c("ABCC9", "PTN","COL4A1","COL4A2","COL4A3","COL4A4", "LAMC3","ADAMTS1","ADAMTS9","COL18A1", "CRISPLD2"),
      `aSMC` = c("ACTA2", "TAGLN"),
      `vSMC` = c("CTNNA3", "SLIT3"),
      `Perivascular Fibroblast-like` = c("LAMA2", "FBLN1","CYP1B1","ABCA6","ABCA8","ABCA9","ABCA10"),
      `Meningeal Fibroblast` = c("KCNMA1","SLC4A4","SLC24A3","SLC7A2","SLC26A7","SLC41A2","SLC35G1","SLC22A23"),
      `Dural Meningeal` = c("LMO4","SLC47A1","SLC26A2"),
      `Archnoid Meningeal` = c("TRPM3", "SLC7A11","TJP1")
    ),
    `Garcia F. J. et al., (2021)` = list(
      `Arteriole` = c("ABCB1","ATP10A","VEGFC", "ARL15"),
      `Capillary` = c("ABCB1","ATP10A","MFSD2A", "SLC7A5"),
      `Venule` = c("ABCB1","ATP10A","TSHZ2", "ADGRG6"),
      `aSMC` = c("TAGLN","ACTA2", "MYH11"),
      `vSMC` = c("TAGLN","MRC1", "CD74"),
      `Pericytes` =c("TAGLN","GRM8", "PDGFRB"),
      `Fib Type I` = c("CEMIP", "ABCA10", "FBLN1"),
      `Fib Type II` = c("CEMIP", "TRPM3", "MYRIP"),
      `Fib Type III` = c("CEMIP", "KCNMA1", "SLC4A4")
    )
  ),
  excitatory = list(
    # Molecular characterization of selectively vulnerable neurons in Alzheimer’s disease
    `Leng K. et al., (2021)` = list(
      `Layers II/III` = c("GLRA3","LAMP5", "CARTPT"),
      `Layers II/III/IV` = c("CUX2"),
      `Layers III/IV` = c("PRSS12"),
      `Layers V/VI` = c("TOX", "ETV1","RPRM","RXFP1"),
      `Layers VI` = c("TLE4", "FOXP2","NTNG2", "OPRK1","NR4A2", "ADRA2A")
    ),
    
    `Cain A. et al., (2022)` = list(
      `Exc.1` = c("CBLN2", "LAMP5", "CUX2", "GLRA3", "CARTPT"),
      `Exc.2` = c("RORB", "FOXP2", "TSHZ2", "DCC"),
      `Exc.3` = c("PLCH1", "RMST", "OTOGL"),
      `Exc.4` = c("CASC15", "TMSB10"),
      `Exc.5` = c("ADGRL4", "TOX"),
      `Exc.6` = c("ETV1", "HTR2C", "ITGA8", "ZNF385D"),
      `Exc.7` = c("THEMIS", "PDZRN4", "CTNNA2"),
      `Exc.8` = c("ROBO2", "RXFP1", "DLC1", "TLE4", "NFIA"),
      `Exc.9` = c("SEMA3E", "TRPM3", "RYR3", "GRIK4"),
      `Exc.10`= c("SYNPR", "POSTN", "RGS12", "NTNG2", "OPRK1", "NR4A2")
    ),
    `Berg J. et al., (2021)` = list(
      `COL22A1` = c("COL22A1","CNR1","ITGA11","MME","SEMA5A","ATP1B2","GNG4","CD163L1","KCNK2","SYT17","FOXO1","TWIST2","CLMP","CDH19","PLSCR4","CRYAB","F11-AS1","ARHGAP4","ARHGAP15","NGB","KLHDC8A","S100A6","CD24","TEX26","PCDH11Y","LY6H","HPCA","NDUFB7","POLR2L"),
      `CARM1P1` = c("PHLDB2","SPHKAP","PLCXD3","CCDC68","ADCY8","CLMN","CNTNAP3","SEMA6A","NBPF19","DCBLD1","CNTNAP3B","COBLL1","NTN4","PTPRZ1","UST","MDFIC","GOLIM4","CALB1","SCRG1","GNG12-AS1","THSD7A","PRSS12","SMAD3","SWAP70","AMPD3","ADGRF5","CARTPT","CNGB1","LUZP2"),
      `Deep FREM3` = c("SEMA5B","LAMP5","HTR4","PVT1","KLF12","MCHR2","ART3","TMTC1","GYPB","NETO2","SLC24A4","SH3RF2","ANKRD30B","MEG8","CPLX2","EPHA6","CBLN4","TLE1","STXBP6","GNAL","BICC1","QKI","KCTD8","RPS6KA2","ST3GAL1","ANO4","LRATD1","SYNJ2","FIGN","CEMIP2","AQP4-AS1")
    )
  ),
  inhibitory = list(
    `Cain A. et al., (2022)` = list(
      `Inh.1` = c("GAD1", "ADARB2", "DPP10", "PVALB"),
      `Inh.2` = c("VIP", "TAC3", "CALB2"),
      `Inh.3` = c("SST", "TRHDE", "RALYL"),
      `Inh.4` = c("RELN", "CXCL14", "CNR1", "PAX6"),
      `Inh.5` = c("FBXL7", "NRG1", "KIT", "FGF13"),
      `Inh.6` = c("TRPS1", "CA8", "PTPRK"),
      `Inh.7` = c("TOX", "NTNG1", "CACNA2D1")
    )
  )
)
