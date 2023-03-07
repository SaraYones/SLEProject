#Prepare Metadata

phenotype=metadata13
phenotype$c3=gsub("c3: (.*)","\\1",phenotype$c3)
phenotype$neutrophil_count=gsub("neutrophil_count: (.*)","\\1",phenotype$neutrophil_count)
phenotype$neutrophil_percent=gsub("neutrophil_percent: (.*)","\\1",phenotype$neutrophil_percent)
phenotype$lymphocyte_count=gsub("lymphocyte_count: (.*)","\\1",phenotype$lymphocyte_count)
phenotype$lymphocyte_percent=gsub("lymphocyte_percent: (.*)","\\1",phenotype$lymphocyte_percent)
phenotype$c4=gsub("c4: (.*)","\\1",phenotype$c4)
phenotype$wbc=gsub("wbc: (.*)","\\1",phenotype$wbc)
phenotype$MC=gsub("monocyte_count: (.*)","\\1",phenotype$monocyte_count)
phenotype$MP=gsub("monocyte_percent: (.*)","\\1",phenotype$monocyte_percent)
phenotype$PC=gsub("platelet_count: (.*)","\\1",phenotype$platelet_count)
phenotype$esr=gsub("esr: (.*)","\\1",phenotype$esr)
phenotype$hgb=gsub("hgb: (.*)","\\1",phenotype$hgb)
phenotype$hct=gsub("hct: (.*)","\\1",phenotype$hct)
phenotype$mcv=gsub("mcv: (.*)","\\1",phenotype$mcv)
phenotype$mch=gsub("mch: (.*)","\\1",phenotype$'mch:')
phenotype$mchc=gsub("mchc: (.*)","\\1",phenotype$mchc)
phenotype$rdw=gsub("rdw: (.*)","\\1",phenotype$rdw)
phenotype$mpv=gsub("mpv: (.*)","\\1",phenotype$mpv)
phenotype$cr=gsub("cr: (.*)","\\1",phenotype$'cr:')
phenotype$alb=gsub("alb: (.*)","\\1",phenotype$alb)
phenotype$ds_dna=gsub("ds_dna: (.*)","\\1",phenotype$ds_dna)
phenotype$ast=gsub("ast: (.*)","\\1",phenotype$ast)
phenotype$ald=gsub("ald: (.*)","\\1",phenotype$ald)
phenotype$alt=gsub("alt: (.*)","\\1",phenotype$alt)
phenotype$ldh=gsub("ldh: (.*)","\\1",phenotype$ldh)
phenotype$sledai=gsub("sledai: (.*)","\\1",phenotype$sledai)
phenotype$platelet_count=gsub("platelet_count: (.*)","\\1",phenotype$platelet_count)

#New phenotypes categorical
phenotype$hematuria=gsub("hematuria: (.*)","\\1",phenotype$'hematuria:')
phenotype$proteinuria=gsub("proteinuria: (.*)","\\1",phenotype$proteinuria)
phenotype$pyuria=gsub("pyuria: (.*)","\\1",phenotype$pyuria)
phenotype$new_rash=gsub("new_rash: (.*)","\\1",phenotype$new_rash)
phenotype$alopecia=gsub("alopecia: (.*)","\\1",phenotype$alopecia)
phenotype$mucosal_ulcers=gsub("mucosal_ulcers: (.*)","\\1",phenotype$mucosal_ulcers)
phenotype$pleurisy=gsub("pleurisy: (.*)","\\1",phenotype$pleurisy)
phenotype$pericarditis=gsub("pericarditis: (.*)","\\1",phenotype$pericarditis)
phenotype$low_complement=gsub("low_complement: (.*)","\\1",phenotype$low_complement)
phenotype$increased_dna_binding=gsub("increased_dna_binding: (.*)","\\1",phenotype$increased_dna_binding)
#phenotype$increased_dna_binding=gsub("increased_dna_binding: (.*)","\\1",phenotype$increased_dna_binding)
phenotype$fever=gsub("fever: (.*)","\\1",phenotype$fever)
phenotype$thrombocytopenia=gsub("thrombocytopenia: (.*)","\\1",phenotype$thrombocytopenia)
phenotype$leukopenia=gsub("leukopenia: (.*)","\\1",phenotype$leukopenia)
phenotype$renal=gsub("renal: (.*)","\\1",phenotype$renal)
phenotype$musculoskeletal=gsub("musculoskeletal: (.*)","\\1",phenotype$musculoskeletal)
phenotype$nephritis_class=gsub("nephritis_class: (.*)","\\1",phenotype$nephritis_class)
phenotype$cyclophosphamide_category=gsub("cyclophosphamide_category: (.*)","\\1",phenotype$cyclophosphamide_category)
phenotype$oral_steroids_category=gsub("oral_steroids_category: (.*)","\\1",phenotype$oral_steroids_category)
phenotype$mycophenolate_category=gsub("mycophenolate_category: (.*)","\\1",phenotype$mycophenolate_category)
phenotype$hydroxychloroquine_category=gsub("hydroxychloroquine_category: (.*)","\\1",phenotype$hydroxychloroquine_category)
phenotype$metotrexate_category=gsub("metotrexate_category: (.*)","\\1",phenotype$metotrexate_category)
phenotype$nsaid_category=gsub("nsaid_category: (.*)","\\1",phenotype$nsaid_category)
phenotype$asa_category=gsub("asa_category: (.*)","\\1",phenotype$asa_category)
phenotype$mdg =gsub("mdg: (.*)","\\1",phenotype$mdg)
#all 0's
#phenotype$seizure =gsub("seizure: (.*)","\\1",phenotype$seizure)
phenotype$psychosis  =gsub("psychosis: (.*)","\\1",phenotype$psychosis)
phenotype$organic_brain_syndrome  =gsub("organic_brain_syndrome: (.*)","\\1",phenotype$organic_brain_syndrome)
phenotype$lupus_headache  =gsub("lupus_headache: (.*)","\\1",phenotype$lupus_headache)
phenotype$cva  =gsub("cva: (.*)","\\1",phenotype$cva)
phenotype$vasculitis =gsub("vasculitis: (.*)","\\1",phenotype$vasculitis)
phenotype$arthritis =gsub("arthritis: (.*)","\\1",phenotype$arthritis)
phenotype$myositis=gsub("myositis: (.*)","\\1",phenotype$myositis)
phenotype$urinary_casts=gsub("urinary_casts: (.*)","\\1",phenotype$urinary_casts)
phenotype$neph_treat_lmm3=gsub("neph_treat_lmm3:(.*)","\\1",phenotype$neph_treat_lmm3)
#phenotype$visual_disturbance  =gsub("visual_disturbance: (.*)","\\1",phenotype$visual_disturbance)
#phenotype$cranial_nerve_disorder  =gsub("cranial_nerve_disorder: (.*)","\\1",phenotype$cranial_nerve_disorder)




indexc3=which(phenotype$c3!="Data Not Available")
indexNC=which(phenotype$neutrophil_count!="Data Not Available")
indexNP=which(phenotype$neutrophil_percent!="Data Not Available")
indexLC=which(phenotype$lymphocyte_count!="Data Not Available")
indexLP=which(phenotype$lymphocyte_percent!="Data Not Available")
indexc4=which(phenotype$c4!="Data Not Available")
indexwbc=which(phenotype$wbc!="Data Not Available")
indexMC=which(!(grepl("Data Not Available",phenotype$MC)))
indexMP=which(!(grepl("Data Not Available",phenotype$MP)))
indexESR=which(phenotype$esr!="Data Not Available")
indexhgb=which(phenotype$hgb!="Data Not Available")
indexhct=which(phenotype$hct!="Data Not Available")
indexmcv=which(phenotype$mcv!="Data Not Available")
indexmch=which(phenotype$mch!="Data Not Available")
indexmchc=which(phenotype$mchc!="Data Not Available")
indexrdw=which(phenotype$rdw!="Data Not Available")
indexmpv=which(phenotype$mpv!="Data Not Available")
indexcr=which(phenotype$cr!="Data Not Available")
indexalb=which(phenotype$alb!="Data Not Available")
indexds_dna=which(phenotype$ds_dna!="Data Not Available"&phenotype$ds_dna!="Verify at TSRH"&phenotype$ds_dna!="(-)"&phenotype$ds_dna!="ND"&!grepl("[0-9]:[0-9]",phenotype$ds_dna))
indexast=which(phenotype$ast!="Data Not Available")
indexald=which(phenotype$ald!="Data Not Available")
indexalt=which(phenotype$alt!="Data Not Available")
indexldh=which(phenotype$ldh!="Data Not Available")
indexsledai=which(phenotype$sledai!="Data Not Available")
indexage=c(1:length(phenotype$age))
indexplatelet_count=which(phenotype$platelet_count!="Data Not Available")

#categorical
indexhydroxychloroquine_category=which(phenotype$hydroxychloroquine_category!="Data Not Available")
indexmetotrexate_category=which(phenotype$metotrexate_category!="Data Not Available")
indexnsaid_category=which(phenotype$nsaid_category!="Data Not Available")
indexasa_category=which(phenotype$asa_category!="Data Not Available")
indexmdg=which(phenotype$mdg!="Data Not Available"& phenotype$mdg!="N/A")
indexmdg=which(phenotype$mdg!="Data Not Available"& phenotype$mdg!="N/A")
indextreatment=which(phenotype$neph_treat_lmm3!=" Data Not Available"& phenotype$neph_treat_lmm3!=" Not Applicable")
indexcommon=c(1:dim(phenotype)[1])

indicesPhenotypes=NULL
indicesPhenotypes=list(as.data.frame(indexc3),as.data.frame(indexc4),as.data.frame(indexNC),as.data.frame(indexNP),as.data.frame(indexLC),as.data.frame(indexLP),as.data.frame(indexwbc),as.data.frame(indexMC),as.data.frame(indexMP),as.data.frame(indexESR),as.data.frame(indexhgb),as.data.frame(indexhct),as.data.frame(indexmcv),as.data.frame(indexmch),as.data.frame(indexmchc),as.data.frame(indexrdw),as.data.frame(indexmpv),as.data.frame(indexcr),as.data.frame(indexalb),as.data.frame(indexds_dna),as.data.frame(indexast),as.data.frame(indexald),as.data.frame(indexalt),as.data.frame(indexldh),as.data.frame(indexsledai),as.data.frame(indexage),as.data.frame(indexplatelet_count))

names(indicesPhenotypes)<-c("c3","c4","neutrophil_count","neutrophil_percent","lymphocyte_count","lymphocyte_percent","wbc","MC","MP","esr","hgb","hct","mcv","mch","mchc","rdw","mpv","cr","alb","ds_dna","ast","ald","alt","ldh","sledai","age","platelet_count")

indicesPhenotypesCategorical=list(as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexhydroxychloroquine_category),
                                  as.data.frame(indexmetotrexate_category)
                                  ,as.data.frame(indexnsaid_category),
                                  as.data.frame(indexasa_category),
                                  as.data.frame(indexmdg),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indexcommon),
                                  as.data.frame(indextreatment))

names(indicesPhenotypesCategorical)<-c("hematuria",
                                       "proteinuria",
                                       "pyuria",
                                       "new_rash",
                                       "alopecia",
                                       "pleurisy",
                                       "mucosal_ulcers",
                                       "pericarditis",
                                       "low_complement",
                                       "increased_dna_binding",
                                       "fever",
                                       "thrombocytopenia",
                                       "leukopenia",
                                       "musculoskeletal",
                                       "cyclophosphamide_category",
                                       "oral_steroids_category",
                                       "mycophenolate_category",
                                       "hydroxychloroquine_category",
                                       "metotrexate_category",
                                       "nsaid_category",
                                       "asa_category",
                                       "mdg",
                                       "psychosis",
                                       "organic_brain_syndrome",
                                       "lupus_headache",
                                       "cva",
                                       "vasculitis",
                                       "arthritis",
                                       "myositis",
                                       "urinary_casts",
                                       "gender",
                                       "race",
                                       "neph_treat_lmm3")


