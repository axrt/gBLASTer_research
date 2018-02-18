# Grant Report Ending 2015
Alexander Tuzhikov  
December 21, 2015  

#Report

##Libraries

```r
source("R/helper.R")
library(gbra)
source("R/pander.lib.R")
source("R/yr.lib.R")
source("R/data.table.lib.R")
source("R/reshape2.lib.R")
source("R/data.table.lib.R")
source("R/scales.lib.R")
report.2015<- function(file=""){
  return(paste("report.2015", file, sep="/"))
}

if(!file.exists(report.2015())){
  dir.create(report.2015())
}
```

##A list of all organisms from the database:


```r
legend<- read.table(file=sixtyseven.genomes.dir("legend.txt"), sep = "\t", stringsAsFactors = FALSE, header = TRUE)
legend %>% 
  mutate(nubmer=1:nrow(legend)) %>% 
  rename(`database id`=id_genomes, `genome name`=name) -> report.legend
pander(report.legend, caption = "Genome List")
```


-------------
 database id 
-------------
      1      

      2      

      3      

      4      

      5      

      6      

      7      

      8      

      9      

     10      

     11      

     12      

     13      

     14      

     15      

     16      

     17      

     18      

     19      

     20      

     21      

     22      

     23      

     24      

     25      

     26      

     27      

     28      

     30      

     31      

     32      

     33      

     34      

     35      

     36      

     37      

     38      

     39      

     40      

     41      

     42      

     43      

     44      

     45      

     46      

     47      

     48      

     49      

     50      

     51      

     52      

     53      

     54      

     55      

     56      

     57      

     58      

     59      

     60      

     62      

     63      

     64      

     65      

     66      

     67      

     68      

     69      
-------------

Table: Genome List (continued below)

 
-------------------------------------------------------------------------------
                                  genome name                                  
-------------------------------------------------------------------------------
            Escherichia_coli_str_K_12_substr_MG1655_complete_genome            

      Bacillus_subtilis_subsp_subtilis_str_168_chromosome_complete_genome      

       Wolbachia_endosymbiont_of_Drosophila_simulans_wHa_complete_genome       

Wolbachia_endosymbiont_of_Culex_quinquefasciatus_Pel_chromosome_complete_genome

             Sphingomonas_wittichii_RW1_chromosome_complete_genome             

                              Drosophila_simulans                              

                             Arabidopsis_thaliana                              

                              Ostreococcus_tauri                               

                       Caenorhabditis_Elegans_Bristol_N2                       

                        Saccharomyces_cerevisiae_S288c                         

                        Cryptosporidium_parvum_Iowa_II                         

                        Schizosaccharomyces_pombe_972h                         

                       Thalassiosira_pseudonana_CCMP1335                       

                          Anopheles_gambiae_str_PEST                           

                        Encephalitozoon_cuniculi_GB_M1                         

                       Methanosarcina_barkeri_str_Fusaro                       

                               Takifugu_rubripes                               

                             Haloquadratum_walsbyi                             

                           Acidilobus_saccharovorans                           

                            Cryptococcus_Neoformans                            

                              Tribolium_castaneum                              

                              Ciona_intestinalis                               

                      Cyanidioschyzon_merolae_strain_10D                       

                        Pyrobaculum_aerophilum_str_IM2                         

                          Methanobacterium_sp_SWAN_1                           

                             Monosiga_brevicollis                              

                             Capsaspora_owczarzaki                             

                       Bartonella_henselae_str_Houston_1                       

                        Schistosoma_mansoni_puerto_rico                        

                   Fusarium_oxysporum_lycopersici_4287.fasta                   

                     Allomyces_macrogynus_ATCC_38327.fasta                     

                       Rubrobacter_xylanophilus_DSM_9941                       

                           Mycoplasma_genitalium_G37                           

                           Mycoplasma_gallisepticum                            

                       Mycobacterium_tuberculosis_H37Rv                        

                             Nema_parisii_ERTm1_V3                             

                            Cyanidioschyzon_merolae                            

                           Plasmodium_falciparum_3D7                           

                         Ichthyophthirius_multifiliis                          

                                Genlisea_aurea                                 

                               Naegleria_fowleri                               

                            Gregarina_niphandrodes                             

                                 Culex_pipiens                                 

                            Archaeon_Loki_Lokiarch                             

                             Prunus_persica_Lowell                             

                        Drosophila_simulans_strain_w501                        

                       Synechococcus_elongatus_PCC_7942                        

                       Prochlorococcus_marinus_CCMP1375                        

                        Thermotoga_neapolitana_DSM_4359                        

                           Thermus_aquaticus_Y51MC23                           

                             Aquifex_aeolicus_VF5                              

                               Oikopleura_dioica                               

                         Dictyostelium_discoideum_AX4                          

                          Entamoeba_histolytica_HM_1                           

                            Rozella_allomycis_CSF55                            

                         Meloidogyne_hapla_strain_VW9                          

                             Pandoravirus_salinus                              

                       Acanthamoeba_polyphaga_mimivirus                        

                             Megavirus_chiliensis                              

                        Rickettsia_prowazekii_Madrid_E                         

                             Meloidogyne_hapla_VW9                             

                          Chondrus_crispus_Stackhouse                          

                         Leishmania_donovani_BPK282A1                          

                              Leishmania_mexicana                              

                             Micromonas_sp_RCC299                              

                              Belgica_antarctica                               

                              Ordospora_colligata                              
-------------------------------------------------------------------------------

Table: Table continues below

 
--------
 nubmer 
--------
   1    

   2    

   3    

   4    

   5    

   6    

   7    

   8    

   9    

   10   

   11   

   12   

   13   

   14   

   15   

   16   

   17   

   18   

   19   

   20   

   21   

   22   

   23   

   24   

   25   

   26   

   27   

   28   

   29   

   30   

   31   

   32   

   33   

   34   

   35   

   36   

   37   

   38   

   39   

   40   

   41   

   42   

   43   

   44   

   45   

   46   

   47   

   48   

   49   

   50   

   51   

   52   

   53   

   54   

   55   

   56   

   57   

   58   

   59   

   60   

   61   

   62   

   63   

   64   

   65   

   66   

   67   
--------

```r
write.table(x=report.legend, file=report.2015("legend.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
```

## Core/non-core distribution per gene

```r
load(sixtyseven.genomes.dir("core.distribution.data.rda"))
expand.df(core.distribution.data) -> core.distribution.data
core.distribution.data.melt<- melt(core.distribution.data, id.vars = c("QUERY_ORF_ID","ID_QUERY_GENOME"))
core.distribution.data.melt %>% group_by(ID_QUERY_GENOME, QUERY_ORF_ID, value) %>% summarize(total.grouping=n()) -> core.distribution.data.melt.summary
#always core
instate<- function(df){
  cores<- df %>% filter(value=="core") %>% select(total.grouping)
  outsiders<- df %>% filter(value=="outsider") %>% select(total.grouping)
  noseps<- df %>% filter(value=="nosep") %>% select(total.grouping)
  if(noseps==61){
    return(data.frame(ID_QUERY_GENOME= df$ID_QUERY_GENOME[1], QUERY_ORF_ID=df$QUERY_ORF_ID[1], state="never separates"))
  }
  if(cores+noseps==61){
    return(data.frame(ID_QUERY_GENOME= df$ID_QUERY_GENOME[1], QUERY_ORF_ID=df$QUERY_ORF_ID[1], state="always core"))
  }
  if(cores+noseps==60){
    return(data.frame(ID_QUERY_GENOME= df$ID_QUERY_GENOME[1], QUERY_ORF_ID=df$QUERY_ORF_ID[1], state="outsider once"))
  }else{
    return(data.frame(ID_QUERY_GENOME= df$ID_QUERY_GENOME[1], QUERY_ORF_ID=df$QUERY_ORF_ID[1], state="outsider multiple"))
  }
}
if(!file.exists(report.2015("df.list.rda"))){
  df.list<- list()
  core.distribution.data.melt.summary %>% 
    group_by(ID_QUERY_GENOME, QUERY_ORF_ID) %>% 
    do({
      dl<-instate(.)
      df.list[[length(df.list)+1]]<<- dl
      return(NA)
    })
  save(df.list, file = report.2015("df.list.rda"))
}else{
  load(file=report.2015("df.list.rda"))
}

#turns out this is the fastes way
if(!file.exists(report.2015("core.distribution.data.melt.summary.total.txt"))){
  for(i in 1:length(df.list)){
    line<- paste(df.list[[i]][[1]], df.list[[i]][[2]], as.character(df.list[[i]][[3]]), sep="\t", collapse = "\n")
    write(line,file=report.2015("core.distribution.data.melt.summary.total.txt"), append=TRUE)
  }
}
core.distribution.data.melt.summary.total<- fread(input = report.2015("core.distribution.data.melt.summary.total.txt"), sep = "\t", data.table = FALSE, col.names = c("ID_QUERY_GENOME", "QUERY_ORF_ID", "state"))

core.distribution.data.melt.summary.total %>%
  group_by(ID_QUERY_GENOME, state) %>% 
  summarize(ammount=n()) %>%
  group_by(ID_QUERY_GENOME) %>%
  do({
    df<-.
    df$ratio<- rescale(.$ammount,from = c(0,sum(.$ammount)), to=c(0, 100))
    return(df)
  }) %>%
  merge(x=., y=legend, by.x="ID_QUERY_GENOME", by.y="id_genomes")-> core.distribution.data.melt.summary.total.ratio
write.table(x = core.distribution.data.melt.summary.total.ratio, file=report.2015("core.distribution.data.melt.summary.total.ratio.txt"), sep="\t", quote = FALSE)
pander(core.distribution.data.melt.summary.total.ratio, caption = "Distribution of gene occurances in core and non-core")
```


-----------------------------------------------------
 ID_QUERY_GENOME        state        ammount   ratio 
----------------- ----------------- --------- -------
        1            always core      1619     75.72 

        1         outsider multiple    364     17.03 

        1           outsider once      155     7.25  

        2            always core      1372     74.36 

        2         outsider multiple    342     18.54 

        2           outsider once      131      7.1  

        3            always core       555     83.21 

        3         outsider multiple    74      11.09 

        3           outsider once      38      5.697 

        4            always core       607     83.38 

        4         outsider multiple    85      11.68 

        4           outsider once      36      4.945 

        5            always core      1998     75.71 

        5         outsider multiple    402     15.23 

        5           outsider once      239     9.056 

        6            always core      1537     76.43 

        6         outsider multiple    353     17.55 

        6           outsider once      121     6.017 

        7            always core      10579    89.71 

        7         outsider multiple    913     7.742 

        7           outsider once      301     2.552 

        8            always core      3415     70.69 

        8         outsider multiple    725     15.01 

        8           outsider once      691     14.3  

        9            always core      7120     85.42 

        9         outsider multiple    709     8.506 

        9           outsider once      506     6.071 

       10            always core      2556     74.71 

       10         outsider multiple    633     18.5  

       10           outsider once      232     6.782 

       11            always core      1167     70.3  

       11         outsider multiple    198     11.93 

       11           outsider once      295     17.77 

       12            always core      2524     75.39 

       12         outsider multiple    505     15.08 

       12           outsider once      319     9.528 

       13            always core      5523     94.47 

       13         outsider multiple    106     1.813 

       13           outsider once      217     3.712 

       14            always core      17988    92.59 

       14         outsider multiple    869     4.473 

       14           outsider once      570     2.934 

       15            always core       798     81.6  

       15         outsider multiple    150     15.34 

       15           outsider once      30      3.067 

       16            always core      1049     71.17 

       16         outsider multiple    340     23.07 

       16           outsider once      85      5.767 

       17            always core      10409    88.75 

       17         outsider multiple    710     6.053 

       17           outsider once      610     5.201 

       18            always core       704     67.43 

       18         outsider multiple    251     24.04 

       18           outsider once      89      8.525 

       19            always core       443     65.53 

       19         outsider multiple    192     28.4  

       19           outsider once      41      6.065 

       20            always core      2232     63.9  

       20         outsider multiple    849     24.31 

       20           outsider once      412     11.8  

       21            always core      12460    92.25 

       21         outsider multiple    628     4.649 

       21           outsider once      419     3.102 

       22            always core      2198     75.09 

       22         outsider multiple    479     16.36 

       22           outsider once      250     8.541 

       24            always core       557     67.6  

       24         outsider multiple    203     24.64 

       24           outsider once      64      7.767 

       25            always core       713     72.02 

       25         outsider multiple    228     23.03 

       25           outsider once      49      4.949 

       26            always core      4803     92.21 

       26         outsider multiple    179     3.436 

       26           outsider once      227     4.358 

       27            always core      5823     93.75 

       27         outsider multiple    187     3.011 

       27           outsider once      201     3.236 

       28            always core       549     64.74 

       28         outsider multiple    242     28.54 

       28           outsider once      57      6.722 

       30            always core      11620    93.24 

       30         outsider multiple    554     4.446 

       30           outsider once      288     2.311 

       31            always core      6912     87.89 

       31         outsider multiple    493     6.269 

       31           outsider once      459     5.837 

       32            always core      10423    89.29 

       32         outsider multiple    222     1.902 

       32           outsider once     1028     8.807 

       33            always core      1447     72.6  

       33         outsider multiple    338     16.96 

       33           outsider once      208     10.44 

       34            always core       159     58.89 

       34         outsider multiple    91      33.7  

       34           outsider once      20      7.407 

       35            always core       199     61.42 

       35         outsider multiple    101     31.17 

       35           outsider once      24      7.407 

       36            always core      1236     69.63 

       36         outsider multiple    379     21.35 

       36           outsider once      160     9.014 

       37            always core       486     56.45 

       37         outsider multiple    267     31.01 

       37           outsider once      108     12.54 

       38            always core      2852     90.4  

       38         outsider multiple    61      1.933 

       38           outsider once      242     7.67  

       39            always core      1846     77.08 

       39         outsider multiple    338     14.11 

       39           outsider once      211     8.81  

       41            always core      5973     82.79 

       41         outsider multiple    979     13.57 

       41           outsider once      263     3.645 

       42            always core      4506     85.03 

       42         outsider multiple    174     3.284 

       42           outsider once      619     11.68 

       43            always core      1983     87.78 

       43         outsider multiple    161     7.127 

       43           outsider once      115     5.091 

       44            always core      46173    96.29 

       44         outsider multiple    817     1.704 

       44           outsider once      960     2.002 

       45            always core      1058     65.51 

       45         outsider multiple    347     21.49 

       45           outsider once      210      13   

       46            always core      16396    91.76 

       46         outsider multiple   1107     6.195 

       46           outsider once      365     2.043 

       48            always core      1030     70.02 

       48         outsider multiple    324     22.03 

       48           outsider once      117     7.954 

       49            always core       645     67.33 

       49         outsider multiple    224     23.38 

       49           outsider once      89      9.29  

       50            always core       773     67.28 

       50         outsider multiple    293     25.5  

       50           outsider once      83      7.224 

       51            always core       96      43.44 

       51         outsider multiple    97      43.89 

       51           outsider once      28      12.67 

       52            always core       626     66.38 

       52         outsider multiple    245     25.98 

       52           outsider once      72      7.635 

       53            always core      5383     85.72 

       53         outsider multiple    522     8.312 

       53           outsider once      375     5.971 

       54            always core      5293     93.52 

       54         outsider multiple    145     2.562 

       54           outsider once      222     3.922 

       55            always core      2664     79.67 

       55         outsider multiple    218     6.519 

       55           outsider once      462     13.82 

       56            always core      3059     88.87 

       56         outsider multiple    217     6.304 

       56           outsider once      166     4.823 

       58            always core       115     45.28 

       58         outsider multiple    95      37.4  

       58           outsider once      44      17.32 

       59            always core       144     60.5  

       59         outsider multiple    78      32.77 

       59           outsider once      16      6.723 

       60            always core       176     58.28 

       60         outsider multiple    104     34.44 

       60           outsider once      22      7.285 

       62            always core       320     57.97 

       62         outsider multiple    177     32.07 

       62           outsider once      55      9.964 

       65            always core      3598     94.58 

       65         outsider multiple    142     3.733 

       65           outsider once      64      1.682 

       66            always core      3775     89.33 

       66         outsider multiple    385     9.11  

       66           outsider once      66      1.562 

       67            always core      5017     78.86 

       67         outsider multiple    729     11.46 

       67           outsider once      616     9.682 

       68            always core      11890    91.61 

       68         outsider multiple    616     4.746 

       68           outsider once      473     3.644 

       69            always core       785     89.51 

       69         outsider multiple    62      7.07  

       69           outsider once      30      3.421 
-----------------------------------------------------

Table: Distribution of gene occurances in core and non-core (continued below)

 
-------------------------------------------------------------------------------
                                     name                                      
-------------------------------------------------------------------------------
            Escherichia_coli_str_K_12_substr_MG1655_complete_genome            

            Escherichia_coli_str_K_12_substr_MG1655_complete_genome            

            Escherichia_coli_str_K_12_substr_MG1655_complete_genome            

      Bacillus_subtilis_subsp_subtilis_str_168_chromosome_complete_genome      

      Bacillus_subtilis_subsp_subtilis_str_168_chromosome_complete_genome      

      Bacillus_subtilis_subsp_subtilis_str_168_chromosome_complete_genome      

       Wolbachia_endosymbiont_of_Drosophila_simulans_wHa_complete_genome       

       Wolbachia_endosymbiont_of_Drosophila_simulans_wHa_complete_genome       

       Wolbachia_endosymbiont_of_Drosophila_simulans_wHa_complete_genome       

Wolbachia_endosymbiont_of_Culex_quinquefasciatus_Pel_chromosome_complete_genome

Wolbachia_endosymbiont_of_Culex_quinquefasciatus_Pel_chromosome_complete_genome

Wolbachia_endosymbiont_of_Culex_quinquefasciatus_Pel_chromosome_complete_genome

             Sphingomonas_wittichii_RW1_chromosome_complete_genome             

             Sphingomonas_wittichii_RW1_chromosome_complete_genome             

             Sphingomonas_wittichii_RW1_chromosome_complete_genome             

                              Drosophila_simulans                              

                              Drosophila_simulans                              

                              Drosophila_simulans                              

                             Arabidopsis_thaliana                              

                             Arabidopsis_thaliana                              

                             Arabidopsis_thaliana                              

                              Ostreococcus_tauri                               

                              Ostreococcus_tauri                               

                              Ostreococcus_tauri                               

                       Caenorhabditis_Elegans_Bristol_N2                       

                       Caenorhabditis_Elegans_Bristol_N2                       

                       Caenorhabditis_Elegans_Bristol_N2                       

                        Saccharomyces_cerevisiae_S288c                         

                        Saccharomyces_cerevisiae_S288c                         

                        Saccharomyces_cerevisiae_S288c                         

                        Cryptosporidium_parvum_Iowa_II                         

                        Cryptosporidium_parvum_Iowa_II                         

                        Cryptosporidium_parvum_Iowa_II                         

                        Schizosaccharomyces_pombe_972h                         

                        Schizosaccharomyces_pombe_972h                         

                        Schizosaccharomyces_pombe_972h                         

                       Thalassiosira_pseudonana_CCMP1335                       

                       Thalassiosira_pseudonana_CCMP1335                       

                       Thalassiosira_pseudonana_CCMP1335                       

                          Anopheles_gambiae_str_PEST                           

                          Anopheles_gambiae_str_PEST                           

                          Anopheles_gambiae_str_PEST                           

                        Encephalitozoon_cuniculi_GB_M1                         

                        Encephalitozoon_cuniculi_GB_M1                         

                        Encephalitozoon_cuniculi_GB_M1                         

                       Methanosarcina_barkeri_str_Fusaro                       

                       Methanosarcina_barkeri_str_Fusaro                       

                       Methanosarcina_barkeri_str_Fusaro                       

                               Takifugu_rubripes                               

                               Takifugu_rubripes                               

                               Takifugu_rubripes                               

                             Haloquadratum_walsbyi                             

                             Haloquadratum_walsbyi                             

                             Haloquadratum_walsbyi                             

                           Acidilobus_saccharovorans                           

                           Acidilobus_saccharovorans                           

                           Acidilobus_saccharovorans                           

                            Cryptococcus_Neoformans                            

                            Cryptococcus_Neoformans                            

                            Cryptococcus_Neoformans                            

                              Tribolium_castaneum                              

                              Tribolium_castaneum                              

                              Tribolium_castaneum                              

                              Ciona_intestinalis                               

                              Ciona_intestinalis                               

                              Ciona_intestinalis                               

                        Pyrobaculum_aerophilum_str_IM2                         

                        Pyrobaculum_aerophilum_str_IM2                         

                        Pyrobaculum_aerophilum_str_IM2                         

                          Methanobacterium_sp_SWAN_1                           

                          Methanobacterium_sp_SWAN_1                           

                          Methanobacterium_sp_SWAN_1                           

                             Monosiga_brevicollis                              

                             Monosiga_brevicollis                              

                             Monosiga_brevicollis                              

                             Capsaspora_owczarzaki                             

                             Capsaspora_owczarzaki                             

                             Capsaspora_owczarzaki                             

                       Bartonella_henselae_str_Houston_1                       

                       Bartonella_henselae_str_Houston_1                       

                       Bartonella_henselae_str_Houston_1                       

                        Schistosoma_mansoni_puerto_rico                        

                        Schistosoma_mansoni_puerto_rico                        

                        Schistosoma_mansoni_puerto_rico                        

                   Fusarium_oxysporum_lycopersici_4287.fasta                   

                   Fusarium_oxysporum_lycopersici_4287.fasta                   

                   Fusarium_oxysporum_lycopersici_4287.fasta                   

                     Allomyces_macrogynus_ATCC_38327.fasta                     

                     Allomyces_macrogynus_ATCC_38327.fasta                     

                     Allomyces_macrogynus_ATCC_38327.fasta                     

                       Rubrobacter_xylanophilus_DSM_9941                       

                       Rubrobacter_xylanophilus_DSM_9941                       

                       Rubrobacter_xylanophilus_DSM_9941                       

                           Mycoplasma_genitalium_G37                           

                           Mycoplasma_genitalium_G37                           

                           Mycoplasma_genitalium_G37                           

                           Mycoplasma_gallisepticum                            

                           Mycoplasma_gallisepticum                            

                           Mycoplasma_gallisepticum                            

                       Mycobacterium_tuberculosis_H37Rv                        

                       Mycobacterium_tuberculosis_H37Rv                        

                       Mycobacterium_tuberculosis_H37Rv                        

                             Nema_parisii_ERTm1_V3                             

                             Nema_parisii_ERTm1_V3                             

                             Nema_parisii_ERTm1_V3                             

                            Cyanidioschyzon_merolae                            

                            Cyanidioschyzon_merolae                            

                            Cyanidioschyzon_merolae                            

                           Plasmodium_falciparum_3D7                           

                           Plasmodium_falciparum_3D7                           

                           Plasmodium_falciparum_3D7                           

                                Genlisea_aurea                                 

                                Genlisea_aurea                                 

                                Genlisea_aurea                                 

                               Naegleria_fowleri                               

                               Naegleria_fowleri                               

                               Naegleria_fowleri                               

                            Gregarina_niphandrodes                             

                            Gregarina_niphandrodes                             

                            Gregarina_niphandrodes                             

                                 Culex_pipiens                                 

                                 Culex_pipiens                                 

                                 Culex_pipiens                                 

                            Archaeon_Loki_Lokiarch                             

                            Archaeon_Loki_Lokiarch                             

                            Archaeon_Loki_Lokiarch                             

                             Prunus_persica_Lowell                             

                             Prunus_persica_Lowell                             

                             Prunus_persica_Lowell                             

                       Synechococcus_elongatus_PCC_7942                        

                       Synechococcus_elongatus_PCC_7942                        

                       Synechococcus_elongatus_PCC_7942                        

                       Prochlorococcus_marinus_CCMP1375                        

                       Prochlorococcus_marinus_CCMP1375                        

                       Prochlorococcus_marinus_CCMP1375                        

                        Thermotoga_neapolitana_DSM_4359                        

                        Thermotoga_neapolitana_DSM_4359                        

                        Thermotoga_neapolitana_DSM_4359                        

                           Thermus_aquaticus_Y51MC23                           

                           Thermus_aquaticus_Y51MC23                           

                           Thermus_aquaticus_Y51MC23                           

                             Aquifex_aeolicus_VF5                              

                             Aquifex_aeolicus_VF5                              

                             Aquifex_aeolicus_VF5                              

                               Oikopleura_dioica                               

                               Oikopleura_dioica                               

                               Oikopleura_dioica                               

                         Dictyostelium_discoideum_AX4                          

                         Dictyostelium_discoideum_AX4                          

                         Dictyostelium_discoideum_AX4                          

                          Entamoeba_histolytica_HM_1                           

                          Entamoeba_histolytica_HM_1                           

                          Entamoeba_histolytica_HM_1                           

                            Rozella_allomycis_CSF55                            

                            Rozella_allomycis_CSF55                            

                            Rozella_allomycis_CSF55                            

                             Pandoravirus_salinus                              

                             Pandoravirus_salinus                              

                             Pandoravirus_salinus                              

                       Acanthamoeba_polyphaga_mimivirus                        

                       Acanthamoeba_polyphaga_mimivirus                        

                       Acanthamoeba_polyphaga_mimivirus                        

                             Megavirus_chiliensis                              

                             Megavirus_chiliensis                              

                             Megavirus_chiliensis                              

                        Rickettsia_prowazekii_Madrid_E                         

                        Rickettsia_prowazekii_Madrid_E                         

                        Rickettsia_prowazekii_Madrid_E                         

                         Leishmania_donovani_BPK282A1                          

                         Leishmania_donovani_BPK282A1                          

                         Leishmania_donovani_BPK282A1                          

                              Leishmania_mexicana                              

                              Leishmania_mexicana                              

                              Leishmania_mexicana                              

                             Micromonas_sp_RCC299                              

                             Micromonas_sp_RCC299                              

                             Micromonas_sp_RCC299                              

                              Belgica_antarctica                               

                              Belgica_antarctica                               

                              Belgica_antarctica                               

                              Ordospora_colligata                              

                              Ordospora_colligata                              

                              Ordospora_colligata                              
-------------------------------------------------------------------------------

##Separation ratio per genome

```r
model.df<- read.table(file=sixtyseven.genomes.dir("model.df.txt"), sep="\t", header=TRUE)
model.df %>% select(query.genome, target.genome, rsquared) %>% 
  mutate(separation=rsquared>=0.73) %>% select(query.genome, separation) %>%
  group_by(query.genome, separation) %>% summarize(separation.sum=n()) %>% group_by(query.genome) %>% do({
    df<-.
    df$separation.ratio<- rescale(df$separation.sum, from=c(0,sum(df$separation.sum)), to=c(0,100))
    return(df)
  }) -> model.df
write.table(x=model.df, file=report.2015("separation.per.genome.txt"), sep="\t", quote = FALSE)
```

##Drosophila vs Wolbachia Plot Update


```r
legend %>% filter(id_genomes!=23) %>% 
  filter(id_genomes!=40) %>% 
  filter(id_genomes!=47) %>%
  filter(id_genomes!=57) %>%
  filter(id_genomes!=63) %>%
  filter(id_genomes!=64) %>%
  select(name, id_genomes) -> legend
load(sixtyseven.genomes.dir("bh.data.normal.rda"))
gp<- gen.pca(df = select.genomes(g.ids = c(3,6), df = bh.data.norm), legend = legend, g.ids = c(3, 6))
```

```
## Loading required package: devtools
## Loading required package: ggbiplot
## Loading required package: ggplot2
## Loading required package: grid
```

```r
gp
```

![](grant.report.end.2015_files/figure-html/dros wold plot-1.png) 

```r
png(filename = report.2015("dros.wolb.pca.png"), width = 700, heigh= 700)
gp
dev.off()
```

```
## png 
##   2
```

##Amoeba and its viruses


```r
gp<- gen.pca(df = select.genomes(g.ids = c(55, 58, 59), df = bh.data.norm), legend = legend, g.ids = c(55, 58, 59))
gp
```

![](grant.report.end.2015_files/figure-html/amoeba virus-1.png) 

```r
png(filename = report.2015("amoeb.virus.pca.png"), width = 700, heigh= 700)
gp
dev.off()
```

```
## png 
##   2
```

##Drosophila separation plot


```r
model.df<-read.table(file = sixtyseven.genomes.dir("model.df.txt"), sep="\t", stringsAsFactors = FALSE, header = TRUE)
#plot the distribution of rsquareds
rsq.plot<- model.df %>% arrange(rsquared) %>% ggplot(data=., mapping=aes(x=rescale_max(as.numeric(factor(file.name, levels=file.name))), y= rsquared)) + 
  geom_line(color="red") + xlab("Model number (scaled)") + ylab("R squared") + theme_bw()
rsq.plot
```

![](grant.report.end.2015_files/figure-html/dros distro-1.png) 

```r
png(file=report.2015("rsq.plot.png"), width = 700, height= 700)
rsq.plot
dev.off()
```

```
## png 
##   2
```

```r
rsq.bar<- core.distribution.data %>% select.genomes(df=., g.ids=6) %>% mutate(id=rownames(.)) %>% melt(id.vars="id") %>% group_by(variable,value) %>%
  summarise(times=n()) %>% do({
    .$times<- rescale(.$times, from = c(0,sum(.$times)),to = c(0,1))
    return(.)
  }) %>% ungroup() %>% mutate(variable=str_replace_all(string=.$variable, pattern ="X", replacement="")) %>%
  merge(x=., y=legend, by.x="variable", by.y="id_genomes") %>% 
  merge(x=., y=filter(model.df,query.genome==6|target.genome==6) %>% mutate(variable=sapply(1:nrow(.),function(x){
    if(.$target.genome[x]==6){
      return(.$query.genome[x])
    }else{
      return(.$target.genome[x])
    }
  }))) %>% arrange(rsquared) %>% 
  mutate(name=str_trim(toupper(str_replace_all(string=.$name, pattern = "\\.fasta|_|complete_genome|chromosome", replacement=" ")))) %>%
  mutate(name=factor(.$name, levels=unique(.$name))) %>%
  ggplot(data=., mapping=aes(x=name, y=times, group=variable, fill=value, order=value)) +geom_bar(stat="identity") + 
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1),plot.margin = unit(c(1, 1, 1, 4), "cm")) +
  xlab("Genome") + ylab("Core/Outsider Ratio") + labs(title="Drosophila Core/Outsider Gene Classification") +
  geom_segment(mapping=aes(x=1, xend=nrow(legend)-1, y=0.1, yend=0.1), arrow=arrow(length=unit(0.5,"cm")),colour="red") +
  scale_fill_manual(values=c("#06960B","#999999","#F7728A"),labels=c("Core Genes","No Clear Separation","Outsider Genes"), guide=guide_legend(title=""))
plot(rsq.bar)
```

![](grant.report.end.2015_files/figure-html/dros distro-2.png) 

```r
png(file=report.2015("core.outsider.barplot.png"),width=1200, height=700)
plot(rsq.bar)
dev.off()
```

```
## png 
##   2
```
