10/1/2021

# Data Wrangling and Integration with dplyr, tidyr, using genome annotations and gene ontology definitions

We often have to integrate and reformat data to integrate into a report
or chart. In my case, I had to provide gene annotation details that are
at deeper level (transcript level) into a plot that was at a higher
level (locus level). Additionally, a column in the annotation file that
had gene ontology terms, was comma separated with no definitions. Those
terms needed to be placed into their own columns with definitions
provided from another file.

Here is a sample of what each files looked like:

**Genome Annotation**:

|  X.pacId | locusName       | transcriptName    | peptideName         | Pfam                            | Panther                   | KOG     | KEGG.ec  | KO     | GO                                                       | Best.hit.arabi.name | arabi.symbol                                                 | arabi.defline                                                              | Best.hit.rice.name | rice.symbol | rice.defline                                                  |
| -------: | :-------------- | :---------------- | :------------------ | :------------------------------ | :------------------------ | :------ | :------- | :----- | :------------------------------------------------------- | :------------------ | :----------------------------------------------------------- | :------------------------------------------------------------------------- | :----------------- | :---------- | :------------------------------------------------------------ |
| 32806936 | Bradi0135s00100 | Bradi0135s00100.1 | Bradi0135s00100.1.p | PF00295                         | PTHR31375,PTHR31375:SF37  |         | 3.2.1.15 | K01213 | <GO:0005975,GO:0004650>                                  | AT3G07820.1         |                                                              | Pectin lyase-like superfamily protein                                      | LOC\_Os02g10300.1  | NA          | polygalacturonase, putative, expressed                        |
| 32805405 | Bradi1g00215    | Bradi1g00215.1    | Bradi1g00215.1.p    | PF00076                         | PTHR10501,PTHR10501:SF24  |         |          |        | <GO:0003676,GO:0017069,GO:0000398>                       | AT1G21320.1         |                                                              | nucleotide binding;nucleic acid binding                                    | LOC\_Os08g43360.1  | NA          | RNA recognition motif containing protein, putative, expressed |
| 32805623 | Bradi1g00272    | Bradi1g00272.1    | Bradi1g00272.1.p    | PF01554                         | PTHR11206,PTHR11206:SF102 | KOG1347 |          |        | <GO:0055085,GO:0016020,GO:0015297,GO:0015238,GO:0006855> | AT1G33110.1         |                                                              | MATE efflux family protein                                                 | LOC\_Os12g03260.1  | NA          | MATE efflux family protein, putative, expressed               |
| 32799309 | Bradi1g00350    | Bradi1g00350.1    | Bradi1g00350.1.p    | PF05000,PF04998,PF04992         | PTHR19376,PTHR19376:SF37  |         | 2.7.7.6  |        | <GO:0006351,GO:0003899,GO:0003677>                       | AT4G35800.1         | NRPB1,RNA\_POL\_II\_LS,RNA\_POL\_II\_LSRNA\_POL\_II\_LS,RPB1 | RNA polymerase II large subunit                                            | LOC\_Os05g05860.1  | NA          | retrotransposon protein, putative, unclassified, expressed    |
| 32793603 | Bradi1g00400    | Bradi1g00400.2    | Bradi1g00400.2.p    | PF02485                         | PTHR31042,PTHR31042:SF25  |         |          |        | <GO:0016020,GO:0008375>                                  | AT1G10280.1         |                                                              | Core-2/I-branching beta-1,6-N-acetylglucosaminyltransferase family protein | LOC\_Os04g20420.1  | NA          | DNA binding protein, putative, expressed                      |
| 32804480 | Bradi1g00607    | Bradi1g00607.1    | Bradi1g00607.1.p    | PF02736,PF00612,PF01843,PF00063 | PTHR13140,PTHR13140:SF382 |         | 3.6.4.1  | K10357 | <GO:0016459,GO:0005524,GO:0003774,GO:0005515>            | AT1G04160.1         | ATXIB,XI-8,XI-B,XIB                                          | myosin XI B                                                                | LOC\_Os03g64290.1  | NA          | myosin, putative, expressed                                   |
| 32804481 | Bradi1g00607    | Bradi1g00607.2    | Bradi1g00607.2.p    | PF00612,PF01843,PF00063         | PTHR13140,PTHR13140:SF382 |         | 3.6.4.1  |        | <GO:0005515,GO:0016459,GO:0005524,GO:0003774>            | AT1G04160.1         | ATXIB,XI-8,XI-B,XIB                                          | myosin XI B                                                                | LOC\_Os03g64290.1  | NA          | myosin, putative, expressed                                   |

**Gene ontology**:

| ID           | Name                                                     | Namespace           | alt\_id                 | Def                                                                                                                                                                                                                                                                                                                                     |
| :----------- | :------------------------------------------------------- | :------------------ | :---------------------- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <GO:0000001> | mitochondrion inheritance                                | biological\_process |                         | The distribution of mitochondria, including the mitochondrial genome, into daughter cells after mitosis or meiosis, mediated by interactions between mitochondria and the cytoskeleton. \[GOC:mcc, <PMID:10873824>, <PMID:11389764>\]                                                                                                   |
| <GO:0000002> | mitochondrial genome maintenance                         | biological\_process |                         | The maintenance of the structure and integrity of the mitochondrial genome; includes replication and segregation of the mitochondrial chromosome. \[GOC:ai, GOC:vw\]                                                                                                                                                                    |
| <GO:0000003> | reproduction                                             | biological\_process | <GO:0019952,GO:0050876> |                                                                                                                                                                                                                                                                                                                                         |
| <GO:0000005> | obsolete ribosomal chaperone activity                    | molecular\_function |                         | OBSOLETE. Assists in the correct assembly of ribosomes or ribosomal subunits in vivo, but is not a component of the assembled ribosome when performing its normal biological function. \[GOC:jl, <PMID:12150913>\]                                                                                                                      |
| <GO:0000006> | high-affinity zinc transmembrane transporter activity    | molecular\_function |                         | Enables the transfer of zinc ions (Zn2+) from one side of a membrane to the other, probably powered by proton motive force. In high-affinity transport the transporter is able to bind the solute even if it is only present at very low concentrations. \[TC:2.A.5.1.1\]                                                               |
| <GO:0000007> | low-affinity zinc ion transmembrane transporter activity | molecular\_function |                         | Enables the transfer of a solute or solutes from one side of a membrane to the other according to the reaction: Zn2+ = Zn2+, probably powered by proton motive force. In low-affinity transport the transporter is able to bind the solute only if it is present at very high concentrations. \[GOC:mtg\_transport, <ISBN:0815340729>\] |

Its important to note that the gene ontology file has already been
modified from its intial format in the [Converting non-tabular data into
tabular data using
Python](https://github.com/patmendoza330/geneontologyconversion)
repository. Feel free to look that over for ways in which I modified
that format into tabular format.

## Issues

So, there were several issues, I needed to:

1.  Collapse the data in the genome annotation file so that only one
    record appeared for each locus but also included a column with all
    unique gene ontology terms for all transcripts that were associated
    with the gene.
2.  Split the comma delimited gene ontology column into as many columns
    were necessary and include the definition for those terms from the
    gene ontology file.

## Running files in respository

All of the script is held in the
[datawrangling.rmd](https://github.com/patmendoza330/annotationwrangling/blob/main/datawrangling.rmd)
r markdown file and files in the
[supporting.files](https://github.com/patmendoza330/annotationwrangling/tree/main/supporting.files)
folder. I would recomend downloading those files and running the
markdown file in RStudio.

## Solution

### Libraries needed

Lets first install all of the needed libraries:

``` r
library(dplyr)
library(tidyr)
```

Briefly, tidyr will allow us to reformat existing data into a format
that we’d like. dplyr will allow us to summarize, join different tables,
and create new fields. dplyr lets us use traditional SQL scripts as well
as some amazing new stuff.

### Gene annotation file

Lets take the genome annotation file first and go over the steps that we
need.

Firs things first, lets load in the table:

``` r
y1 <- read.delim("supporting.files/Bdistachyon_314_v3.1.annotation_info.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
knitr::kable(head(y1), "pipe")
```

|  X.pacId | locusName       | transcriptName    | peptideName         | Pfam                    | Panther                   | KOG | KEGG.ec | KO     | GO           | Best.hit.arabi.name | arabi.symbol | arabi.defline                                   | Best.hit.rice.name | rice.symbol | rice.defline                              |
| -------: | :-------------- | :---------------- | :------------------ | :---------------------- | :------------------------ | :-- | :------ | :----- | :----------- | :------------------ | :----------- | :---------------------------------------------- | :----------------- | :---------- | :---------------------------------------- |
| 32823775 | Bradi0012s00100 | Bradi0012s00100.1 | Bradi0012s00100.1.p | PF03144,PF03143,PF00009 | PTHR23115,PTHR23115:SF157 |     | 3.6.5.3 | K03231 | <GO:0005525> | AT1G07920.1         |              | GTP binding Elongation factor Tu family protein | LOC\_Os03g08010.1  | NA          | elongation factor Tu, putative, expressed |
| 32823776 | Bradi0012s00201 | Bradi0012s00201.1 | Bradi0012s00201.1.p |                         |                           |     |         |        |              |                     |              |                                                 |                    | NA          |                                           |
| 32823777 | Bradi0012s00201 | Bradi0012s00201.2 | Bradi0012s00201.2.p |                         |                           |     |         |        |              |                     |              |                                                 |                    | NA          |                                           |
| 32823779 | Bradi0012s00201 | Bradi0012s00201.3 | Bradi0012s00201.3.p |                         |                           |     |         |        |              |                     |              |                                                 |                    | NA          |                                           |
| 32823780 | Bradi0012s00201 | Bradi0012s00201.4 | Bradi0012s00201.4.p |                         |                           |     |         |        |              |                     |              |                                                 |                    | NA          |                                           |
| 32823778 | Bradi0012s00201 | Bradi0012s00201.5 | Bradi0012s00201.5.p |                         |                           |     |         |        |              |                     |              |                                                 |                    | NA          |                                           |

Next, lets select only the fields that we want:

``` r
y <- y1 %>% select(locusName, GO)
knitr::kable(head(y), "pipe")
```

| locusName       | GO           |
| :-------------- | :----------- |
| Bradi0012s00100 | <GO:0005525> |
| Bradi0012s00201 |              |
| Bradi0012s00201 |              |
| Bradi0012s00201 |              |
| Bradi0012s00201 |              |
| Bradi0012s00201 |              |

``` r
nrow(y)
```

    ## [1] 52972

Ok, so now we have a table that has 52,972 records. We need to collapse
the locusName column and ensure that we have only unique gene ontology
terms:

``` r
y <- y %>%
   group_by(locusName) %>%
   summarise(GO = paste(GO, collapse = ","))
knitr::kable(head(y), "pipe")
```

| locusName       | GO                      |
| :-------------- | :---------------------- |
| Bradi0012s00100 | <GO:0005525>            |
| Bradi0012s00201 | ,,,,                    |
| Bradi0014s00100 | <GO:0043531>            |
| Bradi0135s00100 | <GO:0005975,GO:0004650> |
| Bradi0180s00100 | <GO:0005975,GO:0004650> |
| Bradi1g00200    | <GO:0005515>            |

``` r
nrow(y)
```

    ## [1] 34310

``` r
length(unique(y$locusName))
```

    ## [1] 34310

As you can see, we now have collapsed the table into 34,310 rows.
Additionally, the number of unique locusName is the same 34,310, so
we’re good\! Unfortunately, we have some rows with a bunch of commas,
and also some locusName’s with duplicate GO terms. Lets fix the
duplicate GO terms first.

``` r
y[which(y$locusName=="Bradi1g00517"), ]
```

    ## # A tibble: 1 x 2
    ##   locusName    GO                                                               
    ##   <chr>        <chr>                                                            
    ## 1 Bradi1g00517 GO:0003676,GO:0003676,GO:0003676,GO:0003676,GO:0003676,GO:000367~

This is an example of one of the items that we need to fix. When we
collapsed the table, it pasted each GO term for all of the transcrripts,
but many of the transcripts had the same GO term. We can fix this by
using the sapply function below:

``` r
y$GO <- sapply(strsplit(y$GO, ","), function(x){
  paste(unique(trimws(x)), collapse = ",")
})
y[which(y$locusName=="Bradi1g00517"), ]
```

    ## # A tibble: 1 x 2
    ##   locusName    GO        
    ##   <chr>        <chr>     
    ## 1 Bradi1g00517 GO:0003676

This function takes the GO column and replaces it with unique entries
only. It does this by splitting the values in each row by the comma
delimiter, pulling only unique entries from the revised list, then
collapsing them back into a comma delimited string.

You can see that the Bradi1g00517 locusName now has only one GO term
associated which is exactly what we want.

Here is what the head of the table looks like now:

``` r
knitr::kable(head(y), "pipe")
```

| locusName       | GO                      |
| :-------------- | :---------------------- |
| Bradi0012s00100 | <GO:0005525>            |
| Bradi0012s00201 |                         |
| Bradi0014s00100 | <GO:0043531>            |
| Bradi0135s00100 | <GO:0005975,GO:0004650> |
| Bradi0180s00100 | <GO:0005975,GO:0004650> |
| Bradi1g00200    | <GO:0005515>            |

## Gene Ontology File

Lets load in the file

``` r
GOFile <- read.delim("supporting.files/go.obo.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
knitr::kable(head(GOFile), "pipe")
```

| ID           | Name                                                     | Namespace           | alt\_id                 | Def                                                                                                                                                                                                                                                                                                                                     |
| :----------- | :------------------------------------------------------- | :------------------ | :---------------------- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <GO:0000001> | mitochondrion inheritance                                | biological\_process |                         | The distribution of mitochondria, including the mitochondrial genome, into daughter cells after mitosis or meiosis, mediated by interactions between mitochondria and the cytoskeleton. \[GOC:mcc, <PMID:10873824>, <PMID:11389764>\]                                                                                                   |
| <GO:0000002> | mitochondrial genome maintenance                         | biological\_process |                         | The maintenance of the structure and integrity of the mitochondrial genome; includes replication and segregation of the mitochondrial chromosome. \[GOC:ai, GOC:vw\]                                                                                                                                                                    |
| <GO:0000003> | reproduction                                             | biological\_process | <GO:0019952,GO:0050876> |                                                                                                                                                                                                                                                                                                                                         |
| <GO:0000005> | obsolete ribosomal chaperone activity                    | molecular\_function |                         | OBSOLETE. Assists in the correct assembly of ribosomes or ribosomal subunits in vivo, but is not a component of the assembled ribosome when performing its normal biological function. \[GOC:jl, <PMID:12150913>\]                                                                                                                      |
| <GO:0000006> | high-affinity zinc transmembrane transporter activity    | molecular\_function |                         | Enables the transfer of zinc ions (Zn2+) from one side of a membrane to the other, probably powered by proton motive force. In high-affinity transport the transporter is able to bind the solute even if it is only present at very low concentrations. \[TC:2.A.5.1.1\]                                                               |
| <GO:0000007> | low-affinity zinc ion transmembrane transporter activity | molecular\_function |                         | Enables the transfer of a solute or solutes from one side of a membrane to the other according to the reaction: Zn2+ = Zn2+, probably powered by proton motive force. In low-affinity transport the transporter is able to bind the solute only if it is present at very high concentrations. \[GOC:mtg\_transport, <ISBN:0815340729>\] |

I also want to rename the first column so that it matches our annotation
file:

``` r
names(GOFile)[1] <- 'GO'
knitr::kable(head(GOFile), "pipe")
```

| GO           | Name                                                     | Namespace           | alt\_id                 | Def                                                                                                                                                                                                                                                                                                                                     |
| :----------- | :------------------------------------------------------- | :------------------ | :---------------------- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| <GO:0000001> | mitochondrion inheritance                                | biological\_process |                         | The distribution of mitochondria, including the mitochondrial genome, into daughter cells after mitosis or meiosis, mediated by interactions between mitochondria and the cytoskeleton. \[GOC:mcc, <PMID:10873824>, <PMID:11389764>\]                                                                                                   |
| <GO:0000002> | mitochondrial genome maintenance                         | biological\_process |                         | The maintenance of the structure and integrity of the mitochondrial genome; includes replication and segregation of the mitochondrial chromosome. \[GOC:ai, GOC:vw\]                                                                                                                                                                    |
| <GO:0000003> | reproduction                                             | biological\_process | <GO:0019952,GO:0050876> |                                                                                                                                                                                                                                                                                                                                         |
| <GO:0000005> | obsolete ribosomal chaperone activity                    | molecular\_function |                         | OBSOLETE. Assists in the correct assembly of ribosomes or ribosomal subunits in vivo, but is not a component of the assembled ribosome when performing its normal biological function. \[GOC:jl, <PMID:12150913>\]                                                                                                                      |
| <GO:0000006> | high-affinity zinc transmembrane transporter activity    | molecular\_function |                         | Enables the transfer of zinc ions (Zn2+) from one side of a membrane to the other, probably powered by proton motive force. In high-affinity transport the transporter is able to bind the solute even if it is only present at very low concentrations. \[TC:2.A.5.1.1\]                                                               |
| <GO:0000007> | low-affinity zinc ion transmembrane transporter activity | molecular\_function |                         | Enables the transfer of a solute or solutes from one side of a membrane to the other according to the reaction: Zn2+ = Zn2+, probably powered by proton motive force. In low-affinity transport the transporter is able to bind the solute only if it is present at very high concentrations. \[GOC:mtg\_transport, <ISBN:0815340729>\] |

## Bringing it all together

Ok, so now we’ve made some modifications to the annotation file and
we’ve loaded in the GO terms and their definitions. Now we need a way
to join them together and present a table that has each gene, and all of
their associated GO terms and definitions as respective columns.

First we need to flatten our annotation file so that we have multiple
rows for every GO term:

``` r
gene.GO <- tibble::as_tibble(y)
gene.GO <- gene.GO %>%
  dplyr::select(locusName, GO) %>%
  filter(GO != "") %>%
  arrange(locusName) %>%
  separate_rows(GO, sep = "[,]")
knitr::kable(head(gene.GO), "pipe")
```

| locusName       | GO           |
| :-------------- | :----------- |
| Bradi0012s00100 | <GO:0005525> |
| Bradi0014s00100 | <GO:0043531> |
| Bradi0135s00100 | <GO:0005975> |
| Bradi0135s00100 | <GO:0004650> |
| Bradi0180s00100 | <GO:0005975> |
| Bradi0180s00100 | <GO:0004650> |

As you can see, the separate\_rows function will create duplicate rows
for every GO term in the GO column that is separated by a comma. Now we
can do a join to the gene ontology file:

``` r
gene.GO <- gene.GO %>%
  left_join(GOFile) %>%
  dplyr::select(locusName, GO, Namespace, Name, Def) %>%
  arrange(locusName, GO)
knitr::kable(head(gene.GO), "pipe")
```

| locusName       | GO           | Namespace           | Name                           | Def                                                                                                                                 |
| :-------------- | :----------- | :------------------ | :----------------------------- | :---------------------------------------------------------------------------------------------------------------------------------- |
| Bradi0012s00100 | <GO:0005525> | molecular\_function | GTP binding                    | Interacting selectively and non-covalently with GTP, guanosine triphosphate. \[GOC:ai\]                                             |
| Bradi0014s00100 | <GO:0043531> | molecular\_function | ADP binding                    | Interacting selectively and non-covalently with ADP, adenosine 5’-diphosphate. \[GOC:jl\]                                           |
| Bradi0135s00100 | <GO:0004650> | molecular\_function | polygalacturonase activity     | Catalysis of the random hydrolysis of (1-\>4)-alpha-D-galactosiduronic linkages in pectate and other galacturonans. \[EC:3.2.1.15\] |
| Bradi0135s00100 | <GO:0005975> | biological\_process | carbohydrate metabolic process |                                                                                                                                     |
| Bradi0180s00100 | <GO:0004650> | molecular\_function | polygalacturonase activity     | Catalysis of the random hydrolysis of (1-\>4)-alpha-D-galactosiduronic linkages in pectate and other galacturonans. \[EC:3.2.1.15\] |
| Bradi0180s00100 | <GO:0005975> | biological\_process | carbohydrate metabolic process |                                                                                                                                     |

Now, this is great, and its getting closer to where we want to be, but
we need to have a table with a single row for each locusName and a
column that has the GO term, Namespace, Name, and Def all pasted into an
individual column for each GO term.

Lets start by collapsing the columns in this table by pasting all the
gene ontology related information into one column (separated by a
comma):

``` r
gene.GO$GO <- paste(gene.GO$GO, gene.GO$Namespace, 
                    gene.GO$Name, gene.GO$Def, sep = ", ")
gene.GO <- gene.GO[, 1:2] 
knitr::kable(head(gene.GO), "pipe")
```

| locusName       | GO                                                                                                                                                                                                 |
| :-------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Bradi0012s00100 | <GO:0005525>, molecular\_function, GTP binding, Interacting selectively and non-covalently with GTP, guanosine triphosphate. \[GOC:ai\]                                                            |
| Bradi0014s00100 | <GO:0043531>, molecular\_function, ADP binding, Interacting selectively and non-covalently with ADP, adenosine 5’-diphosphate. \[GOC:jl\]                                                          |
| Bradi0135s00100 | <GO:0004650>, molecular\_function, polygalacturonase activity, Catalysis of the random hydrolysis of (1-\>4)-alpha-D-galactosiduronic linkages in pectate and other galacturonans. \[EC:3.2.1.15\] |
| Bradi0135s00100 | <GO:0005975>, biological\_process, carbohydrate metabolic process,                                                                                                                                 |
| Bradi0180s00100 | <GO:0004650>, molecular\_function, polygalacturonase activity, Catalysis of the random hydrolysis of (1-\>4)-alpha-D-galactosiduronic linkages in pectate and other galacturonans. \[EC:3.2.1.15\] |
| Bradi0180s00100 | <GO:0005975>, biological\_process, carbohydrate metabolic process,                                                                                                                                 |

Now we need to finish up the table by adding an intermediate column that
will label all GO terms incrementally, then using the pivot\_wider
function to put them into columns

``` r
gene.GO <- gene.GO %>% 
  group_by(locusName) %>%
  mutate(Var = paste("GO", 1:n(), sep = "")) %>% 
  arrange(locusName, Var, GO)
knitr::kable(head(gene.GO), "pipe")
```

| locusName       | GO                                                                                                                                                                                                 | Var |
| :-------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :-- |
| Bradi0012s00100 | <GO:0005525>, molecular\_function, GTP binding, Interacting selectively and non-covalently with GTP, guanosine triphosphate. \[GOC:ai\]                                                            | GO1 |
| Bradi0014s00100 | <GO:0043531>, molecular\_function, ADP binding, Interacting selectively and non-covalently with ADP, adenosine 5’-diphosphate. \[GOC:jl\]                                                          | GO1 |
| Bradi0135s00100 | <GO:0004650>, molecular\_function, polygalacturonase activity, Catalysis of the random hydrolysis of (1-\>4)-alpha-D-galactosiduronic linkages in pectate and other galacturonans. \[EC:3.2.1.15\] | GO1 |
| Bradi0135s00100 | <GO:0005975>, biological\_process, carbohydrate metabolic process,                                                                                                                                 | GO2 |
| Bradi0180s00100 | <GO:0004650>, molecular\_function, polygalacturonase activity, Catalysis of the random hydrolysis of (1-\>4)-alpha-D-galactosiduronic linkages in pectate and other galacturonans. \[EC:3.2.1.15\] | GO1 |
| Bradi0180s00100 | <GO:0005975>, biological\_process, carbohydrate metabolic process,                                                                                                                                 | GO2 |

The way that we’ve constructed the Var column creates a column that
stores the number of gene ontology terms for each locusName. When a new
GO term is introduced in the table, the Var column has an incremented
value.

The pivot\_wider function (below) will give us our final desired table.
If a locusName does not have a GO term for a specific column, NA is
entered.

``` r
gene.GO <- gene.GO %>%
  pivot_wider(names_from = Var, values_from = GO)
knitr::kable(head(gene.GO), "pipe")
```

| locusName       | GO1                                                                                                                                                                                                | GO2                                                                | GO3 | GO4 | GO5 | GO6 | GO7 | GO8 | GO9 |
| :-------------- | :------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :----------------------------------------------------------------- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| Bradi0012s00100 | <GO:0005525>, molecular\_function, GTP binding, Interacting selectively and non-covalently with GTP, guanosine triphosphate. \[GOC:ai\]                                                            | NA                                                                 | NA  | NA  | NA  | NA  | NA  | NA  | NA  |
| Bradi0014s00100 | <GO:0043531>, molecular\_function, ADP binding, Interacting selectively and non-covalently with ADP, adenosine 5’-diphosphate. \[GOC:jl\]                                                          | NA                                                                 | NA  | NA  | NA  | NA  | NA  | NA  | NA  |
| Bradi0135s00100 | <GO:0004650>, molecular\_function, polygalacturonase activity, Catalysis of the random hydrolysis of (1-\>4)-alpha-D-galactosiduronic linkages in pectate and other galacturonans. \[EC:3.2.1.15\] | <GO:0005975>, biological\_process, carbohydrate metabolic process, | NA  | NA  | NA  | NA  | NA  | NA  | NA  |
| Bradi0180s00100 | <GO:0004650>, molecular\_function, polygalacturonase activity, Catalysis of the random hydrolysis of (1-\>4)-alpha-D-galactosiduronic linkages in pectate and other galacturonans. \[EC:3.2.1.15\] | <GO:0005975>, biological\_process, carbohydrate metabolic process, | NA  | NA  | NA  | NA  | NA  | NA  | NA  |
| Bradi1g00200    | <GO:0005515>, molecular\_function, protein binding,                                                                                                                                                | NA                                                                 | NA  | NA  | NA  | NA  | NA  | NA  | NA  |
| Bradi1g00210    | <GO:0005515>, molecular\_function, protein binding,                                                                                                                                                | NA                                                                 | NA  | NA  | NA  | NA  | NA  | NA  | NA  |

Now, if we wanted to add this back to our original annotation file, we’d
simply make a join between this new table and our original:

``` r
y1 <- y1 %>%
  left_join(gene.GO, by="locusName")
knitr::kable(head(y1), "pipe")
```

|  X.pacId | locusName       | transcriptName    | peptideName         | Pfam                    | Panther                   | KOG | KEGG.ec | KO     | GO           | Best.hit.arabi.name | arabi.symbol | arabi.defline                                   | Best.hit.rice.name | rice.symbol | rice.defline                              | GO1                                                                                                                                     | GO2 | GO3 | GO4 | GO5 | GO6 | GO7 | GO8 | GO9 |
| -------: | :-------------- | :---------------- | :------------------ | :---------------------- | :------------------------ | :-- | :------ | :----- | :----------- | :------------------ | :----------- | :---------------------------------------------- | :----------------- | :---------- | :---------------------------------------- | :-------------------------------------------------------------------------------------------------------------------------------------- | :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| 32823775 | Bradi0012s00100 | Bradi0012s00100.1 | Bradi0012s00100.1.p | PF03144,PF03143,PF00009 | PTHR23115,PTHR23115:SF157 |     | 3.6.5.3 | K03231 | <GO:0005525> | AT1G07920.1         |              | GTP binding Elongation factor Tu family protein | LOC\_Os03g08010.1  | NA          | elongation factor Tu, putative, expressed | <GO:0005525>, molecular\_function, GTP binding, Interacting selectively and non-covalently with GTP, guanosine triphosphate. \[GOC:ai\] | NA  | NA  | NA  | NA  | NA  | NA  | NA  | NA  |
| 32823776 | Bradi0012s00201 | Bradi0012s00201.1 | Bradi0012s00201.1.p |                         |                           |     |         |        |              |                     |              |                                                 |                    | NA          |                                           | NA                                                                                                                                      | NA  | NA  | NA  | NA  | NA  | NA  | NA  | NA  |
| 32823777 | Bradi0012s00201 | Bradi0012s00201.2 | Bradi0012s00201.2.p |                         |                           |     |         |        |              |                     |              |                                                 |                    | NA          |                                           | NA                                                                                                                                      | NA  | NA  | NA  | NA  | NA  | NA  | NA  | NA  |
| 32823779 | Bradi0012s00201 | Bradi0012s00201.3 | Bradi0012s00201.3.p |                         |                           |     |         |        |              |                     |              |                                                 |                    | NA          |                                           | NA                                                                                                                                      | NA  | NA  | NA  | NA  | NA  | NA  | NA  | NA  |
| 32823780 | Bradi0012s00201 | Bradi0012s00201.4 | Bradi0012s00201.4.p |                         |                           |     |         |        |              |                     |              |                                                 |                    | NA          |                                           | NA                                                                                                                                      | NA  | NA  | NA  | NA  | NA  | NA  | NA  | NA  |
| 32823778 | Bradi0012s00201 | Bradi0012s00201.5 | Bradi0012s00201.5.p |                         |                           |     |         |        |              |                     |              |                                                 |                    | NA          |                                           | NA                                                                                                                                      | NA  | NA  | NA  | NA  | NA  | NA  | NA  | NA  |

# Conclusion

We now have a an annotation table that includes not just the GO terms
without definitions, but all deifinitions for each GO term included in
their own column.

# Citations

Ashburner, et al. (2000, 2000/05/01). Gene Ontology: tool for the
unification of biology. Nature Genetics, 25(1), 25-29.
<https://doi.org/10.1038/75556>

Gene Ontology, C. (2021). The Gene Ontology resource: enriching a GOld
mine. Nucleic acids research, 49(D1), D325-D334.
<https://doi.org/10.1093/nar/gkaa1113>

The International Brachypodium, I. (2010, 02/11/online). Genome
sequencing and analysis of the model grass Brachypodium distachyon
\[Article\]. Nature, 463, 763. <https://doi.org/10.1038/nature08747>
