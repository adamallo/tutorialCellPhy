# Tutorial CellPhy

### Collapsed boxes
You will see collapsed boxes along the tutorial document. Please, try to answer the questions on your own before clicking on them. Otherwise, this tutorial will be of little use for you!

<details><summary>Understood?</summary>
<p>
YES!
</p>
</details>

### Objectives
In this tutorial, we will use the "expert mode" of CellPhy (`cellphy.sh RAXML` command) to do some phylogenetic estimation of single-cell somatic data. We will get familiar with the phylogenetics pipeline, concepts of model selection and hypothesis testing, and with the graphical representation of phylogenies. Our toy experiment's objective will be to test the hypothesis that the different genetic markers of this dataset (SNVs) evolve at different evolutionary rates and generate a phylogeny of the samples with nodal support (bootstrap values).

### Data
We will use the example data distributed with the program and run the software from the `example` directory where sample input files reside. Thus, we will use the respective relative paths for all scripts (e.g., `../cellphy.sh` for the main script). This is a real single-cell SNV dataset thinned down to 500 SNVs.

### Is cellphy working?
Run cellphy's help to make sure the program is working for your computer. Read the output to get familiarized with the program's input options.
```bash
../cellphy.sh -h
```
### Runing cellphy using a genotype matrix (aligment) without modeling single-cell sequencing errors
View the contents of the input file CRC24l.ToySet.phy. For this, you can use a GUI text editor (e.g., atom) or the command line.

We will analize these data using the most general substitution model in CellPhy _GT16+F0_. This model assumes the data is known without error.  
$ ../cellphy.sh RAXML --msa CRC24.ToySet.phy --model GT16+FO --seed 2 --threads 1 --prefix model1

### Runing cellphy using a genotype matrix (alignment) modeling single-cell sequencing errors
Alternatively, we will analyze the data adding a single-cell sequencing error model _GT16+F0+E_.
$ ../cellphy.sh RAXML --msa CRC24.ToySet.phy --model GT16+FO+E --seed 2 --threads 1 --prefix model2

### Explore the outputs
Tree:
Parameters:


### Model selection
What model fits the data better?
_Tip_: look in the output for terms related to model selection, e.g., likelihood, number of parameters, AIC, BIC. Feel free to ask me questions about it, or use your favorite internet search engine to take a quick look at the topic.

<details><summary>Answer</summary>
<p>
Comparing the AIC or BIC of the models will give you a quick and clear answer that the _GT16+F0+E_ fits this data better (lower value). Alternatively, you could use a likelihood ratio test since these models are nested. I would not expect you to know this last bit and is not necessary at this point.
</p>
</details>

### Alternative data that incorporates variant calling uncertainty
Probabilistic variant callers not only output the final variants but also the uncertainty of each call. CellPhy can incorporate these uncertainties in the phylogenetic tree estimation instead of using its single-cell sequencing error model. In order to run CellPhy this way, you will need to provide it with a VCF file that contains genotype likelihood data.

Take a look at the provided input file CRC24.ToySet.vcf.

Run cellphy:

../cellphy.sh RAXML --msa CRC24.ToySet.vcf --model GT16+FO --seed 2 --threads 1 --prefix data2 

<details><summary>Can you compare this model with our previous best-fit?</summary>
<p>
No, you can't, since the data is different. We will continue with the other dataset. You do not need to finish this run.
</p>
</details>
 
### Hypothesis testing
Now that we know what model fits our data best, we will use  model to test the hypothesis that SNVs are evolving at different rates. For this, we will run an alternative model that incorporates site-heterogeneity called _GT16+F0+E+G_.

$ ../cellphy.sh RAXML --msa CRC24.ToySet.phy --model GT16+FO+E+G --seed 2 --threads 1 --prefix model3

<details><summary>Can we reject our null hypothesis that SNVs evolve at the same evolutionary rate?</summary>
<p>
Yes. You can use an LRT or calculate the relative likelihood by comparing AICs.
</p>
</details>


### Bootstrapping the data to estimate the tree support
Finally, we will perform a standard non-parametric bootstrap by re-sampling alignment columns and re-inferring trees for each bootstrap replicate. By default, CellPhy will automatically determine the optimal number of replicates (up to 1000), but here we will manually set the number of replicates to 100 in order to keep runtime reasonable.  

$ ../cellphy.sh RAXML --bootstrap --msa CRC24.ToySet.vcf --model GT16+FO+E+G --seed 2 --threads 1 --bs-trees 100 --prefix bootstrap

Now that we have the bootstrap trees, we can map the BS support values onto the (previously inferred) best-scoring ML tree:
../cellphy.sh RAXML --support -tree model3.raxml.bestTree --bs-trees bootstrap.raxml.bootstraps --prefix final --threads 1

CellPhy will output a tree with support values (in NEWICK format) that can be visualized using any tree viewer software. Open it in figtree.


<details><summary>Is the tree rooted?</summary>
<p>
No. Use the re-rooting option to root it using the healthy/normal tissue sample included.
</p>
</details>

<details><summary>Is the tree strongly supported?</summary>
<p>
No. Remember, we are using a thinned dataset to execute things fast.
</p>
</details>

For you convenience, CellPhy also provides a script called _support-map.R_ for plotting the support tree. Feel free to try to generate a pdf and svg figures.

<!--

# **Mapping mutations onto a phylogenetic tree**
<p style='text-align: justify;'>Cancer genomics studies are, for the most part, interested in understanding when "**driver**" mutations appeared in the malignant cell population. On this basis, we will next show how to map mutations onto a phylogenetic tree using CellPhy. Although the full VCF can be used, users are free to choose which mutations they wish to map onto the inferred phylogeny. In this tutorial, we will focus solely on a set of 15 exonic mutations from the original VCF.   
  
Our input files will therefore consist of a "trimmed" VCF only carrying this subset of exonic mutations, together with the best tree and model estimates from our original CellPhy run (tree search).</p>  

<br/>  

\scriptsize

```{r, engine = 'bash', eval = FALSE}
$ head -n 5 CRC24.MutationsMap 
#Chr	Position	GeneID
2	71042907	CLEC4F
2	142274377	LRP1B
3	33048242	GLB1
4	16764214	LDB2
$ bcftools view -T CRC24.MutationsMap CRC24.ToySet.vcf -O v -o CRC24.ToySet.Exonic.vcf
$ ../cellphy.sh RAXML --mutmap \
    --msa CRC24.ToySet.Exonic.vcf \
    --model CRC24.VCF.GL16.Tree.raxml.bestModel \ 
    --tree CRC24.VCF.GL16.Tree.raxml.bestTree \
    --opt-branches off --prefix CRC24.ToySet.Exonic.Mapped --threads 1
Branch-labeled tree saved to: CRC24.ToySet.Exonic.Mapped.raxml.mutationMapTree
Per-branch mutation list saved to: CRC24.ToySet.Exonic.Mapped.raxml.mutationMapList
```

\normalsize

***

<br/>  


# **Visualizing the results**
<p style='text-align: justify;'>Once it's done, CellPhy will output 2 distinct files:.</p>

* **(A)** _CRC24.ToySet.Exonic.Mapped.raxml.mutationMapTree_  
   &rarr; Newick tree file with indexed branches
* **(B)** _CRC24.ToySet.Exonic.Mapped.raxml.mutationMapList_  
   &rarr; Text file with the number and the list of mutations per branch


<p style='text-align: justify;'>We can now use the _mutation-map.R_ accompanying script to plot the mutations onto the phylogenetic tree. If you run this script wihtout any parameters, it will show a help message:</p>

<br/> 

\scriptsize

```{r, engine = 'bash', eval = FALSE}
$ ../script/mutation-map.R 
CellPhy - Mutation mapping plot - 22.07.2020
Created by: Alexey Kozlov, Joao M Alves, Alexandros Stamatakis & David Posada
Usage: ./mutation-map.R raxml.mutationMapTree raxml.mutationMapList Outgroup Output_prefix [geneIDs]
*Required files:
	-Tree
	-Mutation List
	-Outgroup name (comma-delimited list of taxa or NONE)
	-Output Prefix
*Optional:
	-Gene IDs (tab-delimited)
```

\normalsize

<p style='text-align: justify;'>Now let's run it again but this time with the required parameters:</p>

<br/>  

\scriptsize

```{r, engine = 'bash', eval = FALSE}
$ ../script/mutation-map.R \
    CRC24.ToySet.Exonic.Mapped.raxml.mutationMapTree \
    CRC24.ToySet.Exonic.Mapped.raxml.mutationMapList \
    Healthy \
    CRC24.ToySet.ExonicMutMap
Generating mutation-mapped tree plot...
Done!
```

\normalsize

<p style='text-align: justify;'>If everything went as expected, you should have generated the following figure, in PDF format (_CRC24.ToySet.ExonicMutMap.pdf_), where the mutations are mapped onto the tree branches (**Figure 2**).</p>

![CRC24 phylogenetic tree with 15 exonic mutations mapped (genomic position)](../example/CRC24.ToySet.ExonicMutMap.pdf){ width=100% }

<p style='text-align: justify;'>If you are interested in plotting the gene names instead, you can provide a tab-delimited file (as the one we used to subset our original VCF) linking the genomic position to its gene ID:</p>

<br/>  

\scriptsize

```{r, engine = 'bash', eval = FALSE}
$ head -n 5 CRC24.MutationsMap 
#Chr	Position	GeneID
2	71042907	CLEC4F
2	142274377	LRP1B
3	33048242	GLB1
4	16764214	LDB2
```

\normalsize

<p style='text-align: justify;'>Afterwards, we can run _mutation-map.R_ again, but changing the output prefix so that you don't overwrite the previous results:</p>

<br/>  

\scriptsize
```{r, engine = 'bash', eval = FALSE}
$ ../script/mutation-map.R \
    CRC24.ToySet.Exonic.Mapped.raxml.mutationMapTree \
    CRC24.ToySet.Exonic.Mapped.raxml.mutationMapList \
    Healthy \
    CRC24.ToySet.ExonicMutMap-GeneID \
    CRC24.MutationsMap
Converting positions to GeneID...
Generating mutation-mapped tree plot...
Done!
```
\normalsize

<p style='text-align: justify;'>You will notice that our tree now has the gene names displayed, instead of the genomic positions (Figure 3).</p>

<center>

![CRC24 phylogenetic tree with 15 exonic mutations mapped (gene names)](../example/CRC24.ToySet.ExonicMutMap-GeneID.pdf){ width=100% }

</center>-->
