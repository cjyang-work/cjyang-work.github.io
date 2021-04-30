---
layout: post
title: magicdesign
---

## Introduction
Welcome to :star2: `magicdesign` :star2:!

MAGIC (Multiparental Advanced Generation Inter-Cross) is a highly recombined population of multiple founders. MAGIC population is a versatile genetic resource for quantitative trait locus (QTL) mapping, fine-mapping, GxE dissection, genomic prediction, breeding and many more. Please check [this](https://doi.org/10.1038/s41437-020-0336-6) out for an in-depth review of MAGIC.

`magicdesign` is an R package devoted to creating and testing MAGIC population designs via simulations. Click [here](https://cjyang-sruc.github.io/magicdesign_vignette) for detailed instructions on installation and usage.

<iframe width="560" height="315"
src="https://www.youtube-nocookie.com/embed/94JgxeRFSxc"
title="MAGIC population design"
frameborder="0"
allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture"
allowfullscreen></iframe>

## Frequently asked questions (FAQs)
<details>
  <summary>1. How to use <code>magicdesign</code>?</summary>
  <p>First, make sure you have <code>devtools</code> installed in <code>R</code>. Skip this if you already have <code>devtools</code>.</p>
  <p><code>install.packages("devtools")</code></p>
  <p>Next, install <code>magicdesign</code>.</p>
  <p><code>devtools::install_github("cjyang-sruc/magicdesign")</code></p>
  <p>If the installation is successful, run the following:</p>
  <p><code>library(magicdesign)</code></p>
  <p>Detailed user instructions can be found <a href="https://cjyang-sruc.github.io/files/magicdesign_vignette">here</a>.</p>
  <br>
</details>

<details>
  <summary>2. Which functions in <code>magicdesign</code> do I need?</summary>
  <ul>
    <li><code>magic.eval</code>: create a MAGIC population design and simulate the population.</li>
    <li><code>magic.summary</code>: tabulate the information of all designs.</li>
    <li><code>magic.plot</code>: plot the distributions of recombinant haplotypes and founder genomes from all designs.</li>
    <li><code>magic.ped2plot</code>: plot the pedigree of a MAGIC population design.</li>
  </ul>
  <br>
</details>

<details>
  <summary>3. I have an issue!</summary>
  <p>Please report the issue with error codes to <a href="mailto=:cyang@sruc.ac.uk">cyang@sruc.ac.uk</a>.</p>
  <br>
</details>

<details>
  <summary>4. How do I create a MAGIC population with full design?</summary>
  <code>mpop <- magic.eval(m=8, m=45, reps=c(1,1,2), self=c(0,0,4), balanced=T)</code>
  <br>
</details>

<details>
  <summary>5. How do I create a MAGIC population with partial balanced design?</summary>
  <code>mpop <- magic.eval(m=8, m=1, reps=c(1,1,10), self=c(0,0,4), balanced=T)</code>
  <br>
</details>

<details>
  <summary>6. How do I create a MAGIC population with partial unbalanced design?</summary>
  <code>mpop <- magic.eval(m=8, m=7, reps=c(1,1,10), self=c(0,0,4), balanced=F)</code>
  <br>
</details>

<details>
  <summary>7. How do I create a MAGIC population with basic design?</summary>
  <code>mpop <- magic.eval(m=8, m=0, reps=c(1,4,20), self=c(0,0,4))</code>
  <br>
</details>

<details>
  <summary>8. How do I create a MAGIC population with custom design?</summary>
  <p><code>cped <- cbind(1:106, c(rep(0,4), 1, 3, rep(5,100)), c(rep(0,4), 2, 4, rep(6,100)), c(rep(0,4), rep(1,2), rep(2,100)))</code></p>
  <p><code>mpop <- magic.eval(ped=cped)</code></p>
  <br>
</details>

<details>
  <summary>9. Which R version do I need?</summary>
  <p><code>magicdesign</code> is built in R version 4.0.3 and it should work with any newer version of R. I have not tested <code>magicdesign</code> in older versions so it may or may not work in older R versions.</p>
  <br>
</details>

<details>
  <summary>10. Which are the supported platforms?</summary>
  <p><code>magicdesign</code> is built in Windows, but it should work in Mac or Linux.</p>
  <br>
</details>

<details>
  <summary>11. What are the computer requirements?</summary>
  <p>Nothing specific, but a fast processor (>2 GHz) with decent amount of RAM (>4 Gb) would be nice. <code>magicdesign</code> itself does not require more than one processor core, but I am not sure if its dependency <code>AlphaSimR</code> uses multiple cores.</p>
  <br>
</details>

## Contact us
You can reach me at [cyang@sruc.ac.uk](mailto:cyang@sruc.ac.uk) for any question or suggestion.

## Reference
Yang, C. J., R. N. Edmondson, H.-P. Piepho, W. Powell and I. Mackay, 2021. Crafting for a better MAGIC: systematic design and test for multiparental advanced generation inter-cross population. Biorxiv. [[Link](https://doi.org/10.1101/2021.04.27.441636)]
