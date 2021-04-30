---
layout: post
title: magicdesignee
---

## Introduction
Welcome to :sparkles: `magicdesignee` :sparkles:!

`magicdesignee` is a Shiny app that runs `magicdesign` in its backend. Minimal experience in R is required.

* Local version is available by running `magicdesign::magicdesignee()` in `R`.
* Web version is available [here](https://magicdesign.shinyapps.io/magicdesignee/).

## Frequently asked questions (FAQs)
<details>
  <summary>1. How to use <code>magicdesignee</code> locally?</summary>
  <p>You will need to first install <code>magicdesign</code> from <a href="https://cjyang-sruc/github/magicdesign">here</a>, and then run <code>magicdesign::magicdesignee()</code> in <code>R</code>. The Shiny app should appear in your default browser.</p>
  <br>
</details>

<details>
  <summary>2. How to use <code>magicdesignee</code> on the web?</summary>
  <p>You can access <code>magicdesignee</code> <a href="https://magicdesign.shinyapps.io/magicdesignee">here</a>.</p>
  <br>
</details>

<details>
  <summary>3. Why is the web version not working for me?</summary>
  <p>There are several possibilities: the design is too much for the Shiny app server, the Shiny app has run out of its 25 hours of free monthly allowance, or the Shiny app has encountered a bug. It can be hard to identify the issue, so please <a href="mailto:cyang@sruc.ac.uk">report</a> it to me along with its error messages.</p>
  <br>
</details>

<details>
  <summary>4. How to fix the titles and texts in the plot axes that are hard to read?</summary>
  <p>You can either reduce the width of the browser window or download the plot using the button provided on top of the plot. The downloaded plot should have clearer font sizes.</p>
  <br>
</details>

<details>
  <summary>5. Why does the pedigree plot look incomplete?</summary>
  <p>The replicates from the last crossing generation are not displayed. It is just a mean to minimize clutter and hopefully easier to see the pedigree.</p>
  <br>
</details>

<details>
  <summary>6. Why is the analysis taking so long?</summary>
  <p>The computation time depends on the designs and number of simulations. It can be slow with many crosses and large populations, for example, more than 1,000 total crosses/individuals. It is recommended to start with small number of simulations (like 10) if you are going to use a design with large number of crosses and individuals. If everything looks fine, then just let the Shiny app takes its time if you want to use a higher number of simulations. Please consider using the local version because the web version is slow and can be unstable depending on your network connection.</p>
  <br>
</details>

## Contact us
You can reach me at [cyang@sruc.ac.uk](mailto:cyang@sruc.ac.uk) for any question or suggestion.

## Reference
Yang, C. J., R. Edmondson, H.-P. Piepho, W. Powell and I. Mackay, 2021. Crafting for a better MAGIC: systematic design and test for multiparental advanced generation inter-cross population. Biorxiv. [[Link](https://doi.org/10.1101/2021.04.27.441636)]
