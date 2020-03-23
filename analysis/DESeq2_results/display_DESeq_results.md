DESeq2 Drug Perturbation Results
================
Matthew Berginski
3/23/2020

``` r
DESeq2_summary = read_csv(here('DESeq2_summary.csv'))
```

    ## Parsed with column specification:
    ## cols(
    ##   treatment = col_character(),
    ##   batch = col_double(),
    ##   Genes = col_double(),
    ##   Kinases = col_double(),
    ##   `Dark Kinases` = col_double(),
    ##   `IDG Kinases` = col_character()
    ## )

``` r
kable(DESeq2_summary)
```

<table>

<thead>

<tr>

<th style="text-align:left;">

treatment

</th>

<th style="text-align:right;">

batch

</th>

<th style="text-align:right;">

Genes

</th>

<th style="text-align:right;">

Kinases

</th>

<th style="text-align:right;">

Dark Kinases

</th>

<th style="text-align:left;">

IDG Kinases

</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

Dinaciclib

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

303

</td>

<td style="text-align:right;">

9

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:left;">

ALPK2 | BCKDK | PRKACB

</td>

</tr>

<tr>

<td style="text-align:left;">

Entinostat

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

3508

</td>

<td style="text-align:right;">

132

</td>

<td style="text-align:right;">

30

</td>

<td style="text-align:left;">

ALPK3 | BRSK2 | CAMK1D | CAMK1G | CAMKV | CDK15 | CDK18 | ICK | LMTK3 |
MAP3K15 | MAPK4 | MAST3 | MKNK1 | MKNK2 | OBSCN | PAK6 | PDIK1L |
PIK3C2B | PKMYT1 | PKN3 | POMK | PRAG1 | PRKCQ | PSKH1 | SBK1 | STK17A |
STK32B | TESK2 | TLK1 | WNK2

</td>

</tr>

<tr>

<td style="text-align:left;">

JIB-04 (low dose)

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

13

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:left;">

NA

</td>

</tr>

<tr>

<td style="text-align:left;">

MK2206

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

2546

</td>

<td style="text-align:right;">

94

</td>

<td style="text-align:right;">

16

</td>

<td style="text-align:left;">

ALPK2 | CDK19 | CDKL1 | HIPK3 | MAP3K10 | NEK11 | NRBP2 | PIP5K1C |
PKMYT1 | PSKH1 | RIOK1 | RIOK2 | RIOK3 | SBK3 | STK38L | TESK1

</td>

</tr>

<tr>

<td style="text-align:left;">

Trametinib

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

3102

</td>

<td style="text-align:right;">

105

</td>

<td style="text-align:right;">

23

</td>

<td style="text-align:left;">

ALPK2 | ALPK3 | CDC42BPG | CDK15 | CDK18 | CDKL1 | CDKL5 | LMTK3 | LRRK1
| MAP3K21 | MAST3 | MAST4 | NEK6 | NRBP2 | PAK6 | PIK3C2B | PKMYT1 |
PKN3 | RIOK1 | RNASEL | STK17B | STK38L | TTBK2

</td>

</tr>

<tr>

<td style="text-align:left;">

JIB-04 (high dose)

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

6

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:left;">

NA

</td>

</tr>

<tr>

<td style="text-align:left;">

Navitoclax

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

2296

</td>

<td style="text-align:right;">

74

</td>

<td style="text-align:right;">

19

</td>

<td style="text-align:left;">

BRSK1 | CAMK1G | CDK20 | CLK3 | CLK4 | DYRK3 | DYRK4 | ICK | MAST3 |
NEK9 | NRBP2 | NUAK2 | PIK3C2B | PIP5K1C | PKMYT1 | PKN3 | RIOK3 | RSKR
| TAOK2

</td>

</tr>

<tr>

<td style="text-align:left;">

Paclitaxel

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

111

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:left;">

EEF2K

</td>

</tr>

<tr>

<td style="text-align:left;">

Ruxolibinib

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

472

</td>

<td style="text-align:right;">

18

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

CDC42BPA | CDK10 | EEF2K | NRBP2

</td>

</tr>

<tr>

<td style="text-align:left;">

SGCCP30

</td>

<td style="text-align:right;">

2

</td>

<td style="text-align:right;">

1382

</td>

<td style="text-align:right;">

56

</td>

<td style="text-align:right;">

8

</td>

<td style="text-align:left;">

CDK15 | ICK | MKNK2 | PKMYT1 | PRAG1 | PRKACB | STK19 | TESK2

</td>

</tr>

<tr>

<td style="text-align:left;">

Alisertib

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

22

</td>

<td style="text-align:right;">

1

</td>

<td style="text-align:right;">

0

</td>

<td style="text-align:left;">

NA

</td>

</tr>

<tr>

<td style="text-align:left;">

Bortezomib

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

3252

</td>

<td style="text-align:right;">

92

</td>

<td style="text-align:right;">

25

</td>

<td style="text-align:left;">

ADCK1 | CDC42BPG | CDK15 | CDKL1 | CLK3 | COQ8B | CSNK1G2 | ICK | LMTK2
| MAP3K14 | MAST4 | NEK1 | NEK6 | NEK9 | NRBP2 | PIP5K1A | PKN3 | PRAG1
| RIOK3 | STK17A | STK19 | STK32C | TESK1 | TLK2 | ULK4

</td>

</tr>

<tr>

<td style="text-align:left;">

CoCl2

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

1592

</td>

<td style="text-align:right;">

50

</td>

<td style="text-align:right;">

14

</td>

<td style="text-align:left;">

CDC42BPG | CDK15 | CDK18 | CDK19 | CLK3 | ICK | MKNK2 | NEK3 | NEK9 |
PKMYT1 | RIOK3 | RSKR | STK32A | TESK1

</td>

</tr>

<tr>

<td style="text-align:left;">

FCCP

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

178

</td>

<td style="text-align:right;">

7

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

ALPK2 | CDC42BPG | MARK1 | OBSCN

</td>

</tr>

<tr>

<td style="text-align:left;">

Palbociclib

</td>

<td style="text-align:right;">

3

</td>

<td style="text-align:right;">

941

</td>

<td style="text-align:right;">

40

</td>

<td style="text-align:right;">

4

</td>

<td style="text-align:left;">

NEK3 | PKMYT1 | RSKR | STK36

</td>

</tr>

</tbody>

</table>
