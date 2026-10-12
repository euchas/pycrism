# Parameters implemented in pycrism

Generated from the registry in `src/pycrism/parameters.py` (`pycrism list` prints the same).
Bands are nominal centre wavelengths in nm; numbers in brackets are the Table 2 median-kernel
widths in channels (6.55 nm each). `band depth` means BD = 1 - R_c / (a R_s + b R_l) with the
lever-rule weights of Viviano-Beck et al. (2014); `shoulder height` is the inverted form;
`RB` is a relative band depth and `dropoff` a continuum-removed ratio against a straight
continuum through two distant anchors.

The last column says what each definition was checked against:

- **Table 2**: wavelengths, weights and kernels as printed in Viviano-Beck et al. (2014), Table 2.
- **Table 2, CAT v7.4 hyperspectral form**: the printed table and the CRISM Analysis Toolkit
  differ for hyperspectral data; the CAT form is used and the difference is noted in the docstring.
- **method description only**: follows the published method but uses every valid channel in the
  stated ranges instead of the specific channel lists of Table 2.
- **adapted from the Curiosity-team conventions**: not in Viviano-Beck et al. (2014); band
  positions moved off the VNIR dropped-band gaps, so not identical to the rover parameters.

## VNIR (S detector), Viviano-Beck et al. (2014)

| Name | Kind | Bands [kernels], nm | Use | Checked against |
|---|---|---|---|---|
| `R440` | reflectance | 440 [5] | 0.44 um reflectance; TRU blue, clouds and haze | Table 2 |
| `R530` | reflectance | 530 [5] | 0.53 um reflectance; TRU green | Table 2 |
| `R600` | reflectance | 600 [5] | 0.60 um reflectance; TRU red | Table 2 |
| `R770` | reflectance | 770 [5] | 0.77 um reflectance; VNIR albedo | Table 2 |
| `RBR` | ratio | 770 / 440 [5, 5] | red/blue ratio; nanophase ferric oxide | Table 2 |
| `IRR1` | ratio | 800 / 997 [5, 5] | ice clouds versus dust | Table 2 |
| `BD530_2` | band depth | 440 / 530 / 614 [5, 5, 5] | crystalline hematite | Table 2 |
| `SH600_2` | shoulder height | 533 / 600 / 716 [5, 5, 3] | ferric minerals, textured surfaces | Table 2 |
| `SH770` | shoulder height | 716 / 775 / 860 [3, 5, 5] | ferric minerals | Table 2 |
| `BD640_2` | band depth | 600 / 624 / 760 [5, 3, 5] | maghemite | Table 2 |
| `BD860_2` | band depth | 755 / 860 / 977 [5, 5, 5] | crystalline hematite | Table 2 |
| `BD920_2` | band depth | 807 / 920 / 984 [5, 5, 5] | crystalline ferric minerals, low-Ca pyroxene | Table 2 |
| `RPEAK1` | polynomial peak | fit 440-1000 | wavelength of the VNIR reflectance peak | method description only |
| `BDI1000VIS` | integrated depth | 833-1023 against the fitted peak | 1 um band, VNIR side; olivine, pyroxene, glass | method description only |

## VNIR slopes and ratios (Curiosity conventions)

| Name | Kind | Bands [kernels], nm | Use | Checked against |
|---|---|---|---|---|
| `VIS_SLOPE` | slope | 442 -> 639 [5, 5] | visible slope; oxidation, dust | adapted from the Curiosity-team conventions |
| `S600_840` | slope | 600 -> 840 [5, 5] | red-to-NIR slope; ferric trend | adapted from the Curiosity-team conventions |
| `S750_840` | slope | 750 -> 840 [5, 5] | NIR downturn; red crystalline hematite negative | adapted from the Curiosity-team conventions |
| `RED_BLUE` | ratio | 639 / 527 [5, 5] | red/blue oxidation ratio | adapted from the Curiosity-team conventions |
| `RATIO_600_840` | ratio | 600 / 840 [5, 5] | oxidation ratio | adapted from the Curiosity-team conventions |
| `NIR_RATIO` | ratio | 1013 / 751 [5, 5] | ferric (above 1) versus ferrous (below 1) | adapted from the Curiosity-team conventions |
| `BD535` | band depth | 500 / 535 / 600 [5, 5, 5] | ferric band | adapted from the Curiosity-team conventions |
| `BD867W` | band depth | 751 / 867 / 1012 [5, 5, 5] | crystalline hematite, wide continuum | adapted from the Curiosity-team conventions |
| `BD867N` | band depth | 805 / 867 / 908 [5, 5, 5] | crystalline hematite, narrow continuum | adapted from the Curiosity-team conventions |

## IR (L detector), Viviano-Beck et al. (2014), Table 2 order

| Name | Kind | Bands [kernels], nm | Use | Checked against |
|---|---|---|---|---|
| `BDI1000IR` | integrated depth | 1020-1255, continuum 1330-1875 / 2430-2600 | 1 um band, IR side; Fe2+ silicates | Table 2, CAT v7.4 hyperspectral form |
| `R1330` | reflectance | 1330 [11] | IR albedo | Table 2 |
| `BD1300` | band depth | 1080 / 1320 / 1750 [5, 15, 5] | Fe2+ plagioclase | Table 2 |
| `OLINDEX3` | weighted RB | 11 bands 1080-1470, continuum 1750-2400 [7] | olivine | Table 2 |
| `LCPINDEX2` | weighted RB | 1690-1870, continuum 1560-2450 [7] | low-Ca pyroxene | Table 2 |
| `HCPINDEX2` | weighted RB | 2120-2460, continuum 1810-2530 [7] | high-Ca pyroxene | Table 2 |
| `VAR` | variance | line 1014-2287, residuals 1021-2253 | spectral variance 1.0-2.3 um | Table 2, CAT v7.4 hyperspectral form |
| `ISLOPE1` | slope | 1815 -> 2530 [5, 5], per um | ferric coatings on dark rock | Table 2 |
| `BD1400` | band depth | 1330 / 1395 / 1467 [5, 3, 5] | H2O, structural OH | Table 2 |
| `BD1435` | band depth | 1370 / 1435 / 1470 [3, 1, 3] | CO2 ice, some hydrated minerals | Table 2 |
| `BD1500_2` | band depth | 1367 / 1525 / 1808 [5, 11, 5] | H2O ice | Table 2, CAT v7.4 hyperspectral form |
| `ICER1_2` | ratio of depths | BD1435, BD1500_2 | CO2 versus H2O ice | Table 2 |
| `BD1750_2` | band depth | 1690 / 1750 / 1815 [5, 3, 5] | alunite, gypsum | Table 2 |
| `BD1900_2` | band depth | 1850 / 1930 / 2067 [5, 5, 5] | bound H2O | Table 2, CAT v7.4 hyperspectral form |
| `BD1900r2` | ratio form | 1908-1941 over 1862-1875 and 2112-2126, line 1850-2060 | bound H2O | Table 2 |
| `BDI2000` | integrated depth | 1660-2390, continuum 1330-1875 / 2430-2600 | 2 um band; pyroxene | Table 2, CAT v7.4 hyperspectral form |
| `BD2100_2` | band depth | 1930 / 2132 / 2250 [3, 5, 3] | monohydrated sulfates | Table 2 |
| `BD2165` | band depth | 2120 / 2165 / 2230 [5, 3, 3] | kaolinite group | Table 2 |
| `BD2190` | band depth | 2120 / 2185 / 2250 [5, 3, 3] | beidellite, allophane, imogolite | Table 2 |
| `MIN2200` | min of depths | 2165 and 2210, shoulders 2120 / 2350 [5, 3, 3] | kaolinite group | Table 2 |
| `BD2210_2` | band depth | 2165 / 2210 / 2250 [5, 5, 5] | Al-OH | Table 2 |
| `D2200` | dropoff | 2210, 2230 over 2165, line 1815-2430 | Al-OH | Table 2 |
| `BD2230` | band depth | 2210 / 2235 / 2252 [3, 3, 3] | hydroxylated ferric sulfate | Table 2 |
| `BD2250` | band depth | 2120 / 2245 / 2340 [5, 7, 3] | opal | Table 2 |
| `MIN2250` | min of depths | 2210 and 2265, shoulders 2165 / 2350 [5, 3, 5] | opal doublet | Table 2 |
| `BD2265` | band depth | 2210 / 2265 / 2340 [5, 3, 5] | jarosite, gibbsite | Table 2 |
| `BD2290` | band depth | 2250 / 2290 / 2350 [5, 5, 5] | Mg,Fe-OH; CO2 ice | Table 2 |
| `D2300` | dropoff | 2290, 2320, 2330 over 2120, 2170, 2210, line 1815-2530 | Fe/Mg phyllosilicates | Table 2 |
| `BD2355` | band depth | 2300 / 2355 / 2450 [5, 5, 5] | chlorite, prehnite, carbonate | Table 2 |
| `SINDEX2` | shoulder height | 2120 / 2290 / 2400 [5, 7, 3] | hydrated sulfates | Table 2 |
| `ICER2_2` | RB | 2600, line 2456-2530 [5] | CO2 ice positive, H2O ice negative | Table 2 |
| `MIN2295_2480` | min of depths | 2165 / 2295 / 2364 and 2364 / 2480 / 2570 [5, 5, 5] | Mg carbonate | Table 2 |
| `MIN2345_2537` | min of depths | 2250 / 2345 / 2430 and 2430 / 2537 / 2602 [5, 5, 5] | Ca/Fe carbonate | Table 2 |
| `BD2500_2` | band depth | 2364 / 2480 / 2570 [5, 5, 5] | Mg carbonate | Table 2 |
| `BD3000` | band depth | 3000 against R2530 * R2530 / R2210 [5] | bound H2O | Table 2 |
| `BD3100` | band depth | 3000 / 3120 / 3250 [5, 5, 5] | H2O ice | Table 2 |
| `BD3200` | band depth | 3250 / 3320 / 3390 [5, 5, 5] | CO2 ice | Table 2 |
| `BD3400_2` | band depth | 3250 / 3420 / 3630 [10, 15, 10] | carbonate | Table 2 |
| `CINDEX2` | shoulder height | 3450 / 3610 / 3875 [9, 11, 7] | carbonate | Table 2 |
| `R1080` | reflectance | 1080 [5] | FAL blue | Table 2 |
| `R1506` | reflectance | 1506 [5] | FAL green | Table 2 |
| `R2529` | reflectance | 2529 [5] | FAL red | Table 2 |
| `BD2600` | band depth | 2530 / 2600 / 2630 [5, 5, 5] | H2O vapour | Table 2 |
| `IRR2` | ratio | 2530 / 2210 [5, 5] | ice clouds versus dust | Table 2 |
| `IRR3` | ratio | 3500 / 3390 [7, 7] | ice clouds versus dust | Table 2 |
| `R3920` | reflectance | 3920 [5] | IC2 red | Table 2 |

## Browse recipes (Table 4)

| Recipe | Name | Red | Green | Blue |
|---|---|---|---|---|
| `TRU` | true colour | `R600` | `R530` | `R440` |
| `VNA` | VNIR albedo | `R770` | `R770` | `R770` |
| `FEM` | ferric minerals | `BD530_2` | `SH600_2` | `BDI1000VIS` |
| `FM2` | ferric minerals 2 | `BD530_2` | `BD920_2` | `BDI1000VIS` |
| `FAL` | false colour | `R2529` | `R1506` | `R1080` |
| `MAF` | mafic minerals | `OLINDEX3` | `LCPINDEX2` | `HCPINDEX2` |
| `HYD` | hydrated minerals | `SINDEX2` | `BD2100_2` | `BD1900_2` |
| `PHY` | phyllosilicates | `D2300` | `D2200` | `BD1900r2` |
| `PFM` | Fe/Mg phyllosilicates | `BD2355` | `D2300` | `BD2290` |
| `PAL` | Al phyllosilicates | `BD2210_2` | `BD2190` | `BD2165` |
| `HYS` | hydrated silica | `MIN2250` | `BD2250` | `BD1900r2` |
| `ICE` | ices | `BD1900_2` | `BD1500_2` | `BD1435` |
| `IC2` | ices 2 | `R3920` | `BD1500_2` | `BD1435` |
| `CHL` | chlorides | `ISLOPE1` | `BD3000` | `IRR2` |
| `CAR` | carbonates | `D2300` | `BD2500_2` | `BD1900_2` |
| `CR2` | carbonates 2 | `MIN2295_2480` | `MIN2345_2537` | `CINDEX2` |
