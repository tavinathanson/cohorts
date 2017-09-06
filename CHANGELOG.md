# Change Log

## [0.7.1](https://github.com/hammerlab/cohorts/tree/0.7.1) (2017-09-01)
[Full Changelog](https://github.com/hammerlab/cohorts/compare/0.7.0...0.7.1)

**Fixed bugs:**

- all\_effects caches, but does not return, all effects [\#253](https://github.com/hammerlab/cohorts/issues/253)
- Don't cache after only top-priority effects; post-filter instead [\#252](https://github.com/hammerlab/cohorts/issues/252)

**Merged pull requests:**

- Fix effects filtering for splice sites [\#256](https://github.com/hammerlab/cohorts/pull/256) ([tavinathanson](https://github.com/tavinathanson))

## [0.7.0](https://github.com/hammerlab/cohorts/tree/0.7.0) (2017-08-29)
[Full Changelog](https://github.com/hammerlab/cohorts/compare/0.6.4...0.7.0)

**Fixed bugs:**

- filter\_fn result caching should work for Python 2 and 3 [\#174](https://github.com/hammerlab/cohorts/issues/174)

**Closed issues:**

- ImportError: cannot import name 'AlleleParseError'  [\#233](https://github.com/hammerlab/cohorts/issues/233)
- Invalid RGBA argument error using `plot\_survival` [\#221](https://github.com/hammerlab/cohorts/issues/221)
- support plotting survival curves for more than 2 groups [\#217](https://github.com/hammerlab/cohorts/issues/217)
- Add `cache\_root\_dir` attribute [\#214](https://github.com/hammerlab/cohorts/issues/214)
- Cache sometimes fails if user interrupts writing process [\#170](https://github.com/hammerlab/cohorts/issues/170)
- builds against latest versions of isovar fails [\#69](https://github.com/hammerlab/cohorts/issues/69)

**Merged pull requests:**

- Fix effect priority caching [\#254](https://github.com/hammerlab/cohorts/pull/254) ([tavinathanson](https://github.com/tavinathanson))
- add get-blob method to gcio [\#244](https://github.com/hammerlab/cohorts/pull/244) ([jburos](https://github.com/jburos))
- add `join\_on\_left` parameter to dataframe\_loader [\#242](https://github.com/hammerlab/cohorts/pull/242) ([jburos](https://github.com/jburos))
- Expose gcloud utils [\#241](https://github.com/hammerlab/cohorts/pull/241) ([jburos](https://github.com/jburos))
- Plot surv by strata [\#237](https://github.com/hammerlab/cohorts/pull/237) ([jburos](https://github.com/jburos))
- allow median-vaf-purity with patients having no variants [\#236](https://github.com/hammerlab/cohorts/pull/236) ([jburos](https://github.com/jburos))
- Fix mhcnames to v 0.1.0  [\#235](https://github.com/hammerlab/cohorts/pull/235) ([jburos](https://github.com/jburos))
- "guess" at variant file format, if not obvious from filename [\#232](https://github.com/hammerlab/cohorts/pull/232) ([jburos](https://github.com/jburos))
- Update warning and associated comment [\#230](https://github.com/hammerlab/cohorts/pull/230) ([tavinathanson](https://github.com/tavinathanson))
- Gcloud batch downloads [\#229](https://github.com/hammerlab/cohorts/pull/229) ([jburos](https://github.com/jburos))
- Feature plot\_survival by category \(clean branch\) [\#225](https://github.com/hammerlab/cohorts/pull/225) ([jburos](https://github.com/jburos))
- Fix RGBA error with plot\_survival [\#223](https://github.com/hammerlab/cohorts/pull/223) ([jburos](https://github.com/jburos))
- Warn on vcf errors [\#222](https://github.com/hammerlab/cohorts/pull/222) ([jburos](https://github.com/jburos))
- Google Storage file access add-on for cohorts \(a.k.a. minibucket\) [\#219](https://github.com/hammerlab/cohorts/pull/219) ([armish](https://github.com/armish))

## [0.6.4](https://github.com/hammerlab/cohorts/tree/0.6.4) (2017-06-05)
[Full Changelog](https://github.com/hammerlab/cohorts/compare/0.6.3...0.6.4)

**Closed issues:**

- median\_vaf\_purity\(\) got an unexpected keyword argument 'normalized\_per\_mb' [\#212](https://github.com/hammerlab/cohorts/issues/212)

**Merged pull requests:**

- add cache\_root\_dir attribute & behavior [\#216](https://github.com/hammerlab/cohorts/pull/216) ([jburos](https://github.com/jburos))

## [0.6.3](https://github.com/hammerlab/cohorts/tree/0.6.3) (2017-05-28)
[Full Changelog](https://github.com/hammerlab/cohorts/compare/0.6.2...0.6.3)

**Closed issues:**

- Drop Python 2 support [\#204](https://github.com/hammerlab/cohorts/issues/204)

**Merged pull requests:**

- small fix to address \#212 [\#213](https://github.com/hammerlab/cohorts/pull/213) ([jburos](https://github.com/jburos))

## [0.6.2](https://github.com/hammerlab/cohorts/tree/0.6.2) (2017-05-19)
[Full Changelog](https://github.com/hammerlab/cohorts/compare/0.6.1...0.6.2)

## [0.6.1](https://github.com/hammerlab/cohorts/tree/0.6.1) (2017-05-19)
[Full Changelog](https://github.com/hammerlab/cohorts/compare/0.6.0...0.6.1)

**Merged pull requests:**

- Add new PyPi server [\#210](https://github.com/hammerlab/cohorts/pull/210) ([tavinathanson](https://github.com/tavinathanson))

## [0.6.0](https://github.com/hammerlab/cohorts/tree/0.6.0) (2017-05-18)
[Full Changelog](https://github.com/hammerlab/cohorts/compare/0.5.5...0.6.0)

**Closed issues:**

- Not identifying filtered-cache-name in python 3 [\#195](https://github.com/hammerlab/cohorts/issues/195)

**Merged pull requests:**

- Require Python 3 [\#208](https://github.com/hammerlab/cohorts/pull/208) ([tavinathanson](https://github.com/tavinathanson))
- Add optional show\_progress [\#206](https://github.com/hammerlab/cohorts/pull/206) ([tavinathanson](https://github.com/tavinathanson))

## [0.5.5](https://github.com/hammerlab/cohorts/tree/0.5.5) (2017-05-16)
[Full Changelog](https://github.com/hammerlab/cohorts/compare/0.5.4...0.5.5)

**Closed issues:**

- TypeError: "load\_maf\(\) got an unexpected keyword argument 'encoding'" [\#203](https://github.com/hammerlab/cohorts/issues/203)
- Solve intermittent test failures [\#201](https://github.com/hammerlab/cohorts/issues/201)

**Merged pull requests:**

- Update Pageant logic [\#205](https://github.com/hammerlab/cohorts/pull/205) ([tavinathanson](https://github.com/tavinathanson))

## [0.5.4](https://github.com/hammerlab/cohorts/tree/0.5.4) (2017-04-29)
[Full Changelog](https://github.com/hammerlab/cohorts/compare/0.5.3...0.5.4)

**Merged pull requests:**

- Update encoding and str handling for Py2/3 compatibility [\#199](https://github.com/hammerlab/cohorts/pull/199) ([tavinathanson](https://github.com/tavinathanson))

## [0.5.3](https://github.com/hammerlab/cohorts/tree/0.5.3) (2017-04-05)
[Full Changelog](https://github.com/hammerlab/cohorts/compare/0.5.2...0.5.3)

**Implemented enhancements:**

- Add watermark for default filter\_fn on Cohort initialization [\#183](https://github.com/hammerlab/cohorts/issues/183)

**Fixed bugs:**

- deletion\_count is counting all indels [\#182](https://github.com/hammerlab/cohorts/issues/182)
- exonic\_silent\_snv\_count is confusing [\#181](https://github.com/hammerlab/cohorts/issues/181)

**Closed issues:**

- Don't allow PFS/OS values to be NaN? [\#193](https://github.com/hammerlab/cohorts/issues/193)
- Add frameshift count function [\#191](https://github.com/hammerlab/cohorts/issues/191)
- additional\_data should automatically create Patient attributes [\#190](https://github.com/hammerlab/cohorts/issues/190)
- Minor: old comment [\#184](https://github.com/hammerlab/cohorts/issues/184)
- Number of functions is ballooning [\#173](https://github.com/hammerlab/cohorts/issues/173)
- Support "merged" strelka VCF type [\#167](https://github.com/hammerlab/cohorts/issues/167)
- Using only Strelka tier 1 results in divide by zero sometimes [\#166](https://github.com/hammerlab/cohorts/issues/166)

**Merged pull requests:**

- Add patient \_\_str\_\_ [\#197](https://github.com/hammerlab/cohorts/pull/197) ([tavinathanson](https://github.com/tavinathanson))
- Fixing filter cache issues [\#196](https://github.com/hammerlab/cohorts/pull/196) ([tavinathanson](https://github.com/tavinathanson))
- Add a median VAF filter function and a filter\_fn watermark [\#189](https://github.com/hammerlab/cohorts/pull/189) ([tavinathanson](https://github.com/tavinathanson))
- Refactor and fix functions [\#186](https://github.com/hammerlab/cohorts/pull/186) ([tavinathanson](https://github.com/tavinathanson))

## [0.5.2](https://github.com/hammerlab/cohorts/tree/0.5.2) (2017-01-18)
[Full Changelog](https://github.com/hammerlab/cohorts/compare/0.5.1...0.5.2)

**Merged pull requests:**

- Grab VAFs from MAFs [\#179](https://github.com/hammerlab/cohorts/pull/179) ([tavinathanson](https://github.com/tavinathanson))

## [0.5.1](https://github.com/hammerlab/cohorts/tree/0.5.1) (2017-01-17)
[Full Changelog](https://github.com/hammerlab/cohorts/compare/0.5.0...0.5.1)

**Merged pull requests:**

- Upgrade isovar to the latest version in cohorts [\#177](https://github.com/hammerlab/cohorts/pull/177) ([tavinathanson](https://github.com/tavinathanson))
- Pre-merge exonic counts into develop [\#172](https://github.com/hammerlab/cohorts/pull/172) ([jburos](https://github.com/jburos))

## [0.5.0](https://github.com/hammerlab/cohorts/tree/0.5.0) (2016-12-17)
[Full Changelog](https://github.com/hammerlab/cohorts/compare/0.4.2...0.5.0)

**Implemented enhancements:**

- Create public-facing usage template [\#114](https://github.com/hammerlab/cohorts/issues/114)

**Merged pull requests:**

- Add more granular variant specification and other minor changes [\#165](https://github.com/hammerlab/cohorts/pull/165) ([tavinathanson](https://github.com/tavinathanson))

## [0.4.2](https://github.com/hammerlab/cohorts/tree/0.4.2) (2016-11-30)
[Full Changelog](https://github.com/hammerlab/cohorts/compare/0.4.1...0.4.2)

**Merged pull requests:**

- Random MHC prediction for random cohort [\#164](https://github.com/hammerlab/cohorts/pull/164) ([tavinathanson](https://github.com/tavinathanson))
- Add variant generation to random\_cohorts and more [\#161](https://github.com/hammerlab/cohorts/pull/161) ([tavinathanson](https://github.com/tavinathanson))

## [0.4.1](https://github.com/hammerlab/cohorts/tree/0.4.1) (2016-11-05)
[Full Changelog](https://github.com/hammerlab/cohorts/compare/0.4.0...0.4.1)

**Fixed bugs:**

- Single variant collections do not allow for metadata extraction [\#81](https://github.com/hammerlab/cohorts/issues/81)
- Default to not requiring a single-value on dict [\#153](https://github.com/hammerlab/cohorts/pull/153) ([tavinathanson](https://github.com/tavinathanson))

**Merged pull requests:**

- Add regplot as a plot\_correlation option [\#159](https://github.com/hammerlab/cohorts/pull/159) ([tavinathanson](https://github.com/tavinathanson))
- Simplify functions with decorators [\#158](https://github.com/hammerlab/cohorts/pull/158) ([tavinathanson](https://github.com/tavinathanson))

## [0.4.0](https://github.com/hammerlab/cohorts/tree/0.4.0) (2016-09-28)
[Full Changelog](https://github.com/hammerlab/cohorts/compare/0.3.0...0.4.0)

**Implemented enhancements:**

- Make API more consistent [\#150](https://github.com/hammerlab/cohorts/pull/150) ([tavinathanson](https://github.com/tavinathanson))

**Closed issues:**

- Update README w/ Tavi's bullet list of features and Jacki's TCGA demo [\#147](https://github.com/hammerlab/cohorts/issues/147)
- Show plots in README's usage examples [\#141](https://github.com/hammerlab/cohorts/issues/141)

**Merged pull requests:**

- Update pypi badge every 4 hours [\#151](https://github.com/hammerlab/cohorts/pull/151) ([arahuja](https://github.com/arahuja))
- Use conda to install pypandoc [\#149](https://github.com/hammerlab/cohorts/pull/149) ([arahuja](https://github.com/arahuja))
- Update README.md [\#148](https://github.com/hammerlab/cohorts/pull/148) ([jburos](https://github.com/jburos))

## [0.3.0](https://github.com/hammerlab/cohorts/tree/0.3.0) (2016-09-27)
[Full Changelog](https://github.com/hammerlab/cohorts/compare/0.2.0...0.3.0)

**Implemented enhancements:**

- Break up load.py [\#37](https://github.com/hammerlab/cohorts/issues/37)

**Closed issues:**

- pip install fails on v 0.2.0 [\#144](https://github.com/hammerlab/cohorts/issues/144)
- Plot survival curves with upper/lower bounds [\#65](https://github.com/hammerlab/cohorts/issues/65)

**Merged pull requests:**

- Ensure requirements and README.md is packaged [\#146](https://github.com/hammerlab/cohorts/pull/146) ([arahuja](https://github.com/arahuja))
- upgrade to varcode 0.5.1 and necessary deps [\#143](https://github.com/hammerlab/cohorts/pull/143) ([arahuja](https://github.com/arahuja))
- Add pypi badge [\#142](https://github.com/hammerlab/cohorts/pull/142) ([arahuja](https://github.com/arahuja))
- Fix versioneer settings [\#140](https://github.com/hammerlab/cohorts/pull/140) ([tavinathanson](https://github.com/tavinathanson))

## [0.2.0](https://github.com/hammerlab/cohorts/tree/0.2.0) (2016-09-21)
[Full Changelog](https://github.com/hammerlab/cohorts/compare/0.1.2...0.2.0)

**Implemented enhancements:**

- Split up load.py [\#138](https://github.com/hammerlab/cohorts/pull/138) ([tavinathanson](https://github.com/tavinathanson))
- Add barplot option to plot\_joint and change plot\_joint name to plot\_correlation [\#134](https://github.com/hammerlab/cohorts/pull/134) ([tavinathanson](https://github.com/tavinathanson))
- Support joining on DataFrameLoaders with shared columns [\#132](https://github.com/hammerlab/cohorts/pull/132) ([tavinathanson](https://github.com/tavinathanson))
- Make benefit name configurable [\#131](https://github.com/hammerlab/cohorts/pull/131) ([tavinathanson](https://github.com/tavinathanson))
- Add float\_str to round floats [\#128](https://github.com/hammerlab/cohorts/pull/128) ([tavinathanson](https://github.com/tavinathanson))
- Add default survival colors [\#126](https://github.com/hammerlab/cohorts/pull/126) ([tavinathanson](https://github.com/tavinathanson))
- Pass stat\_func through jointplot [\#125](https://github.com/hammerlab/cohorts/pull/125) ([tavinathanson](https://github.com/tavinathanson))

**Closed issues:**

- Summarize provenance on init [\#80](https://github.com/hammerlab/cohorts/issues/80)

**Merged pull requests:**

- Fix travis run [\#139](https://github.com/hammerlab/cohorts/pull/139) ([tavinathanson](https://github.com/tavinathanson))
- Add quickstart [\#136](https://github.com/hammerlab/cohorts/pull/136) ([jburos](https://github.com/jburos))
- Allow boolean\_value\_map to be passed in to plot\_benefit [\#133](https://github.com/hammerlab/cohorts/pull/133) ([tavinathanson](https://github.com/tavinathanson))
- Fix threshold formatting [\#130](https://github.com/hammerlab/cohorts/pull/130) ([e5c](https://github.com/e5c))
- Fix fishers plot ranges [\#129](https://github.com/hammerlab/cohorts/pull/129) ([arahuja](https://github.com/arahuja))
- Allow configurable labels to survival plots [\#127](https://github.com/hammerlab/cohorts/pull/127) ([arahuja](https://github.com/arahuja))
- Fix survival plots for non-threshold conditions [\#123](https://github.com/hammerlab/cohorts/pull/123) ([arahuja](https://github.com/arahuja))

## [0.1.2](https://github.com/hammerlab/cohorts/tree/0.1.2) (2016-08-16)
[Full Changelog](https://github.com/hammerlab/cohorts/compare/0.1.1...0.1.2)

**Implemented enhancements:**

- Add p-value indicator [\#112](https://github.com/hammerlab/cohorts/pull/112) ([tavinathanson](https://github.com/tavinathanson))
- Make provenance printing optional [\#86](https://github.com/hammerlab/cohorts/pull/86) ([tavinathanson](https://github.com/tavinathanson))
- Filter to expressed variants [\#67](https://github.com/hammerlab/cohorts/pull/67) ([tavinathanson](https://github.com/tavinathanson))

**Fixed bugs:**

- Start plots at 0 on the y-axis [\#110](https://github.com/hammerlab/cohorts/issues/110)
- Summarize\_provenance fails after a cache has been deleted [\#92](https://github.com/hammerlab/cohorts/issues/92)
- VariantCollection metadata is lost after filtering [\#66](https://github.com/hammerlab/cohorts/issues/66)
- Caching issues [\#26](https://github.com/hammerlab/cohorts/issues/26)
- Always use a boolean mask [\#118](https://github.com/hammerlab/cohorts/pull/118) ([tavinathanson](https://github.com/tavinathanson))
- Fix Fisher's exact plot for non-boolean columns [\#117](https://github.com/hammerlab/cohorts/pull/117) ([tavinathanson](https://github.com/tavinathanson))
- Fix travis failure [\#106](https://github.com/hammerlab/cohorts/pull/106) ([tavinathanson](https://github.com/tavinathanson))
- Fix expressed\_missense\_snv\_count for None filter\_fn [\#95](https://github.com/hammerlab/cohorts/pull/95) ([tavinathanson](https://github.com/tavinathanson))

**Closed issues:**

- too easy to over-ride cohort-level defaults [\#120](https://github.com/hammerlab/cohorts/issues/120)
- Statistical significance \* on mann whitney and fisher's exact plot  [\#109](https://github.com/hammerlab/cohorts/issues/109)
- Convert benefit labels [\#108](https://github.com/hammerlab/cohorts/issues/108)
- filter\_fn default works across variants and effects [\#100](https://github.com/hammerlab/cohorts/issues/100)
- strip\_column\_names warning when running tests [\#89](https://github.com/hammerlab/cohorts/issues/89)
- Default filter\_fn back to None [\#77](https://github.com/hammerlab/cohorts/issues/77)
- Spacing of new functions is off [\#73](https://github.com/hammerlab/cohorts/issues/73)
- unit tests for summarize\_provenance [\#71](https://github.com/hammerlab/cohorts/issues/71)
- travis builds failing with isovar==0.0.6 [\#70](https://github.com/hammerlab/cohorts/issues/70)
- Fix NaN situation [\#54](https://github.com/hammerlab/cohorts/issues/54)
- Travis python2.7 build failing [\#53](https://github.com/hammerlab/cohorts/issues/53)
- Handle null counts better [\#16](https://github.com/hammerlab/cohorts/issues/16)

**Merged pull requests:**

- Add more to logrank results object [\#121](https://github.com/hammerlab/cohorts/pull/121) ([tavinathanson](https://github.com/tavinathanson))
- Allow labeling of response/benefit labels [\#111](https://github.com/hammerlab/cohorts/pull/111) ([arahuja](https://github.com/arahuja))
- Return namedtuples in plot functions [\#104](https://github.com/hammerlab/cohorts/pull/104) ([tavinathanson](https://github.com/tavinathanson))
- Remove flier from boxplot on stripboxplot [\#103](https://github.com/hammerlab/cohorts/pull/103) ([arahuja](https://github.com/arahuja))
- Add ax argument for fishers exact plot [\#102](https://github.com/hammerlab/cohorts/pull/102) ([arahuja](https://github.com/arahuja))
- Fix varcode at 0.4.15 until API updates [\#99](https://github.com/hammerlab/cohorts/pull/99) ([arahuja](https://github.com/arahuja))
- add load\_kallisto\(\) [\#98](https://github.com/hammerlab/cohorts/pull/98) ([e5c](https://github.com/e5c))
- Fix warning messages on strip\_column\_names [\#94](https://github.com/hammerlab/cohorts/pull/94) ([jburos](https://github.com/jburos))
- Fix minor bug in summarize\_provenance [\#93](https://github.com/hammerlab/cohorts/pull/93) ([jburos](https://github.com/jburos))
- Add filter\_fn to Cohort [\#90](https://github.com/hammerlab/cohorts/pull/90) ([tavinathanson](https://github.com/tavinathanson))
- first draft of clean\_column\_names utility function [\#85](https://github.com/hammerlab/cohorts/pull/85) ([jburos](https://github.com/jburos))
- Report data sources [\#84](https://github.com/hammerlab/cohorts/pull/84) ([jburos](https://github.com/jburos))
- modify setup.py to reference requirements.txt [\#83](https://github.com/hammerlab/cohorts/pull/83) ([jburos](https://github.com/jburos))
- Add basic bootstrap AUC and coxph to cohorts [\#79](https://github.com/hammerlab/cohorts/pull/79) ([tavinathanson](https://github.com/tavinathanson))
- Use filter\_fn None as default [\#78](https://github.com/hammerlab/cohorts/pull/78) ([tavinathanson](https://github.com/tavinathanson))
- Add responder\_pfs\_equals\_os options [\#75](https://github.com/hammerlab/cohorts/pull/75) ([tavinathanson](https://github.com/tavinathanson))
- Summarize provenance [\#68](https://github.com/hammerlab/cohorts/pull/68) ([jburos](https://github.com/jburos))
- Fix threshold rounding [\#63](https://github.com/hammerlab/cohorts/pull/63) ([arahuja](https://github.com/arahuja))
- Display a rounded threshold [\#62](https://github.com/hammerlab/cohorts/pull/62) ([arahuja](https://github.com/arahuja))
- Updates to filtering [\#61](https://github.com/hammerlab/cohorts/pull/61) ([tavinathanson](https://github.com/tavinathanson))
- Add Pageant CoverageDepth output parsing [\#60](https://github.com/hammerlab/cohorts/pull/60) ([tavinathanson](https://github.com/tavinathanson))
- Add axes argument to roc plots [\#59](https://github.com/hammerlab/cohorts/pull/59) ([arahuja](https://github.com/arahuja))
- Pass along axes argument [\#57](https://github.com/hammerlab/cohorts/pull/57) ([arahuja](https://github.com/arahuja))
- Fix bugs with caching and NaN handling [\#56](https://github.com/hammerlab/cohorts/pull/56) ([tavinathanson](https://github.com/tavinathanson))
- Allow extra keyword arguments in functions [\#55](https://github.com/hammerlab/cohorts/pull/55) ([tavinathanson](https://github.com/tavinathanson))
- Return 0 for patients with no neoantigens [\#52](https://github.com/hammerlab/cohorts/pull/52) ([arahuja](https://github.com/arahuja))
- Update CHANGELOG.md for 0.1.1 release [\#51](https://github.com/hammerlab/cohorts/pull/51) ([arahuja](https://github.com/arahuja))

## [0.1.1](https://github.com/hammerlab/cohorts/tree/0.1.1) (2016-06-14)
[Full Changelog](https://github.com/hammerlab/cohorts/compare/0.1.0...0.1.1)

**Implemented enhancements:**

- Cache a dependency manifest [\#15](https://github.com/hammerlab/cohorts/issues/15)
- Refactor how functions work [\#36](https://github.com/hammerlab/cohorts/pull/36) ([tavinathanson](https://github.com/tavinathanson))
- Pass kwargs through to function [\#34](https://github.com/hammerlab/cohorts/pull/34) ([tavinathanson](https://github.com/tavinathanson))
- Clarify mannwhitney sidedness [\#27](https://github.com/hammerlab/cohorts/pull/27) ([tavinathanson](https://github.com/tavinathanson))

**Fixed bugs:**

- Fix versioneer releasing [\#20](https://github.com/hammerlab/cohorts/issues/20)
- Fix RNA BAM path [\#30](https://github.com/hammerlab/cohorts/pull/30) ([tavinathanson](https://github.com/tavinathanson))

**Closed issues:**

- additional data assumes patient id field [\#44](https://github.com/hammerlab/cohorts/issues/44)
- Is verify\_cache still needed? [\#25](https://github.com/hammerlab/cohorts/issues/25)
- Fix isovar git dependency [\#23](https://github.com/hammerlab/cohorts/issues/23)

**Merged pull requests:**

- Add mutect and strelka specific variant stats parsing [\#50](https://github.com/hammerlab/cohorts/pull/50) ([arahuja](https://github.com/arahuja))
- Add ability to filter variants [\#49](https://github.com/hammerlab/cohorts/pull/49) ([arahuja](https://github.com/arahuja))
- Refactor plot\_benefit [\#48](https://github.com/hammerlab/cohorts/pull/48) ([tavinathanson](https://github.com/tavinathanson))
- Upgrade varcode for frameshift error and effects df [\#47](https://github.com/hammerlab/cohorts/pull/47) ([arahuja](https://github.com/arahuja))
- Update building dataframe with additional data [\#45](https://github.com/hammerlab/cohorts/pull/45) ([arahuja](https://github.com/arahuja))
- Remove verify cache functions [\#42](https://github.com/hammerlab/cohorts/pull/42) ([arahuja](https://github.com/arahuja))
- Move polyphen data dump to Cohorts object [\#41](https://github.com/hammerlab/cohorts/pull/41) ([tavinathanson](https://github.com/tavinathanson))
- Save variant metadata after merging variants [\#39](https://github.com/hammerlab/cohorts/pull/39) ([arahuja](https://github.com/arahuja))
- Add loader for PolyPhen2 annotations [\#38](https://github.com/hammerlab/cohorts/pull/38) ([armish](https://github.com/armish))
- Update requirements [\#33](https://github.com/hammerlab/cohorts/pull/33) ([tavinathanson](https://github.com/tavinathanson))
- Add plotting for ROC curve [\#32](https://github.com/hammerlab/cohorts/pull/32) ([arahuja](https://github.com/arahuja))
- Add docstring for load variant functions [\#31](https://github.com/hammerlab/cohorts/pull/31) ([arahuja](https://github.com/arahuja))
- Update isovar setup link [\#29](https://github.com/hammerlab/cohorts/pull/29) ([arahuja](https://github.com/arahuja))
- Add basic provenance tracking [\#28](https://github.com/hammerlab/cohorts/pull/28) ([tavinathanson](https://github.com/tavinathanson))
- Update README.md [\#24](https://github.com/hammerlab/cohorts/pull/24) ([arahuja](https://github.com/arahuja))
- Fix dataframe loading bug [\#22](https://github.com/hammerlab/cohorts/pull/22) ([tavinathanson](https://github.com/tavinathanson))
- add CHANGELOG.md [\#21](https://github.com/hammerlab/cohorts/pull/21) ([arahuja](https://github.com/arahuja))

## [0.1.0](https://github.com/hammerlab/cohorts/tree/0.1.0) (2016-05-03)
**Implemented enhancements:**

- Refactor cohorts [\#14](https://github.com/hammerlab/cohorts/pull/14) ([tavinathanson](https://github.com/tavinathanson))

**Fixed bugs:**

- Possible issue: sample\_id as int vs str [\#1](https://github.com/hammerlab/cohorts/issues/1)

**Closed issues:**

- Warning on import [\#13](https://github.com/hammerlab/cohorts/issues/13)

**Merged pull requests:**

- Update versioning to use versioneer [\#12](https://github.com/hammerlab/cohorts/pull/12) ([arahuja](https://github.com/arahuja))
- Set version in accesible variable [\#10](https://github.com/hammerlab/cohorts/pull/10) ([arahuja](https://github.com/arahuja))
- Add tests, Travis, etc. [\#9](https://github.com/hammerlab/cohorts/pull/9) ([tavinathanson](https://github.com/tavinathanson))
- Adding comments to plotting functions [\#8](https://github.com/hammerlab/cohorts/pull/8) ([arahuja](https://github.com/arahuja))
- Fix KM plot axis [\#7](https://github.com/hammerlab/cohorts/pull/7) ([arahuja](https://github.com/arahuja))
- Add lifelines as a dependency [\#6](https://github.com/hammerlab/cohorts/pull/6) ([armish](https://github.com/armish))
- Add stripboxplot [\#5](https://github.com/hammerlab/cohorts/pull/5) ([tavinathanson](https://github.com/tavinathanson))
- Handle PFS better [\#4](https://github.com/hammerlab/cohorts/pull/4) ([tavinathanson](https://github.com/tavinathanson))
- Add easy plotting for cohorts [\#3](https://github.com/hammerlab/cohorts/pull/3) ([tavinathanson](https://github.com/tavinathanson))
- Add topiary to the requirements [\#2](https://github.com/hammerlab/cohorts/pull/2) ([armish](https://github.com/armish))



\* *This Change Log was automatically generated by [github_changelog_generator](https://github.com/skywinder/Github-Changelog-Generator)*