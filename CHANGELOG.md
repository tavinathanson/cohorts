# Change Log

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

**Fixed bugs:**

- Start plots at 0 on the y-axis [\#110](https://github.com/hammerlab/cohorts/issues/110)
- Summarize\_provenance fails after a cache has been deleted [\#92](https://github.com/hammerlab/cohorts/issues/92)
- VariantCollection metadata is lost after filtering [\#66](https://github.com/hammerlab/cohorts/issues/66)
- Caching issues [\#26](https://github.com/hammerlab/cohorts/issues/26)
- Always use a boolean mask [\#118](https://github.com/hammerlab/cohorts/pull/118) ([tavinathanson](https://github.com/tavinathanson))
- Fix Fisher's exact plot for non-boolean columns [\#117](https://github.com/hammerlab/cohorts/pull/117) ([tavinathanson](https://github.com/tavinathanson))
- Fix travis failure [\#106](https://github.com/hammerlab/cohorts/pull/106) ([tavinathanson](https://github.com/tavinathanson))

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

## [0.1.1](https://github.com/hammerlab/cohorts/tree/0.1.1) (2016-06-14)
[Full Changelog](https://github.com/hammerlab/cohorts/compare/0.1.0...0.1.1)

**Implemented enhancements:**

- Cache a dependency manifest [\#15](https://github.com/hammerlab/cohorts/issues/15)

**Fixed bugs:**

- Fix versioneer releasing [\#20](https://github.com/hammerlab/cohorts/issues/20)

**Closed issues:**

- additional data assumes patient id field [\#44](https://github.com/hammerlab/cohorts/issues/44)
- Is verify\_cache still needed? [\#25](https://github.com/hammerlab/cohorts/issues/25)
- Fix isovar git dependency [\#23](https://github.com/hammerlab/cohorts/issues/23)

## [0.1.0](https://github.com/hammerlab/cohorts/tree/0.1.0) (2016-05-03)
**Fixed bugs:**

- Possible issue: sample\_id as int vs str [\#1](https://github.com/hammerlab/cohorts/issues/1)

**Closed issues:**

- Warning on import [\#13](https://github.com/hammerlab/cohorts/issues/13)



\* *This Change Log was automatically generated by [github_changelog_generator](https://github.com/skywinder/Github-Changelog-Generator)*