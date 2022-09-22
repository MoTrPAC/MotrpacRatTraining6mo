# MotrpacRatTraining6mo 1.0.1 (2022-09-21)

Fix bugs:
* missing `:` in `transcript_timewise_da()` 
* counts improperly filtered in `transcript_prep_data()` 
* `DESeq2` prefix missing from `counts()` in `transcript_timewise_da()`

Site dev:
* keep examples from printing too much output
* rename `tutorial.Rmd` to `MotrpacRatTraining6mo` to take advantage of "Get Started" `pkgdown` feature
* add `@keywords internal` to functions that are not exported to prevent them from showing up in site index

# MotrpacRatTraining6mo 1.0.0 (2022-09-20)

First release
