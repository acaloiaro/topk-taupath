# Topk-taupath

The topk-taupath algorithm described in "The TopK Tau-Path Screen for Monotone Association": https://arxiv.org/pdf/1509.00549.pdf

Caveat: R FastBCS2 is not yet available
# About

The goal of this repository is to provide implementations of the top-k taupath algorithm in multiple langauges. The key algorithm _FastBCS_ (Fast Backward Conditional Search) will be provided in its original form, in addition to the much faster _FastBCS2_.

Many variable and function names refer to Joe Verducci and Srinath Sampath's original work lain out in their 2013 paper, _Detecting the end of agreement between two long ranked lists_. If a concept or variable name is unclear, their paper is a great resource https://onlinelibrary.wiley.com/doi/abs/10.1002/sam.11205.

# Current Status

As of 08-17-2019, the code here is very much a rough work in progress. Documentation is either poor or non-existent. Code is not well formatted or clearly lain out. Over time, I will add more langauges, implementations, and clean up my work.

The current goal is to finish FastBCS2 R.

# Todo
* Make R _FastBCS_ ✅
* Make R _FastBCS2_ available 🛑
* Make Java _FastBCS_ and _FastBCS2_ available ✅
* Make Go _FastBCS_ and _FastBCS2_ available ✅

